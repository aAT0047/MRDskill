#!/usr/bin/env python3
"""QC a somatic VCF against the local MRD support graph evidence layer.

The script intentionally does not fabricate gene-level biological evidence when
the VCF has no gene or functional annotation. It still produces a technical QC
table and a KG-RAG readiness report.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import html
import json
import math
import re
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from statistics import median
from typing import Any

from mrdskill_collect_case_qc_outputs import collect_case_qc_outputs
from mrdskill_project_paths import artifact_path


MISSING = {
    "",
    ".",
    "na",
    "n/a",
    "none",
    "null",
    "unknown",
    "uncertain",
    "unavailable",
    "not_available",
    "needs_review",
    "not_provided",
    "not_reported",
    "not assessable",
    "not_assessable",
}

CHIP_REVIEW_GENES = {
    "DNMT3A",
    "TET2",
    "ASXL1",
    "PPM1D",
    "TP53",
    "JAK2",
    "SF3B1",
    "SRSF2",
    "U2AF1",
    "IDH1",
    "IDH2",
    "ATM",
    "CHEK2",
    "KRAS",
    "GNAS",
}

TARGETED_OR_RESISTANCE_REVIEW_GENES = {
    "EGFR",
    "KRAS",
    "NRAS",
    "BRAF",
    "ERBB2",
    "ALK",
    "ROS1",
    "RET",
    "MET",
    "NTRK1",
    "NTRK2",
    "NTRK3",
    "ESR1",
    "PIK3CA",
    "BRCA1",
    "BRCA2",
    "IDH1",
    "IDH2",
}

CALLER_DEPENDENT_SCENARIOS = [
    {
        "scenario": "tumor_informed_marker_tracking",
        "why_caller_matters": "Tracking known tumor mutations or personalized panel markers requires variant-level calls, not only a cancer-signal score.",
        "vcf_only_contribution": "Produces a molecular candidate that can be compared with the primary-tumor marker set when available.",
        "required_inputs_for_completion": ["primary_tumor_marker_set", "tumor_concordance_status"],
    },
    {
        "scenario": "clone_dynamics_or_recurrence_origin",
        "why_caller_matters": "Clone, gene, and variant identity are needed to review whether recurrence signal is local, metastatic, treatment-selected, or newly emergent.",
        "vcf_only_contribution": "Preserves variant identity for later trajectory and clone-origin comparison.",
        "required_inputs_for_completion": ["serial_ctdna_samples", "baseline_tumor_clone_map", "recurrence_or_metastasis_sample_if_available"],
    },
    {
        "scenario": "chip_or_hematopoietic_origin_audit",
        "why_caller_matters": "CHIP and hematopoietic-origin review requires variant-level evidence, matched WBC/buffy coat comparison, and VAF-pattern review.",
        "vcf_only_contribution": "Surfaces CHIP-review genes, matched-control trace signal, and variant-level VAF/depth patterns.",
        "required_inputs_for_completion": ["matched_wbc_or_buffy_coat", "germline_or_hematopoietic_filtering_status"],
    },
    {
        "scenario": "low_frequency_serial_confirmation",
        "why_caller_matters": "A single low-support or low-VAF call often cannot drive action; the same variant must be followed across repeat or serial samples.",
        "vcf_only_contribution": "Identifies low-support calls that require repeat sampling or orthogonal confirmation.",
        "required_inputs_for_completion": ["repeat_sample", "serial_ctdna_samples", "assay_lod"],
    },
    {
        "scenario": "targeted_resistance_or_trial_evidence",
        "why_caller_matters": "Therapy selection, resistance monitoring, and trial eligibility often require specific molecular alterations.",
        "vcf_only_contribution": "Surfaces target/resistance genes that need functional annotation and actionability review.",
        "required_inputs_for_completion": ["functional_variant_annotation", "therapy_context", "actionability_database_or_guideline"],
    },
    {
        "scenario": "fragmentomics_caller_discordance_audit",
        "why_caller_matters": "Discordance between fragmentomics cancer-signal scores and caller-derived variants requires auditable source reconciliation.",
        "vcf_only_contribution": "Provides the molecular side of the conflict map: variant identity, QC, KG provenance, and source-trace status.",
        "required_inputs_for_completion": ["fragmentomics_status_or_score", "fragmentomics_assay_metadata"],
    },
    {
        "scenario": "clinical_report_audit_trail",
        "why_caller_matters": "Clinical report review and regulatory audit need traceable molecular evidence, QC rationale, and expert-review flags.",
        "vcf_only_contribution": "Creates variant-level audit evidence: coordinate, gene, VAF, depth, QC flags, KG provenance, and review action.",
        "required_inputs_for_completion": ["expert_signout_policy", "clinical_context_for_final_report"],
    },
]

CALLER_SCENARIO_INFO = {item["scenario"]: item for item in CALLER_DEPENDENT_SCENARIOS}
CALLER_SCENARIO_ORDER = [item["scenario"] for item in CALLER_DEPENDENT_SCENARIOS]

RAG_ROUND_SCENARIO_MAP = {
    "support_retrieval": ["tumor_informed_marker_tracking", "clinical_report_audit_trail"],
    "variant_origin": ["chip_or_hematopoietic_origin_audit", "fragmentomics_caller_discordance_audit"],
    "assay_timepoint_context": [
        "tumor_informed_marker_tracking",
        "low_frequency_serial_confirmation",
        "fragmentomics_caller_discordance_audit",
    ],
    "longitudinal_clonal_context": [
        "clone_dynamics_or_recurrence_origin",
        "low_frequency_serial_confirmation",
        "tumor_informed_marker_tracking",
    ],
    "targeted_resistance_context": ["targeted_resistance_or_trial_evidence"],
    "contradiction_uncertainty": [
        "chip_or_hematopoietic_origin_audit",
        "fragmentomics_caller_discordance_audit",
        "low_frequency_serial_confirmation",
    ],
    "source_trace_protocol": ["clinical_report_audit_trail"],
}

RAG_CATEGORY_TERMS = {
    "tumor_informed": ["tumor informed", "tumor-informed", "primary tumor", "tumour tissue", "tumor tissue", "personalized"],
    "serial_or_clonal": ["serial", "longitudinal", "trajectory", "clearance", "persistence", "recurrence", "relapse", "clonal"],
    "chip_wbc": ["chip", "clonal hematopo", "wbc", "buffy coat", "hematopoietic", "haematopoietic", "germline"],
    "assay_lod": ["lod", "limit of detection", "assay", "panel", "low vaf", "low-vaf"],
    "artifact_or_confounder": ["artifact", "artefact", "inflammation", "necrosis", "radiotherapy", "immunotherapy", "background"],
    "therapy_or_resistance": ["resistance", "targeted", "therapy", "trial", "egfr", "kras", "braf", "met", "esr1", "brca"],
}

SOURCE_TRACE_STATUSES = [
    "traceable_caller_evidence_candidate",
    "traceable_caller_signal_needs_serial_confirmation",
    "hematopoietic_or_CHIP_origin_suspected",
    "technical_artifact_origin_suspected",
    "cross_sectional_signal_needs_contextualization",
    "technical_rescue_required_before_source_assignment",
    "untraceable_missing_annotation_or_context",
    "discordant_signal_origin_needs_reconciliation",
    "low_shedding_or_no_caller_signal_context_needed",
]

EVIDENCE_MODES = [
    "caller_source_trace",
    "signal_only_requires_molecular_trace",
    "discordant_signal_and_caller_trace",
    "vcf_only_source_trace_qc",
    "insufficient_context",
]

NCBI_FEATURE_TABLES = {
    "GRCh37": {
        "assembly": "GRCh37.p13",
        "synonym": "hg19",
        "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_feature_table.txt.gz",
        "filename": "GCF_000001405.25_GRCh37.p13_feature_table.txt.gz",
    },
    "GRCh38": {
        "assembly": "GRCh38.p14",
        "synonym": "hg38",
        "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_feature_table.txt.gz",
        "filename": "GCF_000001405.40_GRCh38.p14_feature_table.txt.gz",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run MRD KG-aware QC for a somatic VCF.")
    parser.add_argument("--vcf", default=None, help="Input VCF path for single-sample QC.")
    parser.add_argument("--case-manifest", default=None, help="CSV/TSV manifest for multi-file longitudinal MRD QC.")
    parser.add_argument("--longitudinal-variant-csv", default=None, help="Output per-variant longitudinal source-trace CSV.")
    parser.add_argument("--write-manifest-template", default=None, help="Write a blank case-manifest template CSV and exit unless --vcf/--case-manifest is also provided.")
    parser.add_argument("--evidence-csv", default="data/evidence/mrd_evidence_scored_full.csv", help="Confidence-scored MRD evidence CSV.")
    parser.add_argument("--variant-csv", default=None, help="Output per-variant QC CSV.")
    parser.add_argument("--summary-json", default=None, help="Output summary JSON.")
    parser.add_argument("--report-md", default=None, help="Output Markdown report.")
    parser.add_argument("--case-report-md", default=None, help="Output EMR-style case-level Markdown report.")
    parser.add_argument("--visual-html", default=None, help="Output HTML evidence-lane visualization.")
    parser.add_argument("--case-id", default=None, help="Case/sample identifier for outputs.")
    parser.add_argument("--ncbi-annotate", dest="ncbi_annotate", action="store_true", default=True, help="Annotate coordinates with NCBI RefSeq gene intervals; enabled by default.")
    parser.add_argument("--no-ncbi-annotate", dest="ncbi_annotate", action="store_false", help="Disable NCBI RefSeq gene-interval annotation.")
    parser.add_argument("--assembly", default="auto", choices=["auto", "GRCh37", "GRCh38"], help="Genome assembly for NCBI interval annotation.")
    parser.add_argument("--annotation-cache-dir", default="data/annotation/ncbi", help="Local cache directory for NCBI feature tables.")
    parser.add_argument("--source-trace-rag", action="store_true", help="Run multi-agent RAG retrieval for caller-signal source tracing.")
    parser.add_argument("--rag-chunks-jsonl", default="data/rag/mrd_latest_unified_rag_chunks.jsonl", help="Unified RAG chunks JSONL.")
    parser.add_argument("--rag-rule-db", default="data/rag/mrd_latest_kg_rule_index.sqlite", help="Structured rule SQLite index.")
    parser.add_argument("--rag-top-k-per-round", type=int, default=3, help="RAG literature hits per source-trace round.")
    parser.add_argument("--rag-rule-top-k-per-round", type=int, default=2, help="Structured rule hits per source-trace round.")
    parser.add_argument("--rag-max-chunks", type=int, default=120000, help="Cap RAG chunks for source-trace retrieval; set 0 for full scan.")
    parser.add_argument("--rag-full-scan", action="store_true", help="Scan the full RAG chunk index; equivalent to --rag-max-chunks 0.")
    parser.add_argument("--rag-mode", choices=["case", "variant", "complete"], default="case", help="RAG depth: case rounds only, variant-specific rounds, or complete GraphRAG-style synthesis.")
    parser.add_argument("--rag-variant-limit", type=int, default=24, help="Maximum variants selected for variant-specific RAG in variant/complete mode.")
    parser.add_argument("--rag-min-score", type=float, default=0.0, help="Minimum denoised RAG score to keep.")
    parser.add_argument("--rag-max-hits-per-variant", type=int, default=3, help="Maximum RAG hits retained per variant after case and variant-specific retrieval.")
    parser.add_argument("--rag-denoise", dest="rag_denoise", action="store_true", default=True, help="Apply source-tracing relevance filters to RAG retrieval.")
    parser.add_argument("--no-rag-denoise", dest="rag_denoise", action="store_false", help="Disable RAG source-tracing relevance filters.")
    parser.add_argument("--rag-json", default=None, help="Output source-trace RAG JSON.")
    parser.add_argument("--rag-md", default=None, help="Output source-trace RAG Markdown.")
    return parser.parse_args()


def open_text_maybe_gzip(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def clean_text(value: Any) -> str:
    if value is None:
        return ""
    return str(value).strip()


def is_missing(value: Any) -> bool:
    return clean_text(value).lower() in MISSING


def to_float(value: Any) -> float | None:
    text = clean_text(value)
    if not text or text == ".":
        return None
    if text.lower() in {"inf", "+inf", "infinity"}:
        return math.inf
    try:
        return float(text)
    except ValueError:
        return None


def to_int(value: Any) -> int | None:
    num = to_float(value)
    if num is None or math.isinf(num):
        return None
    return int(num)


def parse_info_field(info_text: str) -> dict[str, str]:
    info: dict[str, str] = {}
    for item in info_text.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = "true"
    return info


def parse_sample(format_keys: list[str], sample_text: str) -> dict[str, str]:
    values = sample_text.split(":")
    return {key: values[i] if i < len(values) else "" for i, key in enumerate(format_keys)}


def parse_pair(value: Any) -> tuple[int | None, int | None]:
    text = clean_text(value)
    if not text or "," not in text:
        return None, None
    left, right = text.split(",", 1)
    return to_int(left), to_int(right)


def parse_number_list(value: Any) -> list[float]:
    text = clean_text(value)
    if not text or text == ".":
        return []
    values: list[float] = []
    for part in text.split(","):
        parsed = to_float(part)
        if parsed is not None and not math.isinf(parsed):
            values.append(parsed)
    return values


def parse_int_list(value: Any) -> list[int]:
    ints: list[int] = []
    for num in parse_number_list(value):
        ints.append(int(num))
    return ints


def alt_depth_from_ad(value: Any) -> int | None:
    values = parse_int_list(value)
    if len(values) < 2:
        return None
    return max(values[1:])


def alt_strand_from_sb(value: Any) -> tuple[int | None, int | None]:
    values = parse_int_list(value)
    if len(values) < 4:
        return None, None
    return values[-2], values[-1]


def alt_info_number(info: dict[str, str], key: str) -> float | None:
    values = parse_number_list(info.get(key))
    if not values:
        return None
    if len(values) >= 2:
        return values[1]
    return values[0]


def normalize_chrom(chrom: str) -> str:
    text = clean_text(chrom).lower()
    if text.startswith("chr"):
        text = text[3:]
    if text == "m":
        return "mt"
    return text.upper()


def normalize_gene(gene: str) -> str:
    text = clean_text(gene)
    if is_missing(text):
        return ""
    text = re.sub(r"\s+", "", text)
    return text.upper()


def parse_vcf(vcf_path: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    samples: list[str] = []
    info_headers: dict[str, str] = {}
    format_headers: dict[str, str] = {}
    variants: list[dict[str, Any]] = []
    tumor_sample_votes: Counter[str] = Counter()
    normal_sample_votes: Counter[str] = Counter()

    with open_text_maybe_gzip(vcf_path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##GATKCommandLine=<ID=Mutect2"):
                match = re.search(r"--normal-sample\s+([^\s,>]+)", line)
                if match:
                    normal_sample_votes[match.group(1)] += 10
                continue
            if line.startswith("##INFO="):
                match = re.search(r"ID=([^,>]+).*?Description=\"([^\"]*)\"", line)
                if match:
                    info_headers[match.group(1)] = match.group(2)
                continue
            if line.startswith("##FORMAT="):
                match = re.search(r"ID=([^,>]+).*?Description=\"([^\"]*)\"", line)
                if match:
                    format_headers[match.group(1)] = match.group(2)
                continue
            if line.startswith("#CHROM"):
                parts = line.split("\t")
                samples = parts[9:]
                continue
            if line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, vid, ref, alt, qual, flt, info_text = parts[:8]
            fmt = parts[8] if len(parts) > 8 else ""
            sample_texts = parts[9:]
            info = parse_info_field(info_text)
            if info.get("SAMPLE"):
                tumor_sample_votes[info["SAMPLE"]] += 1

            format_keys = fmt.split(":") if fmt else []
            sample_values = {
                samples[i]: parse_sample(format_keys, sample_texts[i])
                for i in range(min(len(samples), len(sample_texts)))
            }

            variants.append(
                {
                    "chrom": chrom,
                    "pos": to_int(pos),
                    "id": vid,
                    "ref": ref,
                    "alt": alt,
                    "qual": to_float(qual),
                    "filter": flt,
                    "info": info,
                    "format_keys": format_keys,
                    "sample_values": sample_values,
                    "raw_line": line,
                }
            )

    tumor_sample = ""
    if tumor_sample_votes:
        tumor_sample = tumor_sample_votes.most_common(1)[0][0]
    elif normal_sample_votes and samples:
        normal_sample = normal_sample_votes.most_common(1)[0][0]
        tumor_candidates = [sample for sample in samples if sample != normal_sample]
        tumor_sample = tumor_candidates[-1] if tumor_candidates else samples[-1]
    elif samples:
        tumor_sample = samples[-1] if len(samples) > 1 else samples[0]
    control_samples = [sample for sample in samples if sample != tumor_sample]

    meta = {
        "samples": samples,
        "tumor_sample": tumor_sample,
        "control_samples": control_samples,
        "info_headers": info_headers,
        "format_headers": format_headers,
        "tumor_sample_votes": dict(tumor_sample_votes),
        "normal_sample_votes": dict(normal_sample_votes),
    }
    return variants, meta


def extract_genes(info: dict[str, str]) -> list[str]:
    genes: set[str] = set()
    direct_keys = [
        "GENE",
        "Gene",
        "gene",
        "SYMBOL",
        "Gene.refGene",
        "GeneDetail.refGene",
        "Hugo_Symbol",
        "GENEINFO",
    ]
    for key in direct_keys:
        value = info.get(key)
        if not value:
            continue
        for token in re.split(r"[,|;:/\s]+", value):
            gene = normalize_gene(token)
            if gene and not gene.isdigit():
                genes.add(gene)

    ann = info.get("ANN")
    if ann:
        for record in ann.split(","):
            fields = record.split("|")
            if len(fields) > 3:
                gene = normalize_gene(fields[3])
                if gene:
                    genes.add(gene)

    # VEP CSQ commonly carries the gene symbol in one of the early fields, but
    # the exact position depends on the header. Use a conservative symbol-like
    # extraction only when no direct annotation was found.
    csq = info.get("CSQ")
    if csq and not genes:
        for record in csq.split(","):
            fields = record.split("|")
            for field in fields[:8]:
                gene = normalize_gene(field)
                if re.match(r"^[A-Z][A-Z0-9.-]{1,15}$", gene) and not gene.startswith("ENSG"):
                    genes.add(gene)

    return sorted(genes)


def ensure_ncbi_feature_table(assembly: str, cache_dir: Path) -> Path:
    spec = NCBI_FEATURE_TABLES[assembly]
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / spec["filename"]
    if path.exists() and path.stat().st_size > 0:
        return path

    tmp = path.with_suffix(path.suffix + ".tmp")
    print(f"Downloading NCBI {spec['assembly']} feature table to {path} ...")
    with urllib.request.urlopen(spec["url"], timeout=120) as response, tmp.open("wb") as out:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(path)
    return path


def parse_gene_id(feature_id: str) -> str:
    text = clean_text(feature_id)
    match = re.search(r"GeneID:(\d+)", text)
    if match:
        return match.group(1)
    return text


def chrom_sort_key(chrom: str) -> tuple[int, str]:
    text = normalize_chrom(chrom)
    if text.isdigit():
        return int(text), text
    order = {"X": 23, "Y": 24, "MT": 25}
    return order.get(text, 99), text


def load_ncbi_gene_intervals(feature_table_path: Path, assembly: str) -> dict[str, list[dict[str, Any]]]:
    cache_json = feature_table_path.with_suffix("").with_suffix(f".{assembly}.gene_intervals.json")
    if cache_json.exists() and cache_json.stat().st_mtime >= feature_table_path.stat().st_mtime:
        data = json.loads(cache_json.read_text(encoding="utf-8"))
        return {chrom: intervals for chrom, intervals in data.items()}

    intervals: dict[str, list[dict[str, Any]]] = defaultdict(list)
    with gzip.open(feature_table_path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("# feature") != "gene":
                continue
            chrom = normalize_chrom(row.get("chromosome", ""))
            if not chrom:
                continue
            start = to_int(row.get("start"))
            end = to_int(row.get("end"))
            if start is None or end is None:
                continue
            if start > end:
                start, end = end, start
            symbol = normalize_gene(row.get("symbol") or row.get("name") or "")
            if not symbol:
                continue
            intervals[chrom].append(
                {
                    "symbol": symbol,
                    "gene_id": parse_gene_id(row.get("GeneID") or row.get("feature_id", "")),
                    "gene_type": clean_text(row.get("class")) or clean_text(row.get("feature_type")),
                    "start": start,
                    "end": end,
                    "strand": clean_text(row.get("strand")),
                    "assembly": assembly,
                }
            )

    for chrom in list(intervals):
        intervals[chrom].sort(key=lambda item: (item["start"], item["end"], item["symbol"]))
    cache_json.write_text(json.dumps(intervals, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    return dict(intervals)


def annotate_coordinate_with_intervals(
    chrom: str,
    pos: int | None,
    intervals_by_chrom: dict[str, list[dict[str, Any]]],
    flank_bp: int = 0,
) -> list[dict[str, Any]]:
    if pos is None:
        return []
    intervals = intervals_by_chrom.get(normalize_chrom(chrom), [])
    matches: list[dict[str, Any]] = []
    query_start = pos - flank_bp
    query_end = pos + flank_bp
    for interval in intervals:
        if interval["start"] > query_end:
            break
        if interval["end"] < query_start:
            continue
        matches.append(interval)
    return matches


def annotate_variants_with_ncbi(
    variants: list[dict[str, Any]],
    requested_assembly: str,
    cache_dir: Path,
) -> dict[str, Any]:
    assemblies = ["GRCh37", "GRCh38"] if requested_assembly == "auto" else [requested_assembly]
    candidate_indexes: dict[str, dict[str, list[dict[str, Any]]]] = {}
    candidate_stats: dict[str, dict[str, Any]] = {}

    for assembly in assemblies:
        feature_table = ensure_ncbi_feature_table(assembly, cache_dir)
        intervals_by_chrom = load_ncbi_gene_intervals(feature_table, assembly)
        exact_count = 0
        near_count = 0
        sentinel_hits: list[dict[str, Any]] = []
        for variant in variants:
            exact = annotate_coordinate_with_intervals(variant["chrom"], variant["pos"], intervals_by_chrom, flank_bp=0)
            near = exact or annotate_coordinate_with_intervals(variant["chrom"], variant["pos"], intervals_by_chrom, flank_bp=5000)
            if exact:
                exact_count += 1
            if near:
                near_count += 1
            if len(sentinel_hits) < 8 and (exact or near):
                sentinel_hits.append(
                    {
                        "variant_id": f"{variant['chrom']}:{variant['pos']}:{variant['ref']}>{variant['alt']}",
                        "exact_genes": [item["symbol"] for item in exact],
                        "nearby_genes": [item["symbol"] for item in near],
                    }
                )
        candidate_indexes[assembly] = intervals_by_chrom
        candidate_stats[assembly] = {
            "exact_gene_hits": exact_count,
            "nearby_gene_hits_5kb": near_count,
            "feature_table": str(feature_table),
            "gene_interval_count": sum(len(items) for items in intervals_by_chrom.values()),
            "sentinel_hits": sentinel_hits,
        }

    selected_assembly = assemblies[0]
    if requested_assembly == "auto" and len(assemblies) > 1:
        selected_assembly = max(
            assemblies,
            key=lambda asm: (
                candidate_stats[asm]["exact_gene_hits"],
                candidate_stats[asm]["nearby_gene_hits_5kb"],
            ),
        )

    intervals_by_chrom = candidate_indexes[selected_assembly]
    per_variant: dict[str, dict[str, Any]] = {}
    for variant in variants:
        variant_id = f"{variant['chrom']}:{variant['pos']}:{variant['ref']}>{variant['alt']}"
        exact = annotate_coordinate_with_intervals(variant["chrom"], variant["pos"], intervals_by_chrom, flank_bp=0)
        nearby = [] if exact else annotate_coordinate_with_intervals(variant["chrom"], variant["pos"], intervals_by_chrom, flank_bp=5000)
        per_variant[variant_id] = {
            "assembly": selected_assembly,
            "genes": exact,
            "nearby_genes_5kb": nearby,
            "annotation_source": "NCBI RefSeq feature_table gene intervals",
            "annotation_status": "exact_gene_interval" if exact else ("nearby_gene_5kb" if nearby else "intergenic_or_unresolved"),
        }

    return {
        "enabled": True,
        "requested_assembly": requested_assembly,
        "selected_assembly": selected_assembly,
        "source": "NCBI RefSeq feature_table",
        "source_url": NCBI_FEATURE_TABLES[selected_assembly]["url"],
        "candidate_stats": candidate_stats,
        "per_variant": per_variant,
    }


def extract_coordinate_key(text: str) -> tuple[str, int] | None:
    value = clean_text(text)
    if is_missing(value):
        return None
    match = re.search(r"(?:chr)?([0-9]{1,2}|X|Y|M|MT)\D{1,8}([0-9][0-9,]{2,})", value, re.IGNORECASE)
    if not match:
        return None
    chrom = normalize_chrom(match.group(1))
    pos = to_int(match.group(2).replace(",", ""))
    if pos is None:
        return None
    return chrom, pos


def add_counter(counter: Counter[str], value: Any, max_items: int = 4) -> None:
    text = clean_text(value)
    if is_missing(text):
        return
    items = re.split(r"[,;/|]+", text)
    for item in items[:max_items]:
        cleaned = clean_text(item)
        if not is_missing(cleaned):
            counter[cleaned] += 1


def load_evidence_context(evidence_csv: Path) -> dict[str, Any]:
    context: dict[str, Any] = {
        "evidence_csv": str(evidence_csv),
        "evidence_rows": 0,
        "confidence_counts": Counter(),
        "graph_inclusion_counts": Counter(),
        "audit_priority_counts": Counter(),
        "gene_counts": Counter(),
        "assay_counts": Counter(),
        "sample_type_counts": Counter(),
        "cancer_type_counts": Counter(),
        "mrd_related_rows": 0,
        "chip_related_rows": 0,
        "false_positive_rows": 0,
        "matched_wbc_rows": 0,
        "germline_filter_rows": 0,
        "chip_filter_rows": 0,
        "coordinate_index": defaultdict(list),
        "gene_index": defaultdict(list),
    }

    if not evidence_csv.exists():
        context["missing"] = True
        return context

    text_fields = [
        "title",
        "quoted_or_paraphrased_evidence",
        "mrd_relevance",
        "relapse_association",
        "recurrence_association",
        "prognosis_association",
        "chip_association",
        "false_positive_risk",
        "assay_type",
        "sample_type",
    ]

    with evidence_csv.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            context["evidence_rows"] += 1
            evidence_id = clean_text(row.get("evidence_id"))
            context["confidence_counts"][clean_text(row.get("confidence_score")) or "not_reported"] += 1
            context["graph_inclusion_counts"][clean_text(row.get("graph_inclusion_status")) or "not_reported"] += 1
            context["audit_priority_counts"][clean_text(row.get("audit_priority")) or "not_reported"] += 1

            gene = normalize_gene(row.get("gene_symbol_normalized") or row.get("gene") or "")
            evidence_brief = {
                "evidence_id": evidence_id,
                "gene": gene,
                "confidence_score": clean_text(row.get("confidence_score")),
                "graph_inclusion_status": clean_text(row.get("graph_inclusion_status")),
                "audit_priority": clean_text(row.get("audit_priority")),
                "title": clean_text(row.get("title")),
                "pmid": clean_text(row.get("pmid")),
                "doi": clean_text(row.get("doi")),
                "cancer_type": clean_text(row.get("cancer_type")),
                "assay_type": clean_text(row.get("assay_type")),
                "sample_type": clean_text(row.get("sample_type")),
                "mrd_relevance": clean_text(row.get("mrd_relevance")),
                "chip_association": clean_text(row.get("chip_association")),
                "false_positive_risk": clean_text(row.get("false_positive_risk")),
                "matched_WBC_used": clean_text(row.get("matched_WBC_used")),
                "CHIP_filtering_used": clean_text(row.get("CHIP_filtering_used")),
                "germline_filtering_used": clean_text(row.get("germline_filtering_used")),
                "evidence_text": clean_text(row.get("quoted_or_paraphrased_evidence")),
            }
            if gene:
                context["gene_counts"][gene] += 1
                if evidence_id and len(context["gene_index"][gene]) < 8:
                    context["gene_index"][gene].append(evidence_brief)

            add_counter(context["assay_counts"], row.get("assay_type"))
            add_counter(context["sample_type_counts"], row.get("sample_type"))
            add_counter(context["cancer_type_counts"], row.get("cancer_type"))

            combined = " ".join(clean_text(row.get(field)) for field in text_fields).lower()
            if any(term in combined for term in ["mrd", "ctdna", "circulating tumour dna", "circulating tumor dna", "minimal residual"]):
                context["mrd_related_rows"] += 1
            if "chip" in combined or "clonal hematopo" in combined:
                context["chip_related_rows"] += 1
            if "false positive" in combined or "artifact" in combined or "artefact" in combined:
                context["false_positive_rows"] += 1

            if not is_missing(row.get("matched_WBC_used")):
                context["matched_wbc_rows"] += 1
            if not is_missing(row.get("germline_filtering_used")):
                context["germline_filter_rows"] += 1
            if not is_missing(row.get("CHIP_filtering_used")):
                context["chip_filter_rows"] += 1

            coord_key = extract_coordinate_key(row.get("genomic_coordinate", ""))
            if coord_key and evidence_id and len(context["coordinate_index"][coord_key]) < 8:
                context["coordinate_index"][coord_key].append(evidence_brief)

    return context


def sortable_counter(counter: Counter[str], limit: int = 15) -> list[dict[str, Any]]:
    return [{"name": key, "count": count} for key, count in counter.most_common(limit)]


def sort_by_numeric(rows: list[dict[str, Any]], key: str, reverse: bool = True) -> list[dict[str, Any]]:
    def value(row: dict[str, Any]) -> float:
        raw = row.get(key)
        if isinstance(raw, (int, float)) and not math.isinf(raw):
            return float(raw)
        parsed = to_float(raw)
        if parsed is None or math.isinf(parsed):
            return float("-inf")
        return parsed

    return sorted(rows, key=value, reverse=reverse)


def make_variant_brief(row: dict[str, Any]) -> dict[str, Any]:
    return {
        "variant_id": row.get("variant_id"),
        "genes": row.get("genes"),
        "type": row.get("type"),
        "tumor_dp": row.get("tumor_dp"),
        "tumor_vd": row.get("tumor_vd"),
        "tumor_af": row.get("tumor_af"),
        "control_af": row.get("control_af"),
        "technical_score": row.get("technical_score"),
        "qc_tier": row.get("qc_tier"),
        "mrd_entry_role": row.get("mrd_entry_role"),
        "kg_match_status": row.get("kg_match_status"),
        "kg_evidence_ids": row.get("kg_evidence_ids"),
        "kg_provenance_summary": row.get("kg_provenance_summary"),
        "kg_evidence_details": row.get("kg_evidence_details"),
        "source_trace_status": row.get("source_trace_status"),
        "evidence_mode": row.get("evidence_mode"),
        "caller_dependent_mrd_scenarios": row.get("caller_dependent_mrd_scenarios"),
        "caller_workflow_summary": row.get("caller_workflow_summary"),
        "caller_required_context_to_complete": row.get("caller_required_context_to_complete"),
        "caller_scenario_completion": row.get("caller_scenario_completion"),
        "recommended_review_action": row.get("recommended_review_action"),
        "variant_nlp_qc_conclusion": row.get("variant_nlp_qc_conclusion"),
        "source_trace_reason_codes": row.get("source_trace_reason_codes"),
        "inner_qc_conclusion": row.get("inner_qc_conclusion"),
        "rag_evidence_status": row.get("rag_evidence_status"),
        "rag_agent_summary": row.get("rag_agent_summary"),
        "rag_support_summary": row.get("rag_support_summary"),
        "rag_supporting_evidence": row.get("rag_supporting_evidence"),
        "rag_uncertainty_checklist": row.get("rag_uncertainty_checklist"),
        "rag_evidence_details": row.get("rag_evidence_details"),
        "graphrag_path_summary": row.get("graphrag_path_summary"),
        "source_trace_synthesis": row.get("source_trace_synthesis"),
        "source_trace_explanation": row.get("source_trace_explanation"),
        "supporting_evidence_summary": row.get("supporting_evidence_summary"),
        "contradictory_evidence_summary": row.get("contradictory_evidence_summary"),
        "missing_context_summary": row.get("missing_context_summary"),
        "next_confirmation_step": row.get("next_confirmation_step"),
        "source_assignment_confidence": row.get("source_assignment_confidence"),
        "qc_flags": row.get("qc_flags"),
    }


def has_flag(row: dict[str, Any], flag: str) -> bool:
    return flag in set(clean_text(row.get("qc_flags")).split(";"))


def split_semicolon(value: Any) -> list[str]:
    return [item for item in clean_text(value).split(";") if item]


def build_caller_workflow_coverage(rows: list[dict[str, Any]], rag_result: dict[str, Any] | None = None) -> list[dict[str, Any]]:
    scenario_counts: Counter[str] = Counter()
    scenario_examples: defaultdict[str, list[str]] = defaultdict(list)
    for row in rows:
        variant_id = clean_text(row.get("variant_id"))
        for scenario in split_semicolon(row.get("caller_dependent_mrd_scenarios")):
            scenario_counts[scenario] += 1
            if variant_id and len(scenario_examples[scenario]) < 6:
                label = variant_id
                genes = clean_text(row.get("genes"))
                if genes:
                    label = f"{label} ({genes.split(';')[0]})"
                scenario_examples[scenario].append(label)

    rag_scenario_counts: Counter[str] = Counter()
    if rag_result:
        for round_result in rag_result.get("rounds", []):
            round_id = clean_text(round_result.get("round_id"))
            mapped = RAG_ROUND_SCENARIO_MAP.get(round_id, [])
            hit_count = len(round_result.get("hits", []) or []) + len(round_result.get("rule_hits", []) or [])
            for scenario in mapped:
                rag_scenario_counts[scenario] += hit_count

    coverage: list[dict[str, Any]] = []
    for scenario in CALLER_SCENARIO_ORDER:
        info = CALLER_SCENARIO_INFO[scenario]
        variant_count = scenario_counts.get(scenario, 0)
        if not variant_count and not rag_scenario_counts.get(scenario, 0):
            status = "not_triggered_by_current_vcf"
        elif variant_count and rag_scenario_counts.get(scenario, 0):
            status = "vcf_layer_plus_rag_context_available"
        elif variant_count:
            status = "vcf_layer_available_rag_not_integrated"
        else:
            status = "rag_context_available_no_current_variant_trigger"
        coverage.append(
            {
                "scenario": scenario,
                "status": status,
                "variant_count": variant_count,
                "example_variants": scenario_examples.get(scenario, []),
                "rag_hit_count": rag_scenario_counts.get(scenario, 0),
                "vcf_only_contribution": info.get("vcf_only_contribution", ""),
                "why_caller_matters": info.get("why_caller_matters", ""),
                "required_inputs_for_completion": info.get("required_inputs_for_completion", []),
            }
        )
    return coverage


def build_molecular_evidence_qc(summary: dict[str, Any], rows: list[dict[str, Any]], gate: dict[str, Any]) -> dict[str, Any]:
    status_counts = Counter(clean_text(row.get("source_trace_status")) or "untraceable_missing_annotation_or_context" for row in rows)
    mode_counts = Counter(clean_text(row.get("evidence_mode")) or "insufficient_context" for row in rows)
    scenario_counts: Counter[str] = Counter()
    action_counts: Counter[str] = Counter()
    for row in rows:
        for scenario in split_semicolon(row.get("caller_dependent_mrd_scenarios")):
            scenario_counts[scenario] += 1
        action = clean_text(row.get("recommended_review_action"))
        if action:
            action_counts[action] += 1

    traceable = status_counts.get("traceable_caller_evidence_candidate", 0)
    serial_needed = status_counts.get("traceable_caller_signal_needs_serial_confirmation", 0)
    chip_like = status_counts.get("hematopoietic_or_CHIP_origin_suspected", 0)
    artifact_like = status_counts.get("technical_artifact_origin_suspected", 0)
    rescue_needed = status_counts.get("technical_rescue_required_before_source_assignment", 0)
    contextualize_needed = status_counts.get("cross_sectional_signal_needs_contextualization", 0)
    total = max(1, len(rows))

    if not rows:
        source_trace_status = "untraceable_missing_annotation_or_context"
        evidence_mode = "insufficient_context"
    elif summary.get("annotated_variant_count", 0) == 0:
        source_trace_status = "untraceable_missing_annotation_or_context"
        evidence_mode = "vcf_only_source_trace_qc"
    elif traceable > 0 and gate.get("entry_permission") == "allow":
        source_trace_status = "traceable_caller_evidence_candidate"
        evidence_mode = "caller_source_trace"
    elif chip_like > 0 and chip_like >= traceable + serial_needed:
        source_trace_status = "hematopoietic_or_CHIP_origin_suspected"
        evidence_mode = "caller_source_trace"
    elif artifact_like + rescue_needed > total * 0.75:
        source_trace_status = "technical_rescue_required_before_source_assignment"
        evidence_mode = "vcf_only_source_trace_qc"
    elif traceable > 0 or serial_needed > 0:
        source_trace_status = "traceable_caller_signal_needs_serial_confirmation"
        evidence_mode = "caller_source_trace"
    else:
        source_trace_status = "untraceable_missing_annotation_or_context"
        evidence_mode = mode_counts.most_common(1)[0][0] if mode_counts else "insufficient_context"

    actions: list[str] = []
    if summary.get("annotated_variant_count", 0) == 0:
        actions.append("add_gene_and_functional_annotation")
    if summary.get("matched_control_signal_count", 0) > 0:
        actions.append("review_matched_control_or_wbc_signal")
    if gate.get("chip_review_variant_count", 0) > 0:
        actions.append("matched_wbc_or_buffy_coat_chip_review")
    if gate.get("low_support_variant_count", 0) > 0:
        actions.append("repeat_sampling_or_orthogonal_confirmation_for_low_support_calls")
    if gate.get("artifact_review_variant_count", 0) > 0:
        actions.append("local_alignment_and_caller_artifact_review")
    if not actions:
        actions.append("expert_mrd_review_with_clinical_timepoint_and_tumor_concordance")

    def row_rank(row: dict[str, Any]) -> tuple[int, float, float]:
        return (
            {
                "traceable_caller_evidence_candidate": 0,
                "traceable_caller_signal_needs_serial_confirmation": 1,
                "hematopoietic_or_CHIP_origin_suspected": 2,
                "untraceable_missing_annotation_or_context": 3,
                "technical_artifact_origin_suspected": 4,
                "cross_sectional_signal_needs_contextualization": 4,
                "technical_rescue_required_before_source_assignment": 5,
            }.get(clean_text(row.get("source_trace_status")), 6),
            -(row.get("technical_score") or 0),
            -(row.get("tumor_af") or 0),
        )

    selected_rows: list[dict[str, Any]] = []
    selected_ids: set[str] = set()

    def add_rows(candidates: list[dict[str, Any]], limit: int) -> None:
        for row in sorted(candidates, key=row_rank):
            vid = clean_text(row.get("variant_id"))
            if vid and vid not in selected_ids:
                selected_rows.append(row)
                selected_ids.add(vid)
            if len([item for item in selected_rows if item in candidates]) >= limit:
                break

    kg_rows = [row for row in rows if clean_text(row.get("kg_match_status")) in {"gene_evidence_match", "coordinate_evidence_match"}]
    chip_rows = [row for row in rows if "chip_gene_overlap" in clean_text(row.get("source_trace_reason_codes"))]
    strong_rows = [row for row in rows if clean_text(row.get("qc_tier")) == "strong_technical_candidate"]
    artifact_rows = [row for row in rows if clean_text(row.get("source_trace_status")) == "technical_artifact_origin_suspected"]
    add_rows(kg_rows, 8)
    add_rows(chip_rows, 4)
    add_rows(strong_rows, 8)
    add_rows(artifact_rows, 4)
    if len(selected_rows) < 16:
        add_rows(rows, 16 - len(selected_rows))

    if source_trace_status == "untraceable_missing_annotation_or_context" and summary.get("annotated_variant_count", 0) == 0:
        conclusion = (
            "This VCF contains caller-derived molecular signals, but their source cannot be traced because gene/functional "
            "annotation and clinical context are incomplete. Use this output as a source-tracing QC packet: add annotation, "
            "matched-control or WBC review, tumor-concordance evidence, and serial sampling context before interpreting MRD origin."
        )
    elif source_trace_status == "traceable_caller_evidence_candidate":
        conclusion = (
            "The caller-derived signal contains technically strong, KG-supported molecular evidence that can be traced at gene/variant level. "
            "The next question is origin confirmation: tumor concordance, timepoint, assay LOD, WBC status, and serial trajectory."
        )
    elif source_trace_status == "traceable_caller_signal_needs_serial_confirmation":
        conclusion = (
            "The caller-derived signal has traceable molecular candidates, but the source explanation is incomplete. "
            "Use repeat sampling, orthogonal confirmation, tumor concordance, and missing-context review to determine whether the signal is persistent, transient, hematopoietic, or technical."
        )
    elif source_trace_status == "hematopoietic_or_CHIP_origin_suspected":
        conclusion = (
            "The dominant source-tracing explanation is hematopoietic/CHIP origin. Matched WBC or buffy-coat review and longitudinal behavior are required "
            "before attributing these caller signals to tumor MRD."
        )
    elif source_trace_status == "technical_rescue_required_before_source_assignment":
        conclusion = (
            "This is a cross-sectional caller-derived signal set with many technical rescue priorities. The KG evidence is still useful for source tracing, "
            "but current data should be interpreted as provisional: prioritize manual caller rescue, local alignment review, orthogonal confirmation, and context collection rather than rejecting the signal."
        )
    elif source_trace_status == "technical_artifact_origin_suspected":
        conclusion = (
            "Several caller-derived signals carry artifact-prone features. The current source explanation remains provisional and should be separated into technical review, KG provenance review, and clinical-context follow-up."
        )
    else:
        conclusion = (
            "The caller-derived signal cannot be traced to a reliable biological source in its current form. Additional annotation, technical review, "
            "tumor/WBC concordance, and disease-context metadata are required."
        )

    return {
        "qc_question": "What is the source-tracing explanation for the caller-derived signal in this MRD setting?",
        "nlp_qc_conclusion": conclusion,
        "source_trace_status": source_trace_status,
        "evidence_mode": evidence_mode,
        "available_source_trace_statuses": SOURCE_TRACE_STATUSES,
        "available_evidence_modes": EVIDENCE_MODES,
        "caller_dependent_mrd_scenarios": CALLER_DEPENDENT_SCENARIOS,
        "caller_workflow_coverage": build_caller_workflow_coverage(rows),
        "source_trace_status_counts": dict(status_counts),
        "evidence_mode_counts": dict(mode_counts),
        "caller_scenario_counts": dict(scenario_counts),
        "recommended_review_actions": actions,
        "recommended_review_action_counts": dict(action_counts),
        "kg_matched_variant_count": len(kg_rows),
        "chip_review_variant_count": len(chip_rows),
        "top_variant_conclusions": [make_variant_brief(row) for row in selected_rows[:20]],
    }


def build_mrd_entry_gate(summary: dict[str, Any], rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Assess source-tracing readiness for VCF-derived molecular calls."""
    primary_rows = [row for row in rows if row.get("mrd_entry_role") == "primary_mrd_adjudication_candidate"]
    conditional_rows = [row for row in rows if row.get("mrd_entry_role") == "conditional_mrd_adjudication_candidate"]
    review_rows = [
        row
        for row in rows
        if row.get("mrd_entry_role")
        in {"review_before_mrd_adjudication", "technical_candidate_needs_biological_context"}
    ]
    technical_rescue_rows = [row for row in rows if row.get("mrd_entry_role") == "do_not_enter_mrd_adjudication_without_manual_rescue"]
    chip_review_rows = [row for row in rows if has_flag(row, "common_chip_gene_review")]
    low_support_rows = [
        row
        for row in rows
        if has_flag(row, "variant_support_low") or has_flag(row, "variant_support_very_low")
    ]
    artifact_rows = [
        row
        for row in rows
        if any(
            has_flag(row, flag)
            for flag in [
                "repeat_or_shift_context",
                "alt_strand_imbalance_review",
                "severe_alt_strand_imbalance",
                "indel_or_complex_variant",
                "mapping_quality_low",
                "ssf_elevated",
            ]
        )
    ]

    reasons: list[str] = []
    blockers: list[str] = []
    guards: list[str] = []

    if summary.get("record_count", 0) == 0:
        blockers.append("No VCF variant records were parsed.")
    if summary.get("matched_control_signal_count", 0) > 0:
        guards.append("Matched/control alternate signal is present and must be reviewed before MRD use.")
    if summary.get("annotated_variant_count", 0) == 0:
        guards.append("No gene annotation is available; direct biological plausibility review is limited.")
    if summary.get("gene_match_variant_count", 0) == 0 and summary.get("coordinate_match_variant_count", 0) == 0:
        guards.append("No direct local KG support was found for current gene/coordinate calls.")
    if not primary_rows and not conditional_rows and not review_rows:
        blockers.append("No variant reached a technical or KG-supported source-tracing candidate role.")
    if len(technical_rescue_rows) == len(rows) and rows:
        blockers.append("All parsed variants need technical rescue before case-specific source assignment.")
    if chip_review_rows:
        guards.append("One or more candidates are in common CHIP-risk genes and need hematopoietic-origin review.")
    if low_support_rows:
        guards.append("Low variant-support records are present; repeat/orthogonal confirmation may be needed.")
    if artifact_rows:
        guards.append("Artifact-prone records are present; local alignment and caller-level review are needed.")

    if primary_rows:
        reasons.append("At least one technically strong, KG-supported candidate is available.")
    if conditional_rows:
        reasons.append("At least one technically strong candidate has gene annotation but limited direct KG support.")
    if review_rows:
        reasons.append("Review-level candidates are present but require expert QC before MRD interpretation.")
    if summary.get("matched_control_signal_count", 0) == 0:
        reasons.append("No matched/control alternate signal was observed in the parsed control fields.")
    if summary.get("annotated_variant_count", 0) > 0:
        reasons.append("At least one variant has gene-level annotation for KG lookup.")

    if blockers:
        gate_class = "NOT_READY_FOR_SOURCE_ASSIGNMENT"
        entry_permission = "technical_rescue_required"
        verdict = "This VCF is not ready for case-specific source assignment until technical rescue and missing context review are completed."
    elif primary_rows and not any("Matched/control" in item for item in guards) and not chip_review_rows:
        gate_class = "READY_FOR_EXPERT_SOURCE_TRACING"
        entry_permission = "allow"
        verdict = "This VCF can enter expert molecular source tracing as a candidate caller-evidence packet."
    else:
        gate_class = "REVIEW_REQUIRED_BEFORE_SOURCE_ASSIGNMENT"
        entry_permission = "conditional"
        verdict = "This VCF can be triaged, but source assignment should remain provisional until the guardrail items are reviewed."

    final_signout_requirements = [
        "cancer_type",
        "disease_stage",
        "sample_type",
        "collection_timepoint",
        "treatment_setting",
        "assay_panel_and_lod",
        "matched_wbc_or_buffy_coat_status",
        "tumor_tissue_confirmation_or_prior_tumor_variants",
        "serial_ctdna_trend_or_repeat_sample_when_low_level",
    ]

    return {
        "entry_question": "Can this VCF-derived molecular-call packet enter expert source tracing?",
        "gate_class": gate_class,
        "entry_permission": entry_permission,
        "verdict": verdict,
        "allowed_scope": "source-tracing readiness triage; not final MRD positivity or treatment recommendation",
        "primary_candidate_count": len(primary_rows),
        "conditional_candidate_count": len(conditional_rows),
        "review_candidate_count": len(review_rows),
        "blocked_variant_count": len(technical_rescue_rows),
        "technical_rescue_variant_count": len(technical_rescue_rows),
        "chip_review_variant_count": len(chip_review_rows),
        "low_support_variant_count": len(low_support_rows),
        "artifact_review_variant_count": len(artifact_rows),
        "reasons": reasons,
        "guardrails": guards,
        "blockers": blockers,
        "missing_for_final_signout": final_signout_requirements,
        "top_entry_candidates": [make_variant_brief(row) for row in sort_by_numeric(primary_rows + conditional_rows + review_rows, "tumor_af")[:12]],
    }


def build_case_assessment(summary: dict[str, Any], rows: list[dict[str, Any]]) -> dict[str, Any]:
    strong_rows = [row for row in rows if row.get("qc_tier") == "strong_technical_candidate"]
    review_rows = [row for row in rows if row.get("qc_tier") == "review_candidate"]
    weak_rows = [row for row in rows if row.get("qc_tier") == "weak_or_artifact_prone_candidate"]
    hold_rows = [row for row in rows if row.get("qc_tier") == "hold_or_reject_without_review"]
    kg_rows = [row for row in rows if row.get("kg_match_status") == "gene_evidence_match"]
    high_af_rows = [row for row in rows if (row.get("tumor_af") or 0) >= 0.10]
    low_support_rows = [
        row
        for row in rows
        if "variant_support_low" in clean_text(row.get("qc_flags"))
        or "variant_support_very_low" in clean_text(row.get("qc_flags"))
    ]
    artifact_rows = [
        row
        for row in rows
        if any(
            flag in clean_text(row.get("qc_flags"))
            for flag in [
                "repeat_or_shift_context",
                "alt_strand_imbalance_review",
                "severe_alt_strand_imbalance",
                "indel_or_complex_variant",
                "mapping_quality_low",
                "ssf_elevated",
            ]
        )
    ]

    total = max(1, len(rows))
    strong_fraction = len(strong_rows) / total
    hold_fraction = len(hold_rows) / total
    review_or_worse_fraction = (len(review_rows) + len(weak_rows) + len(hold_rows)) / total

    if summary.get("matched_control_signal_count", 0) > 0:
        sample_qc_class = "REVIEW_REQUIRED_CONTROL_SIGNAL_PRESENT"
        overall_interpretation = "matched/control sample shows non-zero alternate signal, so case-level use requires manual review."
    elif len(strong_rows) >= 1 and len(hold_rows) <= 5 and strong_fraction >= 0.45:
        sample_qc_class = "PASS_WITH_REVIEW_NOTES"
        overall_interpretation = "VCF has a substantial set of technically strong somatic candidates with a clean matched/control sample; weaker records still need review before MRD use."
    elif len(strong_rows) >= 1:
        sample_qc_class = "BORDERLINE_REVIEW_REQUIRED"
        overall_interpretation = "VCF contains usable technical candidates, but the review burden is high enough to require case-level manual assessment."
    else:
        sample_qc_class = "FAIL_OR_INSUFFICIENT_FOR_MRD_USE"
        overall_interpretation = "VCF does not provide enough strong technical candidates for automated MRD use."

    mrd_signal_class = "not_assessable_without_clinical_context"
    if kg_rows:
        mrd_signal_class = "biologically_supported_candidate_present"
    elif strong_rows:
        mrd_signal_class = "technical_signal_present_biological_context_limited"

    control_class = "matched_control_clean" if summary.get("matched_control_signal_count", 0) == 0 else "matched_control_signal_detected"
    annotation_class = (
        "ncbi_gene_interval_annotated"
        if summary.get("ncbi_exact_gene_interval_count", 0) == summary.get("record_count", -1)
        else "partial_or_missing_gene_annotation"
    )

    top_technical = [make_variant_brief(row) for row in sort_by_numeric(strong_rows, "tumor_af")[:12]]
    top_review = [make_variant_brief(row) for row in sort_by_numeric(review_rows + weak_rows + hold_rows, "tumor_af")[:12]]
    top_kg = [make_variant_brief(row) for row in sort_by_numeric(kg_rows, "tumor_af")[:8]]
    hold_brief = [make_variant_brief(row) for row in sort_by_numeric(hold_rows, "tumor_af")[:8]]

    return {
        "case_qc_class": sample_qc_class,
        "overall_interpretation": overall_interpretation,
        "mrd_signal_class": mrd_signal_class,
        "control_class": control_class,
        "annotation_class": annotation_class,
        "strong_candidate_count": len(strong_rows),
        "review_candidate_count": len(review_rows),
        "weak_candidate_count": len(weak_rows),
        "hold_candidate_count": len(hold_rows),
        "review_or_worse_fraction": round(review_or_worse_fraction, 4),
        "high_af_variant_count": len(high_af_rows),
        "low_support_variant_count": len(low_support_rows),
        "artifact_review_variant_count": len(artifact_rows),
        "kg_supported_variant_count": len(kg_rows),
        "top_technical_candidates": top_technical,
        "top_review_candidates": top_review,
        "kg_supported_candidates": top_kg,
        "hold_candidates": hold_brief,
        "case_limitations": [
            "NCBI interval annotation does not provide coding consequence or protein change.",
            "Cancer type, specimen type, assay/panel design, limit of detection, and clinical timepoint are not provided in the VCF.",
            "MRD positivity cannot be clinically concluded from a single VCF without longitudinal and disease-context metadata.",
        ],
    }


def get_sample_number(sample: dict[str, str], key: str, fallback: Any = None) -> float | None:
    value = sample.get(key)
    if value in (None, ""):
        value = fallback
    return to_float(value)


def compute_variant_qc(
    variant: dict[str, Any],
    meta: dict[str, Any],
    evidence_context: dict[str, Any],
    ncbi_annotation: dict[str, Any] | None = None,
) -> dict[str, Any]:
    info = variant["info"]
    samples = variant["sample_values"]
    tumor_sample = meta["tumor_sample"]
    control_samples = meta["control_samples"]
    tumor = samples.get(tumor_sample, {})
    control_name = control_samples[0] if control_samples else ""
    control = samples.get(control_name, {}) if control_name else {}

    info_dp = to_float(info.get("DP"))
    info_vd = to_float(info.get("VD"))
    info_af = to_float(info.get("AF"))

    tumor_dp = get_sample_number(tumor, "DP", info_dp)
    tumor_vd = get_sample_number(tumor, "VD", info_vd)
    if tumor_vd is None:
        tumor_vd = alt_depth_from_ad(tumor.get("AD"))
    tumor_af = get_sample_number(tumor, "AF", info_af)
    if tumor_af is None:
        af_values = parse_number_list(tumor.get("AF"))
        if af_values:
            tumor_af = max(af_values)
    if tumor_af is None and tumor_dp and tumor_vd is not None and tumor_dp > 0:
        tumor_af = tumor_vd / tumor_dp
    tumor_hiaf = get_sample_number(tumor, "HIAF")
    tumor_mq = get_sample_number(tumor, "MQ")
    if tumor_mq is None:
        tumor_mq = alt_info_number(info, "MMQ")
    tumor_bq = alt_info_number(info, "MBQ")
    tumor_tlod = alt_info_number(info, "TLOD")
    tumor_sn = get_sample_number(tumor, "SN")
    tumor_nm = get_sample_number(tumor, "NM")
    tumor_qual_fmt = get_sample_number(tumor, "QUAL")
    control_dp = get_sample_number(control, "DP")
    control_vd = get_sample_number(control, "VD")
    if control_vd is None:
        control_vd = alt_depth_from_ad(control.get("AD"))
    control_af = get_sample_number(control, "AF")
    if control_af is None:
        control_af_values = parse_number_list(control.get("AF"))
        if control_af_values:
            control_af = max(control_af_values)
    if control_af is None and control_dp and control_vd is not None and control_dp > 0:
        control_af = control_vd / control_dp

    alt_left, alt_right = parse_pair(tumor.get("ALD"))
    if alt_left is None or alt_right is None:
        alt_left, alt_right = alt_strand_from_sb(tumor.get("SB"))
    alt_sum = (alt_left or 0) + (alt_right or 0) if alt_left is not None and alt_right is not None else None
    alt_min_fraction = None
    if alt_sum and alt_sum > 0:
        alt_min_fraction = min(alt_left or 0, alt_right or 0) / alt_sum

    vtype = clean_text(info.get("TYPE")) or ("SNV" if len(variant["ref"]) == 1 and len(variant["alt"]) == 1 else "Complex")
    status = clean_text(info.get("STATUS"))
    msi = to_float(info.get("MSI"))
    msilen = to_float(info.get("MSILEN"))
    shift3 = to_float(info.get("SHIFT3"))
    ssf = to_float(info.get("SSF"))

    flags: list[str] = []
    review_reasons: list[str] = []
    hard_hold = False
    score = 100

    if variant["filter"] != "PASS":
        flags.append("filter_not_pass")
        review_reasons.append("VCF FILTER is not PASS")
        score -= 40
        hard_hold = True
    if status and status != "StrongSomatic":
        flags.append("status_not_strong_somatic")
        review_reasons.append("somatic status is not StrongSomatic")
        score -= 25

    if tumor_dp is not None:
        if tumor_dp < 50:
            flags.append("tumor_depth_low")
            review_reasons.append("tumor depth below 50")
            score -= 25
        elif tumor_dp < 100:
            flags.append("tumor_depth_moderate")
            review_reasons.append("tumor depth between 50 and 99")
            score -= 10
    else:
        flags.append("tumor_depth_missing")
        review_reasons.append("tumor depth missing")
        score -= 15

    if tumor_vd is not None:
        if tumor_vd < 5:
            flags.append("variant_support_very_low")
            review_reasons.append("variant depth below 5")
            score -= 30
        elif tumor_vd < 8:
            flags.append("variant_support_low")
            review_reasons.append("variant depth between 5 and 7")
            score -= 12
    else:
        flags.append("variant_support_missing")
        review_reasons.append("variant depth missing")
        score -= 20

    if tumor_af is not None:
        if tumor_af < 0.02:
            flags.append("vaf_below_2_percent")
            review_reasons.append("tumor AF below 2 percent")
            score -= 30
        elif tumor_af < 0.05:
            flags.append("vaf_low_review")
            review_reasons.append("tumor AF between 2 and 5 percent")
            score -= 10
    else:
        flags.append("vaf_missing")
        review_reasons.append("tumor AF missing")
        score -= 20

    control_has_observed_alt = control_vd is not None and control_vd > 0
    if control_has_observed_alt and control_af is not None and control_af >= 0.02:
        flags.append("matched_control_high_af")
        review_reasons.append("matched/control sample AF is at least 2 percent")
        score -= 40
        hard_hold = True
    elif control_vd is not None and control_vd >= 3:
        flags.append("matched_control_recurrent_alt_reads")
        review_reasons.append("matched/control sample has at least 3 alt reads")
        score -= 30
    elif control_has_observed_alt or (control_vd is None and control_af is not None and control_af > 0):
        flags.append("matched_control_trace_signal")
        review_reasons.append("matched/control sample has non-zero alt signal")
        score -= 15

    if alt_min_fraction is not None and alt_sum is not None and alt_sum >= 5:
        if alt_min_fraction <= 0.10:
            flags.append("severe_alt_strand_imbalance")
            review_reasons.append("alt reads are almost one-strand only")
            score -= 20
        elif alt_min_fraction <= 0.20:
            flags.append("alt_strand_imbalance_review")
            review_reasons.append("alt reads show strand imbalance")
            score -= 8

    repeat_or_indel = False
    if vtype != "SNV":
        flags.append("indel_or_complex_variant")
        review_reasons.append("indel/complex variant needs local alignment review")
        repeat_or_indel = True
        score -= 8
    if (msi is not None and msi >= 3) or (msilen is not None and msilen >= 3) or (shift3 is not None and shift3 >= 3):
        flags.append("repeat_or_shift_context")
        review_reasons.append("repeat/shift context may increase alignment artifact risk")
        repeat_or_indel = True
        score -= 8

    if tumor_mq is not None and tumor_mq < 40:
        flags.append("mapping_quality_low")
        review_reasons.append("mapping quality below 40")
        score -= 15
    if tumor_bq is not None and tumor_bq < 25:
        flags.append("base_quality_low")
        review_reasons.append("Mutect2 alt median base quality below 25")
        score -= 10
    if tumor_tlod is not None and tumor_tlod < 6:
        flags.append("tumor_lod_low")
        review_reasons.append("Mutect2 TLOD below 6")
        score -= 15
    if tumor_sn is not None and tumor_sn < 10:
        flags.append("signal_noise_low")
        review_reasons.append("SN below 10")
        score -= 10
    if tumor_nm is not None and tumor_nm > 3:
        flags.append("mismatch_noise_high")
        review_reasons.append("NM above 3")
        score -= 8
    if ssf is not None and ssf >= 0.05:
        flags.append("ssf_elevated")
        review_reasons.append("SSF at least 0.05")
        score -= 8
    if variant["qual"] is not None and variant["qual"] < 80:
        flags.append("site_qual_low")
        review_reasons.append("site QUAL below 80")
        score -= 8
    if tumor_qual_fmt is not None and tumor_qual_fmt < 20:
        flags.append("sample_qual_low")
        review_reasons.append("sample QUAL below 20")
        score -= 8
    if control_dp is not None and control_dp < 50:
        flags.append("control_depth_low")
        review_reasons.append("matched/control depth below 50")
        score -= 5

    variant_id = f"{variant['chrom']}:{variant['pos']}:{variant['ref']}>{variant['alt']}"
    vcf_genes = extract_genes(info)
    ncbi_record = (ncbi_annotation or {}).get("per_variant", {}).get(variant_id, {})
    ncbi_exact = ncbi_record.get("genes", []) or []
    ncbi_nearby = ncbi_record.get("nearby_genes_5kb", []) or []
    ncbi_genes = sorted({item.get("symbol", "") for item in ncbi_exact if item.get("symbol")})
    nearby_gene_symbols = sorted({item.get("symbol", "") for item in ncbi_nearby if item.get("symbol")})
    genes = sorted(set(vcf_genes) | set(ncbi_genes))
    coordinate_key = None
    if variant["pos"] is not None:
        coordinate_key = (normalize_chrom(variant["chrom"]), int(variant["pos"]))
    coord_matches = evidence_context.get("coordinate_index", {}).get(coordinate_key, []) if coordinate_key else []

    gene_matches: list[dict[str, Any]] = []
    for gene in genes:
        gene_matches.extend(evidence_context.get("gene_index", {}).get(gene, []))

    chip_overlap = sorted(set(genes) & CHIP_REVIEW_GENES)
    if chip_overlap:
        flags.append("common_chip_gene_review")
        review_reasons.append(f"gene overlaps common CHIP/hematopoietic review set: {','.join(chip_overlap)}")
        score -= 3

    if gene_matches:
        kg_match_status = "gene_evidence_match"
    elif coord_matches:
        kg_match_status = "coordinate_evidence_match"
    elif genes:
        kg_match_status = "gene_annotated_no_local_kg_match"
    else:
        kg_match_status = "coordinate_unannotated_no_gene"

    kg_evidence_records = gene_matches[:5] if gene_matches else coord_matches[:5]
    kg_evidence_ids = [item.get("evidence_id", "") for item in kg_evidence_records if item.get("evidence_id")]
    kg_confidences = sorted({clean_text(item.get("confidence_score")) for item in kg_evidence_records if clean_text(item.get("confidence_score"))})
    kg_titles = [clean_text(item.get("title")) for item in kg_evidence_records if clean_text(item.get("title"))]
    if kg_evidence_records:
        kg_provenance_summary = (
            f"Local KG matched {len(kg_evidence_records)} evidence record(s): "
            f"{', '.join(kg_evidence_ids[:4])}. "
            f"Confidence: {', '.join(kg_confidences) if kg_confidences else 'not_reported'}."
        )
    elif genes:
        kg_provenance_summary = (
            "Gene interval annotation is available, but no direct local KG evidence record matched this gene/coordinate."
        )
    else:
        kg_provenance_summary = (
            "No gene annotation is available, so local KG source tracing is blocked at the gene-evidence layer."
        )

    if hard_hold or score < 40:
        tier = "hold_or_reject_without_review"
        action = "do_not_use_for_mrd_without_manual_review"
    elif score < 60:
        tier = "weak_or_artifact_prone_candidate"
        action = "confirm_with_replicate_or_orthogonal_evidence"
    elif score < 80:
        tier = "review_candidate"
        action = "manual_review_before_mrd_use"
    else:
        tier = "strong_technical_candidate"
        action = "retain_as_technical_candidate"

    if repeat_or_indel and tier == "strong_technical_candidate":
        tier = "review_candidate"
        action = "manual_review_before_mrd_use"

    if tier == "hold_or_reject_without_review":
        mrd_entry_role = "do_not_enter_mrd_adjudication_without_manual_rescue"
        mrd_entry_rationale = "technical QC holds this call before MRD adjudication."
    elif kg_match_status in {"gene_evidence_match", "coordinate_evidence_match"} and tier == "strong_technical_candidate":
        mrd_entry_role = "primary_mrd_adjudication_candidate"
        mrd_entry_rationale = "technically strong call with direct local KG support."
    elif kg_match_status in {"gene_evidence_match", "coordinate_evidence_match"}:
        mrd_entry_role = "review_before_mrd_adjudication"
        mrd_entry_rationale = "local KG support exists, but technical or CHIP/artifact guards require expert review."
    elif genes and tier == "strong_technical_candidate":
        mrd_entry_role = "conditional_mrd_adjudication_candidate"
        mrd_entry_rationale = "technically strong annotated call, but direct local KG support is limited."
    elif tier in {"review_candidate", "weak_or_artifact_prone_candidate"}:
        mrd_entry_role = "technical_candidate_needs_biological_context"
        mrd_entry_rationale = "call may be triaged, but needs stronger biological/contextual evidence before MRD use."
    else:
        mrd_entry_role = "do_not_enter_mrd_adjudication_without_manual_rescue"
        mrd_entry_rationale = "call does not meet MRD adjudication entry criteria."

    artifact_flag = any(
        flag in flags
        for flag in [
            "repeat_or_shift_context",
            "alt_strand_imbalance_review",
            "severe_alt_strand_imbalance",
            "indel_or_complex_variant",
            "mapping_quality_low",
            "ssf_elevated",
        ]
    )
    low_support_flag = "variant_support_low" in flags or "variant_support_very_low" in flags
    has_direct_kg_support = kg_match_status in {"gene_evidence_match", "coordinate_evidence_match"}

    caller_scenarios = ["clinical_report_audit_trail"]
    if genes:
        caller_scenarios.extend(["tumor_informed_marker_tracking", "clone_dynamics_or_recurrence_origin"])
    if chip_overlap or control_has_observed_alt:
        caller_scenarios.append("chip_or_hematopoietic_origin_audit")
    if low_support_flag or "vaf_low_review" in flags or "vaf_below_2_percent" in flags:
        caller_scenarios.append("low_frequency_serial_confirmation")
    if sorted(set(genes) & TARGETED_OR_RESISTANCE_REVIEW_GENES):
        caller_scenarios.append("targeted_resistance_or_trial_evidence")
    if artifact_flag or low_support_flag or control_has_observed_alt:
        caller_scenarios.append("fragmentomics_caller_discordance_audit")
    caller_scenarios = sorted(
        set(caller_scenarios),
        key=lambda item: CALLER_SCENARIO_ORDER.index(item) if item in CALLER_SCENARIO_ORDER else len(CALLER_SCENARIO_ORDER),
    )
    required_context_inputs = sorted(
        {
            requirement
            for scenario in caller_scenarios
            for requirement in CALLER_SCENARIO_INFO.get(scenario, {}).get("required_inputs_for_completion", [])
        }
    )
    caller_workflow_summary = (
        "VCF/caller layer completed: variant identity, gene interval annotation when available, VAF/depth, caller QC flags, "
        "KG provenance lookup, and source-trace action. Final biological source assignment still requires: "
        + (", ".join(required_context_inputs) if required_context_inputs else "expert review context")
        + "."
    )
    caller_scenario_completion = [
        {
            "scenario": scenario,
            "vcf_only_contribution": CALLER_SCENARIO_INFO.get(scenario, {}).get("vcf_only_contribution", ""),
            "why_caller_matters": CALLER_SCENARIO_INFO.get(scenario, {}).get("why_caller_matters", ""),
            "required_inputs_for_completion": CALLER_SCENARIO_INFO.get(scenario, {}).get("required_inputs_for_completion", []),
        }
        for scenario in caller_scenarios
    ]

    if not genes:
        evidence_mode = "vcf_only_source_trace_qc"
    elif has_direct_kg_support:
        evidence_mode = "caller_source_trace"
    else:
        evidence_mode = "insufficient_context"

    if tier == "hold_or_reject_without_review":
        source_trace_status = "technical_rescue_required_before_source_assignment"
        recommended_review_action = "manual_caller_rescue_then_source_tracing"
    elif chip_overlap and control_has_observed_alt:
        source_trace_status = "hematopoietic_or_CHIP_origin_suspected"
        recommended_review_action = "matched_wbc_or_buffy_coat_chip_review"
    elif artifact_flag and tier == "weak_or_artifact_prone_candidate":
        source_trace_status = "technical_artifact_origin_suspected"
        recommended_review_action = "local_alignment_and_caller_artifact_review"
    elif not genes and not has_direct_kg_support:
        source_trace_status = "untraceable_missing_annotation_or_context"
        recommended_review_action = "add_gene_functional_annotation_before_source_tracing"
    elif chip_overlap:
        source_trace_status = "traceable_caller_signal_needs_serial_confirmation"
        recommended_review_action = "matched_wbc_review_and_serial_confirmation"
    elif has_direct_kg_support and tier == "strong_technical_candidate" and not control_has_observed_alt:
        source_trace_status = "traceable_caller_evidence_candidate"
        recommended_review_action = "tumor_concordance_timepoint_lod_and_serial_trace_review"
    elif tier in {"strong_technical_candidate", "review_candidate"}:
        source_trace_status = "traceable_caller_signal_needs_serial_confirmation"
        recommended_review_action = "serial_sampling_or_orthogonal_source_confirmation"
    elif artifact_flag or low_support_flag:
        source_trace_status = "cross_sectional_signal_needs_contextualization"
        recommended_review_action = "manual_caller_qc_before_source_assignment"
    else:
        source_trace_status = "untraceable_missing_annotation_or_context"
        recommended_review_action = "add_context_annotation_and_expert_review"

    variant_label = f"{variant_id}"
    if genes:
        variant_label = f"{variant_id} ({','.join(genes[:3])})"
    if source_trace_status == "traceable_caller_evidence_candidate":
        variant_nlp_qc_conclusion = (
            f"{variant_label} is a technically strong caller-derived candidate with direct local KG support; "
            "its source should be traced with tumor concordance, collection timepoint, assay LOD, WBC status, and serial trajectory."
        )
    elif source_trace_status == "traceable_caller_signal_needs_serial_confirmation":
        variant_nlp_qc_conclusion = (
            f"{variant_label} is traceable at caller/gene level, but the source is not resolved; use serial sampling, orthogonal confirmation, "
            "tumor concordance, and WBC review to decide whether it is persistent tumor-derived signal, transient signal, hematopoietic origin, or technical noise."
        )
    elif source_trace_status == "hematopoietic_or_CHIP_origin_suspected":
        variant_nlp_qc_conclusion = (
            f"{variant_label} is best explained first as possible CHIP/hematopoietic-origin signal; matched WBC or buffy-coat review is needed "
            "before assigning it to tumor MRD."
        )
    elif source_trace_status == "technical_artifact_origin_suspected":
        variant_nlp_qc_conclusion = (
            f"{variant_label} is currently best explained as technical/artifact-prone caller signal; prioritize local alignment, strand/context, "
            "depth, and caller-level rescue before biological source assignment."
        )
    elif source_trace_status == "technical_rescue_required_before_source_assignment":
        variant_nlp_qc_conclusion = (
            f"{variant_label} has caller-QC features that require technical rescue before source assignment. This is not a biological rejection; "
            "use local alignment review, replicate/orthogonal confirmation, and KG provenance to decide whether the signal remains traceable."
        )
    elif source_trace_status == "cross_sectional_signal_needs_contextualization":
        variant_nlp_qc_conclusion = (
            f"{variant_label} is a cross-sectional caller signal with unresolved source. It needs technical review plus tumor/WBC/serial context "
            "before assigning the signal to tumor, hematopoietic, transient, or artifact origin."
        )
    else:
        variant_nlp_qc_conclusion = (
            f"{variant_label} cannot be assigned to a biological source with the current context; add gene/function annotation, "
            "tumor concordance, WBC status, and serial sampling information."
        )

    source_trace_reason_codes: list[str] = []
    if tier == "strong_technical_candidate":
        source_trace_reason_codes.append("caller_qc_strong")
    elif tier == "review_candidate":
        source_trace_reason_codes.append("caller_qc_review_required")
    elif tier == "weak_or_artifact_prone_candidate":
        source_trace_reason_codes.append("caller_qc_artifact_prone")
    else:
        source_trace_reason_codes.append("technical_rescue_priority")
    if kg_evidence_records:
        source_trace_reason_codes.append("local_kg_evidence_matched")
    elif genes:
        source_trace_reason_codes.append("gene_annotated_without_local_kg_match")
    else:
        source_trace_reason_codes.append("gene_annotation_missing")
    if chip_overlap:
        source_trace_reason_codes.append("chip_gene_overlap")
    if control_has_observed_alt:
        source_trace_reason_codes.append("matched_control_alt_signal_observed")
    if artifact_flag:
        source_trace_reason_codes.append("technical_artifact_features_present")
    if low_support_flag:
        source_trace_reason_codes.append("low_variant_support")
    source_trace_reason_codes.extend(
        [
            "tumor_tissue_concordance_not_provided",
            "serial_sampling_context_not_provided",
            "assay_lod_context_not_provided",
        ]
    )

    inner_qc_parts = [
        f"technical={tier}",
        f"kg={kg_match_status}",
        f"source={source_trace_status}",
    ]
    if chip_overlap:
        inner_qc_parts.append(f"CHIP_review={','.join(chip_overlap)}")
    if control_has_observed_alt:
        inner_qc_parts.append("matched_control_alt_signal=yes")
    if artifact_flag:
        inner_qc_parts.append("artifact_features=yes")
    if kg_evidence_records:
        inner_qc_parts.append(f"kg_evidence_ids={','.join(kg_evidence_ids[:3])}")
    inner_qc_conclusion = "; ".join(inner_qc_parts) + "."

    return {
        "variant_id": variant_id,
        "chrom": variant["chrom"],
        "pos": variant["pos"],
        "ref": variant["ref"],
        "alt": variant["alt"],
        "type": vtype,
        "filter": variant["filter"],
        "status": status,
        "site_qual": variant["qual"],
        "tumor_sample": tumor_sample,
        "tumor_dp": tumor_dp,
        "tumor_vd": tumor_vd,
        "tumor_af": tumor_af,
        "tumor_hiaf": tumor_hiaf,
        "tumor_mq": tumor_mq,
        "tumor_bq": tumor_bq,
        "tumor_tlod": tumor_tlod,
        "tumor_sn": tumor_sn,
        "tumor_nm": tumor_nm,
        "tumor_alt_forward": alt_left,
        "tumor_alt_reverse": alt_right,
        "tumor_alt_min_strand_fraction": alt_min_fraction,
        "control_sample": control_name,
        "control_dp": control_dp,
        "control_vd": control_vd,
        "control_af": control_af,
        "msi": msi,
        "msilen": msilen,
        "shift3": shift3,
        "ssf": ssf,
        "genes": ";".join(genes) if genes else "",
        "gene_annotation_source": "VCF" if vcf_genes else ("NCBI" if ncbi_genes else ""),
        "ncbi_assembly": ncbi_record.get("assembly", ""),
        "ncbi_annotation_status": ncbi_record.get("annotation_status", ""),
        "ncbi_gene_ids": ";".join(
            item.get("gene_id", "") for item in ncbi_exact if item.get("gene_id")
        ),
        "ncbi_nearby_genes_5kb": ";".join(nearby_gene_symbols),
        "kg_match_status": kg_match_status,
        "kg_gene_match_count": len(gene_matches),
        "kg_coordinate_match_count": len(coord_matches),
        "kg_evidence_ids": ";".join(kg_evidence_ids),
        "kg_provenance_summary": kg_provenance_summary,
        "kg_evidence_details": json.dumps(kg_evidence_records, ensure_ascii=False),
        "technical_score": max(0, min(100, round(score, 1))),
        "qc_tier": tier,
        "qc_action": action,
        "mrd_entry_role": mrd_entry_role,
        "mrd_entry_rationale": mrd_entry_rationale,
        "source_trace_status": source_trace_status,
        "evidence_mode": evidence_mode,
        "caller_dependent_mrd_scenarios": ";".join(caller_scenarios),
        "caller_workflow_summary": caller_workflow_summary,
        "caller_required_context_to_complete": ";".join(required_context_inputs),
        "caller_scenario_completion": json.dumps(caller_scenario_completion, ensure_ascii=False),
        "recommended_review_action": recommended_review_action,
        "variant_nlp_qc_conclusion": variant_nlp_qc_conclusion,
        "source_trace_reason_codes": ";".join(source_trace_reason_codes),
        "inner_qc_conclusion": inner_qc_conclusion,
        "rag_evidence_status": "rag_not_run",
        "rag_agent_summary": "",
        "rag_support_summary": "",
        "rag_supporting_evidence": "",
        "rag_uncertainty_checklist": "",
        "rag_evidence_details": "[]",
        "graphrag_path_summary": "",
        "graphrag_path_nodes": "[]",
        "graphrag_path_edges": "[]",
        "source_trace_synthesis": "",
        "source_trace_explanation": "",
        "supporting_evidence_summary": "",
        "contradictory_evidence_summary": "",
        "missing_context_summary": "",
        "next_confirmation_step": "",
        "source_assignment_confidence": "",
        "qc_flags": ";".join(flags),
        "review_reasons": "; ".join(review_reasons),
    }


def summarize_variants(
    rows: list[dict[str, Any]],
    meta: dict[str, Any],
    evidence_context: dict[str, Any],
    ncbi_annotation: dict[str, Any] | None = None,
) -> dict[str, Any]:
    def nums(key: str) -> list[float]:
        values = []
        for row in rows:
            value = row.get(key)
            if isinstance(value, (int, float)) and not math.isinf(value):
                values.append(float(value))
        return values

    def stats(key: str) -> dict[str, float | None]:
        values = nums(key)
        if not values:
            return {"min": None, "median": None, "max": None}
        return {"min": min(values), "median": median(values), "max": max(values)}

    flag_counts: Counter[str] = Counter()
    for row in rows:
        for flag in clean_text(row.get("qc_flags")).split(";"):
            if flag:
                flag_counts[flag] += 1

    tier_counts = Counter(clean_text(row.get("qc_tier")) for row in rows)
    type_counts = Counter(clean_text(row.get("type")) for row in rows)
    filter_counts = Counter(clean_text(row.get("filter")) for row in rows)
    status_counts = Counter(clean_text(row.get("status")) for row in rows)
    kg_counts = Counter(clean_text(row.get("kg_match_status")) for row in rows)
    annotation_status_counts = Counter(clean_text(row.get("ncbi_annotation_status")) or "not_run" for row in rows)
    gene_annotation_source_counts = Counter(clean_text(row.get("gene_annotation_source")) or "not_available" for row in rows)

    annotated = sum(1 for row in rows if clean_text(row.get("genes")))
    ncbi_exact = sum(1 for row in rows if clean_text(row.get("ncbi_annotation_status")) == "exact_gene_interval")
    ncbi_nearby = sum(1 for row in rows if clean_text(row.get("ncbi_annotation_status")) == "nearby_gene_5kb")
    exact_coord_matches = sum(1 for row in rows if int(row.get("kg_coordinate_match_count") or 0) > 0)
    gene_matches = sum(1 for row in rows if int(row.get("kg_gene_match_count") or 0) > 0)
    control_signal = sum(1 for row in rows if (row.get("control_vd") or 0) > 0)

    if annotated == 0:
        overall_status = "TECHNICAL_QC_COMPLETE_BIOLOGICAL_KG_MATCH_BLOCKED_BY_MISSING_ANNOTATION"
    elif gene_matches == 0 and exact_coord_matches == 0:
        overall_status = "TECHNICAL_QC_COMPLETE_NO_DIRECT_LOCAL_KG_VARIANT_MATCH"
    else:
        overall_status = "TECHNICAL_QC_COMPLETE_WITH_LOCAL_KG_MATCHES"

    result = {
        "overall_status": overall_status,
        "interpretation_scope": "technical_candidate_qc_not_clinical_mrd_interpretation",
        "record_count": len(rows),
        "samples": meta.get("samples", []),
        "tumor_sample": meta.get("tumor_sample", ""),
        "control_samples": meta.get("control_samples", []),
        "type_counts": dict(type_counts),
        "filter_counts": dict(filter_counts),
        "status_counts": dict(status_counts),
        "tier_counts": dict(tier_counts),
        "kg_match_counts": dict(kg_counts),
        "annotated_variant_count": annotated,
        "ncbi_exact_gene_interval_count": ncbi_exact,
        "ncbi_nearby_gene_5kb_count": ncbi_nearby,
        "ncbi_annotation_status_counts": dict(annotation_status_counts),
        "gene_annotation_source_counts": dict(gene_annotation_source_counts),
        "gene_match_variant_count": gene_matches,
        "coordinate_match_variant_count": exact_coord_matches,
        "matched_control_signal_count": control_signal,
        "tumor_dp_stats": stats("tumor_dp"),
        "tumor_vd_stats": stats("tumor_vd"),
        "tumor_af_stats": stats("tumor_af"),
        "control_dp_stats": stats("control_dp"),
        "control_af_stats": stats("control_af"),
        "technical_score_stats": stats("technical_score"),
        "flag_counts": dict(flag_counts.most_common()),
        "kg_context": {
            "evidence_rows": evidence_context.get("evidence_rows", 0),
            "confidence_counts": dict(evidence_context.get("confidence_counts", {})),
            "graph_inclusion_counts": dict(evidence_context.get("graph_inclusion_counts", {})),
            "audit_priority_counts": dict(evidence_context.get("audit_priority_counts", {})),
            "mrd_related_rows": evidence_context.get("mrd_related_rows", 0),
            "chip_related_rows": evidence_context.get("chip_related_rows", 0),
            "false_positive_rows": evidence_context.get("false_positive_rows", 0),
            "matched_wbc_rows": evidence_context.get("matched_wbc_rows", 0),
            "germline_filter_rows": evidence_context.get("germline_filter_rows", 0),
            "chip_filter_rows": evidence_context.get("chip_filter_rows", 0),
            "top_genes": sortable_counter(evidence_context.get("gene_counts", Counter()), 12),
            "top_assays": sortable_counter(evidence_context.get("assay_counts", Counter()), 12),
            "top_sample_types": sortable_counter(evidence_context.get("sample_type_counts", Counter()), 8),
            "top_cancer_types": sortable_counter(evidence_context.get("cancer_type_counts", Counter()), 8),
            "coordinate_index_size": len(evidence_context.get("coordinate_index", {})),
        },
        "ncbi_annotation": {
            key: value
            for key, value in (ncbi_annotation or {"enabled": False}).items()
            if key != "per_variant"
        },
    }
    result["case_assessment"] = build_case_assessment(result, rows)
    result["mrd_clinical_interpretation_entry_gate"] = build_mrd_entry_gate(result, rows)
    result["molecular_evidence_qc"] = build_molecular_evidence_qc(
        result,
        rows,
        result["mrd_clinical_interpretation_entry_gate"],
    )
    return result


def write_variant_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "variant_id",
        "chrom",
        "pos",
        "ref",
        "alt",
        "type",
        "filter",
        "status",
        "site_qual",
        "tumor_sample",
        "tumor_dp",
        "tumor_vd",
        "tumor_af",
        "tumor_hiaf",
        "tumor_mq",
        "tumor_bq",
        "tumor_tlod",
        "tumor_sn",
        "tumor_nm",
        "tumor_alt_forward",
        "tumor_alt_reverse",
        "tumor_alt_min_strand_fraction",
        "control_sample",
        "control_dp",
        "control_vd",
        "control_af",
        "msi",
        "msilen",
        "shift3",
        "ssf",
        "genes",
        "gene_annotation_source",
        "ncbi_assembly",
        "ncbi_annotation_status",
        "ncbi_gene_ids",
        "ncbi_nearby_genes_5kb",
        "kg_match_status",
        "kg_gene_match_count",
        "kg_coordinate_match_count",
        "kg_evidence_ids",
        "kg_provenance_summary",
        "kg_evidence_details",
        "technical_score",
        "qc_tier",
        "qc_action",
        "mrd_entry_role",
        "mrd_entry_rationale",
        "source_trace_status",
        "evidence_mode",
        "caller_dependent_mrd_scenarios",
        "caller_workflow_summary",
        "caller_required_context_to_complete",
        "caller_scenario_completion",
        "recommended_review_action",
        "variant_nlp_qc_conclusion",
        "source_trace_reason_codes",
        "inner_qc_conclusion",
        "rag_evidence_status",
        "rag_agent_summary",
        "rag_support_summary",
        "rag_supporting_evidence",
        "rag_uncertainty_checklist",
        "rag_evidence_details",
        "graphrag_path_summary",
        "graphrag_path_nodes",
        "graphrag_path_edges",
        "source_trace_synthesis",
        "source_trace_explanation",
        "supporting_evidence_summary",
        "contradictory_evidence_summary",
        "missing_context_summary",
        "next_confirmation_step",
        "source_assignment_confidence",
        "qc_flags",
        "review_reasons",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def fmt_stat(stat: dict[str, Any], digits: int = 4) -> str:
    if stat.get("median") is None:
        return "not_available"
    return f"min {stat['min']:.{digits}g}, median {stat['median']:.{digits}g}, max {stat['max']:.{digits}g}"


def counter_line(counter: dict[str, Any]) -> str:
    if not counter:
        return "not_available"
    return ", ".join(f"{key}: {value}" for key, value in counter.items())


def short_text(value: Any, limit: int = 180) -> str:
    text = clean_text(value).replace("|", "/")
    return text if len(text) <= limit else text[: limit - 3] + "..."


def display_genes(row: dict[str, Any]) -> str:
    genes = clean_text(row.get("genes"))
    if genes:
        return genes
    nearby = clean_text(row.get("ncbi_nearby_genes_5kb"))
    if nearby:
        return f"nearby:{nearby}"
    return ""


def write_report(path: Path, summary: dict[str, Any], rows: list[dict[str, Any]], vcf_path: Path, variant_csv: Path, summary_json: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    weakest = sorted(rows, key=lambda row: (row.get("technical_score") if row.get("technical_score") is not None else 999, -float(row.get("tumor_af") or 0)))[:15]
    strongest = sorted(rows, key=lambda row: (-(row.get("technical_score") or 0), -float(row.get("tumor_af") or 0)))[:10]
    kg = summary["kg_context"]
    ncbi = summary.get("ncbi_annotation", {})
    gate = summary["mrd_clinical_interpretation_entry_gate"]
    mqc = summary["molecular_evidence_qc"]

    lines: list[str] = []
    lines.append("# MRD KG-RAG VCF QC Report")
    lines.append("")
    lines.append(f"- Case: `{summary['tumor_sample']}`")
    lines.append(f"- Input VCF: `{vcf_path}`")
    lines.append(f"- Per-variant CSV: `{variant_csv}`")
    lines.append(f"- Summary JSON: `{summary_json}`")
    lines.append(f"- Overall status: `{summary['overall_status']}`")
    lines.append(f"- Scope: `{summary['interpretation_scope']}`")
    lines.append("")
    lines.append("## Executive Summary")
    lines.append("")
    lines.append(f"- VCF records parsed: {summary['record_count']}")
    lines.append(f"- Samples: {', '.join(summary['samples'])}")
    lines.append(f"- Tumor sample inferred from VCF: `{summary['tumor_sample']}`")
    lines.append(f"- Matched/control samples: {', '.join(summary['control_samples']) if summary['control_samples'] else 'not_available'}")
    lines.append(f"- Variant types: {counter_line(summary['type_counts'])}")
    lines.append(f"- FILTER status: {counter_line(summary['filter_counts'])}")
    lines.append(f"- Somatic status: {counter_line(summary['status_counts'])}")
    lines.append(f"- QC tiers: {counter_line(summary['tier_counts'])}")
    lines.append(f"- KG match status: {counter_line(summary['kg_match_counts'])}")
    lines.append(f"- Source trace status: `{mqc['source_trace_status']}`")
    lines.append(f"- Evidence mode: `{mqc['evidence_mode']}`")
    lines.append(f"- Source-tracing readiness: `{gate['gate_class']}` ({gate['entry_permission']})")
    if ncbi.get("enabled"):
        lines.append(f"- NCBI assembly selected: `{ncbi.get('selected_assembly')}` from request `{ncbi.get('requested_assembly')}`")
        lines.append(f"- NCBI exact gene-interval annotations: {summary['ncbi_exact_gene_interval_count']}")
    lines.append("")
    lines.append("## Caller Signal Source-Tracing Explanation")
    lines.append("")
    lines.append(f"- Question: {mqc['qc_question']}")
    lines.append(f"- NLP conclusion: {mqc['nlp_qc_conclusion']}")
    lines.append(f"- Source trace status: `{mqc['source_trace_status']}`")
    lines.append(f"- Evidence mode: `{mqc['evidence_mode']}`")
    lines.append(f"- Recommended review actions: {', '.join(mqc['recommended_review_actions'])}")
    lines.append(f"- Source trace status counts: {counter_line(mqc['source_trace_status_counts'])}")
    lines.append(f"- KG-matched variants surfaced for provenance review: {mqc.get('kg_matched_variant_count', 0)}")
    lines.append(f"- CHIP/WBC-origin review variants surfaced: {mqc.get('chip_review_variant_count', 0)}")
    lines.append("")
    lines.append("### Top Variant-Level Source Explanations")
    lines.append("")
    lines.append("| Variant | Gene | AF | Score | KG evidence | Source trace status | Inner QC | Source explanation |")
    lines.append("|---|---|---:|---:|---|---|---|---|")
    for item in mqc["top_variant_conclusions"][:8]:
        lines.append(
            f"| `{item.get('variant_id')}` | {item.get('genes') or ''} | {fmt_af(item.get('tumor_af'))} | "
            f"{item.get('technical_score') or ''} | {short_text(item.get('kg_evidence_ids') or item.get('kg_match_status'), 90)} | "
            f"`{item.get('source_trace_status') or ''}` | {short_text(item.get('inner_qc_conclusion'), 180)} | "
            f"{short_text(item.get('variant_nlp_qc_conclusion'), 220)} |"
        )
    if summary.get("visual_html"):
        lines.append("")
        lines.append(f"- Visual evidence lanes: `{summary['visual_html']}`")
    rag = summary.get("source_trace_rag", {})
    if rag.get("enabled"):
        lines.append(f"- Multi-agent RAG source-trace packet: `{rag.get('md')}`")
        lines.append(f"- Multi-agent RAG JSON: `{rag.get('json')}`")
        lines.append(f"- RAG mode: `{rag.get('rag_mode', 'case')}`")
        lines.append(f"- RAG integrated variant rows: {rag.get('integrated_variant_rows', 0)}")
        lines.append(f"- RAG direct gene-matched variant rows: {rag.get('direct_gene_matched_variant_rows', 0)}")
        lines.append(f"- Variant-specific RAG rows: {rag.get('variant_specific_rag_rows', 0)}")
        lines.append(f"- GraphRAG path rows: {rag.get('graphrag_path_rows', 0)}")
    lines.append("")
    lines.append("### Caller-Dependent MRD Workflow Coverage")
    lines.append("")
    lines.append("| Workflow | Status | Current VCF/caller layer | To complete final source assignment | Variants | RAG hits |")
    lines.append("|---|---|---|---|---:|---:|")
    for item in mqc.get("caller_workflow_coverage", []):
        requires = item.get("required_inputs_for_completion") or []
        requires_text = ", ".join(requires) if isinstance(requires, list) else clean_text(requires)
        lines.append(
            f"| `{item.get('scenario')}` | `{item.get('status')}` | {short_text(item.get('vcf_only_contribution'), 180)} | "
            f"{short_text(requires_text or 'expert review context', 180)} | {item.get('variant_count', 0)} | {item.get('rag_hit_count', 0)} |"
        )
    lines.append("")
    if mqc.get("rag_integration_summary"):
        rag_summary = mqc["rag_integration_summary"]
        lines.append("### RAG Integration Summary")
        lines.append("")
        lines.append(f"- Variant rows with mapped RAG context: {rag_summary.get('integrated_variant_rows', 0)}")
        lines.append(f"- Variant rows with direct gene-mentioned RAG hits: {rag_summary.get('direct_gene_matched_variant_rows', 0)}")
        lines.append(f"- RAG category hits: {counter_line(rag_summary.get('category_hit_counts', {}))}")
        lines.append(f"- RAG scenario hits: {counter_line(rag_summary.get('scenario_hit_counts', {}))}")
        lines.append(f"- Variant synthesis rows: {rag_summary.get('variant_synthesis_rows', 0)}")
        lines.append(f"- GraphRAG path rows: {rag_summary.get('graphrag_path_rows', 0)}")
        if rag_summary.get("missing_context_items"):
            lines.append("- Missing context retained as unresolved requirements:")
            for item in rag_summary["missing_context_items"][:6]:
                lines.append(f"  - {item}")
        lines.append("")
    if mqc.get("case_source_trace_synthesis"):
        syn = mqc["case_source_trace_synthesis"]
        lines.append("### Final Source-Tracing Synthesis")
        lines.append("")
        lines.append(f"- Dominant explanation: {syn.get('dominant_source_explanation')}")
        lines.append(f"- Workflow readiness: `{syn.get('workflow_readiness')}`")
        eq = syn.get("evidence_quality", {})
        if isinstance(eq, dict):
            lines.append(f"- Evidence quality: case rounds={eq.get('case_round_count')}, variant RAG queries={eq.get('variant_specific_rag_queries')}, synthesized rows={eq.get('variant_rows_with_synthesis')}")
        lines.append("- Requires external metadata: " + ", ".join(syn.get("requires_external_metadata", [])[:8]))
        lines.append("")
    lines.append("## Source-Tracing Readiness")
    lines.append("")
    lines.append(f"- Question: {gate['entry_question']}")
    lines.append(f"- Readiness status: `{gate['gate_class']}`")
    lines.append(f"- Permission: `{gate['entry_permission']}`")
    lines.append(f"- Allowed scope: {gate['allowed_scope']}")
    lines.append(f"- Summary: {gate['verdict']}")
    lines.append(f"- Primary candidates: {gate['primary_candidate_count']}")
    lines.append(f"- Conditional candidates: {gate['conditional_candidate_count']}")
    lines.append(f"- Review candidates: {gate['review_candidate_count']}")
    lines.append(f"- Technical rescue variants: {gate['blocked_variant_count']}")
    if gate["reasons"]:
        lines.append("- Reasons: " + " ".join(gate["reasons"]))
    if gate["guardrails"]:
        lines.append("- Guardrails: " + " ".join(gate["guardrails"]))
    if gate["blockers"]:
        lines.append("- Blockers: " + " ".join(gate["blockers"]))
    lines.append("- Missing for final sign-out: " + ", ".join(gate["missing_for_final_signout"]))
    lines.append("")
    lines.append("## Technical Signal Profile")
    lines.append("")
    lines.append(f"- Tumor depth: {fmt_stat(summary['tumor_dp_stats'])}")
    lines.append(f"- Tumor variant depth: {fmt_stat(summary['tumor_vd_stats'])}")
    lines.append(f"- Tumor AF: {fmt_stat(summary['tumor_af_stats'])}")
    lines.append(f"- Control depth: {fmt_stat(summary['control_dp_stats'])}")
    lines.append(f"- Control AF: {fmt_stat(summary['control_af_stats'])}")
    lines.append(f"- Technical score: {fmt_stat(summary['technical_score_stats'])}")
    lines.append(f"- Matched/control non-zero alt signal records: {summary['matched_control_signal_count']}")
    lines.append("")
    lines.append("## KG-RAG Readiness")
    lines.append("")
    if ncbi.get("enabled"):
        lines.append(f"- NCBI annotation source: {ncbi.get('source')} ({ncbi.get('source_url')})")
        lines.append(f"- NCBI annotation status: {counter_line(summary['ncbi_annotation_status_counts'])}")
        lines.append(f"- Gene annotation source: {counter_line(summary['gene_annotation_source_counts'])}")
    if summary["annotated_variant_count"] == 0:
        lines.append("- No gene or functional annotation fields were detected in this VCF.")
        lines.append("- Direct gene-level biological plausibility checking is blocked unless NCBI or another annotation source supplies a gene interval.")
        lines.append("- The current report uses local KG evidence for assay-level MRD QC principles and records unannotated variants as `coordinate_unannotated_no_gene` unless an exact coordinate match exists.")
    else:
        lines.append(f"- Variants with gene annotation: {summary['annotated_variant_count']}")
        lines.append(f"- Variants with local gene evidence match: {summary['gene_match_variant_count']}")
        lines.append(f"- Variants with exact local coordinate evidence match: {summary['coordinate_match_variant_count']}")
    lines.append(f"- Local evidence rows loaded: {kg['evidence_rows']}")
    lines.append(f"- MRD/ctDNA-related evidence rows: {kg['mrd_related_rows']}")
    lines.append(f"- CHIP-related evidence rows: {kg['chip_related_rows']}")
    lines.append(f"- Artifact/confounder-related evidence rows: {kg['false_positive_rows']}")
    lines.append(f"- Evidence rows mentioning matched WBC use: {kg['matched_wbc_rows']}")
    lines.append(f"- Local coordinate evidence index size: {kg['coordinate_index_size']}")
    lines.append("")
    lines.append("## Main QC Flags")
    lines.append("")
    if summary["flag_counts"]:
        for flag, count in list(summary["flag_counts"].items())[:20]:
            lines.append(f"- {flag}: {count}")
    else:
        lines.append("- No QC flags.")
    lines.append("")
    lines.append("## Lowest-Scoring Variants For Review")
    lines.append("")
    lines.append("| Variant | Gene | Type | DP | VD | AF | Control AF | Score | Tier | Main flags |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|---|---|")
    for row in weakest:
        flags = clean_text(row.get("qc_flags")).replace(";", ", ")
        lines.append(
            f"| `{row['variant_id']}` | {display_genes(row)} | {row['type']} | {row.get('tumor_dp') or ''} | {row.get('tumor_vd') or ''} | "
            f"{row.get('tumor_af') or ''} | {row.get('control_af') or ''} | {row.get('technical_score') or ''} | "
            f"{row.get('qc_tier') or ''} | {flags} |"
        )
    lines.append("")
    lines.append("## Highest-Confidence Technical Candidates")
    lines.append("")
    lines.append("| Variant | Gene | Type | DP | VD | AF | Score | Tier | Entry role |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---|---|")
    for row in strongest:
        lines.append(
            f"| `{row['variant_id']}` | {display_genes(row)} | {row['type']} | {row.get('tumor_dp') or ''} | {row.get('tumor_vd') or ''} | "
            f"{row.get('tumor_af') or ''} | {row.get('technical_score') or ''} | {row.get('qc_tier') or ''} | "
            f"{row.get('mrd_entry_role') or ''} |"
        )
    lines.append("")
    lines.append("## Recommended Next Step")
    lines.append("")
    lines.append("- NCBI interval annotation gives gene-level placement, but it does not predict coding consequence or protein change.")
    lines.append("- For final variant interpretation, add VEP, ANNOVAR, SnpEff, or assay panel annotation when available, then rerun this script.")
    lines.append("- With gene annotation available, KG-RAG can check cancer-gene plausibility, CHIP/hematopoietic origin risk, recurrence/MRD evidence support, and evidence traceability per variant.")
    lines.append("- This report is not a clinical interpretation and should be treated as an automated technical pre-QC plus KG-readiness assessment.")
    lines.append("")

    path.write_text("\n".join(lines), encoding="utf-8")


def fmt_percent(value: Any, digits: int = 1) -> str:
    num = to_float(value)
    if num is None or math.isinf(num):
        return "not_available"
    return f"{num * 100:.{digits}f}%"


def fmt_af(value: Any) -> str:
    num = to_float(value)
    if num is None or math.isinf(num):
        return ""
    return f"{num:.4f}"


def write_case_report(
    path: Path,
    summary: dict[str, Any],
    rows: list[dict[str, Any]],
    vcf_path: Path,
    variant_csv: Path,
    summary_json: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    case = summary["case_assessment"]
    ncbi = summary.get("ncbi_annotation", {})
    kg = summary["kg_context"]
    gate = summary["mrd_clinical_interpretation_entry_gate"]
    mqc = summary["molecular_evidence_qc"]
    generated_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def variant_table(items: list[dict[str, Any]], include_kg: bool = False) -> list[str]:
        if not items:
            return ["No records."]
        header = "| Variant | Gene | Type | AF | VD/DP | Score | Tier | Entry role | Flags |"
        align = "|---|---|---:|---:|---:|---:|---|---|---|"
        if include_kg:
            header = "| Variant | Gene | Type | AF | VD/DP | Score | Entry role | KG evidence |"
            align = "|---|---|---:|---:|---:|---:|---|---|"
        table = [header, align]
        for row in items:
            vd_dp = f"{row.get('tumor_vd')}/{row.get('tumor_dp')}"
            flags = clean_text(row.get("qc_flags")).replace(";", ", ")
            if include_kg:
                table.append(
                    f"| `{row.get('variant_id')}` | {row.get('genes') or ''} | {row.get('type') or ''} | "
                    f"{fmt_af(row.get('tumor_af'))} | {vd_dp} | {row.get('technical_score') or ''} | "
                    f"{row.get('mrd_entry_role') or ''} | {row.get('kg_evidence_ids') or ''} |"
                )
            else:
                table.append(
                    f"| `{row.get('variant_id')}` | {row.get('genes') or ''} | {row.get('type') or ''} | "
                    f"{fmt_af(row.get('tumor_af'))} | {vd_dp} | {row.get('technical_score') or ''} | "
                    f"{row.get('qc_tier') or ''} | {row.get('mrd_entry_role') or ''} | {flags} |"
                )
        return table

    lines: list[str] = []
    lines.append("# MRD VCF Case-Level QC Note")
    lines.append("")
    lines.append("## Patient / Sample")
    lines.append("")
    lines.append(f"- Case ID: `{summary['tumor_sample']}`")
    lines.append("- Patient identifier: not_provided")
    lines.append("- Cancer type: not_provided")
    lines.append("- Specimen type: not_provided")
    lines.append("- Assay / panel: not_provided")
    lines.append(f"- Tumor sample: `{summary['tumor_sample']}`")
    lines.append(f"- Matched/control sample: `{', '.join(summary['control_samples']) if summary['control_samples'] else 'not_available'}`")
    lines.append(f"- Input VCF: `{vcf_path}`")
    lines.append(f"- Generated at: {generated_at}")
    lines.append("")
    lines.append("## Case-Level Impression")
    lines.append("")
    lines.append(f"- Overall QC class: `{case['case_qc_class']}`")
    lines.append(f"- MRD signal class: `{case['mrd_signal_class']}`")
    lines.append(f"- Matched/control assessment: `{case['control_class']}`")
    lines.append(f"- Annotation assessment: `{case['annotation_class']}`")
    lines.append(f"- Source trace status: `{mqc['source_trace_status']}`")
    lines.append(f"- Evidence mode: `{mqc['evidence_mode']}`")
    lines.append(f"- Source-tracing readiness: `{gate['gate_class']}` ({gate['entry_permission']})")
    lines.append(f"- Impression: {mqc['nlp_qc_conclusion']}")
    lines.append("")
    lines.append("## Caller Signal Source-Tracing Explanation")
    lines.append("")
    lines.append(f"- Source trace status: `{mqc['source_trace_status']}`")
    lines.append(f"- Evidence mode: `{mqc['evidence_mode']}`")
    lines.append(f"- NLP QC conclusion: {mqc['nlp_qc_conclusion']}")
    lines.append(f"- Recommended review actions: {', '.join(mqc['recommended_review_actions'])}")
    lines.append(f"- Source trace status counts: {counter_line(mqc['source_trace_status_counts'])}")
    if summary.get("visual_html"):
        lines.append(f"- Visual evidence lanes: `{summary['visual_html']}`")
    rag = summary.get("source_trace_rag", {})
    if rag.get("enabled"):
        lines.append(f"- Multi-agent RAG source-trace packet: `{rag.get('md')}`")
        lines.append(f"- Multi-agent RAG JSON: `{rag.get('json')}`")
        lines.append(f"- RAG mode: `{rag.get('rag_mode', 'case')}`")
        lines.append(f"- RAG integrated variant rows: {rag.get('integrated_variant_rows', 0)}")
        lines.append(f"- RAG direct gene-matched variant rows: {rag.get('direct_gene_matched_variant_rows', 0)}")
        lines.append(f"- Variant-specific RAG rows: {rag.get('variant_specific_rag_rows', 0)}")
        lines.append(f"- GraphRAG path rows: {rag.get('graphrag_path_rows', 0)}")
    lines.append("")
    lines.append("### Caller-Dependent MRD Workflow Coverage")
    lines.append("")
    lines.append("| Workflow | Status | Current VCF/caller layer | To complete final source assignment | Variants | RAG hits |")
    lines.append("|---|---|---|---|---:|---:|")
    for item in mqc.get("caller_workflow_coverage", []):
        requires = item.get("required_inputs_for_completion") or []
        requires_text = ", ".join(requires) if isinstance(requires, list) else clean_text(requires)
        lines.append(
            f"| `{item.get('scenario')}` | `{item.get('status')}` | {short_text(item.get('vcf_only_contribution'), 180)} | "
            f"{short_text(requires_text or 'expert review context', 180)} | {item.get('variant_count', 0)} | {item.get('rag_hit_count', 0)} |"
        )
    lines.append("")
    if mqc.get("rag_integration_summary"):
        rag_summary = mqc["rag_integration_summary"]
        lines.append("### RAG Integration Summary")
        lines.append("")
        lines.append(f"- Variant rows with mapped RAG context: {rag_summary.get('integrated_variant_rows', 0)}")
        lines.append(f"- Variant rows with direct gene-mentioned RAG hits: {rag_summary.get('direct_gene_matched_variant_rows', 0)}")
        lines.append(f"- RAG category hits: {counter_line(rag_summary.get('category_hit_counts', {}))}")
        lines.append(f"- RAG scenario hits: {counter_line(rag_summary.get('scenario_hit_counts', {}))}")
        lines.append(f"- Variant synthesis rows: {rag_summary.get('variant_synthesis_rows', 0)}")
        lines.append(f"- GraphRAG path rows: {rag_summary.get('graphrag_path_rows', 0)}")
        if rag_summary.get("missing_context_items"):
            lines.append("- Missing context retained as unresolved requirements:")
            for item in rag_summary["missing_context_items"][:6]:
                lines.append(f"  - {item}")
        lines.append("")
    if mqc.get("case_source_trace_synthesis"):
        syn = mqc["case_source_trace_synthesis"]
        lines.append("### Final Source-Tracing Synthesis")
        lines.append("")
        lines.append(f"- Dominant explanation: {syn.get('dominant_source_explanation')}")
        lines.append(f"- Workflow readiness: `{syn.get('workflow_readiness')}`")
        eq = syn.get("evidence_quality", {})
        if isinstance(eq, dict):
            lines.append(f"- Evidence quality: case rounds={eq.get('case_round_count')}, variant RAG queries={eq.get('variant_specific_rag_queries')}, synthesized rows={eq.get('variant_rows_with_synthesis')}")
        lines.append("- Requires external metadata: " + ", ".join(syn.get("requires_external_metadata", [])[:8]))
        lines.append("")
    lines.append("### Variant-Level Source Explanations")
    lines.append("")
    lines.append("| Variant | Gene | AF | Score | KG evidence | Source trace status | RAG status | Inner QC | Source explanation |")
    lines.append("|---|---|---:|---:|---|---|---|---|---|")
    for item in mqc["top_variant_conclusions"][:10]:
        lines.append(
            f"| `{item.get('variant_id')}` | {item.get('genes') or ''} | {fmt_af(item.get('tumor_af'))} | "
            f"{item.get('technical_score') or ''} | {short_text(item.get('kg_evidence_ids') or item.get('kg_match_status'), 90)} | "
            f"`{item.get('source_trace_status') or ''}` | `{item.get('rag_evidence_status') or ''}` | {short_text(item.get('inner_qc_conclusion'), 180)} | "
            f"{short_text(item.get('variant_nlp_qc_conclusion'), 220)} |"
        )
    lines.append("")
    lines.append("## Source-Tracing Readiness")
    lines.append("")
    lines.append(f"- Readiness status: `{gate['gate_class']}`")
    lines.append(f"- Permission: `{gate['entry_permission']}`")
    lines.append(f"- Allowed scope: {gate['allowed_scope']}")
    lines.append(f"- Summary: {gate['verdict']}")
    lines.append(f"- Primary candidates: {gate['primary_candidate_count']}")
    lines.append(f"- Conditional candidates: {gate['conditional_candidate_count']}")
    lines.append(f"- Review candidates: {gate['review_candidate_count']}")
    lines.append(f"- Technical rescue variants: {gate['blocked_variant_count']}")
    if gate["reasons"]:
        lines.append("- Reasons: " + " ".join(gate["reasons"]))
    if gate["guardrails"]:
        lines.append("- Guardrails: " + " ".join(gate["guardrails"]))
    if gate["blockers"]:
        lines.append("- Blockers: " + " ".join(gate["blockers"]))
    lines.append("- Missing for final sign-out: " + ", ".join(gate["missing_for_final_signout"]))
    lines.append("")
    lines.append("## Test Overview")
    lines.append("")
    lines.append(f"- Total variants: {summary['record_count']}")
    lines.append(f"- Variant types: {counter_line(summary['type_counts'])}")
    lines.append(f"- PASS records: {summary['filter_counts'].get('PASS', 0)}")
    lines.append(f"- StrongSomatic records: {summary['status_counts'].get('StrongSomatic', 0)}")
    lines.append(f"- Strong technical candidates: {case['strong_candidate_count']}")
    lines.append(f"- Review candidates: {case['review_candidate_count']}")
    lines.append(f"- Weak/artifact-prone candidates: {case['weak_candidate_count']}")
    lines.append(f"- Technical rescue priority records: {case['hold_candidate_count']}")
    lines.append(f"- Review-or-worse fraction: {fmt_percent(case['review_or_worse_fraction'])}")
    lines.append("")
    lines.append("## Technical Adequacy")
    lines.append("")
    lines.append(f"- Tumor depth: {fmt_stat(summary['tumor_dp_stats'])}")
    lines.append(f"- Tumor variant depth: {fmt_stat(summary['tumor_vd_stats'])}")
    lines.append(f"- Tumor AF: {fmt_stat(summary['tumor_af_stats'])}")
    lines.append(f"- Control AF: {fmt_stat(summary['control_af_stats'])}")
    lines.append(f"- Matched/control non-zero alt signal records: {summary['matched_control_signal_count']}")
    lines.append(f"- Low support records: {case['low_support_variant_count']}")
    lines.append(f"- Artifact-review records: {case['artifact_review_variant_count']}")
    lines.append("")
    lines.append("## MRD / KG Evidence Context")
    lines.append("")
    if ncbi.get("enabled"):
        lines.append(f"- NCBI assembly selected: `{ncbi.get('selected_assembly')}`")
        lines.append(f"- NCBI exact gene-interval annotations: {summary['ncbi_exact_gene_interval_count']}")
    lines.append(f"- Local MRD evidence rows loaded: {kg['evidence_rows']}")
    lines.append(f"- MRD/ctDNA-related evidence rows: {kg['mrd_related_rows']}")
    lines.append(f"- CHIP-related evidence rows: {kg['chip_related_rows']}")
    lines.append(f"- Variants with local gene evidence support: {case['kg_supported_variant_count']}")
    lines.append(f"- KG match status: {counter_line(summary['kg_match_counts'])}")
    lines.append("")
    lines.append("## KG-Supported Candidate")
    lines.append("")
    lines.extend(variant_table(case["kg_supported_candidates"], include_kg=True))
    lines.append("")
    lines.append("## Leading Technical Candidates")
    lines.append("")
    lines.extend(variant_table(case["top_technical_candidates"]))
    lines.append("")
    lines.append("## Highest-Priority Review Items")
    lines.append("")
    lines.extend(variant_table(case["top_review_candidates"]))
    lines.append("")
    lines.append("## Technical Rescue Pending Review")
    lines.append("")
    lines.extend(variant_table(case["hold_candidates"]))
    lines.append("")
    lines.append("## Assessment")
    lines.append("")
    if summary["matched_control_signal_count"] > 0:
        lines.append("- The full Mutect/paired VCF contains matched/control alternate-read signal in many records; case-level use requires manual review and should prioritize PASS records with zero observed normal alt depth.")
    else:
        lines.append("- The matched/control sample is clean by observed alternate-read criteria in this parsed VCF.")
    if case["strong_candidate_count"] > 0:
        lines.append("- PASS technical candidates are present; these should be interpreted as technical candidates unless supported by clinical context, assay design, and orthogonal annotation.")
    if case["kg_supported_candidates"]:
        genes = sorted({clean_text(item.get("genes")) for item in case["kg_supported_candidates"] if clean_text(item.get("genes"))})
        lines.append(f"- Current local MRD KG support is strongest for: {', '.join(genes[:8])}.")
    else:
        lines.append("- No variant in the top candidates has direct local MRD KG support.")
    lines.append("- Non-PASS, normal-artifact, germline, multiallelic, strand-bias, slippage, low-depth, or low-TLOD calls should not be used for MRD marker selection without manual review.")
    lines.append("")
    lines.append("## Limitations")
    lines.append("")
    for item in case["case_limitations"]:
        lines.append(f"- {item}")
    lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append(f"- Per-variant QC CSV: `{variant_csv}`")
    lines.append(f"- Summary JSON: `{summary_json}`")
    lines.append("")
    lines.append("This note is an automated research QC summary, not a clinical diagnosis or treatment recommendation.")
    lines.append("")

    path.write_text("\n".join(lines), encoding="utf-8")


def lane_class(value: str) -> str:
    text = clean_text(value).lower()
    if any(token in text for token in ["supported", "strong", "allow", "clean"]):
        return "good"
    if any(token in text for token in ["rescue", "contextualization", "review", "conditional", "missing", "limited"]):
        return "warn"
    if any(token in text for token in ["chip", "artifact", "reject", "block"]):
        return "bad"
    if any(token in text for token in ["repeat"]):
        return "warn"
    return "neutral"


def write_visual_lanes_html(path: Path, summary: dict[str, Any], rows: list[dict[str, Any]], vcf_path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    mqc = summary["molecular_evidence_qc"]
    gate = summary["mrd_clinical_interpretation_entry_gate"]

    def num(value: Any) -> float:
        parsed = to_float(value)
        return parsed if parsed is not None and math.isfinite(parsed) else 0.0

    top_ids = {item.get("variant_id") for item in mqc["top_variant_conclusions"]}
    top_order = {clean_text(item.get("variant_id")): index for index, item in enumerate(mqc["top_variant_conclusions"])}
    display_rows = sorted(
        [row for row in rows if row.get("variant_id") in top_ids],
        key=lambda row: top_order.get(clean_text(row.get("variant_id")), 9999),
    )
    if not display_rows:
        display_rows = sorted(rows, key=lambda row: (-num(row.get("technical_score")), -num(row.get("tumor_af"))))[:12]

    def esc(value: Any) -> str:
        return html.escape(clean_text(value), quote=True)

    def pct(count: int | float, denominator: int | float | None = None) -> str:
        base = denominator if denominator else len(rows)
        if not base:
            return "0.0%"
        return f"{(float(count) / float(base)) * 100:.1f}%"

    aliases = {
        "hold_or_reject_without_review": "Hold / rescue",
        "weak_or_artifact_prone_candidate": "Weak / artifact-prone",
        "review_candidate": "Review candidate",
        "strong_technical_candidate": "Strong candidate",
        "gene_annotated_no_local_kg_match": "No local KG",
        "gene_evidence_match": "KG gene match",
        "coordinate_unannotated_no_gene": "No gene annotation",
        "technical_rescue_required_before_source_assignment": "Technical rescue",
        "cross_sectional_signal_needs_contextualization": "Context needed",
        "traceable_caller_signal_needs_serial_confirmation": "Serial confirmation",
        "technical_artifact_origin_suspected": "Artifact suspected",
        "untraceable_missing_annotation_or_context": "Missing annotation",
        "rag_scenario_context_only": "Scenario RAG",
        "rag_variant_or_direct_context_integrated": "Direct RAG",
        "requires_technical_rescue": "Requires rescue",
        "context_limited_with_scenario_rag": "Context-limited",
        "provisional_kg_supported_context_limited": "KG-supported",
        "provisional_rag_supported_context_limited": "RAG-supported",
        "manual_caller_rescue_then_source_tracing": "Caller rescue",
        "manual_caller_qc_before_source_assignment": "Manual caller QC",
        "serial_sampling_or_orthogonal_source_confirmation": "Serial/orthogonal",
        "local_alignment_and_caller_artifact_review": "Alignment review",
        "add_context_annotation_and_expert_review": "Add context",
        "matched_wbc_review_and_serial_confirmation": "WBC + serial",
        "vcf_layer_plus_rag_context_available": "VCF + mapped RAG context",
        "vcf_layer_available_rag_not_integrated": "VCF only; no mapped RAG",
        "rag_context_available_no_current_variant_trigger": "RAG context only",
        "not_triggered_by_current_vcf": "Not triggered",
        "tumor_informed_marker_tracking": "Tumor-informed",
        "clone_dynamics_or_recurrence_origin": "Clone dynamics",
        "chip_or_hematopoietic_origin_audit": "CHIP/WBC audit",
        "low_frequency_serial_confirmation": "Low-frequency",
        "targeted_resistance_or_trial_evidence": "Target/resistance",
        "fragmentomics_caller_discordance_audit": "Fragmentomics discordance",
        "clinical_report_audit_trail": "Report audit",
        "primary_tumor_marker_set": "tumor markers",
        "tumor_concordance_status": "tumor concordance",
        "serial_ctdna_samples": "serial ctDNA",
        "baseline_tumor_clone_map": "clone map",
        "recurrence_or_metastasis_sample_if_available": "recurrence tissue",
        "matched_wbc_or_buffy_coat": "WBC/buffy",
        "germline_or_hematopoietic_filtering_status": "germline/CH filter",
        "repeat_sample": "repeat sample",
        "assay_lod": "assay LOD",
        "functional_variant_annotation": "functional annotation",
        "therapy_context": "therapy context",
        "actionability_database_or_guideline": "actionability DB",
        "fragmentomics_status_or_score": "fragmentomics score",
        "fragmentomics_assay_metadata": "fragmentomics metadata",
        "expert_signout_policy": "signout policy",
        "clinical_context_for_final_report": "clinical context",
    }

    def label(value: Any) -> str:
        text = clean_text(value)
        if not text:
            return "not available"
        return aliases.get(text, text.replace("_", " "))

    def action_label(row: dict[str, Any]) -> str:
        action = clean_text(row.get("recommended_review_action"))
        if action:
            return label(action)
        step = clean_text(row.get("next_confirmation_step")).lower()
        if "wbc" in step or "buffy" in step:
            return "WBC/buffy review"
        if "serial" in step or "orthogonal" in step:
            return "Serial/orthogonal"
        if "alignment" in step or "caller rescue" in step:
            return "Caller rescue"
        return "Expert review"

    def pill(label: str, value: Any) -> str:
        text = clean_text(value) or "not_available"
        return (
            f'<span class="pill {lane_class(text)}">'
            f'<span class="pill-label">{esc(label)}</span>'
            f'<span class="pill-value">{esc(label_alias(text))}</span>'
            f'</span>'
        )

    def label_alias(value: Any) -> str:
        return label(value)

    def evidence_items(row: dict[str, Any]) -> str:
        raw = clean_text(row.get("kg_evidence_details"))
        if not raw:
            return "<li>No local KG evidence record attached to this variant.</li>"
        try:
            records = json.loads(raw)
        except json.JSONDecodeError:
            records = []
        if not isinstance(records, list) or not records:
            return "<li>No local KG evidence record attached to this variant.</li>"
        items: list[str] = []
        for rec in records[:4]:
            if not isinstance(rec, dict):
                continue
            title = esc(short_text(rec.get("title"), 170))
            evidence_text = esc(short_text(rec.get("evidence_text"), 220))
            ids = (
                f"<code>{esc(rec.get('evidence_id'))}</code> "
                f"confidence={esc(rec.get('confidence_score')) or 'not_reported'} "
                f"PMID={esc(rec.get('pmid')) or 'not_reported'} "
                f"DOI={esc(rec.get('doi')) or 'not_reported'}"
            )
            context = (
                f"cancer={esc(rec.get('cancer_type')) or 'not_reported'}; "
                f"assay={esc(rec.get('assay_type')) or 'not_reported'}; "
                f"sample={esc(rec.get('sample_type')) or 'not_reported'}"
            )
            items.append(f"<li>{ids}<br><b>{title}</b><br><span>{context}</span><br><span>{evidence_text}</span></li>")
        return "".join(items) if items else "<li>No local KG evidence record attached to this variant.</li>"

    def rag_items(row: dict[str, Any]) -> str:
        raw = clean_text(row.get("rag_evidence_details"))
        if not raw:
            return "<li>RAG was not run or no RAG hit was mapped to this variant.</li>"
        try:
            records = json.loads(raw)
        except json.JSONDecodeError:
            records = []
        if not isinstance(records, list) or not records:
            return "<li>No direct RAG hit was mapped to this variant; use the scenario-level RAG context and KG provenance.</li>"
        items: list[str] = []
        for rec in records[:4]:
            if not isinstance(rec, dict):
                continue
            ids = (
                f"<code>{esc(rec.get('round_id'))}</code> "
                f"agent={esc(rec.get('agent')) or 'not_reported'} "
                f"score={esc(rec.get('score')) or 'not_reported'} "
                f"PMID={esc(rec.get('pmid')) or 'not_reported'} "
                f"DOI={esc(rec.get('doi')) or 'not_reported'}"
            )
            categories = ", ".join(rec.get("categories", [])) if isinstance(rec.get("categories"), list) else clean_text(rec.get("categories"))
            items.append(
                f"<li>{ids}<br><b>{esc(short_text(rec.get('title'), 180))}</b>"
                f"<br><span>match={esc(rec.get('match_reason'))}; categories={esc(categories)}</span>"
                f"<br><span>{esc(short_text(rec.get('text'), 240))}</span></li>"
            )
        return "".join(items) if items else "<li>No direct RAG hit was mapped to this variant.</li>"

    def path_items(row: dict[str, Any]) -> str:
        nodes = parse_json_list(row.get("graphrag_path_nodes"))
        edges = parse_json_list(row.get("graphrag_path_edges"))
        if not nodes and not edges:
            return "<li>No GraphRAG path was constructed for this variant.</li>"
        node_lines = []
        for node in nodes[:10]:
            label = esc(node.get("label"))
            name = esc(short_text(node.get("name"), 140))
            node_lines.append(f"<li><code>{label}</code> {name}</li>")
        if edges:
            node_lines.append(f"<li><code>Edges</code> {esc(len(edges))} traceability edge(s) connect variant, gene, evidence, rules, and missing context.</li>")
        return "".join(node_lines)

    def scenario_items(row: dict[str, Any]) -> str:
        raw = clean_text(row.get("caller_scenario_completion"))
        try:
            records = json.loads(raw) if raw else []
        except json.JSONDecodeError:
            records = []
        if not isinstance(records, list) or not records:
            return "<li>No caller-dependent MRD workflow was attached to this row.</li>"
        items = []
        for rec in records[:6]:
            if not isinstance(rec, dict):
                continue
            requires = rec.get("required_inputs_for_completion") or []
            requires_text = ", ".join(requires) if isinstance(requires, list) else clean_text(requires)
            items.append(
                f"<li><code>{esc(rec.get('scenario'))}</code><br>"
                f"<span>{esc(rec.get('vcf_only_contribution'))}</span><br>"
                f"<span>To complete source assignment: {esc(requires_text or 'expert review context')}</span></li>"
            )
        return "".join(items)

    def distribution_panel(title: str, counter: dict[str, Any] | Counter[str], denominator: int | None = None, limit: int = 7) -> str:
        items = [(clean_text(key), int(value)) for key, value in dict(counter).items() if clean_text(key)]
        items.sort(key=lambda item: item[1], reverse=True)
        if not items:
            return f'<section class="viz-panel"><div class="panel-head"><h2>{esc(title)}</h2></div><p class="muted">No data</p></section>'
        max_count = max(value for _, value in items) or 1
        denom = denominator or sum(value for _, value in items) or len(rows)
        bars = []
        for raw_label, value in items[:limit]:
            width = max(2.0, (value / max_count) * 100)
            bars.append(
                f"""
                <div class="bar-row">
                  <div class="bar-label">{esc(label(raw_label))}</div>
                  <div class="bar-track"><span class="bar-fill {lane_class(raw_label)}" style="width:{width:.1f}%"></span></div>
                  <div class="bar-value">{value:,}<small>{pct(value, denom)}</small></div>
                </div>
                """
            )
        return (
            f'<section class="viz-panel">'
            f'<div class="panel-head"><h2>{esc(title)}</h2><span>{sum(value for _, value in items):,}</span></div>'
            f'<div class="bar-list">{"".join(bars)}</div>'
            f'</section>'
        )

    def stacked_bar(counter: dict[str, Any] | Counter[str], denominator: int | None = None) -> str:
        items = [(clean_text(key), int(value)) for key, value in dict(counter).items() if clean_text(key)]
        items.sort(key=lambda item: item[1], reverse=True)
        denom = denominator or sum(value for _, value in items) or 1
        segments = []
        for raw_label, value in items:
            width = max(0.6, (value / denom) * 100)
            segments.append(
                f'<span class="stack-seg {lane_class(raw_label)}" title="{esc(label(raw_label))}: {value:,}" '
                f'style="width:{width:.2f}%"></span>'
            )
        return f'<div class="stacked-bar">{"".join(segments)}</div>'

    def metric(label_text: str, value: Any, sub: str = "", tone: str = "neutral") -> str:
        return (
            f'<div class="metric {tone}">'
            f'<span>{esc(label_text)}</span>'
            f'<b>{esc(value)}</b>'
            f'<small>{esc(sub)}</small>'
            f'</div>'
        )

    def af_histogram() -> str:
        bins = [
            (0.0, 0.02, "<2%"),
            (0.02, 0.05, "2-5%"),
            (0.05, 0.10, "5-10%"),
            (0.10, 0.20, "10-20%"),
            (0.20, math.inf, ">=20%"),
        ]
        counts = Counter()
        for row in rows:
            af = num(row.get("tumor_af"))
            for low, high, name in bins:
                if low <= af < high:
                    counts[name] += 1
                    break
        max_count = max(counts.values()) if counts else 1
        bars = []
        for _, _, name in bins:
            value = counts.get(name, 0)
            height = max(4, int((value / max_count) * 118)) if max_count else 4
            bars.append(
                f"""
                <div class="hist-col">
                  <div class="hist-bar" style="height:{height}px"><span>{value:,}</span></div>
                  <small>{esc(name)}</small>
                </div>
                """
            )
        return (
            '<section class="viz-panel">'
            '<div class="panel-head"><h2>VAF Distribution</h2><span>tumor AF</span></div>'
            f'<div class="histogram">{"".join(bars)}</div>'
            '</section>'
        )

    def af_score_heatmap() -> str:
        af_bins = [
            (0.0, 0.02, "<2%"),
            (0.02, 0.05, "2-5%"),
            (0.05, 0.10, "5-10%"),
            (0.10, 0.20, "10-20%"),
            (0.20, math.inf, ">=20%"),
        ]
        score_bins = [
            (0, 20, "0-20"),
            (20, 40, "20-40"),
            (40, 60, "40-60"),
            (60, 80, "60-80"),
            (80, 101, "80-100"),
        ]
        cells: dict[tuple[str, str], int] = {}
        for row in rows:
            af = num(row.get("tumor_af"))
            score = num(row.get("technical_score"))
            af_name = next(name for low, high, name in af_bins if low <= af < high)
            score_name = next(name for low, high, name in score_bins if low <= score < high)
            cells[(score_name, af_name)] = cells.get((score_name, af_name), 0) + 1
        max_count = max(cells.values()) if cells else 1
        body = []
        body.append('<div class="heat-cell heat-axis"></div>')
        for _, _, name in af_bins:
            body.append(f'<div class="heat-cell heat-axis">{esc(name)}</div>')
        for _, _, score_name in score_bins:
            body.append(f'<div class="heat-cell heat-axis">{esc(score_name)}</div>')
            for _, _, af_name in af_bins:
                value = cells.get((score_name, af_name), 0)
                alpha = 0.08 + (0.78 * (value / max_count if max_count else 0))
                text_class = " heat-strong" if alpha > 0.48 else ""
                body.append(
                    f'<div class="heat-cell{text_class}" style="background:rgba(0,109,119,{alpha:.3f})">'
                    f'{value:,}</div>'
                )
        return (
            '<section class="viz-panel">'
            '<div class="panel-head"><h2>AF x Score Heatmap</h2><span>count</span></div>'
            f'<div class="heatmap">{"".join(body)}</div>'
            '</section>'
        )

    def workflow_matrix() -> str:
        items = mqc.get("caller_workflow_coverage", [])
        if not items:
            return ""
        max_count = max(int(item.get("variant_count") or 0) for item in items) or 1
        rows_html = []
        for item in mqc.get("caller_workflow_coverage", []):
            requires = item.get("required_inputs_for_completion") or []
            requires_list = requires if isinstance(requires, list) else split_semicolon(requires)
            requires_chips = "".join(f'<span class="mini-chip">{esc(label(req))}</span>' for req in requires_list[:4])
            examples = item.get("example_variants") or []
            examples_text = ", ".join(examples[:4]) if isinstance(examples, list) else clean_text(examples)
            variant_count = int(item.get("variant_count") or 0)
            width = max(2, (variant_count / max_count) * 100)
            rows_html.append(
                f"""
                <tr>
                  <td><b>{esc(label(item.get('scenario')))}</b><small>{esc(label(item.get('status')))}</small></td>
                  <td>
                    <div class="table-bar"><span style="width:{width:.1f}%"></span></div>
                    <small>{variant_count:,} variants</small>
                  </td>
                  <td><b>{esc(item.get('rag_hit_count'))}</b><small>mapped RAG/rule hits</small></td>
                  <td><div class="chip-strip">{requires_chips}</div></td>
                  <td>
                    <details>
                      <summary>Details</summary>
                      <p>{esc(item.get('vcf_only_contribution'))}</p>
                      <p>{esc(item.get('why_caller_matters'))}</p>
                      <p>{esc(examples_text)}</p>
                    </details>
                  </td>
                </tr>
                """
            )
        return (
            '<section class="workflow-section">'
            '<div class="section-head"><h2>Caller Workflow Coverage</h2><span>VCF triggers + mapped context</span></div>'
            '<div class="table-wrap"><table class="workflow-table">'
            '<thead><tr><th>Workflow</th><th>Coverage</th><th>Mapped context</th><th>Missing inputs</th><th></th></tr></thead>'
            f'<tbody>{"".join(rows_html)}</tbody>'
            '</table></div>'
            '</section>'
        )

    def scenario_conclusion_panel() -> str:
        items = mqc.get("caller_workflow_coverage", [])
        if not items:
            return ""
        conclusion_by_scenario = {
            "tumor_informed_marker_tracking": "VCF/caller candidates trigger this workflow, but KG-backed tumor-informed MRD assignment requires re-mined tumor marker and concordance evidence",
            "clone_dynamics_or_recurrence_origin": "VCF/caller candidates preserve variant identity, but KG-backed clone-origin interpretation requires serial ctDNA and baseline tumor-clone evidence",
            "chip_or_hematopoietic_origin_audit": "CHIP/WBC review is triggered, but KG-backed downgrading requires extracted matched WBC/buffy coat or germline/hematopoietic filtering evidence",
            "low_frequency_serial_confirmation": "low-frequency review is triggered, but KG-backed interpretation requires extracted repeat/serial, orthogonal confirmation, and assay LOD evidence",
            "targeted_resistance_or_trial_evidence": "target/resistance review is triggered, but KG-backed actionability requires extracted functional annotation, therapy context, and guideline/database evidence",
            "fragmentomics_caller_discordance_audit": "caller-side discordance review is triggered, but KG-backed reconciliation requires extracted fragmentomics score and assay metadata evidence",
            "clinical_report_audit_trail": "the audit trail is structured from VCF fields, but KG-backed report-facing source assignment requires extracted signout-policy and clinical-context evidence",
        }
        fragments = []
        for item in items[:7]:
            scenario = clean_text(item.get("scenario"))
            variant_count = int(item.get("variant_count") or 0)
            rag_count = int(item.get("rag_hit_count") or 0)
            conclusion = conclusion_by_scenario.get(
                scenario,
                "current VCF/KG/RAG evidence is useful for workflow triage but not sufficient for KG-backed biological source assignment",
            )
            fragments.append(
                f"<span><b>{esc(label(scenario))}</b>: {esc(conclusion)} "
                f"<small>({variant_count:,} VCF-triggered variants; {rag_count} mapped RAG/rule hits)</small>.</span>"
            )
        return (
            '<section class="workflow-section scenario-conclusion">'
            '<div class="section-head"><h2>MRD 7-Scenario Workflow Interpretation</h2><span>workflow-derived; rerun needed for KG-backed view</span></div>'
            f'<p class="scenario-paragraph">{" ".join(fragments)}</p>'
            '</section>'
        )

    def row_signal(row: dict[str, Any], field: str, positive: set[str] | None = None) -> str:
        text = clean_text(row.get(field))
        active = bool(text)
        if positive is not None:
            active = text in positive
        css = lane_class(text) if active else "neutral"
        return f'<span class="dot {css}" title="{esc(label(text or "not_available"))}"></span>'

    def af_bar(row: dict[str, Any]) -> str:
        af = num(row.get("tumor_af"))
        width = min(100, max(2, (af / 0.20) * 100))
        return (
            f'<div class="af-cell"><div class="table-bar"><span style="width:{width:.1f}%"></span></div>'
            f'<small>{esc(fmt_af(row.get("tumor_af")))}</small></div>'
        )

    def score_bar(row: dict[str, Any]) -> str:
        score = num(row.get("technical_score"))
        width = min(100, max(2, score))
        return (
            f'<div class="score-cell"><div class="table-bar score"><span style="width:{width:.1f}%"></span></div>'
            f'<small>{score:.0f}</small></div>'
        )

    def variant_matrix() -> str:
        matrix_rows = []
        for row in display_rows[:14]:
            gene = clean_text(row.get("genes")) or "gene_not_annotated"
            matrix_rows.append(
                f"""
                <tr>
                  <td><code>{esc(row.get('variant_id'))}</code><small>{esc(short_text(gene, 42))}</small></td>
                  <td>{af_bar(row)}</td>
                  <td><b>{esc(row.get('tumor_vd'))}/{esc(row.get('tumor_dp'))}</b><small>VD/DP</small></td>
                  <td>{score_bar(row)}</td>
                  <td>{esc(label(row.get('qc_tier')))}</td>
                  <td>{esc(label(row.get('source_trace_status')))}</td>
                  <td>
                    <div class="dot-row">
                      {row_signal(row, 'kg_match_status', {'gene_evidence_match'})}
                      {row_signal(row, 'rag_evidence_status', {'rag_variant_or_direct_context_integrated'})}
                      {row_signal(row, 'source_assignment_confidence')}
                    </div>
                  </td>
                  <td>{esc(action_label(row))}</td>
                </tr>
                """
            )
        return (
            '<section class="workflow-section">'
            '<div class="section-head"><h2>Top Variant Matrix</h2><span>AF, support, source trace</span></div>'
            '<div class="table-wrap"><table class="variant-table">'
            '<thead><tr><th>Variant / gene</th><th>AF</th><th>Depth</th><th>Score</th><th>QC</th><th>Source</th><th>KG/RAG/conf.</th><th>Next</th></tr></thead>'
            f'<tbody>{"".join(matrix_rows)}</tbody>'
            '</table></div>'
            '</section>'
        )

    def evidence_drawers() -> str:
        cards: list[str] = []
        for row in display_rows[:12]:
            gene = clean_text(row.get("genes")) or "gene_not_annotated"
            lane = "".join(
                [
                    pill("Caller QC", row.get("qc_tier")),
                    pill("KG", row.get("kg_evidence_ids") or row.get("kg_match_status")),
                    pill("Risk", row.get("qc_flags") or "no_major_flag"),
                    pill("Source trace", row.get("source_trace_status")),
                    pill("RAG", row.get("rag_evidence_status")),
                    pill("Action", row.get("recommended_review_action")),
                ]
            )
            cards.append(
                f"""
                <details class="variant-drawer status-{lane_class(row.get('source_trace_status'))}">
                  <summary>
                    <span><code>{esc(row.get('variant_id'))}</code> {esc(short_text(gene, 70))}</span>
                    <span>{esc(fmt_af(row.get('tumor_af')))} | {esc(action_label(row))}</span>
                  </summary>
                  <div class="lane">{lane}</div>
                  <div class="drawer-grid">
                    <section class="drawer-panel">
                      <h3>Source Trace</h3>
                      <p>{esc(row.get('source_trace_explanation') or row.get('variant_nlp_qc_conclusion'))}</p>
                      <p><b>Confidence:</b> {esc(label(row.get('source_assignment_confidence')))}</p>
                      <p><b>Next:</b> {esc(row.get('next_confirmation_step'))}</p>
                    </section>
                    <section class="drawer-panel">
                      <h3>RAG / GraphRAG</h3>
                      <p>{esc(row.get('rag_support_summary'))}</p>
                      <p>{esc(row.get('graphrag_path_summary'))}</p>
                    </section>
                    <section class="drawer-panel">
                      <h3>Unresolved Context</h3>
                      <p>{esc(row.get('missing_context_summary'))}</p>
                      <p>{esc(row.get('rag_uncertainty_checklist'))}</p>
                    </section>
                    <section class="drawer-panel">
                      <h3>QC Flags</h3>
                      <p>{esc(row.get('inner_qc_conclusion'))}</p>
                      <p>{esc(row.get('source_trace_reason_codes'))}</p>
                    </section>
                  </div>
                  <div class="evidence-grid">
                    <details class="kg-details">
                      <summary>RAG agent evidence</summary>
                      <ul>{rag_items(row)}</ul>
                    </details>
                    <details class="kg-details">
                      <summary>GraphRAG path provenance</summary>
                      <ul>{path_items(row)}</ul>
                    </details>
                    <details class="kg-details">
                      <summary>Caller-dependent workflow coverage</summary>
                      <ul>{scenario_items(row)}</ul>
                    </details>
                    <details class="kg-details">
                      <summary>KG provenance evidence</summary>
                      <ul>{evidence_items(row)}</ul>
                    </details>
                  </div>
                </details>
                """
            )
        return (
            '<section class="workflow-section">'
            '<div class="section-head"><h2>Evidence Drawers</h2><span>collapsed by default</span></div>'
            f'{"".join(cards)}'
            '</section>'
        )

    total_rows = len(rows)
    tier_counts = summary.get("tier_counts") or Counter(clean_text(row.get("qc_tier")) for row in rows)
    kg_counts = summary.get("kg_match_counts") or Counter(clean_text(row.get("kg_match_status")) for row in rows)
    source_counts = mqc.get("source_trace_status_counts") or Counter(clean_text(row.get("source_trace_status")) for row in rows)
    confidence_counts = Counter(clean_text(row.get("source_assignment_confidence")) or "not_available" for row in rows)
    rag_counts = Counter(clean_text(row.get("rag_evidence_status")) or "not_available" for row in rows)
    chrom_counts = Counter(clean_text(row.get("chrom")) or "NA" for row in rows)
    type_counts = summary.get("type_counts") or Counter(clean_text(row.get("type")) for row in rows)
    flag_counts = Counter(flag for row in rows for flag in split_semicolon(row.get("qc_flags")) if flag)
    af_values = [num(row.get("tumor_af")) for row in rows if clean_text(row.get("tumor_af"))]
    score_values = [num(row.get("technical_score")) for row in rows if clean_text(row.get("technical_score"))]
    tech_rescue_count = int(dict(source_counts).get("technical_rescue_required_before_source_assignment", 0))
    kg_matched_count = int(dict(kg_counts).get("gene_evidence_match", 0))
    rag_direct_count = int(rag_counts.get("rag_variant_or_direct_context_integrated", 0))
    strong_count = int(dict(tier_counts).get("strong_technical_candidate", 0))
    dominant_source = (mqc.get("case_source_trace_synthesis") or {}).get("dominant_source_explanation", "not_available")

    metrics = "".join(
        [
            metric("VCF records", f"{total_rows:,}", "caller-derived rows", "neutral"),
            metric("Technical rescue", f"{tech_rescue_count:,}", pct(tech_rescue_count, total_rows), "warn"),
            metric("Strong candidates", f"{strong_count:,}", pct(strong_count, total_rows), "good"),
            metric("KG gene match", f"{kg_matched_count:,}", pct(kg_matched_count, total_rows), "neutral"),
            metric("CHIP/WBC review", f"{int(mqc.get('chip_review_variant_count') or 0):,}", "variants surfaced", "bad"),
            metric("Direct RAG rows", f"{rag_direct_count:,}", f"{int(summary.get('source_trace_rag', {}).get('variant_specific_rag_rows', 0) or 0):,} variant queries", "neutral"),
            metric("Median AF", fmt_af(median(af_values) if af_values else 0), f"max {fmt_af(max(af_values) if af_values else 0)}", "neutral"),
            metric("Median score", f"{median(score_values):.0f}" if score_values else "0", f"max {max(score_values):.0f}" if score_values else "max 0", "neutral"),
        ]
    )

    chart_panels = "".join(
        [
            distribution_panel("QC Tier", tier_counts, total_rows),
            distribution_panel("Source Trace", source_counts, total_rows),
            distribution_panel("Source Confidence", confidence_counts, total_rows),
            distribution_panel("KG Evidence", kg_counts, total_rows),
            distribution_panel("RAG Status", rag_counts, total_rows),
            af_histogram(),
            af_score_heatmap(),
            distribution_panel("QC Flags", flag_counts, total_rows, limit=8),
            distribution_panel("Chromosome Load", chrom_counts, total_rows, limit=10),
        ]
    )

    type_stack = stacked_bar(type_counts, total_rows)
    source_stack = stacked_bar(source_counts, total_rows)

    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>MRD-KGQC Source-Tracing Dashboard</title>
  <style>
    :root {{
      color-scheme: light;
      --bg: #f4f6f3;
      --panel: #ffffff;
      --ink: #20262d;
      --muted: #65717d;
      --soft: #f8faf8;
      --line: #d7ded8;
      --line-strong: #b8c4bd;
      --accent: #006d77;
      --accent-soft: #e1f2f1;
      --plum: #7a3e58;
      --plum-soft: #f3e7ed;
      --gold: #9a6a00;
      --gold-soft: #fff2cf;
      --good: #dff1e4;
      --good-ink: #235a33;
      --warn: #fff2cf;
      --warn-ink: #7a5200;
      --bad: #f7dfe5;
      --bad-ink: #8a2441;
      --neutral: #edf1f3;
      --neutral-ink: #394650;
      --olive: #667b3a;
      --olive-soft: #edf3df;
    }}
    * {{
      box-sizing: border-box;
    }}
    body {{
      margin: 0;
      font-family: "Segoe UI", Arial, Helvetica, sans-serif;
      color: var(--ink);
      background: var(--bg);
      font-size: 14px;
    }}
    main {{
      max-width: 1240px;
      margin: 0 auto;
      padding: 24px;
    }}
    header {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 22px;
      margin-bottom: 18px;
      box-shadow: 0 10px 28px rgba(32, 38, 45, 0.06);
    }}
    h1, h2, h3 {{
      margin: 0;
      letter-spacing: 0;
    }}
    h1 {{
      font-size: 26px;
      line-height: 1.15;
    }}
    h2 {{
      font-size: 18px;
      line-height: 1.25;
    }}
    h3 {{
      font-size: 14px;
      line-height: 1.25;
    }}
    header p, .variant-head p, .conclusion, .muted {{
      color: var(--muted);
      line-height: 1.5;
    }}
    .topline {{
      display: flex;
      justify-content: space-between;
      gap: 18px;
      align-items: flex-start;
      margin-bottom: 14px;
    }}
    .topline p {{
      max-width: 760px;
      margin: 8px 0 0;
    }}
    .case-badge {{
      background: var(--accent-soft);
      color: var(--accent);
      border: 1px solid #bddbd8;
      border-radius: 8px;
      padding: 8px 10px;
      font-weight: 700;
      white-space: nowrap;
    }}
    .mini-strip {{
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 8px;
      margin-top: 10px;
    }}
    .mini-strip > div {{
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 10px;
      background: var(--soft);
    }}
    .mini-strip b {{
      display: block;
      margin-bottom: 8px;
    }}
    .workflow-section {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 16px;
      margin-bottom: 16px;
      box-shadow: 0 8px 22px rgba(32, 38, 45, 0.045);
    }}
    .section-head, .panel-head {{
      display: flex;
      justify-content: space-between;
      gap: 12px;
      align-items: baseline;
      margin-bottom: 12px;
    }}
    .section-head span, .panel-head span {{
      color: var(--muted);
      font-size: 12px;
      font-weight: 700;
    }}
    .kg-details {{
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 11px 12px;
      background: #fff;
      color: var(--ink);
    }}
    .kg-details summary {{
      cursor: pointer;
      font-weight: bold;
      color: var(--ink);
      list-style-position: outside;
    }}
    .kg-details li {{
      margin-bottom: 10px;
      line-height: 1.45;
    }}
    .kg-details span {{
      color: var(--muted);
    }}
    .metric-grid {{
      display: grid;
      grid-template-columns: repeat(4, minmax(0, 1fr));
      gap: 10px;
      margin-top: 16px;
    }}
    .metric {{
      min-height: 92px;
      border-radius: 8px;
      padding: 13px;
      border: 1px solid var(--line);
      background: var(--soft);
      display: flex;
      flex-direction: column;
      justify-content: space-between;
    }}
    .metric span {{
      color: var(--muted);
      font-weight: 700;
      font-size: 12px;
      text-transform: uppercase;
    }}
    .metric b {{
      font-size: 28px;
      line-height: 1;
    }}
    .metric small {{
      color: var(--muted);
      font-size: 12px;
    }}
    .metric.good {{ border-color: #a8cfb4; }}
    .metric.warn {{ border-color: #e3c66f; }}
    .metric.bad {{ border-color: #ddb1bf; }}
    .viz-grid {{
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 10px;
      margin: 12px 0 18px;
    }}
    .viz-panel {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 14px;
      box-shadow: 0 8px 22px rgba(32, 38, 45, 0.035);
    }}
    .bar-list {{
      display: grid;
      gap: 9px;
    }}
    .bar-row {{
      display: grid;
      grid-template-columns: minmax(150px, 1.3fr) minmax(120px, 2fr) 82px;
      gap: 10px;
      align-items: center;
    }}
    .bar-label {{
      font-weight: 700;
      font-size: 12px;
      overflow-wrap: anywhere;
    }}
    .bar-track, .table-bar {{
      height: 10px;
      border-radius: 999px;
      background: #e8eeeb;
      overflow: hidden;
      border: 1px solid #dce4df;
    }}
    .bar-fill, .table-bar span {{
      display: block;
      height: 100%;
      border-radius: inherit;
      background: var(--accent);
    }}
    .bar-fill.good {{ background: #3d8b54; }}
    .bar-fill.warn {{ background: var(--gold); }}
    .bar-fill.bad {{ background: var(--plum); }}
    .bar-fill.neutral {{ background: var(--accent); }}
    .bar-value {{
      text-align: right;
      font-weight: 800;
      font-variant-numeric: tabular-nums;
    }}
    .bar-value small {{
      display: block;
      color: var(--muted);
      font-weight: 600;
      margin-top: 2px;
    }}
    .stacked-bar {{
      display: flex;
      height: 16px;
      overflow: hidden;
      border-radius: 999px;
      border: 1px solid var(--line);
      background: #e8eeeb;
    }}
    .stack-seg.good {{ background: #3d8b54; }}
    .stack-seg.warn {{ background: var(--gold); }}
    .stack-seg.bad {{ background: var(--plum); }}
    .stack-seg.neutral {{ background: var(--accent); }}
    .histogram {{
      height: 160px;
      display: grid;
      grid-template-columns: repeat(5, 1fr);
      gap: 12px;
      align-items: end;
      padding-top: 6px;
    }}
    .hist-col {{
      display: grid;
      gap: 7px;
      align-items: end;
      text-align: center;
      color: var(--muted);
      font-weight: 700;
    }}
    .hist-bar {{
      min-height: 4px;
      background: var(--olive);
      border-radius: 6px 6px 2px 2px;
      display: flex;
      align-items: flex-start;
      justify-content: center;
      color: #fff;
      font-size: 11px;
      padding-top: 4px;
      overflow: hidden;
    }}
    .heatmap {{
      display: grid;
      grid-template-columns: 64px repeat(5, 1fr);
      gap: 4px;
    }}
    .heat-cell {{
      min-height: 32px;
      border-radius: 6px;
      display: flex;
      align-items: center;
      justify-content: center;
      font-weight: 800;
      color: var(--ink);
      font-variant-numeric: tabular-nums;
      font-size: 12px;
    }}
    .heat-axis {{
      background: var(--soft);
      color: var(--muted);
      border: 1px solid var(--line);
      font-weight: 700;
    }}
    .heat-strong {{
      color: #fff;
    }}
    .table-wrap {{
      overflow-x: auto;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      min-width: 920px;
    }}
    th {{
      color: var(--muted);
      text-transform: uppercase;
      font-size: 11px;
      letter-spacing: 0;
      text-align: left;
      padding: 0 10px 8px;
    }}
    td {{
      border-top: 1px solid var(--line);
      padding: 11px 10px;
      vertical-align: middle;
      font-size: 13px;
    }}
    td small {{
      display: block;
      color: var(--muted);
      margin-top: 4px;
    }}
    .variant-table code {{
      font-weight: 800;
      color: var(--ink);
      white-space: nowrap;
    }}
    .table-bar.score span {{
      background: var(--gold);
    }}
    .af-cell, .score-cell {{
      min-width: 120px;
    }}
    .dot-row {{
      display: flex;
      gap: 6px;
    }}
    .dot {{
      width: 13px;
      height: 13px;
      border-radius: 999px;
      display: inline-block;
      border: 1px solid rgba(32, 38, 45, 0.16);
    }}
    .dot.good {{ background: #3d8b54; }}
    .dot.warn {{ background: var(--gold); }}
    .dot.bad {{ background: var(--plum); }}
    .dot.neutral {{ background: #b9c4be; }}
    .chip-strip {{
      display: flex;
      gap: 5px;
      flex-wrap: wrap;
    }}
    .mini-chip {{
      display: inline-flex;
      align-items: center;
      min-height: 24px;
      border-radius: 999px;
      padding: 3px 8px;
      background: var(--neutral);
      color: var(--neutral-ink);
      font-size: 12px;
      font-weight: 700;
    }}
    .scenario-conclusion {{
      border-left: 4px solid var(--accent);
    }}
    .scenario-paragraph {{
      margin: 0;
      color: var(--ink);
      line-height: 1.65;
    }}
    .scenario-paragraph span {{
      display: inline;
    }}
    .scenario-paragraph small {{
      color: var(--muted);
      font-weight: 700;
    }}
    .variant-drawer {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      margin-top: 10px;
      position: relative;
      overflow: hidden;
    }}
    .variant-drawer::before {{
      content: "";
      position: absolute;
      inset: 0 auto 0 0;
      width: 5px;
      background: var(--line-strong);
    }}
    .variant-drawer.status-good::before {{ background: #3d8b54; }}
    .variant-drawer.status-warn::before {{ background: #b77900; }}
    .variant-drawer.status-bad::before {{ background: var(--plum); }}
    .variant-drawer.status-neutral::before {{ background: var(--accent); }}
    .variant-drawer > * {{
      position: relative;
    }}
    .variant-drawer > summary {{
      display: flex;
      justify-content: space-between;
      gap: 12px;
      padding: 13px 14px 13px 18px;
      cursor: pointer;
      font-weight: 800;
      list-style-position: outside;
    }}
    .variant-drawer > summary span:last-child {{
      color: var(--muted);
      font-weight: 700;
      white-space: nowrap;
    }}
    .drawer-grid {{
      display: grid;
      grid-template-columns: repeat(4, minmax(0, 1fr));
      gap: 10px;
      padding: 0 14px 14px 18px;
    }}
    .drawer-panel {{
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 11px;
      background: var(--soft);
      line-height: 1.45;
      overflow-wrap: anywhere;
    }}
    .drawer-panel p {{
      margin: 8px 0 0;
      color: var(--muted);
    }}
    .lane {{
      display: grid;
      grid-template-columns: repeat(6, minmax(110px, 1fr));
      gap: 8px;
      margin: 0 14px 12px 18px;
      position: relative;
    }}
    .pill {{
      min-height: 64px;
      border-radius: 8px;
      padding: 9px 10px;
      border: 1px solid rgba(32, 38, 45, 0.09);
      overflow-wrap: anywhere;
      font-size: 12px;
      line-height: 1.35;
      display: flex;
      flex-direction: column;
      justify-content: space-between;
    }}
    .pill-label {{
      font-size: 11px;
      text-transform: uppercase;
      letter-spacing: 0;
      font-weight: 700;
      opacity: 0.78;
    }}
    .pill-value {{
      margin-top: 6px;
      font-weight: 600;
    }}
    .evidence-grid {{
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 10px;
      margin-top: 12px;
    }}
    .full-span {{
      grid-column: 1 / -1;
    }}
    .good {{ background: var(--good); color: var(--good-ink); }}
    .warn {{ background: var(--warn); color: var(--warn-ink); }}
    .bad {{ background: var(--bad); color: var(--bad-ink); }}
    .neutral {{ background: var(--neutral); color: var(--neutral-ink); }}
    @media (max-width: 760px) {{
      main {{ padding: 14px; }}
      .topline {{ display: block; }}
      .case-badge {{ display: inline-block; margin-top: 10px; white-space: normal; }}
      .metric-grid, .viz-grid, .mini-strip, .drawer-grid {{ grid-template-columns: 1fr; }}
      .lane {{ grid-template-columns: 1fr; margin-right: 10px; }}
      .evidence-grid {{ grid-template-columns: 1fr; }}
      .variant-head {{ display: block; }}
      .status-code {{ display: inline-block; margin-top: 8px; white-space: normal; }}
      .variant-drawer > summary {{ display: block; }}
      .variant-drawer > summary span:last-child {{ display: block; margin-top: 6px; white-space: normal; }}
      .bar-row {{ grid-template-columns: 1fr; gap: 5px; }}
      .bar-value {{ text-align: left; }}
    }}
  </style>
</head>
<body>
<main>
  <header>
    <div class="topline">
      <div>
        <h1>MRD-KGQC Source-Tracing Dashboard</h1>
        <p>{esc(short_text(dominant_source, 190))}</p>
      </div>
      <div class="case-badge">{esc(label(mqc.get('source_trace_status')))}</div>
    </div>
    <div class="mini-strip">
      <div><b>Source mix</b>{source_stack}</div>
      <div><b>Variant types</b>{type_stack}</div>
    </div>
    <div class="metric-grid">{metrics}</div>
  </header>
  <div class="viz-grid">{chart_panels}</div>
  {workflow_matrix()}
  {scenario_conclusion_panel()}
  {variant_matrix()}
  {evidence_drawers()}
  <section class="workflow-section">
    <div class="section-head"><h2>Case Metadata</h2><span>review boundary</span></div>
    <div class="chip-strip">
      <span class="mini-chip">mode: {esc(label(mqc.get('evidence_mode')))}</span>
      <span class="mini-chip">gate: {esc(label(gate.get('gate_class')))}</span>
      <span class="mini-chip">permission: {esc(label(gate.get('entry_permission')))}</span>
      <span class="mini-chip">GraphRAG paths: {esc(summary.get('source_trace_rag', {}).get('graphrag_path_rows', 0))}</span>
      <span class="mini-chip">VCF: {esc(vcf_path.name or vcf_path)}</span>
    </div>
  </section>
</main>
</body>
</html>
"""
    path.write_text(html_text, encoding="utf-8")


def build_source_trace_rag_payload(summary: dict[str, Any], rows: list[dict[str, Any]]) -> dict[str, str]:
    mqc = summary["molecular_evidence_qc"]
    top_items = mqc.get("top_variant_conclusions", [])
    genes: list[str] = []
    variants: list[str] = []
    statuses: list[str] = []
    actions: list[str] = []
    for item in top_items:
        variants.append(clean_text(item.get("variant_id")))
        for gene in split_semicolon(item.get("genes")):
            genes.append(gene)
        statuses.append(clean_text(item.get("source_trace_status")))
        actions.append(clean_text(item.get("recommended_review_action")))
    if not genes:
        gene_counter: Counter[str] = Counter()
        for row in rows:
            for gene in split_semicolon(row.get("genes")):
                gene_counter[gene] += 1
        genes = [gene for gene, _ in gene_counter.most_common(12)]

    query = (
        "MRD caller-derived signal source tracing. Explain whether variant-level caller signals are traceable to "
        "tumor-informed MRD, clonal recurrence, CHIP/hematopoietic origin, germline, technical artifact, transient low-support signal, "
        "or unresolved missing context. Retrieve evidence for matched WBC/buffy coat review, serial sampling, tumor concordance, assay LOD, "
        "low VAF, and local alignment artifact review."
    )

    return {
        "query": query,
        "source_trace_status": clean_text(mqc.get("source_trace_status")),
        "evidence_mode": clean_text(mqc.get("evidence_mode")),
        "nlp_qc_conclusion": clean_text(mqc.get("nlp_qc_conclusion")),
        "gene": " ".join(sorted(set(genes))[:20]),
        "variant": " ".join(item for item in variants[:12] if item),
        "mrd_status": "not_classified_source_tracing_only",
        "matched_wbc_status": "needs_review" if summary.get("matched_control_signal_count", 0) > 0 else "not_provided",
        "tumor_tissue_status": "not_provided",
        "assay_type": "caller-derived VCF",
        "sample_type": "not_provided",
        "timepoint": "not_provided",
        "source_trace_actions": " ".join(sorted(set(actions))[:12]),
    }


def write_source_trace_rag_markdown(path: Path, result: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# MRD-KGQC Caller Signal Source-Tracing RAG Packet",
        "",
        f"- Generated: {result['generated_at']}",
        f"- Case text: {result['case_text']}",
        f"- Chunks JSONL: `{result['chunks_jsonl']}`",
        f"- Rule DB: `{result['rule_db']}`",
        "",
        "## Source-Trace Question",
        "",
        result["case"].get("query", ""),
        "",
        "## Current VCF Source-Trace Summary",
        "",
        f"- Source trace status: `{result['case'].get('source_trace_status', '')}`",
        f"- Evidence mode: `{result['case'].get('evidence_mode', '')}`",
        f"- RAG mode: `{result.get('rag_mode', 'case')}`",
        f"- Genes queried: `{result['case'].get('gene', '')}`",
        f"- Top variants queried: `{result['case'].get('variant', '')}`",
        f"- NLP source explanation: {result['case'].get('nlp_qc_conclusion', '')}",
        "",
        "## How To Use",
        "",
        "This packet retrieves evidence and counter-evidence for source tracing. It does not classify MRD positivity or recommend treatment.",
        "",
    ]
    for round_result in result["rounds"]:
        lines.extend(["", f"## {round_result['round_id']} ({round_result['agent']})", "", round_result["purpose"], ""])
        lines.append(f"Query: `{round_result['query']}`")
        lines.append("")
        for index, hit in enumerate(round_result.get("hits", []), start=1):
            meta = hit["metadata"]
            lines.append(f"{index}. score={hit['score']} `{hit['chunk_id']}`")
            lines.append(f"   Title: {hit['title']}")
            lines.append(f"   DOI: {meta.get('doi')} | PMID: {meta.get('pmid')} | confidence: {meta.get('confidence_score')}")
            lines.append(f"   Text: {short_text(hit['text'], 700)}")
        if round_result.get("rule_hits"):
            lines.extend(["", "Structured rule hits:", ""])
            for index, row in enumerate(round_result["rule_hits"], start=1):
                lines.append(
                    f"{index}. `{row.get('rule_id')}` {row.get('rule_type')} {row.get('cancer_type')} "
                    f"{row.get('feature_a_label') or row.get('gene_a') or row.get('feature_a_id')} / "
                    f"{row.get('feature_b_label') or row.get('gene_b') or row.get('feature_b_id')} "
                    f"pair={row.get('pair_count')}/{row.get('sample_count')}"
                )
    if result.get("variant_rag_results"):
        lines.extend(["", "## Variant-Specific RAG", ""])
        for variant_result in result.get("variant_rag_results", [])[:24]:
            lines.extend(
                [
                    f"### {variant_result.get('variant_id')}",
                    "",
                    f"- Genes: `{', '.join(variant_result.get('genes', []))}`",
                    f"- Selection: {variant_result.get('selection_score')} ({', '.join(variant_result.get('selection_reasons', []))})",
                    f"- Query: `{variant_result.get('query')}`",
                    "",
                ]
            )
            for round_result in variant_result.get("rounds", [])[:3]:
                hit_count = len(round_result.get("hits", []) or [])
                rule_count = len(round_result.get("rule_hits", []) or [])
                lines.append(f"- `{round_result.get('round_id')}`: {hit_count} RAG hit(s), {rule_count} rule hit(s)")
            lines.append("")
    if result.get("case_source_trace_synthesis"):
        syn = result["case_source_trace_synthesis"]
        lines.extend(["", "## Final Deterministic Source-Tracing Synthesis", ""])
        lines.append(f"- Dominant explanation: {syn.get('dominant_source_explanation')}")
        lines.append(f"- Workflow readiness: `{syn.get('workflow_readiness')}`")
        eq = syn.get("evidence_quality", {})
        if isinstance(eq, dict):
            lines.append(f"- Evidence quality: case rounds={eq.get('case_round_count')}, variant RAG queries={eq.get('variant_specific_rag_queries')}, synthesized rows={eq.get('variant_rows_with_synthesis')}")
        if syn.get("unresolved_context"):
            lines.append("- Unresolved context:")
            for item in syn.get("unresolved_context", [])[:8]:
                lines.append(f"  - {item}")
        if syn.get("selected_variant_summaries"):
            lines.append("- Selected variant source explanations:")
            for item in syn.get("selected_variant_summaries", [])[:8]:
                lines.append(
                    f"  - `{item.get('variant_id')}` {item.get('genes') or ''}: "
                    f"{short_text(item.get('source_explanation'), 220)}"
                )
    lines.extend(
        [
            "",
            "## Source-Tracing Answer Shape For Expert Review",
            "",
            "- source_trace_answer",
            "- most_likely_source_explanation",
            "- evidence_supporting_tumor_origin",
            "- evidence_supporting_CHIP_or_hematopoietic_origin",
            "- evidence_supporting_technical_artifact",
            "- unresolved_missing_context",
            "- recommended_next_review_step",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def classify_rag_hit(hit: dict[str, Any]) -> list[str]:
    text = " ".join(
        [
            clean_text(hit.get("title")),
            clean_text(hit.get("text")),
            json.dumps(hit.get("metadata", {}), ensure_ascii=False),
        ]
    ).lower()
    categories = []
    for category, terms in RAG_CATEGORY_TERMS.items():
        if any(term in text for term in terms):
            categories.append(category)
    return categories or ["general_source_trace"]


def compact_rag_hit(
    round_result: dict[str, Any],
    hit: dict[str, Any],
    match_reason: str,
    scope: str = "case",
    variant_id: str = "",
) -> dict[str, Any]:
    meta = hit.get("metadata", {})
    if not isinstance(meta, dict):
        meta = {}
    return {
        "scope": scope,
        "variant_id": variant_id,
        "round_id": clean_text(round_result.get("round_id")),
        "agent": clean_text(round_result.get("agent")),
        "match_reason": match_reason,
        "score": hit.get("score"),
        "chunk_id": clean_text(hit.get("chunk_id")),
        "evidence_id": clean_text(meta.get("evidence_id")),
        "title": clean_text(hit.get("title")),
        "doi": clean_text(meta.get("doi")),
        "pmid": clean_text(meta.get("pmid")),
        "confidence_score": clean_text(meta.get("confidence_score")),
        "categories": classify_rag_hit(hit),
        "text": short_text(hit.get("text"), 360),
    }


def variant_rag_priority(row: dict[str, Any]) -> tuple[float, list[str]]:
    reasons: list[str] = []
    score = 0.0
    if clean_text(row.get("kg_match_status")) in {"gene_evidence_match", "coordinate_evidence_match"}:
        score += 120
        reasons.append("kg_matched")
    source_status = clean_text(row.get("source_trace_status"))
    if source_status == "traceable_caller_evidence_candidate":
        score += 55
        reasons.append("traceable_candidate")
    if source_status == "traceable_caller_signal_needs_serial_confirmation":
        score += 35
        reasons.append("serial_confirmation_needed")
    if source_status == "hematopoietic_or_CHIP_origin_suspected":
        score += 80
        reasons.append("chip_or_hematopoietic_review")
    if source_status in {"technical_artifact_origin_suspected", "technical_rescue_required_before_source_assignment"}:
        score += 45
        reasons.append("technical_rescue_or_artifact_review")
    if "chip_gene_overlap" in clean_text(row.get("source_trace_reason_codes")):
        score += 70
        reasons.append("chip_gene_overlap")
    if "targeted_resistance_or_trial_evidence" in split_semicolon(row.get("caller_dependent_mrd_scenarios")):
        score += 55
        reasons.append("targeted_or_resistance_gene")
    if "low_variant_support" in clean_text(row.get("source_trace_reason_codes")):
        score += 25
        reasons.append("low_frequency_confirmation")
    af = row.get("tumor_af") or 0
    technical_score = row.get("technical_score") or 0
    try:
        score += min(float(af) * 100, 25)
    except (TypeError, ValueError):
        pass
    try:
        score += min(float(technical_score), 100) / 10
    except (TypeError, ValueError):
        pass
    return score, reasons or ["representative_caller_signal"]


def select_variant_rag_rows(rows: list[dict[str, Any]], limit: int) -> list[dict[str, Any]]:
    if limit <= 0:
        return []
    selected: list[dict[str, Any]] = []
    seen: set[str] = set()
    for row in sorted(rows, key=lambda item: variant_rag_priority(item)[0], reverse=True):
        variant_id = clean_text(row.get("variant_id"))
        if not variant_id or variant_id in seen:
            continue
        selected.append(row)
        seen.add(variant_id)
        if len(selected) >= limit:
            break
    return selected


def build_variant_rag_payload(row: dict[str, Any], summary: dict[str, Any]) -> dict[str, str]:
    genes = clean_text(row.get("genes"))
    kg_ids = clean_text(row.get("kg_evidence_ids"))
    reason_codes = clean_text(row.get("source_trace_reason_codes"))
    query = (
        "MRD caller-derived variant source tracing. Retrieve evidence that helps explain whether this exact caller signal "
        "is tumor-derived MRD/clonal recurrence, CHIP or hematopoietic origin, germline, technical artifact, transient "
        "low-support signal, or unresolved context."
    )
    return {
        "query": query,
        "variant_id": clean_text(row.get("variant_id")),
        "gene": genes,
        "variant_type": clean_text(row.get("type")),
        "tumor_af": fmt_af(row.get("tumor_af")),
        "tumor_depth": clean_text(row.get("tumor_dp")),
        "tumor_variant_depth": clean_text(row.get("tumor_vd")),
        "qc_tier": clean_text(row.get("qc_tier")),
        "source_trace_status": clean_text(row.get("source_trace_status")),
        "reason_codes": reason_codes,
        "caller_scenarios": clean_text(row.get("caller_dependent_mrd_scenarios")),
        "kg_evidence_ids": kg_ids,
        "kg_match_status": clean_text(row.get("kg_match_status")),
        "matched_wbc_status": "needs_review" if summary.get("matched_control_signal_count", 0) > 0 else "not_provided",
        "tumor_tissue_status": "not_provided",
        "serial_sampling_status": "not_provided",
        "assay_lod_status": "not_provided",
        "needed_evidence": "matched WBC buffy coat tumor tissue concordance serial ctDNA low VAF assay LOD local alignment artifact germline",
    }


def payload_text(payload: dict[str, str]) -> str:
    return " | ".join(f"{key}: {value}" for key, value in payload.items() if clean_text(value))


def collect_genes_from_rows(rows: list[dict[str, Any]]) -> set[str]:
    genes: set[str] = set()
    for row in rows:
        genes.update(gene.upper() for gene in split_semicolon(row.get("genes")))
    return genes


def build_variant_rag_results(
    rows: list[dict[str, Any]],
    summary: dict[str, Any],
    args: argparse.Namespace,
    chunks_path: Path,
) -> list[dict[str, Any]]:
    if args.rag_mode not in {"variant", "complete"}:
        return []
    from mrdskill_query_mrd_kg_rag import search_rule_db
    from mrdskill_run_adjudication_rag import retrieve_rounds_streaming

    selected_rows = select_variant_rag_rows(rows, args.rag_variant_limit)
    selected_genes = {
        gene.upper()
        for row in selected_rows
        for gene in split_semicolon(row.get("genes"))
    }
    results: list[dict[str, Any]] = []
    max_chunks = 0 if args.rag_full_scan else args.rag_max_chunks
    overall_gene_set = collect_genes_from_rows(rows)
    overall_genes = sorted(overall_gene_set)
    for row in selected_rows:
        priority_score, reasons = variant_rag_priority(row)
        payload = build_variant_rag_payload(row, summary)
        query_text = payload_text(payload)
        genes = [gene.upper() for gene in split_semicolon(row.get("genes"))]
        rounds = retrieve_rounds_streaming(
            chunks_path,
            query_text,
            args.rag_top_k_per_round,
            max_chunks,
            boost_terms=genes or overall_genes[:20],
            min_score=args.rag_min_score,
            denoise=args.rag_denoise,
        )
        for round_result in rounds:
            round_result["rule_hits"] = search_rule_db(
                artifact_path(args.rag_rule_db),
                round_result["query"],
                args.rag_rule_top_k_per_round,
                allowed_genes=genes,
                blocked_genes=overall_gene_set - selected_genes,
                require_allowed_gene=bool(genes),
            )
        results.append(
            {
                "variant_id": clean_text(row.get("variant_id")),
                "genes": genes,
                "selection_score": round(priority_score, 3),
                "selection_reasons": reasons,
                "query_payload": payload,
                "query": query_text,
                "rounds": rounds,
            }
        )
    return results


def flatten_variant_rag_hits(result: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    hits_by_variant: defaultdict[str, list[dict[str, Any]]] = defaultdict(list)
    for variant_result in result.get("variant_rag_results", []) or []:
        variant_id = clean_text(variant_result.get("variant_id"))
        genes = set(clean_text(gene).upper() for gene in variant_result.get("genes", []) if clean_text(gene))
        for round_result in variant_result.get("rounds", []) or []:
            for hit in round_result.get("hits", []) or []:
                haystack = " ".join(
                    [
                        clean_text(hit.get("title")),
                        clean_text(hit.get("text")),
                        json.dumps(hit.get("metadata", {}), ensure_ascii=False),
                    ]
                ).upper()
                reasons = ["variant_specific_rag_hit"]
                if genes and any(re.search(rf"\b{re.escape(gene)}\b", haystack) for gene in genes):
                    reasons.append("gene_mentioned_by_rag")
                hits_by_variant[variant_id].append(
                    compact_rag_hit(
                        round_result,
                        hit,
                        ",".join(sorted(set(reasons))),
                        scope="variant",
                        variant_id=variant_id,
                    )
                )
    return hits_by_variant


def flatten_variant_rule_hits(result: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    rules_by_variant: defaultdict[str, list[dict[str, Any]]] = defaultdict(list)
    for variant_result in result.get("variant_rag_results", []) or []:
        variant_id = clean_text(variant_result.get("variant_id"))
        seen: set[str] = set()
        for round_result in variant_result.get("rounds", []) or []:
            for row in round_result.get("rule_hits", []) or []:
                key = clean_text(row.get("rule_id")) or f"{row.get('source_id')}:{row.get('row_num')}"
                if key and key in seen:
                    continue
                seen.add(key)
                compact = dict(row)
                compact["round_id"] = clean_text(round_result.get("round_id"))
                rules_by_variant[variant_id].append(compact)
    return rules_by_variant


def parse_json_list(value: Any) -> list[dict[str, Any]]:
    text = clean_text(value)
    if not text:
        return []
    try:
        parsed = json.loads(text)
    except json.JSONDecodeError:
        return []
    if not isinstance(parsed, list):
        return []
    return [item for item in parsed if isinstance(item, dict)]


def build_graphrag_path(
    row: dict[str, Any],
    rag_hits: list[dict[str, Any]],
    rule_hits: list[dict[str, Any]],
    missing_context: list[str],
) -> dict[str, Any]:
    variant_id = clean_text(row.get("variant_id"))
    genes = split_semicolon(row.get("genes"))
    nodes: list[dict[str, Any]] = [
        {"id": f"Variant:{variant_id}", "label": "Variant", "name": variant_id, "properties": {"type": row.get("type"), "tumor_af": row.get("tumor_af")}},
        {"id": f"CallerQC:{variant_id}", "label": "CallerQC", "name": clean_text(row.get("qc_tier")), "properties": {"score": row.get("technical_score"), "flags": row.get("qc_flags")}},
        {"id": f"SourceTrace:{variant_id}", "label": "SourceTraceStatus", "name": clean_text(row.get("source_trace_status")), "properties": {"action": row.get("recommended_review_action")}},
    ]
    edges: list[dict[str, Any]] = [
        {"source": f"Variant:{variant_id}", "target": f"CallerQC:{variant_id}", "type": "HAS_CALLER_QC"},
        {"source": f"Variant:{variant_id}", "target": f"SourceTrace:{variant_id}", "type": "HAS_SOURCE_TRACE_STATUS"},
    ]
    for gene in genes:
        nodes.append({"id": f"Gene:{gene}", "label": "Gene", "name": gene, "properties": {}})
        edges.append({"source": f"Variant:{variant_id}", "target": f"Gene:{gene}", "type": "HAS_GENE_ANNOTATION"})
    kg_records = parse_json_list(row.get("kg_evidence_details"))
    for rec in kg_records[:5]:
        evidence_id = clean_text(rec.get("evidence_id")) or clean_text(rec.get("title"))[:40]
        node_id = f"KGEvidence:{evidence_id}"
        nodes.append(
            {
                "id": node_id,
                "label": "KGEvidence",
                "name": evidence_id,
                "properties": {
                    "confidence_score": rec.get("confidence_score"),
                    "graph_inclusion_status": rec.get("graph_inclusion_status"),
                    "doi": rec.get("doi"),
                    "pmid": rec.get("pmid"),
                    "title": rec.get("title"),
                },
            }
        )
        edges.append({"source": f"Variant:{variant_id}", "target": node_id, "type": "MATCHED_LOCAL_KG_EVIDENCE"})
    for hit in rag_hits[:5]:
        evidence_id = clean_text(hit.get("evidence_id")) or clean_text(hit.get("chunk_id"))
        node_id = f"RAGEvidence:{evidence_id}"
        nodes.append(
            {
                "id": node_id,
                "label": "RAGEvidence",
                "name": evidence_id,
                "properties": {
                    "scope": hit.get("scope"),
                    "round_id": hit.get("round_id"),
                    "categories": hit.get("categories"),
                    "confidence_score": hit.get("confidence_score"),
                    "doi": hit.get("doi"),
                    "pmid": hit.get("pmid"),
                    "title": hit.get("title"),
                },
            }
        )
        edges.append({"source": node_id, "target": f"SourceTrace:{variant_id}", "type": "SUPPORTS_SOURCE_TRACE_CONTEXT"})
    for rule in rule_hits[:4]:
        rule_id = clean_text(rule.get("rule_id")) or f"{rule.get('source_id')}:{rule.get('row_num')}"
        node_id = f"Rule:{rule_id}"
        nodes.append(
            {
                "id": node_id,
                "label": "Rule",
                "name": rule_id,
                "properties": {
                    "rule_type": rule.get("rule_type"),
                    "rule_status": rule.get("rule_status"),
                    "cancer_type": rule.get("cancer_type"),
                    "supporting_dois": rule.get("supporting_dois"),
                },
            }
        )
        target = f"Gene:{genes[0]}" if genes else f"Variant:{variant_id}"
        edges.append({"source": node_id, "target": target, "type": "STRUCTURED_RULE_CONTEXT"})
    for item in missing_context[:6]:
        node_id = f"MissingContext:{variant_id}:{item}"
        nodes.append({"id": node_id, "label": "MissingContext", "name": item, "properties": {}})
        edges.append({"source": node_id, "target": f"SourceTrace:{variant_id}", "type": "BLOCKS_FINAL_SOURCE_ASSIGNMENT"})

    category_counts: Counter[str] = Counter()
    for hit in rag_hits:
        for category in hit.get("categories", []) if isinstance(hit.get("categories"), list) else []:
            category_counts[category] += 1
    path_parts = [
        f"Variant {variant_id}",
        f"Gene {','.join(genes) if genes else 'not_annotated'}",
        f"KG evidence {len(kg_records)}",
        f"RAG evidence {len(rag_hits)}",
        f"Rule context {len(rule_hits)}",
        f"Missing context {len(missing_context)}",
    ]
    if category_counts:
        path_parts.append("RAG categories " + counter_line(dict(category_counts)))
    return {
        "summary": " -> ".join(path_parts),
        "nodes": nodes,
        "edges": edges,
    }


def source_explanation_for_status(row: dict[str, Any]) -> str:
    variant_id = clean_text(row.get("variant_id"))
    genes = clean_text(row.get("genes")) or "gene_not_annotated"
    status = clean_text(row.get("source_trace_status"))
    if status == "traceable_caller_evidence_candidate":
        return (
            f"{variant_id} ({genes}) is a caller-derived molecular evidence candidate with local KG provenance. "
            "Tumor-informed MRD or clonal recurrence remains plausible, but final source assignment needs tumor concordance, WBC status, assay LOD, and serial behavior."
        )
    if status == "traceable_caller_signal_needs_serial_confirmation":
        return (
            f"{variant_id} ({genes}) is traceable at caller/gene level but remains source-limited. "
            "Serial or orthogonal confirmation is needed to separate persistent tumor-derived signal from transient low-support, hematopoietic, germline, or technical origin."
        )
    if status == "hematopoietic_or_CHIP_origin_suspected":
        return (
            f"{variant_id} ({genes}) should first be audited as possible CHIP/hematopoietic-origin signal. "
            "Matched WBC or buffy-coat comparison is required before assigning it to tumor-derived MRD."
        )
    if status == "technical_artifact_origin_suspected":
        return (
            f"{variant_id} ({genes}) is artifact-prone in the caller layer. "
            "Local alignment, strand/context, depth, and caller rescue are prerequisites before biological source assignment."
        )
    if status == "technical_rescue_required_before_source_assignment":
        return (
            f"{variant_id} ({genes}) requires technical rescue before source assignment. "
            "Current evidence should be treated as provisional caller-derived signal, not as a biological rejection."
        )
    if status == "cross_sectional_signal_needs_contextualization":
        return (
            f"{variant_id} ({genes}) is a cross-sectional caller signal with unresolved source. "
            "Tumor/WBC/serial metadata are needed before assigning tumor, hematopoietic, transient, germline, or artifact origin."
        )
    return (
        f"{variant_id} ({genes}) cannot be assigned to a biological source with the current context. "
        "Additional annotation and source-tracing metadata are required."
    )


def next_step_for_row(row: dict[str, Any]) -> str:
    status = clean_text(row.get("source_trace_status"))
    reasons = clean_text(row.get("source_trace_reason_codes"))
    action = clean_text(row.get("recommended_review_action"))
    if status == "technical_rescue_required_before_source_assignment" or "technical_artifact_features_present" in reasons:
        return "Perform manual caller rescue and local alignment/artifact review, then repeat source tracing."
    if status == "hematopoietic_or_CHIP_origin_suspected" or "chip_gene_overlap" in reasons:
        return "Review matched WBC/buffy coat or germline/hematopoietic filtering for the same variant."
    if "low_variant_support" in reasons:
        return "Confirm with repeat/serial plasma sampling or orthogonal assay, using assay LOD as the interpretation boundary."
    if clean_text(row.get("kg_evidence_ids")):
        return "Compare against primary tumor markers and follow serial ctDNA trajectory for tumor-informed source assignment."
    return action.replace("_", " ") if action else "Collect tumor tissue, WBC, serial sampling, and assay metadata before source assignment."


def synthesize_variant_source_trace(
    row: dict[str, Any],
    rag_hits: list[dict[str, Any]],
    rule_hits: list[dict[str, Any]],
    checklist_missing: list[str],
) -> dict[str, Any]:
    kg_ids = split_semicolon(row.get("kg_evidence_ids"))
    missing_context = split_semicolon(row.get("caller_required_context_to_complete"))
    for item in checklist_missing:
        field = item.split(":", 1)[0]
        if field and field not in missing_context:
            missing_context.append(field)
    categories: Counter[str] = Counter()
    for hit in rag_hits:
        for category in hit.get("categories", []) if isinstance(hit.get("categories"), list) else []:
            categories[category] += 1
    direct_rag_hits = [
        hit
        for hit in rag_hits
        if clean_text(hit.get("scope")) == "variant" or "gene_mentioned_by_rag" in clean_text(hit.get("match_reason"))
    ]
    support_bits = []
    if kg_ids:
        support_bits.append(f"local KG evidence: {', '.join(kg_ids[:5])}")
    if direct_rag_hits:
        support_bits.append(f"direct/variant RAG evidence categories: {counter_line(dict(categories)) if categories else 'general_source_trace'}")
    elif rag_hits:
        support_bits.append(f"scenario-level RAG context categories: {counter_line(dict(categories)) if categories else 'general_source_trace'}")
    if rule_hits:
        support_bits.append(f"structured rule context: {', '.join(clean_text(rule.get('rule_id')) for rule in rule_hits[:3] if clean_text(rule.get('rule_id')))}")
    if not support_bits:
        support_bits.append("caller QC and scenario context only; no direct KG/RAG evidence attached")

    contradiction_bits = []
    status = clean_text(row.get("source_trace_status"))
    if status in {"traceable_caller_evidence_candidate", "traceable_caller_signal_needs_serial_confirmation"}:
        if categories.get("chip_wbc"):
            contradiction_bits.append("CHIP/WBC-origin evidence remains an alternative explanation")
        if categories.get("artifact_or_confounder"):
            contradiction_bits.append("artifact/confounder evidence remains an alternative explanation")
    if not clean_text(row.get("kg_evidence_ids")):
        contradiction_bits.append("no direct local KG evidence matched this exact gene/coordinate")
    if "matched_control_alt_signal_observed" in clean_text(row.get("source_trace_reason_codes")):
        contradiction_bits.append("matched/control alternate signal requires hematopoietic or technical review")
    if not contradiction_bits:
        contradiction_bits.append("no explicit contradictory source evidence was mapped; missing context still limits assignment")

    if status == "technical_rescue_required_before_source_assignment":
        confidence = "requires_technical_rescue"
    elif status == "hematopoietic_or_CHIP_origin_suspected":
        confidence = "requires_wbc_or_buffy_coat_review"
    elif kg_ids and direct_rag_hits:
        confidence = "provisional_kg_rag_supported"
    elif kg_ids:
        confidence = "provisional_kg_supported_context_limited"
    elif direct_rag_hits:
        confidence = "provisional_rag_supported_context_limited"
    elif rag_hits:
        confidence = "context_limited_with_scenario_rag"
    else:
        confidence = "context_limited"

    synthesis = {
        "variant_id": clean_text(row.get("variant_id")),
        "source_trace_category": status or "unresolved_missing_context",
        "confidence_level": confidence,
        "source_explanation": source_explanation_for_status(row),
        "supporting_evidence_summary": "; ".join(support_bits),
        "contradictory_evidence_summary": "; ".join(contradiction_bits),
        "kg_path_summary": "",
        "missing_context": missing_context[:12],
        "next_confirmation_step": next_step_for_row(row),
        "rag_hit_count": len(rag_hits),
        "direct_rag_hit_count": len(direct_rag_hits),
        "rule_hit_count": len(rule_hits),
    }
    return synthesis


def build_case_source_trace_synthesis(
    summary: dict[str, Any],
    rows: list[dict[str, Any]],
    result: dict[str, Any],
) -> dict[str, Any]:
    status_counts = Counter(clean_text(row.get("source_trace_status")) or "unresolved_missing_context" for row in rows)
    confidence_counts = Counter(clean_text(row.get("source_assignment_confidence")) or "not_synthesized" for row in rows)
    variant_results = result.get("variant_rag_results", []) or []
    missing_context_items = result.get("integrated_source_trace_summary", {}).get("missing_context_items", [])
    dominant_status = status_counts.most_common(1)[0][0] if status_counts else "unresolved_missing_context"
    if dominant_status == "technical_rescue_required_before_source_assignment":
        dominant = "technical rescue is the dominant prerequisite before biological source assignment"
    elif dominant_status == "hematopoietic_or_CHIP_origin_suspected":
        dominant = "CHIP/hematopoietic origin audit is the dominant source-tracing concern"
    elif dominant_status in {"traceable_caller_evidence_candidate", "traceable_caller_signal_needs_serial_confirmation"}:
        dominant = "traceable caller-derived molecular evidence is present, but final source assignment remains context-limited"
    else:
        dominant = "source assignment remains unresolved because context and direct evidence are incomplete"
    selected = []
    for row in rows:
        if clean_text(row.get("source_trace_synthesis")):
            selected.append(
                {
                    "variant_id": row.get("variant_id"),
                    "genes": row.get("genes"),
                    "source_trace_category": row.get("source_trace_status"),
                    "confidence_level": row.get("source_assignment_confidence"),
                    "source_explanation": row.get("source_trace_explanation"),
                    "next_confirmation_step": row.get("next_confirmation_step"),
                }
            )
        if len(selected) >= 24:
            break
    return {
        "case_source_trace_answer": summary.get("molecular_evidence_qc", {}).get("nlp_qc_conclusion", ""),
        "dominant_source_explanation": dominant,
        "workflow_readiness": summary.get("mrd_clinical_interpretation_entry_gate", {}).get("gate_class", ""),
        "evidence_quality": {
            "case_round_count": len(result.get("rounds", []) or []),
            "variant_specific_rag_queries": len(variant_results),
            "variant_rows_with_synthesis": sum(1 for row in rows if clean_text(row.get("source_trace_synthesis"))),
            "source_trace_status_counts": dict(status_counts),
            "source_assignment_confidence_counts": dict(confidence_counts),
        },
        "unresolved_context": missing_context_items,
        "caller_layer_completed": [
            "variant identity",
            "gene interval annotation when available",
            "VAF/depth/caller QC",
            "local KG provenance lookup",
            "RAG evidence mapping",
            "GraphRAG-style path summary",
        ],
        "requires_external_metadata": [
            "matched WBC or buffy coat",
            "primary tumor marker set or tumor tissue concordance",
            "serial ctDNA samples",
            "assay LOD and panel design",
            "clinical timepoint and treatment setting",
        ],
        "selected_variant_summaries": selected,
    }


def integrate_source_trace_rag_into_rows(
    rows: list[dict[str, Any]],
    result: dict[str, Any],
    max_hits_per_variant: int = 5,
    complete_mode: bool = False,
) -> dict[str, Any]:
    """Attach multi-agent RAG evidence summaries back to variant-level QC rows."""
    round_hits: list[tuple[dict[str, Any], dict[str, Any]]] = []
    round_counts: Counter[str] = Counter()
    category_counts: Counter[str] = Counter()
    scenario_counts: Counter[str] = Counter()
    for round_result in result.get("rounds", []):
        round_id = clean_text(round_result.get("round_id"))
        for hit in round_result.get("hits", []) or []:
            round_hits.append((round_result, hit))
            round_counts[round_id] += 1
            for category in classify_rag_hit(hit):
                category_counts[category] += 1
            for scenario in RAG_ROUND_SCENARIO_MAP.get(round_id, []):
                scenario_counts[scenario] += 1
        rule_hits = round_result.get("rule_hits", []) or []
        for scenario in RAG_ROUND_SCENARIO_MAP.get(round_id, []):
            scenario_counts[scenario] += len(rule_hits)

    checklist_missing = [
        f"{item.get('field')}: {item.get('reason')}"
        for item in result.get("missing_information_checklist", [])
        if clean_text(item.get("status")) == "missing"
    ]
    agent_rounds = [
        f"{clean_text(round_result.get('round_id'))}/{clean_text(round_result.get('agent'))}"
        for round_result in result.get("rounds", [])
    ]
    agent_summary = (
        f"RAG rounds integrated: {', '.join(agent_rounds)}."
        if agent_rounds
        else "RAG was requested but no retrieval rounds were available."
    )

    variant_specific_hits = flatten_variant_rag_hits(result)
    variant_specific_rules = flatten_variant_rule_hits(result)

    integrated_rows = 0
    direct_gene_rows = 0
    direct_or_variant_rows = 0
    synthesized_rows = 0
    kg_path_rows = 0
    for row in rows:
        variant_id = clean_text(row.get("variant_id"))
        genes = [gene.upper() for gene in split_semicolon(row.get("genes"))]
        scenarios = split_semicolon(row.get("caller_dependent_mrd_scenarios"))
        row_status = clean_text(row.get("source_trace_status"))
        reasons = clean_text(row.get("source_trace_reason_codes"))
        matched_hits: list[dict[str, Any]] = list(variant_specific_hits.get(variant_id, []))

        for round_result, hit in round_hits:
            haystack = " ".join(
                [
                    clean_text(hit.get("title")),
                    clean_text(hit.get("text")),
                    json.dumps(hit.get("metadata", {}), ensure_ascii=False),
                ]
            ).upper()
            round_id = clean_text(round_result.get("round_id"))
            mapped_scenarios = set(RAG_ROUND_SCENARIO_MAP.get(round_id, []))
            direct_gene_match = bool(genes and any(re.search(rf"\b{re.escape(gene)}\b", haystack) for gene in genes))
            scenario_match = bool(mapped_scenarios & set(scenarios))
            status_match = any(token in haystack.lower() for token in row_status.split("_") if len(token) > 5)
            chip_context_match = "chip_gene_overlap" in reasons and any(
                token in haystack.lower()
                for token in ["chip", "clonal hematopo", "wbc", "buffy", "hematopoietic", "germline"]
            )
            serial_context_match = "low_variant_support" in reasons and any(
                token in haystack.lower()
                for token in ["serial", "longitudinal", "repeat", "trajectory", "persistence", "clearance"]
            )
            artifact_context_match = "technical_artifact_features_present" in reasons and any(
                token in haystack.lower()
                for token in ["artifact", "artefact", "alignment", "confound", "background"]
            )

            match_reasons: list[str] = []
            if direct_gene_match:
                match_reasons.append("gene_mentioned_by_rag")
            if scenario_match:
                match_reasons.append("caller_scenario_round_match")
            if status_match:
                match_reasons.append("source_trace_status_context_match")
            if chip_context_match:
                match_reasons.append("chip_wbc_context_match")
            if serial_context_match:
                match_reasons.append("serial_confirmation_context_match")
            if artifact_context_match:
                match_reasons.append("artifact_context_match")

            if match_reasons:
                matched_hits.append(
                    compact_rag_hit(
                        round_result,
                        hit,
                        ",".join(sorted(set(match_reasons))),
                        scope="case",
                        variant_id=variant_id,
                    )
                )

        deduped: list[dict[str, Any]] = []
        seen_keys: set[str] = set()
        for hit in sorted(matched_hits, key=lambda item: float(item.get("score") or 0), reverse=True):
            key = hit.get("chunk_id") or f"{hit.get('round_id')}::{hit.get('title')}"
            if key in seen_keys:
                continue
            seen_keys.add(key)
            deduped.append(hit)
            if len(deduped) >= max_hits_per_variant:
                break

        direct_gene_hit_count = sum(1 for hit in deduped if "gene_mentioned_by_rag" in clean_text(hit.get("match_reason")))
        variant_specific_hit_count = sum(1 for hit in deduped if clean_text(hit.get("scope")) == "variant")
        direct_or_variant_hit_count = direct_gene_hit_count + variant_specific_hit_count
        if deduped:
            integrated_rows += 1
        if direct_gene_hit_count:
            direct_gene_rows += 1
        if direct_or_variant_hit_count:
            direct_or_variant_rows += 1

        missing_context = split_semicolon(row.get("caller_required_context_to_complete"))
        uncertainty = []
        if missing_context:
            uncertainty.append("needs: " + ", ".join(missing_context[:8]))
        if checklist_missing:
            uncertainty.append("case missing: " + "; ".join(short_text(item, 120) for item in checklist_missing[:4]))
        if not deduped:
            uncertainty.append("no direct RAG hit mapped to this variant; use KG provenance and scenario-level RAG context")

        scenario_labels = [
            scenario
            for scenario in scenarios
            if scenario_counts.get(scenario, 0) > 0
        ]
        if direct_or_variant_hit_count:
            row["rag_evidence_status"] = "rag_variant_or_direct_context_integrated"
        elif deduped:
            row["rag_evidence_status"] = "rag_scenario_context_only"
        else:
            row["rag_evidence_status"] = "rag_scenario_context_only"
        row["rag_agent_summary"] = agent_summary
        row["rag_support_summary"] = (
            f"RAG mapped {len(deduped)} evidence hit(s) to this variant; "
            f"direct/variant hits={direct_or_variant_hit_count}; scenario support={', '.join(scenario_labels) if scenario_labels else 'general_source_trace'}."
        )
        row["rag_supporting_evidence"] = "; ".join(
            f"{hit.get('round_id')}:{short_text(hit.get('title'), 110)}"
            for hit in deduped[:3]
        )
        row["rag_uncertainty_checklist"] = " | ".join(uncertainty)
        row["rag_evidence_details"] = json.dumps(deduped, ensure_ascii=False)

        row_rule_hits = variant_specific_rules.get(variant_id, [])
        missing_context = split_semicolon(row.get("caller_required_context_to_complete"))
        path = build_graphrag_path(row, deduped, row_rule_hits, missing_context)
        row["graphrag_path_summary"] = path["summary"]
        row["graphrag_path_nodes"] = json.dumps(path["nodes"], ensure_ascii=False)
        row["graphrag_path_edges"] = json.dumps(path["edges"], ensure_ascii=False)
        if path["nodes"]:
            kg_path_rows += 1

        if complete_mode or deduped or row_rule_hits or clean_text(row.get("kg_evidence_ids")):
            synthesis = synthesize_variant_source_trace(row, deduped, row_rule_hits, checklist_missing)
            synthesis["kg_path_summary"] = path["summary"]
            row["source_trace_synthesis"] = json.dumps(synthesis, ensure_ascii=False)
            row["source_trace_explanation"] = synthesis["source_explanation"]
            row["supporting_evidence_summary"] = synthesis["supporting_evidence_summary"]
            row["contradictory_evidence_summary"] = synthesis["contradictory_evidence_summary"]
            row["missing_context_summary"] = ";".join(synthesis["missing_context"])
            row["next_confirmation_step"] = synthesis["next_confirmation_step"]
            row["source_assignment_confidence"] = synthesis["confidence_level"]
            synthesized_rows += 1

    summary = {
        "integrated_variant_rows": integrated_rows,
        "direct_gene_matched_variant_rows": direct_gene_rows,
        "direct_or_variant_context_rows": direct_or_variant_rows,
        "variant_specific_rag_rows": len(result.get("variant_rag_results", []) or []),
        "variant_synthesis_rows": synthesized_rows,
        "graphrag_path_rows": kg_path_rows,
        "round_hit_counts": dict(round_counts),
        "category_hit_counts": dict(category_counts),
        "scenario_hit_counts": dict(scenario_counts),
        "missing_context_items": checklist_missing,
        "caller_workflow_coverage": build_caller_workflow_coverage(rows, result),
    }
    result["integrated_source_trace_summary"] = summary
    return summary


def run_source_trace_rag(
    summary: dict[str, Any],
    rows: list[dict[str, Any]],
    args: argparse.Namespace,
    rag_json: Path,
    rag_md: Path,
) -> dict[str, Any]:
    from mrdskill_query_mrd_kg_rag import existing_chunks_path, search_rule_db
    from mrdskill_run_adjudication_rag import checklist, retrieve_rounds_streaming

    payload = build_source_trace_rag_payload(summary, rows)
    base_text = payload_text(payload)
    chunks_path = existing_chunks_path(args.rag_chunks_jsonl)
    max_chunks = 0 if args.rag_full_scan else args.rag_max_chunks
    all_genes = sorted(collect_genes_from_rows(rows))
    rounds = retrieve_rounds_streaming(
        chunks_path,
        base_text,
        args.rag_top_k_per_round,
        max_chunks,
        boost_terms=all_genes[:50],
        min_score=args.rag_min_score,
        denoise=args.rag_denoise,
    )
    for round_result in rounds:
        round_result["rule_hits"] = search_rule_db(
            artifact_path(args.rag_rule_db),
            round_result["query"],
            args.rag_rule_top_k_per_round,
            allowed_genes=all_genes,
            require_allowed_gene=bool(all_genes),
        )
    variant_rag_results = build_variant_rag_results(rows, summary, args, chunks_path)

    result = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "chunks_jsonl": str(chunks_path),
        "rule_db": str(artifact_path(args.rag_rule_db)),
        "rag_mode": args.rag_mode,
        "rag_profile": {
            "top_k_per_round": args.rag_top_k_per_round,
            "rule_top_k_per_round": args.rag_rule_top_k_per_round,
            "max_chunks": max_chunks,
            "full_scan": bool(args.rag_full_scan or args.rag_max_chunks == 0),
            "variant_limit": args.rag_variant_limit,
            "min_score": args.rag_min_score,
            "denoise": args.rag_denoise,
            "max_hits_per_variant": args.rag_max_hits_per_variant,
        },
        "case": payload,
        "case_text": base_text,
        "missing_information_checklist": checklist(payload),
        "rounds": rounds,
        "variant_rag_results": variant_rag_results,
    }
    integration_summary = integrate_source_trace_rag_into_rows(
        rows,
        result,
        max_hits_per_variant=args.rag_max_hits_per_variant,
        complete_mode=args.rag_mode == "complete",
    )
    result["case_source_trace_synthesis"] = build_case_source_trace_synthesis(summary, rows, result)
    mqc = summary.get("molecular_evidence_qc", {})
    if mqc:
        top_ids = [clean_text(item.get("variant_id")) for item in mqc.get("top_variant_conclusions", [])]
        row_by_id = {clean_text(row.get("variant_id")): row for row in rows}
        if top_ids:
            mqc["top_variant_conclusions"] = [
                make_variant_brief(row_by_id[variant_id])
                for variant_id in top_ids
                if variant_id in row_by_id
            ]
        mqc["caller_workflow_coverage"] = integration_summary["caller_workflow_coverage"]
        mqc["rag_integration_summary"] = {
            key: value
            for key, value in integration_summary.items()
            if key != "caller_workflow_coverage"
        }
        mqc["case_source_trace_synthesis"] = result["case_source_trace_synthesis"]
    rag_json.parent.mkdir(parents=True, exist_ok=True)
    rag_json.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    write_source_trace_rag_markdown(rag_md, result)
    return {
        "enabled": True,
        "json": str(rag_json),
        "md": str(rag_md),
        "round_count": len(rounds),
        "rag_mode": args.rag_mode,
        "integrated_variant_rows": integration_summary.get("integrated_variant_rows", 0),
        "direct_gene_matched_variant_rows": integration_summary.get("direct_gene_matched_variant_rows", 0),
        "direct_or_variant_context_rows": integration_summary.get("direct_or_variant_context_rows", 0),
        "variant_specific_rag_rows": integration_summary.get("variant_specific_rag_rows", 0),
        "variant_synthesis_rows": integration_summary.get("variant_synthesis_rows", 0),
        "graphrag_path_rows": integration_summary.get("graphrag_path_rows", 0),
        "scenario_hit_counts": integration_summary.get("scenario_hit_counts", {}),
        "chunks_scanned": rounds[0].get("chunks_scanned", 0) if rounds else 0,
        "top_k_per_round": args.rag_top_k_per_round,
        "rule_top_k_per_round": args.rag_rule_top_k_per_round,
    }


MANIFEST_TEMPLATE_COLUMNS = [
    "case_id",
    "sample_id",
    "sample_type",
    "timepoint",
    "days_from_event",
    "event",
    "vcf_path",
    "assay_type",
    "assay_lod",
    "matched_tumor_id",
    "matched_wbc_id",
    "collection_date",
    "notes",
]


LONGITUDINAL_VARIANT_FIELDS = [
    "case_id",
    "canonical_variant_id",
    "representative_variant_id",
    "genes",
    "variant_type",
    "source_assignment_status",
    "source_assignment_confidence",
    "technical_readiness",
    "tumor_concordance_status",
    "wbc_concordance_status",
    "plasma_series_status",
    "plasma_af_trend",
    "sample_observation_count",
    "tumor_observation_count",
    "wbc_observation_count",
    "plasma_observation_count",
    "first_plasma_timepoint",
    "last_plasma_timepoint",
    "first_plasma_af",
    "last_plasma_af",
    "min_plasma_af",
    "max_plasma_af",
    "near_lod_observation_count",
    "assay_lod_values",
    "chip_review_genes",
    "kg_match_statuses",
    "kg_evidence_ids",
    "qc_tier_counts",
    "source_trace_status_counts",
    "caller_scenarios",
    "recommended_next_step",
    "observed_samples",
    "missing_or_not_reported_samples",
    "observation_trace",
]


LONGITUDINAL_OBSERVATION_FIELDS = [
    "case_id",
    "manifest_sample_id",
    "manifest_sample_type",
    "sample_type_group",
    "timepoint",
    "days_from_event",
    "event",
    "assay_type",
    "assay_lod",
    "matched_tumor_id",
    "matched_wbc_id",
    "source_vcf",
    "canonical_variant_id",
] + [
    "variant_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "type",
    "filter",
    "status",
    "tumor_dp",
    "tumor_vd",
    "tumor_af",
    "control_dp",
    "control_vd",
    "control_af",
    "genes",
    "kg_match_status",
    "kg_evidence_ids",
    "technical_score",
    "qc_tier",
    "source_trace_status",
    "evidence_mode",
    "caller_dependent_mrd_scenarios",
    "recommended_review_action",
    "source_trace_reason_codes",
    "qc_flags",
    "review_reasons",
]


def write_case_manifest_template(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_TEMPLATE_COLUMNS)
        writer.writeheader()
        writer.writerow(
            {
                "case_id": "CASE001",
                "sample_id": "CASE001_TUMOR",
                "sample_type": "tumor_tissue",
                "timepoint": "baseline",
                "days_from_event": "0",
                "event": "surgery_or_baseline",
                "vcf_path": "inputs/CASE001_tumor.vcf",
                "assay_type": "tumor_ngs",
                "assay_lod": "",
                "matched_tumor_id": "",
                "matched_wbc_id": "CASE001_WBC",
                "collection_date": "",
                "notes": "Primary tumor or baseline tissue marker set.",
            }
        )
        writer.writerow(
            {
                "case_id": "CASE001",
                "sample_id": "CASE001_WBC",
                "sample_type": "matched_wbc",
                "timepoint": "baseline",
                "days_from_event": "0",
                "event": "germline_or_chip_control",
                "vcf_path": "inputs/CASE001_wbc.vcf",
                "assay_type": "matched_control_ngs",
                "assay_lod": "",
                "matched_tumor_id": "CASE001_TUMOR",
                "matched_wbc_id": "",
                "collection_date": "",
                "notes": "Matched WBC/buffy coat for CHIP, hematopoietic, and germline review.",
            }
        )
        writer.writerow(
            {
                "case_id": "CASE001",
                "sample_id": "CASE001_P01",
                "sample_type": "plasma_ctdna",
                "timepoint": "postop_4w",
                "days_from_event": "28",
                "event": "postoperative_surveillance",
                "vcf_path": "inputs/CASE001_plasma_t1.vcf",
                "assay_type": "ctdna_vcf",
                "assay_lod": "0.001",
                "matched_tumor_id": "CASE001_TUMOR",
                "matched_wbc_id": "CASE001_WBC",
                "collection_date": "",
                "notes": "Serial plasma/cfDNA timepoint.",
            }
        )


def read_case_manifest(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        raise FileNotFoundError(f"Case manifest not found: {path}")
    delimiter = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        rows = []
        for index, row in enumerate(reader, start=1):
            cleaned = {clean_text(key): clean_text(value) for key, value in row.items() if key is not None}
            if not any(cleaned.values()):
                continue
            cleaned["_manifest_index"] = index
            rows.append(cleaned)
    if not rows:
        raise ValueError(f"Case manifest has no data rows: {path}")
    missing = [column for column in ["sample_id", "vcf_path"] if not clean_text(rows[0].get(column))]
    if missing:
        raise ValueError(f"Case manifest requires columns with values: {', '.join(missing)}")
    return rows


def resolve_manifest_vcf_path(raw_path: str, manifest_path: Path) -> Path:
    path = Path(raw_path)
    if path.is_absolute():
        return path
    relative = manifest_path.parent / path
    if relative.exists():
        return relative
    return artifact_path(path)


def sample_type_group(sample_type: str) -> str:
    text = clean_text(sample_type).lower()
    if any(term in text for term in ["tumor", "tumour", "tissue", "primary", "metastasis", "metastatic"]):
        return "tumor"
    if any(term in text for term in ["wbc", "buffy", "leukocyte", "leucocyte", "normal", "germline", "hematopoietic", "haematopoietic"]):
        return "wbc"
    if any(term in text for term in ["plasma", "cfdna", "ctdna", "liquid", "serum"]):
        return "plasma"
    return "other"


def canonical_variant_id_from_row(row: dict[str, Any]) -> str:
    return f"{normalize_chrom(clean_text(row.get('chrom')))}:{row.get('pos')}:{clean_text(row.get('ref'))}>{clean_text(row.get('alt'))}"


def numeric_or_none(value: Any) -> float | None:
    num = to_float(value)
    if num is None or math.isinf(num):
        return None
    return num


def sort_observations(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    def key(row: dict[str, Any]) -> tuple[int, float, str, str]:
        day = numeric_or_none(row.get("days_from_event"))
        if day is None:
            return (1, 0.0, clean_text(row.get("timepoint")), clean_text(row.get("manifest_sample_id")))
        return (0, day, clean_text(row.get("timepoint")), clean_text(row.get("manifest_sample_id")))

    return sorted(rows, key=key)


def summarize_observation_trace(rows: list[dict[str, Any]]) -> str:
    parts = []
    for row in sort_observations(rows):
        sample = clean_text(row.get("manifest_sample_id"))
        sample_type = clean_text(row.get("manifest_sample_type"))
        timepoint = clean_text(row.get("timepoint")) or clean_text(row.get("days_from_event"))
        af = row.get("tumor_af")
        vd = row.get("tumor_vd")
        tier = clean_text(row.get("qc_tier"))
        parts.append(f"{sample}:{sample_type}:{timepoint}:AF={af if af is not None else ''}:VD={vd if vd is not None else ''}:QC={tier}")
    return " | ".join(parts)


def sample_label(row: dict[str, Any]) -> str:
    bits = [clean_text(row.get("manifest_sample_id"))]
    if clean_text(row.get("manifest_sample_type")):
        bits.append(clean_text(row.get("manifest_sample_type")))
    if clean_text(row.get("timepoint")):
        bits.append(clean_text(row.get("timepoint")))
    return ":".join(bits)


def collect_split_values(rows: list[dict[str, Any]], key: str) -> list[str]:
    values: set[str] = set()
    for row in rows:
        for token in re.split(r"[;,|]+", clean_text(row.get(key))):
            token = clean_text(token)
            if token:
                values.add(token)
    return sorted(values)


def assay_lod_status(rows: list[dict[str, Any]]) -> tuple[int, list[str]]:
    near_lod = 0
    lod_values: set[str] = set()
    for row in rows:
        lod = numeric_or_none(row.get("assay_lod"))
        af = numeric_or_none(row.get("tumor_af"))
        if lod is None:
            continue
        lod_values.add(f"{lod:g}")
        if af is not None and af <= lod * 2:
            near_lod += 1
    return near_lod, sorted(lod_values)


def plasma_trend(plasma_rows: list[dict[str, Any]]) -> tuple[str, str, dict[str, Any]]:
    if not plasma_rows:
        return "no_plasma_signal_observed_or_reported", "not_assessable", {}
    ordered = sort_observations(plasma_rows)
    af_rows = [row for row in ordered if numeric_or_none(row.get("tumor_af")) is not None]
    if len(ordered) == 1:
        row = ordered[0]
        return "single_plasma_timepoint_observed", "single_timepoint", {
            "first": row,
            "last": row,
            "values": [numeric_or_none(row.get("tumor_af"))],
        }
    if len(af_rows) < 2:
        return "serial_plasma_observed_vaf_not_assessable", "not_assessable", {
            "first": ordered[0],
            "last": ordered[-1],
            "values": [numeric_or_none(row.get("tumor_af")) for row in ordered],
        }
    first = af_rows[0]
    last = af_rows[-1]
    first_af = numeric_or_none(first.get("tumor_af")) or 0.0
    last_af = numeric_or_none(last.get("tumor_af")) or 0.0
    delta_guard = max(0.001, abs(first_af) * 0.2)
    if last_af > first_af + delta_guard:
        trend = "rising"
    elif last_af < first_af - delta_guard:
        trend = "falling"
    else:
        trend = "persistent_or_stable"
    values = [numeric_or_none(row.get("tumor_af")) for row in af_rows]
    return "serial_plasma_signal_observed", trend, {"first": first, "last": last, "values": values}


def build_longitudinal_variant_rows(
    observation_rows: list[dict[str, Any]],
    manifest_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    all_sample_ids = [clean_text(row.get("sample_id")) for row in manifest_rows if clean_text(row.get("sample_id"))]
    manifest_groups = Counter(sample_type_group(row.get("sample_type", "")) for row in manifest_rows)
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in observation_rows:
        grouped[clean_text(row.get("canonical_variant_id"))].append(row)

    output: list[dict[str, Any]] = []
    for canonical_id, rows in grouped.items():
        tumor_rows = [row for row in rows if row.get("sample_type_group") == "tumor"]
        wbc_rows = [row for row in rows if row.get("sample_type_group") == "wbc"]
        plasma_rows = [row for row in rows if row.get("sample_type_group") == "plasma"]
        other_rows = [row for row in rows if row.get("sample_type_group") == "other"]

        genes = collect_split_values(rows, "genes")
        chip_genes = sorted(set(genes) & CHIP_REVIEW_GENES)
        near_lod, lod_values = assay_lod_status(rows)
        plasma_status, af_trend, trend_meta = plasma_trend(plasma_rows)
        plasma_values = [value for value in trend_meta.get("values", []) if value is not None]

        tumor_concordance = (
            "tumor_concordant"
            if tumor_rows
            else ("tumor_context_available_variant_not_reported" if manifest_groups.get("tumor") else "tumor_not_provided")
        )
        wbc_concordance = (
            "wbc_concordant_chip_or_germline_possible"
            if wbc_rows
            else ("wbc_context_available_variant_not_reported" if manifest_groups.get("wbc") else "wbc_not_provided")
        )
        technical_hold = any(clean_text(row.get("qc_tier")) == "hold_or_reject_without_review" for row in rows)
        artifact_or_rescue = any(
            clean_text(row.get("source_trace_status"))
            in {"technical_rescue_required_before_source_assignment", "technical_artifact_origin_suspected"}
            for row in rows
        )
        if technical_hold or artifact_or_rescue:
            technical_readiness = "technical_rescue_or_local_alignment_review_required"
        elif near_lod:
            technical_readiness = "near_lod_or_low_signal_confirmation_required"
        else:
            technical_readiness = "technical_context_review_ready"

        if wbc_rows and tumor_rows:
            source_status = "tumor_and_wbc_shared_signal_unresolved"
        elif wbc_rows:
            source_status = "wbc_concordant_chip_or_germline_candidate"
        elif tumor_rows and plasma_rows and len(plasma_rows) > 1:
            source_status = "tumor_concordant_serial_mrd_trace_candidate"
        elif tumor_rows and plasma_rows:
            source_status = "tumor_concordant_single_timepoint_mrd_trace_candidate"
        elif plasma_rows and len(plasma_rows) > 1:
            source_status = "plasma_only_serial_signal_unresolved"
        elif plasma_rows:
            source_status = "plasma_only_single_timepoint_unresolved"
        elif tumor_rows:
            source_status = "tumor_reference_marker_only"
        elif other_rows:
            source_status = "nonstandard_sample_context_unresolved"
        else:
            source_status = "unresolved_missing_context"

        if source_status.startswith("tumor_concordant_serial") and technical_readiness == "technical_context_review_ready":
            source_confidence = "moderate_to_high_for_expert_review"
        elif source_status.startswith("tumor_concordant") and not wbc_rows:
            source_confidence = "moderate_requires_serial_or_lod_context"
        elif wbc_rows or technical_readiness != "technical_context_review_ready":
            source_confidence = "low_until_control_and_technical_review"
        elif plasma_rows:
            source_confidence = "low_to_moderate_unresolved"
        else:
            source_confidence = "context_only_not_assignable"

        if technical_readiness.startswith("technical_rescue"):
            next_step = "manual_caller_rescue_and_local_alignment_review"
        elif not manifest_groups.get("wbc"):
            next_step = "add_matched_wbc_or_buffy_coat_vcf"
        elif not manifest_groups.get("tumor"):
            next_step = "add_primary_tumor_or_baseline_marker_vcf"
        elif len(plasma_rows) < 2:
            next_step = "add_serial_plasma_ctdna_timepoints_or_repeat_sample"
        elif near_lod:
            next_step = "review_assay_lod_and_orthogonal_confirmation"
        else:
            next_step = "expert_source_tracing_ready"

        observed_samples = sorted({clean_text(row.get("manifest_sample_id")) for row in rows if clean_text(row.get("manifest_sample_id"))})
        missing_samples = [sample for sample in all_sample_ids if sample and sample not in observed_samples]
        first_plasma = trend_meta.get("first", {}) or {}
        last_plasma = trend_meta.get("last", {}) or {}
        tier_counts = Counter(clean_text(row.get("qc_tier")) for row in rows)
        source_counts = Counter(clean_text(row.get("source_trace_status")) for row in rows)

        output.append(
            {
                "case_id": clean_text(rows[0].get("case_id")),
                "canonical_variant_id": canonical_id,
                "representative_variant_id": clean_text(rows[0].get("variant_id")),
                "genes": ";".join(genes),
                "variant_type": ";".join(sorted({clean_text(row.get("type")) for row in rows if clean_text(row.get("type"))})),
                "source_assignment_status": source_status,
                "source_assignment_confidence": source_confidence,
                "technical_readiness": technical_readiness,
                "tumor_concordance_status": tumor_concordance,
                "wbc_concordance_status": wbc_concordance,
                "plasma_series_status": plasma_status,
                "plasma_af_trend": af_trend,
                "sample_observation_count": len(rows),
                "tumor_observation_count": len(tumor_rows),
                "wbc_observation_count": len(wbc_rows),
                "plasma_observation_count": len(plasma_rows),
                "first_plasma_timepoint": clean_text(first_plasma.get("timepoint")),
                "last_plasma_timepoint": clean_text(last_plasma.get("timepoint")),
                "first_plasma_af": first_plasma.get("tumor_af", ""),
                "last_plasma_af": last_plasma.get("tumor_af", ""),
                "min_plasma_af": min(plasma_values) if plasma_values else "",
                "max_plasma_af": max(plasma_values) if plasma_values else "",
                "near_lod_observation_count": near_lod,
                "assay_lod_values": ";".join(lod_values),
                "chip_review_genes": ";".join(chip_genes),
                "kg_match_statuses": ";".join(sorted({clean_text(row.get("kg_match_status")) for row in rows if clean_text(row.get("kg_match_status"))})),
                "kg_evidence_ids": ";".join(collect_split_values(rows, "kg_evidence_ids")),
                "qc_tier_counts": json.dumps(dict(tier_counts), ensure_ascii=False),
                "source_trace_status_counts": json.dumps(dict(source_counts), ensure_ascii=False),
                "caller_scenarios": ";".join(collect_split_values(rows, "caller_dependent_mrd_scenarios")),
                "recommended_next_step": next_step,
                "observed_samples": ";".join(observed_samples),
                "missing_or_not_reported_samples": ";".join(missing_samples),
                "observation_trace": summarize_observation_trace(rows),
            }
        )

    return sorted(
        output,
        key=lambda row: (
            clean_text(row.get("source_assignment_status")),
            -int(row.get("sample_observation_count") or 0),
            clean_text(row.get("canonical_variant_id")),
        ),
    )


def write_csv_rows(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def build_longitudinal_summary(
    case_id: str,
    manifest_path: Path,
    manifest_rows: list[dict[str, Any]],
    sample_summaries: list[dict[str, Any]],
    observation_rows: list[dict[str, Any]],
    longitudinal_rows: list[dict[str, Any]],
    missing_inputs: list[dict[str, str]],
) -> dict[str, Any]:
    source_counts = Counter(clean_text(row.get("source_assignment_status")) for row in longitudinal_rows)
    confidence_counts = Counter(clean_text(row.get("source_assignment_confidence")) for row in longitudinal_rows)
    sample_type_counts = Counter(sample_type_group(row.get("sample_type", "")) for row in manifest_rows)
    scenario_counts: Counter[str] = Counter()
    for row in longitudinal_rows:
        for scenario in clean_text(row.get("caller_scenarios")).split(";"):
            if scenario:
                scenario_counts[scenario] += 1

    tumor_ready = sum(1 for row in longitudinal_rows if clean_text(row.get("tumor_concordance_status")) == "tumor_concordant")
    wbc_flagged = sum(1 for row in longitudinal_rows if clean_text(row.get("wbc_concordance_status")).startswith("wbc_concordant"))
    serial_ready = sum(1 for row in longitudinal_rows if clean_text(row.get("plasma_series_status")) == "serial_plasma_signal_observed")
    near_lod = sum(1 for row in longitudinal_rows if int(row.get("near_lod_observation_count") or 0) > 0)

    if sample_type_counts.get("tumor") and sample_type_counts.get("wbc") and sample_type_counts.get("plasma", 0) >= 2:
        readiness = "LONGITUDINAL_SOURCE_TRACE_READY_FOR_EXPERT_REVIEW"
    elif len(manifest_rows) > 1:
        readiness = "PARTIAL_LONGITUDINAL_CONTEXT_AVAILABLE"
    else:
        readiness = "SINGLE_FILE_MANIFEST_SMOKE_TEST_ONLY"

    return {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "mode": "case_manifest_longitudinal_mrd_qc",
        "case_id": case_id,
        "case_manifest": str(manifest_path),
        "readiness": readiness,
        "interpretation_scope": "longitudinal source-tracing QC; not final MRD positivity or treatment recommendation",
        "manifest_sample_count": len(manifest_rows),
        "parsed_sample_count": len(sample_summaries),
        "missing_input_count": len(missing_inputs),
        "observation_row_count": len(observation_rows),
        "longitudinal_variant_count": len(longitudinal_rows),
        "sample_type_counts": dict(sample_type_counts),
        "source_assignment_counts": dict(source_counts),
        "source_assignment_confidence_counts": dict(confidence_counts),
        "caller_scenario_counts": dict(scenario_counts),
        "tumor_concordant_variant_count": tumor_ready,
        "wbc_concordant_variant_count": wbc_flagged,
        "serial_plasma_signal_variant_count": serial_ready,
        "near_lod_variant_count": near_lod,
        "sample_summaries": sample_summaries,
        "missing_inputs": missing_inputs,
        "required_manifest_columns": ["sample_id", "vcf_path"],
        "recommended_manifest_columns": MANIFEST_TEMPLATE_COLUMNS,
        "experiment_readiness_note": (
            "Use tumor tissue, matched WBC/buffy coat, and at least two serial plasma/cfDNA VCFs to test source-tracing lift. "
            "A one-row manifest only validates IO and baseline technical QC."
        ),
    }


def write_longitudinal_report(
    path: Path,
    summary: dict[str, Any],
    longitudinal_rows: list[dict[str, Any]],
    observation_csv: Path,
    longitudinal_csv: Path,
    summary_json: Path,
    html_path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# MRD Longitudinal Source-Tracing QC",
        "",
        f"- Case ID: `{summary['case_id']}`",
        f"- Generated: {summary['generated']}",
        f"- Readiness: `{summary['readiness']}`",
        f"- Manifest samples: {summary['manifest_sample_count']} parsed={summary['parsed_sample_count']} missing={summary['missing_input_count']}",
        f"- Observation rows: {summary['observation_row_count']}",
        f"- Longitudinal variants: {summary['longitudinal_variant_count']}",
        f"- Sample types: {counter_line(summary.get('sample_type_counts', {}))}",
        f"- Source assignments: {counter_line(summary.get('source_assignment_counts', {}))}",
        f"- Confidence: {counter_line(summary.get('source_assignment_confidence_counts', {}))}",
        "",
        "## Current Innovation Hypothesis",
        "",
        "- Build MRD QC as source tracing, not as a single positivity classifier: each caller signal is audited as tumor-concordant, WBC/CHIP-concordant, serially persistent, near-LOD, artifact-prone, or unresolved.",
        "- Separate three evidence layers: raw VCF/caller evidence, longitudinal patient context, and KG/RAG literature support. This avoids claiming KG support for workflow-only scenarios.",
        "- Use a manifest-driven tumor/WBC/plasma matrix so the same variant can be followed across baseline tissue, hematopoietic control, and serial cfDNA timepoints.",
        "- Treat missing context as a measurable QC output. Missing tumor, WBC, serial samples, and LOD become explicit experimental variables for ablation.",
        "",
        "## Experiments To Run",
        "",
        "1. Single VCF vs manifest longitudinal QC: compare unresolved/technical-rescue calls before and after adding tumor, WBC, and serial plasma files.",
        "2. WBC ablation: remove matched WBC rows from the manifest and quantify how many CHIP/hematopoietic-origin candidates become unresolved.",
        "3. Tumor concordance ablation: remove tumor tissue/baseline marker rows and quantify loss of tumor-informed source assignment.",
        "4. Serial sampling experiment: compare one plasma timepoint versus two or more timepoints for persistence, rising/falling VAF, and transient signal calls.",
        "5. LOD/low-VAF experiment: vary assay_lod in the manifest and count near-LOD calls requiring orthogonal confirmation.",
        "6. KG source-trace evidence experiment: merge the new source-trace evidence delta into GraphRAG and measure scenario-specific retrieval coverage.",
        "",
        "## Top Variant Rows",
        "",
    ]
    for row in longitudinal_rows[:20]:
        lines.append(
            f"- `{row.get('canonical_variant_id')}` {row.get('genes') or 'no_gene'}: "
            f"{row.get('source_assignment_status')} / {row.get('source_assignment_confidence')} / next={row.get('recommended_next_step')}"
        )
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            f"- Observation CSV: `{observation_csv}`",
            f"- Longitudinal variant CSV: `{longitudinal_csv}`",
            f"- Summary JSON: `{summary_json}`",
            f"- HTML dashboard: `{html_path}`",
            "",
            "This report is automated research QC, not a clinical diagnosis or treatment recommendation.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_longitudinal_html(
    path: Path,
    summary: dict[str, Any],
    longitudinal_rows: list[dict[str, Any]],
    observation_csv: Path,
    longitudinal_csv: Path,
) -> None:
    def esc(value: Any) -> str:
        return html.escape(clean_text(value))

    def bars(counter: dict[str, Any]) -> str:
        if not counter:
            return "<p class=\"muted\">not_available</p>"
        total = max(1, sum(int(value or 0) for value in counter.values()))
        items = []
        for key, value in sorted(counter.items(), key=lambda item: -int(item[1] or 0)):
            width = max(2, int((int(value or 0) / total) * 100))
            items.append(
                f"<div class=\"bar-row\"><span>{esc(key)}</span><strong>{int(value or 0)}</strong>"
                f"<div class=\"bar\"><i style=\"width:{width}%\"></i></div></div>"
            )
        return "\n".join(items)

    top_rows = longitudinal_rows[:60]
    row_html = "\n".join(
        "<tr>"
        f"<td>{esc(row.get('canonical_variant_id'))}</td>"
        f"<td>{esc(row.get('genes'))}</td>"
        f"<td>{esc(row.get('source_assignment_status'))}</td>"
        f"<td>{esc(row.get('tumor_concordance_status'))}</td>"
        f"<td>{esc(row.get('wbc_concordance_status'))}</td>"
        f"<td>{esc(row.get('plasma_series_status'))}</td>"
        f"<td>{esc(row.get('plasma_af_trend'))}</td>"
        f"<td>{esc(row.get('recommended_next_step'))}</td>"
        "</tr>"
        for row in top_rows
    )
    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>MRD Longitudinal Source-Tracing QC</title>
  <style>
    :root {{ color-scheme: light; --ink:#14212b; --muted:#607080; --line:#d8e0e7; --blue:#2a6fbb; --green:#1b8a5a; --amber:#b87516; --red:#b94040; --bg:#f7f9fb; }}
    * {{ box-sizing:border-box; }}
    body {{ margin:0; font-family:Segoe UI, Arial, sans-serif; color:var(--ink); background:var(--bg); }}
    header {{ padding:28px 34px 18px; background:#ffffff; border-bottom:1px solid var(--line); }}
    h1 {{ margin:0 0 8px; font-size:28px; letter-spacing:0; }}
    h2 {{ margin:0 0 14px; font-size:18px; letter-spacing:0; }}
    main {{ padding:22px 34px 38px; }}
    .muted {{ color:var(--muted); }}
    .kpis {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(190px,1fr)); gap:12px; margin:18px 0; }}
    .kpi {{ background:#fff; border:1px solid var(--line); border-radius:8px; padding:14px; }}
    .kpi b {{ display:block; font-size:28px; margin-top:6px; }}
    .grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(320px,1fr)); gap:16px; align-items:start; }}
    section {{ background:#fff; border:1px solid var(--line); border-radius:8px; padding:16px; margin-bottom:16px; }}
    .bar-row {{ display:grid; grid-template-columns:minmax(150px,1fr) 48px; gap:8px; align-items:center; margin:9px 0; }}
    .bar-row span {{ min-width:0; overflow-wrap:anywhere; }}
    .bar-row strong {{ text-align:right; }}
    .bar {{ grid-column:1 / -1; height:8px; background:#edf2f6; border-radius:999px; overflow:hidden; }}
    .bar i {{ display:block; height:100%; background:linear-gradient(90deg,var(--blue),var(--green)); }}
    table {{ width:100%; border-collapse:collapse; font-size:13px; }}
    th,td {{ padding:8px 9px; border-bottom:1px solid var(--line); text-align:left; vertical-align:top; }}
    th {{ background:#f0f4f7; position:sticky; top:0; }}
    .table-wrap {{ overflow:auto; max-height:620px; border:1px solid var(--line); border-radius:8px; }}
    .pill {{ display:inline-block; padding:4px 8px; border:1px solid var(--line); border-radius:999px; background:#fff; margin-right:6px; color:var(--muted); }}
    a {{ color:var(--blue); }}
  </style>
</head>
<body>
  <header>
    <h1>MRD Longitudinal Source-Tracing QC</h1>
    <div class="muted">Case {esc(summary.get('case_id'))} · {esc(summary.get('generated'))} · {esc(summary.get('readiness'))}</div>
  </header>
  <main>
    <div class="kpis">
      <div class="kpi">Manifest samples <b>{int(summary.get('manifest_sample_count') or 0)}</b></div>
      <div class="kpi">Observation rows <b>{int(summary.get('observation_row_count') or 0):,}</b></div>
      <div class="kpi">Longitudinal variants <b>{int(summary.get('longitudinal_variant_count') or 0):,}</b></div>
      <div class="kpi">Serial plasma variants <b>{int(summary.get('serial_plasma_signal_variant_count') or 0):,}</b></div>
      <div class="kpi">WBC-concordant variants <b>{int(summary.get('wbc_concordant_variant_count') or 0):,}</b></div>
      <div class="kpi">Near LOD variants <b>{int(summary.get('near_lod_variant_count') or 0):,}</b></div>
    </div>
    <section>
      <span class="pill">Observation CSV: {esc(observation_csv)}</span>
      <span class="pill">Longitudinal CSV: {esc(longitudinal_csv)}</span>
    </section>
    <div class="grid">
      <section>
        <h2>Source Assignment Counts</h2>
        {bars(summary.get('source_assignment_counts', {}))}
      </section>
      <section>
        <h2>Scenario Coverage</h2>
        {bars(summary.get('caller_scenario_counts', {}))}
      </section>
      <section>
        <h2>Sample Types</h2>
        {bars(summary.get('sample_type_counts', {}))}
      </section>
    </div>
    <section>
      <h2>Variant Matrix</h2>
      <div class="table-wrap">
        <table>
          <thead><tr><th>Variant</th><th>Genes</th><th>Source</th><th>Tumor</th><th>WBC</th><th>Plasma</th><th>Trend</th><th>Next step</th></tr></thead>
          <tbody>{row_html}</tbody>
        </table>
      </div>
    </section>
  </main>
</body>
</html>
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(html_text, encoding="utf-8")


def run_longitudinal_case_qc(args: argparse.Namespace) -> None:
    manifest_path = artifact_path(args.case_manifest)
    manifest_rows = read_case_manifest(manifest_path)
    case_id = args.case_id or clean_text(manifest_rows[0].get("case_id")) or manifest_path.stem
    evidence_context = load_evidence_context(artifact_path(args.evidence_csv))
    observation_rows: list[dict[str, Any]] = []
    sample_summaries: list[dict[str, Any]] = []
    missing_inputs: list[dict[str, str]] = []

    for manifest_row in manifest_rows:
        vcf_path = resolve_manifest_vcf_path(clean_text(manifest_row.get("vcf_path")), manifest_path)
        sample_id = clean_text(manifest_row.get("sample_id")) or vcf_path.stem
        if not vcf_path.exists():
            missing_inputs.append({"sample_id": sample_id, "vcf_path": str(vcf_path), "reason": "vcf_not_found"})
            continue
        variants, meta = parse_vcf(vcf_path)
        ncbi_annotation = None
        if args.ncbi_annotate:
            ncbi_annotation = annotate_variants_with_ncbi(variants, args.assembly, artifact_path(args.annotation_cache_dir))
        rows = [compute_variant_qc(variant, meta, evidence_context, ncbi_annotation) for variant in variants]
        sample_group = sample_type_group(manifest_row.get("sample_type", ""))
        sample_case_id = clean_text(manifest_row.get("case_id")) or case_id
        for row in rows:
            row.update(
                {
                    "case_id": sample_case_id,
                    "manifest_sample_id": sample_id,
                    "manifest_sample_type": clean_text(manifest_row.get("sample_type")),
                    "sample_type_group": sample_group,
                    "timepoint": clean_text(manifest_row.get("timepoint")),
                    "days_from_event": clean_text(manifest_row.get("days_from_event")),
                    "event": clean_text(manifest_row.get("event")),
                    "assay_type": clean_text(manifest_row.get("assay_type")),
                    "assay_lod": clean_text(manifest_row.get("assay_lod")),
                    "matched_tumor_id": clean_text(manifest_row.get("matched_tumor_id")),
                    "matched_wbc_id": clean_text(manifest_row.get("matched_wbc_id")),
                    "source_vcf": str(vcf_path),
                    "canonical_variant_id": canonical_variant_id_from_row(row),
                }
            )
        observation_rows.extend(rows)
        sample_summary = summarize_variants(rows, meta, evidence_context, ncbi_annotation)
        sample_summaries.append(
            {
                "sample_id": sample_id,
                "sample_type": clean_text(manifest_row.get("sample_type")),
                "sample_type_group": sample_group,
                "timepoint": clean_text(manifest_row.get("timepoint")),
                "days_from_event": clean_text(manifest_row.get("days_from_event")),
                "vcf_path": str(vcf_path),
                "record_count": sample_summary.get("record_count", 0),
                "annotated_variant_count": sample_summary.get("annotated_variant_count", 0),
                "gene_match_variant_count": sample_summary.get("gene_match_variant_count", 0),
                "matched_control_signal_count": sample_summary.get("matched_control_signal_count", 0),
                "tier_counts": sample_summary.get("tier_counts", {}),
                "source_trace_status": sample_summary.get("molecular_evidence_qc", {}).get("source_trace_status"),
            }
        )

    longitudinal_rows = build_longitudinal_variant_rows(observation_rows, manifest_rows)
    default_prefix = safe_case_component(case_id)
    observation_csv = artifact_path(args.variant_csv) if args.variant_csv else artifact_path(Path("data/qc") / f"{default_prefix}.longitudinal_observations.csv")
    longitudinal_csv = artifact_path(args.longitudinal_variant_csv) if args.longitudinal_variant_csv else artifact_path(Path("data/qc") / f"{default_prefix}.longitudinal_variants.csv")
    summary_json = artifact_path(args.summary_json) if args.summary_json else artifact_path(Path("data/qc") / f"{default_prefix}.longitudinal_summary.json")
    report_md = artifact_path(args.report_md) if args.report_md else artifact_path(Path("docs/qc") / f"{default_prefix}_longitudinal_mrd_qc.md")
    visual_html = artifact_path(args.visual_html) if args.visual_html else artifact_path(Path("docs/qc") / f"{default_prefix}_longitudinal_mrd_qc.html")

    summary = build_longitudinal_summary(
        case_id,
        manifest_path,
        manifest_rows,
        sample_summaries,
        observation_rows,
        longitudinal_rows,
        missing_inputs,
    )
    summary["observation_csv"] = str(observation_csv)
    summary["longitudinal_variant_csv"] = str(longitudinal_csv)
    summary["summary_json"] = str(summary_json)
    summary["report_md"] = str(report_md)
    summary["visual_html"] = str(visual_html)
    summary["source_trace_rag"] = {
        "enabled": False,
        "reason": "case_manifest_mode_uses_longitudinal_matrix; run per-sample or merged RAG after KG source-trace delta is integrated",
    }

    write_csv_rows(observation_csv, observation_rows, LONGITUDINAL_OBSERVATION_FIELDS)
    write_csv_rows(longitudinal_csv, longitudinal_rows, LONGITUDINAL_VARIANT_FIELDS)
    summary_json.parent.mkdir(parents=True, exist_ok=True)
    summary_json.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    write_longitudinal_report(report_md, summary, longitudinal_rows, observation_csv, longitudinal_csv, summary_json, visual_html)
    write_longitudinal_html(visual_html, summary, longitudinal_rows, observation_csv, longitudinal_csv)

    print(
        json.dumps(
            {
                "mode": summary["mode"],
                "readiness": summary["readiness"],
                "case_id": case_id,
                "manifest_sample_count": summary["manifest_sample_count"],
                "parsed_sample_count": summary["parsed_sample_count"],
                "observation_row_count": summary["observation_row_count"],
                "longitudinal_variant_count": summary["longitudinal_variant_count"],
                "source_assignment_counts": summary["source_assignment_counts"],
                "observation_csv": str(observation_csv),
                "longitudinal_variant_csv": str(longitudinal_csv),
                "summary_json": str(summary_json),
                "report_md": str(report_md),
                "visual_html": str(visual_html),
            },
            ensure_ascii=False,
            indent=2,
        )
    )


def safe_case_component(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", clean_text(value))
    text = re.sub(r"_+", "_", text).strip("._ ")
    return text or "longitudinal_case"


def run_single_vcf_qc(args: argparse.Namespace) -> None:
    if not args.vcf:
        raise SystemExit("--vcf is required unless --case-manifest is provided.")
    vcf_path = Path(args.vcf)
    evidence_csv = artifact_path(args.evidence_csv)
    variants, meta = parse_vcf(vcf_path)
    ncbi_annotation = None
    if args.ncbi_annotate:
        ncbi_annotation = annotate_variants_with_ncbi(variants, args.assembly, artifact_path(args.annotation_cache_dir))
    evidence_context = load_evidence_context(evidence_csv)
    rows = [compute_variant_qc(variant, meta, evidence_context, ncbi_annotation) for variant in variants]

    case_id = args.case_id or meta.get("tumor_sample") or vcf_path.stem
    variant_csv = artifact_path(args.variant_csv) if args.variant_csv else artifact_path(Path("data/qc") / f"{case_id}.qc_variants.csv")
    summary_json = artifact_path(args.summary_json) if args.summary_json else artifact_path(Path("data/qc") / f"{case_id}.qc_summary.json")
    report_md = artifact_path(args.report_md) if args.report_md else artifact_path(Path("docs/qc") / f"{case_id}_qc_report.md")
    case_report_md = artifact_path(args.case_report_md) if args.case_report_md else artifact_path(Path("docs/qc") / f"{case_id}_case_qc_note.md")
    visual_html = artifact_path(args.visual_html) if args.visual_html else artifact_path(Path("docs/qc") / f"{case_id}_source_trace_lanes.html")
    rag_json = artifact_path(args.rag_json) if args.rag_json else artifact_path(Path("docs/qc") / f"{case_id}_source_trace_rag.json")
    rag_md = artifact_path(args.rag_md) if args.rag_md else artifact_path(Path("docs/qc") / f"{case_id}_source_trace_rag.md")

    summary = summarize_variants(rows, meta, evidence_context, ncbi_annotation)
    summary["input_vcf"] = str(vcf_path)
    summary["variant_csv"] = str(variant_csv)
    summary["summary_json"] = str(summary_json)
    summary["report_md"] = str(report_md)
    summary["case_report_md"] = str(case_report_md)
    summary["visual_html"] = str(visual_html)
    summary["source_trace_rag"] = {"enabled": False, "json": str(rag_json), "md": str(rag_md)}
    if args.source_trace_rag:
        summary["source_trace_rag"] = run_source_trace_rag(summary, rows, args, rag_json, rag_md)

    write_variant_csv(variant_csv, rows)
    summary_json.parent.mkdir(parents=True, exist_ok=True)
    summary_json.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    write_report(report_md, summary, rows, vcf_path, variant_csv, summary_json)
    write_case_report(case_report_md, summary, rows, vcf_path, variant_csv, summary_json)
    write_visual_lanes_html(visual_html, summary, rows, vcf_path)
    case_manifest = collect_case_qc_outputs(
        case_id=case_id,
        sample_id=meta.get("tumor_sample", ""),
        summary_json=summary_json,
        variant_csv=variant_csv,
        qc_report_md=report_md,
        case_report_md=case_report_md,
        include_decision=False,
    )

    print(json.dumps({
        "overall_status": summary["overall_status"],
        "case_qc_class": summary["case_assessment"]["case_qc_class"],
        "mrd_signal_class": summary["case_assessment"]["mrd_signal_class"],
        "source_trace_status": summary["molecular_evidence_qc"]["source_trace_status"],
        "evidence_mode": summary["molecular_evidence_qc"]["evidence_mode"],
        "nlp_qc_conclusion": summary["molecular_evidence_qc"]["nlp_qc_conclusion"],
        "mrd_entry_gate": summary["mrd_clinical_interpretation_entry_gate"]["gate_class"],
        "mrd_entry_permission": summary["mrd_clinical_interpretation_entry_gate"]["entry_permission"],
        "record_count": summary["record_count"],
        "tier_counts": summary["tier_counts"],
        "kg_match_counts": summary["kg_match_counts"],
        "ncbi_selected_assembly": summary.get("ncbi_annotation", {}).get("selected_assembly"),
        "ncbi_exact_gene_interval_count": summary.get("ncbi_exact_gene_interval_count"),
        "variant_csv": str(variant_csv),
        "summary_json": str(summary_json),
        "report_md": str(report_md),
        "case_report_md": str(case_report_md),
        "visual_html": str(visual_html),
        "source_trace_rag": summary.get("source_trace_rag"),
        "case_qc_dir": case_manifest["qc_dir"],
    }, ensure_ascii=False, indent=2))


def main() -> None:
    args = parse_args()
    if args.write_manifest_template:
        template_path = artifact_path(args.write_manifest_template)
        write_case_manifest_template(template_path)
        if not args.vcf and not args.case_manifest:
            print(json.dumps({"manifest_template": str(template_path)}, ensure_ascii=False, indent=2))
            return
    if args.case_manifest:
        run_longitudinal_case_qc(args)
        return
    run_single_vcf_qc(args)


if __name__ == "__main__":
    main()
