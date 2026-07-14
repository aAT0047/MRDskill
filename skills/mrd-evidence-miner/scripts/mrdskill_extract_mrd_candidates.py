#!/usr/bin/env python3
"""Extract traceable MRD evidence candidates from kg_ready_no_pdf text files.

This script is intentionally conservative. It creates candidate evidence rows
from title/abstract metadata and nearby text snippets; it does not claim final
confidence or replace expert/LLM review.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable


RAW_COLUMNS = [
    "evidence_id",
    "paper_id",
    "title",
    "year",
    "journal",
    "doi",
    "pmid",
    "pmcid",
    "cancer_type",
    "cancer_subtype",
    "disease_stage",
    "treatment_setting",
    "gene",
    "gene_symbol_normalized",
    "variant",
    "variant_type",
    "mutation_hotspot",
    "genomic_coordinate",
    "pathway",
    "biomarker_type",
    "assay_type",
    "assay_platform",
    "tumor_informed_or_tumor_naive",
    "sample_type",
    "VAF_cutoff",
    "limit_of_detection",
    "assay_sensitivity",
    "assay_specificity",
    "matched_WBC_used",
    "germline_filtering_used",
    "CHIP_filtering_used",
    "mrd_relevance",
    "relapse_association",
    "recurrence_association",
    "prognosis_association",
    "treatment_monitoring_value",
    "resistance_association",
    "recurrence_lead_time",
    "survival_endpoint",
    "hazard_ratio",
    "confidence_interval",
    "p_value",
    "cohort_type",
    "sample_size",
    "validation_status",
    "evidence_source_type_preliminary",
    "source_section",
    "source_table",
    "source_figure",
    "supplementary_file",
    "supplementary_table",
    "quoted_or_paraphrased_evidence",
    "limitation",
    "extraction_notes",
]


MRD_TERMS = [
    r"\bMRD\b",
    r"minimal residual disease",
    r"measurable residual disease",
    r"molecular residual disease",
    r"\bctDNA\b",
    r"circulating tumou?r DNA",
    r"liquid biopsy",
    r"molecular relapse",
]

OUTCOME_TERMS = [
    "recurrence",
    "relapse",
    "progression-free survival",
    "disease-free survival",
    "overall survival",
    "prognosis",
    "treatment response",
    "clearance",
    "persistence",
    "reappearance",
    "surveillance",
]

CHIP_GENES = {
    "DNMT3A",
    "TET2",
    "ASXL1",
    "JAK2",
    "TP53",
    "PPM1D",
    "SF3B1",
    "SRSF2",
    "U2AF1",
    "IDH1",
    "IDH2",
    "ATM",
    "CHEK2",
    "GNAS",
    "KRAS",
}

GENE_RE = re.compile(
    r"\b(?:DNMT3A|TET2|ASXL1|JAK2|TP53|PPM1D|SF3B1|SRSF2|U2AF1|IDH1|IDH2|ATM|CHEK2|GNAS|KRAS|BRAF|EGFR|ALK|ROS1|RET|MET|ERBB2|HER2|PIK3CA|APC|NRAS|NPM1|FLT3|KIT|WT1|RUNX1|CEBPA|IGH|IGK|IGL|ESR1|BRCA1|BRCA2|PALB2|PTEN|RB1|SMAD4|CTNNB1|FBXW7|ARID1A|CDKN2A|CDK4|CDK6|MDM2|MYC|CCND1|FGFR1|FGFR2|FGFR3|NTRK1|NTRK2|NTRK3|ERBB3|ERBB4|MAP2K1|MAP2K2|AKT1|AKT2|AKT3|MTOR|NF1|NF2|VHL|TERT|NOTCH1|NOTCH2|KMT2A|NUP98|MLL|PML|RARA)\b"
)
VARIANT_RE = re.compile(
    r"\b(?:p\.[A-Za-z]{1,3}\d+[A-Za-z]{1,3}|p\.[A-Z]\d+[A-Z*]|[A-Z]\d+[A-Z*]|c\.\d+(?:[+-]\d+)?[ACGT]>[ACGT]|c\.\d+_\d+(?:del|ins|dup)[A-Z]*|V600E|G12[ACDRV]|G13[CD]|Q61[HKLR]|L858R|T790M|E746_A750del|exon\s+\d+|fusion|breakpoint|methylation|fragmentomics|chr(?:omosome)?\s*\d+[:_]\d+|[cg]\.\d+[ACGT]>[ACGT])\b",
    re.IGNORECASE,
)
HR_RE = re.compile(r"\bHR\s*[:=]?\s*([0-9]+(?:\.[0-9]+)?)", re.IGNORECASE)
CI_RE = re.compile(
    r"(?:95\s*%\s*)?C[Il]\s*[:=]?\s*([0-9]+(?:\.[0-9]+)?\s*[-–]\s*[0-9]+(?:\.[0-9]+)?)",
    re.IGNORECASE,
)
P_RE = re.compile(r"\bp\s*(?:=|<|<=)\s*([0-9]*\.?[0-9]+)", re.IGNORECASE)
N_RE = re.compile(
    r"\b(?:n\s*=\s*|patients?|cohort of|samples?)\s*(\d{2,5})\b",
    re.IGNORECASE,
)
PMID_FROM_URL_RE = re.compile(r"pubmed\.ncbi\.nlm\.nih\.gov/(\d+)")
KWRAW_RE = re.compile(r"(KWRAW-\d+)")


CANCER_PATTERNS = [
    ("colorectal cancer", r"colorectal|colon|rectal"),
    ("non-small cell lung cancer", r"non-small cell lung|NSCLC|lung adenocarcinoma"),
    ("breast cancer", r"breast"),
    ("acute myeloid leukemia", r"acute myeloid|AML"),
    ("myelodysplastic syndrome", r"myelodysplastic|MDS"),
    ("lymphoma", r"lymphoma|DLBCL|diffuse large B-cell"),
    ("pancreatic cancer", r"pancreatic"),
    ("urothelial carcinoma", r"urothelial|bladder"),
    ("gastric cancer", r"gastric|stomach"),
    ("ovarian cancer", r"ovarian"),
    ("melanoma", r"melanoma"),
    ("hepatocellular carcinoma", r"hepatocellular|HCC"),
    ("esophageal cancer", r"esophageal|oesophageal"),
    ("head and neck cancer", r"head and neck"),
    ("multiple myeloma", r"multiple myeloma|myeloma"),
]


ASSAY_PATTERNS = [
    ("ddPCR", r"\bddPCR\b|droplet digital"),
    ("targeted_NGS", r"targeted (?:next-generation sequencing|NGS)|targeted sequencing"),
    ("ultra_deep_sequencing", r"ultra[- ]deep"),
    ("WES", r"whole[- ]exome|WES"),
    ("WGS", r"whole[- ]genome|WGS"),
    ("methylation_assay", r"methylation"),
    ("personalized_panel", r"personalized|patient-specific"),
    ("tumor_informed_assay", r"tumou?r-informed|tumou?r informed"),
    ("tumor_naive_assay", r"tumou?r-na[iï]ve|tumou?r naive"),
    ("fragmentomics_assay", r"fragmentomics|fragmentation"),
    ("hybrid_capture_panel", r"hybrid capture"),
    ("amplicon_panel", r"amplicon"),
    ("IG_NGS", r"IG[- ]NGS|immunoglobulin gene"),
]


def normalize_text(text: str) -> str:
    text = text.replace("\r", "\n")
    text = re.sub(r"[ \t]+", " ", text)
    text = re.sub(r"\n{3,}", "\n\n", text)
    return text.strip()


def first_match(patterns: Iterable[tuple[str, str]], text: str) -> str:
    for label, pattern in patterns:
        if re.search(pattern, text, re.IGNORECASE):
            return label
    return "not_reported"


def all_gene_matches(text: str) -> list[str]:
    return sorted(set(match.group(0).upper() for match in GENE_RE.finditer(text)))


def first_regex(pattern: re.Pattern[str], text: str) -> str:
    match = pattern.search(text)
    if not match:
        return "not_reported"
    return match.group(1).strip()


def load_candidate_table(path: Path) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    if not path.exists():
        return rows
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            record_id = row.get("record_id", "")
            if record_id:
                rows[record_id] = row
    return rows


def load_converted_manifest(path: Path, kg_ready_root: Path) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    if not path.exists():
        return rows
    converted_dir = kg_ready_root / "downloaded_pdfs" / "converted_from_pdf"
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            source_filename = row.get("source_filename", "")
            match = KWRAW_RE.search(source_filename)
            if not match:
                continue
            paper_id = match.group(1)
            stem = row.get("output_stem", "")
            text_path = converted_dir / f"{stem}.txt"
            json_path = converted_dir / f"{stem}.json"
            html_path = converted_dir / f"{stem}.html"
            rows[paper_id] = {
                **row,
                "kg_ready_text_path": str(text_path),
                "kg_ready_json_path": str(json_path),
                "kg_ready_html_path": str(html_path),
            }
    return rows


def iter_mrd_snippets(text: str, context_chars: int) -> Iterable[tuple[str, str]]:
    term_re = re.compile("|".join(MRD_TERMS), re.IGNORECASE)
    seen: set[str] = set()
    for match in term_re.finditer(text):
        start = max(0, match.start() - context_chars)
        end = min(len(text), match.end() + context_chars)
        snippet = normalize_text(text[start:end])
        key = snippet[:180]
        if key in seen:
            continue
        seen.add(key)
        section = infer_section(text[: match.start()])
        yield section, snippet


def infer_section(prefix: str) -> str:
    tail = prefix[-4000:]
    headings = [
        ("Results", r"\n\s*Results\s*\n"),
        ("Methods", r"\n\s*(?:Methods|Materials and methods)\s*\n"),
        ("Discussion", r"\n\s*Discussion\s*\n"),
        ("Abstract", r"\n\s*Abstract\s*\n"),
        ("Introduction", r"\n\s*Introduction\s*\n"),
    ]
    best = ("unknown", -1)
    for label, pattern in headings:
        matches = list(re.finditer(pattern, tail, re.IGNORECASE))
        if matches and matches[-1].start() > best[1]:
            best = (label, matches[-1].start())
    return best[0]


def preliminary_source_type(section: str, text: str) -> str:
    if section in {"Supplementary_Table", "Supplementary_Figure", "Supplementary_Dataset"}:
        return "supplementary_result"
    if re.search(r"\bvalidation\b|prospective|clinical trial|multicenter|randomi[sz]ed", text, re.IGNORECASE):
        return "clinical_validation"
    if section in {"Results"}:
        return "primary_result"
    if section in {"Methods"}:
        return "methods_assay_detail"
    if section in {"Discussion", "Introduction"}:
        return "discussion_or_background_claim"
    if "systematic review" in text.lower() or "meta-analysis" in text.lower():
        return "meta_analysis"
    return "uncertain"


def mrd_relevance(text: str) -> str:
    lowered = text.lower()
    if "chip" in lowered or "clonal hematopoiesis" in lowered or "haematopoiesis" in lowered:
        return "CHIP_related"
    if "recurrence" in lowered or "surveillance" in lowered:
        return "recurrence_marker"
    if "relapse" in lowered:
        return "relapse_driver"
    if "prognos" in lowered or "survival" in lowered:
        return "prognostic_marker"
    if "response" in lowered or "clearance" in lowered or "persistence" in lowered:
        return "treatment_monitoring_marker"
    if "methylation" in lowered:
        return "methylation_marker"
    if "fragment" in lowered:
        return "fragmentomics_marker"
    if GENE_RE.search(text) and VARIANT_RE.search(text) and re.search(r"\bctdna\b|plasma|tumou?r|tumor|tissue|vaf|variant allele", text, re.IGNORECASE):
        return "tumor_specific_marker"
    return "uncertain"


def tumor_informed_status(text: str) -> str:
    if re.search(r"tumou?r-informed|tumou?r informed|personalized|patient-specific", text, re.IGNORECASE):
        return "tumor_informed"
    if re.search(r"tumou?r-na[iï]ve|tumou?r naive", text, re.IGNORECASE):
        return "tumor_naive"
    return "not_reported"


def yes_no_from_terms(text: str, terms: Iterable[str]) -> str:
    lowered = text.lower()
    return "yes" if any(term.lower() in lowered for term in terms) else "not_reported"


def build_row(
    idx: int,
    paper_id: str,
    metadata: dict[str, str],
    section: str,
    snippet: str,
    max_quote_chars: int,
) -> dict[str, str]:
    combined = " ".join(
        [
            metadata.get("title", ""),
            metadata.get("abstract_if_available", ""),
            metadata.get("likely_topic_tags", ""),
            snippet,
        ]
    )
    genes = all_gene_matches(combined)
    gene = ";".join(genes) if genes else "not_reported"
    chip_overlap = sorted(set(genes) & CHIP_GENES)
    variant_match = VARIANT_RE.search(combined)
    variant = variant_match.group(0) if variant_match else "not_reported"

    row = {column: "not_reported" for column in RAW_COLUMNS}
    row.update(
        {
            "evidence_id": f"EV-{idx:07d}",
            "paper_id": paper_id,
            "title": metadata.get("title") or "not_reported",
            "year": metadata.get("year") or "not_reported",
            "journal": metadata.get("journal") or "not_reported",
            "doi": metadata.get("doi") or "not_reported",
            "pmid": metadata.get("pmid") or first_regex(PMID_FROM_URL_RE, metadata.get("landing_page_url", "")),
            "pmcid": metadata.get("pmcid") or "not_reported",
            "cancer_type": first_match(CANCER_PATTERNS, combined),
            "cancer_subtype": "not_reported",
            "disease_stage": first_stage(combined),
            "treatment_setting": first_treatment_setting(combined),
            "gene": gene,
            "gene_symbol_normalized": gene,
            "variant": variant,
            "variant_type": infer_variant_type(variant, combined),
            "biomarker_type": infer_biomarker_type(combined, bool(chip_overlap)),
            "assay_type": first_match(ASSAY_PATTERNS, combined),
            "assay_platform": "not_reported",
            "tumor_informed_or_tumor_naive": tumor_informed_status(combined),
            "sample_type": infer_sample_type(combined),
            "VAF_cutoff": first_vaf(combined),
            "limit_of_detection": first_lod(combined),
            "matched_WBC_used": yes_no_from_terms(combined, ["matched WBC", "matched white blood", "buffy coat"]),
            "germline_filtering_used": yes_no_from_terms(combined, ["germline filtering", "germline variant"]),
            "CHIP_filtering_used": yes_no_from_terms(combined, ["CHIP", "clonal hematopoiesis", "clonal haematopoiesis"]),
            "mrd_relevance": mrd_relevance(combined),
            "relapse_association": "reported" if re.search(r"relapse", combined, re.IGNORECASE) else "not_reported",
            "recurrence_association": "reported" if re.search(r"recurrence", combined, re.IGNORECASE) else "not_reported",
            "prognosis_association": "reported"
            if re.search(r"prognos|survival|PFS|DFS|OS", combined, re.IGNORECASE)
            else "not_reported",
            "treatment_monitoring_value": "reported"
            if re.search(r"treatment response|clearance|persistence|monitoring", combined, re.IGNORECASE)
            else "not_reported",
            "recurrence_lead_time": first_lead_time(combined),
            "survival_endpoint": infer_survival_endpoint(combined),
            "hazard_ratio": first_regex(HR_RE, combined),
            "confidence_interval": first_regex(CI_RE, combined),
            "p_value": first_regex(P_RE, combined),
            "sample_size": first_regex(N_RE, combined),
            "cohort_type": infer_cohort_type(combined, metadata),
            "validation_status": infer_validation(combined),
            "evidence_source_type_preliminary": preliminary_source_type(section, combined),
            "source_section": section,
            "quoted_or_paraphrased_evidence": snippet[:max_quote_chars],
            "limitation": infer_limitation(section, combined, genes),
            "extraction_notes": "python_candidate_extraction; requires domain review",
        }
    )
    return row


def first_stage(text: str) -> str:
    match = re.search(r"\bstage\s+(?:I{1,3}V?|[1-4])\b", text, re.IGNORECASE)
    return match.group(0) if match else "not_reported"


def first_treatment_setting(text: str) -> str:
    for label, pattern in [
        ("postoperative", r"post[- ]?operative|after surgery|postsurgery"),
        ("adjuvant", r"adjuvant"),
        ("neoadjuvant", r"neoadjuvant"),
        ("post-transplant", r"post[- ]?transplant"),
        ("metastatic", r"metastatic|advanced disease"),
        ("surveillance", r"surveillance"),
    ]:
        if re.search(pattern, text, re.IGNORECASE):
            return label
    return "not_reported"


def infer_variant_type(variant: str, text: str) -> str:
    lowered = text.lower()
    if variant == "not_reported":
        if "methylation" in lowered:
            return "methylation"
        if "fragment" in lowered:
            return "fragmentomics"
        return "not_reported"
    if "fusion" in variant.lower() or "breakpoint" in variant.lower():
        return "fusion"
    if "methylation" in variant.lower():
        return "methylation"
    if re.search(r"[ACGT]>", variant, re.IGNORECASE):
        return "SNV"
    if re.search(r"[A-Z]\d+[A-Z]|p\.", variant, re.IGNORECASE):
        return "protein_change"
    return "uncertain"


def infer_biomarker_type(text: str, has_chip_gene: bool) -> str:
    lowered = text.lower()
    if "chip" in lowered or "clonal hematopoiesis" in lowered or "clonal haematopoiesis" in lowered:
        return "CHIP_related"
    if "recurrence" in lowered:
        return "recurrence_marker"
    if "prognos" in lowered or "survival" in lowered:
        return "prognostic_marker"
    if "response" in lowered or "clearance" in lowered:
        return "treatment_response_marker"
    if "ctdna" in lowered or "circulating tumour dna" in lowered or "circulating tumor dna" in lowered:
        return "tumor-derived_ctDNA"
    if has_chip_gene:
        return "technical_noise_candidate"
    return "uncertain"


def infer_sample_type(text: str) -> str:
    for label, pattern in [
        ("plasma", r"\bplasma\b"),
        ("serum", r"\bserum\b"),
        ("bone_marrow", r"bone marrow"),
        ("peripheral_blood", r"peripheral blood"),
        ("tissue", r"\btissue\b|tumou?r tissue"),
        ("matched_WBC", r"matched WBC|white blood cell"),
        ("buffy_coat", r"buffy coat"),
    ]:
        if re.search(pattern, text, re.IGNORECASE):
            return label
    return "not_reported"


def first_vaf(text: str) -> str:
    match = re.search(r"(?:VAF|variant allele fraction)[^\n]{0,40}?([0-9]*\.?[0-9]+\s*%)", text, re.IGNORECASE)
    return match.group(1) if match else "not_reported"


def first_lod(text: str) -> str:
    match = re.search(
        r"(?:limit of detection|LOD)[^\n]{0,40}?([0-9]*\.?[0-9]+\s*%|10\^-?\d+|10-\d+)",
        text,
        re.IGNORECASE,
    )
    return match.group(1) if match else "not_reported"


def first_lead_time(text: str) -> str:
    match = re.search(
        r"(\d+(?:\.\d+)?\s*(?:days?|weeks?|months?))[^.\n]{0,80}(?:before|prior to|preced)",
        text,
        re.IGNORECASE,
    )
    return match.group(1) if match else "not_reported"


def infer_survival_endpoint(text: str) -> str:
    endpoints = []
    for label, pattern in [
        ("PFS", r"\bPFS\b|progression-free survival"),
        ("DFS", r"\bDFS\b|disease-free survival"),
        ("RFS", r"\bRFS\b|recurrence-free survival"),
        ("OS", r"\bOS\b|overall survival"),
    ]:
        if re.search(pattern, text, re.IGNORECASE):
            endpoints.append(label)
    return ";".join(endpoints) if endpoints else "not_reported"


def infer_cohort_type(text: str, metadata: dict[str, str]) -> str:
    combined = " ".join([text, metadata.get("publication_type", ""), metadata.get("article_type_guess", "")])
    for label, pattern in [
        ("prospective", r"prospective"),
        ("randomized_controlled_trial", r"randomi[sz]ed|clinical trial"),
        ("multicenter", r"multicenter|multi-centre|multi-center"),
        ("retrospective", r"retrospective"),
        ("review", r"review|perspective"),
    ]:
        if re.search(pattern, combined, re.IGNORECASE):
            return label
    return "not_reported"


def infer_validation(text: str) -> str:
    if re.search(r"validation cohort|independent validation|validated", text, re.IGNORECASE):
        return "validated"
    if re.search(r"clinical trial|prospective|randomi[sz]ed", text, re.IGNORECASE):
        return "clinical_study"
    return "not_reported"


def infer_limitation(section: str, text: str, genes: list[str]) -> str:
    limitations = []
    if section == "unknown":
        limitations.append("source_section_uncertain")
    if not genes:
        limitations.append("no_gene_detected")
    if not any(term in text.lower() for term in OUTCOME_TERMS):
        limitations.append("outcome_not_explicit_in_snippet")
    if "review" in text.lower() or "perspective" in text.lower():
        limitations.append("may_be_review_or_background")
    return ";".join(limitations) if limitations else "not_reported"


def write_csv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_jsonl(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Extract candidate MRD evidence rows from kg_ready_no_pdf corpus."
    )
    parser.add_argument("--corpus-dir", required=True, help="Path to kg_ready_no_pdf.")
    parser.add_argument("--out-csv", default="data/evidence/mrd_evidence_raw.csv")
    parser.add_argument("--out-jsonl", default="data/evidence/mrd_evidence_raw.jsonl")
    parser.add_argument("--max-papers", type=int, default=0, help="0 means all papers.")
    parser.add_argument("--max-snippets-per-paper", type=int, default=6)
    parser.add_argument("--context-chars", type=int, default=900)
    parser.add_argument("--max-quote-chars", type=int, default=1200)
    args = parser.parse_args()

    corpus_dir = Path(args.corpus_dir).resolve()
    manifests = corpus_dir / "manifests"
    candidate_rows = load_candidate_table(manifests / "keyword_research_candidate_table.csv")
    converted_rows = load_converted_manifest(
        manifests / "downloaded_pdfs_converted_manifest.csv",
        corpus_dir,
    )

    evidence_rows: list[dict[str, str]] = []
    per_paper_counts: defaultdict[str, int] = defaultdict(int)
    paper_ids = sorted(converted_rows)
    if args.max_papers > 0:
        paper_ids = paper_ids[: args.max_papers]

    for paper_id in paper_ids:
        converted = converted_rows[paper_id]
        text_path = Path(converted["kg_ready_text_path"])
        if not text_path.exists():
            continue
        metadata = candidate_rows.get(paper_id, {})
        text = normalize_text(text_path.read_text(encoding="utf-8", errors="replace"))
        if not text:
            continue

        for section, snippet in iter_mrd_snippets(text, args.context_chars):
            if per_paper_counts[paper_id] >= args.max_snippets_per_paper:
                break
            row = build_row(
                len(evidence_rows) + 1,
                paper_id,
                metadata,
                section,
                snippet,
                args.max_quote_chars,
            )
            evidence_rows.append(row)
            per_paper_counts[paper_id] += 1

    write_csv(Path(args.out_csv), evidence_rows, RAW_COLUMNS)
    write_jsonl(Path(args.out_jsonl), evidence_rows)
    print(f"Wrote raw evidence CSV: {Path(args.out_csv).resolve()} ({len(evidence_rows)} rows)")
    print(f"Wrote raw evidence JSONL: {Path(args.out_jsonl).resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
