#!/usr/bin/env python3
"""Build an MRD source-tracing benchmark from Razavi et al. 2019 supplements.

The benchmark is centered on cfDNA variants and labels each call by whether the
same patient-coordinate-allele is present in matched WBC and/or tumor tissue.
It is intended for evaluating MRD-specific caller-evidence QC skills, not for
making clinical calls.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd

from mrdskill_project_paths import artifact_root

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


ARTIFACTS = artifact_root()
DEFAULT_INPUT = ARTIFACTS / "inputs" / "razavi_2019_supplementary_tables.xlsx"
DEFAULT_OUT_CSV = ARTIFACTS / "data" / "benchmarks" / "razavi_2019_mrd_source_trace_skill_benchmark.csv"
DEFAULT_OUT_JSON = ARTIFACTS / "data" / "benchmarks" / "razavi_2019_mrd_source_trace_skill_benchmark_summary.json"
DEFAULT_OUT_MD = ARTIFACTS / "docs" / "benchmarks" / "razavi_2019_mrd_source_trace_skill_benchmark.md"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build Razavi 2019 MRD source-trace skill benchmark.")
    parser.add_argument("--input-xlsx", default=str(DEFAULT_INPUT), help="Razavi supplementary tables XLSX.")
    parser.add_argument("--out-csv", default=str(DEFAULT_OUT_CSV), help="Output benchmark CSV.")
    parser.add_argument("--out-summary-json", default=str(DEFAULT_OUT_JSON), help="Output summary JSON.")
    parser.add_argument("--out-md", default=str(DEFAULT_OUT_MD), help="Output Markdown report.")
    return parser.parse_args()


def clean(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    return str(value).strip()


def to_int(value: Any) -> int | None:
    text = clean(value)
    if not text:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def to_float(value: Any) -> float | None:
    text = clean(value)
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def normalize_chrom(value: Any) -> str:
    text = clean(value)
    if text.lower().startswith("chr"):
        text = text[3:]
    return text.upper()


def variant_key(row: dict[str, Any]) -> tuple[str, str, int, str, str]:
    patient = clean(row.get("Tumor_Sample_Barcode"))
    chrom = normalize_chrom(row.get("Chromosome"))
    pos = to_int(row.get("Start_Position")) or -1
    ref = clean(row.get("Reference_Allele")).upper()
    alt = clean(row.get("Tumor_Seq_Allele2")).upper()
    return patient, chrom, pos, ref, alt


def variant_id_from_key(key: tuple[str, str, int, str, str]) -> str:
    _, chrom, pos, ref, alt = key
    return f"{chrom}:{pos}:{ref}>{alt}"


def load_sheet(path: Path, sheet_name: str) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name=sheet_name, header=1)
    return df.dropna(how="all").copy()


def first_by_key(df: pd.DataFrame) -> dict[tuple[str, str, int, str, str], dict[str, Any]]:
    result: dict[tuple[str, str, int, str, str], dict[str, Any]] = {}
    for row in df.to_dict(orient="records"):
        key = variant_key(row)
        if key[0] and key[2] > 0 and key[3] and key[4] and key not in result:
            result[key] = row
    return result


def af_from_counts(alt: Any, depth: Any) -> float | None:
    alt_i = to_int(alt)
    depth_i = to_int(depth)
    if alt_i is None or depth_i is None or depth_i <= 0:
        return None
    return alt_i / depth_i


def baseline_label_from_chip_gene(gene: str) -> str:
    return "chip_or_hematopoietic" if clean(gene).upper() in CHIP_REVIEW_GENES else "not_chip_gene"


def baseline_label_from_vaf(cfdna_af: float | None) -> str:
    if cfdna_af is None:
        return "not_assessable"
    if cfdna_af < 0.005:
        return "low_vaf_possible_noise_or_chip"
    if cfdna_af >= 0.02:
        return "higher_vaf_possible_tumor"
    return "intermediate_vaf_unresolved"


def build_rows(cf_df: pd.DataFrame, wbc_df: pd.DataFrame, tumor_df: pd.DataFrame) -> list[dict[str, Any]]:
    wbc_by_key = first_by_key(wbc_df)
    tumor_by_key = first_by_key(tumor_df)
    rows: list[dict[str, Any]] = []

    for raw in cf_df.to_dict(orient="records"):
        key = variant_key(raw)
        if not key[0] or key[2] <= 0 or not key[3] or not key[4]:
            continue
        wbc = wbc_by_key.get(key)
        tumor = tumor_by_key.get(key)
        gene = clean(raw.get("Hugo_Symbol")).upper()
        cfdna_af = af_from_counts(raw.get("t_alt_count"), raw.get("t_depth"))
        wbc_af = af_from_counts(wbc.get("n_alt_count"), wbc.get("n_depth")) if wbc else None
        tumor_af = af_from_counts(tumor.get("t_alt_count"), tumor.get("t_depth")) if tumor else None
        has_wbc = wbc is not None
        has_tumor = tumor is not None

        if has_tumor and has_wbc:
            label = "tumor_and_wbc_shared_signal_unresolved"
            binary_origin = "ambiguous_shared"
        elif has_tumor:
            label = "tumor_concordant_cfDNA_signal"
            binary_origin = "tumor_derived"
        elif has_wbc:
            label = "wbc_concordant_chip_or_hematopoietic_signal"
            binary_origin = "chip_or_hematopoietic"
        else:
            label = "cfDNA_only_unresolved_signal"
            binary_origin = "unresolved_cfDNA_only"

        low_support = (to_int(raw.get("t_alt_count")) or 0) < 3
        rows.append(
            {
                "benchmark_id": "razavi_2019_mrd_source_trace_skill",
                "paper_id": "KWRAW-000862",
                "title": "High-intensity sequencing reveals the sources of plasma circulating cell-free DNA variants",
                "doi": "10.1038/s41591-019-0652-7",
                "patient_id": key[0],
                "tissue": clean(raw.get("Matched_Norm_Sample_Barcode")),
                "gene": gene,
                "chrom": key[1],
                "pos": key[2],
                "ref": key[3],
                "alt": key[4],
                "variant_id": variant_id_from_key(key),
                "sample_type": "cfDNA",
                "assembly": clean(raw.get("NCBI_Build")),
                "variant_classification": clean(raw.get("Variant_Classification")),
                "variant_type": clean(raw.get("Variant_Type")),
                "hgvsc": clean(raw.get("HGVSc")),
                "hgvsp": clean(raw.get("HGVSp")),
                "cfdna_depth": to_int(raw.get("t_depth")),
                "cfdna_ref_count": to_int(raw.get("t_ref_count")),
                "cfdna_alt_count": to_int(raw.get("t_alt_count")),
                "cfdna_af": cfdna_af,
                "wbc_matched": int(has_wbc),
                "wbc_depth": to_int(wbc.get("n_depth")) if wbc else "",
                "wbc_ref_count": to_int(wbc.get("n_ref_count")) if wbc else "",
                "wbc_alt_count": to_int(wbc.get("n_alt_count")) if wbc else "",
                "wbc_af": wbc_af if wbc_af is not None else "",
                "tumor_tissue_matched": int(has_tumor),
                "tumor_depth": to_int(tumor.get("t_depth")) if tumor else "",
                "tumor_ref_count": to_int(tumor.get("t_ref_count")) if tumor else "",
                "tumor_alt_count": to_int(tumor.get("t_alt_count")) if tumor else "",
                "tumor_af": tumor_af if tumor_af is not None else "",
                "source_trace_gold_label": label,
                "binary_origin_label": binary_origin,
                "chip_gene_overlap": int(gene in CHIP_REVIEW_GENES),
                "low_alt_support": int(low_support),
                "skill_input_available_context": "cfDNA_variant;matched_wbc_table;tumor_tissue_table;gene;depth;allele_counts",
                "expected_skill_output": label,
                "recommended_review_action": (
                    "matched_wbc_review"
                    if has_wbc and not has_tumor
                    else "tumor_concordance_review"
                    if has_tumor and not has_wbc
                    else "expert_review_shared_tumor_wbc_signal"
                    if has_tumor and has_wbc
                    else "unresolved_add_context_or_repeat_sampling"
                ),
                "baseline_chip_gene_label": baseline_label_from_chip_gene(gene),
                "baseline_vaf_label": baseline_label_from_vaf(cfdna_af),
                "source_tables": "cfDNA=Supplementary Table 11;WBC=Supplementary Table 12;tumor=Supplementary Table 13",
            }
        )
    return rows


def binary_metrics(rows: list[dict[str, Any]], prediction_key: str, positive_prediction: str) -> dict[str, Any]:
    evaluable = [
        row
        for row in rows
        if row["binary_origin_label"] in {"tumor_derived", "chip_or_hematopoietic"}
    ]
    tp = sum(
        1
        for row in evaluable
        if row[prediction_key] == positive_prediction and row["binary_origin_label"] == "chip_or_hematopoietic"
    )
    fp = sum(
        1
        for row in evaluable
        if row[prediction_key] == positive_prediction and row["binary_origin_label"] == "tumor_derived"
    )
    fn = sum(
        1
        for row in evaluable
        if row[prediction_key] != positive_prediction and row["binary_origin_label"] == "chip_or_hematopoietic"
    )
    tn = sum(
        1
        for row in evaluable
        if row[prediction_key] != positive_prediction and row["binary_origin_label"] == "tumor_derived"
    )
    precision = tp / (tp + fp) if tp + fp else None
    recall = tp / (tp + fn) if tp + fn else None
    specificity = tn / (tn + fp) if tn + fp else None
    accuracy = (tp + tn) / len(evaluable) if evaluable else None
    return {
        "evaluable_rows": len(evaluable),
        "tp_chip": tp,
        "fp_chip": fp,
        "fn_chip": fn,
        "tn_tumor": tn,
        "chip_precision": precision,
        "chip_recall": recall,
        "tumor_specificity": specificity,
        "accuracy": accuracy,
    }


def build_summary(rows: list[dict[str, Any]], cf_df: pd.DataFrame, wbc_df: pd.DataFrame, tumor_df: pd.DataFrame) -> dict[str, Any]:
    label_counts = Counter(row["source_trace_gold_label"] for row in rows)
    binary_counts = Counter(row["binary_origin_label"] for row in rows)
    tissue_counts = Counter(row["patient_id"].split("-")[0] for row in rows)
    gene_counts = Counter(row["gene"] for row in rows)
    chip_gene_metrics = binary_metrics(rows, "baseline_chip_gene_label", "chip_or_hematopoietic")
    vaf_metrics = binary_metrics(rows, "baseline_vaf_label", "low_vaf_possible_noise_or_chip")
    return {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "benchmark_id": "razavi_2019_mrd_source_trace_skill",
        "paper_id": "KWRAW-000862",
        "source_xlsx": str(DEFAULT_INPUT),
        "task": "Given a cfDNA caller variant, trace whether it is tumor-concordant, WBC/CHIP-concordant, shared, or unresolved using matched WBC and tumor tissue tables.",
        "cfDNA_rows": len(cf_df),
        "wbc_rows": len(wbc_df),
        "tumor_rows": len(tumor_df),
        "benchmark_rows": len(rows),
        "unique_patients": len({row["patient_id"] for row in rows}),
        "unique_genes": len({row["gene"] for row in rows if row["gene"]}),
        "source_trace_gold_label_counts": dict(label_counts),
        "binary_origin_label_counts": dict(binary_counts),
        "top_genes": dict(gene_counts.most_common(20)),
        "chip_gene_overlap_rows": sum(int(row["chip_gene_overlap"]) for row in rows),
        "low_alt_support_rows": sum(int(row["low_alt_support"]) for row in rows),
        "baseline_metrics": {
            "chip_gene_blacklist_for_chip_origin": chip_gene_metrics,
            "low_vaf_rule_for_chip_or_noise": vaf_metrics,
        },
        "intended_experiment": {
            "primary_question": "Do MRD-specific source-tracing skills resolve cfDNA caller variants better than generic VCF QC, CHIP gene blacklist, or VAF-only rules?",
            "baselines": ["generic VCF QC fields only", "CHIP gene blacklist", "VAF-only rule", "LLM-only", "MRD source-tracing skill"],
            "metrics": [
                "source_trace_label_accuracy",
                "CHIP/hematopoietic recall",
                "tumor-derived precision",
                "unresolved/shared flag rate",
                "evidence_traceability_rate",
                "unsupported_claim_rate",
            ],
        },
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def fmt_float(value: Any) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.3f}"
    return str(value)


def write_md(path: Path, summary: dict[str, Any], out_csv: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    metrics = summary["baseline_metrics"]["chip_gene_blacklist_for_chip_origin"]
    vaf_metrics = summary["baseline_metrics"]["low_vaf_rule_for_chip_or_noise"]
    lines = [
        "# Razavi 2019 MRD Source-Trace Skill Benchmark",
        "",
        f"- Generated: {summary['generated']}",
        f"- Benchmark ID: `{summary['benchmark_id']}`",
        f"- Source: `{summary['source_xlsx']}`",
        f"- Output CSV: `{out_csv}`",
        "",
        "## Task",
        "",
        summary["task"],
        "",
        "## Data Shape",
        "",
        f"- cfDNA rows: {summary['cfDNA_rows']:,}",
        f"- WBC rows: {summary['wbc_rows']:,}",
        f"- Tumor tissue rows: {summary['tumor_rows']:,}",
        f"- Benchmark cfDNA rows: {summary['benchmark_rows']:,}",
        f"- Unique patients: {summary['unique_patients']:,}",
        f"- Unique genes: {summary['unique_genes']:,}",
        "",
        "## Gold Labels",
        "",
    ]
    for label, count in summary["source_trace_gold_label_counts"].items():
        lines.append(f"- `{label}`: {count:,}")
    lines.extend(
        [
            "",
            "## Baseline Gap Check",
            "",
            "CHIP gene blacklist as a binary CHIP-origin baseline:",
            f"- Evaluable tumor-vs-WBC rows: {metrics['evaluable_rows']:,}",
            f"- CHIP precision: {fmt_float(metrics['chip_precision'])}",
            f"- CHIP recall: {fmt_float(metrics['chip_recall'])}",
            f"- Tumor specificity: {fmt_float(metrics['tumor_specificity'])}",
            f"- Accuracy: {fmt_float(metrics['accuracy'])}",
            "",
            "Low-VAF rule as a CHIP/noise proxy:",
            f"- CHIP precision: {fmt_float(vaf_metrics['chip_precision'])}",
            f"- CHIP recall: {fmt_float(vaf_metrics['chip_recall'])}",
            f"- Tumor specificity: {fmt_float(vaf_metrics['tumor_specificity'])}",
            f"- Accuracy: {fmt_float(vaf_metrics['accuracy'])}",
            "",
            "## Interpretation",
            "",
            "This benchmark directly tests the paper's central claim: generic caller or VCF QC fields are not enough for MRD source tracing. "
            "The same cfDNA variant must be adjudicated against matched WBC and tumor tissue context before it can be treated as tumor-derived MRD evidence, CHIP/hematopoietic signal, or unresolved.",
            "",
            "## Next Experiment Step",
            "",
            "After review, run MRD-KGQC/source-tracing skill predictions on this CSV and compare against the matched WBC/tumor-derived gold labels.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    input_xlsx = Path(args.input_xlsx)
    out_csv = Path(args.out_csv)
    out_json = Path(args.out_summary_json)
    out_md = Path(args.out_md)

    cf_df = load_sheet(input_xlsx, "Supplementary Table 11")
    wbc_df = load_sheet(input_xlsx, "Supplementary Table 12")
    tumor_df = load_sheet(input_xlsx, "Supplementary Table 13")
    rows = build_rows(cf_df, wbc_df, tumor_df)
    summary = build_summary(rows, cf_df, wbc_df, tumor_df)
    summary["source_xlsx"] = str(input_xlsx)
    write_csv(out_csv, rows)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    write_md(out_md, summary, out_csv)
    print(json.dumps({
        "benchmark_rows": len(rows),
        "unique_patients": summary["unique_patients"],
        "source_trace_gold_label_counts": summary["source_trace_gold_label_counts"],
        "out_csv": str(out_csv),
        "out_summary_json": str(out_json),
        "out_md": str(out_md),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
