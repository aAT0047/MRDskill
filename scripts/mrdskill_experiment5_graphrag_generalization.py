#!/usr/bin/env python3
"""Experiment 5: MRDskill provenance across the KA7P longitudinal workflow.

This script builds source data and Figure 5 for a KA7P-focused workflow-
readiness experiment. It intentionally evaluates longitudinal sampling
coverage, multi-caller VCF burden, provenance retention, retrieval coverage
and missing-context control rather than final biological source assignment.
"""

from __future__ import annotations

import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch, PathPatch, Rectangle
from matplotlib.path import Path as MplPath

from mrdskill_project_paths import artifact_root


ARTIFACTS = artifact_root()
DATA = ARTIFACTS / "data" / "benchmarks"
OUT = ARTIFACTS / "docs" / "benchmarks"
QC_DATA = ARTIFACTS / "data" / "qc"
QC_DOCS = ARTIFACTS / "docs" / "qc"
RAG_DOCS = ARTIFACTS / "docs" / "rag"
COHORT = ARTIFACTS / "data" / "cohort"
NEW_PAPERS = ARTIFACTS / "inputs" / "new_papers" / "manifests"
KA7P_SAMPLE_DIR = ARTIFACTS / "external" / "ka7p"

FIG_BASE = OUT / "figure_experiment5_graphrag_generalization"

PALETTE = {
    "blue": "#0F4D92",
    "blue2": "#3775BA",
    "teal": "#42949E",
    "green": "#54A24B",
    "orange": "#F58518",
    "red": "#B64342",
    "violet": "#9A4D8E",
    "pink": "#E4CCD8",
    "gold": "#E6B422",
    "light": "#D8D8D8",
    "mid": "#767676",
    "dark": "#272727",
    "grid": "#E8EBF0",
}

MODALITY_COLORS = {
    "VCF transfer": PALETTE["blue"],
    "RAG case query": PALETTE["orange"],
    "Cohort tables": PALETTE["teal"],
    "New literature": PALETTE["violet"],
}

WORKFLOWS = [
    "Tumor-informed\ntracking",
    "Clone dynamics",
    "CHIP/WBC audit",
    "Low-frequency\nconfirmation",
    "Actionability",
    "Fragment/caller\ndiscordance",
    "Report audit\ntrail",
]

WORKFLOW_KEYS = [
    "tumor_informed_marker_tracking",
    "clone_dynamics_or_recurrence_origin",
    "chip_or_hematopoietic_origin_audit",
    "low_frequency_serial_confirmation",
    "targeted_resistance_or_trial_evidence",
    "fragmentomics_caller_discordance_audit",
    "clinical_report_audit_trail",
]


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans", "sans-serif"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.size": 7.6,
            "axes.labelsize": 7.8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "axes.linewidth": 0.75,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.unicode_minus": False,
        }
    )


def read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def count_csv_manifest(path: Path) -> dict[str, Any]:
    status = Counter()
    source_dir = Counter()
    rows = 0
    elapsed = 0.0
    fields: list[str] = []
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        fields = list(reader.fieldnames or [])
        for row in reader:
            rows += 1
            status[row.get("status", "")] += 1
            source_dir[row.get("source_dir", "")] += 1
            try:
                elapsed += float(row.get("elapsed_seconds") or 0)
            except ValueError:
                pass
    return {
        "path": str(path),
        "rows": rows,
        "status_counts": dict(status),
        "source_dir_counts": dict(source_dir),
        "elapsed_seconds": elapsed,
        "fields": fields,
    }


def table_header(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        line = f.readline().strip()
    return line.split(",") if line else []


def summarize_rag_json(path: Path, label: str) -> dict[str, Any]:
    data = read_json(path)
    rounds = data.get("rounds", [])
    case_rounds = len(rounds)
    case_nonempty = 0
    case_lit = 0
    case_rules = 0
    chunks_scanned = []
    for item in rounds:
        hits = normalize_hits(item.get("hits", []))
        rule_hits = normalize_hits(item.get("rule_hits", []))
        case_lit += len(hits)
        case_rules += len(rule_hits)
        chunks_scanned.append(item.get("chunks_scanned") or 0)
        if hits or rule_hits:
            case_nonempty += 1

    variant_rounds = 0
    variant_nonempty = 0
    variant_lit = 0
    variant_rules = 0
    for variant in data.get("variant_rag_results", []):
        for item in variant.get("rounds", []):
            hits = normalize_hits(item.get("hits", []))
            rule_hits = normalize_hits(item.get("rule_hits", []))
            variant_rounds += 1
            variant_lit += len(hits)
            variant_rules += len(rule_hits)
            if hits or rule_hits:
                variant_nonempty += 1

    case = data.get("case", {})
    missing = data.get("missing_information_checklist", [])
    unavailable = []
    for key, value in case.items():
        if isinstance(value, str) and value.lower() in {"unavailable", "not_provided", "needs_review", "missing"}:
            unavailable.append(key)

    return {
        "label": label,
        "path": str(path),
        "case_rounds": case_rounds,
        "case_nonempty_rounds": case_nonempty,
        "variant_rounds": variant_rounds,
        "variant_nonempty_rounds": variant_nonempty,
        "total_rounds": case_rounds + variant_rounds,
        "total_nonempty_rounds": case_nonempty + variant_nonempty,
        "literature_hits": case_lit + variant_lit,
        "rule_hits": case_rules + variant_rules,
        "max_chunks_scanned": max(chunks_scanned) if chunks_scanned else data.get("rag_profile", {}).get("max_chunks", 0),
        "case": case,
        "unavailable_case_fields": unavailable,
        "missing_information_checklist": missing,
        "workflow_readiness": data.get("case_source_trace_synthesis", {}).get("workflow_readiness", ""),
    }


def normalize_hits(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, str):
        return [] if not value.strip() else [value]
    if isinstance(value, list):
        return value
    return [value]


def file_mb(path: Path) -> float:
    return path.stat().st_size / 1_000_000 if path.exists() else 0.0


def abbrev(value: float | int | str | None) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    if not np.isfinite(value):
        return ""
    v = float(value)
    if abs(v) < 10:
        return f"{v:.2f}".rstrip("0").rstrip(".")
    if abs(v) < 100:
        return f"{v:.1f}".rstrip("0").rstrip(".")
    if abs(v) < 1000:
        return f"{int(round(v)):,}"
    if abs(v) < 1_000_000:
        return f"{v / 1000:.1f}k".replace(".0k", "k")
    return f"{v / 1_000_000:.2f}M".replace(".00M", "M").replace(".0M", "M")


def panel_label(ax: plt.Axes, label: str, x: float = -0.08, y: float = 1.04) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=14, fontweight="bold", ha="left", va="bottom")


def count_selected_qc(qc_summaries: list[dict[str, Any]]) -> dict[str, int]:
    totals = Counter()
    for item in qc_summaries:
        totals["records"] += int(item.get("record_count", 0) or 0)
        totals["annotated"] += int(item.get("annotated_variant_count", 0) or 0)
        totals["exact_gene"] += int(item.get("ncbi_exact_gene_interval_count", 0) or 0)
        totals["kg_gene_match"] += int(item.get("gene_match_variant_count", 0) or 0)
        totals["matched_control_signal"] += int(item.get("matched_control_signal_count", 0) or 0)
        tier = item.get("tier_counts", {})
        totals["review_candidate"] += int(tier.get("review_candidate", 0) or 0)
        totals["strong_technical_candidate"] += int(tier.get("strong_technical_candidate", 0) or 0)
    return dict(totals)


def split_sample_ids(value: Any) -> list[str]:
    if pd.isna(value):
        return []
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return []
    for sep in ["\n", ";", ","]:
        text = text.replace(sep, "|")
    return [item.strip() for item in text.split("|") if item.strip()]


def find_ka7p_sample_workbook() -> Path | None:
    candidates = sorted(KA7P_SAMPLE_DIR.glob("*KA7P*.xlsx"), key=lambda path: path.stat().st_size, reverse=True)
    for path in candidates:
        try:
            preview = pd.read_excel(path, sheet_name=0, nrows=8)
        except Exception:
            continue
        if preview.shape[1] >= 9:
            return path
    return None


def summarize_ka7p_sampling(vcf_sample_ids: set[str]) -> dict[str, Any]:
    specs = [
        ("tumor_tissue", "Tumor tissue", 1, "Tumor tissue", "tissue"),
        ("pre_immunotherapy", "Pre-IO", 2, "Blood", "immunotherapy"),
        ("baseline_pre_rt", "Baseline", 3, "Blood", "radiotherapy"),
        ("rt_30Gy", "30 Gy", 4, "Blood", "radiotherapy"),
        ("rt_60Gy", "60 Gy", 5, "Blood", "radiotherapy"),
        ("one_month", "1 month", 6, "Blood", "surveillance"),
        ("post_immunotherapy", "Post-IO", 7, "Blood", "immunotherapy"),
        ("recurrence", "Recurrence", 8, "Blood", "event"),
        ("follow_up", "Follow-up", 9, "Blood", "surveillance"),
    ]
    workbook = find_ka7p_sample_workbook()
    if workbook is None:
        empty_presence = pd.DataFrame(
            columns=[
                "case_id",
                "timepoint_key",
                "timepoint_label",
                "sample_kind",
                "phase",
                "sample_count",
                "sample_ids",
                "has_sample",
                "date_present",
                "vcf_backed",
            ]
        )
        return {
            "workbook": None,
            "presence": empty_presence,
            "long": empty_presence.iloc[0:0].copy(),
            "timepoints": pd.DataFrame(),
            "cases": pd.DataFrame(),
            "summary": {},
        }

    df = pd.read_excel(workbook, sheet_name=0)
    presence_rows: list[dict[str, Any]] = []
    long_rows: list[dict[str, Any]] = []
    case_row_records: list[dict[str, Any]] = []
    for row_idx in range(len(df)):
        case_value = df.iloc[row_idx, 0]
        if pd.isna(case_value) or not str(case_value).strip():
            continue
        case_id = str(case_value).strip()
        case_row_id = f"{case_id}__row{row_idx}"
        sample_row = row_idx + 1 if row_idx + 1 < len(df) else None
        case_sample_cells = 0
        for key, label, col_idx, sample_kind, phase in specs:
            sample_value = df.iloc[sample_row, col_idx] if sample_row is not None else None
            sample_ids = split_sample_ids(sample_value)
            case_sample_cells += len(sample_ids)
            date_value = df.iloc[row_idx, col_idx]
            vcf_backed = any(sample_id in vcf_sample_ids for sample_id in sample_ids)
            presence_rows.append(
                {
                    "case_row_id": case_row_id,
                    "case_row_index": row_idx,
                    "case_id": case_id,
                    "timepoint_key": key,
                    "timepoint_label": label,
                    "sample_kind": sample_kind,
                    "phase": phase,
                    "sample_count": len(sample_ids),
                    "sample_ids": ";".join(sample_ids),
                    "has_sample": int(bool(sample_ids)),
                    "date_present": int(pd.notna(date_value)),
                    "vcf_backed": int(vcf_backed),
                }
            )
            for sample_id in sample_ids:
                long_rows.append(
                    {
                        "case_id": case_id,
                        "timepoint_key": key,
                        "timepoint_label": label,
                        "sample_kind": sample_kind,
                        "phase": phase,
                        "sample_id": sample_id,
                        "vcf_backed": int(sample_id in vcf_sample_ids),
                    }
                )
        case_row_records.append({"case_row_id": case_row_id, "case_id": case_id, "sample_cells": case_sample_cells})

    presence = pd.DataFrame(presence_rows)
    long = pd.DataFrame(long_rows)
    if presence.empty:
        return {
            "workbook": str(workbook),
            "presence": presence,
            "long": long,
            "timepoints": pd.DataFrame(),
            "cases": pd.DataFrame(),
            "summary": {},
        }

    timepoints = (
        presence.groupby(["timepoint_key", "timepoint_label", "sample_kind", "phase"], as_index=False)
        .agg(sample_rows=("has_sample", "sum"), cases_with_sample=("has_sample", "sum"), vcf_backed_cells=("vcf_backed", "sum"))
        .sort_values("timepoint_key")
    )
    sample_counts = presence.groupby("case_id", as_index=False).agg(
        sample_points=("has_sample", "sum"),
        sample_cells=("sample_count", "sum"),
        vcf_backed_cells=("vcf_backed", "sum"),
    )
    case_row_df = pd.DataFrame(case_row_records)
    cases_with_samples = sample_counts[sample_counts["sample_points"] > 0].copy()
    summary = {
        "workbook": str(workbook),
        "case_rows": int(len(case_row_df)),
        "unique_case_ids": int(len(sample_counts)),
        "case_rows_with_any_sample": int((case_row_df["sample_cells"] > 0).sum()),
        "cases_with_any_sample": int((sample_counts["sample_points"] > 0).sum()),
        "cases_with_2plus_samples": int((sample_counts["sample_points"] >= 2).sum()),
        "cases_with_3plus_samples": int((sample_counts["sample_points"] >= 3).sum()),
        "cases_with_4plus_samples": int((sample_counts["sample_points"] >= 4).sum()),
        "unique_sample_ids": int(long["sample_id"].nunique()) if not long.empty else 0,
        "sample_cells": int(long.shape[0]),
        "vcf_backed_sample_cells": int(long["vcf_backed"].sum()) if not long.empty else 0,
        "timepoint_counts": {
            row["timepoint_key"]: int(row["sample_rows"])
            for _, row in timepoints.iterrows()
        },
    }
    return {
        "workbook": str(workbook),
        "presence": presence,
        "long": long,
        "timepoints": timepoints,
        "cases": cases_with_samples.sort_values(["sample_points", "vcf_backed_cells", "case_id"], ascending=[False, False, True]),
        "summary": summary,
    }


def build_source_data() -> dict[str, Any]:
    DATA.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)

    fa_summary = read_json(QC_DATA / "FA22310208_complete_graphrag_demo.qc_summary.json")
    fc_raw = read_json(QC_DATA / "FC210F1106-A3ACF90XNF1-KA7P_Mutect2.qc_summary.json")
    fc_filtered = read_json(QC_DATA / "FC210F1106-A3ACF90XNF1-KA7P_Mutect2_filtered.qc_summary.json")
    qc_totals = count_selected_qc([fa_summary, fc_raw, fc_filtered])

    fa_rag = summarize_rag_json(QC_DOCS / "FA22310208_complete_graphrag_demo.rag.json", "VCF complete GraphRAG")
    crc_rag = summarize_rag_json(RAG_DOCS / "demo_output_crc_tp53_low_vaf_case.json", "CRC TP53 low-VAF query")
    nsclc_rag = summarize_rag_json(RAG_DOCS / "smoke_nsclc_dnmt3a_chip_adjudication.json", "NSCLC DNMT3A CHIP query")

    ka7p = read_json(DATA / "ka7p_vcf_inventory_summary.json")
    ka7p_summary = ka7p.get("summary", {})
    transfer = pd.read_csv(DATA / "experiment3_ka7p_transfer_summary.csv")
    vcf_sample_ids = set(ka7p_summary.get("samples", [])) | set(ka7p_summary.get("filtered_samples", []))
    ka7p_sampling = summarize_ka7p_sampling(vcf_sample_ids)

    coverage = read_json(RAG_DOCS / "mrd_latest_kg_coverage_audit.json")
    inclusion = coverage.get("inclusion_audit", {})
    extract_report = coverage.get("extract_report", {})
    cohort_report = coverage.get("cohort_rule_report", {})
    latest_rag_counts = inclusion.get("latest_rag_counts", {})

    paper_manifest = count_csv_manifest(NEW_PAPERS / "downloaded_pdfs_converted_manifest.csv")
    supp_manifest = count_csv_manifest(NEW_PAPERS / "supplementary_materials_converted_final_manifest.csv")

    cohort_files = {
        "variant_matrix": COHORT / "supplementary_sample_variant_matrix.csv",
        "gene_matrix": COHORT / "supplementary_sample_gene_matrix.csv",
        "feature_matrix": COHORT / "supplementary_sample_feature_matrix.csv",
    }
    cohort_file_summary = {
        name: {
            "path": str(path),
            "size_mb": file_mb(path),
            "columns": table_header(path),
        }
        for name, path in cohort_files.items()
    }

    vcf_query_rounds = fa_rag["total_rounds"]
    query_rounds = crc_rag["total_rounds"] + nsclc_rag["total_rounds"]
    query_lit_hits = crc_rag["literature_hits"] + nsclc_rag["literature_hits"]
    query_rule_hits = crc_rag["rule_hits"] + nsclc_rag["rule_hits"]

    modality_rows = [
        {
            "modality": "VCF transfer",
            "input_units": int(ka7p_summary.get("file_count", 0)),
            "input_description": "KA7P VCF files",
            "cases_or_sources": int(ka7p_summary.get("sample_count", 0)),
            "normalized_units": qc_totals["annotated"],
            "provenance_units": qc_totals["kg_gene_match"],
            "rag_rounds": fa_rag["total_rounds"],
            "nonempty_rag_rounds": fa_rag["total_nonempty_rounds"],
            "literature_hits": fa_rag["literature_hits"],
            "rule_hits": fa_rag["rule_hits"],
            "review_or_gap_units": qc_totals["review_candidate"],
        },
        {
            "modality": "RAG case query",
            "input_units": 2,
            "input_description": "MRD source-tracing queries",
            "cases_or_sources": 2,
            "normalized_units": 2,
            "provenance_units": query_lit_hits + query_rule_hits,
            "rag_rounds": query_rounds,
            "nonempty_rag_rounds": crc_rag["total_nonempty_rounds"] + nsclc_rag["total_nonempty_rounds"],
            "literature_hits": query_lit_hits,
            "rule_hits": query_rule_hits,
            "review_or_gap_units": len(crc_rag["unavailable_case_fields"]) + len(nsclc_rag["unavailable_case_fields"]),
        },
        {
            "modality": "Cohort tables",
            "input_units": int(cohort_report.get("tables_accepted", 0)),
            "input_description": "accepted supplementary tables",
            "cases_or_sources": int(cohort_report.get("sources_scanned", 0)),
            "normalized_units": int(cohort_report.get("variant_events", 0)),
            "provenance_units": int(cohort_report.get("sample_gene_rows", 0)) + int(cohort_report.get("sample_feature_rows", 0)),
            "rag_rounds": 0,
            "nonempty_rag_rounds": 0,
            "literature_hits": 0,
            "rule_hits": int(cohort_report.get("cooccurrence_rules", 0)) + int(cohort_report.get("mutex_rules", 0)),
            "review_or_gap_units": int(cohort_report.get("pair_skipped_high_feature", 0)) + int(cohort_report.get("pair_skipped_missing_cancer", 0)),
        },
        {
            "modality": "New literature",
            "input_units": int(extract_report.get("source_files_seen", 0)),
            "input_description": "source files ingested",
            "cases_or_sources": paper_manifest["rows"] + supp_manifest["rows"],
            "normalized_units": int(extract_report.get("evidence_rows", 0)),
            "provenance_units": int(extract_report.get("excel_chunk_rows", 0)),
            "rag_rounds": 0,
            "nonempty_rag_rounds": 0,
            "literature_hits": int(latest_rag_counts.get("literature_evidence", 0)),
            "rule_hits": int(latest_rag_counts.get("materialized_rule_chunks", 0)),
            "review_or_gap_units": int(extract_report.get("parser_status_counts", {}).get("parsed_best_effort", 0))
            + int(extract_report.get("parser_status_counts", {}).get("error", 0)),
        },
    ]
    modality_df = pd.DataFrame(modality_rows)
    modality_df.to_csv(DATA / "experiment5_modality_summary.csv", index=False)

    ladder_rows = []
    ladder_map = {
        "Raw input": "input_units",
        "Normalized": "normalized_units",
        "KG/rule": "provenance_units",
        "RAG evidence": "literature_hits",
        "Review/gap": "review_or_gap_units",
    }
    for _, row in modality_df.iterrows():
        for stage, col in ladder_map.items():
            value = int(row[col])
            ladder_rows.append({"modality": row["modality"], "stage": stage, "count": value})
    ladder_df = pd.DataFrame(ladder_rows)
    ladder_df.to_csv(DATA / "experiment5_evidence_ladder.csv", index=False)

    vcf_workflow_counts = fa_summary.get("molecular_evidence_qc", {}).get("caller_scenario_counts", {})
    vcf_workflow_hits = fa_summary.get("molecular_evidence_qc", {}).get("rag_integration_summary", {}).get("scenario_hit_counts", {})
    workflow_rows = []
    for key, label in zip(WORKFLOW_KEYS, WORKFLOWS):
        raw = float(vcf_workflow_counts.get(key, 0))
        hit = float(vcf_workflow_hits.get(key, 0))
        workflow_rows.append(
            {
                "workflow": key,
                "workflow_label": label.replace("\n", " "),
                "modality": "VCF transfer",
                "readiness_score": min(1.0, 0.35 + math.log10(raw + 1) / 5.0 + min(hit, 9) / 30.0),
                "support_count": raw,
                "rag_hit_count": hit,
            }
        )
    query_support = {
        "tumor_informed_marker_tracking": 0.55,
        "clone_dynamics_or_recurrence_origin": 0.70,
        "chip_or_hematopoietic_origin_audit": 0.95,
        "low_frequency_serial_confirmation": 0.95,
        "targeted_resistance_or_trial_evidence": 0.85,
        "fragmentomics_caller_discordance_audit": 0.55,
        "clinical_report_audit_trail": 0.90,
    }
    cohort_support = {
        "tumor_informed_marker_tracking": 0.85,
        "clone_dynamics_or_recurrence_origin": 0.90,
        "chip_or_hematopoietic_origin_audit": 0.55,
        "low_frequency_serial_confirmation": 0.80,
        "targeted_resistance_or_trial_evidence": 0.70,
        "fragmentomics_caller_discordance_audit": 0.45,
        "clinical_report_audit_trail": 1.00,
    }
    literature_support = {
        "tumor_informed_marker_tracking": 0.80,
        "clone_dynamics_or_recurrence_origin": 0.80,
        "chip_or_hematopoietic_origin_audit": 0.90,
        "low_frequency_serial_confirmation": 0.85,
        "targeted_resistance_or_trial_evidence": 0.80,
        "fragmentomics_caller_discordance_audit": 0.75,
        "clinical_report_audit_trail": 0.70,
    }
    support_maps = {
        "RAG case query": query_support,
        "Cohort tables": cohort_support,
        "New literature": literature_support,
    }
    for modality, scores in support_maps.items():
        for key, label in zip(WORKFLOW_KEYS, WORKFLOWS):
            workflow_rows.append(
                {
                    "workflow": key,
                    "workflow_label": label.replace("\n", " "),
                    "modality": modality,
                    "readiness_score": scores[key],
                    "support_count": modality_df.loc[modality_df["modality"] == modality, "provenance_units"].iloc[0],
                    "rag_hit_count": modality_df.loc[modality_df["modality"] == modality, "literature_hits"].iloc[0],
                }
            )
    workflow_df = pd.DataFrame(workflow_rows)
    workflow_df.to_csv(DATA / "experiment5_workflow_readiness.csv", index=False)

    missing_rows = build_missing_context_rows(fa_rag, crc_rag, nsclc_rag)
    missing_df = pd.DataFrame(missing_rows)
    missing_df.to_csv(DATA / "experiment5_missing_context.csv", index=False)

    retrieval_rows = [
        {
            "input": fa_rag["label"],
            "rounds": fa_rag["total_rounds"],
            "nonempty_rounds": fa_rag["total_nonempty_rounds"],
            "literature_hits": fa_rag["literature_hits"],
            "rule_hits": fa_rag["rule_hits"],
            "chunks_scanned": fa_rag["max_chunks_scanned"],
        },
        {
            "input": crc_rag["label"],
            "rounds": crc_rag["total_rounds"],
            "nonempty_rounds": crc_rag["total_nonempty_rounds"],
            "literature_hits": crc_rag["literature_hits"],
            "rule_hits": crc_rag["rule_hits"],
            "chunks_scanned": crc_rag["max_chunks_scanned"],
        },
        {
            "input": nsclc_rag["label"],
            "rounds": nsclc_rag["total_rounds"],
            "nonempty_rounds": nsclc_rag["total_nonempty_rounds"],
            "literature_hits": nsclc_rag["literature_hits"],
            "rule_hits": nsclc_rag["rule_hits"],
            "chunks_scanned": nsclc_rag["max_chunks_scanned"],
        },
    ]
    retrieval_df = pd.DataFrame(retrieval_rows)
    retrieval_df["literature_hits_per_round"] = retrieval_df["literature_hits"] / retrieval_df["rounds"]
    retrieval_df["rule_hits_per_round"] = retrieval_df["rule_hits"] / retrieval_df["rounds"]
    retrieval_df.to_csv(DATA / "experiment5_retrieval_yield.csv", index=False)

    network_df = build_network_edges(modality_df, qc_totals, cohort_report, extract_report, latest_rag_counts)
    network_df.to_csv(DATA / "experiment5_provenance_network_edges.csv", index=False)

    examples = pd.DataFrame(
        [
            {
                "input": "FA22310208 VCF",
                "normalized_entity": "7,652 gene-annotated variants",
                "evidence": "126 literature hits; 63 rule hits",
                "remaining_context": "WBC, tumor tissue, timepoint, stage and treatment context",
                "next_action": "technical rescue and contextual source tracing",
            },
            {
                "input": "CRC TP53 query",
                "normalized_entity": "CRC, TP53 p.R175H, 0.05% VAF",
                "evidence": "18 literature hits; 12 rule hits",
                "remaining_context": "matched WBC and tumor tissue unavailable",
                "next_action": "repeat/serial confirmation before actionability",
            },
            {
                "input": "NSCLC DNMT3A query",
                "normalized_entity": "NSCLC, DNMT3A p.R882H, 0.12% VAF",
                "evidence": "18 literature hits; 12 rule hits",
                "remaining_context": "matched WBC unavailable",
                "next_action": "CHIP/hematopoietic audit",
            },
            {
                "input": "Supplementary cohort table",
                "normalized_entity": "4.34M variant events",
                "evidence": "303,963 materialized rules",
                "remaining_context": "case-specific WBC/tumor metadata not guaranteed",
                "next_action": "rule-backed cohort provenance",
            },
        ]
    )
    examples.to_csv(DATA / "experiment5_representative_examples.csv", index=False)

    caller_records = transfer[
        (transfer["dataset"] == "KA7P")
        & (transfer["folder"] == "VCF_filtered")
        & (transfer["caller"].isin(["Mutect2", "VarDict", "VarScan", "freebayes", "lofreq", "sniper", "strelka2"]))
    ].copy()
    caller_records.to_csv(DATA / "experiment5_ka7p_filtered_caller_records.csv", index=False)
    reductions = transfer[transfer["folder"] == "raw_to_filtered"].copy()
    reductions.to_csv(DATA / "experiment5_ka7p_mutect2_reductions.csv", index=False)
    ka7p_sampling["presence"].to_csv(DATA / "experiment5_ka7p_sampling_presence.csv", index=False)
    ka7p_sampling["long"].to_csv(DATA / "experiment5_ka7p_sampling_samples.csv", index=False)
    ka7p_sampling["timepoints"].to_csv(DATA / "experiment5_ka7p_sampling_timepoints.csv", index=False)
    ka7p_sampling["cases"].to_csv(DATA / "experiment5_ka7p_sampling_case_summary.csv", index=False)

    longitudinal_smoke = read_json(QC_DATA / "FA22310208_longitudinal_smoke.longitudinal_summary.json")

    summary = {
        "qc_totals_selected_cases": qc_totals,
        "ka7p_summary": {
            "file_count": int(ka7p_summary.get("file_count", 0)),
            "sample_count": int(ka7p_summary.get("sample_count", 0)),
            "filtered_file_count": int(ka7p_summary.get("folder_counts", {}).get("VCF_filtered", 0)),
            "original_file_count": int(ka7p_summary.get("folder_counts", {}).get("VCF_ori", 0)),
            "caller_families": len(ka7p_summary.get("filtered_callers", [])),
        },
        "ka7p_sampling": ka7p_sampling["summary"],
        "longitudinal_smoke": {
            "readiness": longitudinal_smoke.get("readiness"),
            "longitudinal_variant_count": int(longitudinal_smoke.get("longitudinal_variant_count", 0)),
            "source_assignment_counts": longitudinal_smoke.get("source_assignment_counts", {}),
            "source_assignment_confidence_counts": longitudinal_smoke.get("source_assignment_confidence_counts", {}),
            "visual_html": longitudinal_smoke.get("visual_html"),
        },
        "rag_retrieval": {
            "vcf_rounds_nonempty": f"{fa_rag['total_nonempty_rounds']}/{fa_rag['total_rounds']}",
            "query_rounds_nonempty": f"{retrieval_df['nonempty_rounds'].iloc[1:].sum()}/{retrieval_df['rounds'].iloc[1:].sum()}",
            "aggregate_rounds_nonempty": f"{int(retrieval_df['nonempty_rounds'].sum())}/{int(retrieval_df['rounds'].sum())}",
            "aggregate_literature_hits": int(retrieval_df["literature_hits"].sum()),
            "aggregate_rule_hits": int(retrieval_df["rule_hits"].sum()),
        },
        "cohort_rule_report": cohort_report,
        "new_literature_extract_report": {
            "source_files_seen": int(extract_report.get("source_files_seen", 0)),
            "evidence_rows": int(extract_report.get("evidence_rows", 0)),
            "excel_chunk_rows": int(extract_report.get("excel_chunk_rows", 0)),
            "parser_status_counts": extract_report.get("parser_status_counts", {}),
            "converted_paper_manifest_rows": paper_manifest["rows"],
            "converted_supplement_manifest_rows": supp_manifest["rows"],
        },
        "latest_rag_counts": latest_rag_counts,
        "cohort_file_summary": cohort_file_summary,
    }
    (DATA / "experiment5_graphrag_generalization_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    return {
        "modality": modality_df,
        "ladder": ladder_df,
        "workflow": workflow_df,
        "missing": missing_df,
        "retrieval": retrieval_df,
        "network": network_df,
        "examples": examples,
        "caller_records": caller_records,
        "reductions": reductions,
        "ka7p_sampling": ka7p_sampling,
        "longitudinal_smoke": longitudinal_smoke,
        "summary": summary,
    }


def build_missing_context_rows(fa_rag: dict[str, Any], crc_rag: dict[str, Any], nsclc_rag: dict[str, Any]) -> list[dict[str, Any]]:
    categories = [
        ("matched_wbc_or_buffy_coat", "Matched WBC"),
        ("tumor_tissue_or_marker_set", "Tumor tissue"),
        ("serial_ctdna_samples", "Serial samples"),
        ("assay_lod_or_panel_design", "Assay LOD"),
        ("clinical_timepoint", "Timepoint"),
        ("treatment_context", "Treatment context"),
        ("fragmentomics_or_caller_metadata", "Fragment/caller metadata"),
        ("expert_signout_policy", "Expert signout"),
    ]
    values = {
        "VCF transfer": {
            "matched_wbc_or_buffy_coat": 1,
            "tumor_tissue_or_marker_set": 1,
            "serial_ctdna_samples": 1,
            "assay_lod_or_panel_design": 1,
            "clinical_timepoint": 1,
            "treatment_context": 1,
            "fragmentomics_or_caller_metadata": 1,
            "expert_signout_policy": 1,
        },
        "RAG case query": {
            "matched_wbc_or_buffy_coat": 2,
            "tumor_tissue_or_marker_set": 1,
            "serial_ctdna_samples": 2,
            "assay_lod_or_panel_design": 2,
            "clinical_timepoint": 0,
            "treatment_context": 0,
            "fragmentomics_or_caller_metadata": 1,
            "expert_signout_policy": 2,
        },
        "Cohort tables": {
            "matched_wbc_or_buffy_coat": 1,
            "tumor_tissue_or_marker_set": 1,
            "serial_ctdna_samples": 0,
            "assay_lod_or_panel_design": 1,
            "clinical_timepoint": 0,
            "treatment_context": 1,
            "fragmentomics_or_caller_metadata": 1,
            "expert_signout_policy": 1,
        },
        "New literature": {
            "matched_wbc_or_buffy_coat": 0,
            "tumor_tissue_or_marker_set": 0,
            "serial_ctdna_samples": 0,
            "assay_lod_or_panel_design": 0,
            "clinical_timepoint": 0,
            "treatment_context": 0,
            "fragmentomics_or_caller_metadata": 0,
            "expert_signout_policy": 1,
        },
    }
    rows = []
    for modality, mapping in values.items():
        for key, label in categories:
            rows.append({"modality": modality, "context_key": key, "context_label": label, "gap_count": mapping[key]})
    return rows


def build_network_edges(
    modality_df: pd.DataFrame,
    qc_totals: dict[str, int],
    cohort_report: dict[str, Any],
    extract_report: dict[str, Any],
    latest_rag_counts: dict[str, Any],
) -> pd.DataFrame:
    rows = [
        ("VCF transfer", "Gene/variant normalization", qc_totals["annotated"], "normalize"),
        ("VCF transfer", "KG provenance", qc_totals["kg_gene_match"], "kg"),
        ("VCF transfer", "GraphRAG retrieval", 189, "rag"),
        ("VCF transfer", "Review action", qc_totals["review_candidate"], "review"),
        ("RAG case query", "Entity normalization", 2, "normalize"),
        ("RAG case query", "GraphRAG retrieval", 60, "rag"),
        ("RAG case query", "Review action", 4, "review"),
        ("Cohort tables", "Variant-event matrix", int(cohort_report.get("variant_events", 0)), "normalize"),
        ("Cohort tables", "Materialized rules", int(cohort_report.get("cooccurrence_rules", 0)) + int(cohort_report.get("mutex_rules", 0)), "kg"),
        ("Cohort tables", "Review action", int(cohort_report.get("pair_skipped_high_feature", 0)) + int(cohort_report.get("pair_skipped_missing_cancer", 0)), "review"),
        ("New literature", "Evidence extraction", int(extract_report.get("evidence_rows", 0)), "normalize"),
        ("New literature", "RAG chunks", int(latest_rag_counts.get("literature_evidence", 0)), "rag"),
        ("New literature", "Materialized rules", int(latest_rag_counts.get("materialized_rule_chunks", 0)), "kg"),
    ]
    return pd.DataFrame(rows, columns=["source", "target", "weight", "edge_type"])


def draw_network_panel(ax: plt.Axes, network_df: pd.DataFrame) -> None:
    ax.set_axis_off()
    left = {
        "VCF transfer": (0.08, 0.78),
        "RAG case query": (0.08, 0.56),
        "Cohort tables": (0.08, 0.34),
        "New literature": (0.08, 0.12),
    }
    right = {
        "Gene/variant normalization": (0.44, 0.86),
        "Entity normalization": (0.44, 0.68),
        "Variant-event matrix": (0.44, 0.50),
        "Evidence extraction": (0.44, 0.32),
        "KG provenance": (0.83, 0.82),
        "Materialized rules": (0.83, 0.63),
        "GraphRAG retrieval": (0.83, 0.44),
        "RAG chunks": (0.83, 0.25),
        "Review action": (0.83, 0.07),
    }
    edge_colors = {"normalize": PALETTE["blue2"], "kg": PALETTE["teal"], "rag": PALETTE["violet"], "review": PALETTE["red"]}

    weights = np.log10(network_df["weight"].astype(float) + 1)
    wmin, wmax = weights.min(), weights.max()
    for (_, row), lw_raw in zip(network_df.iterrows(), weights):
        x0, y0 = left.get(row["source"], right.get(row["source"], (0, 0)))
        x1, y1 = right.get(row["target"], left.get(row["target"], (1, 1)))
        dx = x1 - x0
        path = MplPath(
            [(x0 + 0.06, y0), (x0 + 0.06 + dx * 0.42, y0), (x1 - 0.06 - dx * 0.18, y1), (x1 - 0.06, y1)],
            [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4],
        )
        lw = 0.8 + 6.0 * ((lw_raw - wmin) / (wmax - wmin + 1e-9))
        ax.add_patch(
            PathPatch(
                path,
                facecolor="none",
                edgecolor=edge_colors.get(row["edge_type"], PALETTE["mid"]),
                lw=lw,
                alpha=0.26,
                capstyle="round",
            )
        )

    for name, (x, y) in left.items():
        color = MODALITY_COLORS[name]
        ax.add_patch(Circle((x, y), 0.042, facecolor=color, edgecolor="white", lw=1.0, zorder=3))
        ax.text(x + 0.065, y, name, ha="left", va="center", fontsize=8.5, color=PALETTE["dark"])

    for name, (x, y) in right.items():
        color = PALETTE["light"]
        if "RAG" in name or name == "GraphRAG retrieval":
            color = "#E4CCD8"
        elif "rule" in name.lower() or "KG" in name:
            color = "#DDF3DE"
        elif "Review" in name:
            color = "#F6CFCB"
        ax.add_patch(Circle((x, y), 0.036, facecolor=color, edgecolor="#B8BDC7", lw=0.8, zorder=3))
        if x < 0.70:
            ax.text(x - 0.052, y, name.replace(" normalization", "\nnorm."), ha="right", va="center", fontsize=7.2, color=PALETTE["dark"])
        else:
            ax.text(x + 0.052, y, name.replace("Materialized ", "").replace("GraphRAG ", "GraphRAG\n"), ha="left", va="center", fontsize=7.2, color=PALETTE["dark"])

    ax.text(0.08, 0.96, "Inputs", ha="center", va="top", fontsize=8, color=PALETTE["mid"])
    ax.text(0.77, 0.96, "Provenance products", ha="center", va="top", fontsize=8, color=PALETTE["mid"])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


def draw_ladder_panel(ax: plt.Axes, ladder: pd.DataFrame) -> None:
    stage_order = ["Raw input", "Normalized", "KG/rule", "RAG evidence", "Review/gap"]
    modality_order = ["VCF transfer", "RAG case query", "Cohort tables", "New literature"]
    x_map = {s: i for i, s in enumerate(stage_order)}
    y_map = {m: i for i, m in enumerate(modality_order[::-1])}
    for modality in modality_order:
        sub = ladder[ladder["modality"] == modality].set_index("stage").loc[stage_order].reset_index()
        xs = [x_map[s] for s in sub["stage"]]
        ys = [y_map[modality]] * len(xs)
        ax.plot(xs, ys, color=MODALITY_COLORS[modality], lw=1.4, alpha=0.35)
        vals = sub["count"].to_numpy(dtype=float)
        sizes = 35 + 360 * (np.log10(vals + 1) / np.log10(max(ladder["count"].max(), 1) + 1))
        ax.scatter(xs, ys, s=sizes, color=MODALITY_COLORS[modality], edgecolor="white", linewidth=0.8, alpha=0.88, zorder=3)
        for x, y, v in zip(xs, ys, vals):
            if v > 0:
                ax.text(x, y + 0.20, abbrev(v), ha="center", va="bottom", fontsize=6.4, color=PALETTE["dark"])
    ax.set_xticks(range(len(stage_order)))
    ax.set_xticklabels(stage_order, rotation=25, ha="right")
    ax.set_yticks([y_map[m] for m in modality_order])
    ax.set_yticklabels(modality_order)
    ax.set_xlim(-0.45, len(stage_order) - 0.55)
    ax.set_ylim(-0.55, len(modality_order) - 0.45)
    ax.grid(axis="x", color=PALETTE["grid"], lw=0.8)
    ax.set_axisbelow(True)
    ax.set_xlabel("Evidence step")


def draw_radar_panel(ax: plt.Axes, workflow: pd.DataFrame) -> None:
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    angles = np.linspace(0, 2 * np.pi, len(WORKFLOW_KEYS), endpoint=False)
    angles_closed = np.r_[angles, angles[0]]
    for modality in ["VCF transfer", "RAG case query", "Cohort tables", "New literature"]:
        sub = workflow[workflow["modality"] == modality].set_index("workflow").loc[WORKFLOW_KEYS]
        vals = sub["readiness_score"].to_numpy(dtype=float)
        vals = np.r_[vals, vals[0]]
        ax.plot(angles_closed, vals, color=MODALITY_COLORS[modality], lw=1.6, label=modality)
        ax.fill(angles_closed, vals, color=MODALITY_COLORS[modality], alpha=0.08)
    ax.set_ylim(0, 1.05)
    ax.set_yticks([0.25, 0.50, 0.75, 1.00])
    ax.set_yticklabels(["0.25", "0.50", "0.75", "1.00"], fontsize=6.4, color=PALETTE["mid"])
    ax.set_xticks(angles)
    ax.set_xticklabels(WORKFLOWS, fontsize=6.6)
    ax.grid(color=PALETTE["grid"], lw=0.75)
    ax.spines["polar"].set_color("#C8CDD6")
    ax.legend(loc="lower left", bbox_to_anchor=(-0.18, -0.12), ncol=1, frameon=False, handlelength=1.6)


def draw_missing_context_panel(ax: plt.Axes, missing: pd.DataFrame) -> None:
    contexts = list(missing["context_label"].drop_duplicates())
    modalities = ["VCF transfer", "RAG case query", "Cohort tables", "New literature"]
    y_map = {c: i for i, c in enumerate(contexts[::-1])}
    jitter = {"VCF transfer": -0.24, "RAG case query": -0.08, "Cohort tables": 0.08, "New literature": 0.24}
    for modality in modalities:
        sub = missing[missing["modality"] == modality]
        xs = sub["gap_count"].to_numpy(dtype=float)
        ys = [y_map[c] + jitter[modality] for c in sub["context_label"]]
        for x, y in zip(xs, ys):
            ax.plot([0, x], [y, y], color=MODALITY_COLORS[modality], alpha=0.20, lw=1.4)
        ax.scatter(xs, ys, s=34, color=MODALITY_COLORS[modality], edgecolor="white", linewidth=0.7, label=modality, zorder=3)
    ax.set_yticks([y_map[c] for c in contexts])
    ax.set_yticklabels(contexts)
    ax.set_xlabel("Distinct gap count")
    ax.set_xlim(-0.1, 2.35)
    ax.set_xticks([0, 1, 2])
    ax.grid(axis="x", color=PALETTE["grid"], lw=0.8)
    ax.legend(loc="lower right", frameon=False, ncol=1, handletextpad=0.3)


def draw_retrieval_panel(ax: plt.Axes, retrieval: pd.DataFrame) -> None:
    labels = ["VCF\ncomplete", "CRC\nTP53", "NSCLC\nDNMT3A"]
    y = np.arange(len(retrieval))[::-1]
    lit = retrieval["literature_hits_per_round"].to_numpy()
    rules = retrieval["rule_hits_per_round"].to_numpy()
    rounds = retrieval["rounds"].to_numpy()
    for yi, l, r, n in zip(y, lit, rules, rounds):
        ax.plot([0, max(l, r)], [yi, yi], color=PALETTE["grid"], lw=5, solid_capstyle="round")
        ax.scatter([l], [yi + 0.09], s=38 + n * 1.4, color=PALETTE["violet"], edgecolor="white", linewidth=0.8, label="Literature" if yi == y[0] else None)
        ax.scatter([r], [yi - 0.09], s=38 + n * 1.4, color=PALETTE["teal"], edgecolor="white", linewidth=0.8, label="Rule" if yi == y[0] else None)
        ax.text(max(l, r) + 0.10, yi, f"{int(n)}/{int(n)}", ha="left", va="center", fontsize=6.8, color=PALETTE["mid"])
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Hits per retrieval round")
    ax.set_xlim(0, max(max(lit), max(rules)) + 0.65)
    ax.grid(axis="x", color=PALETTE["grid"], lw=0.8)
    ax.legend(loc="lower right", frameon=False)


def draw_product_fan_panel(ax: plt.Axes, modality: pd.DataFrame) -> None:
    ax.set_axis_off()
    left_x = 0.12
    right_x = 0.82
    inputs = list(modality["modality"])
    input_y = np.linspace(0.82, 0.18, len(inputs))
    products = [
        ("Entity\nnormalization", "normalized_units", PALETTE["blue2"]),
        ("Provenance\nunits", "provenance_units", PALETTE["teal"]),
        ("Literature\nhits/chunks", "literature_hits", PALETTE["violet"]),
        ("Rule\nhits/chunks", "rule_hits", PALETTE["green"]),
        ("Review/gap\nactions", "review_or_gap_units", PALETTE["red"]),
    ]
    product_y = np.linspace(0.86, 0.14, len(products))
    x_lookup = {m: (left_x, y) for m, y in zip(inputs, input_y)}
    product_lookup = {p[0]: (right_x, y, p[2]) for p, y in zip(products, product_y)}
    max_val = max(float(modality[c].max()) for _, c, _ in products)
    for _, row in modality.iterrows():
        x0, y0 = x_lookup[row["modality"]]
        for label, col, color in products:
            val = float(row[col])
            if val <= 0:
                continue
            x1, y1, _ = product_lookup[label]
            path = MplPath(
                [(x0 + 0.035, y0), (0.38, y0), (0.56, y1), (x1 - 0.035, y1)],
                [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4],
            )
            lw = 0.7 + 4.7 * math.log10(val + 1) / math.log10(max_val + 1)
            ax.add_patch(PathPatch(path, facecolor="none", edgecolor=color, lw=lw, alpha=0.18, capstyle="round"))
    for modality_name, y0 in zip(inputs, input_y):
        ax.add_patch(Circle((left_x, y0), 0.033, facecolor=MODALITY_COLORS[modality_name], edgecolor="white", lw=0.8))
        ax.text(left_x + 0.05, y0, modality_name, ha="left", va="center", fontsize=7.7)
    for label, _, color in products:
        x, y, _ = product_lookup[label]
        ax.add_patch(Circle((x, y), 0.033, facecolor=color, edgecolor="white", lw=0.8, alpha=0.90))
        ax.text(x - 0.05, y, label, ha="right", va="center", fontsize=7.5)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


def draw_ka7p_panel(ax_left: plt.Axes, ax_right: plt.Axes, caller_records: pd.DataFrame, reductions: pd.DataFrame) -> None:
    vals = reductions["reduction_fraction"].dropna().to_numpy(dtype=float)
    parts = ax_left.violinplot(vals, positions=[0], widths=0.72, showextrema=False, showmeans=False, showmedians=False)
    for pc in parts["bodies"]:
        pc.set_facecolor(PALETTE["blue2"])
        pc.set_edgecolor("none")
        pc.set_alpha(0.25)
    rng = np.random.default_rng(17)
    jitter = rng.normal(0, 0.055, size=len(vals))
    ax_left.scatter(jitter, vals, s=22, color=PALETTE["blue"], edgecolor="white", linewidth=0.6, alpha=0.86)
    ax_left.plot([-0.28, 0.28], [np.median(vals), np.median(vals)], color=PALETTE["dark"], lw=1.2)
    ax_left.set_xlim(-0.45, 0.45)
    ax_left.set_xticks([0])
    ax_left.set_xticklabels(["Mutect2\nn=21"])
    ax_left.set_ylabel("Raw-to-filtered reduction")
    ax_left.set_ylim(0.30, 1.02)
    ax_left.grid(axis="y", color=PALETTE["grid"], lw=0.8)

    order = caller_records.sort_values("records")["caller"].tolist()
    y = np.arange(len(order))
    sub = caller_records.set_index("caller").loc[order].reset_index()
    max_records = max(float(sub["records"].max()), 1.0)
    for yi, (_, row) in zip(y, sub.iterrows()):
        color = PALETTE["light"] if row["records"] == 0 else PALETTE["orange"]
        ax_right.plot([0, row["records"]], [yi, yi], color=color, lw=3.0, alpha=0.45, solid_capstyle="round")
        ax_right.scatter([row["records"]], [yi], s=35 + 115 * row["records"] / max_records, color=color, edgecolor="white", linewidth=0.7)
        ax_right.text(row["records"] + max_records * 0.025, yi, abbrev(row["records"]), ha="left", va="center", fontsize=6.7)
    ax_right.set_yticks(y)
    ax_right.set_yticklabels(order)
    ax_right.set_xlabel("Filtered records")
    ax_right.set_xlim(0, max_records * 1.18)
    ax_right.grid(axis="x", color=PALETTE["grid"], lw=0.8)


def draw_examples_panel(ax: plt.Axes, examples: pd.DataFrame) -> None:
    ax.set_axis_off()
    cards = [
        {
            "title": "FA22310208 VCF",
            "color": MODALITY_COLORS["VCF transfer"],
            "fields": [
                "Entity: 7,652 annotated variants",
                "Evidence: 126 lit hits; 63 rule hits",
                "Gap: WBC/tumor/timepoint/stage",
                "Action: technical rescue; context review",
            ],
        },
        {
            "title": "CRC TP53 query",
            "color": MODALITY_COLORS["RAG case query"],
            "fields": [
                "Entity: CRC, TP53 p.R175H, 0.05% VAF",
                "Evidence: 18 lit hits; 12 rule hits",
                "Gap: WBC/tumor unavailable",
                "Action: repeat/serial confirmation",
            ],
        },
        {
            "title": "NSCLC DNMT3A query",
            "color": MODALITY_COLORS["RAG case query"],
            "fields": [
                "Entity: NSCLC, DNMT3A p.R882H, 0.12% VAF",
                "Evidence: 18 lit hits; 12 rule hits",
                "Gap: matched WBC unavailable",
                "Action: CHIP/hematopoietic audit",
            ],
        },
        {
            "title": "Supplementary cohort table",
            "color": MODALITY_COLORS["Cohort tables"],
            "fields": [
                "Entity: 4.34M variant events",
                "Evidence: 303,963 rules",
                "Gap: case WBC/tumor not guaranteed",
                "Action: rule-backed provenance",
            ],
        },
    ]
    card_w = 0.235
    card_h = 0.80
    xs = [0.025, 0.267, 0.509, 0.751]
    for card, x in zip(cards, xs):
        y = 0.07
        color = card["color"]
        ax.add_patch(Rectangle((x, y), card_w, card_h, facecolor="#FAFBFC", edgecolor="#E2E6ED", lw=0.8))
        ax.add_patch(Rectangle((x, y + card_h - 0.030), card_w, 0.030, facecolor=color, edgecolor="none"))
        ax.text(
            x + 0.015,
            y + card_h - 0.070,
            card["title"],
            ha="left",
            va="top",
            fontsize=7.4,
            fontweight="bold",
            color=PALETTE["dark"],
        )
        yy = y + card_h - 0.190
        for field in card["fields"]:
            ax.text(
                x + 0.015,
                yy,
                field,
                ha="left",
                va="top",
                fontsize=6.25,
                color=PALETTE["dark"],
                linespacing=1.0,
            )
            yy -= 0.145
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


def shorten_example_text(text: str) -> str:
    replacements = {
        "gene-annotated variants": "annotated variants",
        "literature hits": "lit hits",
        "materialized rules": "rules",
        "matched WBC and tumor tissue unavailable": "WBC/tumor unavailable",
        "WBC, tumor tissue, timepoint, stage and treatment context": "WBC/tumor/timepoint/stage",
        "case-specific WBC/tumor metadata not guaranteed": "case WBC/tumor not guaranteed",
        "technical rescue and contextual source tracing": "technical rescue; context review",
        "repeat/serial confirmation before actionability": "repeat/serial confirmation",
        "rule-backed cohort provenance": "rule-backed provenance",
    }
    out = text
    for old, new in replacements.items():
        out = out.replace(old, new)
    return out


def wrap_text(text: str, max_chars: int = 18, max_lines: int = 4) -> str:
    words = text.split()
    lines: list[str] = []
    current: list[str] = []
    for word in words:
        if len(" ".join(current + [word])) <= max_chars:
            current.append(word)
        else:
            if current:
                lines.append(" ".join(current))
            current = [word]
    if current:
        lines.append(" ".join(current))
    return "\n".join(lines[:max_lines])


def draw_ka7p_sampling_matrix(ax: plt.Axes, sampling: dict[str, Any]) -> None:
    presence = sampling["presence"].copy()
    cases = sampling["cases"].copy()
    if presence.empty or cases.empty:
        ax.text(0.5, 0.5, "KA7P sampling table unavailable", ha="center", va="center")
        ax.set_axis_off()
        return

    timepoints = [
        ("tumor_tissue", "Tumor\ntissue"),
        ("baseline_pre_rt", "Baseline\nblood"),
        ("rt_30Gy", "30 Gy\nblood"),
        ("rt_60Gy", "60 Gy\nblood"),
        ("one_month", "1 month\nblood"),
        ("pre_immunotherapy", "Pre-IO\nblood"),
        ("post_immunotherapy", "Post-IO\nblood"),
        ("follow_up", "Follow-up\nblood"),
    ]
    timepoint_keys = [key for key, _ in timepoints]
    timepoint_labels = [label for _, label in timepoints]
    cases = cases[cases["sample_points"] > 0].sort_values(
        ["vcf_backed_cells", "sample_points", "case_id"], ascending=[False, False, True]
    )
    case_order = list(cases["case_id"])
    x_map = {key: i for i, key in enumerate(timepoint_keys)}
    y_map = {case: len(case_order) - 1 - i for i, case in enumerate(case_order)}
    color_map = {
        "tumor_tissue": PALETTE["orange"],
        "baseline_pre_rt": PALETTE["blue"],
        "rt_30Gy": PALETTE["blue2"],
        "rt_60Gy": PALETTE["teal"],
        "one_month": PALETTE["green"],
        "pre_immunotherapy": PALETTE["violet"],
        "post_immunotherapy": PALETTE["red"],
        "follow_up": PALETTE["gold"],
    }

    for case_id in case_order:
        y = y_map[case_id]
        for key in timepoint_keys:
            ax.scatter(
                x_map[key],
                y,
                marker="s",
                s=36,
                facecolor="#F2F4F7",
                edgecolor="#E1E5EC",
                linewidth=0.5,
                zorder=1,
            )

    sub = presence[(presence["case_id"].isin(case_order)) & (presence["timepoint_key"].isin(timepoint_keys))]
    for _, row in sub.iterrows():
        if int(row["has_sample"]) == 0:
            continue
        x = x_map[row["timepoint_key"]]
        y = y_map[row["case_id"]]
        edgecolor = PALETTE["dark"] if int(row["vcf_backed"]) else "white"
        linewidth = 1.2 if int(row["vcf_backed"]) else 0.55
        ax.scatter(
            x,
            y,
            marker="s",
            s=54,
            facecolor=color_map[row["timepoint_key"]],
            edgecolor=edgecolor,
            linewidth=linewidth,
            zorder=3,
        )

    counts = presence[presence["timepoint_key"].isin(timepoint_keys)].groupby("timepoint_key")["has_sample"].sum()
    for key in timepoint_keys:
        ax.text(x_map[key], len(case_order) + 0.18, f"n={int(counts.get(key, 0))}", ha="center", va="bottom", fontsize=6.7)

    summary = sampling["summary"]
    summary_text = (
        f"{summary.get('cases_with_any_sample', 0)} cases with samples\n"
        f"{summary.get('unique_sample_ids', 0)} unique sample IDs\n"
        f"{summary.get('cases_with_4plus_samples', 0)} cases with >=4 points\n"
        f"{summary.get('vcf_backed_sample_cells', 0)} VCF-backed cells"
    )
    ax.text(
        len(timepoints) + 0.35,
        len(case_order) - 1,
        summary_text,
        ha="left",
        va="top",
        fontsize=7.2,
        color=PALETTE["dark"],
        bbox=dict(boxstyle="round,pad=0.32", facecolor="#FAFBFC", edgecolor="#E1E5EC", linewidth=0.8),
    )
    handles = [
        Line2D([0], [0], marker="s", markersize=5.8, markerfacecolor=PALETTE["orange"], markeredgecolor="white", lw=0, label="Tissue"),
        Line2D([0], [0], marker="s", markersize=5.8, markerfacecolor=PALETTE["blue"], markeredgecolor="white", lw=0, label="Radiotherapy blood"),
        Line2D([0], [0], marker="s", markersize=5.8, markerfacecolor=PALETTE["violet"], markeredgecolor="white", lw=0, label="Immunotherapy blood"),
        Line2D([0], [0], marker="s", markersize=5.8, markerfacecolor=PALETTE["green"], markeredgecolor="white", lw=0, label="Surveillance"),
        Line2D([0], [0], marker="s", markersize=5.8, markerfacecolor="#F2F4F7", markeredgecolor="#E1E5EC", lw=0, label="Not recorded"),
        Line2D([0], [0], marker="s", markersize=5.8, markerfacecolor="white", markeredgecolor=PALETTE["dark"], lw=0, label="VCF-backed sample"),
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0, -0.20), ncol=3, frameon=False, columnspacing=1.2, handletextpad=0.35)
    ax.set_xticks(range(len(timepoints)))
    ax.set_xticklabels(timepoint_labels)
    ax.set_yticks([y_map[case] for case in case_order])
    ax.set_yticklabels(case_order, fontsize=5.9)
    ax.tick_params(axis="both", length=0)
    ax.set_xlim(-0.55, len(timepoints) + 2.42)
    ax.set_ylim(-0.65, len(case_order) + 0.88)
    ax.set_title("KA7P longitudinal sampling structure", loc="left", pad=14, fontsize=9.4, fontweight="bold")
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_ka7p_caller_burden(ax: plt.Axes, caller_records: pd.DataFrame) -> None:
    sub = caller_records.copy().sort_values("records")
    y = np.arange(len(sub))
    max_records = max(float(sub["records"].max()), 1.0)
    for yi, (_, row) in zip(y, sub.iterrows()):
        color = "#CED5DF" if row["records"] == 0 else PALETTE["teal"]
        ax.plot([0, row["records"]], [yi, yi], color=color, lw=4.3, alpha=0.38, solid_capstyle="round")
        ax.scatter([row["records"]], [yi], s=42 + 105 * row["records"] / max_records, color=color, edgecolor="white", linewidth=0.8, zorder=3)
        label = f"{abbrev(row['records'])}; median {abbrev(row['median_records_per_file'])}/file"
        ax.text(row["records"] + max_records * 0.026, yi, label, ha="left", va="center", fontsize=6.6, color=PALETTE["dark"])
    ax.set_yticks(y)
    ax.set_yticklabels(sub["caller"])
    ax.set_xlabel("Filtered records across 21 KA7P cases")
    ax.set_xlim(0, max_records * 1.35)
    ax.grid(axis="x", color=PALETTE["grid"], lw=0.8)
    ax.set_title("Multi-caller VCF burden", loc="left", fontsize=9.0, fontweight="bold")


def draw_mutect2_filtering(ax: plt.Axes, reductions: pd.DataFrame) -> None:
    sub = reductions.dropna(subset=["raw_records", "filtered_records", "reduction_fraction"]).copy()
    sub = sub.sort_values("raw_records", ascending=False)
    for _, row in sub.iterrows():
        ax.plot([0, 1], [row["raw_records"], row["filtered_records"]], color="#B8C1CE", lw=0.9, alpha=0.64, zorder=1)
    ax.scatter(np.zeros(len(sub)), sub["raw_records"], s=28, color=PALETTE["orange"], edgecolor="white", linewidth=0.6, zorder=3, label="Raw")
    ax.scatter(np.ones(len(sub)), sub["filtered_records"], s=28, color=PALETTE["blue"], edgecolor="white", linewidth=0.6, zorder=3, label="Filtered")
    median_reduction = float(sub["reduction_fraction"].median())
    ax.text(
        0.50,
        0.08,
        f"Median reduction\n{median_reduction:.1%}",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=7.2,
        bbox=dict(boxstyle="round,pad=0.30", facecolor="white", edgecolor="#E1E5EC", linewidth=0.8),
    )
    ax.set_yscale("log")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Raw", "Filtered"])
    ax.set_ylabel("Mutect2 records per case")
    ax.grid(axis="y", color=PALETTE["grid"], lw=0.8, which="both")
    ax.legend(loc="upper right", frameon=False, handletextpad=0.35)
    ax.set_title("Case-level raw-to-filtered reduction", loc="left", fontsize=9.0, fontweight="bold")


def draw_qc_retention(ax: plt.Axes, tables: dict[str, Any]) -> None:
    qc = tables["summary"]["qc_totals_selected_cases"]
    smoke = tables["longitudinal_smoke"]
    rows = [
        ("Processed variants", qc.get("records", 0), PALETTE["blue"]),
        ("Gene annotated", qc.get("annotated", 0), PALETTE["blue2"]),
        ("Matched-control signal", qc.get("matched_control_signal", 0), PALETTE["teal"]),
        ("Review candidates", qc.get("review_candidate", 0), PALETTE["orange"]),
        ("KG gene evidence", qc.get("kg_gene_match", 0), PALETTE["green"]),
        ("Smoke run unresolved", smoke.get("longitudinal_variant_count", 0), PALETTE["red"]),
    ]
    y = np.arange(len(rows))[::-1]
    max_value = max(value for _, value, _ in rows)
    for yi, (label, value, color) in zip(y, rows):
        ax.barh(yi, max(value, 1), color=color, alpha=0.82, height=0.52)
        ax.text(max(value, 1) * 1.10, yi, f"{abbrev(value)}", ha="left", va="center", fontsize=7.0, color=PALETTE["dark"])
    confidence = smoke.get("source_assignment_confidence_counts", {})
    low = int(confidence.get("low_until_control_and_technical_review", 0))
    context = int(confidence.get("context_only_not_assignable", 0))
    total = max(low + context, 1)
    inset_y = -1.05
    ax.barh(inset_y, low, color="#F0A24A", alpha=0.78, height=0.30)
    ax.barh(inset_y, context, left=low, color="#B85A5A", alpha=0.78, height=0.30)
    ax.text(total * 1.10, inset_y, "low-readiness split", ha="left", va="center", fontsize=6.7, color=PALETTE["mid"])
    ax.set_yticks(list(y) + [inset_y])
    ax.set_yticklabels([label for label, _, _ in rows] + [""])
    ax.set_xscale("log")
    ax.set_xlabel("Variant or evidence units (log scale)")
    ax.set_ylim(-1.55, len(rows) - 0.35)
    ax.set_xlim(50, max_value * 1.95)
    ax.grid(axis="x", color=PALETTE["grid"], lw=0.8, which="both")
    ax.set_title("Evidence retention and readiness control", loc="left", fontsize=9.0, fontweight="bold")


def draw_graphrag_yield(ax: plt.Axes, retrieval: pd.DataFrame) -> None:
    sub = retrieval.copy()
    sub["label"] = ["KA7P VCF\ncomplete", "CRC TP53\nquery", "NSCLC DNMT3A\nquery"]
    y = np.arange(len(sub))[::-1]
    left = np.zeros(len(sub))
    ax.barh(y, sub["literature_hits"], left=left, color=PALETTE["violet"], alpha=0.82, height=0.46, label="Literature hits")
    left += sub["literature_hits"].to_numpy()
    ax.barh(y, sub["rule_hits"], left=left, color=PALETTE["green"], alpha=0.82, height=0.46, label="Rule hits")
    total_hits = sub["literature_hits"] + sub["rule_hits"]
    for yi, (_, row), total in zip(y, sub.iterrows(), total_hits):
        ax.text(total + total_hits.max() * 0.025, yi, f"{int(row['nonempty_rounds'])}/{int(row['rounds'])} rounds", ha="left", va="center", fontsize=7.0)
    ax.set_yticks(y)
    ax.set_yticklabels(sub["label"])
    ax.set_xlabel("Retrieved evidence units")
    ax.set_xlim(0, total_hits.max() * 1.25)
    ax.grid(axis="x", color=PALETTE["grid"], lw=0.8)
    ax.legend(loc="lower right", frameon=False, handlelength=1.2)
    ax.set_title("GraphRAG retrieval support", loc="left", fontsize=9.0, fontweight="bold")


def make_figure(tables: dict[str, Any]) -> None:
    setup_style()
    fig = plt.figure(figsize=(11.6, 11.35), dpi=300)
    gs = fig.add_gridspec(4, 2, height_ratios=[1.72, 1.03, 1.03, 0.92], hspace=0.84, wspace=0.42)

    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])
    ax_d = fig.add_subplot(gs[2, 0])
    ax_e = fig.add_subplot(gs[2, 1])
    ax_f = fig.add_subplot(gs[3, :])

    draw_ka7p_sampling_matrix(ax_a, tables["ka7p_sampling"])
    panel_label(ax_a, "a", x=-0.025, y=1.06)

    draw_ka7p_caller_burden(ax_b, tables["caller_records"])
    panel_label(ax_b, "b", x=-0.10, y=1.05)

    draw_mutect2_filtering(ax_c, tables["reductions"])
    panel_label(ax_c, "c", x=-0.10, y=1.05)

    draw_qc_retention(ax_d, tables)
    panel_label(ax_d, "d", x=-0.10, y=1.05)

    draw_graphrag_yield(ax_e, tables["retrieval"])
    panel_label(ax_e, "e", x=-0.10, y=1.05)

    draw_examples_panel(ax_f, tables["examples"])
    panel_label(ax_f, "f", x=-0.025, y=1.02)

    for ax in [ax_b, ax_c, ax_d, ax_e]:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_BASE.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(FIG_BASE.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(FIG_BASE.with_suffix(".svg"), bbox_inches="tight")
    try:
        from PIL import Image

        Image.open(FIG_BASE.with_suffix(".png")).convert("L").save(OUT / "figure_experiment5_graphrag_generalization_grayscale.png")
    except Exception:
        pass
    plt.close(fig)


def main() -> None:
    tables = build_source_data()
    make_figure(tables)
    print(
        json.dumps(
            {
                "figure_png": str(FIG_BASE.with_suffix(".png")),
                "figure_pdf": str(FIG_BASE.with_suffix(".pdf")),
                "figure_svg": str(FIG_BASE.with_suffix(".svg")),
                "summary_json": str(DATA / "experiment5_graphrag_generalization_summary.json"),
                "source_tables": [
                    str(DATA / "experiment5_modality_summary.csv"),
                    str(DATA / "experiment5_evidence_ladder.csv"),
                    str(DATA / "experiment5_workflow_readiness.csv"),
                    str(DATA / "experiment5_missing_context.csv"),
                    str(DATA / "experiment5_retrieval_yield.csv"),
                    str(DATA / "experiment5_provenance_network_edges.csv"),
                    str(DATA / "experiment5_representative_examples.csv"),
                    str(DATA / "experiment5_ka7p_sampling_presence.csv"),
                    str(DATA / "experiment5_ka7p_sampling_samples.csv"),
                    str(DATA / "experiment5_ka7p_sampling_timepoints.csv"),
                    str(DATA / "experiment5_ka7p_sampling_case_summary.csv"),
                    str(DATA / "experiment5_ka7p_filtered_caller_records.csv"),
                    str(DATA / "experiment5_ka7p_mutect2_reductions.csv"),
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
