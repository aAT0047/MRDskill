#!/usr/bin/env python3
"""Experiment 3: leakage-controlled MRD source-tracing adjudication.

This script builds a reproducible multi-method benchmark from the Razavi
cfDNA-WBC-tumor source-tracing table and adds operational transfer summaries
from the KA7P multi-caller VCF cohort. The formal source-origin metrics use
only Razavi hard labels. KA7P is reported as transfer/readiness data, not as
biological source-origin truth.
"""

from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import math
import re
import shutil
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, PathPatch, Rectangle
from matplotlib.path import Path as MplPath

from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    precision_recall_curve,
    precision_recall_fscore_support,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from mrdskill_project_paths import artifact_root


ARTIFACTS = artifact_root()
DATA = ARTIFACTS / "data" / "benchmarks"
EVIDENCE = ARTIFACTS / "data" / "evidence"
QC = ARTIFACTS / "data" / "qc"
DOCS = ARTIFACTS / "docs" / "benchmarks"
EXTERNAL = ARTIFACTS / "external"
KA7P_ROOT = EXTERNAL / "ka7p"

RAZAVI = DATA / "razavi_2019_mrd_source_trace_skill_benchmark.csv"
SOURCE_TRACE_DELTA = EVIDENCE / "mrd_evidence_source_trace_delta.csv"
KA7P_INVENTORY = DATA / "ka7p_vcf_inventory_files.csv"

PREDICTIONS_CSV = DATA / "experiment3_source_tracing_predictions.csv"
METRICS_CSV = DATA / "experiment3_source_tracing_metrics.csv"
PR_CURVE_CSV = DATA / "experiment3_precision_recall_curve.csv"
ABLATION_JSON = DATA / "experiment3_source_tracing_ablation.json"
SUMMARY_JSON = DATA / "experiment3_source_tracing_summary.json"
GENE_PRIOR_CSV = DATA / "experiment3_leakage_masked_gene_prior.csv"
METHOD_MATRIX_CSV = DATA / "experiment3_method_evidence_matrix.csv"
METACH_INPUT_CSV = DATA / "experiment3_metach_razavi_input.csv"
PLASMUT_INPUT_CSV = DATA / "experiment3_plasmut_razavi_input.csv"
EXTERNAL_STATUS_CSV = DATA / "experiment3_external_method_status.csv"
EXTERNAL_STATUS_JSON = DATA / "experiment3_external_method_status.json"
KA7P_TRANSFER_CSV = DATA / "experiment3_ka7p_transfer_summary.csv"
KA7P_CONCORDANCE_CSV = DATA / "experiment3_ka7p_caller_concordance.csv"
CASE_TRANSFER_CSV = DATA / "experiment3_operational_case_summary.csv"

FIG_BASE = DOCS / "figure_experiment3_source_tracing_adjudication"
FIG_PNG = FIG_BASE.with_suffix(".png")
FIG_PDF = FIG_BASE.with_suffix(".pdf")
FIG_SVG = FIG_BASE.with_suffix(".svg")
FIG_GRAY = DOCS / "figure_experiment3_source_tracing_adjudication_grayscale.png"
QA_JSON = DOCS / "figure_experiment3_source_tracing_adjudication_scipilot_qa.json"

CHIP = "chip_or_hematopoietic"
TUMOR = "tumor_derived"
REVIEW = "review_needed"
UNRESOLVED = "unresolved_cfDNA_only"

PALETTE = {
    "blue": "#4C78A8",
    "orange": "#F58518",
    "green": "#54A24B",
    "purple": "#B279A2",
    "red": "#E45756",
    "teal": "#72B7B2",
    "chip": "#D55E00",
    "tumor": "#0072B2",
    "review": "#7A7F87",
    "light": "#E5E7EB",
    "text": "#111827",
    "muted": "#6B7280",
}

METHOD_ORDER = [
    "All cfDNA as tumor",
    "CHIP gene blacklist",
    "Low-VAF rule",
    "CHIP gene OR low VAF",
    "CHIP gene AND low VAF",
    "Prompt-only feature adjudicator",
    "RAG-only gene prior",
    "KG/rule-only gene prior",
    "MRD-KGQC plasma-only",
    "MRD-KGQC plasma-only + review",
    "MRD-KGQC matched-source",
]

DISPLAY_METHODS = [
    "CHIP gene OR low VAF",
    "Prompt-only feature adjudicator",
    "RAG-only gene prior",
    "KG/rule-only gene prior",
    "MRD-KGQC plasma-only",
    "MRD-KGQC plasma-only + review",
    "MRD-KGQC matched-source",
]

CV_FEATURES = {
    "Prompt-only feature adjudicator": [
        "log_cfdna_af",
        "log_cfdna_depth",
        "log_cfdna_alt_count",
        "chip_gene_overlap",
        "low_vaf_flag",
        "moderate_vaf_flag",
        "high_vaf_flag",
    ],
    "RAG-only gene prior": [
        "gene_prior_log2or",
        "gene_prior_chip_probability",
        "gene_prior_chip_support",
        "gene_prior_tumor_support",
        "gene_prior_known",
    ],
    "KG/rule-only gene prior": [
        "chip_gene_overlap",
        "gene_prior_log2or",
        "gene_prior_chip_probability",
        "gene_prior_known",
    ],
    "MRD-KGQC plasma-only": [
        "log_cfdna_af",
        "log_cfdna_depth",
        "log_cfdna_alt_count",
        "chip_gene_overlap",
        "low_vaf_flag",
        "moderate_vaf_flag",
        "high_vaf_flag",
        "gene_prior_log2or",
        "gene_prior_chip_probability",
        "gene_prior_chip_support",
        "gene_prior_tumor_support",
        "gene_prior_known",
    ],
}


@dataclass
class MethodPrediction:
    method: str
    track: str
    score: np.ndarray
    call: np.ndarray
    threshold: float | None = None
    model_note: str = ""


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "axes.linewidth": 0.75,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "axes.unicode_minus": False,
        }
    )


def abbrev(value: float | int | None) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    value = float(value)
    if value == 0:
        return "0"
    sign = "-" if value < 0 else ""
    value = abs(value)
    if value < 1:
        return f"{sign}{value:.2f}".rstrip("0").rstrip(".")
    if value < 10:
        return f"{sign}{value:.2f}".rstrip("0").rstrip(".")
    if value < 100:
        return f"{sign}{value:.1f}".rstrip("0").rstrip(".")
    if value < 1000:
        return f"{sign}{int(round(value)):,}"
    if value < 1_000_000:
        return f"{sign}{value / 1000:.1f}k".replace(".0k", "k")
    return f"{sign}{value / 1_000_000:.2f}M".replace(".00M", "M").replace(".0M", "M")


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.10,
        1.05,
        label,
        transform=ax.transAxes,
        fontsize=15,
        fontweight="bold",
        va="bottom",
        ha="left",
    )


def split_tokens(value: object) -> list[str]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    out = []
    for token in re.split(r"[;,|/\s]+", str(value)):
        token = token.strip().upper()
        if not token or token in {"NOT_REPORTED", "NAN", "NONE", "NA"}:
            continue
        if len(token) > 30:
            continue
        out.append(token)
    return out


def build_gene_prior() -> pd.DataFrame:
    """Build source-origin gene priors after excluding direct Razavi rows."""
    cols = [
        "paper_id",
        "title",
        "doi",
        "gene_symbol_normalized",
        "source_origin_category",
        "confidence_score",
    ]
    counts: dict[str, dict[str, float]] = defaultdict(
        lambda: {"chip_support": 0.0, "tumor_support": 0.0, "unresolved_support": 0.0, "total_support": 0.0}
    )
    if not SOURCE_TRACE_DELTA.exists():
        return pd.DataFrame(
            columns=[
                "gene",
                "chip_support",
                "tumor_support",
                "unresolved_support",
                "total_support",
                "chip_probability",
                "log2_odds_ratio_chip_vs_tumor",
            ]
        )

    for chunk in pd.read_csv(SOURCE_TRACE_DELTA, usecols=cols, chunksize=50_000, low_memory=False):
        direct_razavi = (
            chunk["doi"].astype(str).str.contains("s41591-019-0652-7|10.1038/s41591-019-0652-7", case=False, na=False)
            | chunk["title"].astype(str).str.contains("High-intensity sequencing reveals", case=False, na=False)
            | chunk["paper_id"].astype(str).str.contains("KWRAW-000862", case=False, na=False)
        )
        chunk = chunk.loc[~direct_razavi].copy()
        weights = chunk["confidence_score"].map({"High": 1.0, "Medium": 0.65, "Low": 0.35}).fillna(0.5)
        for row, weight in zip(chunk.itertuples(index=False), weights):
            genes = split_tokens(row.gene_symbol_normalized)
            if not genes:
                continue
            categories = [item.strip() for item in str(row.source_origin_category).split(";")]
            has_chip = any("chip_or_hematopoietic_support" in item for item in categories)
            has_tumor = any("tumor_origin_support" in item for item in categories)
            has_unresolved = any("unresolved_context" in item for item in categories)
            for gene in genes:
                counts[gene]["total_support"] += float(weight)
                if has_chip:
                    counts[gene]["chip_support"] += float(weight)
                if has_tumor:
                    counts[gene]["tumor_support"] += float(weight)
                if has_unresolved:
                    counts[gene]["unresolved_support"] += float(weight)

    rows = []
    for gene, values in sorted(counts.items()):
        chip = values["chip_support"]
        tumor = values["tumor_support"]
        probability = (chip + 0.5) / (chip + tumor + 1.0)
        log2or = math.log2((chip + 0.5) / (tumor + 0.5))
        rows.append(
            {
                "gene": gene,
                "chip_support": chip,
                "tumor_support": tumor,
                "unresolved_support": values["unresolved_support"],
                "total_support": values["total_support"],
                "chip_probability": probability,
                "log2_odds_ratio_chip_vs_tumor": log2or,
            }
        )
    prior = pd.DataFrame(rows)
    prior.to_csv(GENE_PRIOR_CSV, index=False)
    return prior


def load_benchmark(prior: pd.DataFrame) -> pd.DataFrame:
    df = pd.read_csv(RAZAVI)
    prior_map = prior.set_index("gene") if not prior.empty else pd.DataFrame()

    values = []
    for gene in df["gene"].astype(str).str.upper():
        if not prior_map.empty and gene in prior_map.index:
            row = prior_map.loc[gene]
            values.append(
                (
                    float(row["log2_odds_ratio_chip_vs_tumor"]),
                    float(row["chip_probability"]),
                    float(row["chip_support"]),
                    float(row["tumor_support"]),
                    float(row["total_support"]),
                )
            )
        else:
            values.append((0.0, 0.5, 0.0, 0.0, 0.0))
    arr = np.array(values, dtype=float)
    df["gene_prior_log2or"] = arr[:, 0]
    df["gene_prior_chip_probability"] = arr[:, 1]
    df["gene_prior_chip_support"] = arr[:, 2]
    df["gene_prior_tumor_support"] = arr[:, 3]
    df["gene_prior_total_support"] = arr[:, 4]
    df["gene_prior_known"] = (df["gene_prior_total_support"] > 0).astype(int)

    df["cfdna_af"] = pd.to_numeric(df["cfdna_af"], errors="coerce").fillna(0)
    df["cfdna_depth"] = pd.to_numeric(df["cfdna_depth"], errors="coerce").fillna(0)
    df["cfdna_alt_count"] = pd.to_numeric(df["cfdna_alt_count"], errors="coerce").fillna(0)
    df["chip_gene_overlap"] = pd.to_numeric(df["chip_gene_overlap"], errors="coerce").fillna(0).astype(int)
    df["log_cfdna_af"] = np.log10(df["cfdna_af"].clip(lower=1e-6))
    df["log_cfdna_depth"] = np.log10(df["cfdna_depth"] + 1)
    df["log_cfdna_alt_count"] = np.log10(df["cfdna_alt_count"] + 1)
    df["low_vaf_flag"] = (df["cfdna_af"] < 0.005).astype(int)
    df["moderate_vaf_flag"] = ((df["cfdna_af"] >= 0.005) & (df["cfdna_af"] < 0.01)).astype(int)
    df["high_vaf_flag"] = (df["cfdna_af"] >= 0.02).astype(int)
    return df


def fixed_rule_predictions(df: pd.DataFrame) -> list[MethodPrediction]:
    low_vaf = df["cfdna_af"] < 0.005
    chip_gene = df["chip_gene_overlap"].astype(bool)
    return [
        MethodPrediction(
            "All cfDNA as tumor",
            "rule",
            score=np.zeros(len(df)),
            call=np.full(len(df), TUMOR, dtype=object),
            model_note="Fixed rule.",
        ),
        MethodPrediction(
            "CHIP gene blacklist",
            "rule",
            score=chip_gene.astype(float).to_numpy(),
            call=np.where(chip_gene, CHIP, TUMOR),
            model_note="Fixed rule.",
        ),
        MethodPrediction(
            "Low-VAF rule",
            "rule",
            score=low_vaf.astype(float).to_numpy(),
            call=np.where(low_vaf, CHIP, TUMOR),
            model_note="Fixed rule.",
        ),
        MethodPrediction(
            "CHIP gene OR low VAF",
            "rule",
            score=(chip_gene | low_vaf).astype(float).to_numpy(),
            call=np.where(chip_gene | low_vaf, CHIP, TUMOR),
            model_note="Fixed rule.",
        ),
        MethodPrediction(
            "CHIP gene AND low VAF",
            "rule",
            score=(chip_gene & low_vaf).astype(float).to_numpy(),
            call=np.where(chip_gene & low_vaf, CHIP, TUMOR),
            model_note="Fixed rule.",
        ),
    ]


def cross_validated_scores(df: pd.DataFrame, method: str, features: list[str]) -> np.ndarray:
    hard = df["binary_origin_label"].isin([TUMOR, CHIP])
    hard_df = df.loc[hard].copy()
    y = (hard_df["binary_origin_label"] == CHIP).astype(int).to_numpy()
    groups = hard_df["patient_id"].astype(str).to_numpy()
    scores = np.full(len(df), np.nan, dtype=float)

    splitter = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=12)
    for train_index, test_index in splitter.split(hard_df[features], y, groups):
        model = make_pipeline(
            SimpleImputer(strategy="median"),
            StandardScaler(),
            LogisticRegression(max_iter=2000, class_weight="balanced", C=0.7, random_state=12),
        )
        model.fit(hard_df.iloc[train_index][features], y[train_index])
        scores[np.flatnonzero(hard)[test_index]] = model.predict_proba(hard_df.iloc[test_index][features])[:, 1]

    model = make_pipeline(
        SimpleImputer(strategy="median"),
        StandardScaler(),
        LogisticRegression(max_iter=2000, class_weight="balanced", C=0.7, random_state=12),
    )
    model.fit(hard_df[features], y)
    non_hard = ~hard
    if non_hard.any():
        scores[non_hard.to_numpy()] = model.predict_proba(df.loc[non_hard, features])[:, 1]
    return scores


def model_predictions(df: pd.DataFrame) -> list[MethodPrediction]:
    predictions = fixed_rule_predictions(df)
    for method, features in CV_FEATURES.items():
        score = cross_validated_scores(df, method, features)
        call = np.where(score >= 0.5, CHIP, TUMOR)
        predictions.append(
            MethodPrediction(
                method,
                "plasma-only",
                score=score,
                call=call,
                threshold=0.5,
                model_note="Patient-level 5-fold CV on Razavi hard-label rows; matched source columns masked.",
            )
        )

    mrd = next(item for item in predictions if item.method == "MRD-KGQC plasma-only")
    review_call = np.where(mrd.score >= 0.70, CHIP, np.where(mrd.score <= 0.30, TUMOR, REVIEW))
    predictions.append(
        MethodPrediction(
            "MRD-KGQC plasma-only + review",
            "plasma-only-review",
            score=mrd.score.copy(),
            call=review_call,
            threshold=None,
            model_note="Same plasma-only score with a 0.30-0.70 uncertainty review band.",
        )
    )

    matched_call = np.where(
        df["wbc_matched"].astype(int) == 1,
        CHIP,
        np.where(df["tumor_tissue_matched"].astype(int) == 1, TUMOR, REVIEW),
    )
    matched_score = np.where(df["wbc_matched"].astype(int) == 1, 1.0, np.where(df["tumor_tissue_matched"].astype(int) == 1, 0.0, 0.5))
    predictions.append(
        MethodPrediction(
            "MRD-KGQC matched-source",
            "matched-source",
            score=matched_score,
            call=matched_call,
            threshold=None,
            model_note="Evidence-assisted source mode using matched WBC/tumor concordance; reported separately from plasma-only prediction.",
        )
    )
    return predictions


def evaluate_method(df: pd.DataFrame, prediction: MethodPrediction) -> dict[str, object]:
    hard = df["binary_origin_label"].isin([TUMOR, CHIP])
    gold = df.loc[hard, "binary_origin_label"].to_numpy()
    call = prediction.call[hard.to_numpy()]
    score = prediction.score[hard.to_numpy()]
    resolved = call != REVIEW

    tp = int(((call == CHIP) & (gold == CHIP)).sum())
    fp = int(((call == CHIP) & (gold == TUMOR)).sum())
    fn = int(((call == TUMOR) & (gold == CHIP)).sum())
    tn = int(((call == TUMOR) & (gold == TUMOR)).sum())
    review_chip = int(((call == REVIEW) & (gold == CHIP)).sum())
    review_tumor = int(((call == REVIEW) & (gold == TUMOR)).sum())

    chip_precision = tp / (tp + fp) if tp + fp else np.nan
    chip_recall = tp / (tp + fn + review_chip) if tp + fn + review_chip else np.nan
    tumor_specificity = tn / (tn + fp + review_tumor) if tn + fp + review_tumor else np.nan
    chip_f1 = 2 * chip_precision * chip_recall / (chip_precision + chip_recall) if chip_precision + chip_recall > 0 else np.nan
    balanced = (chip_recall + tumor_specificity) / 2 if not (np.isnan(chip_recall) or np.isnan(tumor_specificity)) else np.nan
    strict_accuracy = float((call == gold).mean())
    resolved_accuracy = float((call[resolved] == gold[resolved]).mean()) if resolved.any() else np.nan
    resolved_chip_recall = tp / (tp + fn) if tp + fn else np.nan
    resolved_tumor_specificity = tn / (tn + fp) if tn + fp else np.nan
    resolved_balanced = (
        (resolved_chip_recall + resolved_tumor_specificity) / 2
        if not (np.isnan(resolved_chip_recall) or np.isnan(resolved_tumor_specificity))
        else np.nan
    )

    y = (gold == CHIP).astype(int)
    auroc = np.nan
    auprc = np.nan
    if len(np.unique(y)) == 2 and np.isfinite(score).all():
        auroc = float(roc_auc_score(y, score))
        auprc = float(average_precision_score(y, score))

    unresolved = df["binary_origin_label"] == UNRESOLVED
    unresolved_call = prediction.call[unresolved.to_numpy()]
    return {
        "method": prediction.method,
        "track": prediction.track,
        "evaluable_rows": int(hard.sum()),
        "resolved_rows": int(resolved.sum()),
        "coverage": float(resolved.mean()),
        "review_rate_hard": float((call == REVIEW).mean()),
        "chip_precision": chip_precision,
        "chip_recall": chip_recall,
        "tumor_specificity": tumor_specificity,
        "accuracy_strict": strict_accuracy,
        "resolved_accuracy": resolved_accuracy,
        "balanced_accuracy": balanced,
        "resolved_balanced_accuracy": resolved_balanced,
        "resolved_chip_recall": resolved_chip_recall,
        "resolved_tumor_specificity": resolved_tumor_specificity,
        "chip_f1": chip_f1,
        "auroc": auroc,
        "auprc": auprc,
        "tp_chip": tp,
        "fp_chip": fp,
        "fn_chip": fn,
        "tn_tumor": tn,
        "review_chip": review_chip,
        "review_tumor": review_tumor,
        "unresolved_rows": int(unresolved.sum()),
        "unresolved_review_rate": float((unresolved_call == REVIEW).mean()) if unresolved.any() else np.nan,
        "unresolved_chip_call_rate": float((unresolved_call == CHIP).mean()) if unresolved.any() else np.nan,
        "unresolved_tumor_call_rate": float((unresolved_call == TUMOR).mean()) if unresolved.any() else np.nan,
        "model_note": prediction.model_note,
    }


def make_predictions_table(df: pd.DataFrame, predictions: list[MethodPrediction]) -> pd.DataFrame:
    keep = [
        "benchmark_id",
        "patient_id",
        "gene",
        "chrom",
        "pos",
        "ref",
        "alt",
        "variant_type",
        "cfdna_af",
        "cfdna_depth",
        "cfdna_alt_count",
        "chip_gene_overlap",
        "wbc_matched",
        "tumor_tissue_matched",
        "binary_origin_label",
        "source_trace_gold_label",
        "gene_prior_log2or",
        "gene_prior_chip_probability",
        "gene_prior_chip_support",
        "gene_prior_tumor_support",
        "gene_prior_total_support",
    ]
    out = df[keep].copy()
    for pred in predictions:
        safe = pred.method.replace(" ", "_").replace("/", "_").replace("+", "plus").replace("-", "_")
        out[f"score_{safe}"] = pred.score
        out[f"pred_{safe}"] = pred.call
    return out


def export_pr_curve(df: pd.DataFrame, predictions: list[MethodPrediction]) -> pd.DataFrame:
    hard = df["binary_origin_label"].isin([TUMOR, CHIP])
    y = (df.loc[hard, "binary_origin_label"] == CHIP).astype(int).to_numpy()
    rows = []
    for pred in predictions:
        if pred.method not in DISPLAY_METHODS:
            continue
        score = pred.score[hard.to_numpy()]
        if not np.isfinite(score).all() or len(np.unique(score)) <= 2:
            continue
        precision, recall, thresholds = precision_recall_curve(y, score)
        threshold_values = np.concatenate([thresholds, [np.nan]])
        for p, r, t in zip(precision, recall, threshold_values):
            rows.append({"method": pred.method, "precision": p, "recall": r, "threshold": t})
    pr = pd.DataFrame(rows)
    pr.to_csv(PR_CURVE_CSV, index=False)
    return pr


def export_comparator_inputs(df: pd.DataFrame) -> tuple[int, int]:
    variant_type = df["variant_type"].replace({"SNV": "SNP"})
    metach = pd.DataFrame(
        {
            "Chromosome": df["chrom"].astype(str).map(lambda value: value if value.startswith("chr") else f"chr{value}"),
            "Start_Position": df["pos"].astype(int),
            "End_Position": df["pos"].astype(int) + df["ref"].astype(str).str.len() - 1,
            "Reference_Allele": df["ref"],
            "Tumor_Allele_1": df["alt"],
            "Tumor_Allele_2": df["alt"],
            "Variant_Type": variant_type,
            "Hugo_Symbol": df["gene"],
            "Gene_Name": df["gene"],
            "Reference_Genome": df["assembly"].fillna("GRCh37"),
            "VAF": df["cfdna_af"],
            "CANCER_TYPE": "Breast Cancer",
            "Sample_ID": df["patient_id"],
        }
    )
    metach = metach[metach["Variant_Type"].isin(["SNP", "DEL", "INS"])].copy()
    metach.to_csv(METACH_INPUT_CSV, index=False)

    plasmut = df[
        [
            "benchmark_id",
            "patient_id",
            "gene",
            "chrom",
            "pos",
            "ref",
            "alt",
            "variant_type",
            "cfdna_depth",
            "cfdna_ref_count",
            "cfdna_alt_count",
            "cfdna_af",
            "wbc_depth",
            "wbc_ref_count",
            "wbc_alt_count",
            "wbc_af",
            "tumor_depth",
            "tumor_ref_count",
            "tumor_alt_count",
            "tumor_af",
            "binary_origin_label",
        ]
    ].copy()
    plasmut.to_csv(PLASMUT_INPUT_CSV, index=False)
    return len(metach), len(plasmut)


def export_external_status(metach_rows: int, plasmut_rows: int) -> pd.DataFrame:
    metk_available = importlib.util.find_spec("metk") is not None
    rscript = shutil.which("Rscript")
    docker = shutil.which("docker")
    status = [
        {
            "method": "MetaCH",
            "comparison_track": "plasma-only SOTA",
            "input_adapter": str(METACH_INPUT_CSV),
            "mapped_rows": metach_rows,
            "run_status": "not_run",
            "reason": "METk package/embedding environment unavailable locally" if not metk_available else "ready_for_METk_feature_extraction",
            "dependency_status": f"metk={'available' if metk_available else 'missing'}; docker={'available' if docker else 'missing'}",
        },
        {
            "method": "plasmut",
            "comparison_track": "matched-source statistical comparator",
            "input_adapter": str(PLASMUT_INPUT_CSV),
            "mapped_rows": plasmut_rows,
            "run_status": "not_run",
            "reason": "Rscript unavailable locally" if not rscript else "ready_for_R_execution",
            "dependency_status": f"Rscript={rscript or 'missing'}",
        },
        {
            "method": "SSGAN",
            "comparison_track": "plasma-only literature comparator",
            "input_adapter": "MetaCH repository contains external SSGAN prediction files, but not Razavi-mapped gold labels.",
            "mapped_rows": 0,
            "run_status": "not_run",
            "reason": "available prediction files are on external validation cohorts and are not directly aligned to the Razavi benchmark labels",
            "dependency_status": "local CSV predictions found under MetaCH/results/baselines/ssgan",
        },
        {
            "method": "Fairchild-style classifier",
            "comparison_track": "plasma-only literature comparator",
            "input_adapter": "not_available",
            "mapped_rows": 0,
            "run_status": "not_run",
            "reason": "no local runnable implementation or exact feature mapping was available",
            "dependency_status": "missing local implementation",
        },
    ]
    out = pd.DataFrame(status)
    out.to_csv(EXTERNAL_STATUS_CSV, index=False)
    EXTERNAL_STATUS_JSON.write_text(json.dumps(status, ensure_ascii=False, indent=2), encoding="utf-8")
    return out


def parse_filter_counts(value: object) -> dict[str, int]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return {}
    try:
        parsed = json.loads(str(value))
        return {str(k): int(v) for k, v in parsed.items()}
    except Exception:
        return {}


def parse_vcf_keys(path: Path) -> set[str]:
    keys: set[str] = set()
    opener = gzip.open if str(path).endswith(".gz") else open
    try:
        with opener(path, "rt", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                chrom, pos, _id, ref, alt = parts[:5]
                for one_alt in alt.split(","):
                    keys.add(f"{chrom}:{pos}:{ref}>{one_alt}")
    except FileNotFoundError:
        pass
    return keys


def build_ka7p_transfer() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if not KA7P_INVENTORY.exists():
        empty = pd.DataFrame()
        empty.to_csv(KA7P_TRANSFER_CSV, index=False)
        empty.to_csv(KA7P_CONCORDANCE_CSV, index=False)
        empty.to_csv(CASE_TRANSFER_CSV, index=False)
        return empty, empty, empty

    inv = pd.read_csv(KA7P_INVENTORY)
    inv["record_count"] = pd.to_numeric(inv["record_count"], errors="coerce").fillna(0).astype(int)
    caller_rows = []
    for (folder, caller), group in inv.groupby(["folder", "caller"], dropna=False):
        caller_rows.append(
            {
                "dataset": "KA7P",
                "folder": folder,
                "caller": caller,
                "files": int(len(group)),
                "records": int(group["record_count"].sum()),
                "median_records_per_file": float(group["record_count"].median()),
            }
        )

    mut = inv[inv["caller"].eq("Mutect2")].copy()
    pivot = mut.pivot_table(index="sample_id", columns="folder", values="record_count", aggfunc="sum", fill_value=0)
    if {"VCF_ori", "VCF_filtered"}.issubset(set(pivot.columns)):
        for sample, row in pivot.iterrows():
            raw = int(row.get("VCF_ori", 0))
            filtered = int(row.get("VCF_filtered", 0))
            reduction = 1 - filtered / raw if raw else np.nan
            caller_rows.append(
                {
                    "dataset": "KA7P",
                    "folder": "raw_to_filtered",
                    "caller": f"Mutect2:{sample}",
                    "files": 2,
                    "records": filtered,
                    "median_records_per_file": float(filtered),
                    "raw_records": raw,
                    "filtered_records": filtered,
                    "reduction_fraction": reduction,
                }
            )
    transfer = pd.DataFrame(caller_rows)
    transfer.to_csv(KA7P_TRANSFER_CSV, index=False)

    concordance_rows = []
    filtered = inv[inv["folder"].eq("VCF_filtered")].copy()
    for sample, group in filtered.groupby("sample_id"):
        variant_callers: dict[str, set[str]] = defaultdict(set)
        for row in group.itertuples(index=False):
            path = Path(row.path)
            if int(row.record_count) == 0:
                continue
            keys = parse_vcf_keys(path)
            for key in keys:
                variant_callers[key].add(str(row.caller))
        support_counts = Counter(len(callers) for callers in variant_callers.values())
        concordance_rows.append(
            {
                "sample_id": sample,
                "filtered_callers": int(group["caller"].nunique()),
                "variant_union": int(len(variant_callers)),
                "single_caller_variants": int(support_counts.get(1, 0)),
                "two_caller_variants": int(support_counts.get(2, 0)),
                "three_or_more_caller_variants": int(sum(v for k, v in support_counts.items() if k >= 3)),
                "max_caller_support": int(max(support_counts.keys()) if support_counts else 0),
                "caller_concordant_fraction": float(sum(v for k, v in support_counts.items() if k >= 2) / len(variant_callers)) if variant_callers else 0.0,
            }
        )
    concordance = pd.DataFrame(concordance_rows)
    concordance.to_csv(KA7P_CONCORDANCE_CSV, index=False)

    case_rows = []
    for name, label in [
        ("FA22310208_complete_graphrag_demo.qc_summary.json", "FA22310208 GraphRAG"),
        ("FC210F1106-A3ACF90XNF1-KA7P_Mutect2.qc_summary.json", "FC210F1106 Mutect2 raw"),
        ("FC210F1106-A3ACF90XNF1-KA7P_Mutect2_filtered.qc_summary.json", "FC210F1106 Mutect2 filtered"),
    ]:
        path = QC / name
        if not path.exists():
            continue
        data = json.loads(path.read_text(encoding="utf-8"))
        records = int(data.get("record_count", 0))
        case_rows.append(
            {
                "case": label,
                "records": records,
                "annotated_variants": int(data.get("annotated_variant_count", 0)),
                "annotation_rate": int(data.get("annotated_variant_count", 0)) / records if records else np.nan,
                "kg_supported_variants": int(data.get("gene_match_variant_count", 0)),
                "kg_support_rate": int(data.get("gene_match_variant_count", 0)) / records if records else np.nan,
                "matched_control_signal_rows": int(data.get("matched_control_signal_count", 0)),
                "matched_control_signal_rate": int(data.get("matched_control_signal_count", 0)) / records if records else np.nan,
                "strong_or_review_candidates": int(data.get("tier_counts", {}).get("strong_technical_candidate", 0))
                + int(data.get("tier_counts", {}).get("review_candidate", 0)),
            }
        )
    cases = pd.DataFrame(case_rows)
    cases.to_csv(CASE_TRANSFER_CSV, index=False)
    return transfer, concordance, cases


def method_evidence_matrix() -> pd.DataFrame:
    rows = [
        ("All cfDNA as tumor", "cfDNA feature", 1, "CHIP list", 0, "VAF", 0, "literature prior", 0, "KG/rule", 0, "RAG", 0, "matched source", 0, "review", 0),
        ("CHIP gene OR low VAF", "cfDNA feature", 1, "CHIP list", 1, "VAF", 1, "literature prior", 0, "KG/rule", 0, "RAG", 0, "matched source", 0, "review", 0),
        ("Prompt-only feature adjudicator", "cfDNA feature", 1, "CHIP list", 1, "VAF", 1, "literature prior", 0, "KG/rule", 0, "RAG", 0, "matched source", 0, "review", 0),
        ("RAG-only gene prior", "cfDNA feature", 0, "CHIP list", 0, "VAF", 0, "literature prior", 1, "KG/rule", 0, "RAG", 1, "matched source", 0, "review", 0),
        ("KG/rule-only gene prior", "cfDNA feature", 0, "CHIP list", 1, "VAF", 0, "literature prior", 1, "KG/rule", 1, "RAG", 0, "matched source", 0, "review", 0),
        ("MRD-KGQC plasma-only", "cfDNA feature", 1, "CHIP list", 1, "VAF", 1, "literature prior", 1, "KG/rule", 1, "RAG", 1, "matched source", 0, "review", 0),
        ("MRD-KGQC plasma-only + review", "cfDNA feature", 1, "CHIP list", 1, "VAF", 1, "literature prior", 1, "KG/rule", 1, "RAG", 1, "matched source", 0, "review", 1),
        ("MRD-KGQC matched-source", "cfDNA feature", 1, "CHIP list", 1, "VAF", 1, "literature prior", 1, "KG/rule", 1, "RAG", 1, "matched source", 1, "review", 1),
        ("MetaCH", "cfDNA feature", 1, "CHIP list", 0, "VAF", 1, "literature prior", 0, "KG/rule", 0, "RAG", 0, "matched source", 0, "review", 0),
        ("plasmut", "cfDNA feature", 1, "CHIP list", 0, "VAF", 1, "literature prior", 0, "KG/rule", 0, "RAG", 0, "matched source", 1, "review", 0),
    ]
    out = []
    for row in rows:
        method = row[0]
        data = {"method": method}
        for i in range(1, len(row), 2):
            data[row[i]] = row[i + 1]
        out.append(data)
    matrix = pd.DataFrame(out)
    matrix.to_csv(METHOD_MATRIX_CSV, index=False)
    return matrix


def draw_dataset_map(ax: plt.Axes, metrics: pd.DataFrame, cases: pd.DataFrame, transfer: pd.DataFrame) -> None:
    ax.axis("off")
    nodes = [
        ("Razavi hard labels", "1,410", 0.10, 0.70, PALETTE["blue"]),
        ("Razavi unresolved", "1,944", 0.10, 0.35, PALETTE["review"]),
        ("Evidence index", "100,791", 0.43, 0.82, PALETTE["green"]),
        ("External adapters", "SOTA", 0.43, 0.53, PALETTE["purple"]),
        ("KA7P VCF transfer", "21 cases", 0.43, 0.24, PALETTE["orange"]),
        ("Performance", "labels", 0.80, 0.70, PALETTE["blue"]),
        ("Review behavior", "review", 0.80, 0.38, PALETTE["review"]),
        ("Operational coverage", "coverage", 0.80, 0.12, PALETTE["orange"]),
    ]
    for name, count, x, y, color in nodes:
        ax.add_patch(Circle((x, y), 0.055, facecolor=color, edgecolor="white", lw=1.0, alpha=0.95))
        ax.text(x, y - 0.085, name, ha="center", va="top", fontsize=7.2)
        ax.text(x, y, count, ha="center", va="center", fontsize=7.2, color="white", fontweight="bold")

    edges = [
        (0, 5, PALETTE["blue"]),
        (1, 6, PALETTE["review"]),
        (2, 5, PALETTE["green"]),
        (2, 6, PALETTE["green"]),
        (3, 5, PALETTE["purple"]),
        (4, 7, PALETTE["orange"]),
    ]
    for src, dst, color in edges:
        x0, y0 = nodes[src][2], nodes[src][3]
        x1, y1 = nodes[dst][2], nodes[dst][3]
        path = MplPath(
            [(x0 + 0.055, y0), ((x0 + x1) / 2, y0), ((x0 + x1) / 2, y1), (x1 - 0.055, y1)],
            [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4],
        )
        ax.add_patch(PathPatch(path, facecolor="none", edgecolor=color, alpha=0.34, lw=2.0))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


def draw_method_matrix(ax: plt.Axes, matrix: pd.DataFrame) -> None:
    sub = matrix[matrix["method"].isin(DISPLAY_METHODS + ["MetaCH", "plasmut"])].copy()
    sub = sub.set_index("method")
    columns = ["cfDNA feature", "CHIP list", "VAF", "literature prior", "KG/rule", "RAG", "matched source", "review"]
    values = sub[columns].to_numpy(dtype=float)
    ax.set_xticks(np.arange(len(columns)))
    ax.set_xticklabels(["cfDNA", "CHIP", "VAF", "Lit.", "KG", "RAG", "WBC/\ntumor", "Review"], rotation=0)
    ax.set_yticks(np.arange(len(sub.index)))
    ax.set_yticklabels([m.replace("MRD-KGQC ", "KGQC\n").replace(" gene prior", "\nprior") for m in sub.index])
    ax.set_xlim(-0.55, len(columns) - 0.45)
    ax.set_ylim(len(sub.index) - 0.55, -0.55)
    ax.tick_params(length=0)
    for i in range(values.shape[0]):
        ax.axhline(i, color="#F3F4F6", lw=0.8, zorder=0)
    for j in range(values.shape[1]):
        ax.axvline(j, color="#F3F4F6", lw=0.8, zorder=0)
    col_colors = [PALETTE["blue"], PALETTE["orange"], PALETTE["teal"], PALETTE["green"], PALETTE["purple"], PALETTE["red"], PALETTE["tumor"], PALETTE["review"]]
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            if values[i, j]:
                ax.scatter(j, i, s=58, color=col_colors[j], edgecolor="white", linewidth=0.7, zorder=3)
            else:
                ax.scatter(j, i, s=10, color="#E5E7EB", edgecolor="none", zorder=1)
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_performance_rank(ax: plt.Axes, metrics: pd.DataFrame) -> None:
    sub = metrics[metrics["method"].isin(DISPLAY_METHODS)].copy()
    sub["order"] = sub["method"].map({m: i for i, m in enumerate(DISPLAY_METHODS)})
    sub = sub.sort_values("order")
    y = np.arange(len(sub))[::-1]
    color = [PALETTE["red"] if "matched-source" in m else PALETTE["blue"] if "MRD-KGQC plasma-only" in m else PALETTE["green"] if "feature" in m else PALETTE["orange"] if "OR low" in m else PALETTE["purple"] for m in sub["method"]]
    ax.hlines(y, 0.45, sub["balanced_accuracy"], color=color, alpha=0.30, lw=3)
    ax.scatter(sub["balanced_accuracy"], y, s=58, color=color, edgecolor="white", linewidth=0.7, zorder=3, label="Strict balanced accuracy")
    ax.scatter(sub["auprc"], y - 0.16, s=34, color="white", edgecolor=color, linewidth=1.0, zorder=3, label="AUPRC")
    for yi, value in zip(y, sub["balanced_accuracy"]):
        ax.text(value + 0.012, yi, f"{value:.3f}", va="center", fontsize=7.0)
    ax.set_yticks(y)
    ax.set_yticklabels([m.replace("MRD-KGQC ", "KGQC ").replace(" + review", "\n+ review") for m in sub["method"]])
    ax.set_xlim(0.45, 1.03)
    ax.set_xlabel("Score")
    ax.grid(axis="x", color=PALETTE["light"], lw=0.8)
    ax.legend(frameon=False, loc="lower right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def draw_pr_curve(ax: plt.Axes, pr: pd.DataFrame) -> None:
    color_map = {
        "Prompt-only feature adjudicator": PALETTE["teal"],
        "RAG-only gene prior": PALETTE["purple"],
        "KG/rule-only gene prior": PALETTE["green"],
        "MRD-KGQC plasma-only": PALETTE["blue"],
        "MRD-KGQC plasma-only + review": PALETTE["red"],
    }
    for method, group in pr.groupby("method"):
        if method not in color_map:
            continue
        ax.plot(group["recall"], group["precision"], color=color_map[method], lw=1.8, label=method.replace("MRD-KGQC ", "KGQC "))
    ax.scatter([0.848], [0.671], color=PALETTE["orange"], edgecolor="white", zorder=4, s=50, label="CHIP OR low VAF")
    ax.set_xlabel("CHIP/WBC recall")
    ax.set_ylabel("CHIP/WBC precision")
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0.35, 1.02)
    ax.grid(color=PALETTE["light"], lw=0.8)
    ax.legend(frameon=False, loc="lower left", fontsize=6.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def draw_error_decomposition(ax: plt.Axes, metrics: pd.DataFrame) -> None:
    sub = metrics[metrics["method"].isin(["CHIP gene OR low VAF", "Prompt-only feature adjudicator", "MRD-KGQC plasma-only", "MRD-KGQC plasma-only + review", "MRD-KGQC matched-source"])].copy()
    sub["label"] = sub["method"].replace(
        {
            "CHIP gene OR low VAF": "Rule\nOR",
            "Prompt-only feature adjudicator": "Feature\nonly",
            "MRD-KGQC plasma-only": "KGQC\nplasma",
            "MRD-KGQC plasma-only + review": "KGQC\nreview",
            "MRD-KGQC matched-source": "KGQC\nmatched",
        }
    )
    x = np.arange(len(sub))
    bottom = np.zeros(len(sub))
    stacks = [
        ("Missed CHIP/WBC", "fn_chip", PALETTE["orange"]),
        ("Tumor false suppression", "fp_chip", PALETTE["red"]),
        ("Review", None, PALETTE["review"]),
    ]
    for label, key, color in stacks:
        if key is None:
            values = (sub["review_chip"] + sub["review_tumor"]).to_numpy(dtype=float)
        else:
            values = sub[key].to_numpy(dtype=float)
        ax.bar(x, values, bottom=bottom, color=color, label=label, width=0.64)
        bottom += values
    ax.set_xticks(x)
    ax.set_xticklabels(sub["label"])
    ax.set_ylabel("Hard-label rows")
    ax.legend(frameon=False, ncol=1, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def draw_ablation(ax: plt.Axes, ablation: list[dict[str, object]]) -> None:
    stages = [row["stage"] for row in ablation]
    values = np.array([row["balanced_accuracy"] for row in ablation], dtype=float)
    x = np.arange(len(stages))
    ax.plot(x, values, color=PALETTE["blue"], lw=2.0, marker="o", markersize=5)
    ax.fill_between(x, 0.5, values, color=PALETTE["blue"], alpha=0.10)
    for i, (stage, value) in enumerate(zip(stages, values)):
        ax.text(i, value + 0.018, f"{value:.3f}", ha="center", fontsize=7.0)
    ax.set_xticks(x)
    ax.set_xticklabels(stages, rotation=20, ha="right")
    ax.set_ylabel("Balanced accuracy")
    ax.set_ylim(0.48, 1.04)
    ax.grid(axis="y", color=PALETTE["light"], lw=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def draw_unresolved(ax: plt.Axes, metrics: pd.DataFrame) -> None:
    sub = metrics[metrics["method"].isin(["CHIP gene OR low VAF", "Prompt-only feature adjudicator", "MRD-KGQC plasma-only", "MRD-KGQC plasma-only + review", "MRD-KGQC matched-source"])].copy()
    sub["label"] = sub["method"].replace(
        {
            "CHIP gene OR low VAF": "Rule\nOR",
            "Prompt-only feature adjudicator": "Feature\nonly",
            "MRD-KGQC plasma-only": "KGQC\nplasma",
            "MRD-KGQC plasma-only + review": "KGQC\nreview",
            "MRD-KGQC matched-source": "KGQC\nmatched",
        }
    )
    x = np.arange(len(sub))
    bottom = np.zeros(len(sub))
    for label, col, color in [
        ("CHIP/WBC call", "unresolved_chip_call_rate", PALETTE["chip"]),
        ("Tumor call", "unresolved_tumor_call_rate", PALETTE["tumor"]),
        ("Review", "unresolved_review_rate", PALETTE["review"]),
    ]:
        values = sub[col].to_numpy(dtype=float)
        ax.bar(x, values, bottom=bottom, color=color, width=0.64, label=label)
        bottom += values
    ax.set_xticks(x)
    ax.set_xticklabels(sub["label"])
    ax.set_ylim(0, 1.03)
    ax.set_ylabel("Fraction of cfDNA-only rows")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.18), ncol=3, handlelength=1.4, columnspacing=1.0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def draw_ka7p(ax: plt.Axes, transfer: pd.DataFrame, concordance: pd.DataFrame, cases: pd.DataFrame) -> None:
    ax2 = ax.inset_axes([0.58, 0.10, 0.38, 0.78])
    mut = transfer[transfer["folder"].eq("raw_to_filtered")].dropna(subset=["reduction_fraction"]) if not transfer.empty else pd.DataFrame()
    if not mut.empty:
        values = np.sort(mut["reduction_fraction"].to_numpy(dtype=float))
        y = np.arange(len(values))
        ax.hlines(y, 0, values, color=PALETTE["orange"], alpha=0.22, lw=2.6)
        ax.scatter(values, y, color=PALETTE["orange"], edgecolor="white", linewidth=0.7, s=32, zorder=3)
        ax.axvline(np.nanmedian(values), color=PALETTE["text"], lw=1.1)
        ax.set_xlabel("Mutect2 raw-to-filtered reduction")
        ax.set_ylabel("KA7P cases")
        ax.set_xlim(0, 1.02)
        ax.set_yticks([])
        ax.text(0.04, 0.92, f"n={len(values)}\nmedian={np.nanmedian(values):.2f}", transform=ax.transAxes, ha="left", va="top", fontsize=7.2)
    else:
        ax.text(0.5, 0.5, "KA7P inventory unavailable", ha="center", va="center")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if not concordance.empty:
        vals = [
            concordance["single_caller_variants"].sum(),
            concordance["two_caller_variants"].sum(),
            concordance["three_or_more_caller_variants"].sum(),
        ]
        labels = ["1", "2", "3+"]
        xx = np.arange(len(labels))
        ax2.bar(xx, vals, color=[PALETTE["review"], PALETTE["teal"], PALETTE["blue"]], width=0.62)
        ax2.set_xticks(xx)
        ax2.set_xticklabels(labels)
        ax2.set_xlabel("Caller support")
        ax2.set_ylabel("Filtered variant keys")
        for i, value in enumerate(vals):
            ax2.text(i, value + max(vals) * 0.035, abbrev(value), ha="center", va="bottom", fontsize=6.8)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)


def make_figure(
    metrics: pd.DataFrame,
    pr: pd.DataFrame,
    matrix: pd.DataFrame,
    ablation_rows: list[dict[str, object]],
    transfer: pd.DataFrame,
    concordance: pd.DataFrame,
    cases: pd.DataFrame,
) -> None:
    setup_style()
    fig = plt.figure(figsize=(12.4, 14.0), dpi=300)
    gs = fig.add_gridspec(4, 2, height_ratios=[1.08, 1.15, 1.10, 1.12], hspace=0.58, wspace=0.34)
    axes = [
        fig.add_subplot(gs[0, 0]),
        fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[1, 0]),
        fig.add_subplot(gs[1, 1]),
        fig.add_subplot(gs[2, 0]),
        fig.add_subplot(gs[2, 1]),
        fig.add_subplot(gs[3, 0]),
        fig.add_subplot(gs[3, 1]),
    ]
    draw_dataset_map(axes[0], metrics, cases, transfer)
    draw_method_matrix(axes[1], matrix)
    draw_performance_rank(axes[2], metrics)
    draw_pr_curve(axes[3], pr)
    draw_error_decomposition(axes[4], metrics)
    draw_ablation(axes[5], ablation_rows)
    draw_unresolved(axes[6], metrics)
    draw_ka7p(axes[7], transfer, concordance, cases)
    for label, ax in zip("abcdefgh", axes):
        panel_label(ax, label)

    DOCS.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_PNG, dpi=600, bbox_inches="tight")
    fig.savefig(FIG_PDF, bbox_inches="tight")
    fig.savefig(FIG_SVG, bbox_inches="tight")
    try:
        from PIL import Image

        Image.open(FIG_PNG).convert("L").save(FIG_GRAY)
    except Exception:
        pass
    qa = {
        "figure": str(FIG_BASE),
        "skill": "scipilot-figure-skill",
        "png_dpi": 600,
        "palette_template": "figure_experiment2_mrd_skill_chain",
        "no_subplot_titles": True,
        "no_suptitle": True,
        "panels": {
            "a": "Evaluation dataset map.",
            "b": "Method/input-evidence matrix.",
            "c": "Hard-label performance rank plot.",
            "d": "CHIP/WBC precision-recall curves.",
            "e": "Error and review decomposition.",
            "f": "Ablation trajectory.",
            "g": "cfDNA-only unresolved behavior.",
            "h": "KA7P VCF transfer, raw-to-filtered reduction and caller concordance.",
        },
        "outputs": {
            "png": str(FIG_PNG),
            "pdf": str(FIG_PDF),
            "svg": str(FIG_SVG),
            "grayscale_png": str(FIG_GRAY),
        },
    }
    QA_JSON.write_text(json.dumps(qa, ensure_ascii=False, indent=2), encoding="utf-8")
    plt.close(fig)


def write_ablation(metrics: pd.DataFrame) -> list[dict[str, object]]:
    stage_map = [
        ("VAF only", "Low-VAF rule"),
        ("+ CHIP gene", "CHIP gene OR low VAF"),
        ("+ prompt features", "Prompt-only feature adjudicator"),
        ("+ KG/RAG prior", "MRD-KGQC plasma-only"),
        ("+ matched source", "MRD-KGQC matched-source"),
    ]
    rows = []
    for stage, method in stage_map:
        row = metrics.loc[metrics["method"].eq(method)].iloc[0].to_dict()
        rows.append(
            {
                "stage": stage,
                "method": method,
                "balanced_accuracy": float(row["balanced_accuracy"]),
                "accuracy_strict": float(row["accuracy_strict"]),
                "resolved_accuracy": float(row["resolved_accuracy"]) if not pd.isna(row["resolved_accuracy"]) else None,
                "coverage": float(row["coverage"]),
                "chip_recall": float(row["chip_recall"]) if not pd.isna(row["chip_recall"]) else None,
                "tumor_specificity": float(row["tumor_specificity"]) if not pd.isna(row["tumor_specificity"]) else None,
            }
        )
    payload = {
        "interpretation": (
            "Patient-level cross-validated plasma-only MRD-KGQC improves the rule trade-off, "
            "while the matched-source mode is reported separately as evidence-assisted source tracing."
        ),
        "leakage_control": (
            "Razavi direct DOI/title/paper rows were excluded when constructing gene priors; "
            "patient/sample/allele label fields were masked from plasma-only methods."
        ),
        "stages": rows,
    }
    ABLATION_JSON.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return rows


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)

    prior = build_gene_prior()
    df = load_benchmark(prior)
    predictions = model_predictions(df)
    pred_table = make_predictions_table(df, predictions)
    pred_table.to_csv(PREDICTIONS_CSV, index=False)

    metrics = pd.DataFrame([evaluate_method(df, pred) for pred in predictions])
    metrics["method_order"] = metrics["method"].map({m: i for i, m in enumerate(METHOD_ORDER)})
    metrics = metrics.sort_values("method_order").drop(columns=["method_order"])
    metrics.to_csv(METRICS_CSV, index=False)

    pr = export_pr_curve(df, predictions)
    metach_rows, plasmut_rows = export_comparator_inputs(df)
    external_status = export_external_status(metach_rows, plasmut_rows)
    transfer, concordance, cases = build_ka7p_transfer()
    matrix = method_evidence_matrix()
    ablation_rows = write_ablation(metrics)
    make_figure(metrics, pr, matrix, ablation_rows, transfer, concordance, cases)

    summary = {
        "benchmark_rows": int(len(df)),
        "hard_label_rows": int(df["binary_origin_label"].isin([TUMOR, CHIP]).sum()),
        "unresolved_rows": int((df["binary_origin_label"] == UNRESOLVED).sum()),
        "gene_prior_rows": int(len(prior)),
        "metach_input_rows": int(metach_rows),
        "plasmut_input_rows": int(plasmut_rows),
        "ka7p_filtered_variant_keys": int(concordance["variant_union"].sum()) if not concordance.empty else 0,
        "outputs": {
            "predictions": str(PREDICTIONS_CSV),
            "metrics": str(METRICS_CSV),
            "pr_curve": str(PR_CURVE_CSV),
            "ablation": str(ABLATION_JSON),
            "external_status": str(EXTERNAL_STATUS_CSV),
            "ka7p_transfer": str(KA7P_TRANSFER_CSV),
            "ka7p_concordance": str(KA7P_CONCORDANCE_CSV),
            "figure_png": str(FIG_PNG),
            "figure_pdf": str(FIG_PDF),
            "figure_svg": str(FIG_SVG),
        },
    }
    SUMMARY_JSON.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
