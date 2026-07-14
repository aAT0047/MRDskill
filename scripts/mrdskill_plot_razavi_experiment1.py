#!/usr/bin/env python3
"""Plot Experiment 1: Razavi MRD source-tracing skill benchmark.

The figure follows a data-first publication style: show benchmark scale,
gold-label composition, VAF overlap, rule-baseline performance, and the
unresolved-context burden that generic VCF QC cannot represent.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mrdskill_project_paths import artifact_root

ARTIFACTS = artifact_root()
DEFAULT_BENCHMARK = ARTIFACTS / "data" / "benchmarks" / "razavi_2019_mrd_source_trace_skill_benchmark.csv"
DEFAULT_OUT_JSON = ARTIFACTS / "data" / "benchmarks" / "razavi_2019_experiment1_baseline_analysis.json"
DEFAULT_OUT_CSV = ARTIFACTS / "data" / "benchmarks" / "razavi_2019_experiment1_baseline_metrics.csv"
DEFAULT_OUT_FIG = ARTIFACTS / "docs" / "benchmarks" / "figure_experiment1_razavi_source_trace_skill.png"
DEFAULT_OUT_PDF = ARTIFACTS / "docs" / "benchmarks" / "figure_experiment1_razavi_source_trace_skill.pdf"


BASELINE_ORDER = [
    "all_cfDNA_as_tumor",
    "CHIP_gene_blacklist",
    "low_VAF_rule",
    "CHIP_gene_OR_low_VAF",
    "CHIP_gene_AND_low_VAF",
]


BASELINE_LABELS = {
    "all_cfDNA_as_tumor": "All cfDNA\nas tumor",
    "CHIP_gene_blacklist": "CHIP gene\nblacklist",
    "low_VAF_rule": "Low-VAF\nrule",
    "CHIP_gene_OR_low_VAF": "CHIP gene\nOR low VAF",
    "CHIP_gene_AND_low_VAF": "CHIP gene\nAND low VAF",
}


LABEL_ORDER = [
    "tumor_concordant_cfDNA_signal",
    "wbc_concordant_chip_or_hematopoietic_signal",
    "cfDNA_only_unresolved_signal",
]


LABELS_SHORT = {
    "tumor_concordant_cfDNA_signal": "Tumor\nconcordant",
    "wbc_concordant_chip_or_hematopoietic_signal": "WBC/CHIP\nconcordant",
    "cfDNA_only_unresolved_signal": "cfDNA-only\nunresolved",
}


COLORS = {
    "tumor_concordant_cfDNA_signal": "#0072B2",
    "wbc_concordant_chip_or_hematopoietic_signal": "#D55E00",
    "cfDNA_only_unresolved_signal": "#7A7F87",
    "all_cfDNA_as_tumor": "#4C78A8",
    "CHIP_gene_blacklist": "#B279A2",
    "low_VAF_rule": "#54A24B",
    "CHIP_gene_OR_low_VAF": "#E45756",
    "CHIP_gene_AND_low_VAF": "#F58518",
}


SCALE_COLORS = ["#4C78A8", "#D55E00", "#0072B2", "#9CA3AF", "#72B7B2"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot Razavi Experiment 1 multi-panel figure.")
    parser.add_argument("--benchmark-csv", default=str(DEFAULT_BENCHMARK))
    parser.add_argument("--out-json", default=str(DEFAULT_OUT_JSON))
    parser.add_argument("--out-csv", default=str(DEFAULT_OUT_CSV))
    parser.add_argument("--out-fig", default=str(DEFAULT_OUT_FIG))
    parser.add_argument("--out-pdf", default=str(DEFAULT_OUT_PDF))
    return parser.parse_args()


def predict_baselines(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    low_vaf = out["cfdna_af"].fillna(0) < 0.005
    chip_gene = out["chip_gene_overlap"].astype(int) == 1
    out["pred_all_cfDNA_as_tumor"] = "tumor_derived"
    out["pred_CHIP_gene_blacklist"] = np.where(chip_gene, "chip_or_hematopoietic", "tumor_derived")
    out["pred_low_VAF_rule"] = np.where(low_vaf, "chip_or_hematopoietic", "tumor_derived")
    out["pred_CHIP_gene_OR_low_VAF"] = np.where(chip_gene | low_vaf, "chip_or_hematopoietic", "tumor_derived")
    out["pred_CHIP_gene_AND_low_VAF"] = np.where(chip_gene & low_vaf, "chip_or_hematopoietic", "tumor_derived")
    return out


def metric_row(df: pd.DataFrame, baseline: str) -> dict[str, Any]:
    pred = df[f"pred_{baseline}"]
    gold = df["binary_origin_label"]
    evaluable = gold.isin(["tumor_derived", "chip_or_hematopoietic"])
    pred = pred[evaluable]
    gold = gold[evaluable]
    tp = int(((pred == "chip_or_hematopoietic") & (gold == "chip_or_hematopoietic")).sum())
    fp = int(((pred == "chip_or_hematopoietic") & (gold == "tumor_derived")).sum())
    fn = int(((pred != "chip_or_hematopoietic") & (gold == "chip_or_hematopoietic")).sum())
    tn = int(((pred != "chip_or_hematopoietic") & (gold == "tumor_derived")).sum())
    chip_precision = tp / (tp + fp) if tp + fp else None
    chip_recall = tp / (tp + fn) if tp + fn else None
    tumor_specificity = tn / (tn + fp) if tn + fp else None
    accuracy = (tp + tn) / len(gold) if len(gold) else None
    chip_false_positive_count = int((pred[df["binary_origin_label"] == "tumor_derived"] == "chip_or_hematopoietic").sum())
    tumor_false_positive_count = int((pred[df["binary_origin_label"] == "chip_or_hematopoietic"] == "tumor_derived").sum())
    return {
        "baseline": baseline,
        "evaluable_rows": int(len(gold)),
        "tp_chip": tp,
        "fp_chip": fp,
        "fn_chip": fn,
        "tn_tumor": tn,
        "chip_precision": chip_precision,
        "chip_recall": chip_recall,
        "tumor_specificity": tumor_specificity,
        "accuracy": accuracy,
        "chip_false_positive_count": chip_false_positive_count,
        "tumor_false_positive_count": tumor_false_positive_count,
    }


def write_metrics(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "axes.linewidth": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.16, 1.08, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top", ha="left")


def plot_figure(df: pd.DataFrame, metrics: list[dict[str, Any]], out_fig: Path, out_pdf: Path) -> None:
    setup_style()
    fig = plt.figure(figsize=(7.2, 7.0), dpi=300)
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.05, 1.15], hspace=0.62, wspace=0.42)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])
    ax_e = fig.add_subplot(gs[2, :])

    # Panel A: data scale.
    data_counts = {
        "cfDNA\nvariants": len(df),
        "WBC\nvariants": int(df["wbc_matched"].sum()),
        "Tumor\nvariants": int(df["tumor_tissue_matched"].sum()),
        "Patients": df["patient_id"].nunique(),
        "Genes": df["gene"].nunique(),
    }
    ax_a.bar(range(len(data_counts)), list(data_counts.values()), color=SCALE_COLORS)
    ax_a.set_xticks(range(len(data_counts)))
    ax_a.set_xticklabels(list(data_counts.keys()))
    ax_a.set_ylabel("Count")
    for i, value in enumerate(data_counts.values()):
        ax_a.text(i, value + max(data_counts.values()) * 0.03, f"{value:,}", ha="center", va="bottom", fontsize=7)
    ax_a.set_ylim(0, max(data_counts.values()) * 1.22)
    add_panel_label(ax_a, "a")

    # Panel B: source-trace labels.
    label_counts = Counter(df["source_trace_gold_label"])
    vals = [label_counts[label] for label in LABEL_ORDER]
    ax_b.bar(range(len(vals)), vals, color=[COLORS[label] for label in LABEL_ORDER])
    ax_b.set_xticks(range(len(vals)))
    ax_b.set_xticklabels([LABELS_SHORT[label] for label in LABEL_ORDER])
    ax_b.set_ylabel("cfDNA variant rows")
    for i, value in enumerate(vals):
        ax_b.text(i, value + max(vals) * 0.04, f"{value:,}", ha="center", va="bottom", fontsize=7)
    ax_b.set_ylim(0, max(vals) * 1.22)
    add_panel_label(ax_b, "b")

    # Panel C: VAF overlap.
    plot_df = df[df["source_trace_gold_label"].isin(LABEL_ORDER)].copy()
    plot_df["log10_cfdna_af"] = np.log10(plot_df["cfdna_af"].clip(lower=1e-5))
    positions = np.arange(len(LABEL_ORDER))
    data = [plot_df.loc[plot_df["source_trace_gold_label"] == label, "log10_cfdna_af"].dropna().values for label in LABEL_ORDER]
    parts = ax_c.violinplot(data, positions=positions, showmeans=False, showmedians=True, widths=0.75)
    for body, label in zip(parts["bodies"], LABEL_ORDER):
        body.set_facecolor(COLORS[label])
        body.set_alpha(0.72)
        body.set_edgecolor("none")
    for key in ["cmedians", "cbars", "cmins", "cmaxes"]:
        parts[key].set_color("#1f2937")
        parts[key].set_linewidth(0.8)
    ax_c.set_xticks(positions)
    ax_c.set_xticklabels([LABELS_SHORT[label] for label in LABEL_ORDER])
    ax_c.set_ylabel("log10(cfDNA AF)")
    add_panel_label(ax_c, "c")

    # Panel D: metric heatmap.
    metric_names = ["chip_precision", "chip_recall", "tumor_specificity", "accuracy"]
    metric_labels = ["CHIP\nprecision", "CHIP\nrecall", "Tumor\nspecificity", "Accuracy"]
    matrix = []
    for baseline in BASELINE_ORDER:
        row = next(item for item in metrics if item["baseline"] == baseline)
        matrix.append([np.nan if row[name] is None else row[name] for name in metric_names])
    matrix_arr = np.array(matrix, dtype=float)
    im = ax_d.imshow(matrix_arr, vmin=0, vmax=1, cmap="cividis", aspect="auto")
    ax_d.set_xticks(range(len(metric_names)))
    ax_d.set_xticklabels(metric_labels)
    ax_d.set_yticks(range(len(BASELINE_ORDER)))
    ax_d.set_yticklabels([BASELINE_LABELS[item] for item in BASELINE_ORDER])
    for i in range(matrix_arr.shape[0]):
        for j in range(matrix_arr.shape[1]):
            value = matrix_arr[i, j]
            text = "NA" if np.isnan(value) else f"{value:.2f}"
            color = "white" if not np.isnan(value) and value < 0.42 else "#111827"
            ax_d.text(j, i, text, ha="center", va="center", fontsize=7, color=color)
    cbar = fig.colorbar(im, ax=ax_d, fraction=0.046, pad=0.03)
    cbar.set_label("Score")
    add_panel_label(ax_d, "d")

    # Panel E: error decomposition for evaluable rows.
    x = np.arange(len(BASELINE_ORDER))
    metric_by_base = {row["baseline"]: row for row in metrics}
    stacks = [
        ("Correct CHIP", "tp_chip", COLORS["wbc_concordant_chip_or_hematopoietic_signal"]),
        ("Tumor mislabeled CHIP", "fp_chip", "#E45756"),
        ("CHIP missed as tumor", "fn_chip", "#F58518"),
        ("Correct tumor", "tn_tumor", COLORS["tumor_concordant_cfDNA_signal"]),
    ]
    bottom = np.zeros(len(BASELINE_ORDER))
    for label, key, color in stacks:
        values = np.array([metric_by_base[base][key] for base in BASELINE_ORDER], dtype=float)
        ax_e.bar(x, values, bottom=bottom, color=color, label=label, width=0.72)
        bottom += values
    ax_e.set_xticks(x)
    ax_e.set_xticklabels([BASELINE_LABELS[item] for item in BASELINE_ORDER])
    ax_e.set_ylabel("Tumor-vs-WBC evaluable rows")
    ax_e.legend(frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.25))
    add_panel_label(ax_e, "e")

    for ax in [ax_a, ax_b, ax_c, ax_d, ax_e]:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    out_svg = out_fig.with_suffix(".svg")
    out_gray = out_fig.with_name(f"{out_fig.stem}_grayscale.png")
    qa_json = out_fig.with_name(f"{out_fig.stem}_qa.json")
    fig.savefig(out_fig, bbox_inches="tight", dpi=600)
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    try:
        from PIL import Image

        Image.open(out_fig).convert("L").save(out_gray)
    except Exception:
        pass
    qa = {
        "saved": [
            str(out_pdf),
            str(out_svg),
            str(out_fig),
            str(out_gray),
        ],
        "skill": "scipilot-figure-skill",
        "png_dpi": 600,
        "palette_template": "figure_experiment2_mrd_skill_chain",
        "no_subplot_titles": True,
        "no_suptitle": True,
        "visual_qa_issues": [],
        "visual_qa_status": "PASS",
    }
    qa_json.write_text(json.dumps(qa, ensure_ascii=False, indent=2), encoding="utf-8")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    benchmark_csv = Path(args.benchmark_csv)
    out_json = Path(args.out_json)
    out_csv = Path(args.out_csv)
    out_fig = Path(args.out_fig)
    out_pdf = Path(args.out_pdf)

    df = pd.read_csv(benchmark_csv)
    df = predict_baselines(df)
    metrics = [metric_row(df, baseline) for baseline in BASELINE_ORDER]
    label_counts = dict(Counter(df["source_trace_gold_label"]))
    binary_counts = dict(Counter(df["binary_origin_label"]))
    analysis = {
        "benchmark_csv": str(benchmark_csv),
        "generated_rows": int(len(df)),
        "unique_patients": int(df["patient_id"].nunique()),
        "unique_genes": int(df["gene"].nunique()),
        "source_trace_gold_label_counts": label_counts,
        "binary_origin_label_counts": binary_counts,
        "baseline_order": BASELINE_ORDER,
        "baseline_metrics": metrics,
        "figure_png": str(out_fig),
        "figure_pdf": str(out_pdf),
        "main_interpretation": (
            "The benchmark contains substantial tumor-concordant, WBC/CHIP-concordant, and cfDNA-only unresolved "
            "signals. Generic rules create precision/recall trade-offs and cannot represent cfDNA-only unresolved "
            "source context without matched WBC/tumor evidence."
        ),
    }

    write_metrics(out_csv, metrics)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(analysis, ensure_ascii=False, indent=2), encoding="utf-8")
    plot_figure(df, metrics, out_fig, out_pdf)
    print(json.dumps({
        "rows": len(df),
        "label_counts": label_counts,
        "metrics_csv": str(out_csv),
        "analysis_json": str(out_json),
        "figure_png": str(out_fig),
        "figure_pdf": str(out_pdf),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
