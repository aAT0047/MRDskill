#!/usr/bin/env python3
"""Plot Experiment 2: MRD evidence skill-chain and source-tracing index.

The figure is organized as an eight-panel, publication-style composite with
non-grid visual encodings:

a, source-origin evidence fan
b, source-origin to confidence alluvial flow
c, source-type by confidence bubble plot
d, CHIP/WBC-versus-tumor butterfly volcano
e, multi-skill input/output association fan
f, cancer-type by evidence-type association bubble map
g, gene-pair association arc diagram
h, pipeline output ladder across build, update, index, and demo stages

All subplot titles are omitted by design.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import cm, colormaps
from matplotlib import colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, PathPatch, Rectangle, Wedge
from matplotlib.path import Path as MplPath

from mrdskill_project_paths import artifact_root


ARTIFACTS = artifact_root()
DATA = ARTIFACTS / "data" / "benchmarks"
OUT = ARTIFACTS / "docs" / "benchmarks"

FIG_PNG = OUT / "figure_experiment2_mrd_skill_chain.png"
FIG_PDF = OUT / "figure_experiment2_mrd_skill_chain.pdf"
FIG_SVG = OUT / "figure_experiment2_mrd_skill_chain.svg"
FIG_GRAY = OUT / "figure_experiment2_mrd_skill_chain_grayscale.png"
QA_JSON = OUT / "figure_experiment2_mrd_skill_chain_scipilot_qa.json"


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8,
            "axes.labelsize": 8,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "axes.unicode_minus": False,
        }
    )


def abbrev(value: float | int | str | None) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    if isinstance(value, str):
        return value
    v = float(value)
    if v == 0:
        return "0"
    if v < 0:
        return f"-{abbrev(abs(v))}"
    if v < 10:
        return f"{v:.2f}".rstrip("0").rstrip(".")
    if v < 100:
        return f"{v:.1f}".rstrip("0").rstrip(".")
    if v < 1000:
        return f"{int(round(v)):,}"
    if v < 1_000_000:
        return f"{v/1000:.1f}k".replace(".0k", "k")
    return f"{v/1_000_000:.2f}M".replace(".00M", "M").replace(".0M", "M")


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.10, 1.05, label, transform=ax.transAxes, fontsize=16, fontweight="bold", va="bottom", ha="left")


def shared_heatmap_cmap() -> tuple[mcolors.Colormap, mcolors.LogNorm]:
    cmap = colormaps.get_cmap("cividis").copy()
    cmap.set_bad("#f1f3f5")
    norm = mcolors.LogNorm(vmin=1, vmax=1_000_000)
    return cmap, norm


def prepare_heatmap(values: np.ndarray) -> np.ndarray:
    arr = values.astype(float).copy()
    arr[arr <= 0] = np.nan
    return arr


def draw_heatmap(
    ax: plt.Axes,
    values: np.ndarray,
    row_labels: Sequence[str],
    col_labels: Sequence[str],
    *,
    xlabel: str,
    ylabel: str,
    cmap: mcolors.Colormap,
    norm: mcolors.LogNorm,
    annotate_zero: bool = True,
) -> None:
    show = prepare_heatmap(values)
    im = ax.imshow(show, cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=25, ha="right")
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.tick_params(length=0)

    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            raw = values[i, j]
            if np.isnan(raw):
                continue
            txt = "0" if raw == 0 and annotate_zero else abbrev(raw)
            color = "#1f2937"
            if raw > 0:
                normed = norm(raw)
                color = "white" if normed > 0.62 else "#1f2937"
            ax.text(j, i, txt, ha="center", va="center", fontsize=6.7, color=color)
    return im


def draw_heatmap_simple(
    ax: plt.Axes,
    values: np.ndarray,
    row_labels: Sequence[str],
    col_labels: Sequence[str],
    *,
    xlabel: str,
    ylabel: str,
    cmap: mcolors.Colormap,
    norm: mcolors.LogNorm,
) -> None:
    im = draw_heatmap(ax, values, row_labels, col_labels, xlabel=xlabel, ylabel=ylabel, cmap=cmap, norm=norm)
    return im


def load_tables() -> dict[str, pd.DataFrame]:
    return {
        "landscape": pd.read_csv(DATA / "experiment2_evidence_landscape_table.csv"),
        "origin_evidence": pd.read_csv(DATA / "experiment2_origin_evidence_matrix.csv"),
        "origin_confidence": pd.read_csv(DATA / "experiment2_origin_confidence_matrix.csv"),
        "source_type_confidence": pd.read_csv(DATA / "experiment2_source_type_confidence_matrix.csv"),
        "volcano": pd.read_csv(DATA / "experiment2_gene_origin_enrichment_volcano.csv"),
        "cancer_counts": pd.read_csv(DATA / "experiment2_cancer_evidence_association_counts.csv"),
        "cancer_log2oe": pd.read_csv(DATA / "experiment2_cancer_evidence_association_log2oe.csv"),
        "rule_edges": pd.read_csv(DATA / "experiment2_rule_fan_edges.csv"),
        "skill_chain": pd.read_csv(DATA / "experiment2_skill_chain_table.csv"),
        "metrics": json.loads((DATA / "experiment2_evidence_landscape_metrics.json").read_text(encoding="utf-8")),
        "skill_metrics": json.loads((DATA / "experiment2_skill_chain_metrics.json").read_text(encoding="utf-8")),
    }


def skill_fan(ax: plt.Axes, metrics: dict[str, object]) -> None:
    ax.set_axis_off()
    skill_color = {
        "Harvest": "#4C78A8",
        "Miner": "#F58518",
        "Rater": "#54A24B",
        "KG builder": "#B279A2",
        "GraphRAG": "#E45756",
    }

    skills = [
        ("Harvest", 0.88),
        ("Miner", 0.68),
        ("Rater", 0.48),
        ("KG builder", 0.28),
        ("GraphRAG", 0.10),
    ]
    outputs = [
        ("Files 6,255", 0.90, 6255),
        ("Evidence 61,619", 0.79, 61619),
        ("Source-tracing 100,791", 0.68, 100791),
        ("High conf. 39,111", 0.57, 39111),
        ("KG nodes 173,430", 0.46, 173430),
        ("Rule chunks 575,502", 0.34, 575502),
        ("RAG chunks 728,126", 0.22, 728126),
        ("Rounds 63/63", 0.11, 63),
        ("Hits 126/63", 0.00, 189),
    ]

    edges = [
        ("Harvest", "Files 6,255", 6255),
        ("Harvest", "Evidence 61,619", 61619),
        ("Miner", "Evidence 61,619", 61619),
        ("Miner", "Source-tracing 100,791", 100791),
        ("Rater", "High conf. 39,111", 39111),
        ("KG builder", "KG nodes 173,430", 173430),
        ("KG builder", "Rule chunks 575,502", 575502),
        ("GraphRAG", "RAG chunks 728,126", 728126),
        ("GraphRAG", "Rounds 63/63", 63),
        ("GraphRAG", "Hits 126/63", 189),
    ]

    sx = 0.08
    ox = 0.78
    skill_pos = {name: y for name, y in skills}
    out_pos = {name: (ox, y, count) for name, y, count in outputs}

    for name, y in skills:
        circ = Circle((sx, y), 0.035, facecolor=skill_color[name], edgecolor="none", alpha=0.96)
        ax.add_patch(circ)
        ax.text(sx + 0.045, y, name, va="center", ha="left", fontsize=8.7)

    for name, y, count in outputs:
        circ = Circle((ox, y), 0.030, facecolor="#dfe3ea", edgecolor="#c7ced9", lw=0.7)
        ax.add_patch(circ)
        ax.text(ox - 0.02, y, name, va="center", ha="right", fontsize=8.0)

    strengths = np.array([math.log10(c + 1) for _, _, c in edges], dtype=float)
    smin, smax = strengths.min(), strengths.max()
    for (src, dst, count), strength in zip(edges, strengths):
        y0 = skill_pos[src]
        x0 = sx + 0.035
        y1 = out_pos[dst][1]
        x1 = ox - 0.030
        dx = x1 - x0
        ctrl1 = (x0 + dx * 0.42, y0)
        ctrl2 = (x0 + dx * 0.58, y1)
        path = MplPath(
            [(x0, y0), ctrl1, ctrl2, (x1, y1)],
            [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4],
        )
        color = skill_color[src]
        lw = 0.8 + 5.0 * ((strength - smin) / (smax - smin + 1e-9))
        ax.add_patch(PathPatch(path, facecolor="none", edgecolor=color, lw=lw, alpha=0.25, capstyle="round"))

    ax.set_xlim(0.0, 0.98)
    ax.set_ylim(-0.05, 0.97)
    ax.text(0.04, 0.97, "Skills", ha="left", va="top", fontsize=8.5, color="#6b7280")
    ax.text(0.78, 0.97, "Outputs", ha="right", va="top", fontsize=8.5, color="#6b7280")


ORIGIN_LABELS = [
    "Unresolved",
    "CHIP/WBC",
    "Tumor",
    "Germline",
    "Low support",
    "Artifact",
]

ORIGIN_COLORS = {
    "Unresolved": "#7A7F87",
    "CHIP/WBC": "#D55E00",
    "Tumor": "#0072B2",
    "Germline": "#009E73",
    "Low support": "#CC79A7",
    "Artifact": "#6B7280",
}

EVIDENCE_LABELS = [
    "WBC filter",
    "Serial tracking",
    "Fragment analysis",
    "Tumor tissue",
    "Germline filter",
    "Resistance/action",
    "Not reported",
]

EVIDENCE_COLORS = {
    "WBC filter": "#4C78A8",
    "Serial tracking": "#72B7B2",
    "Fragment analysis": "#F58518",
    "Tumor tissue": "#E45756",
    "Germline filter": "#54A24B",
    "Resistance/action": "#B279A2",
    "Not reported": "#9CA3AF",
}


def radial_evidence_fan(ax: plt.Axes, origin_evidence: pd.DataFrame) -> None:
    ax.set_aspect("equal")
    ax.axis("off")
    mat = origin_evidence.iloc[:, 1:].to_numpy(dtype=float)
    origin_totals = mat.sum(axis=1)
    total = origin_totals.sum()
    start = 92.0
    gap = 2.5
    inner_r0, inner_r1 = 0.34, 0.58
    outer_r0, outer_r1 = 0.62, 0.95

    for i, origin in enumerate(ORIGIN_LABELS):
        span = max(4.0, 360.0 * origin_totals[i] / total - gap)
        theta1, theta2 = start, start + span
        ax.add_patch(
            Wedge(
                (0, 0),
                inner_r1,
                theta1,
                theta2,
                width=inner_r1 - inner_r0,
                facecolor=ORIGIN_COLORS[origin],
                edgecolor="white",
                lw=1.0,
                alpha=0.92,
            )
        )
        sub_start = theta1
        subtotal = mat[i].sum()
        for j, evidence in enumerate(EVIDENCE_LABELS):
            value = mat[i, j]
            if subtotal <= 0 or value <= 0:
                continue
            sub_span = span * value / subtotal
            ax.add_patch(
                Wedge(
                    (0, 0),
                    outer_r1,
                    sub_start,
                    sub_start + sub_span,
                    width=outer_r1 - outer_r0,
                    facecolor=EVIDENCE_COLORS[evidence],
                    edgecolor="white",
                    lw=0.7,
                    alpha=0.82,
                )
            )
            sub_start += sub_span
        if origin in {"Unresolved", "CHIP/WBC", "Tumor"}:
            mid = math.radians((theta1 + theta2) / 2)
            lx, ly = 1.14 * math.cos(mid), 1.14 * math.sin(mid)
            ha = "left" if lx >= 0 else "right"
            ax.text(lx, ly, f"{origin}\n{abbrev(origin_totals[i])}", ha=ha, va="center", fontsize=7.4)
        start = theta2 + gap

    ax.add_patch(Circle((0, 0), inner_r0 - 0.035, facecolor="white", edgecolor="#d1d5db", lw=0.8))
    ax.text(0, 0.02, abbrev(total), ha="center", va="center", fontsize=9, fontweight="bold")
    ax.text(0, -0.09, "origin-evidence\nlinks", ha="center", va="center", fontsize=6.8, color="#6b7280")
    legend_items = [
        ("WBC", "WBC filter"),
        ("Serial", "Serial tracking"),
        ("Frag.", "Fragment analysis"),
        ("Tumor", "Tumor tissue"),
        ("Germline", "Germline filter"),
        ("Resist.", "Resistance/action"),
        ("Unrep.", "Not reported"),
    ]
    for k, (label, key) in enumerate(legend_items):
        x = -1.24 + (k % 4) * 0.64
        y = -1.23 - (k // 4) * 0.13
        ax.add_patch(Rectangle((x, y), 0.06, 0.06, facecolor=EVIDENCE_COLORS[key], edgecolor="none"))
        ax.text(x + 0.08, y + 0.03, label, va="center", ha="left", fontsize=6.6)
    ax.set_xlim(-1.35, 1.35)
    ax.set_ylim(-1.35, 1.25)


def alluvial_confidence_panel(ax: plt.Axes, origin_conf: pd.DataFrame) -> None:
    ax.axis("off")
    mat = origin_conf.iloc[:, 1:].to_numpy(dtype=float)
    conf_labels = ["Low", "Medium", "High"]
    x_left, x_right = 0.18, 0.82
    scale = 0.72 / mat.sum()
    y0_left = 0.86
    left_segments = {}
    for i, origin in enumerate(ORIGIN_LABELS):
        height = mat[i].sum() * scale
        left_segments[i] = (y0_left - height, y0_left)
        ax.add_patch(Rectangle((x_left - 0.035, y0_left - height), 0.07, height, facecolor=ORIGIN_COLORS[origin], edgecolor="white", lw=0.7))
        if origin in {"Unresolved", "CHIP/WBC", "Tumor"}:
            ax.text(x_left - 0.055, y0_left - height / 2, origin, ha="right", va="center", fontsize=7.2)
        y0_left -= height + 0.012

    col_totals = mat.sum(axis=0)
    y0_right = 0.86
    right_segments = {}
    for j, conf in enumerate(conf_labels):
        height = col_totals[j] * scale
        right_segments[j] = (y0_right - height, y0_right)
        ax.add_patch(Rectangle((x_right - 0.035, y0_right - height), 0.07, height, facecolor="#d9dde5", edgecolor="white", lw=0.7))
        ax.text(x_right + 0.055, y0_right - height / 2, f"{conf}\n{abbrev(col_totals[j])}", ha="left", va="center", fontsize=7.2)
        y0_right -= height + 0.035

    left_cursor = {i: left_segments[i][0] for i in range(mat.shape[0])}
    right_cursor = {j: right_segments[j][0] for j in range(mat.shape[1])}
    for i, origin in enumerate(ORIGIN_LABELS):
        for j in range(mat.shape[1]):
            value = mat[i, j]
            if value <= 0:
                continue
            h = value * scale
            yl0, yl1 = left_cursor[i], left_cursor[i] + h
            yr0, yr1 = right_cursor[j], right_cursor[j] + h
            left_cursor[i] += h
            right_cursor[j] += h
            c = ORIGIN_COLORS[origin]
            verts = [
                (x_left + 0.035, yl0),
                (0.42, yl0),
                (0.58, yr0),
                (x_right - 0.035, yr0),
                (x_right - 0.035, yr1),
                (0.58, yr1),
                (0.42, yl1),
                (x_left + 0.035, yl1),
                (x_left + 0.035, yl0),
            ]
            codes = [
                MplPath.MOVETO,
                MplPath.CURVE4,
                MplPath.CURVE4,
                MplPath.CURVE4,
                MplPath.LINETO,
                MplPath.CURVE4,
                MplPath.CURVE4,
                MplPath.CURVE4,
                MplPath.CLOSEPOLY,
            ]
            ax.add_patch(PathPatch(MplPath(verts, codes), facecolor=c, edgecolor="none", alpha=0.24))
    ax.text(x_left, 0.97, "Origin", ha="center", va="top", fontsize=8.0, color="#6b7280")
    ax.text(x_right, 0.97, "Confidence", ha="center", va="top", fontsize=8.0, color="#6b7280")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.03, 1.0)


def confidence_bubble_panel(ax: plt.Axes, source_conf: pd.DataFrame) -> None:
    labels = [
        "clinical\nvalidation",
        "discussion/\nbackground",
        "meta\nanalysis",
        "methods/\nassay detail",
        "primary\nresult",
        "supplementary\nresult",
        "uncertain",
    ]
    cols = ["Low", "Medium", "High"]
    mat = source_conf.iloc[:, 1:].to_numpy(dtype=float)
    xs = np.arange(len(cols))
    ys = np.arange(len(labels))
    X, Y = np.meshgrid(xs, ys)
    max_count = mat.max()
    sizes = 25 + 850 * np.sqrt(mat / max_count)
    colors = {"Low": "#9CA3AF", "Medium": "#E69F00", "High": "#009E73"}
    for j, conf in enumerate(cols):
        ax.scatter(
            X[:, j],
            Y[:, j],
            s=sizes[:, j],
            color=colors[conf],
            edgecolors="white",
            linewidths=0.7,
            alpha=0.78,
            label=conf,
        )
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if mat[i, j] >= 8000:
                ax.text(j, i, abbrev(mat[i, j]), ha="center", va="center", fontsize=6.8, color="#111827")
    ax.set_xticks(xs)
    ax.set_xticklabels(cols)
    ax.set_yticks(ys)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Confidence")
    ax.set_ylabel("Evidence source type")
    ax.set_ylim(len(labels) - 0.5, -0.5)
    ax.grid(axis="x", color="#e5e7eb", lw=0.8)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)


def butterfly_volcano_panel(ax: plt.Axes, volcano: pd.DataFrame) -> None:
    df = volcano.copy()
    df["direction"] = np.where(df["log2_odds_ratio_chip_vs_tumor"] >= 0, "CHIP/WBC", "Tumor")
    top_chip = df[df["direction"] == "CHIP/WBC"].nsmallest(8, "q_value").copy()
    top_tumor = df[df["direction"] == "Tumor"].nsmallest(8, "q_value").copy()
    plot_df = pd.concat([top_tumor, top_chip], ignore_index=True)
    plot_df["signed_fdr"] = plot_df["neg_log10_q"].clip(upper=160) * np.where(plot_df["direction"] == "CHIP/WBC", 1, -1)
    plot_df = plot_df.sort_values(["direction", "signed_fdr"])
    y = np.arange(len(plot_df))
    colors = np.where(plot_df["direction"] == "CHIP/WBC", "#D55E00", "#0072B2")
    sizes = 30 + 160 * (np.log10(plot_df["total_count"] + 1) / np.log10(volcano["total_count"].max() + 1))
    ax.axvline(0, color="#9ca3af", lw=1.0)
    for yi, x, c in zip(y, plot_df["signed_fdr"], colors):
        ax.plot([0, x], [yi, yi], color=c, lw=2.0, alpha=0.45)
    ax.scatter(plot_df["signed_fdr"], y, s=sizes, c=colors, edgecolors="white", linewidths=0.7, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["gene"])
    ax.set_xlabel("Signed -log10 FDR (capped at 160)")
    ax.set_xlim(-165, 165)
    ax.set_ylim(-0.8, len(plot_df) - 0.2)
    ax.text(0.02, 0.97, "Tumor", transform=ax.transAxes, ha="left", va="top", color="#0072B2", fontsize=8.0)
    ax.text(0.98, 0.97, "CHIP/WBC", transform=ax.transAxes, ha="right", va="top", color="#D55E00", fontsize=8.0)
    ax.grid(axis="x", color="#e5e7eb", lw=0.8)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)


def volcano_panel(ax: plt.Axes, volcano: pd.DataFrame) -> None:
    df = volcano.copy()
    df["direction"] = np.where(df["log2_odds_ratio_chip_vs_tumor"] >= 0, "CHIP/WBC-enriched", "Tumor-enriched")
    df["sig"] = df["q_value"] < 0.05
    chip_color = "#D55E00"
    tumor_color = "#0072B2"
    neutral = "#b8c0cc"
    colors = np.where(df["sig"], np.where(df["direction"] == "CHIP/WBC-enriched", chip_color, tumor_color), neutral)
    sizes = 18 + 140 * (np.log10(df["total_count"] + 1) / np.log10(df["total_count"].max() + 1))
    ax.scatter(
        df["log2_odds_ratio_chip_vs_tumor"],
        df["neg_log10_q"],
        s=sizes,
        c=colors,
        alpha=0.82,
        edgecolors="white",
        linewidths=0.45,
    )
    ax.axvline(0, color="#9aa3af", lw=0.9)
    ax.axhline(-math.log10(0.05), color="#9aa3af", lw=0.9, ls="--")
    ax.set_xlabel("log2 odds ratio (CHIP/WBC vs tumor)")
    ax.set_ylabel("-log10 FDR")
    ax.set_xlim(min(-6.3, df["log2_odds_ratio_chip_vs_tumor"].min() - 0.5), max(10.5, df["log2_odds_ratio_chip_vs_tumor"].max() + 0.5))
    ax.set_ylim(0, max(70, float(df["neg_log10_q"].max()) + 5))
    ax.text(0.98, 0.96, f"{len(df)} genes\n{(df['q_value'] < 0.05).sum()} q<0.05", transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color="#374151")
    top_chip = df[df["direction"] == "CHIP/WBC-enriched"].nsmallest(4, "q_value")
    top_tumor = df[df["direction"] == "Tumor-enriched"].nsmallest(4, "q_value")
    for _, row in pd.concat([top_chip, top_tumor]).iterrows():
        x = row["log2_odds_ratio_chip_vs_tumor"]
        y = row["neg_log10_q"]
        dx = 0.15 if x >= 0 else -0.15
        dy = 0.7 if y < 60 else -1.3
        ha = "left" if x >= 0 else "right"
        ax.annotate(
            row["gene"],
            (x, y),
            xytext=(x + dx, y + dy),
            textcoords="data",
            fontsize=7.5,
            ha=ha,
            va="center",
            arrowprops=dict(arrowstyle="-", lw=0.5, color="#9aa3af"),
        )
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=chip_color, markeredgecolor="white", markersize=6, label="CHIP/WBC-enriched"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=tumor_color, markeredgecolor="white", markersize=6, label="Tumor-enriched"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=neutral, markeredgecolor="white", markersize=6, label="Not significant"),
    ]
    ax.legend(handles=handles, frameon=False, loc="lower right", ncol=1, handletextpad=0.4, columnspacing=0.8)


def bubble_panel(ax: plt.Axes, counts: pd.DataFrame, log2oe: pd.DataFrame) -> None:
    evidence_order = [
        "matched_wbc_filtering",
        "serial_longitudinal_tracking",
        "fragmentomics_discordance",
        "tumor_tissue_concordance",
        "germline_filtering",
        "target_resistance_actionability",
    ]
    evidence_labels = [
        "WBC\nfilter",
        "Serial\ntracking",
        "Fragment\nanalysis",
        "Tumor\ntissue",
        "Germline\nfilter",
        "Resistance/\naction",
    ]
    cancer_order = [
        "non-small cell lung cancer",
        "breast cancer",
        "colorectal cancer",
        "acute myeloid leukemia",
        "urothelial carcinoma",
        "lymphoma",
        "pancreatic cancer",
        "melanoma",
    ]
    cancer_labels = ["NSCLC", "Breast", "CRC", "AML", "Urothelial", "Lymphoma", "Pancreatic", "Melanoma"]

    counts = counts.set_index("cancer_type").loc[cancer_order, evidence_order]
    log2oe = log2oe.set_index("cancer_type").loc[cancer_order, evidence_order]
    max_count = float(np.nanmax(counts.values))
    xs = np.arange(len(evidence_order))
    ys = np.arange(len(cancer_order))
    X, Y = np.meshgrid(xs, ys)
    sizes = 40 + 620 * (counts.values / max_count)
    sc = ax.scatter(
        X.ravel(),
        Y.ravel(),
        s=sizes.ravel(),
        c=log2oe.values.ravel(),
        cmap="RdBu_r",
        vmin=-3,
        vmax=3,
        edgecolors="white",
        linewidths=0.45,
    )
    ax.set_xticks(xs)
    ax.set_xticklabels(evidence_labels, rotation=0)
    ax.set_yticks(ys)
    ax.set_yticklabels(cancer_labels)
    ax.set_xlabel("Evidence type")
    ax.set_ylabel("Cancer type")
    ax.set_xlim(-0.6, len(evidence_order) - 0.4)
    ax.set_ylim(len(cancer_order) - 0.5, -0.5)
    ax.grid(True, color="#e5e7eb", lw=0.8)
    ax.set_axisbelow(True)
    cb = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label("log2 observed/expected")


def arc_panel(ax: plt.Axes, edges: pd.DataFrame) -> None:
    ax.set_facecolor("white")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color("#d1d5db")
    ax.spines["bottom"].set_linewidth(0.8)

    genes = sorted(
        set(edges["u"]).union(edges["v"]),
        key=lambda g: (-float(edges[(edges["u"] == g) | (edges["v"] == g)]["score"].sum()), g),
    )
    positions = np.linspace(0.08, 0.92, len(genes))
    x_map = {gene: x for gene, x in zip(genes, positions)}

    subset = pd.concat(
        [
            edges[edges["rule_class"] == "cooccurrence"].head(7),
            edges[edges["rule_class"] == "mutual_exclusion"].head(7),
        ],
        ignore_index=True,
    )
    scores = np.log10(subset["score"].astype(float) + 1)
    smin, smax = float(scores.min()), float(scores.max())
    colors = {"cooccurrence": "#E67E22", "mutual_exclusion": "#2E86DE"}
    baseline = 0.0
    ax.axhline(baseline, color="#9ca3af", lw=0.8)

    for _, row in subset.iterrows():
        x0 = x_map[row["u"]]
        x1 = x_map[row["v"]]
        if x0 > x1:
            x0, x1 = x1, x0
        dx = x1 - x0
        score = float(np.log10(float(row["score"]) + 1))
        height = 0.11 + 0.24 * ((score - smin) / (smax - smin + 1e-9))
        sign = 1 if row["rule_class"] == "cooccurrence" else -1
        y = baseline + sign * height
        c1 = (x0 + dx * 0.20, y)
        c2 = (x1 - dx * 0.20, y)
        path = MplPath(
            [(x0, baseline), c1, c2, (x1, baseline)],
            [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4],
        )
        lw = 0.7 + 3.9 * ((score - smin) / (smax - smin + 1e-9))
        ax.add_patch(PathPatch(path, facecolor="none", edgecolor=colors[row["rule_class"]], lw=lw, alpha=0.28))

    for gene in genes:
        x = x_map[gene]
        ax.scatter([x], [baseline], s=18, color="#111827", zorder=3)
        ax.text(x, -0.045, gene, ha="center", va="top", rotation=35, fontsize=7.5)

    legend_handles = [
        Line2D([0], [0], color=colors["cooccurrence"], lw=2.0, label="Co-occurrence"),
        Line2D([0], [0], color=colors["mutual_exclusion"], lw=2.0, label="Mutual exclusion"),
    ]
    ax.legend(handles=legend_handles, frameon=False, loc="upper right", ncol=1, handlelength=1.8)
    ax.set_xlim(0.03, 0.97)
    ax.set_ylim(-0.38, 0.42)
    ax.set_xticks([])
    ax.set_yticks([])


def pipeline_panel(ax: plt.Axes, metrics: dict[str, object]) -> None:
    stage_cols = ["Full build", "Incremental update", "Latest index", "Demo"]
    row_labels = [
        "Source files\nprocessed",
        "Scored evidence\nrows",
        "Original graph\nnodes",
        "Original graph\nrelationships",
        "Incremental source\ndocuments",
        "Incremental scored\nevidence rows",
        "Augmented KG\nnode IDs",
        "Augmented KG\nrelationship IDs",
        "Literature-evidence\nchunks",
        "Materialized\nrule chunks",
        "Total RAG\nchunks",
        "Non-empty retrieval\nrounds",
    ]
    values = np.array(
        [
            [6255, np.nan, np.nan, np.nan],
            [61619, 152633, np.nan, np.nan],
            [173430, np.nan, np.nan, np.nan],
            [904807, np.nan, np.nan, np.nan],
            [np.nan, 21051, np.nan, np.nan],
            [np.nan, 152633, np.nan, np.nan],
            [np.nan, 1133875, np.nan, np.nan],
            [np.nan, 5364370, np.nan, np.nan],
            [np.nan, np.nan, 152622, np.nan],
            [np.nan, np.nan, 575502, np.nan],
            [np.nan, np.nan, 728126, np.nan],
            [np.nan, np.nan, np.nan, 63],
        ],
        dtype=float,
    )
    show = values.copy()
    show[show <= 0] = np.nan
    cmap = colormaps.get_cmap("cividis").copy()
    cmap.set_bad("#f1f3f5")
    norm = mcolors.LogNorm(vmin=1, vmax=5_364_370)
    im = ax.imshow(show, cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
    ax.set_xticks(np.arange(len(stage_cols)))
    ax.set_xticklabels(stage_cols)
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels)
    ax.set_xlabel("Workflow stage")
    ax.set_ylabel("Output layer")
    ax.tick_params(length=0)

    anno = [
        ["6,255", "", "", ""],
        ["61,619", "152,633", "", ""],
        ["173,430", "", "", ""],
        ["904,807", "", "", ""],
        ["", "21,051", "", ""],
        ["", "152,633", "", ""],
        ["", "1.13M", "", ""],
        ["", "5.36M", "", ""],
        ["", "", "152,622", ""],
        ["", "", "575,502", ""],
        ["", "", "728,126", ""],
        ["", "", "", "63/63"],
    ]
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            if np.isnan(values[i, j]):
                continue
            color = "white" if norm(values[i, j]) > 0.62 else "#1f2937"
            ax.text(j, i, anno[i][j], ha="center", va="center", fontsize=6.7, color=color)
    ax.text(0.99, -0.18, "Query demo returned 126 literature hits and 63 rule hits", transform=ax.transAxes, ha="right", va="top", fontsize=7.2, color="#6b7280")
    return im


def output_ladder_panel(ax: plt.Axes) -> None:
    rows = [
        ("Full build", "Source files", 6_255),
        ("Full build", "Scored evidence", 61_619),
        ("Full build", "Graph nodes", 173_430),
        ("Full build", "Graph relationships", 904_807),
        ("Incremental update", "Source documents", 21_051),
        ("Incremental update", "Scored evidence", 152_633),
        ("Incremental update", "Augmented node IDs", 1_133_875),
        ("Incremental update", "Augmented relationship IDs", 5_364_370),
        ("Latest index", "Literature-evidence chunks", 152_622),
        ("Latest index", "Rule chunks", 575_502),
        ("Latest index", "Total RAG chunks", 728_126),
        ("Source-tracing demo", "Retrieval rounds", 63),
        ("Source-tracing demo", "Literature hits", 126),
        ("Source-tracing demo", "Rule hits", 63),
    ]
    stage_colors = {
        "Full build": "#4C78A8",
        "Incremental update": "#F58518",
        "Latest index": "#54A24B",
        "Source-tracing demo": "#E45756",
    }
    y = np.arange(len(rows))[::-1]
    x = np.log10([r[2] for r in rows])
    labels = [r[1] for r in rows]
    stages = [r[0] for r in rows]
    counts = [r[2] for r in rows]
    for yi, xi, stage, count in zip(y, x, stages, counts):
        ax.plot([0, xi], [yi, yi], color=stage_colors[stage], alpha=0.23, lw=3)
        ax.scatter([xi], [yi], s=42 + 110 * (xi / max(x)), color=stage_colors[stage], edgecolor="white", linewidth=0.8, zorder=3)
        ax.text(xi + 0.07, yi, abbrev(count), va="center", ha="left", fontsize=7.2)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("log10 count")
    ax.set_xlim(0, max(x) + 0.65)
    ax.grid(axis="x", color="#e5e7eb", lw=0.8)
    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=c, markeredgecolor="white", markersize=7, label=s) for s, c in stage_colors.items()]
    ax.legend(handles=handles, frameon=False, ncol=2, loc="lower right", handletextpad=0.3, columnspacing=1.0)
    ax.text(0.01, -0.20, "The source-tracing demo returned evidence in 63/63 retrieval rounds.", transform=ax.transAxes, ha="left", va="top", fontsize=7.2, color="#6b7280")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)


def summarize_and_write(metrics: dict[str, object], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(metrics, ensure_ascii=False, indent=2), encoding="utf-8")


def make_figure(tables: dict[str, pd.DataFrame]) -> None:
    setup_style()

    origin_evidence = tables["origin_evidence"]
    origin_conf = tables["origin_confidence"]
    source_conf = tables["source_type_confidence"]
    volcano = tables["volcano"]
    cancer_counts = tables["cancer_counts"]
    cancer_log2oe = tables["cancer_log2oe"]
    rule_edges = tables["rule_edges"]

    fig = plt.figure(figsize=(12.4, 15.0), dpi=300)
    gs = fig.add_gridspec(
        5,
        2,
        height_ratios=[1.22, 1.10, 1.34, 1.12, 1.24],
        hspace=0.62,
        wspace=0.35,
    )

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])
    ax_e = fig.add_subplot(gs[2, :])
    ax_f = fig.add_subplot(gs[3, 0])
    ax_g = fig.add_subplot(gs[3, 1])
    ax_h = fig.add_subplot(gs[4, :])

    radial_evidence_fan(ax_a, origin_evidence)
    panel_label(ax_a, "a")

    alluvial_confidence_panel(ax_b, origin_conf)
    panel_label(ax_b, "b")

    confidence_bubble_panel(ax_c, source_conf)
    panel_label(ax_c, "c")

    butterfly_volcano_panel(ax_d, volcano)
    panel_label(ax_d, "d")

    skill_fan(ax_e, tables["metrics"])
    panel_label(ax_e, "e")

    bubble_panel(ax_f, cancer_counts, cancer_log2oe)
    panel_label(ax_f, "f")

    arc_panel(ax_g, rule_edges)
    panel_label(ax_g, "g")

    output_ladder_panel(ax_h)
    panel_label(ax_h, "h")

    for ax in [ax_d, ax_f]:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_PNG, dpi=600, bbox_inches="tight")
    fig.savefig(FIG_PDF, bbox_inches="tight")
    fig.savefig(FIG_SVG, bbox_inches="tight")

    try:
        from PIL import Image

        Image.open(FIG_PNG).convert("L").save(FIG_GRAY)
    except Exception:
        pass

    qa = {
        "figure": str(FIG_PNG.with_suffix("")),
        "skill": "scipilot-figure-skill",
        "no_subplot_titles": True,
        "no_suptitle": True,
        "png_dpi": 600,
        "panels": {
            "a": "Source-origin evidence fan.",
            "b": "Source-origin to confidence alluvial flow.",
            "c": "Evidence-source by confidence bubble plot.",
            "d": "CHIP/WBC-versus-tumor butterfly volcano.",
            "e": "Multi-skill input/output association fan.",
            "f": "Cancer-type by evidence-type association bubble map.",
            "g": "Top gene-pair association arc diagram.",
            "h": "Pipeline output ladder across full build, incremental update, latest index, and demo.",
        },
        "outputs": {
            "png": str(FIG_PNG),
            "pdf": str(FIG_PDF),
            "svg": str(FIG_SVG),
            "grayscale_png": str(FIG_GRAY),
        },
    }
    summarize_and_write(qa, QA_JSON)
    plt.close(fig)


def main() -> None:
    tables = load_tables()
    make_figure(tables)
    print(json.dumps(
        {
            "figure_png": str(FIG_PNG),
            "figure_pdf": str(FIG_PDF),
            "figure_svg": str(FIG_SVG),
            "figure_gray": str(FIG_GRAY),
            "qa_json": str(QA_JSON),
        },
        ensure_ascii=False,
        indent=2,
    ))


if __name__ == "__main__":
    main()
