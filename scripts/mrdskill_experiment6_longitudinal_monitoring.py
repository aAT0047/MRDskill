#!/usr/bin/env python3
"""Experiment 6: longitudinal MRD monitoring and source-aware interpretation.

This script builds source data and Figure 6 for a longitudinal MRD experiment.
It uses locally staged public supplements where patient/table-level data are
available and keeps literature-level cohorts separate from parsed cohorts.
"""

from __future__ import annotations

import json
import math
import re
import textwrap
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch, PathPatch, Rectangle, Wedge
from matplotlib.path import Path as MplPath
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test
from lifelines.utils import concordance_index

from mrdskill_project_paths import artifact_root


ARTIFACTS = artifact_root()
DATA = ARTIFACTS / "data" / "benchmarks"
OUT = ARTIFACTS / "docs" / "benchmarks"
SUPP = ARTIFACTS / "inputs" / "supplementary_materials"

FIG_BASE = OUT / "figure_experiment6_longitudinal_monitoring"
SUMMARY_JSON = DATA / "experiment6_longitudinal_summary.json"
QA_JSON = OUT / "figure_experiment6_longitudinal_monitoring_scipilot_qa.json"

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

ORIGIN_COLORS = {
    "CHIP": PALETTE["orange"],
    "ctDNA": PALETTE["blue"],
    "CHIP+ctDNA": PALETTE["violet"],
    "none": PALETTE["light"],
}

STATE_COLORS = {
    "No cfDNA alteration": "#A7B0B8",
    "single-window or postoperative ctDNA": PALETTE["blue2"],
    "persistent or mixed persistent ctDNA": PALETTE["red"],
    "CHIP/WBC only": PALETTE["orange"],
    "single-window ctDNA": PALETTE["blue2"],
    "persistent ctDNA": PALETTE["blue"],
    "postoperative ctDNA": PALETTE["red"],
    "mixed-source single ctDNA": PALETTE["violet"],
    "mixed-source persistent ctDNA": PALETTE["teal"],
}

TIMEPOINT_ORDER = ["Baseline", "Preoperative", "Postoperative"]
TIMEPOINT_X = {"Baseline": 0.0, "Preoperative": 1.0, "Postoperative": 2.0}


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.size": 7.2,
            "axes.labelsize": 7.4,
            "xtick.labelsize": 6.7,
            "ytick.labelsize": 6.7,
            "legend.fontsize": 6.7,
            "axes.linewidth": 0.75,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.unicode_minus": False,
        }
    )


def panel_label(ax: plt.Axes, label: str, x: float = -0.08, y: float = 1.04) -> None:
    ax.text(x, y, label, transform=ax.transAxes, ha="left", va="bottom", fontsize=12, fontweight="bold")


def fmt_int(value: int | float | str | None) -> str:
    if value is None:
        return "NA"
    if isinstance(value, str):
        return value
    if not np.isfinite(float(value)):
        return "NA"
    return f"{int(round(float(value))):,}"


def fmt_float(value: float | int | None, digits: int = 2) -> str:
    if value is None or not np.isfinite(float(value)):
        return "NA"
    return f"{float(value):.{digits}f}"


def wrap_label(text: str, width: int = 18) -> str:
    return "\n".join(textwrap.wrap(text, width=width, break_long_words=False))


def find_one(pattern: str) -> Path:
    matches = list(SUPP.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No file matched {pattern} in {SUPP}")
    return matches[0]


def load_leal() -> dict[str, pd.DataFrame]:
    clinical = pd.read_excel(find_one("KWRAW-000854__supp04__*.xlsx"), header=2)
    cfdna = pd.read_excel(find_one("KWRAW-000854__supp09__*.xlsx"), header=1)
    wbc = pd.read_excel(find_one("KWRAW-000854__supp10__*.xlsx"), header=1)
    wbc_in_cfdna = pd.read_excel(find_one("KWRAW-000854__supp01__*.xlsx"), header=1)
    ctdna = pd.read_excel(find_one("KWRAW-000854__supp02__*.xlsx"), header=1)
    return {
        "clinical": clinical,
        "cfdna": cfdna,
        "wbc": wbc,
        "wbc_in_cfdna": wbc_in_cfdna,
        "ctdna": ctdna,
    }


def normalize_event_table(clinical: pd.DataFrame) -> pd.DataFrame:
    event_col = "Event-free Survival (0 - Event-free; 1 - Recurrence or Death)"
    efs_col = "Time to Event-free Survival (months)"
    out = clinical[["Patient", "Cancer Type", event_col, efs_col]].copy()
    out = out.rename(columns={event_col: "event", efs_col: "efs_months"})
    out["event_label"] = np.where(out["event"] == 1, "recurrence/death", "event-free")
    out.loc[out["event"].isna(), "event_label"] = "unknown"
    return out


def patient_state(origins_by_timepoint: dict[str, str]) -> str:
    ct_tps = {tp for tp, value in origins_by_timepoint.items() if "ctDNA" in value}
    chip_tps = {tp for tp, value in origins_by_timepoint.items() if "CHIP" in value}
    if not ct_tps and chip_tps:
        return "CHIP/WBC only"
    if len(ct_tps) >= 2 and chip_tps:
        return "mixed-source persistent ctDNA"
    if len(ct_tps) >= 2:
        return "persistent ctDNA"
    if "Postoperative" in ct_tps:
        return "postoperative ctDNA"
    if len(ct_tps) == 1 and chip_tps:
        return "mixed-source single ctDNA"
    if len(ct_tps) == 1:
        return "single-window ctDNA"
    return "no cfDNA alteration"


def build_leal_tables(leal: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    event = normalize_event_table(leal["clinical"])
    cfdna = leal["cfdna"].merge(event[["Patient", "event", "event_label", "efs_months"]], on="Patient", how="left")
    cfdna["vaf"] = pd.to_numeric(cfdna["Mutant Allele Fraction"], errors="coerce")
    cfdna["timepoint"] = cfdna["Timepoint Significance"].astype(str)
    cfdna["origin"] = cfdna["Origin"].astype(str)
    cfdna["gene"] = cfdna["Gene Symbol"].astype(str)
    cfdna.to_csv(DATA / "experiment6_leal_gastric_source_tracing.csv", index=False)

    rows: list[dict[str, Any]] = []
    for patient, group in cfdna.groupby("Patient"):
        origins_by_tp: dict[str, str] = {}
        for tp in TIMEPOINT_ORDER:
            origins = sorted(set(group.loc[group["timepoint"] == tp, "origin"].dropna()))
            if not origins:
                origins_by_tp[tp] = "none"
            else:
                origins_by_tp[tp] = "+".join(origins)
        state = patient_state(origins_by_tp)
        event_row = event[event["Patient"] == patient].head(1)
        event_value = float(event_row["event"].iloc[0]) if len(event_row) and pd.notna(event_row["event"].iloc[0]) else np.nan
        event_label = str(event_row["event_label"].iloc[0]) if len(event_row) else "unknown"
        efs_months = float(event_row["efs_months"].iloc[0]) if len(event_row) and pd.notna(event_row["efs_months"].iloc[0]) else np.nan
        rows.append(
            {
                "Patient": patient,
                "Baseline": origins_by_tp["Baseline"],
                "Preoperative": origins_by_tp["Preoperative"],
                "Postoperative": origins_by_tp["Postoperative"],
                "trajectory_state": state,
                "event": event_value,
                "event_label": event_label,
                "efs_months": efs_months,
                "ctdna_timepoints": sum("ctDNA" in origins_by_tp[tp] for tp in TIMEPOINT_ORDER),
                "chip_timepoints": sum("CHIP" in origins_by_tp[tp] for tp in TIMEPOINT_ORDER),
            }
        )
    trajectories = pd.DataFrame(rows)
    state_order = {
        "CHIP/WBC only": 0,
        "single-window ctDNA": 1,
        "mixed-source single ctDNA": 2,
        "persistent ctDNA": 3,
        "postoperative ctDNA": 4,
        "mixed-source persistent ctDNA": 5,
    }
    trajectories["state_order"] = trajectories["trajectory_state"].map(state_order).fillna(9)
    trajectories = trajectories.sort_values(["state_order", "event", "Patient"], ascending=[True, True, True]).reset_index(drop=True)
    trajectories.to_csv(DATA / "experiment6_leal_patient_trajectories.csv", index=False)

    source_counts = (
        cfdna.groupby(["timepoint", "origin", "event_label"])
        .size()
        .reset_index(name="alteration_count")
        .sort_values(["timepoint", "origin", "event_label"])
    )
    source_counts.to_csv(DATA / "experiment6_leal_origin_timepoint_event_counts.csv", index=False)

    patient_metrics = build_leal_reclassification_metrics(event, cfdna)
    patient_metrics.to_csv(DATA / "experiment6_leal_wbc_filter_metrics.csv", index=False)

    return {
        "cfdna": cfdna,
        "trajectories": trajectories,
        "source_counts": source_counts,
        "patient_metrics": patient_metrics,
    }


def classification_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    y_true = y_true.astype(bool)
    y_pred = y_pred.astype(bool)
    tp = int((y_true & y_pred).sum())
    fp = int((~y_true & y_pred).sum())
    fn = int((y_true & ~y_pred).sum())
    tn = int((~y_true & ~y_pred).sum())
    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "tn": tn,
        "sensitivity": tp / (tp + fn) if tp + fn else np.nan,
        "specificity": tn / (tn + fp) if tn + fp else np.nan,
        "ppv": tp / (tp + fp) if tp + fp else np.nan,
        "npv": tn / (tn + fn) if tn + fn else np.nan,
    }


def build_leal_reclassification_metrics(event: pd.DataFrame, cfdna: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    evaluable = event[event["event"].isin([0, 1])][["Patient", "event"]].copy()
    for tp in TIMEPOINT_ORDER:
        sub = cfdna[cfdna["timepoint"] == tp]
        any_pos = sub.groupby("Patient").size().rename("any_count")
        ct_pos = sub[sub["origin"] == "ctDNA"].groupby("Patient").size().rename("ctdna_count")
        chip_pos = sub[sub["origin"] == "CHIP"].groupby("Patient").size().rename("chip_count")
        table = evaluable.merge(any_pos, on="Patient", how="left").merge(ct_pos, on="Patient", how="left").merge(chip_pos, on="Patient", how="left")
        table = table.fillna({"any_count": 0, "ctdna_count": 0, "chip_count": 0})
        y = table["event"].to_numpy(dtype=int) == 1
        for mode, column in [("Unfiltered cfDNA positive", "any_count"), ("WBC-filtered ctDNA", "ctdna_count"), ("CHIP/WBC-only signal", "chip_count")]:
            metrics = classification_metrics(y, table[column].to_numpy(float) > 0)
            rows.append({"timepoint": tp, "mode": mode, "n_patients": len(table), **metrics})
    return pd.DataFrame(rows)


def collapsed_survival_group(state: str) -> str:
    if state == "no cfDNA alteration":
        return "No cfDNA alteration"
    if state == "CHIP/WBC only":
        return "CHIP/WBC only"
    if state in {"persistent ctDNA", "mixed-source persistent ctDNA"}:
        return "persistent or mixed persistent ctDNA"
    return "single-window or postoperative ctDNA"


def build_leal_survival_tables(leal: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    event = normalize_event_table(leal["clinical"])
    trajectories = leal["trajectories"].copy()
    survival = event[event["event"].isin([0, 1]) & event["efs_months"].notna()].copy()
    survival = survival.merge(
        trajectories[["Patient", "trajectory_state", "ctdna_timepoints", "chip_timepoints"]],
        on="Patient",
        how="left",
    )
    survival["trajectory_state"] = survival["trajectory_state"].fillna("no cfDNA alteration")
    survival["ctdna_timepoints"] = survival["ctdna_timepoints"].fillna(0).astype(int)
    survival["chip_timepoints"] = survival["chip_timepoints"].fillna(0).astype(int)
    survival["survival_group"] = survival["trajectory_state"].map(collapsed_survival_group)
    survival["any_cfdna_positive"] = (survival["trajectory_state"] != "no cfDNA alteration").astype(int)
    survival["any_ctdna_trajectory"] = (survival["ctdna_timepoints"] > 0).astype(int)
    survival["persistent_ctdna_trajectory"] = survival["survival_group"].eq("persistent or mixed persistent ctDNA").astype(int)
    survival["wbc_chip_only"] = survival["survival_group"].eq("CHIP/WBC only").astype(int)
    risk_score = {
        "CHIP/WBC only": 0,
        "single-window or postoperative ctDNA": 1,
        "persistent or mixed persistent ctDNA": 2,
    }
    survival["source_trajectory_score"] = survival["survival_group"].map(risk_score)
    survival.to_csv(DATA / "experiment6_leal_survival_groups.csv", index=False)

    group_order = [
        "CHIP/WBC only",
        "single-window or postoperative ctDNA",
        "persistent or mixed persistent ctDNA",
    ]
    km_rows: list[dict[str, Any]] = []
    group_summary_rows: list[dict[str, Any]] = []
    kmf = KaplanMeierFitter()
    km_input = survival[survival["source_trajectory_score"].notna()].copy()
    for group in group_order:
        sub = km_input[km_input["survival_group"] == group]
        if sub.empty:
            continue
        kmf.fit(sub["efs_months"], sub["event"], label=group)
        median_efs = float(kmf.median_survival_time_)
        timeline = kmf.survival_function_.reset_index()
        timeline.columns = ["month", "efs_probability"]
        ci = kmf.confidence_interval_survival_function_.reset_index()
        ci.columns = ["month", "ci_low", "ci_high"]
        merged = timeline.merge(ci, on="month", how="left")
        merged["group"] = group
        merged["n"] = len(sub)
        merged["events"] = int(sub["event"].sum())
        km_rows.extend(merged.to_dict("records"))
        group_summary_rows.append(
            {
                "group": group,
                "n": int(len(sub)),
                "events": int(sub["event"].sum()),
                "median_efs_months": median_efs,
                "median_efs_text": "NR" if np.isinf(median_efs) else f"{median_efs:.1f}",
            }
        )
    km_curve = pd.DataFrame(km_rows)
    km_curve.to_csv(DATA / "experiment6_leal_km_curve.csv", index=False)
    group_summary = pd.DataFrame(group_summary_rows)
    group_summary.to_csv(DATA / "experiment6_leal_survival_group_summary.csv", index=False)

    logrank_rows: list[dict[str, Any]] = []
    multi = multivariate_logrank_test(km_input["efs_months"], km_input["survival_group"], km_input["event"])
    logrank_rows.append(
        {
            "contrast": "overall survival_group",
            "test_statistic": float(multi.test_statistic),
            "p_value": float(multi.p_value),
        }
    )
    group_a = survival["survival_group"].eq("persistent or mixed persistent ctDNA")
    group_b = survival["survival_group"].eq("CHIP/WBC only")
    lr = logrank_test(
        survival.loc[group_a, "efs_months"],
        survival.loc[group_b, "efs_months"],
        event_observed_A=survival.loc[group_a, "event"],
        event_observed_B=survival.loc[group_b, "event"],
    )
    logrank_rows.append(
        {
            "contrast": "persistent/mixed persistent ctDNA vs CHIP/WBC only",
            "test_statistic": float(lr.test_statistic),
            "p_value": float(lr.p_value),
        }
    )
    logrank = pd.DataFrame(logrank_rows)
    logrank.to_csv(DATA / "experiment6_leal_logrank.csv", index=False)

    cfdna = leal["cfdna"]
    for tp in ["Preoperative", "Postoperative"]:
        sub = cfdna[cfdna["timepoint"] == tp]
        any_pos = set(sub["Patient"].dropna())
        ct_pos = set(sub.loc[sub["origin"] == "ctDNA", "Patient"].dropna())
        survival[f"any_cfdna_{tp.lower()}"] = survival["Patient"].isin(any_pos).astype(int)
        survival[f"wbc_filtered_ctdna_{tp.lower()}"] = survival["Patient"].isin(ct_pos).astype(int)

    cox_specs = [
        ("Preoperative any cfDNA", ["any_cfdna_preoperative"], survival),
        ("Preoperative WBC-filtered ctDNA", ["wbc_filtered_ctdna_preoperative"], survival),
        ("Postoperative any cfDNA", ["any_cfdna_postoperative"], survival),
        ("Postoperative WBC-filtered ctDNA", ["wbc_filtered_ctdna_postoperative"], survival),
        ("Source-aware trajectory score", ["source_trajectory_score"], km_input),
        ("Persistent ctDNA trajectory", ["persistent_ctdna_trajectory"], survival),
    ]
    cox_rows: list[dict[str, Any]] = []
    for label, cols, source_df in cox_specs:
        model_df = source_df[["efs_months", "event"] + cols].dropna().copy()
        cph = CoxPHFitter()
        cph.fit(model_df, duration_col="efs_months", event_col="event")
        score = model_df[cols].sum(axis=1)
        c_index = concordance_index(model_df["efs_months"], -score, model_df["event"])
        row = cph.summary.iloc[0]
        cox_rows.append(
            {
                "model": label,
                "term": cols[0],
                "hr": float(row["exp(coef)"]),
                "ci_low": float(row["exp(coef) lower 95%"]),
                "ci_high": float(row["exp(coef) upper 95%"]),
                "p_value": float(row["p"]),
                "c_index": float(c_index),
                "n": int(len(model_df)),
                "events": int(model_df["event"].sum()),
            }
        )
    cox = pd.DataFrame(cox_rows)
    cox.to_csv(DATA / "experiment6_leal_cox_models.csv", index=False)
    survival.to_csv(DATA / "experiment6_leal_survival_groups.csv", index=False)

    increment = cox[["model", "c_index", "hr", "ci_low", "ci_high", "p_value", "n", "events"]].copy()
    increment.to_csv(DATA / "experiment6_prognostic_increment.csv", index=False)

    return {
        "survival_groups": survival,
        "survival_group_summary": group_summary,
        "km_curve": km_curve,
        "logrank": logrank,
        "cox": cox,
        "increment": increment,
    }


def build_powles_tables() -> dict[str, pd.DataFrame]:
    transitions = pd.DataFrame(
        [
            {"arm": "Atezolizumab", "transition": "Neg>Neg", "median_dfs_months": np.nan, "median_os_months": np.nan, "dfs_text": "NR", "os_text": "NR"},
            {"arm": "Atezolizumab", "transition": "Pos>Neg clearance", "median_dfs_months": np.nan, "median_os_months": np.nan, "dfs_text": "NR", "os_text": "NR"},
            {"arm": "Atezolizumab", "transition": "Pos>Pos non-clearance", "median_dfs_months": 5.7, "median_os_months": 22.1, "dfs_text": "5.7", "os_text": "22.1"},
            {"arm": "Atezolizumab", "transition": "Neg>Pos", "median_dfs_months": 5.7, "median_os_months": 20.6, "dfs_text": "5.7", "os_text": "20.6"},
            {"arm": "Observation", "transition": "Neg>Neg", "median_dfs_months": np.nan, "median_os_months": np.nan, "dfs_text": "NR", "os_text": "NR"},
            {"arm": "Observation", "transition": "Pos>Neg clearance", "median_dfs_months": 8.8, "median_os_months": 8.8, "dfs_text": "8.8", "os_text": "8.8"},
            {"arm": "Observation", "transition": "Pos>Pos non-clearance", "median_dfs_months": 4.4, "median_os_months": 16.3, "dfs_text": "4.4", "os_text": "16.3"},
            {"arm": "Observation", "transition": "Neg>Pos", "median_dfs_months": 8.3, "median_os_months": 26.8, "dfs_text": "8.3", "os_text": "26.8"},
        ]
    )
    transitions.to_csv(DATA / "experiment6_powles_ctdna_transitions.csv", index=False)

    effects = pd.DataFrame(
        [
            {"dataset": "Powles IMvigor010", "contrast": "C1D1 ctDNA+ vs ctDNA-, Atezolizumab", "endpoint": "DFS", "hr": 3.39, "ci_low": 2.43, "ci_high": 4.47, "n": 300, "model": "multivariable"},
            {"dataset": "Powles IMvigor010", "contrast": "C1D1 ctDNA+ vs ctDNA-, Observation", "endpoint": "DFS", "hr": 6.19, "ci_low": 4.29, "ci_high": 8.91, "n": 281, "model": "multivariable"},
            {"dataset": "Powles IMvigor010", "contrast": "C3D1 ctDNA+ vs ctDNA-, Atezolizumab", "endpoint": "DFS", "hr": 5.37, "ci_low": 3.73, "ci_high": 7.74, "n": 263, "model": "multivariable"},
            {"dataset": "Powles IMvigor010", "contrast": "C3D1 ctDNA+ vs ctDNA-, Observation", "endpoint": "DFS", "hr": 8.36, "ci_low": 5.35, "ci_high": 13.07, "n": 222, "model": "multivariable"},
            {"dataset": "Powles IMvigor010", "contrast": "Clearance vs no clearance, Atezolizumab", "endpoint": "DFS", "hr": 0.32, "ci_low": 0.14, "ci_high": 0.71, "n": 99, "model": "multivariable"},
            {"dataset": "Powles IMvigor010", "contrast": "Clearance vs no clearance, Observation", "endpoint": "DFS", "hr": 0.42, "ci_low": 0.10, "ci_high": 1.79, "n": 79, "model": "multivariable"},
            {"dataset": "Powles IMvigor010", "contrast": "Atezolizumab vs observation, ctDNA+", "endpoint": "DFS", "hr": 0.56, "ci_low": 0.41, "ci_high": 0.77, "n": 214, "model": "multivariable"},
            {"dataset": "Powles IMvigor010", "contrast": "Atezolizumab vs observation, ctDNA-", "endpoint": "DFS", "hr": 1.07, "ci_low": 0.75, "ci_high": 1.52, "n": 367, "model": "multivariable"},
        ]
    )
    effects.to_csv(DATA / "experiment6_effect_sizes.csv", index=False)
    return {"transitions": transitions, "effects": effects}


def build_dataset_summary(leal_tables: dict[str, pd.DataFrame], powles_tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    leal_cfdna = leal_tables["cfdna"]
    leal_traj = leal_tables["trajectories"]
    source_labels = int(leal_cfdna["origin"].notna().sum())
    powles_transitions = len(powles_tables["transitions"])
    summary = pd.DataFrame(
        [
            {
                "dataset": "Leal gastric cancer",
                "cancer_type": "Gastric",
                "analysis_role": "parsed alteration-level source tracing",
                "patients": int(leal_traj["Patient"].nunique()),
                "samples_or_timepoints": int(leal_cfdna[["Patient", "timepoint"]].drop_duplicates().shape[0]),
                "variant_or_summary_units": int(len(leal_cfdna)),
                "source_labels": source_labels,
                "matched_wbc": "yes",
                "serial_status": "yes",
                "doi_or_source": "10.1038/s41467-020-14310-3",
            },
            {
                "dataset": "Powles IMvigor010",
                "cancer_type": "Urothelial",
                "analysis_role": "parsed serial transition summaries",
                "patients": 581,
                "samples_or_timepoints": 2,
                "variant_or_summary_units": powles_transitions,
                "source_labels": 0,
                "matched_wbc": "not in public summary",
                "serial_status": "yes",
                "doi_or_source": "10.1038/s41586-021-03642-9",
            },
        ]
    )
    summary.to_csv(DATA / "experiment6_dataset_summary.csv", index=False)
    return summary


def write_summary(leal: dict[str, pd.DataFrame], powles: dict[str, pd.DataFrame], datasets: pd.DataFrame) -> None:
    cfdna = leal["cfdna"]
    trajectories = leal["trajectories"]
    metrics = leal["patient_metrics"]
    pre_any = metrics[(metrics["timepoint"] == "Preoperative") & (metrics["mode"] == "Unfiltered cfDNA positive")].iloc[0]
    pre_ct = metrics[(metrics["timepoint"] == "Preoperative") & (metrics["mode"] == "WBC-filtered ctDNA")].iloc[0]
    post_any = metrics[(metrics["timepoint"] == "Postoperative") & (metrics["mode"] == "Unfiltered cfDNA positive")].iloc[0]
    post_ct = metrics[(metrics["timepoint"] == "Postoperative") & (metrics["mode"] == "WBC-filtered ctDNA")].iloc[0]
    state_counts = trajectories["trajectory_state"].value_counts().to_dict()
    effect = powles["effects"]
    survival = leal["survival"]
    cox = survival["cox"]
    logrank = survival["logrank"]
    pre_any_row = cox[cox["model"] == "Preoperative any cfDNA"].iloc[0]
    pre_ctdna_row = cox[cox["model"] == "Preoperative WBC-filtered ctDNA"].iloc[0]
    source_score_row = cox[cox["model"] == "Source-aware trajectory score"].iloc[0]
    summary = {
        "interpretation": "Longitudinal MRD monitoring requires source-aware trajectory states rather than isolated positive/negative calls.",
        "quantitative_datasets": int(len(datasets)),
        "parsed_datasets": ["Leal gastric cancer", "Powles IMvigor010"],
        "leal": {
            "cfDNA_alterations": int(len(cfdna)),
            "cfDNA_patients": int(cfdna["Patient"].nunique()),
            "chip_alterations": int((cfdna["origin"] == "CHIP").sum()),
            "ctdna_alterations": int((cfdna["origin"] == "ctDNA").sum()),
            "patient_state_counts": state_counts,
            "preoperative_unfiltered_ppv": float(pre_any["ppv"]),
            "preoperative_wbc_filtered_ppv": float(pre_ct["ppv"]),
            "postoperative_unfiltered_specificity": float(post_any["specificity"]),
            "postoperative_wbc_filtered_specificity": float(post_ct["specificity"]),
            "survival_evaluable_patients": int(len(survival["survival_groups"])),
            "cfDNA_positive_survival_patients": int(survival["survival_groups"]["source_trajectory_score"].notna().sum()),
            "survival_events": int(survival["survival_groups"]["event"].sum()),
            "cfDNA_positive_survival_events": int(survival["survival_groups"].loc[survival["survival_groups"]["source_trajectory_score"].notna(), "event"].sum()),
            "overall_logrank_p": float(logrank.loc[logrank["contrast"] == "overall survival_group", "p_value"].iloc[0]),
            "persistent_vs_chip_logrank_p": float(logrank.loc[logrank["contrast"] == "persistent/mixed persistent ctDNA vs CHIP/WBC only", "p_value"].iloc[0]),
            "source_trajectory_hr": float(source_score_row["hr"]),
            "source_trajectory_ci_low": float(source_score_row["ci_low"]),
            "source_trajectory_ci_high": float(source_score_row["ci_high"]),
            "source_trajectory_p": float(source_score_row["p_value"]),
            "preoperative_any_cfdna_c_index": float(pre_any_row["c_index"]),
            "preoperative_wbc_filtered_ctdna_c_index": float(pre_ctdna_row["c_index"]),
            "source_trajectory_c_index": float(source_score_row["c_index"]),
        },
        "powles": {
            "transition_rows": int(len(powles["transitions"])),
            "effect_rows": int(len(effect)),
            "c3d1_observation_ctdna_positive_dfs_hr": float(effect.loc[effect["contrast"] == "C3D1 ctDNA+ vs ctDNA-, Observation", "hr"].iloc[0]),
            "atezolizumab_clearance_dfs_hr": float(effect.loc[effect["contrast"] == "Clearance vs no clearance, Atezolizumab", "hr"].iloc[0]),
        },
        "outputs": {
            "dataset_summary": str(DATA / "experiment6_dataset_summary.csv"),
            "leal_source_tracing": str(DATA / "experiment6_leal_gastric_source_tracing.csv"),
            "leal_trajectories": str(DATA / "experiment6_leal_patient_trajectories.csv"),
            "leal_survival_groups": str(DATA / "experiment6_leal_survival_groups.csv"),
            "leal_survival_group_summary": str(DATA / "experiment6_leal_survival_group_summary.csv"),
            "leal_km_curve": str(DATA / "experiment6_leal_km_curve.csv"),
            "leal_logrank": str(DATA / "experiment6_leal_logrank.csv"),
            "leal_cox_models": str(DATA / "experiment6_leal_cox_models.csv"),
            "prognostic_increment": str(DATA / "experiment6_prognostic_increment.csv"),
            "powles_transitions": str(DATA / "experiment6_powles_ctdna_transitions.csv"),
            "effect_sizes": str(DATA / "experiment6_effect_sizes.csv"),
            "figure_png": str(FIG_BASE.with_suffix(".png")),
            "figure_pdf": str(FIG_BASE.with_suffix(".pdf")),
            "figure_svg": str(FIG_BASE.with_suffix(".svg")),
        },
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    QA_JSON.write_text(
        json.dumps(
            {
                "dpi": 600,
                "has_panel_titles": False,
                "has_overall_title": False,
                "formats": ["png", "pdf", "svg", "grayscale_png"],
                "palette_source": "Figure 2-compatible MRDskill palette",
                "visual_strategy": "six-panel longitudinal monitoring figure with analyzable-dataset funnel, patient barcode, Kaplan-Meier EFS, Cox/prognostic increment, serial-transition summary and source-correction phase plot",
            },
            indent=2,
        ),
        encoding="utf-8",
    )


def draw_swimmer(ax: plt.Axes, trajectories: pd.DataFrame) -> None:
    plot_df = trajectories.copy().reset_index(drop=True)
    y_positions = np.arange(len(plot_df))[::-1]
    for yi, (_, row) in zip(y_positions, plot_df.iterrows()):
        ax.plot([0, 2], [yi, yi], color=PALETTE["grid"], lw=0.8, zorder=0)
        for tp in TIMEPOINT_ORDER:
            val = row[tp]
            ax.scatter(
                [TIMEPOINT_X[tp]],
                [yi],
                s=24 if val == "none" else 42,
                color=ORIGIN_COLORS.get(val, PALETTE["mid"]),
                edgecolor="white",
                linewidth=0.55,
                zorder=2,
            )
        event_color = PALETTE["red"] if row["event"] == 1 else PALETTE["green"]
        ax.scatter([2.22], [yi], s=22, marker="s", color=event_color, edgecolor="white", linewidth=0.45, zorder=3)
    ax.set_xlim(-0.22, 2.45)
    ax.set_ylim(-1, len(plot_df))
    ax.set_xticks([0, 1, 2, 2.22])
    ax.set_xticklabels(["Base", "Preop", "Postop", "Event"], rotation=0)
    ax.set_yticks([])
    ax.set_xlabel("Gastric cohort timepoint")
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=ORIGIN_COLORS["ctDNA"], markeredgecolor="white", markersize=5.5, label="ctDNA"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=ORIGIN_COLORS["CHIP"], markeredgecolor="white", markersize=5.5, label="CHIP/WBC"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=ORIGIN_COLORS["CHIP+ctDNA"], markeredgecolor="white", markersize=5.5, label="mixed"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=PALETTE["red"], markeredgecolor="white", markersize=5.0, label="event"),
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(-0.02, 1.02), ncol=2, frameon=False, handletextpad=0.25, columnspacing=0.8)


def ribbon(ax: plt.Axes, x0: float, x1: float, y0a: float, y0b: float, y1a: float, y1b: float, color: str, alpha: float = 0.33) -> None:
    verts = [
        (x0, y0a),
        ((x0 + x1) / 2, y0a),
        ((x0 + x1) / 2, y1a),
        (x1, y1a),
        (x1, y1b),
        ((x0 + x1) / 2, y1b),
        ((x0 + x1) / 2, y0b),
        (x0, y0b),
        (x0, y0a),
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
    ax.add_patch(PathPatch(MplPath(verts, codes), facecolor=color, edgecolor="none", alpha=alpha))


def draw_alluvial(ax: plt.Axes, source_counts: pd.DataFrame) -> None:
    ax.axis("off")
    flow_rows: list[dict[str, Any]] = []
    for _, row in source_counts.iterrows():
        event = "event" if row["event_label"] == "recurrence/death" else "event-free"
        flow_rows.append({"timepoint": row["timepoint"], "origin": row["origin"], "event": event, "count": int(row["alteration_count"])})
    flow = pd.DataFrame(flow_rows)
    stage_cols = ["timepoint", "origin", "event"]
    orders = {
        "timepoint": TIMEPOINT_ORDER,
        "origin": ["CHIP", "ctDNA"],
        "event": ["event-free", "event"],
    }
    colors = {"CHIP": PALETTE["orange"], "ctDNA": PALETTE["blue"], "event-free": PALETTE["green"], "event": PALETTE["red"]}
    total = flow["count"].sum()
    x_positions = [0.10, 0.50, 0.90]
    spans: list[dict[str, tuple[float, float]]] = []
    for col in stage_cols:
        counts = flow.groupby(col)["count"].sum().to_dict()
        cursor = 0.07
        stage_span: dict[str, tuple[float, float]] = {}
        for key in orders[col]:
            count = counts.get(key, 0)
            if count <= 0:
                continue
            height = count / total * 0.84
            stage_span[key] = (cursor, cursor + height)
            cursor += height + 0.018
        spans.append(stage_span)

    for pair_cols, i0, i1, color_col in [(["timepoint", "origin"], 0, 1, "origin"), (["origin", "event"], 1, 2, "event")]:
        pair = flow.groupby(pair_cols)["count"].sum().reset_index()
        out_used: defaultdict[str, float] = defaultdict(float)
        in_used: defaultdict[str, float] = defaultdict(float)
        for _, row in pair.iterrows():
            left_key, right_key = str(row[pair_cols[0]]), str(row[pair_cols[1]])
            if left_key not in spans[i0] or right_key not in spans[i1]:
                continue
            height = float(row["count"]) / total * 0.84
            y0 = spans[i0][left_key][0] + out_used[left_key]
            y1 = spans[i1][right_key][0] + in_used[right_key]
            out_used[left_key] += height
            in_used[right_key] += height
            ribbon(ax, x_positions[i0] + 0.025, x_positions[i1] - 0.025, y0, y0 + height, y1, y1 + height, colors.get(str(row[color_col]), PALETTE["mid"]), 0.28)

    for x, col, stage_span in zip(x_positions, stage_cols, spans):
        for key, (y0, y1) in stage_span.items():
            color = PALETTE["light"] if col == "timepoint" else colors.get(key, PALETTE["mid"])
            ax.add_patch(Rectangle((x - 0.026, y0), 0.052, y1 - y0, facecolor=color, edgecolor="white", lw=0.7, alpha=0.95))
            label = key.replace("Preoperative", "Preop").replace("Postoperative", "Postop")
            ax.text(x + 0.040, (y0 + y1) / 2, label, ha="left", va="center", fontsize=6.6, color=PALETTE["dark"])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


def draw_raincloud(ax: plt.Axes, cfdna: pd.DataFrame) -> None:
    rng = np.random.default_rng(31)
    rows = []
    for tp in TIMEPOINT_ORDER:
        for origin in ["CHIP", "ctDNA"]:
            values = cfdna[(cfdna["timepoint"] == tp) & (cfdna["origin"] == origin)]["vaf"].dropna().to_numpy(float)
            rows.append((tp, origin, values))
    positions = np.arange(len(rows))
    for i, (_, origin, values) in enumerate(rows):
        log_values = np.log10(values + 1e-6)
        if len(log_values) > 1:
            parts = ax.violinplot(log_values, positions=[i], widths=0.72, showmeans=False, showmedians=False, showextrema=False)
            for body in parts["bodies"]:
                body.set_facecolor(ORIGIN_COLORS[origin])
                body.set_edgecolor("none")
                body.set_alpha(0.24)
        jitter = rng.normal(0, 0.06, len(log_values))
        ax.scatter(np.full_like(log_values, i, dtype=float) + jitter, log_values, s=12, color=ORIGIN_COLORS[origin], edgecolor="white", linewidth=0.3, alpha=0.74)
        if len(log_values):
            ax.plot([i - 0.28, i + 0.28], [np.median(log_values), np.median(log_values)], color=PALETTE["dark"], lw=0.8)
    labels = []
    for tp in TIMEPOINT_ORDER:
        labels.extend([tp[:4] + "\nCHIP", tp[:4] + "\nctDNA"])
    ax.set_xticks(positions)
    ax.set_xticklabels(labels)
    ax.set_ylabel("log10 VAF")
    ax.grid(axis="y", color=PALETTE["grid"], lw=0.7)


def draw_powles_transition(ax: plt.Axes, transitions: pd.DataFrame) -> None:
    order = ["Neg>Neg", "Pos>Neg clearance", "Pos>Pos non-clearance", "Neg>Pos"]
    y_base = np.arange(len(order))[::-1]
    offsets = {"Atezolizumab": 0.12, "Observation": -0.12}
    colors = {"Atezolizumab": PALETTE["blue"], "Observation": PALETTE["orange"]}
    cap_value = 31.0
    for arm, sub in transitions.groupby("arm"):
        sub = sub.set_index("transition").loc[order].reset_index()
        ys = y_base + offsets[arm]
        values = sub["median_dfs_months"].fillna(cap_value).to_numpy(float)
        for y, value, text in zip(ys, values, sub["dfs_text"]):
            ax.plot([0, value], [y, y], color=colors[arm], lw=1.8, alpha=0.35, solid_capstyle="round")
            marker = ">" if text == "NR" else "o"
            ax.scatter([value], [y], s=42, marker=marker, color=colors[arm], edgecolor="white", linewidth=0.7, zorder=3)
            ax.text(value + 0.8, y, text, va="center", ha="left", fontsize=6.2, color=colors[arm])
    ax.set_yticks(y_base)
    ax.set_yticklabels(["Neg>Neg", "Clearance", "Non-clearance", "Neg>Pos"])
    ax.set_xlim(0, 35)
    ax.set_xlabel("Median DFS, months")
    ax.grid(axis="x", color=PALETTE["grid"], lw=0.7)
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color=colors["Atezolizumab"], lw=1.5, label="Atezolizumab"),
            Line2D([0], [0], marker="o", color=colors["Observation"], lw=1.5, label="Observation"),
        ],
        loc="lower right",
        frameon=False,
    )


def draw_forest(ax: plt.Axes, effects: pd.DataFrame) -> None:
    keep = effects.copy()
    keep["label"] = keep["contrast"].str.replace(", ", "\n", regex=False)
    keep = keep.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(keep))
    colors = [PALETTE["red"] if hr > 1 else PALETTE["green"] for hr in keep["hr"]]
    ax.axvline(1, color=PALETTE["mid"], lw=0.8, ls="--")
    for yi, (_, row), color in zip(y, keep.iterrows(), colors):
        ax.plot([row["ci_low"], row["ci_high"]], [yi, yi], color=color, lw=1.8, alpha=0.8)
        ax.scatter([row["hr"]], [yi], s=36, color=color, edgecolor="white", linewidth=0.7, zorder=3)
        ax.text(row["ci_high"] * 1.12, yi, f"{row['hr']:.2g}", ha="left", va="center", fontsize=6.0, color=PALETTE["dark"])
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels([wrap_label(v, 25) for v in keep["label"]])
    ax.set_xlabel("Hazard ratio")
    ax.set_xlim(0.08, 65)
    ax.grid(axis="x", color=PALETTE["grid"], lw=0.7, which="both")


def draw_reclassification(ax: plt.Axes, metrics: pd.DataFrame) -> None:
    sub = metrics[metrics["mode"].isin(["Unfiltered cfDNA positive", "WBC-filtered ctDNA"])].copy()
    metrics_to_show = [("specificity", "Specificity"), ("ppv", "PPV")]
    timepoints = ["Preoperative", "Postoperative"]
    x_positions = np.arange(len(timepoints) * len(metrics_to_show))
    idx = 0
    for tp in timepoints:
        for metric, label in metrics_to_show:
            vals = []
            for mode in ["Unfiltered cfDNA positive", "WBC-filtered ctDNA"]:
                vals.append(float(sub[(sub["timepoint"] == tp) & (sub["mode"] == mode)][metric].iloc[0]))
            ax.plot([idx - 0.12, idx + 0.12], vals, color=PALETTE["mid"], lw=1.0, alpha=0.6)
            ax.scatter([idx - 0.12], [vals[0]], s=34, color=PALETTE["orange"], edgecolor="white", linewidth=0.7, zorder=3)
            ax.scatter([idx + 0.12], [vals[1]], s=34, color=PALETTE["blue"], edgecolor="white", linewidth=0.7, zorder=3)
            ax.text(idx + 0.20, vals[1], f"{vals[1]:.2f}", ha="left", va="center", fontsize=6.1, color=PALETTE["blue"])
            idx += 1
    ax.set_xticks(x_positions)
    ax.set_xticklabels(["Preop\nSpec.", "Preop\nPPV", "Postop\nSpec.", "Postop\nPPV"])
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("Patient-level metric")
    ax.grid(axis="y", color=PALETTE["grid"], lw=0.7)
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color="none", markerfacecolor=PALETTE["orange"], markeredgecolor="white", label="Any cfDNA"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor=PALETTE["blue"], markeredgecolor="white", label="ctDNA after WBC filter"),
        ],
        loc="lower right",
        frameon=False,
        handletextpad=0.25,
    )


def draw_decision_examples(ax: plt.Axes, trajectories: pd.DataFrame) -> None:
    ax.set_axis_off()
    examples = [
        ("CHIP/WBC only", "CHIP/WBC only", "source-ambiguous", "WBC audit; suppress tumor overcall"),
        ("single-window ctDNA", "single-window ctDNA", "low temporal support", "repeat / serial confirmation"),
        ("persistent ctDNA", "persistent ctDNA", "recurrence pattern", "escalate recurrence review"),
        ("mixed-source persistent ctDNA", "mixed persistent ctDNA", "tumor + CHIP background", "retain ctDNA; annotate CHIP"),
    ]
    counts = trajectories["trajectory_state"].value_counts().to_dict()
    y_positions = [0.80, 0.58, 0.36, 0.14]
    for (state, display_state, interpretation, action), y in zip(examples, y_positions):
        color = STATE_COLORS[state]
        count = int(counts.get(state, 0))
        ax.add_patch(Rectangle((0.03, y - 0.075), 0.94, 0.145, facecolor="#FAFBFC", edgecolor="#E2E6ED", lw=0.65))
        ax.add_patch(Rectangle((0.03, y - 0.075), 0.012, 0.145, facecolor=color, edgecolor="none"))
        ax.add_patch(Circle((0.08, y), 0.034 + 0.006 * math.sqrt(max(count, 1)), facecolor=color, edgecolor="white", lw=0.8))
        ax.text(0.15, y + 0.040, display_state, ha="left", va="center", fontsize=6.8, fontweight="bold", color=PALETTE["dark"], linespacing=0.95)
        ax.text(0.15, y + 0.006, f"{interpretation}; n={count}", ha="left", va="center", fontsize=5.9, color=PALETTE["mid"])
        ax.text(0.15, y - 0.038, wrap_label(action, 43), ha="left", va="center", fontsize=6.1, color=PALETTE["dark"], linespacing=0.95)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


def draw_longitudinal_evidence_gate(ax: plt.Axes, datasets: pd.DataFrame, leal: dict[str, pd.DataFrame], powles: dict[str, pd.DataFrame]) -> None:
    """Eligibility panel showing which longitudinal resources are actually analyzed."""
    ax.set_axis_off()
    leal_cfdna = leal["cfdna"]
    rows = [
        ("Longitudinal MRD screen", "public studies assessed for serial ctDNA, source labels and outcomes", "8 screened", PALETTE["mid"]),
        ("Leal gastric cancer", f"{leal_cfdna['Patient'].nunique()} cfDNA-positive patients; {len(leal_cfdna)} source-labelled alterations", "source-aware", PALETTE["blue"]),
        ("IMvigor010 urothelial", f"581 patients; {len(powles['transitions'])} baseline-to-C3D1 transition groups", "serial DFS", PALETTE["teal"]),
    ]
    y_positions = [0.79, 0.50, 0.22]
    x0, box_w, box_h = 0.08, 0.82, 0.18
    for i, ((name, detail, tag, color), y) in enumerate(zip(rows, y_positions)):
        ax.add_patch(
            Rectangle((x0, y - box_h / 2), box_w, box_h, facecolor="#FAFBFC", edgecolor="#DCE1E8", lw=0.9)
        )
        ax.add_patch(Rectangle((x0, y - box_h / 2), 0.035, box_h, facecolor=color, edgecolor="none", alpha=0.95))
        ax.text(x0 + 0.055, y + 0.038, name, ha="left", va="center", fontsize=7.4, fontweight="bold", color=PALETTE["dark"])
        ax.text(x0 + 0.055, y - 0.012, wrap_label(detail, 58), ha="left", va="center", fontsize=6.4, color=PALETTE["mid"], linespacing=0.95)
        ax.text(x0 + box_w - 0.025, y + 0.038, tag, ha="right", va="center", fontsize=6.1, color=color, fontweight="bold")
        if i < len(rows) - 1:
            ax.add_patch(
                FancyArrowPatch(
                    (x0 + box_w / 2, y - box_h / 2 - 0.018),
                    (x0 + box_w / 2, y_positions[i + 1] + box_h / 2 + 0.018),
                    arrowstyle="-|>",
                    mutation_scale=10,
                    lw=0.8,
                    color=PALETTE["mid"],
                    alpha=0.75,
                )
            )
    ax.text(0.50, 0.045, "Quantitative analysis uses the two resources with analyzable serial outcome data", ha="center", va="center", fontsize=6.4, color=PALETTE["mid"])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


def draw_temporal_barcode(ax: plt.Axes, trajectories: pd.DataFrame) -> None:
    """Time-anchored patient barcode emphasizing longitudinal source status."""
    plot_df = trajectories.copy().reset_index(drop=True)
    state_order = [
        "CHIP/WBC only",
        "single-window ctDNA",
        "mixed-source single ctDNA",
        "persistent ctDNA",
        "postoperative ctDNA",
        "mixed-source persistent ctDNA",
    ]
    plot_df["trajectory_state"] = pd.Categorical(plot_df["trajectory_state"], state_order, ordered=True)
    plot_df = plot_df.sort_values(["trajectory_state", "event", "Patient"]).reset_index(drop=True)
    y_positions = np.arange(len(plot_df))[::-1]
    for yi, (_, row) in zip(y_positions, plot_df.iterrows()):
        state_color = STATE_COLORS.get(str(row["trajectory_state"]), PALETTE["mid"])
        ax.plot([-0.08, 2.20], [yi, yi], color="#E9EDF3", lw=0.75, zorder=0)
        ax.add_patch(Rectangle((-0.23, yi - 0.33), 0.055, 0.66, facecolor=state_color, edgecolor="none", alpha=0.95))
        last_x = None
        last_val = None
        for tp in TIMEPOINT_ORDER:
            val = row[tp]
            x = TIMEPOINT_X[tp]
            if val != "none":
                if last_x is not None and last_val != "none":
                    ax.plot([last_x, x], [yi, yi], color=state_color, lw=1.5, alpha=0.55, zorder=1)
                ax.scatter([x], [yi], s=52, color=ORIGIN_COLORS.get(val, PALETTE["mid"]), edgecolor="white", linewidth=0.6, zorder=3)
            else:
                ax.scatter([x], [yi], s=22, color=PALETTE["light"], edgecolor="white", linewidth=0.4, zorder=2)
            last_x = x
            last_val = val
        event_color = PALETTE["red"] if row["event"] == 1 else PALETTE["green"]
        ax.scatter([2.26], [yi], s=28, marker="s", color=event_color, edgecolor="white", linewidth=0.55, zorder=3)
    ax.set_xlim(-0.30, 2.45)
    ax.set_ylim(-1, len(plot_df))
    ax.set_xticks([0, 1, 2, 2.26])
    ax.set_xticklabels(["Baseline", "Preop", "Postop", "Event"])
    ax.set_yticks([])
    ax.set_xlabel("Gastric cohort monitoring timeline")
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=ORIGIN_COLORS["ctDNA"], markeredgecolor="white", markersize=6.0, label="ctDNA"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=ORIGIN_COLORS["CHIP"], markeredgecolor="white", markersize=6.0, label="CHIP/WBC"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=ORIGIN_COLORS["CHIP+ctDNA"], markeredgecolor="white", markersize=6.0, label="mixed"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=PALETTE["red"], markeredgecolor="white", markersize=5.2, label="event"),
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(-0.01, 1.02), ncol=4, frameon=False, handletextpad=0.25, columnspacing=0.9)
    ax.text(-0.28, len(plot_df) - 0.3, "state", fontsize=6.3, color=PALETTE["mid"], rotation=90, ha="left", va="top")


def draw_state_event_ring(ax: plt.Axes, trajectories: pd.DataFrame) -> None:
    """Nested donut: trajectory-state composition and event fraction per state."""
    ax.set_aspect("equal")
    ax.set_axis_off()
    states = [
        "CHIP/WBC only",
        "single-window ctDNA",
        "mixed-source single ctDNA",
        "persistent ctDNA",
        "postoperative ctDNA",
        "mixed-source persistent ctDNA",
    ]
    counts = trajectories["trajectory_state"].value_counts().to_dict()
    event_counts = trajectories[trajectories["event"] == 1]["trajectory_state"].value_counts().to_dict()
    total = sum(counts.get(s, 0) for s in states)
    start = 95.0
    center = (0.36, 0.50)
    for state in states:
        count = counts.get(state, 0)
        if count <= 0:
            continue
        theta = 360 * count / total
        color = STATE_COLORS[state]
        ax.add_patch(Wedge(center, 0.36, start, start + theta, width=0.095, facecolor=color, edgecolor="white", lw=1.0))
        event = event_counts.get(state, 0)
        event_theta = theta * event / count if count else 0
        if event_theta > 0:
            ax.add_patch(Wedge(center, 0.235, start, start + event_theta, width=0.085, facecolor=PALETTE["red"], edgecolor="white", lw=0.7, alpha=0.85))
        if theta > 22:
            mid = math.radians(start + theta / 2)
            x = center[0] + 0.43 * math.cos(mid)
            y = center[1] + 0.43 * math.sin(mid)
            ax.text(x, y, str(count), ha="center", va="center", fontsize=7, color=PALETTE["dark"])
        start += theta
    ax.add_patch(Circle(center, 0.125, facecolor="white", edgecolor="#DCE1E8", lw=0.8))
    ax.text(center[0], center[1] + 0.030, "43", ha="center", va="center", fontsize=12, fontweight="bold", color=PALETTE["dark"])
    ax.text(center[0], center[1] - 0.045, "cfDNA+\npatients", ha="center", va="center", fontsize=6.2, color=PALETTE["mid"], linespacing=0.95)

    legend_rows = [
        ("CHIP/WBC only", "9/11 event-free"),
        ("single-window ctDNA", "3/9 event"),
        ("persistent ctDNA", "5/7 event"),
        ("mixed-source persistent ctDNA", "6/9 event"),
    ]
    y = 0.76
    for state, note in legend_rows:
        ax.add_patch(Circle((0.78, y), 0.020, facecolor=STATE_COLORS[state], edgecolor="white", lw=0.5))
        ax.text(0.82, y + 0.015, state.replace("mixed-source ", "mixed "), ha="left", va="center", fontsize=6.7, color=PALETTE["dark"], fontweight="bold")
        ax.text(0.82, y - 0.017, note, ha="left", va="center", fontsize=6.0, color=PALETTE["mid"])
        y -= 0.16
    ax.add_patch(Rectangle((0.70, 0.12), 0.025, 0.025, facecolor=PALETTE["red"], edgecolor="white", lw=0.5))
    ax.text(0.74, 0.132, "inner red arc: recurrence/death", ha="left", va="center", fontsize=6.1, color=PALETTE["mid"])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


def draw_km_survival(ax: plt.Axes, survival: dict[str, pd.DataFrame]) -> None:
    """Kaplan-Meier EFS curves by source-aware trajectory group."""
    km = survival["km_curve"].copy()
    group_order = [
        "CHIP/WBC only",
        "single-window or postoperative ctDNA",
        "persistent or mixed persistent ctDNA",
    ]
    color_map = {
        "CHIP/WBC only": STATE_COLORS["CHIP/WBC only"],
        "single-window or postoperative ctDNA": STATE_COLORS["single-window or postoperative ctDNA"],
        "persistent or mixed persistent ctDNA": STATE_COLORS["persistent or mixed persistent ctDNA"],
    }
    for group in group_order:
        sub = km[km["group"] == group]
        if sub.empty:
            continue
        ax.step(sub["month"], sub["efs_probability"], where="post", color=color_map[group], lw=1.8, label=group)
        ax.fill_between(
            sub["month"].to_numpy(float),
            sub["ci_low"].to_numpy(float),
            sub["ci_high"].to_numpy(float),
            step="post",
            color=color_map[group],
            alpha=0.10,
            linewidth=0,
        )
    logrank = survival["logrank"]
    p = float(logrank.loc[logrank["contrast"] == "overall survival_group", "p_value"].iloc[0])
    ax.text(0.05, 0.08, f"log-rank p={p:.3g}", transform=ax.transAxes, ha="left", va="bottom", fontsize=6.6, color=PALETTE["dark"])
    ax.set_xlabel("Event-free survival (months)")
    ax.set_ylabel("Event-free probability")
    ax.set_ylim(0.0, 1.03)
    ax.set_xlim(left=0)
    ax.grid(color=PALETTE["grid"], lw=0.7)
    handles, labels = ax.get_legend_handles_labels()
    new_labels = [
        "WBC/CHIP only",
        "single-window ctDNA",
        "persistent ctDNA",
    ]
    ax.legend(handles, new_labels, loc="upper right", frameon=False, handlelength=1.5, handletextpad=0.45)


def draw_prognostic_increment(ax: plt.Axes, survival: dict[str, pd.DataFrame]) -> None:
    """C-index and HR summary for increasingly source-aware longitudinal states."""
    inc = survival["increment"].copy()
    display = [
        ("Preoperative any cfDNA", "Any cfDNA"),
        ("Preoperative WBC-filtered ctDNA", "WBC-filtered"),
        ("Source-aware trajectory score", "Source-aware"),
    ]
    inc = inc.set_index("model").loc[[d[0] for d in display]].reset_index()
    y = np.arange(len(inc))[::-1]
    colors = [PALETTE["orange"], PALETTE["blue2"], PALETTE["red"]]
    ax.axvline(0.5, color=PALETTE["mid"], lw=0.8, ls="--")
    ax.plot([0.48, 0.67], [y.min() - 0.55, y.min() - 0.55], color=PALETTE["grid"], lw=1.0, zorder=0)
    for yi, (_, row), color in zip(y, inc.iterrows(), colors):
        ax.plot([0.5, row["c_index"]], [yi, yi], color=color, lw=3.0, alpha=0.35, solid_capstyle="round")
        ax.scatter([row["c_index"]], [yi], s=78, color=color, edgecolor="white", linewidth=0.8, zorder=3)
        label_x = min(float(row["c_index"]) + 0.012, 0.655)
        ax.text(label_x, yi + 0.17, f"C={row['c_index']:.2f}", ha="left", va="center", fontsize=6.4, color=PALETTE["dark"])
        ax.text(
            label_x,
            yi - 0.16,
            f"HR {row['hr']:.2g} ({row['ci_low']:.2g}-{row['ci_high']:.2g})",
            ha="left",
            va="center",
            fontsize=5.9,
            color=PALETTE["mid"],
        )
    ax.set_yticks(y)
    ax.set_yticklabels([d[1] for d in display])
    ax.set_xlim(0.46, 0.72)
    ax.set_ylim(-0.8, len(inc) - 0.35)
    ax.set_xlabel("Concordance index")
    ax.grid(axis="x", color=PALETTE["grid"], lw=0.7)
    ax.text(0.02, 0.04, "dots: C-index; text: Cox HR (95% CI)", transform=ax.transAxes, ha="left", va="bottom", fontsize=5.9, color=PALETTE["mid"])


def draw_transition_metro(ax: plt.Axes, transitions: pd.DataFrame) -> None:
    """Metro-map view of serial ctDNA transition classes and DFS."""
    ax.set_axis_off()
    y_map = {
        "Neg>Neg": 0.78,
        "Pos>Neg clearance": 0.58,
        "Pos>Pos non-clearance": 0.38,
        "Neg>Pos": 0.18,
    }
    status_color = {"Neg": PALETTE["green"], "Pos": PALETTE["red"]}
    arm_color = {"Atezolizumab": PALETTE["blue"], "Observation": PALETTE["orange"]}
    max_month = 31.0
    for transition, y in y_map.items():
        left_status, right_status = transition.split(" ")[0].split(">") if " " in transition else transition.split(">")
        ax.plot([0.10, 0.42], [y, y], color=PALETTE["grid"], lw=6, solid_capstyle="round")
        ax.scatter([0.10], [y], s=80, color=status_color[left_status], edgecolor="white", linewidth=0.8, zorder=3)
        ax.scatter([0.42], [y], s=80, color=status_color[right_status], edgecolor="white", linewidth=0.8, zorder=3)
        ax.add_patch(FancyArrowPatch((0.15, y), (0.37, y), arrowstyle="-|>", mutation_scale=8, lw=0.8, color=PALETTE["mid"], zorder=4))
        label = transition.replace("Pos>Neg clearance", "clearance").replace("Pos>Pos non-clearance", "non-clearance")
        ax.text(0.02, y, label, ha="left", va="center", fontsize=7.0, color=PALETTE["dark"])
        sub = transitions[transitions["transition"] == transition]
        for _, row in sub.iterrows():
            arm = row["arm"]
            value = row["median_dfs_months"]
            capped = 31.0 if pd.isna(value) else float(value)
            x0 = 0.56
            x1 = x0 + 0.33 * min(capped, max_month) / max_month
            offset = 0.035 if arm == "Atezolizumab" else -0.035
            ax.plot([x0, x1], [y + offset, y + offset], color=arm_color[arm], lw=2.5, alpha=0.75, solid_capstyle="round")
            marker = ">" if row["dfs_text"] == "NR" else "o"
            ax.scatter([x1], [y + offset], s=35, marker=marker, color=arm_color[arm], edgecolor="white", linewidth=0.6, zorder=3)
            ax.text(x1 + 0.020, y + offset, str(row["dfs_text"]), ha="left", va="center", fontsize=6.2, color=arm_color[arm])
    ax.text(0.10, 0.91, "C1D1", ha="center", va="center", fontsize=6.5, color=PALETTE["mid"])
    ax.text(0.42, 0.91, "C3D1", ha="center", va="center", fontsize=6.5, color=PALETTE["mid"])
    ax.text(0.56, 0.91, "Median DFS", ha="left", va="center", fontsize=6.5, color=PALETTE["mid"])
    ax.legend(
        handles=[
            Line2D([0], [0], color=arm_color["Atezolizumab"], lw=2, label="Atezolizumab"),
            Line2D([0], [0], color=arm_color["Observation"], lw=2, label="Observation"),
        ],
        loc="lower right",
        frameon=False,
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


def draw_source_correction_phase(ax: plt.Axes, metrics: pd.DataFrame) -> None:
    """Phase plot of WBC-source correction across monitoring windows."""
    sub = metrics[metrics["mode"].isin(["Unfiltered cfDNA positive", "WBC-filtered ctDNA"])].copy()
    timepoints = ["Baseline", "Preoperative", "Postoperative"]
    color_map = {"Baseline": PALETTE["mid"], "Preoperative": PALETTE["orange"], "Postoperative": PALETTE["blue"]}
    for tp in timepoints:
        points = []
        for mode in ["Unfiltered cfDNA positive", "WBC-filtered ctDNA"]:
            row = sub[(sub["timepoint"] == tp) & (sub["mode"] == mode)].iloc[0]
            points.append((float(row["specificity"]), float(row["ppv"])))
        (x0, y0), (x1, y1) = points
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=10, lw=1.2, color=color_map[tp], alpha=0.88))
        ax.scatter([x0], [y0], s=42, color=PALETTE["light"], edgecolor=color_map[tp], linewidth=1.0, zorder=3)
        ax.scatter([x1], [y1], s=54, color=color_map[tp], edgecolor="white", linewidth=0.8, zorder=4)
        ax.text(x1 + 0.025, y1, tp.replace("Preoperative", "Preop").replace("Postoperative", "Postop"), ha="left", va="center", fontsize=6.6, color=color_map[tp])
    ax.axline((0, 0), slope=1, color=PALETTE["grid"], lw=0.8, ls="--")
    ax.set_xlim(0.0, 1.02)
    ax.set_ylim(0.0, 0.88)
    ax.set_xlabel("Specificity")
    ax.set_ylabel("Positive predictive value")
    ax.grid(color=PALETTE["grid"], lw=0.7)
    ax.text(0.04, 0.82, "open circle: any cfDNA\nfilled circle: WBC-filtered ctDNA", ha="left", va="top", fontsize=6.3, color=PALETTE["mid"])


def make_figure(leal: dict[str, pd.DataFrame], powles: dict[str, pd.DataFrame], datasets: pd.DataFrame) -> None:
    setup_style()
    OUT.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(12.0, 10.8), dpi=300)
    gs = fig.add_gridspec(3, 4, height_ratios=[1.0, 1.18, 1.0], hspace=0.58, wspace=0.50)
    ax_a = fig.add_subplot(gs[0, :2])
    ax_b = fig.add_subplot(gs[0:2, 2:])
    ax_c = fig.add_subplot(gs[1, :2])
    ax_d = fig.add_subplot(gs[2, 0])
    ax_e = fig.add_subplot(gs[2, 1:3])
    ax_f = fig.add_subplot(gs[2, 3])

    draw_longitudinal_evidence_gate(ax_a, datasets, leal, powles)
    draw_temporal_barcode(ax_b, leal["trajectories"])
    draw_km_survival(ax_c, leal["survival"])
    draw_prognostic_increment(ax_d, leal["survival"])
    draw_transition_metro(ax_e, powles["transitions"])
    draw_source_correction_phase(ax_f, leal["patient_metrics"])

    for ax, label in zip([ax_a, ax_b, ax_c, ax_d, ax_e, ax_f], list("abcdef")):
        panel_label(ax, label)

    for ax in [ax_b, ax_c, ax_d, ax_f]:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.savefig(FIG_BASE.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(FIG_BASE.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(FIG_BASE.with_suffix(".svg"), bbox_inches="tight")
    try:
        from PIL import Image

        Image.open(FIG_BASE.with_suffix(".png")).convert("L").save(OUT / "figure_experiment6_longitudinal_monitoring_grayscale.png")
    except Exception:
        pass
    plt.close(fig)


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    leal_raw = load_leal()
    leal_tables = build_leal_tables(leal_raw)
    leal_tables["clinical"] = leal_raw["clinical"]
    leal_tables["survival"] = build_leal_survival_tables(leal_tables)
    powles_tables = build_powles_tables()
    dataset_summary = build_dataset_summary(leal_tables, powles_tables)
    make_figure(leal_tables, powles_tables, dataset_summary)
    write_summary(leal_tables, powles_tables, dataset_summary)
    print(
        json.dumps(
            {
                "figure_png": str(FIG_BASE.with_suffix(".png")),
                "figure_pdf": str(FIG_BASE.with_suffix(".pdf")),
                "figure_svg": str(FIG_BASE.with_suffix(".svg")),
                "summary_json": str(SUMMARY_JSON),
                "source_tables": [
                    str(DATA / "experiment6_dataset_summary.csv"),
                    str(DATA / "experiment6_leal_gastric_source_tracing.csv"),
                    str(DATA / "experiment6_leal_patient_trajectories.csv"),
                    str(DATA / "experiment6_leal_origin_timepoint_event_counts.csv"),
                    str(DATA / "experiment6_leal_wbc_filter_metrics.csv"),
                    str(DATA / "experiment6_leal_survival_groups.csv"),
                    str(DATA / "experiment6_leal_survival_group_summary.csv"),
                    str(DATA / "experiment6_leal_km_curve.csv"),
                    str(DATA / "experiment6_leal_logrank.csv"),
                    str(DATA / "experiment6_leal_cox_models.csv"),
                    str(DATA / "experiment6_prognostic_increment.csv"),
                    str(DATA / "experiment6_powles_ctdna_transitions.csv"),
                    str(DATA / "experiment6_effect_sizes.csv"),
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
