#!/usr/bin/env python3
"""Collect MRDskill QC outputs into one folder per sample.

The collector copies files instead of moving them so legacy data/qc and docs/qc
paths remain valid for existing reports and notebooks.
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path
from typing import Any

from mrdskill_project_paths import artifact_path, qc_case_dir, safe_path_component


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Mirror QC outputs into artifacts/cases/<sample>/qc.")
    parser.add_argument("--case-id", required=True, help="QC case id/prefix, for example sample_Mutect2_filtered.")
    parser.add_argument("--sample-id", default="", help="Biological sample folder name. Defaults to tumor_sample in summary JSON or case id.")
    parser.add_argument("--summary-json", default="", help="Summary JSON. Defaults to data/qc/<case-id>.qc_summary.json.")
    parser.add_argument("--variant-csv", default="", help="Variant QC CSV. Defaults to data/qc/<case-id>.qc_variants.csv.")
    parser.add_argument("--qc-report-md", default="", help="Detailed QC report. Defaults to docs/qc/<case-id>_qc_report.md.")
    parser.add_argument("--case-report-md", default="", help="Case-level QC note. Defaults to docs/qc/<case-id>_case_qc_note.md.")
    parser.add_argument("--decision-html", default="", help="Decision HTML. Defaults to data/qc/<case-id>_decision_qc.html.")
    parser.add_argument("--decision-md", default="", help="Decision Markdown. Defaults to data/qc/<case-id>_decision_qc.md.")
    parser.add_argument("--shortlist-csv", default="", help="Decision shortlist CSV. Defaults to data/qc/<case-id>_decision_shortlist.csv.")
    parser.add_argument("--manifest-json", default="", help="Optional output manifest path.")
    return parser.parse_args()


def read_summary(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def default_path(value: str, fallback: str) -> Path:
    return artifact_path(value) if value else artifact_path(fallback)


def copy_if_exists(source: Path, target: Path, copied: list[dict[str, str]], missing: list[str]) -> None:
    if not source.exists():
        missing.append(str(source))
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    copied.append({"source": str(source), "target": str(target), "bytes": str(target.stat().st_size)})


def collect_case_qc_outputs(
    case_id: str,
    sample_id: str = "",
    summary_json: str | Path = "",
    variant_csv: str | Path = "",
    qc_report_md: str | Path = "",
    case_report_md: str | Path = "",
    decision_html: str | Path = "",
    decision_md: str | Path = "",
    shortlist_csv: str | Path = "",
    manifest_json: str | Path = "",
    include_decision: bool = True,
) -> dict[str, Any]:
    """Copy available QC artifacts into artifacts/cases/<sample>/qc/<case-id>."""
    safe_case_id = safe_path_component(case_id)
    summary_path = default_path(str(summary_json), f"data/qc/{safe_case_id}.qc_summary.json")
    summary = read_summary(summary_path)
    resolved_sample_id = sample_id or summary.get("tumor_sample") or safe_case_id
    case_dir = qc_case_dir(resolved_sample_id, safe_case_id)
    qc_dir = case_dir / "qc"
    detail_dir = qc_dir / safe_case_id
    detail_dir.mkdir(parents=True, exist_ok=True)

    files = {
        "summary_json": (summary_path, detail_dir / f"{safe_case_id}.qc_summary.json"),
        "variant_csv": (
            default_path(str(variant_csv), f"data/qc/{safe_case_id}.qc_variants.csv"),
            detail_dir / f"{safe_case_id}.qc_variants.csv",
        ),
        "qc_report_md": (
            default_path(str(qc_report_md), f"docs/qc/{safe_case_id}_qc_report.md"),
            detail_dir / f"{safe_case_id}_qc_report.md",
        ),
        "case_report_md": (
            default_path(str(case_report_md), f"docs/qc/{safe_case_id}_case_qc_note.md"),
            detail_dir / f"{safe_case_id}_case_qc_note.md",
        ),
    }
    if include_decision:
        files.update({
            "decision_html": (
            default_path(str(decision_html), f"data/qc/{safe_case_id}_decision_qc.html"),
            detail_dir / f"{safe_case_id}_decision_qc.html",
            ),
            "decision_md": (
            default_path(str(decision_md), f"data/qc/{safe_case_id}_decision_qc.md"),
            detail_dir / f"{safe_case_id}_decision_qc.md",
            ),
            "shortlist_csv": (
            default_path(str(shortlist_csv), f"data/qc/{safe_case_id}_decision_shortlist.csv"),
            detail_dir / f"{safe_case_id}_decision_shortlist.csv",
            ),
        })

    copied: list[dict[str, str]] = []
    missing: list[str] = []
    for _, (source, target) in files.items():
        copy_if_exists(source, target, copied, missing)

    manifest = {
        "case_id": safe_case_id,
        "sample_id": str(resolved_sample_id),
        "case_dir": str(case_dir),
        "qc_dir": str(detail_dir),
        "input_vcf": summary.get("input_vcf", ""),
        "summary": {
            "record_count": summary.get("record_count"),
            "case_qc_class": summary.get("case_assessment", {}).get("case_qc_class"),
            "mrd_signal_class": summary.get("case_assessment", {}).get("mrd_signal_class"),
            "tier_counts": summary.get("tier_counts", {}),
            "kg_match_counts": summary.get("kg_match_counts", {}),
        },
        "copied": copied,
        "missing": missing,
    }

    manifest_path = artifact_path(manifest_json) if manifest_json else detail_dir / f"{safe_case_id}_qc_manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    manifest["manifest_json"] = str(manifest_path)
    return manifest


def main() -> None:
    args = parse_args()
    manifest = collect_case_qc_outputs(
        case_id=args.case_id,
        sample_id=args.sample_id,
        summary_json=args.summary_json,
        variant_csv=args.variant_csv,
        qc_report_md=args.qc_report_md,
        case_report_md=args.case_report_md,
        decision_html=args.decision_html,
        decision_md=args.decision_md,
        shortlist_csv=args.shortlist_csv,
        manifest_json=args.manifest_json,
    )

    print(json.dumps({"case_dir": manifest["case_dir"], "qc_dir": manifest["qc_dir"], "copied": len(manifest["copied"]), "missing": len(manifest["missing"]), "manifest_json": manifest["manifest_json"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
