#!/usr/bin/env python3
"""Freshly rebuild MRDskill KG and RAG artifacts from retained evidence inputs.

The rebuild is intentionally archive-first instead of destructive: previous
graph, RAG, cohort, and rule outputs are moved into a timestamped archive before
new artifacts are generated from the stable evidence/import layers.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

from mrdskill_project_paths import PROJECT_CODE_ROOT, artifact_path, artifact_root, script_path, skill_script_path


ARCHIVE_TARGETS = [
    "data/graph",
    "data/rag",
    "data/qc_rules",
    "data/cohort",
    "docs/rag",
    "docs/graph",
]

STATIC_DOCS = [
    "docs/rag/incremental_literature_ingest.md",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Archive and rebuild the latest MRDskill KG/RAG from evidence inputs.")
    parser.add_argument("--archive-label", default="", help="Optional label appended to the archive folder name.")
    parser.add_argument("--extensions", default=".xlsx,.xls,.csv,.tsv")
    parser.add_argument("--max-excel-sources", type=int, default=0, help="0 means no Excel source cap.")
    parser.add_argument("--max-candidate-sheets-total", type=int, default=0, help="0 means no candidate sheet cap.")
    parser.add_argument("--max-source-size-mb", type=float, default=0.0, help="0 means no size cap.")
    parser.add_argument("--max-sheets-per-source", type=int, default=0, help="0 means no sheet cap.")
    parser.add_argument("--max-rows-per-sheet", type=int, default=200000)
    parser.add_argument("--min-samples", type=int, default=20)
    parser.add_argument("--min-pair-support", type=int, default=3)
    parser.add_argument("--max-mutex-rules", type=int, default=2500)
    parser.add_argument("--skip-visualization", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--report-json", default="docs/rag/mrd_fresh_rebuild_report.json")
    return parser.parse_args()


def run_step(command: list[str], dry_run: bool = False) -> dict[str, Any]:
    started = datetime.now().isoformat(timespec="seconds")
    if dry_run:
        return {
            "command": command,
            "returncode": 0,
            "started": started,
            "completed": datetime.now().isoformat(timespec="seconds"),
            "output_tail": "dry-run",
        }
    proc = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return {
        "command": command,
        "returncode": proc.returncode,
        "started": started,
        "completed": datetime.now().isoformat(timespec="seconds"),
        "output_tail": proc.stdout[-8000:],
    }


def ensure_inputs() -> None:
    required = [
        artifact_path("data/evidence/mrd_evidence_scored_full.csv"),
        artifact_path("data/imported/full_corpus_sources.csv"),
        artifact_path("data/imported/full_corpus_excel_chunks.csv"),
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Required rebuild inputs are missing: " + "; ".join(missing))


def archive_existing(label: str, dry_run: bool = False) -> dict[str, Any]:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    suffix = f"_{label}" if label else ""
    archive_dir = artifact_root() / "archive" / f"kg_rag_rebuild_{timestamp}{suffix}"
    moved: list[dict[str, str]] = []
    for relative in ARCHIVE_TARGETS:
        source = artifact_path(relative)
        if not source.exists():
            continue
        destination = archive_dir / relative
        moved.append({"from": str(source), "to": str(destination)})
        if dry_run:
            continue
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(source), str(destination))
    return {"archive_dir": str(archive_dir), "moved": moved}


def build_commands(args: argparse.Namespace) -> list[list[str]]:
    python = sys.executable
    kg_builder = skill_script_path("mrd-kg-builder", "mrdskill_build_mrd_kg_csv.py")
    commands = [
        [
            python,
            str(kg_builder),
            "--input-csv",
            str(artifact_path("data/evidence/mrd_evidence_scored_full.csv")),
            "--nodes-csv",
            str(artifact_path("data/graph/nodes_full.csv")),
            "--relationships-csv",
            str(artifact_path("data/graph/relationships_full.csv")),
            "--review-tasks-csv",
            str(artifact_path("data/graph/review_tasks_full.csv")),
            "--evidence-summaries-jsonl",
            str(artifact_path("data/rag/evidence_summaries_full.jsonl")),
            "--sources-csv",
            str(artifact_path("data/imported/full_corpus_sources.csv")),
        ],
        [
            python,
            str(script_path("mrdskill_build_supplementary_cohort_rules.py")),
            "--extensions",
            args.extensions,
            "--max-excel-sources",
            str(args.max_excel_sources),
            "--max-candidate-sheets-total",
            str(args.max_candidate_sheets_total),
            "--max-source-size-mb",
            str(args.max_source_size_mb),
            "--max-sheets-per-source",
            str(args.max_sheets_per_source),
            "--max-rows-per-sheet",
            str(args.max_rows_per_sheet),
            "--min-samples",
            str(args.min_samples),
            "--min-pair-support",
            str(args.min_pair_support),
        ],
        [
            python,
            str(script_path("mrdskill_augment_kg_with_cohort_rules.py")),
            "--max-mutex-rules",
            str(args.max_mutex_rules),
        ],
        [python, str(script_path("mrdskill_build_adjudication_layer.py"))],
        [python, str(script_path("mrdskill_build_unified_rag_index.py"))],
        [python, str(script_path("mrdskill_build_latest_kg_index.py"))],
        [
            python,
            str(script_path("mrdskill_audit_kg_coverage_fast.py")),
            "--out-json",
            "docs/rag/mrd_latest_kg_coverage_audit.json",
            "--out-md",
            "docs/rag/mrd_latest_kg_coverage_audit.md",
        ],
    ]
    if not args.skip_visualization:
        commands.append([python, str(script_path("mrdskill_build_mrd_kg_visualization.py"))])
    return commands


def sync_static_docs(dry_run: bool = False) -> list[dict[str, str]]:
    copied: list[dict[str, str]] = []
    for relative in STATIC_DOCS:
        source = PROJECT_CODE_ROOT / relative
        destination = artifact_path(relative)
        if not source.exists():
            continue
        copied.append({"from": str(source), "to": str(destination)})
        if dry_run:
            continue
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
    return copied


def main() -> int:
    args = parse_args()
    ensure_inputs()
    report: dict[str, Any] = {
        "started": datetime.now().isoformat(timespec="seconds"),
        "mode": "fresh_rebuild",
        "dry_run": args.dry_run,
        "archive": archive_existing(args.archive_label, args.dry_run),
        "steps": [],
    }
    for command in build_commands(args):
        result = run_step(command, args.dry_run)
        report["steps"].append(result)
        print(result["output_tail"])
        if result["returncode"] != 0:
            report["status"] = "failed"
            break
    else:
        report["status"] = "completed"
    report["static_docs"] = sync_static_docs(args.dry_run)
    report["completed"] = datetime.now().isoformat(timespec="seconds")

    report_path = artifact_path(args.report_json)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"status": report["status"], "report_json": str(report_path)}, ensure_ascii=False, indent=2))
    return 0 if report["status"] == "completed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
