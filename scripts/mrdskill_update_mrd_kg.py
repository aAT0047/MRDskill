#!/usr/bin/env python3
"""Orchestrate MRD KG/RAG rebuilds and incremental literature ingests."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

from mrdskill_project_paths import artifact_path, script_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Update MRD KG, cohort rules, latest rule DB, visualization, and unified RAG index.")
    parser.add_argument("--mode", choices=["rebuild", "incremental"], default="rebuild")
    parser.add_argument("--input", default="", help="Prepared new literature/supplement folder for incremental mode.")
    parser.add_argument("--extensions", default=".xlsx,.xls,.csv,.tsv", help="Supplementary extensions to scan for cohort rules.")
    parser.add_argument("--min-samples", type=int, default=5)
    parser.add_argument("--min-pair-support", type=int, default=2)
    parser.add_argument("--max-mutex-rules", type=int, default=2500)
    parser.add_argument("--max-files", type=int, default=0, help="Incremental ingest cap; 0 means all files.")
    parser.add_argument("--max-rows-per-sheet", type=int, default=0, help="0 means all rows.")
    parser.add_argument("--skip-visualization", action="store_true")
    parser.add_argument("--rebuild-cohort-rules", action="store_true", help="Re-scan all supplementary sources for cohort rules. Default preserves existing rule CSVs.")
    parser.add_argument("--dry-run", action="store_true", help="Print and record the planned ingest/rebuild steps without changing artifacts.")
    parser.add_argument("--report-json", default="docs/rag/mrd_kg_update_report.json")
    return parser.parse_args()


def run_step(command: list[str]) -> dict[str, Any]:
    started = datetime.now().isoformat(timespec="seconds")
    completed = None
    proc = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    completed = datetime.now().isoformat(timespec="seconds")
    return {
        "command": command,
        "returncode": proc.returncode,
        "started": started,
        "completed": completed,
        "output_tail": proc.stdout[-4000:],
    }


def main() -> None:
    args = parse_args()
    if args.mode == "incremental" and not args.input:
        print("incremental mode requires --input with a prepared literature/supplement folder", file=sys.stderr)
        sys.exit(2)

    steps: list[list[str]] = []
    if args.mode == "incremental":
        steps.append([
            sys.executable,
            str(script_path("mrdskill_ingest_new_literature.py")),
            "--input",
            args.input,
            "--max-files",
            str(args.max_files),
            "--max-rows-per-sheet",
            str(args.max_rows_per_sheet),
        ])
        if args.dry_run:
            steps[-1].append("--dry-run")

    if args.rebuild_cohort_rules:
        rebuild_command = [
            sys.executable,
            str(script_path("mrdskill_rebuild_latest_kg_rag.py")),
            "--archive-label",
            args.mode,
            "--extensions",
            args.extensions,
            "--max-rows-per-sheet",
            str(args.max_rows_per_sheet),
            "--min-samples",
            str(args.min_samples),
            "--min-pair-support",
            str(args.min_pair_support),
            "--max-mutex-rules",
            str(args.max_mutex_rules),
        ]
        if args.skip_visualization:
            rebuild_command.append("--skip-visualization")
        if args.dry_run:
            rebuild_command.append("--dry-run")
        steps.append(rebuild_command)
    else:
        steps.extend(
            [
                [sys.executable, str(script_path("mrdskill_augment_kg_with_cohort_rules.py")), "--max-mutex-rules", str(args.max_mutex_rules)],
                [sys.executable, str(script_path("mrdskill_build_adjudication_layer.py"))],
                [sys.executable, str(script_path("mrdskill_build_latest_kg_index.py"))],
                [
                    sys.executable,
                    str(script_path("mrdskill_audit_kg_coverage_fast.py")),
                    "--out-json",
                    "docs/rag/mrd_latest_kg_coverage_audit.json",
                    "--out-md",
                    "docs/rag/mrd_latest_kg_coverage_audit.md",
                ],
            ]
        )
        if not args.skip_visualization:
            steps.append([sys.executable, str(script_path("mrdskill_build_mrd_kg_visualization.py"))])

    report: dict[str, Any] = {
        "mode": args.mode,
        "input": args.input,
        "dry_run": args.dry_run,
        "started": datetime.now().isoformat(timespec="seconds"),
        "steps": [],
    }
    for command in steps:
        result = run_step(command)
        report["steps"].append(result)
        print(result["output_tail"])
        if result["returncode"] != 0:
            report["status"] = "failed"
            break
    else:
        report["status"] = "completed"
    report["completed"] = datetime.now().isoformat(timespec="seconds")

    report_path = artifact_path(args.report_json)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"status": report["status"], "report_json": str(report_path)}, ensure_ascii=False, indent=2))
    if report["status"] != "completed":
        sys.exit(1)


if __name__ == "__main__":
    main()
