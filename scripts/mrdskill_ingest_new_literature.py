#!/usr/bin/env python3
"""Ingest a new local literature/supplement folder into MRDskill evidence inputs.

This is the incremental-data entry point. It expects a prepared corpus folder
with the same broad layout as `kg_ready_no_pdf` and keeps every run in an
archive-first staging directory before merging de-duplicated rows into the full
evidence/import layers used by KG/RAG rebuilds.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Callable

from mrdskill_project_paths import artifact_path, artifact_root, safe_path_component, skill_script_path


EVIDENCE_MINER = skill_script_path("mrd-evidence-miner", "mrdskill_extract_full_corpus_evidence.py")
CONFIDENCE_RATER = skill_script_path("mrd-confidence-rater", "mrdskill_score_mrd_evidence.py")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Incrementally ingest new MRD literature and supplementary materials.")
    parser.add_argument("--input", required=True, help="Prepared literature/supplement corpus folder.")
    parser.add_argument("--run-id", default="", help="Stable run id. Defaults to current timestamp.")
    parser.add_argument("--max-files", type=int, default=0, help="0 means all files.")
    parser.add_argument("--max-snippets-per-document", type=int, default=6)
    parser.add_argument("--max-table-evidence-per-document", type=int, default=200)
    parser.add_argument("--max-rows-per-sheet", type=int, default=0, help="0 means all rows.")
    parser.add_argument("--excel-workers", type=int, default=4)
    parser.add_argument("--no-progress", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--full-raw-csv", default="data/evidence/mrd_evidence_raw_full.csv")
    parser.add_argument("--full-scored-csv", default="data/evidence/mrd_evidence_scored_full.csv")
    parser.add_argument("--full-scored-jsonl", default="data/evidence/mrd_evidence_scored_full.jsonl")
    parser.add_argument("--full-sources-csv", default="data/imported/full_corpus_sources.csv")
    parser.add_argument("--full-excel-chunks-csv", default="data/imported/full_corpus_excel_chunks.csv")
    parser.add_argument("--report-json", default="", help="Optional report path. Defaults inside the ingest run directory.")
    return parser.parse_args()


def clean(value: Any) -> str:
    return "" if value is None else str(value).strip()


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.exists():
        return [], []
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_jsonl(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False) + "\n")


def union_fields(*field_lists: list[str], rows: list[dict[str, str]] | None = None) -> list[str]:
    fields: list[str] = []
    seen: set[str] = set()
    for field_list in field_lists:
        for field in field_list:
            if field and field not in seen:
                fields.append(field)
                seen.add(field)
    for row in rows or []:
        for field in row:
            if field and field not in seen:
                fields.append(field)
                seen.add(field)
    return fields


def evidence_content_key(row: dict[str, str]) -> str:
    parts = [
        clean(row.get("paper_id")),
        clean(row.get("doi")),
        clean(row.get("pmid")),
        clean(row.get("source_document_id")),
        clean(row.get("source_file")),
        clean(row.get("source_sheet")),
        clean(row.get("source_row")),
        clean(row.get("source_column")),
        clean(row.get("source_section")),
        clean(row.get("gene_symbol_normalized") or row.get("gene")),
        clean(row.get("variant")),
        clean(row.get("genomic_coordinate")),
        clean(row.get("quoted_or_paraphrased_evidence"))[:900],
    ]
    return "\u241f".join(parts).lower()


def stable_incremental_evidence_id(row: dict[str, str]) -> str:
    digest = hashlib.sha1(evidence_content_key(row).encode("utf-8", errors="ignore")).hexdigest()[:16]
    return f"EV-INCR-{digest}"


def source_key(row: dict[str, str]) -> str:
    return "\u241f".join(
        [
            clean(row.get("source_document_id")),
            clean(row.get("paper_id")),
            clean(row.get("source_file")),
            clean(row.get("extension")),
            clean(row.get("size_bytes")),
        ]
    ).lower()


def chunk_key(row: dict[str, str]) -> str:
    return "\u241f".join(
        [
            clean(row.get("chunk_id")),
            clean(row.get("source_document_id")),
            clean(row.get("source_sheet")),
            clean(row.get("row_start")),
            clean(row.get("row_end")),
        ]
    ).lower()


def prepare_evidence_rows(rows: list[dict[str, str]], run_id: str, input_dir: Path, existing_ids: set[str]) -> list[dict[str, str]]:
    prepared: list[dict[str, str]] = []
    for row in rows:
        item = dict(row)
        evidence_id = stable_incremental_evidence_id(item)
        item["evidence_id"] = evidence_id
        item["ingest_run_id"] = run_id
        item["ingest_input_dir"] = str(input_dir)
        existing_ids.add(evidence_id)
        prepared.append(item)
    return prepared


def merge_csv_unique(
    target: Path,
    staged: Path,
    key_func: Callable[[dict[str, str]], str],
    transform_rows: Callable[[list[dict[str, str]], set[str]], list[dict[str, str]]] | None = None,
) -> dict[str, Any]:
    target_fields, target_rows = read_csv(target)
    staged_fields, staged_rows = read_csv(staged)
    existing_keys = {key_func(row) for row in target_rows}
    existing_ids = {clean(row.get("evidence_id")) for row in target_rows if clean(row.get("evidence_id"))}
    incoming_rows = transform_rows(staged_rows, existing_ids) if transform_rows else staged_rows
    added: list[dict[str, str]] = []
    skipped = 0
    for row in incoming_rows:
        key = key_func(row)
        if key in existing_keys:
            skipped += 1
            continue
        existing_keys.add(key)
        added.append(row)
    if added:
        all_rows = target_rows + added
        fields = union_fields(target_fields, staged_fields, rows=all_rows)
        write_csv(target, all_rows, fields)
    return {
        "target": str(target),
        "staged": str(staged),
        "existing_rows": len(target_rows),
        "staged_rows": len(staged_rows),
        "added_rows": len(added),
        "skipped_duplicate_rows": skipped,
        "final_rows": len(target_rows) + len(added),
    }


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


def main() -> int:
    args = parse_args()
    input_dir = Path(args.input).resolve()
    if not input_dir.exists() and not args.dry_run:
        raise FileNotFoundError(f"Input folder not found: {input_dir}")
    if not EVIDENCE_MINER.exists():
        raise FileNotFoundError(f"Evidence miner script not found: {EVIDENCE_MINER}")
    if not CONFIDENCE_RATER.exists():
        raise FileNotFoundError(f"Confidence rater script not found: {CONFIDENCE_RATER}")

    run_id = safe_path_component(args.run_id or datetime.now().strftime("%Y%m%d_%H%M%S"))
    run_dir = artifact_root() / "imports" / "literature_ingest" / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    staged_raw = run_dir / "mrd_evidence_raw_increment.csv"
    staged_raw_jsonl = run_dir / "mrd_evidence_raw_increment.jsonl"
    staged_sources = run_dir / "full_corpus_sources_increment.csv"
    staged_chunks = run_dir / "full_corpus_excel_chunks_increment.csv"
    staged_extract_report = run_dir / "extract_report.json"
    staged_scored = run_dir / "mrd_evidence_scored_increment.csv"
    staged_scored_jsonl = run_dir / "mrd_evidence_scored_increment.jsonl"

    extract_command = [
        sys.executable,
        str(EVIDENCE_MINER),
        "--corpus-dir",
        str(input_dir),
        "--out-csv",
        str(staged_raw),
        "--out-jsonl",
        str(staged_raw_jsonl),
        "--sources-csv",
        str(staged_sources),
        "--chunk-log-csv",
        str(staged_chunks),
        "--report-json",
        str(staged_extract_report),
        "--max-files",
        str(args.max_files),
        "--max-snippets-per-document",
        str(args.max_snippets_per_document),
        "--max-table-evidence-per-document",
        str(args.max_table_evidence_per_document),
        "--max-rows-per-sheet",
        str(args.max_rows_per_sheet),
        "--excel-workers",
        str(args.excel_workers),
    ]
    if args.no_progress:
        extract_command.append("--no-progress")

    score_command = [
        sys.executable,
        str(CONFIDENCE_RATER),
        "--input-csv",
        str(staged_raw),
        "--out-csv",
        str(staged_scored),
        "--out-jsonl",
        str(staged_scored_jsonl),
    ]

    report: dict[str, Any] = {
        "mode": "incremental_literature_ingest",
        "run_id": run_id,
        "input_dir": str(input_dir),
        "run_dir": str(run_dir),
        "started": datetime.now().isoformat(timespec="seconds"),
        "steps": [],
        "merge": {},
    }

    for command in [extract_command, score_command]:
        result = run_step(command, args.dry_run)
        report["steps"].append(result)
        print(result["output_tail"])
        if result["returncode"] != 0:
            report["status"] = "failed"
            break
    else:
        if args.dry_run:
            report["status"] = "completed_dry_run"
        else:
            full_raw = artifact_path(args.full_raw_csv)
            full_scored = artifact_path(args.full_scored_csv)
            full_scored_jsonl = artifact_path(args.full_scored_jsonl)
            full_sources = artifact_path(args.full_sources_csv)
            full_chunks = artifact_path(args.full_excel_chunks_csv)

            report["merge"]["raw_evidence"] = merge_csv_unique(
                full_raw,
                staged_raw,
                evidence_content_key,
                lambda rows, existing_ids: prepare_evidence_rows(rows, run_id, input_dir, existing_ids),
            )
            report["merge"]["scored_evidence"] = merge_csv_unique(
                full_scored,
                staged_scored,
                evidence_content_key,
                lambda rows, existing_ids: prepare_evidence_rows(rows, run_id, input_dir, existing_ids),
            )
            scored_fields, scored_rows = read_csv(full_scored)
            write_jsonl(full_scored_jsonl, scored_rows)
            report["merge"]["scored_jsonl"] = {"target": str(full_scored_jsonl), "rows": len(scored_rows)}
            report["merge"]["sources"] = merge_csv_unique(full_sources, staged_sources, source_key)
            report["merge"]["excel_chunks"] = merge_csv_unique(full_chunks, staged_chunks, chunk_key)
            report["status"] = "completed"

    report["completed"] = datetime.now().isoformat(timespec="seconds")
    report_path = artifact_path(args.report_json) if args.report_json else run_dir / "ingest_report.json"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"status": report["status"], "run_id": run_id, "run_dir": str(run_dir), "report_json": str(report_path)}, ensure_ascii=False, indent=2))
    return 0 if report["status"] in {"completed", "completed_dry_run"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
