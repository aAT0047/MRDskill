#!/usr/bin/env python3
"""Inventory an MRD KG-ready corpus without moving source files."""

from __future__ import annotations

import argparse
import csv
import json
import os
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


def count_extensions(path: Path) -> dict[str, object]:
    counts: Counter[str] = Counter()
    if not path.exists():
        return {"exists": False, "extensions": {}, "scan_errors": []}

    scan_errors: list[dict[str, str]] = []

    def walk(directory: Path) -> None:
        try:
            with os.scandir(directory) as entries:
                for entry in entries:
                    try:
                        if entry.is_dir(follow_symlinks=False):
                            walk(Path(entry.path))
                        elif entry.is_file(follow_symlinks=False):
                            suffix = os.path.splitext(entry.name)[1].lower() or "<no_ext>"
                            counts[suffix] += 1
                    except OSError as exc:
                        scan_errors.append({"path": entry.path, "error": repr(exc)})
        except OSError as exc:
            scan_errors.append({"path": str(directory), "error": repr(exc)})

    walk(path)
    return {
        "exists": True,
        "extensions": dict(sorted(counts.items())),
        "total_files": sum(counts.values()),
        "scan_errors": scan_errors,
    }


def read_csv_summary(path: Path) -> dict[str, object]:
    if not path.exists():
        return {"exists": False}

    rows = 0
    columns: list[str] = []
    status_counts: Counter[str] = Counter()
    sample_rows: list[dict[str, str]] = []

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        columns = list(reader.fieldnames or [])
        for row in reader:
            rows += 1
            if "status" in row:
                status_counts[row.get("status", "") or "<blank>"] += 1
            if len(sample_rows) < 3:
                sample_rows.append({key: row.get(key, "") for key in columns[:12]})

    return {
        "exists": True,
        "rows": rows,
        "columns": columns,
        "status_counts": dict(status_counts),
        "sample_rows": sample_rows,
    }


def summarize_candidate_table(path: Path) -> dict[str, object]:
    summary = read_csv_summary(path)
    if not summary.get("exists"):
        return summary

    year_counts: Counter[str] = Counter()
    journal_counts: Counter[str] = Counter()
    topic_counts: Counter[str] = Counter()
    research_flag_counts: Counter[str] = Counter()

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            year_counts[row.get("year", "") or "not_reported"] += 1
            journal_counts[row.get("journal", "") or "not_reported"] += 1
            research_flag_counts[row.get("research_article_flag", "") or "not_reported"] += 1
            for tag in (row.get("likely_topic_tags", "") or "").split(";"):
                tag = tag.strip()
                if tag:
                    topic_counts[tag] += 1

    summary.update(
        {
            "year_counts": dict(year_counts.most_common()),
            "top_journals": dict(journal_counts.most_common(25)),
            "top_topic_tags": dict(topic_counts.most_common(50)),
            "research_article_flag_counts": dict(research_flag_counts),
        }
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create a JSON inventory for an MRD kg_ready_no_pdf corpus."
    )
    parser.add_argument("--corpus-dir", required=True, help="Path to kg_ready_no_pdf.")
    parser.add_argument(
        "--out-json",
        default="data/imported/kg_ready_inventory.json",
        help="Output JSON path inside the active project.",
    )
    args = parser.parse_args()

    corpus_dir = Path(args.corpus_dir).resolve()
    out_json = Path(args.out_json).resolve()

    if not corpus_dir.exists():
        raise SystemExit(f"Corpus directory not found: {corpus_dir}")

    manifests = corpus_dir / "manifests"
    inventory = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "corpus_dir": str(corpus_dir),
        "layout": {
            "downloaded_pdfs_converted_from_pdf": count_extensions(
                corpus_dir / "downloaded_pdfs" / "converted_from_pdf"
            ),
            "downloaded_pdfs_publisher_fulltext_html_xml": count_extensions(
                corpus_dir / "downloaded_pdfs" / "publisher_fulltext_html_xml"
            ),
            "supplementary_materials_converted_from_pdf": count_extensions(
                corpus_dir / "supplementary_materials" / "converted_from_pdf"
            ),
            "supplementary_materials_native_non_pdf": count_extensions(
                corpus_dir / "supplementary_materials" / "native_non_pdf"
            ),
        },
        "manifests": {
            "downloaded_pdfs_converted_manifest": read_csv_summary(
                manifests / "downloaded_pdfs_converted_manifest.csv"
            ),
            "supplementary_materials_converted_manifest": read_csv_summary(
                manifests / "supplementary_materials_converted_manifest.csv"
            ),
            "keyword_research_candidate_table": summarize_candidate_table(
                manifests / "keyword_research_candidate_table.csv"
            ),
        },
    }

    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(
        json.dumps(inventory, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(f"Wrote inventory: {out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
