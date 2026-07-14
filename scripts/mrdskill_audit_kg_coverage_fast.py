#!/usr/bin/env python3
"""Fast post-build coverage audit for very large MRDskill KG/RAG artifacts.

The deep audit scans GB-scale CSV/JSONL files and can be slow after a full
supplementary rebuild. This fast audit reuses build reports, inclusion audit
outputs, SQLite metadata, and file headers so the post-build gate remains
practical before testing.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sqlite3
from datetime import datetime
from pathlib import Path
from typing import Any

from mrdskill_project_paths import artifact_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fast audit of latest MRDskill KG/RAG outputs.")
    parser.add_argument("--extract-report-json", default="imports/literature_ingest/20260617_144701/extract_report.json")
    parser.add_argument("--ingest-report-json", default="imports/literature_ingest/20260617_144701/ingest_report.json")
    parser.add_argument("--fresh-rebuild-report-json", default="docs/rag/mrd_fresh_rebuild_report.json")
    parser.add_argument("--update-report-json", default="docs/rag/finalpaper_2010_2026_update_report.json")
    parser.add_argument("--inclusion-audit-json", default="docs/rag/mrd_latest_kg_inclusion_audit.json")
    parser.add_argument("--cohort-report-md", default="docs/qc/supplementary_cohort_rule_build_report.md")
    parser.add_argument("--nodes-augmented", default="data/graph/nodes_augmented.csv")
    parser.add_argument("--relationships-augmented", default="data/graph/relationships_augmented.csv")
    parser.add_argument("--adjudication-nodes", default="data/graph/adjudication_nodes.csv")
    parser.add_argument("--adjudication-relationships", default="data/graph/adjudication_relationships.csv")
    parser.add_argument("--templates-json", default="data/rag/mrd_adjudication_templates.json")
    parser.add_argument("--latest-rag-jsonl", default="data/rag/mrd_latest_unified_rag_chunks.jsonl")
    parser.add_argument("--rule-db", default="data/rag/mrd_latest_kg_rule_index.sqlite")
    parser.add_argument("--out-json", default="docs/rag/mrd_latest_kg_coverage_audit_fast.json")
    parser.add_argument("--out-md", default="docs/rag/mrd_latest_kg_coverage_audit_fast.md")
    return parser.parse_args()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def file_info(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"exists": False, "path": str(path), "bytes": 0, "modified": ""}
    stat = path.stat()
    return {
        "exists": True,
        "path": str(path),
        "bytes": stat.st_size,
        "modified": datetime.fromtimestamp(stat.st_mtime).isoformat(timespec="seconds"),
    }


def csv_header(path: Path) -> list[str]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.reader(handle)
        return next(reader, [])


def jsonl_first(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.strip():
                try:
                    return json.loads(line)
                except json.JSONDecodeError:
                    return {"json_decode_error": True, "line_prefix": line[:200]}
    return {}


def parse_cohort_report(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"exists": False}
    text = path.read_text(encoding="utf-8", errors="replace")
    labels = {
        "sources_scanned": "Sources scanned",
        "tables_accepted": "Tables accepted",
        "variant_events": "Variant events",
        "sample_gene_rows": "Sample-gene rows",
        "sample_feature_rows": "Sample-feature rows",
        "baseline_rows": "Cancer-gene baselines",
        "cooccurrence_rules": "Co-occurrence rules",
        "mutex_rules": "Mutual-exclusion candidate rules",
        "pair_sample_units_total": "Pair-rule sample units total",
        "pair_sample_units_eligible": "Pair-rule sample units eligible",
        "pair_skipped_high_feature": "Pair-rule sample units skipped for high feature load",
        "pair_skipped_missing_cancer": "Pair-rule sample units skipped for missing cancer context",
    }
    result: dict[str, Any] = {"exists": True, "path": str(path)}
    for key, label in labels.items():
        match = re.search(rf"^- {re.escape(label)}:\s*([0-9]+)", text, flags=re.MULTILINE)
        result[key] = int(match.group(1)) if match else None
    top_tables = []
    for match in re.finditer(r"^- ([0-9]+) rows \| ([^|]+) \| ([^|]+) \| `([^`]+)`", text, flags=re.MULTILINE):
        top_tables.append(
            {
                "event_count": int(match.group(1)),
                "paper_id": match.group(2).strip(),
                "source_sheet": match.group(3).strip(),
                "source_file": match.group(4).strip(),
            }
        )
    result["top_tables"] = top_tables[:20]
    return result


def sqlite_summary(path: Path) -> dict[str, Any]:
    info = file_info(path)
    if not info["exists"]:
        return info
    conn = sqlite3.connect(path)
    try:
        info["total_rules"] = conn.execute("SELECT COUNT(*) FROM rules").fetchone()[0]
        info["rule_sources"] = [
            {"source_id": row[0], "row_count": row[1], "path": row[2]}
            for row in conn.execute("SELECT source_id, row_count, path FROM rule_sources ORDER BY source_id")
        ]
        support_counts = {}
        for field in [
            "supporting_paper_ids",
            "supporting_journals",
            "supporting_titles",
            "supporting_dois",
            "supporting_pmids",
            "supporting_source_document_ids",
            "supporting_event_ids",
        ]:
            support_counts[field] = conn.execute(f"SELECT COUNT(*) FROM rules WHERE {field} IS NOT NULL AND {field} != ''").fetchone()[0]
        info["support_field_counts"] = support_counts
        info["support_fields_populated"] = any(count > 0 for count in support_counts.values())
    finally:
        conn.close()
    return info


def successful_steps(report: dict[str, Any]) -> dict[str, Any]:
    steps = report.get("steps", [])
    return {
        "status": report.get("status", ""),
        "step_count": len(steps),
        "successful_steps": sum(1 for step in steps if step.get("returncode") == 0),
        "failed_steps": [
            {
                "command": Path(step.get("command", ["", "unknown"])[1]).name if len(step.get("command", [])) > 1 else "unknown",
                "returncode": step.get("returncode"),
                "started": step.get("started"),
                "completed": step.get("completed"),
                "output_tail": (step.get("output_tail") or "")[-800:],
            }
            for step in steps
            if step.get("returncode") != 0
        ],
    }


def build_audit(args: argparse.Namespace) -> dict[str, Any]:
    paths = {key: artifact_path(value) for key, value in vars(args).items() if key.startswith("out_") is False}
    extract_report = load_json(paths["extract_report_json"])
    ingest_report = load_json(paths["ingest_report_json"])
    fresh_report = load_json(paths["fresh_rebuild_report_json"])
    update_report = load_json(paths["update_report_json"])
    inclusion_audit = load_json(paths["inclusion_audit_json"])
    cohort_report = parse_cohort_report(paths["cohort_report_md"])
    template_payload = load_json(paths["templates_json"])
    rule_db = sqlite_summary(paths["rule_db"])

    files = {
        "nodes_augmented": {**file_info(paths["nodes_augmented"]), "header": csv_header(paths["nodes_augmented"])},
        "relationships_augmented": {**file_info(paths["relationships_augmented"]), "header": csv_header(paths["relationships_augmented"])},
        "adjudication_nodes": {**file_info(paths["adjudication_nodes"]), "header": csv_header(paths["adjudication_nodes"])},
        "adjudication_relationships": {**file_info(paths["adjudication_relationships"]), "header": csv_header(paths["adjudication_relationships"])},
        "latest_rag_jsonl": {**file_info(paths["latest_rag_jsonl"]), "first_chunk": jsonl_first(paths["latest_rag_jsonl"])},
        "templates_json": file_info(paths["templates_json"]),
    }

    expected_files_ok = all(info.get("exists") and info.get("bytes", 0) > 0 for info in files.values())
    adjudication_ok = (
        len(template_payload.get("agents", [])) >= 6
        and len(template_payload.get("question_templates", [])) >= 9
        and len(template_payload.get("decisions", [])) >= 7
        and files["adjudication_nodes"]["exists"]
        and files["adjudication_relationships"]["exists"]
    )
    inclusion_ok = bool(inclusion_audit.get("first_kg_node_subset", {}).get("is_subset")) and bool(
        inclusion_audit.get("first_kg_relationship_subset", {}).get("is_subset")
    )
    latest_rag_total = inclusion_audit.get("latest_rag_counts", {}).get("total", 0)
    rule_db_ok = bool(rule_db.get("exists")) and int(rule_db.get("total_rules") or 0) > 0
    extract_ok = int(extract_report.get("source_files_seen") or 0) > 0 and int(extract_report.get("evidence_rows") or 0) > 0
    ingest_ok = ingest_report.get("status") == "completed" and int(ingest_report.get("merge", {}).get("scored_evidence", {}).get("final_rows") or 0) > 0
    cohort_ok = int(cohort_report.get("variant_events") or 0) > 0 and int(cohort_report.get("cooccurrence_rules") or 0) > 0

    gates = {
        "extract_completed_with_evidence": extract_ok,
        "ingest_merge_completed": ingest_ok,
        "all_expected_output_files_exist": expected_files_ok,
        "first_kg_included_in_augmented_graph": inclusion_ok,
        "latest_rag_has_chunks": int(latest_rag_total or 0) > 0,
        "sqlite_rule_index_has_rows": rule_db_ok,
        "adjudication_layer_present": adjudication_ok,
        "cohort_rule_report_has_counts": cohort_ok,
    }
    overall = "pass" if all(gates.values()) else "fail"

    return {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "mode": "fast_coverage_audit",
        "paths": {key: str(value) for key, value in paths.items()},
        "extract_report": extract_report,
        "ingest_report_summary": {
            "status": ingest_report.get("status"),
            "started": ingest_report.get("started"),
            "completed": ingest_report.get("completed"),
            "scored_evidence_final_rows": ingest_report.get("merge", {}).get("scored_evidence", {}).get("final_rows"),
            "source_documents_final_rows": ingest_report.get("merge", {}).get("sources", {}).get("final_rows"),
            "excel_chunks_final_rows": ingest_report.get("merge", {}).get("excel_chunks", {}).get("final_rows"),
        },
        "fresh_rebuild_report_summary": successful_steps(fresh_report),
        "update_report_summary": successful_steps(update_report),
        "inclusion_audit": inclusion_audit,
        "cohort_rule_report": cohort_report,
        "adjudication_template_counts": {
            "agents": len(template_payload.get("agents", [])),
            "question_templates": len(template_payload.get("question_templates", [])),
            "decisions": len(template_payload.get("decisions", [])),
            "uncertainty_flags": len(template_payload.get("uncertainty_flags", [])),
            "workflow_rounds": len(template_payload.get("workflow_rounds", [])),
        },
        "rule_db": rule_db,
        "files": files,
        "gates": gates,
        "overall_status": overall,
        "notes": [
            "This fast audit intentionally avoids full-table scanning of the 16GB supplementary variant matrix.",
            "Use mrdskill_audit_kg_coverage.py only for deep audits on smaller/capped artifacts or after optimizing it for streaming summaries.",
            "The old fresh/update reports may remain status=failed because their final deep audit step failed, while the KG/RAG build steps succeeded.",
        ],
    }


def write_outputs(audit: dict[str, Any], json_path: Path, md_path: Path) -> None:
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")

    incl = audit["inclusion_audit"]
    cohort = audit["cohort_rule_report"]
    ingest = audit["ingest_report_summary"]
    tmpl = audit["adjudication_template_counts"]
    rule_db = audit["rule_db"]
    lines = [
        "# MRDskill Fast KG/RAG Coverage Audit",
        "",
        f"- Generated at: `{audit['generated_at']}`",
        f"- Overall status: `{audit['overall_status']}`",
        "",
        "## Gates",
        "",
    ]
    for key, value in audit["gates"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(
        [
            "",
            "## Literature Ingest",
            "",
            f"- Extracted source files: {audit['extract_report'].get('source_files_seen')}",
            f"- Extracted evidence rows: {audit['extract_report'].get('evidence_rows')}",
            f"- Excel chunk rows: {audit['extract_report'].get('excel_chunk_rows')}",
            f"- Final scored evidence rows: {ingest.get('scored_evidence_final_rows')}",
            f"- Final source document rows: {ingest.get('source_documents_final_rows')}",
            f"- Final Excel chunk rows: {ingest.get('excel_chunks_final_rows')}",
            "",
            "## KG/RAG Inclusion",
            "",
            f"- `nodes_full` subset of augmented: {incl.get('first_kg_node_subset', {}).get('is_subset')}",
            f"- `relationships_full` subset of augmented: {incl.get('first_kg_relationship_subset', {}).get('is_subset')}",
            f"- Latest RAG total chunks: {incl.get('latest_rag_counts', {}).get('total')}",
            f"- Literature evidence chunks: {incl.get('latest_rag_counts', {}).get('literature_evidence')}",
            f"- Materialized rule chunks: {incl.get('latest_rag_counts', {}).get('materialized_rule_chunks')}",
            "",
            "## Supplementary Rules",
            "",
            f"- Variant events: {cohort.get('variant_events')}",
            f"- Sample-gene rows: {cohort.get('sample_gene_rows')}",
            f"- Sample-feature rows: {cohort.get('sample_feature_rows')}",
            f"- Baseline rows: {cohort.get('baseline_rows')}",
            f"- Co-occurrence rules: {cohort.get('cooccurrence_rules')}",
            f"- Mutual-exclusion rules: {cohort.get('mutex_rules')}",
            f"- SQLite rule rows: {rule_db.get('total_rules')}",
            f"- SQLite support fields populated: {rule_db.get('support_fields_populated')}",
            "",
            "## Adjudication Layer",
            "",
            f"- Agent roles: {tmpl.get('agents')}",
            f"- Question templates: {tmpl.get('question_templates')}",
            f"- Decision categories: {tmpl.get('decisions')}",
            f"- Uncertainty flags: {tmpl.get('uncertainty_flags')}",
            f"- Workflow rounds: {tmpl.get('workflow_rounds')}",
            "",
            "## Notes",
            "",
        ]
    )
    for note in audit["notes"]:
        lines.append(f"- {note}")
    md_path.parent.mkdir(parents=True, exist_ok=True)
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    audit = build_audit(args)
    write_outputs(audit, artifact_path(args.out_json), artifact_path(args.out_md))
    print(json.dumps({"overall_status": audit["overall_status"], "gates": audit["gates"], "out_json": str(artifact_path(args.out_json)), "out_md": str(artifact_path(args.out_md))}, ensure_ascii=False, indent=2))
    if audit["overall_status"] != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
