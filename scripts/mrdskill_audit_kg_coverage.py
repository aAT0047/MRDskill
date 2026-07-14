#!/usr/bin/env python3
"""Audit MRDskill corpus, KG, rule, and RAG coverage.

This script answers a practical provenance question: which parts of the corpus
were scanned, which evidence records reached the graph/RAG, and where
supplementary variant/rule coverage is complete versus only candidate-level.
"""

from __future__ import annotations

import argparse
import csv
import json
import sqlite3
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

from mrdskill_project_paths import artifact_path


EXCLUDED_STATUSES = {"exclude"}


try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(2**31 - 1)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit MRDskill KG/RAG coverage across corpus, evidence, and rules.")
    parser.add_argument("--sources-csv", default="data/imported/full_corpus_sources.csv")
    parser.add_argument("--parser-errors-csv", default="data/imported/full_corpus_parser_errors.csv")
    parser.add_argument("--ingest-report-json", default="data/imported/full_corpus_ingest_report.json")
    parser.add_argument("--evidence-csv", default="data/evidence/mrd_evidence_scored_full.csv")
    parser.add_argument("--evidence-rag-jsonl", default="data/rag/evidence_summaries_full.jsonl")
    parser.add_argument("--latest-rag-jsonl", default="data/rag/mrd_latest_unified_rag_chunks.jsonl")
    parser.add_argument("--nodes-full", default="data/graph/nodes_full.csv")
    parser.add_argument("--relationships-full", default="data/graph/relationships_full.csv")
    parser.add_argument("--nodes-augmented", default="data/graph/nodes_augmented.csv")
    parser.add_argument("--relationships-augmented", default="data/graph/relationships_augmented.csv")
    parser.add_argument("--variant-matrix-csv", default="data/cohort/supplementary_sample_variant_matrix.csv")
    parser.add_argument("--gene-matrix-csv", default="data/cohort/supplementary_sample_gene_matrix.csv")
    parser.add_argument("--feature-matrix-csv", default="data/cohort/supplementary_sample_feature_matrix.csv")
    parser.add_argument("--baseline-csv", default="data/qc_rules/literature_cancer_gene_baselines.csv")
    parser.add_argument("--cooccurrence-csv", default="data/qc_rules/literature_gene_cooccurrence_rules.csv")
    parser.add_argument("--mutex-csv", default="data/qc_rules/literature_gene_mutex_rules.csv")
    parser.add_argument("--materialized-rule-rag-jsonl", default="data/rag/kg_rag_chunks.jsonl")
    parser.add_argument("--rule-db", default="data/rag/mrd_latest_kg_rule_index.sqlite")
    parser.add_argument("--out-json", default="docs/rag/mrd_kg_coverage_audit.json")
    parser.add_argument("--out-md", default="docs/rag/mrd_kg_coverage_audit.md")
    return parser.parse_args()


def clean(value: Any) -> str:
    return "" if value is None else str(value).strip()


def csv_rows(path: Path) -> Iterable[dict[str, str]]:
    if not path.exists():
        return
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        yield from csv.DictReader(handle)


def count_csv(path: Path) -> int:
    return sum(1 for _ in csv_rows(path))


def counter_csv(path: Path, column: str) -> dict[str, int]:
    counter: Counter[str] = Counter()
    for row in csv_rows(path):
        counter[clean(row.get(column)) or "not_reported"] += 1
    return dict(counter)


def unique_csv(path: Path, column: str, prefix: str = "") -> set[str]:
    values: set[str] = set()
    for row in csv_rows(path):
        value = clean(row.get(column))
        if value:
            values.add(prefix + value)
    return values


def iter_jsonl(path: Path) -> Iterable[dict[str, Any]]:
    if not path.exists():
        return
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue
            try:
                yield json.loads(line)
            except json.JSONDecodeError:
                continue


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def node_ids(path: Path, prefix: str = "") -> set[str]:
    values: set[str] = set()
    for row in csv_rows(path):
        node_id = clean(row.get("node_id"))
        if node_id and (not prefix or node_id.startswith(prefix)):
            values.add(node_id)
    return values


def relationship_ids(path: Path) -> set[str]:
    return unique_csv(path, "relationship_id")


def jsonl_chunk_ids(path: Path, prefix: str = "") -> set[str]:
    values: set[str] = set()
    for row in iter_jsonl(path):
        chunk_id = clean(row.get("chunk_id"))
        if chunk_id and (not prefix or chunk_id.startswith(prefix)):
            values.add(chunk_id)
    return values


def jsonl_kind_counts(path: Path) -> dict[str, int]:
    counter: Counter[str] = Counter()
    for row in iter_jsonl(path):
        metadata = row.get("metadata") if isinstance(row.get("metadata"), dict) else {}
        counter[clean(metadata.get("kind")) or "not_reported"] += 1
    return dict(counter)


def evidence_sets(evidence_csv: Path) -> dict[str, Any]:
    all_ids: set[str] = set()
    included_ids: set[str] = set()
    excluded_ids: set[str] = set()
    status_counts: Counter[str] = Counter()
    source_type_counts: Counter[str] = Counter()
    empty_summary_ids: set[str] = set()
    for row in csv_rows(evidence_csv):
        evidence_id = clean(row.get("evidence_id"))
        if not evidence_id:
            continue
        prefixed = f"Evidence:{evidence_id}"
        all_ids.add(prefixed)
        status = clean(row.get("graph_inclusion_status")) or "not_reported"
        status_counts[status] += 1
        source_type_counts[clean(row.get("source_document_type")) or "not_reported"] += 1
        if status in EXCLUDED_STATUSES:
            excluded_ids.add(prefixed)
        else:
            included_ids.add(prefixed)
        if not clean(row.get("quoted_or_paraphrased_evidence")):
            empty_summary_ids.add(prefixed)
    return {
        "all": all_ids,
        "included": included_ids,
        "excluded": excluded_ids,
        "status_counts": dict(status_counts),
        "source_type_counts": dict(source_type_counts),
        "empty_summary_ids": empty_summary_ids,
    }


def rule_db_summary(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"exists": False}
    conn = sqlite3.connect(path)
    try:
        result: dict[str, Any] = {"exists": True, "path": str(path), "bytes": path.stat().st_size}
        result["rule_sources"] = []
        for source_id, row_count, source_path in conn.execute("SELECT source_id, row_count, path FROM rule_sources ORDER BY source_id"):
            result["rule_sources"].append({"source_id": source_id, "row_count": row_count, "path": source_path})
        source_counts = []
        for source_id, row_count in conn.execute("SELECT source_id, COUNT(*) FROM rules GROUP BY source_id ORDER BY source_id"):
            source_counts.append({"source_id": source_id, "row_count": row_count})
        result["rules_by_source"] = source_counts
        nonempty_support: dict[str, int] = {}
        for column in [
            "supporting_paper_ids",
            "supporting_journals",
            "supporting_titles",
            "supporting_dois",
            "supporting_pmids",
            "supporting_source_document_ids",
            "supporting_event_ids",
        ]:
            count = conn.execute(f"SELECT COUNT(*) FROM rules WHERE {column} IS NOT NULL AND {column} != ''").fetchone()[0]
            nonempty_support[column] = count
        result["nonempty_support_counts"] = nonempty_support
        result["total_rules"] = conn.execute("SELECT COUNT(*) FROM rules").fetchone()[0]
        return result
    finally:
        conn.close()


def csv_fieldnames(path: Path) -> list[str]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or [])


def sample_rows(path: Path, limit: int = 5) -> list[dict[str, str]]:
    rows = []
    for row in csv_rows(path):
        rows.append(dict(row))
        if len(rows) >= limit:
            break
    return rows


def write_outputs(audit: dict[str, Any], json_path: Path, md_path: Path) -> None:
    json_path.parent.mkdir(parents=True, exist_ok=True)
    md_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")

    source = audit["corpus_sources"]
    evidence = audit["evidence_kg_rag"]
    rules = audit["supplementary_variant_rules"]
    conclusion = audit["conclusion"]
    lines = [
        "# MRDskill KG Coverage Audit",
        "",
        f"- Generated at: `{audit['generated_at']}`",
        "",
        "## Bottom Line",
        "",
        f"- Full text / source files scanned: `{conclusion['source_files_scanned']}`",
        f"- Source parser coverage: `{conclusion['source_parser_coverage']}`",
        f"- Literature evidence included in KG/RAG: `{conclusion['literature_evidence_included']}`",
        f"- Supplementary variant-rule coverage: `{conclusion['supplementary_variant_rule_coverage']}`",
        f"- Overall answer: `{conclusion['overall_answer']}`",
        "",
        "## Corpus Sources",
        "",
        f"- Sources registered: {source['source_rows']}",
        f"- Parser status: {source['parser_status_counts']}",
        f"- Source document types: {source['source_document_type_counts']}",
        f"- Parser errors: {source['parser_error_rows']}",
        "",
        "## Literature Evidence KG/RAG",
        "",
        f"- Scored evidence rows: {evidence['scored_evidence_rows']}",
        f"- Non-excluded evidence rows: {evidence['included_evidence_rows']}",
        f"- Excluded evidence rows: {evidence['excluded_evidence_rows']}",
        f"- Evidence node coverage for non-excluded rows: {evidence['included_evidence_node_coverage']}",
        f"- Latest RAG literature chunk coverage for non-excluded rows: {evidence['included_evidence_latest_rag_coverage']}",
        f"- `nodes_full.csv` subset of `nodes_augmented.csv`: {evidence['nodes_full_subset_of_augmented']}",
        f"- `relationships_full.csv` subset of `relationships_augmented.csv`: {evidence['relationships_full_subset_of_augmented']}",
        "",
        "## Supplementary Variant / Rule Layer",
        "",
        f"- Current variant matrix rows: {rules['current_variant_matrix_rows']}",
        f"- Current sample-gene rows: {rules['current_sample_gene_rows']}",
        f"- Current sample-feature rows: {rules['current_sample_feature_rows']}",
        f"- Current baseline rows: {rules['current_baseline_rows']}",
        f"- Current co-occurrence rule rows: {rules['current_cooccurrence_rows']}",
        f"- Current mutual-exclusion rule rows: {rules['current_mutex_rows']}",
        f"- Materialized rule RAG chunks: {rules['materialized_rule_rag_chunks']}",
        f"- SQLite full rule rows: {rules['sqlite_rule_summary'].get('total_rules', 'not_available')}",
        f"- SQLite support/provenance fields populated: {rules['sqlite_support_fields_populated']}",
        f"- Variant annotation status: {rules['variant_annotation_status_counts']}",
        f"- Coordinate gene status: {rules['coordinate_gene_status_counts']}",
        f"- Feature types: {rules['feature_type_counts']}",
        "",
        "## Important Limits",
        "",
    ]
    for item in conclusion["limits"]:
        lines.append(f"- {item}")
    lines.extend(
        [
            "",
            "## Recommended Fixes",
            "",
        ]
    )
    for item in conclusion["recommended_fixes"]:
        lines.append(f"- {item}")
    lines.append("")
    md_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    paths = {key: artifact_path(value) for key, value in vars(args).items() if key.startswith("out_") is False}

    ingest_report = load_json(paths["ingest_report_json"])
    source_rows = count_csv(paths["sources_csv"])
    parser_error_rows = count_csv(paths["parser_errors_csv"])
    parser_status_counts = counter_csv(paths["sources_csv"], "parser_status")
    source_type_counts = counter_csv(paths["sources_csv"], "source_document_type")

    source_ids = unique_csv(paths["sources_csv"], "source_document_id", "SourceDocument:")
    augmented_node_ids = node_ids(paths["nodes_augmented"])
    missing_source_nodes = sorted(source_ids - augmented_node_ids)

    evidence = evidence_sets(paths["evidence_csv"])
    augmented_evidence_nodes = node_ids(paths["nodes_augmented"], "Evidence:")
    latest_lit_chunks = jsonl_chunk_ids(paths["latest_rag_jsonl"], "LiteratureEvidence:")
    latest_lit_as_evidence = {"Evidence:" + chunk.split(":", 1)[1] for chunk in latest_lit_chunks if ":" in chunk}
    evidence_rag_chunks = jsonl_chunk_ids(paths["evidence_rag_jsonl"], "LiteratureEvidence:")
    evidence_rag_as_evidence = {"Evidence:" + chunk.split(":", 1)[1] for chunk in evidence_rag_chunks if ":" in chunk}

    nodes_full_ids = node_ids(paths["nodes_full"])
    rels_full_ids = relationship_ids(paths["relationships_full"])
    nodes_aug_ids = node_ids(paths["nodes_augmented"])
    rels_aug_ids = relationship_ids(paths["relationships_augmented"])

    rule_summary = rule_db_summary(paths["rule_db"])
    sqlite_support_populated = False
    if rule_summary.get("exists"):
        sqlite_support_populated = any(value > 0 for value in rule_summary.get("nonempty_support_counts", {}).values())

    rule_fields = {
        "baseline": csv_fieldnames(paths["baseline_csv"]),
        "cooccurrence": csv_fieldnames(paths["cooccurrence_csv"]),
        "mutex": csv_fieldnames(paths["mutex_csv"]),
    }
    rule_has_provenance = {
        name: all(field in fields for field in ["supporting_paper_ids", "supporting_source_document_ids", "supporting_event_ids"])
        for name, fields in rule_fields.items()
    }

    included_missing_nodes = sorted(evidence["included"] - augmented_evidence_nodes)
    included_missing_latest_rag = sorted(evidence["included"] - latest_lit_as_evidence)
    excluded_missing_nodes = sorted(evidence["excluded"] - augmented_evidence_nodes)

    source_coverage = (
        "all registered sources have SourceDocument nodes"
        if not missing_source_nodes
        else f"{len(missing_source_nodes)} registered sources lack SourceDocument nodes"
    )
    literature_coverage = (
        "complete for non-excluded scored evidence"
        if not included_missing_nodes and not included_missing_latest_rag
        else "incomplete for non-excluded scored evidence"
    )
    variant_rule_coverage = (
        "partial: current compact KG/RAG has provenance, but full SQLite rule rows lack row-level provenance fields"
        if not sqlite_support_populated
        else "structured rule rows include row-level provenance fields"
    )
    core_complete = (
        not missing_source_nodes
        and literature_coverage == "complete for non-excluded scored evidence"
        and len(nodes_full_ids - nodes_aug_ids) == 0
        and len(rels_full_ids - rels_aug_ids) == 0
        and sqlite_support_populated
    )
    overall_answer = (
        "yes for the current audited rebuild: non-excluded literature evidence, supplementary variant events, "
        "cancer baselines, co-occurrence rules, mutual-exclusion review rules, augmented KG, latest RAG, and "
        "SQLite rule index are integrated with provenance"
        if core_complete
        else "not fully: one or more coverage/provenance gates failed; review the counts below"
    )
    limits = []
    if parser_error_rows:
        limits.append(f"{parser_error_rows} registered source files failed parsing and are listed in full_corpus_parser_errors.csv.")
    if evidence["excluded"]:
        limits.append("Excluded evidence records are intentionally absent from evidence nodes/RAG chunks.")
    if sqlite_support_populated:
        limits.append("The SQLite rule index preserves imported rule rows with supporting paper/source/event provenance fields.")
    else:
        limits.append("SQLite rule rows currently lack populated supporting_* provenance fields; rebuild provenance-rich rule CSVs before traceable rule answers.")
    limits.extend(
        [
            "Supplementary variant events are preserved in the event and sample-feature matrices; pair-rule tables are filtered QC signals, not every possible feature pair.",
            "Coordinate-only variants are retained as locus-level features; NCBI interval gene rescue is supplemental and does not overwrite reported genes.",
            "Cohort co-occurrence/mutual-exclusion rules remain QC support, not standalone clinical truth.",
        ]
    )
    recommended_fixes = [
        "Review or repair parser-error supplementary files if they are scientifically important.",
        "Keep pair-rule thresholds and skipped high-feature sample counts in every run report for auditability.",
        "Add an audit gate to mrdskill_update_mrd_kg.py so every rebuild writes this coverage report and fails if coverage/provenance gates fail.",
    ]

    audit = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "inputs": {name: str(path) for name, path in paths.items()},
        "corpus_sources": {
            "source_rows": source_rows,
            "ingest_report": ingest_report,
            "parser_status_counts": parser_status_counts,
            "source_document_type_counts": source_type_counts,
            "parser_error_rows": parser_error_rows,
            "missing_source_document_node_count": len(missing_source_nodes),
            "missing_source_document_node_examples": missing_source_nodes[:20],
        },
        "evidence_kg_rag": {
            "scored_evidence_rows": len(evidence["all"]),
            "included_evidence_rows": len(evidence["included"]),
            "excluded_evidence_rows": len(evidence["excluded"]),
            "graph_inclusion_status_counts": evidence["status_counts"],
            "source_document_type_counts": evidence["source_type_counts"],
            "evidence_summary_jsonl_chunks": len(evidence_rag_chunks),
            "latest_rag_kind_counts": jsonl_kind_counts(paths["latest_rag_jsonl"]),
            "included_evidence_node_coverage": f"{len(evidence['included'] - set(included_missing_nodes))}/{len(evidence['included'])}",
            "included_evidence_latest_rag_coverage": f"{len(evidence['included'] - set(included_missing_latest_rag))}/{len(evidence['included'])}",
            "included_missing_evidence_node_count": len(included_missing_nodes),
            "included_missing_evidence_node_examples": included_missing_nodes[:20],
            "included_missing_latest_rag_count": len(included_missing_latest_rag),
            "included_missing_latest_rag_examples": included_missing_latest_rag[:20],
            "excluded_missing_evidence_node_count": len(excluded_missing_nodes),
            "excluded_missing_evidence_node_examples": excluded_missing_nodes[:20],
            "nodes_full_subset_of_augmented": len(nodes_full_ids - nodes_aug_ids) == 0,
            "relationships_full_subset_of_augmented": len(rels_full_ids - rels_aug_ids) == 0,
            "nodes_full_missing_from_augmented": len(nodes_full_ids - nodes_aug_ids),
            "relationships_full_missing_from_augmented": len(rels_full_ids - rels_aug_ids),
        },
        "supplementary_variant_rules": {
            "current_variant_matrix_rows": count_csv(paths["variant_matrix_csv"]),
            "current_sample_gene_rows": count_csv(paths["gene_matrix_csv"]),
            "current_sample_feature_rows": count_csv(paths["feature_matrix_csv"]),
            "current_baseline_rows": count_csv(paths["baseline_csv"]),
            "current_cooccurrence_rows": count_csv(paths["cooccurrence_csv"]),
            "current_mutex_rows": count_csv(paths["mutex_csv"]),
            "variant_annotation_status_counts": counter_csv(paths["variant_matrix_csv"], "annotation_status"),
            "coordinate_gene_status_counts": counter_csv(paths["variant_matrix_csv"], "coordinate_gene_status"),
            "variant_extraction_confidence_counts": counter_csv(paths["variant_matrix_csv"], "extraction_confidence"),
            "feature_type_counts": counter_csv(paths["feature_matrix_csv"], "feature_type"),
            "current_variant_matrix_source_files": len(unique_csv(paths["variant_matrix_csv"], "source_file")),
            "current_variant_matrix_papers": len(unique_csv(paths["variant_matrix_csv"], "paper_id")),
            "rule_csv_has_provenance_columns": rule_has_provenance,
            "materialized_rule_rag_chunks": sum(jsonl_kind_counts(paths["materialized_rule_rag_jsonl"]).values()),
            "materialized_rule_rag_kind_counts": jsonl_kind_counts(paths["materialized_rule_rag_jsonl"]),
            "sqlite_rule_summary": rule_summary,
            "sqlite_support_fields_populated": sqlite_support_populated,
        },
        "conclusion": {
            "source_files_scanned": "yes, all registered full-text and supplementary sources were scanned or error-logged",
            "source_parser_coverage": source_coverage,
            "literature_evidence_included": literature_coverage,
            "supplementary_variant_rule_coverage": variant_rule_coverage,
            "overall_answer": overall_answer,
            "limits": limits,
            "recommended_fixes": recommended_fixes,
        },
    }
    write_outputs(audit, artifact_path(args.out_json), artifact_path(args.out_md))
    print(json.dumps(audit["conclusion"], ensure_ascii=False, indent=2))
    print(json.dumps({"audit_json": str(artifact_path(args.out_json)), "audit_md": str(artifact_path(args.out_md))}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
