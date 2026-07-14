#!/usr/bin/env python3
"""Build the latest MRDskill KG-RAG inclusion index.

The large supplementary rule tables are too big to turn every row into a JSONL
RAG chunk. This script preserves every rule row in a structured SQLite index,
then writes a compact RAG manifest that points retrieval to the full tables.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import sqlite3
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

from mrdskill_project_paths import artifact_path, artifact_root


RULE_FIELDS = [
    "rule_id",
    "rule_type",
    "rule_status",
    "cancer_type",
    "gene_a",
    "gene_b",
    "feature_a_id",
    "feature_b_id",
    "feature_a_type",
    "feature_b_type",
    "feature_a_label",
    "feature_b_label",
    "sample_count",
    "gene_a_count",
    "gene_b_count",
    "pair_count",
    "gene_a_frequency",
    "gene_b_frequency",
    "pair_frequency",
    "lift",
    "odds_ratio",
    "fisher_p_value",
    "source_layer",
    "supporting_paper_ids",
    "supporting_journals",
    "supporting_titles",
    "supporting_dois",
    "supporting_pmids",
    "supporting_source_document_ids",
    "supporting_event_ids",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build latest KG inclusion audit and structured rule index.")
    parser.add_argument("--nodes-full", default="data/graph/nodes_full.csv")
    parser.add_argument("--relationships-full", default="data/graph/relationships_full.csv")
    parser.add_argument("--nodes-augmented", default="data/graph/nodes_augmented.csv")
    parser.add_argument("--relationships-augmented", default="data/graph/relationships_augmented.csv")
    parser.add_argument("--evidence-jsonl", default="data/rag/evidence_summaries_full.jsonl")
    parser.add_argument("--existing-rule-jsonl", default="data/rag/kg_rag_chunks.jsonl")
    parser.add_argument("--out-db", default="data/rag/mrd_latest_kg_rule_index.sqlite")
    parser.add_argument("--out-jsonl", default="data/rag/mrd_latest_unified_rag_chunks.jsonl")
    parser.add_argument("--out-audit-json", default="docs/rag/mrd_latest_kg_inclusion_audit.json")
    parser.add_argument("--out-audit-md", default="docs/rag/mrd_latest_kg_inclusion_audit.md")
    parser.add_argument("--batch-size", type=int, default=50000)
    parser.add_argument(
        "--rule-source",
        action="append",
        default=[],
        help="Optional rule CSV source_id=relative/path.csv. Defaults include current and legacy rule tables.",
    )
    return parser.parse_args()


def default_rule_sources() -> list[tuple[str, Path]]:
    return [
        ("current_cooccurrence", artifact_path("data/qc_rules/literature_gene_cooccurrence_rules.csv")),
        ("current_mutex", artifact_path("data/qc_rules/literature_gene_mutex_rules.csv")),
    ]


def parse_rule_source(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise ValueError("--rule-source must be source_id=path")
    source_id, path = value.split("=", 1)
    return source_id.strip(), artifact_path(path.strip())


def count_lines(path: Path) -> int:
    count = 0
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        for _ in handle:
            count += 1
    return count


def first_column_set(path: Path, column: str) -> set[str]:
    values: set[str] = set()
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            value = (row.get(column) or "").strip()
            if value:
                values.add(value)
    return values


def subset_audit(source: Path, target: Path, key_column: str) -> dict[str, Any]:
    source_ids = first_column_set(source, key_column)
    target_ids = first_column_set(target, key_column)
    missing = sorted(source_ids - target_ids)
    return {
        "source": str(source),
        "target": str(target),
        "key_column": key_column,
        "source_unique_ids": len(source_ids),
        "target_unique_ids": len(target_ids),
        "missing_count": len(missing),
        "missing_examples": missing[:20],
        "is_subset": len(missing) == 0,
    }


def connect_db(path: Path) -> sqlite3.Connection:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")
    conn.execute("PRAGMA temp_store=MEMORY")
    conn.execute("PRAGMA locking_mode=EXCLUSIVE")
    conn.execute("PRAGMA cache_size=-1000000")
    return conn


def create_schema(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        CREATE TABLE rule_sources (
          source_id TEXT PRIMARY KEY,
          path TEXT NOT NULL,
          bytes INTEGER NOT NULL,
          last_write_time TEXT NOT NULL,
          row_count INTEGER NOT NULL,
          imported_at TEXT NOT NULL
        );

        CREATE TABLE rules (
          source_id TEXT NOT NULL,
          row_num INTEGER NOT NULL,
          rule_id TEXT,
          rule_type TEXT,
          rule_status TEXT,
          cancer_type TEXT,
          gene_a TEXT,
          gene_b TEXT,
          feature_a_id TEXT,
          feature_b_id TEXT,
          feature_a_type TEXT,
          feature_b_type TEXT,
          feature_a_label TEXT,
          feature_b_label TEXT,
          sample_count TEXT,
          gene_a_count TEXT,
          gene_b_count TEXT,
          pair_count TEXT,
          gene_a_frequency TEXT,
          gene_b_frequency TEXT,
          pair_frequency TEXT,
          lift TEXT,
          odds_ratio TEXT,
          fisher_p_value TEXT,
          source_layer TEXT,
          supporting_paper_ids TEXT,
          supporting_journals TEXT,
          supporting_titles TEXT,
          supporting_dois TEXT,
          supporting_pmids TEXT,
          supporting_source_document_ids TEXT,
          supporting_event_ids TEXT,
          PRIMARY KEY (source_id, row_num)
        );
        """
    )


def row_values(source_id: str, row_num: int, row: dict[str, str]) -> tuple[Any, ...]:
    return tuple([source_id, row_num] + [(row.get(field) or "").strip() for field in RULE_FIELDS])


def import_rule_csv(conn: sqlite3.Connection, source_id: str, path: Path, batch_size: int) -> dict[str, Any]:
    if not path.exists():
        return {"source_id": source_id, "path": str(path), "exists": False, "row_count": 0}
    insert_sql = f"""
        INSERT INTO rules (
          source_id, row_num, {", ".join(RULE_FIELDS)}
        ) VALUES ({", ".join(["?"] * (2 + len(RULE_FIELDS)))})
    """
    row_count = 0
    rule_type_counts: dict[str, int] = {}
    batch: list[tuple[Any, ...]] = []
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle)
        for row_count, row in enumerate(reader, start=1):
            rule_type = (row.get("rule_type") or "not_reported").strip()
            rule_type_counts[rule_type] = rule_type_counts.get(rule_type, 0) + 1
            batch.append(row_values(source_id, row_count, row))
            if len(batch) >= batch_size:
                conn.executemany(insert_sql, batch)
                conn.commit()
                batch.clear()
    if batch:
        conn.executemany(insert_sql, batch)
        conn.commit()
    stat = path.stat()
    conn.execute(
        "INSERT INTO rule_sources VALUES (?, ?, ?, ?, ?, ?)",
        (source_id, str(path), stat.st_size, datetime.fromtimestamp(stat.st_mtime).isoformat(timespec="seconds"), row_count, datetime.now().isoformat(timespec="seconds")),
    )
    conn.commit()
    return {
        "source_id": source_id,
        "path": str(path),
        "exists": True,
        "bytes": stat.st_size,
        "row_count": row_count,
        "rule_type_counts": rule_type_counts,
    }


def create_indexes(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        CREATE INDEX idx_rules_rule_id ON rules(rule_id);
        CREATE INDEX idx_rules_gene_a ON rules(gene_a);
        CREATE INDEX idx_rules_gene_b ON rules(gene_b);
        CREATE INDEX idx_rules_feature_a ON rules(feature_a_id);
        CREATE INDEX idx_rules_feature_b ON rules(feature_b_id);
        CREATE INDEX idx_rules_cancer ON rules(cancer_type);
        CREATE INDEX idx_rules_pair ON rules(cancer_type, gene_a, gene_b, rule_type);
        CREATE INDEX idx_rules_feature_pair ON rules(cancer_type, feature_a_id, feature_b_id, rule_type);
        CREATE INDEX idx_rules_status ON rules(rule_status);
        """
    )
    conn.commit()


def support_field_counts(conn: sqlite3.Connection) -> dict[str, int]:
    counts: dict[str, int] = {}
    for field in [
        "supporting_paper_ids",
        "supporting_journals",
        "supporting_titles",
        "supporting_dois",
        "supporting_pmids",
        "supporting_source_document_ids",
        "supporting_event_ids",
    ]:
        counts[field] = conn.execute(f"SELECT COUNT(*) FROM rules WHERE {field} IS NOT NULL AND {field} != ''").fetchone()[0]
    return counts


def iter_jsonl(path: Path) -> Iterable[dict[str, Any]]:
    if not path.exists():
        return
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.strip():
                try:
                    yield json.loads(line)
                except json.JSONDecodeError:
                    continue


def evidence_chunk(row: dict[str, Any]) -> dict[str, Any]:
    evidence_id = str(row.get("evidence_id") or "").strip()
    title = str(row.get("title") or evidence_id).strip()
    text = str(row.get("summary") or "").strip()
    return {
        "chunk_id": f"LiteratureEvidence:{evidence_id}",
        "title": f"Literature evidence {evidence_id}: {title}",
        "text": text,
        "metadata": {
            "kind": "literature_evidence",
            "evidence_id": evidence_id,
            "paper_id": row.get("paper_id", ""),
            "title": title,
            "confidence_score": row.get("confidence_score", ""),
            "graph_inclusion_status": row.get("graph_inclusion_status", ""),
            "audit_priority": row.get("audit_priority", ""),
            "doi": row.get("doi", ""),
            "pmid": row.get("pmid", ""),
        },
    }


def write_latest_rag(path: Path, evidence_jsonl: Path, existing_rule_jsonl: Path, rule_sources: list[dict[str, Any]], db_path: Path) -> dict[str, int]:
    path.parent.mkdir(parents=True, exist_ok=True)
    counts = {"literature_evidence": 0, "materialized_rule_chunks": 0, "rule_source_manifest": 0, "total": 0}
    seen_chunks: set[str] = set()
    with path.open("w", encoding="utf-8") as out:
        for row in iter_jsonl(evidence_jsonl) or []:
            chunk = evidence_chunk(row)
            out.write(json.dumps(chunk, ensure_ascii=False, sort_keys=True) + "\n")
            seen_chunks.add(chunk["chunk_id"])
            counts["literature_evidence"] += 1
            counts["total"] += 1
        for chunk in iter_jsonl(existing_rule_jsonl) or []:
            chunk_id = str(chunk.get("chunk_id") or "")
            if not chunk_id or chunk_id in seen_chunks:
                continue
            out.write(json.dumps(chunk, ensure_ascii=False, sort_keys=True) + "\n")
            seen_chunks.add(chunk_id)
            counts["materialized_rule_chunks"] += 1
            counts["total"] += 1
        for source in rule_sources:
            chunk = {
                "chunk_id": f"RuleSource:{source['source_id']}",
                "title": f"Full supplementary rule table {source['source_id']}",
                "text": (
                    f"Full rule table {source['source_id']} contains {source.get('row_count', 0)} rows. "
                    f"Rows are available in the structured SQLite rule index at {db_path}. "
                    f"Source CSV: {source.get('path')}. This manifest preserves full-table retrieval without loading every row as a text chunk."
                ),
                "metadata": {"kind": "rule_source_manifest", "sqlite_rule_db": str(db_path), **source},
            }
            out.write(json.dumps(chunk, ensure_ascii=False, sort_keys=True) + "\n")
            counts["rule_source_manifest"] += 1
            counts["total"] += 1
    return counts


def write_audit(audit_json: Path, audit_md: Path, audit: dict[str, Any]) -> None:
    audit_json.parent.mkdir(parents=True, exist_ok=True)
    audit_json.write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    lines = [
        "# MRDskill Latest KG Inclusion Audit",
        "",
        f"- Generated at: `{audit['generated_at']}`",
        f"- Structured rule DB: `{audit['rule_db']}`",
        f"- Latest RAG JSONL: `{audit['latest_rag_jsonl']}`",
        "",
        "## First KG Inclusion",
        "",
        f"- `nodes_full.csv` subset of `nodes_augmented.csv`: `{audit['first_kg_node_subset']['is_subset']}`; missing: {audit['first_kg_node_subset']['missing_count']}",
        f"- `relationships_full.csv` subset of `relationships_augmented.csv`: `{audit['first_kg_relationship_subset']['is_subset']}`; missing: {audit['first_kg_relationship_subset']['missing_count']}",
        "",
        "## Rule Sources Loaded",
        "",
    ]
    for source in audit["rule_sources"]:
        lines.append(f"- `{source['source_id']}`: {source.get('row_count', 0)} rows from `{source.get('path')}`")
    lines.extend(
        [
            "",
            "## Latest RAG Counts",
            "",
            f"- Literature evidence chunks: {audit['latest_rag_counts']['literature_evidence']}",
            f"- Existing materialized rule chunks: {audit['latest_rag_counts']['materialized_rule_chunks']}",
            f"- Full rule source manifest chunks: {audit['latest_rag_counts']['rule_source_manifest']}",
            f"- Total latest RAG chunks: {audit['latest_rag_counts']['total']}",
            "",
            "## Rule Provenance Fields",
            "",
            f"- SQLite rule rows with any populated `supporting_*` field: `{audit.get('rule_support_fields_populated', False)}`",
        ]
    )
    for field, count in audit.get("rule_support_field_counts", {}).items():
        lines.append(f"- `{field}` populated rows: {count}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- The first literature KG is included in the augmented graph when both subset checks are true.",
            "- The second large rule build is included as full structured rows in the SQLite rule index.",
            "- The latest JSONL RAG keeps compact manifest chunks and delegates row-level rule retrieval to SQLite.",
            "- If `supporting_*` field counts are zero, full SQLite rule rows are statistical rows without row-level paper/source/event provenance; use the materialized rule chunks or rebuild provenance-rich rule CSVs for traceable rule answers.",
        ]
    )
    audit_md.parent.mkdir(parents=True, exist_ok=True)
    audit_md.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    rule_sources = [parse_rule_source(value) for value in args.rule_source] if args.rule_source else default_rule_sources()
    db_path = artifact_path(args.out_db)
    temp_db = db_path.with_suffix(".tmp.sqlite")
    if temp_db.exists():
        temp_db.unlink()
    conn = connect_db(temp_db)
    create_schema(conn)
    imported_sources = []
    for source_id, path in rule_sources:
        imported_sources.append(import_rule_csv(conn, source_id, path, args.batch_size))
    create_indexes(conn)
    support_counts = support_field_counts(conn)
    conn.close()
    if db_path.exists():
        backup = db_path.with_suffix(f".previous_{datetime.now().strftime('%Y%m%d_%H%M%S')}.sqlite")
        shutil.move(str(db_path), str(backup))
    shutil.move(str(temp_db), str(db_path))

    nodes_audit = subset_audit(artifact_path(args.nodes_full), artifact_path(args.nodes_augmented), "node_id")
    rels_audit = subset_audit(artifact_path(args.relationships_full), artifact_path(args.relationships_augmented), "relationship_id")
    rag_counts = write_latest_rag(artifact_path(args.out_jsonl), artifact_path(args.evidence_jsonl), artifact_path(args.existing_rule_jsonl), imported_sources, db_path)

    audit = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "rule_db": str(db_path),
        "latest_rag_jsonl": str(artifact_path(args.out_jsonl)),
        "first_kg_node_subset": nodes_audit,
        "first_kg_relationship_subset": rels_audit,
        "rule_sources": imported_sources,
        "latest_rag_counts": rag_counts,
        "rule_support_field_counts": support_counts,
        "rule_support_fields_populated": any(count > 0 for count in support_counts.values()),
    }
    write_audit(artifact_path(args.out_audit_json), artifact_path(args.out_audit_md), audit)
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
