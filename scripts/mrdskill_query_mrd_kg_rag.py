#!/usr/bin/env python3
"""Simple local retriever over MRD KG-RAG chunks."""

from __future__ import annotations

import argparse
import heapq
import json
import math
import re
from collections import Counter
from pathlib import Path
from typing import Any
import sys
import sqlite3

from mrdskill_project_paths import artifact_path


BLOCKED_RULE_GENE_TERMS = {
    "ASSAY",
    "ALIGNMENT",
    "AND",
    "ACTIONABLE",
    "ARTIFACT",
    "BACKGROUND",
    "BIOMARKER",
    "BLOOD",
    "BUFFY",
    "BY",
    "CANCER",
    "CELL",
    "CELLS",
    "CFDNA",
    "CHIP",
    "CIRCULATING",
    "CLINICAL",
    "CLONAL",
    "COAT",
    "CONTEXT",
    "CONTROL",
    "CTDNA",
    "DATA",
    "DETECTION",
    "DISEASE",
    "DNA",
    "ELIGIBILITY",
    "EXPRESSION",
    "FOR",
    "FROM",
    "GENE",
    "GENES",
    "GERMLINE",
    "HEMATOPOIESIS",
    "HEMATOPOIETIC",
    "HIGH",
    "IN",
    "INFORMED",
    "LOW",
    "MATCHED",
    "MONITORING",
    "MRD",
    "MOLECULAR",
    "OF",
    "ON",
    "OR",
    "ORIGIN",
    "PATHWAY",
    "PATIENT",
    "PATIENTS",
    "PLASMA",
    "PROTEIN",
    "PROTEINS",
    "RECURRENCE",
    "REGULATION",
    "RELAPSE",
    "RESIDUAL",
    "RESPONSE",
    "RESULT",
    "RESULTS",
    "REVIEW",
    "RNA",
    "SAMPLE",
    "SAMPLES",
    "SIGNAL",
    "SOURCE",
    "STUDY",
    "SURVIVAL",
    "TECHNICAL",
    "THE",
    "THERAPY",
    "TO",
    "TRACE",
    "TRACING",
    "TRACKING",
    "TRANSDUCTION",
    "TREATMENT",
    "TRIAL",
    "TUMOR",
    "TUMOUR",
    "USED",
    "USING",
    "VARIANT",
    "VARIANTS",
    "VAF",
    "WBC",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Query local MRD KG-RAG chunks.")
    parser.add_argument("--query", required=True)
    parser.add_argument("--chunks-jsonl", default="data/rag/mrd_latest_unified_rag_chunks.jsonl")
    parser.add_argument("--rule-db", default="data/rag/mrd_latest_kg_rule_index.sqlite")
    parser.add_argument("--top-k", type=int, default=8)
    parser.add_argument("--rule-top-k", type=int, default=5)
    parser.add_argument("--max-chunks", type=int, default=0, help="Optional cap for quick smoke tests; 0 scans all chunks.")
    parser.add_argument("--max-text-chars", type=int, default=1200)
    return parser.parse_args()


def tokenize(text: str) -> list[str]:
    return re.findall(r"[A-Za-z0-9_.:-]+", text.upper())


def load_chunks(path: Path) -> list[dict[str, Any]]:
    rows = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                rows.append(json.loads(line))
    return rows


def iter_chunks(path: Path, max_chunks: int = 0):
    seen = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue
            try:
                yield json.loads(line)
            except json.JSONDecodeError:
                continue
            seen += 1
            if max_chunks and seen >= max_chunks:
                break


def safe_print(text: Any = "") -> None:
    value = "" if text is None else str(text)
    value = value.replace("\r\n", "\n").replace("\r", "\n")
    try:
        print(value)
    except UnicodeEncodeError:
        sys.stdout.buffer.write(value.encode("utf-8", "replace") + b"\n")


def truncate(text: Any, limit: int) -> str:
    value = "" if text is None else str(text)
    value = value.replace("\r\n", "\n").replace("\r", "\n")
    return value if limit <= 0 or len(value) <= limit else value[: limit - 3] + "..."


def score_chunk(query_terms: Counter[str], chunk: dict[str, Any]) -> float:
    text = " ".join([chunk.get("title", ""), chunk.get("text", ""), json.dumps(chunk.get("metadata", {}), ensure_ascii=False)])
    terms = Counter(tokenize(text))
    if not terms:
        return 0.0
    score = 0.0
    for term, q_count in query_terms.items():
        tf = terms.get(term, 0)
        if tf:
            score += (1 + math.log(tf)) * q_count
    phrase = " ".join(query_terms.keys())
    if phrase and phrase in text.upper():
        score += 5
    return score


def existing_chunks_path(path: str) -> Path:
    resolved = artifact_path(path)
    if resolved.exists():
        return resolved
    fallback = artifact_path("data/rag/mrd_unified_rag_chunks.jsonl")
    return fallback


def normalize_gene_terms(values: Any) -> set[str]:
    if values is None:
        return set()
    if isinstance(values, str):
        raw_values = re.split(r"[,;/|\s]+", values)
    else:
        raw_values = list(values)
    genes = set()
    for value in raw_values:
        text = str(value).strip().upper()
        if not text:
            continue
        if text.startswith("GENE:"):
            text = text.split(":", 1)[1]
        if re.fullmatch(r"[A-Z][A-Z0-9.-]{1,15}", text) and text not in BLOCKED_RULE_GENE_TERMS:
            genes.add(text)
    return genes


def row_gene_values(row: dict[str, Any]) -> set[str]:
    values = {
        str(row.get("gene_a") or "").upper(),
        str(row.get("gene_b") or "").upper(),
        str(row.get("feature_a_label") or "").upper(),
        str(row.get("feature_b_label") or "").upper(),
    }
    for key in ("feature_a_id", "feature_b_id"):
        text = str(row.get(key) or "").upper()
        if text.startswith("GENE:"):
            values.add(text.split(":", 1)[1])
    return normalize_gene_terms(values)


def is_blocked_rule_row(row: dict[str, Any], blocked_genes: set[str]) -> bool:
    values = {
        str(row.get("gene_a") or "").upper(),
        str(row.get("gene_b") or "").upper(),
        str(row.get("feature_a_label") or "").upper(),
        str(row.get("feature_b_label") or "").upper(),
    }
    return bool(values & blocked_genes)


def search_rule_db(
    path: Path,
    query: str,
    limit: int,
    allowed_genes: Any = None,
    blocked_genes: Any = None,
    allowed_cancer_types: Any = None,
    require_allowed_gene: bool = False,
) -> list[dict[str, Any]]:
    if not path.exists() or limit <= 0:
        return []
    terms = [term for term in tokenize(query) if len(term) >= 2]
    rule_ids = [term for term in terms if term.startswith(("COOCCUR-", "MUTEX-"))]
    blocked_gene_set = set(BLOCKED_RULE_GENE_TERMS) | normalize_gene_terms(blocked_genes)
    allowed_gene_set = normalize_gene_terms(allowed_genes)
    allowed_cancers = {str(item).strip().lower() for item in (allowed_cancer_types or []) if str(item).strip()}
    gene_terms = [
        term
        for term in terms
        if re.fullmatch(r"[A-Z][A-Z0-9.-]{1,15}", term) and term not in blocked_gene_set
    ]
    if allowed_gene_set:
        gene_terms = [term for term in gene_terms if term in allowed_gene_set]
    feature_terms = []
    for term in terms:
        if term.startswith("LOCUS:"):
            feature_terms.append(term)
        elif re.fullmatch(r"(?:CHR)?[0-9XYMT]+:[0-9][0-9,]*(?::[ACGTN-]+>[ACGTN-]+)?", term):
            feature_terms.append("LOCUS:" + re.sub(r"^CHR", "", term, flags=re.I).replace(",", ""))
    rows: list[dict[str, Any]] = []
    conn = sqlite3.connect(path)
    conn.row_factory = sqlite3.Row
    try:
        if rule_ids:
            placeholders = ",".join(["?"] * len(rule_ids))
            sql = f"SELECT * FROM rules WHERE UPPER(rule_id) IN ({placeholders}) LIMIT ?"
            rows = [dict(row) for row in conn.execute(sql, [*rule_ids, limit])]
        elif gene_terms or feature_terms:
            clauses = []
            params: list[str] = []
            if gene_terms:
                placeholders = ",".join(["?"] * len(gene_terms))
                clauses.append(f"gene_a IN ({placeholders}) OR gene_b IN ({placeholders})")
                params.extend(gene_terms)
                params.extend(gene_terms)
            if feature_terms:
                placeholders = ",".join(["?"] * len(feature_terms))
                clauses.append(f"feature_a_id IN ({placeholders}) OR feature_b_id IN ({placeholders})")
                params.extend(feature_terms)
                params.extend(feature_terms)
            sql = f"SELECT * FROM rules WHERE {' OR '.join(f'({clause})' for clause in clauses)} ORDER BY CAST(pair_count AS REAL) DESC LIMIT ?"
            rows = [dict(row) for row in conn.execute(sql, [*params, max(limit * 20, limit)])]
            query_terms = Counter(terms)
            scored = []
            for row in rows:
                row_genes = row_gene_values(row)
                if is_blocked_rule_row(row, blocked_gene_set):
                    continue
                if allowed_gene_set and not (row_genes & allowed_gene_set):
                    continue
                if require_allowed_gene and allowed_gene_set and not (row_genes & allowed_gene_set):
                    continue
                if allowed_cancers and str(row.get("cancer_type") or "").strip().lower() not in allowed_cancers:
                    continue
                text = " ".join(
                    str(row.get(key, ""))
                    for key in (
                        "rule_id",
                        "rule_type",
                        "rule_status",
                        "cancer_type",
                        "gene_a",
                        "gene_b",
                        "feature_a_id",
                        "feature_b_id",
                        "feature_a_label",
                        "feature_b_label",
                        "source_layer",
                    )
                )
                token_counts = Counter(tokenize(text))
                score = sum(token_counts.get(term, 0) * count for term, count in query_terms.items())
                if allowed_gene_set and row_genes & allowed_gene_set:
                    score += 8 + len(row_genes & allowed_gene_set)
                if str(row.get("rule_status") or "").lower() == "high_confidence_rule":
                    score += 2
                if str(row.get("rule_type") or "").lower() in {"gene_cooccurrence", "gene_mutex"} and not row_genes:
                    score -= 5
                try:
                    score += min(float(row.get("pair_count") or 0), 100) / 100
                except ValueError:
                    pass
                scored.append((score, row))
            rows = [row for score, row in sorted(scored, key=lambda item: item[0], reverse=True)[:limit] if score > 0]
    finally:
        conn.close()
    return rows[:limit]


def main() -> None:
    args = parse_args()
    chunks_path = existing_chunks_path(args.chunks_jsonl)
    query_terms = Counter(tokenize(args.query))
    heap: list[tuple[float, int, dict[str, Any]]] = []
    scanned = 0
    seq = 0
    for chunk in iter_chunks(chunks_path, args.max_chunks):
        scanned += 1
        score = score_chunk(query_terms, chunk)
        if score <= 0:
            continue
        seq += 1
        item = (float(score), seq, chunk)
        if len(heap) < args.top_k:
            heapq.heappush(heap, item)
        elif score > heap[0][0]:
            heapq.heapreplace(heap, item)
    hits = [(score, chunk) for score, _, chunk in sorted(heap, key=lambda item: item[0], reverse=True)]
    safe_print(f"Query: {args.query}")
    safe_print(f"Hits: {len(hits)} / {scanned} chunks scanned")
    safe_print()
    for i, (score, chunk) in enumerate(hits, start=1):
        meta = chunk.get("metadata", {})
        paper = meta.get("supporting_paper_ids") or meta.get("paper_id") or ""
        journal = meta.get("supporting_journals") or meta.get("journal") or ""
        doi = meta.get("supporting_dois") or meta.get("doi") or ""
        pmid = meta.get("supporting_pmids") or meta.get("pmid") or ""
        safe_print(f"[{i}] score={score:.3f} {chunk.get('title')}")
        safe_print(truncate(chunk.get("text"), args.max_text_chars))
        safe_print(f"paper={paper} | journal={journal} | doi={doi} | pmid={pmid}")
        safe_print(f"source_documents={meta.get('supporting_source_document_ids','')}")
        safe_print()
    rule_hits = search_rule_db(artifact_path(args.rule_db), args.query, args.rule_top_k)
    if rule_hits:
        safe_print("Structured rule hits:")
        for i, row in enumerate(rule_hits, start=1):
            safe_print(
                f"[R{i}] {row.get('rule_id')} | {row.get('rule_type')} | {row.get('rule_status')} | "
                f"{row.get('cancer_type')} | "
                f"{row.get('feature_a_label') or row.get('gene_a') or row.get('feature_a_id')}-"
                f"{row.get('feature_b_label') or row.get('gene_b') or row.get('feature_b_id')} | "
                f"pair={row.get('pair_count')}/{row.get('sample_count')} | lift={row.get('lift')} | "
                f"p={row.get('fisher_p_value')} | source={row.get('source_id')}"
            )
        safe_print()


if __name__ == "__main__":
    main()
