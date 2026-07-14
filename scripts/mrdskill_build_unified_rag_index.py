#!/usr/bin/env python3
"""Build a unified RAG index from literature evidence and cohort-rule chunks."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from mrdskill_project_paths import artifact_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Merge MRD literature and cohort-rule RAG chunks.")
    parser.add_argument("--evidence-jsonl", default="data/rag/evidence_summaries_full.jsonl")
    parser.add_argument("--rule-jsonl", default="data/rag/kg_rag_chunks.jsonl")
    parser.add_argument("--out-jsonl", default="data/rag/mrd_unified_rag_chunks.jsonl")
    parser.add_argument("--out-report-md", default="docs/rag/mrd_unified_rag_usage.md")
    return parser.parse_args()


def clean(value: Any) -> str:
    return "" if value is None else str(value).strip()


def iter_jsonl(path: Path):
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
    evidence_id = clean(row.get("evidence_id"))
    title = clean(row.get("title")) or evidence_id
    text = clean(row.get("summary"))
    confidence = clean(row.get("confidence_score"))
    status = clean(row.get("graph_inclusion_status"))
    doi = clean(row.get("doi"))
    pmid = clean(row.get("pmid"))
    return {
        "chunk_id": f"LiteratureEvidence:{evidence_id}",
        "title": f"Literature evidence {evidence_id}: {title}",
        "text": f"{title}\n{text}\nConfidence: {confidence}; graph status: {status}; DOI: {doi}; PMID: {pmid}.",
        "metadata": {
            "kind": "literature_evidence",
            "evidence_id": evidence_id,
            "paper_id": clean(row.get("paper_id")),
            "title": title,
            "confidence_score": confidence,
            "graph_inclusion_status": status,
            "audit_priority": clean(row.get("audit_priority")),
            "doi": doi,
            "pmid": pmid,
        },
    }


def rule_chunk(row: dict[str, Any]) -> dict[str, Any]:
    metadata = row.get("metadata", {})
    metadata = {**metadata, "kind": metadata.get("kind", "cohort_rule")}
    return {
        "chunk_id": clean(row.get("chunk_id")),
        "title": clean(row.get("title")),
        "text": clean(row.get("text")),
        "metadata": metadata,
    }


def write_report(path: Path, counts: dict[str, int], out_jsonl: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = f"""# MRD Unified RAG Index

Unified index: `{out_jsonl}`

## Counts

- Literature evidence chunks: {counts['literature_evidence']}
- Cohort/rule chunks: {counts['cohort_rule']}
- Total chunks: {counts['total']}

## Query

```powershell
python scripts\\mrdskill_query_mrd_kg_rag.py --chunks-jsonl {out_jsonl} --query "Bladder ASXL1 TET2 cooccur"
python scripts\\mrdskill_query_mrd_kg_rag.py --chunks-jsonl {out_jsonl} --query "TP53 ctDNA MRD recurrence"
```

## Intended Use

Retrieve relevant chunks, then answer with source traceability. Rule chunks are cohort-derived and should be treated as QC support, not standalone clinical truth.
"""
    path.write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    out_path = artifact_path(args.out_jsonl)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    counts = {"literature_evidence": 0, "cohort_rule": 0, "total": 0}
    with out_path.open("w", encoding="utf-8") as out:
        for row in iter_jsonl(artifact_path(args.evidence_jsonl)) or []:
            chunk = evidence_chunk(row)
            out.write(json.dumps(chunk, ensure_ascii=False, sort_keys=True) + "\n")
            counts["literature_evidence"] += 1
            counts["total"] += 1
        for row in iter_jsonl(artifact_path(args.rule_jsonl)) or []:
            chunk = rule_chunk(row)
            out.write(json.dumps(chunk, ensure_ascii=False, sort_keys=True) + "\n")
            counts["cohort_rule"] += 1
            counts["total"] += 1
    write_report(artifact_path(args.out_report_md), counts, out_path)
    print(json.dumps(counts, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
