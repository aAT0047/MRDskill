#!/usr/bin/env python3
"""Add supplementary cohort rules to the MRD KG and create GraphRAG chunks."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path
from typing import Any
import hashlib

from mrdskill_project_paths import artifact_path


NODE_FIELDS = ["node_id", "node_label", "name", "normalized_name", "type", "source", "ontology_id", "properties_json"]
REL_FIELDS = [
    "relationship_id", "start_node_id", "end_node_id", "relationship_type", "confidence_score", "evidence_strength_score",
    "evidence_ids", "source_pmids", "source_dois", "graph_inclusion_status", "audit_priority", "properties_json",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Augment MRD KG with supplementary cohort baseline and association rules.")
    parser.add_argument("--base-nodes", default="data/graph/nodes_full.csv")
    parser.add_argument("--base-relationships", default="data/graph/relationships_full.csv")
    parser.add_argument("--baseline-csv", default="data/qc_rules/literature_cancer_gene_baselines.csv")
    parser.add_argument("--cooccurrence-csv", default="data/qc_rules/literature_gene_cooccurrence_rules.csv")
    parser.add_argument("--mutex-csv", default="data/qc_rules/literature_gene_mutex_rules.csv")
    parser.add_argument("--out-nodes-extension", default="data/graph/cohort_rule_nodes.csv")
    parser.add_argument("--out-relationships-extension", default="data/graph/cohort_rule_relationships.csv")
    parser.add_argument("--out-nodes-augmented", default="data/graph/nodes_augmented.csv")
    parser.add_argument("--out-relationships-augmented", default="data/graph/relationships_augmented.csv")
    parser.add_argument("--out-rag-jsonl", default="data/rag/kg_rag_chunks.jsonl")
    parser.add_argument("--out-report-md", default="docs/rag/kg_rag_usage.md")
    parser.add_argument("--max-mutex-rules", type=int, default=2500, help="Cap mutex rules added to KG/RAG to keep first-pass graph compact.")
    return parser.parse_args()


def clean(value: Any) -> str:
    return "" if value is None else str(value).strip()


def norm(value: Any) -> str:
    return clean(value).upper()


def hash_id(prefix: str, *parts: Any) -> str:
    text = "|".join(clean(part) for part in parts)
    return f"{prefix}:{hashlib.sha1(text.encode('utf-8', 'ignore')).hexdigest()[:16]}"


def props(**kwargs: Any) -> str:
    return json.dumps({k: v for k, v in kwargs.items() if clean(v) != ""}, ensure_ascii=False, sort_keys=True)


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def node(node_id: str, label: str, name: str, typ: str, source: str, properties: str = "{}", ontology_id: str = "not_reported") -> dict[str, Any]:
    return {
        "node_id": node_id,
        "node_label": label,
        "name": name,
        "normalized_name": name.upper(),
        "type": typ,
        "source": source,
        "ontology_id": ontology_id,
        "properties_json": properties,
    }


def rel(
    start: str,
    end: str,
    rtype: str,
    confidence: str,
    strength: Any,
    evidence_ids: str,
    pmids: str,
    dois: str,
    inclusion: str,
    audit: str,
    properties: str = "{}",
) -> dict[str, Any]:
    return {
        "relationship_id": hash_id("REL", start, end, rtype, evidence_ids, properties),
        "start_node_id": start,
        "end_node_id": end,
        "relationship_type": rtype,
        "confidence_score": confidence,
        "evidence_strength_score": strength,
        "evidence_ids": evidence_ids,
        "source_pmids": pmids,
        "source_dois": dois,
        "graph_inclusion_status": inclusion,
        "audit_priority": audit,
        "properties_json": properties,
    }


def rule_confidence(rule_status: str, p_value: str = "") -> tuple[str, int, str]:
    if rule_status == "high_confidence_rule":
        return "High", 5, "no_review_needed"
    try:
        p = float(p_value)
    except ValueError:
        p = 1.0
    if p < 0.05:
        return "Medium", 4, "medium_review"
    return "Low", 2, "high_review"


def add_unique(nodes: dict[str, dict[str, Any]], row: dict[str, Any]) -> None:
    nodes.setdefault(row["node_id"], row)


def make_rag_chunk(chunk_id: str, title: str, text: str, metadata: dict[str, Any]) -> dict[str, Any]:
    return {"chunk_id": chunk_id, "title": title, "text": text, "metadata": metadata}


def feature_node_id(feature_id: str, feature_type: str, fallback_label: str) -> str:
    if feature_type == "gene":
        return f"Gene:{norm(fallback_label or feature_id.replace('GENE:', ''))}"
    safe = clean(feature_id).replace("LOCUS:", "")
    return f"LocusFeature:{safe}" if safe else f"LocusFeature:{norm(fallback_label)}"


def add_feature_node(nodes: dict[str, dict[str, Any]], feature_id: str, feature_type: str, feature_label: str) -> str:
    label = clean(feature_label) or clean(feature_id)
    node_id = feature_node_id(feature_id, feature_type, label)
    if feature_type == "gene":
        add_unique(nodes, node(node_id, "Gene", norm(label), "gene", "supplementary_cohort_rules"))
    else:
        add_unique(
            nodes,
            node(
                node_id,
                "LocusFeature",
                label,
                "locus_feature",
                "supplementary_cohort_rules",
                props(feature_id=feature_id, feature_type=feature_type),
            ),
        )
    return node_id


def build_extension(baselines: list[dict[str, str]], co_rules: list[dict[str, str]], mutex_rules: list[dict[str, str]], max_mutex: int) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    nodes: dict[str, dict[str, Any]] = {}
    rels: list[dict[str, Any]] = []
    chunks: list[dict[str, Any]] = []

    source_node_id = "CohortLayer:supplementary_material_cohort"
    add_unique(nodes, node(source_node_id, "SupplementaryCohortLayer", "supplementary_material_cohort", "cohort_rule_layer", "supplementary_cohort_rules"))

    for row in baselines:
        cancer = clean(row.get("cancer_type")) or "not_reported"
        feature_id = clean(row.get("feature_id")) or (f"GENE:{norm(row.get('gene'))}" if clean(row.get("gene")) else "")
        feature_type = clean(row.get("feature_type")) or ("gene" if clean(row.get("gene")) else "feature")
        feature_label = clean(row.get("feature_label")) or clean(row.get("gene")) or feature_id
        if not feature_id:
            continue
        cancer_id = f"Cancer:{cancer.lower()}"
        feature_node = add_feature_node(nodes, feature_id, feature_type, feature_label)
        baseline_id = f"CancerGeneBaseline:{row['baseline_id']}"
        paper_id = clean(row.get("supporting_paper_ids"))
        paper_node_id = f"Paper:{paper_id}" if paper_id and ";" not in paper_id else ""
        add_unique(nodes, node(cancer_id, "Cancer", cancer, "cancer", "supplementary_cohort_rules"))
        add_unique(
            nodes,
            node(
                baseline_id,
                "CancerGeneBaseline",
                f"{cancer}:{feature_label}",
                "cancer_gene_baseline",
                "supplementary_cohort_rules",
                props(
                    feature_id=feature_id,
                    feature_type=feature_type,
                    feature_label=feature_label,
                    sample_count=row.get("sample_count"),
                    gene_sample_count=row.get("gene_sample_count"),
                    gene_frequency=row.get("gene_frequency"),
                    rule_status=row.get("rule_status"),
                    supporting_paper_ids=row.get("supporting_paper_ids"),
                    supporting_journals=row.get("supporting_journals"),
                    supporting_titles=row.get("supporting_titles"),
                    supporting_dois=row.get("supporting_dois"),
                    supporting_pmids=row.get("supporting_pmids"),
                    supporting_source_document_ids=row.get("supporting_source_document_ids"),
                ),
            ),
        )
        rels.append(rel(cancer_id, baseline_id, "HAS_GENE_BASELINE", "Medium", 3, row.get("supporting_event_ids", ""), row.get("supporting_pmids", ""), row.get("supporting_dois", ""), "include_as_cohort_rule", "medium_review"))
        rels.append(rel(baseline_id, feature_node, "BASELINE_FOR_FEATURE", "Medium", 3, row.get("supporting_event_ids", ""), row.get("supporting_pmids", ""), row.get("supporting_dois", ""), "include_as_cohort_rule", "medium_review"))
        if feature_type == "gene":
            rels.append(rel(baseline_id, feature_node, "BASELINE_FOR_GENE", "Medium", 3, row.get("supporting_event_ids", ""), row.get("supporting_pmids", ""), row.get("supporting_dois", ""), "include_as_cohort_rule", "medium_review"))
        rels.append(rel(baseline_id, source_node_id, "DERIVED_FROM_COHORT_LAYER", "Medium", 3, row.get("supporting_event_ids", ""), row.get("supporting_pmids", ""), row.get("supporting_dois", ""), "include_as_cohort_rule", "medium_review"))
        if paper_node_id:
            add_unique(nodes, node(paper_node_id, "Paper", paper_id, "paper", "supplementary_cohort_rules", props(title=row.get("supporting_titles"), journal=row.get("supporting_journals"), doi=row.get("supporting_dois"), pmid=row.get("supporting_pmids"))))
            rels.append(rel(baseline_id, paper_node_id, "SUPPORTED_BY_PAPER", "Medium", 3, row.get("supporting_event_ids", ""), row.get("supporting_pmids", ""), row.get("supporting_dois", ""), "include_as_cohort_rule", "medium_review"))

        text = (
            f"In {cancer}, {feature_label} ({feature_type}) appears in {row.get('gene_sample_count')} of {row.get('sample_count')} supplementary-material samples "
            f"(frequency {row.get('gene_frequency')}). Source: {row.get('supporting_journals')} | {row.get('supporting_titles')} | "
            f"DOI {row.get('supporting_dois')} | PMID {row.get('supporting_pmids')}."
        )
        chunks.append(make_rag_chunk(baseline_id, f"{cancer} {feature_label} baseline", text, {"kind": "cancer_gene_baseline", **row}))

    all_rules = list(co_rules) + list(mutex_rules[:max_mutex])
    for row in all_rules:
        cancer = clean(row.get("cancer_type")) or "not_reported"
        feature_a = clean(row.get("feature_a_id")) or (f"GENE:{norm(row.get('gene_a'))}" if clean(row.get("gene_a")) else "")
        feature_b = clean(row.get("feature_b_id")) or (f"GENE:{norm(row.get('gene_b'))}" if clean(row.get("gene_b")) else "")
        feature_a_type = clean(row.get("feature_a_type")) or ("gene" if clean(row.get("gene_a")) else "feature")
        feature_b_type = clean(row.get("feature_b_type")) or ("gene" if clean(row.get("gene_b")) else "feature")
        label_a = clean(row.get("feature_a_label")) or clean(row.get("gene_a")) or feature_a
        label_b = clean(row.get("feature_b_label")) or clean(row.get("gene_b")) or feature_b
        if not (feature_a and feature_b):
            continue
        rule_id = f"GeneAssociationRule:{row['rule_id']}"
        cancer_id = f"Cancer:{cancer.lower()}"
        feature_a_node = add_feature_node(nodes, feature_a, feature_a_type, label_a)
        feature_b_node = add_feature_node(nodes, feature_b, feature_b_type, label_b)
        rtype = clean(row.get("rule_type"))
        rel_type = "FEATURE_PAIR_CO_OCCURS_IN_CANCER" if "cooccurrence" in rtype else "FEATURE_PAIR_MUTUALLY_EXCLUSIVE_IN_CANCER"
        if feature_a_type == "gene" and feature_b_type == "gene":
            rel_type = "GENE_PAIR_CO_OCCURS_IN_CANCER" if "cooccurrence" in rtype else "GENE_PAIR_MUTUALLY_EXCLUSIVE_IN_CANCER"
        confidence, strength, audit = rule_confidence(clean(row.get("rule_status")), clean(row.get("fisher_p_value")))
        add_unique(nodes, node(cancer_id, "Cancer", cancer, "cancer", "supplementary_cohort_rules"))
        add_unique(
            nodes,
            node(
                rule_id,
                "GeneAssociationRule",
                f"{cancer}:{label_a}-{label_b}:{rtype}",
                "gene_association_rule",
                "supplementary_cohort_rules",
                props(**row),
            ),
        )
        rels.append(rel(cancer_id, rule_id, "HAS_GENE_ASSOCIATION_RULE", confidence, strength, row.get("supporting_event_ids", ""), row.get("supporting_pmids", ""), row.get("supporting_dois", ""), "include_as_cohort_rule", audit))
        rels.append(rel(rule_id, feature_a_node, "RULE_INVOLVES_FEATURE", confidence, strength, row.get("supporting_event_ids", ""), row.get("supporting_pmids", ""), row.get("supporting_dois", ""), "include_as_cohort_rule", audit, props(role="feature_a", feature_type=feature_a_type)))
        rels.append(rel(rule_id, feature_b_node, "RULE_INVOLVES_FEATURE", confidence, strength, row.get("supporting_event_ids", ""), row.get("supporting_pmids", ""), row.get("supporting_dois", ""), "include_as_cohort_rule", audit, props(role="feature_b", feature_type=feature_b_type)))
        if feature_a_type == "gene":
            rels.append(rel(rule_id, feature_a_node, "RULE_INVOLVES_GENE", confidence, strength, row.get("supporting_event_ids", ""), row.get("supporting_pmids", ""), row.get("supporting_dois", ""), "include_as_cohort_rule", audit, props(role="gene_a")))
        if feature_b_type == "gene":
            rels.append(rel(rule_id, feature_b_node, "RULE_INVOLVES_GENE", confidence, strength, row.get("supporting_event_ids", ""), row.get("supporting_pmids", ""), row.get("supporting_dois", ""), "include_as_cohort_rule", audit, props(role="gene_b")))
        rels.append(rel(feature_a_node, feature_b_node, rel_type, confidence, strength, row.get("supporting_event_ids", ""), row.get("supporting_pmids", ""), row.get("supporting_dois", ""), "include_as_cohort_rule", audit, props(rule_id=row.get("rule_id"), cancer_type=cancer, pair_count=row.get("pair_count"), lift=row.get("lift"), odds_ratio=row.get("odds_ratio"), fisher_p_value=row.get("fisher_p_value"))))
        rels.append(rel(rule_id, source_node_id, "DERIVED_FROM_COHORT_LAYER", confidence, strength, row.get("supporting_event_ids", ""), row.get("supporting_pmids", ""), row.get("supporting_dois", ""), "include_as_cohort_rule", audit))

        association_word = "co-occurs with" if "cooccurrence" in rtype else "is a mutual-exclusion candidate with"
        text = (
            f"In {cancer}, {label_a} ({feature_a_type}) {association_word} {label_b} ({feature_b_type}). Pair count {row.get('pair_count')} among {row.get('sample_count')} samples; "
            f"lift {row.get('lift')}, odds ratio {row.get('odds_ratio')}, Fisher p={row.get('fisher_p_value')}. "
            f"Rule status: {row.get('rule_status')}. Source: {row.get('supporting_journals')} | {row.get('supporting_titles')} | "
            f"DOI {row.get('supporting_dois')} | PMID {row.get('supporting_pmids')}."
        )
        chunks.append(make_rag_chunk(rule_id, f"{cancer} {label_a}-{label_b} {rtype}", text, {"kind": "gene_association_rule", **row}))

    return list(nodes.values()), rels, chunks


def append_csv(base: Path, extension: Path, out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    if base.exists():
        shutil.copyfile(base, out)
        mode = "a"
        write_header = False
    else:
        mode = "w"
        write_header = True
    with out.open(mode, encoding="utf-8", newline="") as out_handle, extension.open("r", encoding="utf-8", newline="") as ext_handle:
        reader = csv.reader(ext_handle)
        writer = csv.writer(out_handle)
        header = next(reader, None)
        if write_header and header:
            writer.writerow(header)
        for row in reader:
            writer.writerow(row)


def write_jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n")


def write_report(path: Path, stats: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = f"""# KG-RAG Usage

## Generated Files

- Extension nodes: `data/graph/cohort_rule_nodes.csv`
- Extension relationships: `data/graph/cohort_rule_relationships.csv`
- Augmented nodes: `data/graph/nodes_augmented.csv`
- Augmented relationships: `data/graph/relationships_augmented.csv`
- GraphRAG chunks: `data/rag/kg_rag_chunks.jsonl`

## Counts

- Cohort/rule nodes: {stats['nodes']}
- Cohort/rule relationships: {stats['relationships']}
- RAG chunks: {stats['chunks']}

## How To Use As RAG

1. Retrieve relevant chunks by cancer type, gene, variant, or question terms.
2. Pull the supporting node/edge metadata from the chunk metadata.
3. Answer with traceability: paper, journal, DOI/PMID, source document, event IDs.
4. Treat `candidate_rule` as review evidence, not automatic clinical truth.

Example:

```powershell
python scripts\\mrdskill_query_mrd_kg_rag.py --query "Bladder ASXL1 TET2 cooccur"
```

"""
    path.write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    baselines = read_csv(artifact_path(args.baseline_csv))
    co_rules = read_csv(artifact_path(args.cooccurrence_csv))
    mutex_rules = read_csv(artifact_path(args.mutex_csv))
    nodes, rels, chunks = build_extension(baselines, co_rules, mutex_rules, args.max_mutex_rules)
    write_csv(artifact_path(args.out_nodes_extension), nodes, NODE_FIELDS)
    write_csv(artifact_path(args.out_relationships_extension), rels, REL_FIELDS)
    append_csv(artifact_path(args.base_nodes), artifact_path(args.out_nodes_extension), artifact_path(args.out_nodes_augmented))
    append_csv(artifact_path(args.base_relationships), artifact_path(args.out_relationships_extension), artifact_path(args.out_relationships_augmented))
    write_jsonl(artifact_path(args.out_rag_jsonl), chunks)
    stats = {"nodes": len(nodes), "relationships": len(rels), "chunks": len(chunks)}
    write_report(artifact_path(args.out_report_md), stats)
    print(json.dumps(stats, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
