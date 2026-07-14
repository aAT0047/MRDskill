#!/usr/bin/env python3
"""Build graph-ready CSV files from scored MRD evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


NODE_COLUMNS = [
    "node_id",
    "node_label",
    "name",
    "normalized_name",
    "type",
    "source",
    "ontology_id",
    "properties_json",
]

REL_COLUMNS = [
    "relationship_id",
    "start_node_id",
    "end_node_id",
    "relationship_type",
    "confidence_score",
    "evidence_strength_score",
    "evidence_ids",
    "source_pmids",
    "source_dois",
    "graph_inclusion_status",
    "audit_priority",
    "properties_json",
]

REVIEW_COLUMNS = [
    "review_task_id",
    "evidence_id",
    "paper_id",
    "audit_priority",
    "reason",
    "title",
    "doi",
    "pmid",
    "source_traceability_status",
    "biological_plausibility_status",
    "graph_inclusion_status",
]

ALLOWED_INCLUSION = {
    "auto_include_high_confidence",
    "include_as_medium_confidence",
    "include_as_low_confidence_annotation",
    "candidate_only",
    "needs_human_review",
}


def stable_id(prefix: str, *parts: str) -> str:
    raw = "|".join(part.strip().lower() for part in parts if part is not None)
    digest = hashlib.sha1(raw.encode("utf-8")).hexdigest()[:14]
    return f"{prefix}:{digest}"


def clean(value: str | None) -> str:
    value = (value or "").strip()
    return value if value else "not_reported"


def is_reported(value: str | None) -> bool:
    return clean(value) not in {"not_reported", "unavailable", "uncertain", "not_assessable"}


def split_values(value: str | None) -> list[str]:
    values = []
    for item in re.split(r"[;|]", clean(value)):
        item = item.strip()
        if is_reported(item):
            values.append(item)
    return sorted(set(values))


def add_node(
    nodes: dict[str, dict[str, str]],
    node_id: str,
    label: str,
    name: str,
    normalized_name: str | None = None,
    node_type: str = "entity",
    source: str = "mrd_evidence_scored",
    ontology_id: str = "not_reported",
    properties: dict[str, object] | None = None,
) -> str:
    if node_id not in nodes:
        nodes[node_id] = {
            "node_id": node_id,
            "node_label": label,
            "name": name,
            "normalized_name": normalized_name or name,
            "type": node_type,
            "source": source,
            "ontology_id": ontology_id,
            "properties_json": json.dumps(properties or {}, ensure_ascii=False, sort_keys=True),
        }
    return node_id


def add_relationship(
    relationships: dict[str, dict[str, str]],
    start: str,
    end: str,
    rel_type: str,
    row: dict[str, str],
    properties: dict[str, object] | None = None,
) -> str:
    evidence_id = clean(row.get("evidence_id"))
    rel_id = stable_id("REL", start, rel_type, end, evidence_id)
    relationships[rel_id] = {
        "relationship_id": rel_id,
        "start_node_id": start,
        "end_node_id": end,
        "relationship_type": rel_type,
        "confidence_score": clean(row.get("confidence_score")),
        "evidence_strength_score": clean(row.get("evidence_strength_score")),
        "evidence_ids": evidence_id,
        "source_pmids": clean(row.get("pmid")),
        "source_dois": clean(row.get("doi")),
        "graph_inclusion_status": clean(row.get("graph_inclusion_status")),
        "audit_priority": clean(row.get("audit_priority")),
        "properties_json": json.dumps(properties or {}, ensure_ascii=False, sort_keys=True),
    }
    return rel_id


def clinical_outcome_nodes(row: dict[str, str]) -> list[tuple[str, str]]:
    outcomes: list[tuple[str, str]] = []
    if row.get("recurrence_association") == "reported" or row.get("mrd_relevance") == "recurrence_marker":
        outcomes.append(("recurrence", "recurrence"))
    if row.get("relapse_association") == "reported" or row.get("mrd_relevance") == "relapse_driver":
        outcomes.append(("relapse", "relapse"))
    if row.get("prognosis_association") == "reported" or is_reported(row.get("survival_endpoint")):
        outcomes.append(("prognosis", clean(row.get("survival_endpoint"))))
    if row.get("treatment_monitoring_value") == "reported":
        outcomes.append(("treatment_response", "treatment_response"))
    return sorted(set(outcomes))


def build_graph(rows: list[dict[str, str]]) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    nodes: dict[str, dict[str, str]] = {}
    relationships: dict[str, dict[str, str]] = {}
    review_tasks: list[dict[str, str]] = []
    evidence_summaries: list[dict[str, str]] = []

    for row in rows:
        inclusion = clean(row.get("graph_inclusion_status"))
        if inclusion == "exclude":
            continue
        if inclusion not in ALLOWED_INCLUSION:
            inclusion = "candidate_only"

        evidence_id = clean(row.get("evidence_id"))
        evidence_node = add_node(
            nodes,
            f"Evidence:{evidence_id}",
            "Evidence",
            evidence_id,
            node_type="evidence",
            properties={
                "confidence_score": clean(row.get("confidence_score")),
                "evidence_source_type_final": clean(row.get("evidence_source_type_final")),
                "graph_inclusion_status": inclusion,
                "audit_priority": clean(row.get("audit_priority")),
                "source_traceability_status": clean(row.get("source_traceability_status")),
                "confidence_reason": clean(row.get("confidence_reason")),
            },
        )

        source_document_node = None
        if is_reported(row.get("source_document_id")):
            source_document_node = add_node(
                nodes,
                f"SourceDocument:{clean(row.get('source_document_id'))}",
                "SourceDocument",
                clean(row.get("source_file")),
                node_type="evidence_source",
                properties={
                    "source_document_type": clean(row.get("source_document_type")),
                    "source_path": clean(row.get("source_path")),
                    "source_sheet": clean(row.get("source_sheet")),
                    "source_row": clean(row.get("source_row")),
                    "source_column": clean(row.get("source_column")),
                    "source_parser": clean(row.get("source_parser")),
                    "source_parser_status": clean(row.get("source_parser_status")),
                },
            )
            add_relationship(relationships, evidence_node, source_document_node, "FROM_SOURCE_DOCUMENT", row)

        paper_id = clean(row.get("paper_id"))
        paper_node = add_node(
            nodes,
            f"Paper:{paper_id}",
            "Paper",
            clean(row.get("title")),
            node_type="evidence_source",
            properties={
                "paper_id": paper_id,
                "year": clean(row.get("year")),
                "journal": clean(row.get("journal")),
                "doi": clean(row.get("doi")),
                "pmid": clean(row.get("pmid")),
                "pmcid": clean(row.get("pmcid")),
            },
        )
        add_relationship(relationships, evidence_node, paper_node, "FROM_PAPER", row)
        if source_document_node:
            add_relationship(relationships, source_document_node, paper_node, "BELONGS_TO_PAPER", row)

        source_location_node = add_node(
            nodes,
            stable_id(
                "SourceLocation",
                evidence_id,
                clean(row.get("source_section")),
                clean(row.get("source_table")),
                clean(row.get("source_figure")),
                clean(row.get("supplementary_table")),
            ),
            "SourceLocation",
            clean(row.get("source_section")),
            node_type="evidence_source",
            properties={
                "source_section": clean(row.get("source_section")),
                "source_table": clean(row.get("source_table")),
                "source_figure": clean(row.get("source_figure")),
                "supplementary_file": clean(row.get("supplementary_file")),
                "supplementary_table": clean(row.get("supplementary_table")),
                "quoted_or_paraphrased_evidence": clean(row.get("quoted_or_paraphrased_evidence")),
            },
        )
        add_relationship(relationships, evidence_node, source_location_node, "LOCATED_IN", row)

        cancer_node = None
        if is_reported(row.get("cancer_type")):
            cancer_name = clean(row.get("cancer_type"))
            cancer_node = add_node(
                nodes,
                stable_id("Cancer", cancer_name),
                "Cancer",
                cancer_name,
                node_type="biomedical_entity",
            )
            add_relationship(relationships, evidence_node, cancer_node, "SUPPORTS", row)

        assay_node = None
        if is_reported(row.get("assay_type")):
            assay_name = clean(row.get("assay_type"))
            assay_node = add_node(
                nodes,
                stable_id("Assay", assay_name, clean(row.get("tumor_informed_or_tumor_naive"))),
                "Assay",
                assay_name,
                node_type="biomedical_entity",
                properties={
                    "tumor_informed_or_tumor_naive": clean(row.get("tumor_informed_or_tumor_naive")),
                    "sample_type": clean(row.get("sample_type")),
                    "limit_of_detection": clean(row.get("limit_of_detection")),
                    "VAF_cutoff": clean(row.get("VAF_cutoff")),
                    "matched_WBC_used": clean(row.get("matched_WBC_used")),
                },
            )
            add_relationship(relationships, evidence_node, assay_node, "SUPPORTS", row)

        mrd_signal_node = None
        if is_reported(row.get("mrd_relevance")):
            signal_name = clean(row.get("mrd_relevance"))
            signal_key = "|".join(
                [
                    clean(row.get("cancer_type")),
                    clean(row.get("gene_symbol_normalized")),
                    clean(row.get("variant")),
                    signal_name,
                    clean(row.get("assay_type")),
                ]
            )
            mrd_signal_node = add_node(
                nodes,
                stable_id("MRDSignal", signal_key),
                "MRDSignal",
                signal_name,
                node_type="qc_interpretation",
                properties={
                    "biomarker_type": clean(row.get("biomarker_type")),
                    "biological_plausibility_status": clean(row.get("biological_plausibility_status")),
                    "clinical_relevance_status": clean(row.get("clinical_relevance_status")),
                },
            )
            add_relationship(relationships, evidence_node, mrd_signal_node, "SUPPORTS", row)
            if cancer_node:
                add_relationship(relationships, mrd_signal_node, cancer_node, "OBSERVED_IN", row)
            if assay_node:
                add_relationship(relationships, mrd_signal_node, assay_node, "DETECTED_BY", row)

        for gene in split_values(row.get("gene_symbol_normalized")):
            gene_node = add_node(
                nodes,
                stable_id("Gene", gene),
                "Gene",
                gene,
                node_type="biomedical_entity",
                properties={"known_CHIP_gene": gene in {"DNMT3A", "TET2", "ASXL1", "JAK2", "TP53", "PPM1D"}},
            )
            add_relationship(relationships, evidence_node, gene_node, "SUPPORTS", row)
            if cancer_node:
                add_relationship(relationships, cancer_node, gene_node, "ASSOCIATED_WITH", row)
            if assay_node:
                add_relationship(relationships, assay_node, gene_node, "DETECTS", row)
            if mrd_signal_node:
                add_relationship(relationships, gene_node, mrd_signal_node, "HAS_MRD_RELEVANCE", row)

            if is_reported(row.get("variant")):
                variant_name = clean(row.get("variant"))
                variant_node = add_node(
                    nodes,
                    stable_id("Variant", gene, variant_name),
                    "Variant",
                    variant_name,
                    node_type="biomedical_entity",
                    properties={
                        "gene_symbol": gene,
                        "variant_type": clean(row.get("variant_type")),
                        "mutation_hotspot": clean(row.get("mutation_hotspot")),
                        "genomic_coordinate": clean(row.get("genomic_coordinate")),
                    },
                )
                add_relationship(relationships, gene_node, variant_node, "HAS_VARIANT", row)
                add_relationship(relationships, evidence_node, variant_node, "SUPPORTS", row)
                if cancer_node:
                    add_relationship(relationships, variant_node, cancer_node, "OBSERVED_IN", row)
                if assay_node:
                    add_relationship(relationships, assay_node, variant_node, "DETECTS_VARIANT", row)
                if mrd_signal_node:
                    add_relationship(relationships, variant_node, mrd_signal_node, "HAS_MRD_RELEVANCE", row)

        if mrd_signal_node:
            for outcome_type, outcome_name in clinical_outcome_nodes(row):
                outcome_node = add_node(
                    nodes,
                    stable_id("ClinicalOutcome", outcome_type, outcome_name),
                    "ClinicalOutcome",
                    outcome_name,
                    node_type="biomedical_entity",
                    properties={"outcome_type": outcome_type},
                )
                rel_type = {
                    "recurrence": "ASSOCIATED_WITH_RECURRENCE",
                    "relapse": "ASSOCIATED_WITH_RELAPSE",
                    "prognosis": "ASSOCIATED_WITH_PROGNOSIS",
                    "treatment_response": "USED_FOR_TREATMENT_MONITORING",
                }.get(outcome_type, "INDICATES")
                add_relationship(relationships, mrd_signal_node, outcome_node, rel_type, row)

        if clean(row.get("audit_priority")) in {"urgent_review", "high_review", "medium_review"} or inclusion in {
            "candidate_only",
            "needs_human_review",
        }:
            review_tasks.append(
                {
                    "review_task_id": stable_id("ReviewTask", evidence_id),
                    "evidence_id": evidence_id,
                    "paper_id": paper_id,
                    "audit_priority": clean(row.get("audit_priority")),
                    "reason": clean(row.get("downgrade_reason")),
                    "title": clean(row.get("title")),
                    "doi": clean(row.get("doi")),
                    "pmid": clean(row.get("pmid")),
                    "source_traceability_status": clean(row.get("source_traceability_status")),
                    "biological_plausibility_status": clean(row.get("biological_plausibility_status")),
                    "graph_inclusion_status": inclusion,
                }
            )

        evidence_summaries.append(
            {
                "evidence_id": evidence_id,
                "paper_id": paper_id,
                "title": clean(row.get("title")),
                "summary": clean(row.get("quoted_or_paraphrased_evidence")),
                "confidence_score": clean(row.get("confidence_score")),
                "graph_inclusion_status": inclusion,
                "audit_priority": clean(row.get("audit_priority")),
                "doi": clean(row.get("doi")),
                "pmid": clean(row.get("pmid")),
            }
        )

    return list(nodes.values()), list(relationships.values()), review_tasks, evidence_summaries


def add_source_documents(
    nodes: list[dict[str, str]],
    relationships: list[dict[str, str]],
    sources: list[dict[str, str]],
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    node_map = {node["node_id"]: node for node in nodes}
    rel_map = {rel["relationship_id"]: rel for rel in relationships}
    for source in sources:
        source_document_id = clean(source.get("source_document_id"))
        if not is_reported(source_document_id):
            continue
        source_node_id = f"SourceDocument:{source_document_id}"
        if source_node_id not in node_map:
            add_node(
                node_map,
                source_node_id,
                "SourceDocument",
                clean(source.get("source_file")),
                node_type="evidence_source",
                properties={
                    "paper_id": clean(source.get("paper_id")),
                    "source_document_type": clean(source.get("source_document_type")),
                    "source_path": clean(source.get("source_path")),
                    "extension": clean(source.get("extension")),
                    "size_bytes": clean(source.get("size_bytes")),
                    "parser": clean(source.get("parser")),
                    "parser_status": clean(source.get("parser_status")),
                    "skip_reason": clean(source.get("skip_reason")),
                    "evidence_count": clean(source.get("evidence_count")),
                },
            )
        paper_id = clean(source.get("paper_id"))
        if is_reported(paper_id):
            paper_node_id = f"Paper:{paper_id}"
            if paper_node_id not in node_map:
                add_node(
                    node_map,
                    paper_node_id,
                    "Paper",
                    paper_id,
                    node_type="evidence_source",
                    properties={"paper_id": paper_id},
                )
            rel_id = stable_id("REL", source_node_id, "BELONGS_TO_PAPER", paper_node_id, "source_manifest")
            if rel_id not in rel_map:
                rel_map[rel_id] = {
                    "relationship_id": rel_id,
                    "start_node_id": source_node_id,
                    "end_node_id": paper_node_id,
                    "relationship_type": "BELONGS_TO_PAPER",
                    "confidence_score": "not_applicable",
                    "evidence_strength_score": "not_applicable",
                    "evidence_ids": "not_applicable",
                    "source_pmids": "not_reported",
                    "source_dois": "not_reported",
                    "graph_inclusion_status": "source_manifest",
                    "audit_priority": "not_applicable",
                    "properties_json": json.dumps(
                        {
                            "parser_status": clean(source.get("parser_status")),
                            "evidence_count": clean(source.get("evidence_count")),
                        },
                        ensure_ascii=False,
                        sort_keys=True,
                    ),
                }
    return list(node_map.values()), list(rel_map.values())


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def write_csv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_jsonl(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description="Build graph-ready MRD CSV files.")
    parser.add_argument("--input-csv", default="data/evidence/mrd_evidence_scored.csv")
    parser.add_argument("--nodes-csv", default="data/graph/nodes.csv")
    parser.add_argument("--relationships-csv", default="data/graph/relationships.csv")
    parser.add_argument("--review-tasks-csv", default="data/graph/review_tasks.csv")
    parser.add_argument("--evidence-summaries-jsonl", default="data/rag/evidence_summaries.jsonl")
    parser.add_argument("--sources-csv", default="", help="Optional full source-document manifest CSV.")
    args = parser.parse_args()

    rows = read_csv(Path(args.input_csv))
    nodes, relationships, review_tasks, summaries = build_graph(rows)
    if args.sources_csv:
        source_path = Path(args.sources_csv)
        if source_path.exists():
            sources = read_csv(source_path)
            nodes, relationships = add_source_documents(nodes, relationships, sources)
    write_csv(Path(args.nodes_csv), NODE_COLUMNS, nodes)
    write_csv(Path(args.relationships_csv), REL_COLUMNS, relationships)
    write_csv(Path(args.review_tasks_csv), REVIEW_COLUMNS, review_tasks)
    write_jsonl(Path(args.evidence_summaries_jsonl), summaries)
    print(f"Wrote nodes: {Path(args.nodes_csv).resolve()} ({len(nodes)} rows)")
    print(f"Wrote relationships: {Path(args.relationships_csv).resolve()} ({len(relationships)} rows)")
    print(f"Wrote review tasks: {Path(args.review_tasks_csv).resolve()} ({len(review_tasks)} rows)")
    print(f"Wrote evidence summaries: {Path(args.evidence_summaries_jsonl).resolve()} ({len(summaries)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
