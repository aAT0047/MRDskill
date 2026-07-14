#!/usr/bin/env python3
"""Build a browser-friendly visualization for the full MRD support KG."""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from mrdskill_project_paths import artifact_path, artifact_root


MISSING = {"", "not_reported", "unavailable", "uncertain", "not_assessable", "none", "nan"}

CONFIDENCE_ORDER = {"High": 3, "Medium": 2, "Low": 1, "Exclude": 0}

NODE_COLORS = {
    "Cancer": "#2f6f73",
    "Gene": "#6d5bd0",
    "Variant": "#9a6a2f",
    "Assay": "#2d6cdf",
    "MRDSignal": "#b33f62",
    "ClinicalOutcome": "#3f7d3f",
    "GeneAssociationRule": "#c23b7a",
    "CancerGeneBaseline": "#138a72",
    "LocusFeature": "#8b5a2b",
}

EDGE_COLORS = {
    "High": "#25905a",
    "Medium": "#b8860b",
    "Low": "#b94a48",
    "Exclude": "#8a8f98",
}

EDGE_TYPE_COLORS = {
    "GENE_COOCCURS_WITH": "#188a54",
    "GENE_MUTUALLY_EXCLUSIVE_WITH": "#c23b7a",
    "CANCER_HAS_BASELINE": "#138a72",
    "RULE_IN_CANCER": "#9b4d96",
    "RULE_INVOLVES_GENE": "#6d5bd0",
}

DEFAULT_EDGE_LIMITS = {
    "ASSOCIATED_WITH": 260,
    "HAS_VARIANT": 260,
    "OBSERVED_IN": 240,
    "HAS_MRD_RELEVANCE": 420,
    "HAS_MRD_RELEVANCE_VARIANT": 180,
    "ASSOCIATED_WITH_RECURRENCE": 120,
    "ASSOCIATED_WITH_RELAPSE": 120,
    "ASSOCIATED_WITH_PROGNOSIS": 160,
    "USED_FOR_TREATMENT_MONITORING": 120,
    "DETECTS": 220,
    "DETECTS_VARIANT": 160,
    "DETECTED_BY": 220,
    "CANCER_HAS_BASELINE": 180,
    "RULE_IN_CANCER": 120,
    "RULE_INVOLVES_GENE": 220,
    "GENE_COOCCURS_WITH": 320,
    "GENE_MUTUALLY_EXCLUSIVE_WITH": 220,
}


def clean(value: str | None) -> str:
    value = (value or "").strip()
    return value if value else "not_reported"


def is_reported(value: str | None) -> bool:
    return clean(value).lower() not in MISSING


def split_terms(value: str | None, max_items: int = 8) -> list[str]:
    if not is_reported(value):
        return []
    value = clean(value)
    parts = [item.strip() for item in re.split(r"[;|]", value) if item.strip()]
    seen: set[str] = set()
    terms: list[str] = []
    for part in parts:
        if part.lower() in MISSING:
            continue
        if part not in seen:
            seen.add(part)
            terms.append(part)
        if len(terms) >= max_items:
            break
    return terms


def clipped(value: str | None, limit: int = 160) -> str:
    value = clean(value)
    return value if len(value) <= limit else value[: limit - 3] + "..."


def csv_row_count(path: Path) -> int:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return sum(1 for _ in csv.DictReader(handle))


def jsonl_count(path: Path) -> int:
    with path.open("r", encoding="utf-8") as handle:
        return sum(1 for _ in handle)


def node_key(label: str, name: str) -> str:
    normalized = re.sub(r"\s+", " ", clean(name)).strip().lower()
    return f"{label}:{normalized}"


def relation_label(edge_type: str) -> str:
    labels = {
        "ASSOCIATED_WITH": "Cancer-Gene",
        "HAS_VARIANT": "Gene-Variant",
        "OBSERVED_IN": "Variant-Cancer",
        "HAS_MRD_RELEVANCE": "MRD relevance",
        "HAS_MRD_RELEVANCE_VARIANT": "Variant MRD",
        "ASSOCIATED_WITH_RECURRENCE": "Recurrence",
        "ASSOCIATED_WITH_RELAPSE": "Relapse",
        "ASSOCIATED_WITH_PROGNOSIS": "Prognosis",
        "USED_FOR_TREATMENT_MONITORING": "Monitoring",
        "DETECTS": "Assay-Gene",
        "DETECTS_VARIANT": "Assay-Variant",
        "DETECTED_BY": "MRD-Assay",
        "CANCER_HAS_BASELINE": "Cancer baseline",
        "RULE_IN_CANCER": "Rule cancer",
        "RULE_INVOLVES_GENE": "Rule gene",
        "GENE_COOCCURS_WITH": "Gene co-occurrence",
        "GENE_MUTUALLY_EXCLUSIVE_WITH": "Gene mutual exclusion",
    }
    return labels.get(edge_type, edge_type.replace("_", " ").title())


def dominant_confidence(counter: Counter[str]) -> str:
    if not counter:
        return "Low"
    return max(counter, key=lambda key: (CONFIDENCE_ORDER.get(key, -1), counter[key]))


def edge_score(edge: dict[str, Any]) -> float:
    counts: Counter[str] = edge["confidence_counts"]
    audit: Counter[str] = edge["audit_counts"]
    return (
        edge["count"]
        + 2.2 * counts.get("High", 0)
        + 1.2 * counts.get("Medium", 0)
        + 0.4 * counts.get("Low", 0)
        + 0.35 * audit.get("urgent_review", 0)
    )


def safe_int(value: str | None) -> int:
    try:
        return int(float(clean(value)))
    except ValueError:
        return 0


def add_node(
    nodes: dict[str, dict[str, Any]],
    label: str,
    name: str,
    evidence_id: str,
    confidence: str,
    audit: str,
    tags: set[str] | None = None,
) -> str:
    key = node_key(label, name)
    if key not in nodes:
        nodes[key] = {
            "id": key,
            "label": label,
            "name": clean(name),
            "count": 0,
            "confidence_counts": Counter(),
            "audit_counts": Counter(),
            "sample_evidence": [],
            "tags": set(tags or []),
        }
    node = nodes[key]
    node["count"] += 1
    node["confidence_counts"][confidence] += 1
    node["audit_counts"][audit] += 1
    node["tags"].update(tags or set())
    if len(node["sample_evidence"]) < 4 and evidence_id:
        node["sample_evidence"].append(evidence_id)
    return key


def add_edge(
    edges: dict[tuple[str, str, str], dict[str, Any]],
    source: str,
    target: str,
    edge_type: str,
    row: dict[str, str],
    confidence: str,
    audit: str,
    tags: set[str] | None = None,
) -> None:
    key = (source, target, edge_type)
    if key not in edges:
        edges[key] = {
            "source": source,
            "target": target,
            "type": edge_type,
            "label": relation_label(edge_type),
            "count": 0,
            "confidence_counts": Counter(),
            "audit_counts": Counter(),
            "samples": [],
            "tags": set(tags or []),
        }
    edge = edges[key]
    edge["count"] += 1
    edge["confidence_counts"][confidence] += 1
    edge["audit_counts"][audit] += 1
    edge["tags"].update(tags or set())
    if len(edge["samples"]) < 4:
        edge["samples"].append(
            {
                "evidence_id": clean(row.get("evidence_id")),
                "title": clipped(row.get("title"), 120),
                "confidence_score": confidence,
                "audit_priority": audit,
                "source_file": clipped(row.get("source_file"), 100),
                "source_sheet": clipped(row.get("source_sheet"), 60),
                "source_row": clean(row.get("source_row")),
                "doi": clean(row.get("doi")),
                "pmid": clean(row.get("pmid")),
            }
        )


def add_rule_sample(edge: dict[str, Any], row: dict[str, str]) -> None:
    if len(edge["samples"]) >= 4:
        return
    edge["samples"].append(
        {
            "evidence_id": clean(row.get("rule_id") or row.get("baseline_id")),
            "title": clipped(row.get("supporting_titles"), 120),
            "confidence_score": "High" if clean(row.get("rule_status")) == "high_confidence_rule" else "Medium",
            "audit_priority": "rule_review" if "mutex" in clean(row.get("rule_type")).lower() else "no_review_needed",
            "source_file": clipped(row.get("supporting_source_document_ids"), 100),
            "source_sheet": "",
            "source_row": clipped(row.get("supporting_event_ids"), 60),
            "doi": clean(row.get("supporting_dois")),
            "pmid": clean(row.get("supporting_pmids")),
        }
    )


def add_rule_edge(
    edges: dict[tuple[str, str, str], dict[str, Any]],
    source: str,
    target: str,
    edge_type: str,
    row: dict[str, str],
    count: int,
    confidence: str,
    tags: set[str],
) -> None:
    key = (source, target, edge_type)
    if key not in edges:
        edges[key] = {
            "source": source,
            "target": target,
            "type": edge_type,
            "label": relation_label(edge_type),
            "count": 0,
            "confidence_counts": Counter(),
            "audit_counts": Counter(),
            "samples": [],
            "tags": set(tags),
        }
    edge = edges[key]
    edge["count"] += max(1, count)
    edge["confidence_counts"][confidence] += max(1, count)
    edge["audit_counts"]["rule_review" if "mutex" in edge_type.lower() else "no_review_needed"] += 1
    edge["tags"].update(tags)
    add_rule_sample(edge, row)


def add_rule_node(
    nodes: dict[str, dict[str, Any]],
    label: str,
    name: str,
    row_id: str,
    confidence: str,
    tags: set[str],
    count: int = 1,
) -> str:
    key = node_key(label, name)
    if key not in nodes:
        nodes[key] = {
            "id": key,
            "label": label,
            "name": clean(name),
            "count": 0,
            "confidence_counts": Counter(),
            "audit_counts": Counter(),
            "sample_evidence": [],
            "tags": set(tags),
        }
    node = nodes[key]
    node["count"] += max(1, count)
    node["confidence_counts"][confidence] += max(1, count)
    node["audit_counts"]["rule_review" if "mutex" in tags else "no_review_needed"] += 1
    node["tags"].update(tags)
    if len(node["sample_evidence"]) < 4 and row_id:
        node["sample_evidence"].append(row_id)
    return key


def make_outcomes(row: dict[str, str]) -> list[tuple[str, str]]:
    outcomes: list[tuple[str, str]] = []
    mrd = clean(row.get("mrd_relevance"))
    if clean(row.get("recurrence_association")) == "reported" or mrd == "recurrence_marker":
        outcomes.append(("ClinicalOutcome", "recurrence"))
    if clean(row.get("relapse_association")) == "reported" or mrd == "relapse_driver":
        outcomes.append(("ClinicalOutcome", "relapse"))
    if clean(row.get("prognosis_association")) == "reported" or is_reported(row.get("survival_endpoint")):
        endpoint = clean(row.get("survival_endpoint"))
        outcomes.append(("ClinicalOutcome", endpoint if is_reported(endpoint) else "prognosis"))
    if clean(row.get("treatment_monitoring_value")) == "reported":
        outcomes.append(("ClinicalOutcome", "treatment_response"))
    dedup: dict[str, tuple[str, str]] = {}
    for label, name in outcomes:
        dedup[node_key(label, name)] = (label, name)
    return list(dedup.values())


def collect_summary(
    scored_csv: Path,
    sources_csv: Path,
    chunks_csv: Path,
    nodes_csv: Path,
    relationships_csv: Path,
    review_csv: Path,
    rag_jsonl: Path,
    raw_csv: Path,
) -> dict[str, Any]:
    confidence = Counter()
    inclusion = Counter()
    audit = Counter()
    source_type = Counter()
    evidence_type = Counter()
    clinical = Counter()
    bio = Counter()
    genes = Counter()
    cancers = Counter()
    assays = Counter()
    supplement_evidence = 0
    sheet_row_evidence = 0

    with scored_csv.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            confidence[clean(row.get("confidence_score"))] += 1
            inclusion[clean(row.get("graph_inclusion_status"))] += 1
            audit[clean(row.get("audit_priority"))] += 1
            evidence_type[clean(row.get("evidence_source_type_final"))] += 1
            clinical[clean(row.get("clinical_relevance_status"))] += 1
            bio[clean(row.get("biological_plausibility_status"))] += 1
            for gene in split_terms(row.get("gene_symbol_normalized") or row.get("gene"), 12):
                genes[gene] += 1
            if is_reported(row.get("cancer_type")):
                cancers[clean(row.get("cancer_type"))] += 1
            if is_reported(row.get("assay_type")):
                assays[clean(row.get("assay_type"))] += 1
            source_document_type = clean(row.get("source_document_type"))
            if source_document_type.startswith("supplementary") or clean(row.get("source_section")).startswith("Supplementary"):
                supplement_evidence += 1
            if is_reported(row.get("source_sheet")) or is_reported(row.get("source_row")):
                sheet_row_evidence += 1

    parser_status = Counter()
    source_docs = Counter()
    capped_sources = 0
    with sources_csv.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            parser_status[clean(row.get("parser_status"))] += 1
            source_docs[clean(row.get("source_document_type"))] += 1
            source_type[clean(row.get("source_document_type"))] += 1
            if clean(row.get("evidence_count")) == "500":
                capped_sources += 1

    chunk_status = Counter()
    rows_seen = 0
    candidate_rows = 0
    with chunks_csv.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            chunk_status[clean(row.get("parser_status"))] += 1
            rows_seen += safe_int(row.get("rows_seen"))
            candidate_rows += safe_int(row.get("candidate_rows"))

    node_labels = Counter()
    with nodes_csv.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            node_labels[clean(row.get("node_label"))] += 1

    rel_types = Counter()
    rel_confidence = Counter()
    with relationships_csv.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            rel_types[clean(row.get("relationship_type"))] += 1
            rel_confidence[clean(row.get("confidence_score"))] += 1

    review_priority = Counter()
    with review_csv.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            review_priority[clean(row.get("audit_priority"))] += 1

    return {
        "counts": {
            "raw_evidence": csv_row_count(raw_csv),
            "scored_evidence": csv_row_count(scored_csv),
            "source_files": csv_row_count(sources_csv),
            "excel_chunks": csv_row_count(chunks_csv),
            "graph_nodes": csv_row_count(nodes_csv),
            "graph_relationships": csv_row_count(relationships_csv),
            "review_tasks": csv_row_count(review_csv),
            "rag_summaries": jsonl_count(rag_jsonl),
            "supplementary_evidence": supplement_evidence,
            "sheet_row_evidence": sheet_row_evidence,
            "excel_rows_seen": rows_seen,
            "excel_candidate_rows": candidate_rows,
            "capped_sources": capped_sources,
        },
        "confidence": confidence,
        "inclusion": inclusion,
        "audit": audit,
        "evidence_type": evidence_type,
        "clinical_relevance": clinical,
        "biological_plausibility": bio,
        "source_document_type": source_docs,
        "parser_status": parser_status,
        "chunk_status": chunk_status,
        "node_labels": node_labels,
        "relationship_types": rel_types,
        "relationship_confidence": rel_confidence,
        "review_priority": review_priority,
        "top_genes": genes,
        "top_cancers": cancers,
        "top_assays": assays,
    }


def rule_rank(row: dict[str, str]) -> float:
    try:
        pair_count = float(clean(row.get("pair_count")))
    except ValueError:
        pair_count = 0
    try:
        sample_count = float(clean(row.get("sample_count")))
    except ValueError:
        sample_count = 0
    try:
        lift = float(clean(row.get("lift")))
    except ValueError:
        lift = 0
    p_text = clean(row.get("fisher_p_value"))
    try:
        p_value = float(p_text)
    except ValueError:
        p_value = 1
    status_bonus = 120 if clean(row.get("rule_status")) == "high_confidence_rule" else 0
    return status_bonus + pair_count * 2 + min(lift, 100) + max(0, -math.log10(max(p_value, 1e-300))) + sample_count / 100000


def baseline_rank(row: dict[str, str]) -> float:
    try:
        count = float(clean(row.get("gene_sample_count")))
    except ValueError:
        count = 0
    try:
        freq = float(clean(row.get("gene_frequency")))
    except ValueError:
        freq = 0
    return count + freq * 500


def top_rule_rows(path: Path, limit: int, kind: str) -> list[dict[str, str]]:
    if not path.exists():
        return []
    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        for row in csv.DictReader(handle):
            if not clean(row.get("gene_a")) or not clean(row.get("gene_b")):
                continue
            if kind == "cooccurrence" and clean(row.get("rule_status")) != "high_confidence_rule":
                continue
            rows.append(row)
    return sorted(rows, key=rule_rank, reverse=True)[:limit]


def top_baseline_rows(path: Path, limit: int) -> list[dict[str, str]]:
    if not path.exists():
        return []
    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        for row in csv.DictReader(handle):
            if clean(row.get("feature_type")) != "gene" or not clean(row.get("gene")):
                continue
            rows.append(row)
    return sorted(rows, key=baseline_rank, reverse=True)[:limit]


def add_rule_layers(
    node_stats: dict[str, dict[str, Any]],
    edge_stats: dict[tuple[str, str, str], dict[str, Any]],
    baseline_csv: Path,
    cooccurrence_csv: Path,
    mutex_csv: Path,
) -> dict[str, int]:
    added = {"baseline_rows": 0, "cooccurrence_rows": 0, "mutex_rows": 0}
    for row in top_baseline_rows(baseline_csv, 260):
        cancer = clean(row.get("cancer_type"))
        gene = clean(row.get("gene"))
        row_id = clean(row.get("baseline_id"))
        if not cancer or not gene:
            continue
        try:
            count = safe_int(row.get("gene_sample_count"))
        except ValueError:
            count = 1
        tags = {"rules", "baseline"}
        cancer_id = add_rule_node(node_stats, "Cancer", cancer, row_id, "Medium", tags, count)
        gene_id = add_rule_node(node_stats, "Gene", gene, row_id, "Medium", tags, count)
        add_rule_edge(edge_stats, cancer_id, gene_id, "CANCER_HAS_BASELINE", row, count, "Medium", tags)
        added["baseline_rows"] += 1

    for row in top_rule_rows(cooccurrence_csv, 340, "cooccurrence"):
        gene_a = clean(row.get("gene_a"))
        gene_b = clean(row.get("gene_b"))
        cancer = clean(row.get("cancer_type"))
        row_id = clean(row.get("rule_id"))
        count = safe_int(row.get("pair_count"))
        tags = {"rules", "cooccurrence"}
        confidence = "High" if clean(row.get("rule_status")) == "high_confidence_rule" else "Medium"
        gene_a_id = add_rule_node(node_stats, "Gene", gene_a, row_id, confidence, tags, count)
        gene_b_id = add_rule_node(node_stats, "Gene", gene_b, row_id, confidence, tags, count)
        rule_name = f"{cancer}: {gene_a}-{gene_b} cooccur"
        rule_id = add_rule_node(node_stats, "GeneAssociationRule", rule_name, row_id, confidence, tags, count)
        cancer_id = add_rule_node(node_stats, "Cancer", cancer, row_id, confidence, tags, count)
        add_rule_edge(edge_stats, gene_a_id, gene_b_id, "GENE_COOCCURS_WITH", row, count, confidence, tags)
        add_rule_edge(edge_stats, rule_id, gene_a_id, "RULE_INVOLVES_GENE", row, count, confidence, tags)
        add_rule_edge(edge_stats, rule_id, gene_b_id, "RULE_INVOLVES_GENE", row, count, confidence, tags)
        add_rule_edge(edge_stats, cancer_id, rule_id, "RULE_IN_CANCER", row, count, confidence, tags)
        added["cooccurrence_rows"] += 1

    for row in top_rule_rows(mutex_csv, 220, "mutex"):
        gene_a = clean(row.get("gene_a"))
        gene_b = clean(row.get("gene_b"))
        cancer = clean(row.get("cancer_type"))
        row_id = clean(row.get("rule_id"))
        count = max(1, safe_int(row.get("gene_a_count")) + safe_int(row.get("gene_b_count")))
        tags = {"rules", "mutex"}
        gene_a_id = add_rule_node(node_stats, "Gene", gene_a, row_id, "Medium", tags, count)
        gene_b_id = add_rule_node(node_stats, "Gene", gene_b, row_id, "Medium", tags, count)
        rule_name = f"{cancer}: {gene_a}-{gene_b} mutex"
        rule_id = add_rule_node(node_stats, "GeneAssociationRule", rule_name, row_id, "Medium", tags, count)
        cancer_id = add_rule_node(node_stats, "Cancer", cancer, row_id, "Medium", tags, count)
        add_rule_edge(edge_stats, gene_a_id, gene_b_id, "GENE_MUTUALLY_EXCLUSIVE_WITH", row, count, "Medium", tags)
        add_rule_edge(edge_stats, rule_id, gene_a_id, "RULE_INVOLVES_GENE", row, count, "Medium", tags)
        add_rule_edge(edge_stats, rule_id, gene_b_id, "RULE_INVOLVES_GENE", row, count, "Medium", tags)
        add_rule_edge(edge_stats, cancer_id, rule_id, "RULE_IN_CANCER", row, count, "Medium", tags)
        added["mutex_rows"] += 1
    return added


def build_network(
    scored_csv: Path,
    baseline_csv: Path,
    cooccurrence_csv: Path,
    mutex_csv: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    node_stats: dict[str, dict[str, Any]] = {}
    edge_stats: dict[tuple[str, str, str], dict[str, Any]] = {}
    gene_counts = Counter()
    cancer_counts = Counter()
    assay_counts = Counter()
    variant_counts = Counter()
    mrd_counts = Counter()
    outcome_counts = Counter()

    rows: list[dict[str, str]] = []
    with scored_csv.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            if clean(row.get("graph_inclusion_status")) == "exclude":
                continue
            rows.append(row)
            for gene in split_terms(row.get("gene_symbol_normalized") or row.get("gene"), 10):
                gene_counts[gene] += 1
            if is_reported(row.get("cancer_type")):
                cancer_counts[clean(row.get("cancer_type"))] += 1
            if is_reported(row.get("assay_type")):
                assay_counts[clean(row.get("assay_type"))] += 1
            for variant in split_terms(row.get("variant"), 5):
                variant_counts[variant] += 1
            if is_reported(row.get("mrd_relevance")) and clean(row.get("mrd_relevance")) != "uncertain":
                mrd_counts[clean(row.get("mrd_relevance"))] += 1
            for _, outcome in make_outcomes(row):
                outcome_counts[outcome] += 1

    selected_genes = {gene for gene, _ in gene_counts.most_common(70)}
    selected_cancers = {cancer for cancer, _ in cancer_counts.most_common(18)}
    selected_assays = {assay for assay, _ in assay_counts.most_common(16)}
    selected_mrd = {mrd for mrd, _ in mrd_counts.most_common(14)}
    selected_outcomes = {outcome for outcome, _ in outcome_counts.most_common(10)}
    selected_variants = {variant for variant, count in variant_counts.most_common(95) if count >= 2}

    for row in rows:
        evidence_id = clean(row.get("evidence_id"))
        confidence = clean(row.get("confidence_score"))
        audit = clean(row.get("audit_priority"))
        tags: set[str] = set()
        if audit == "urgent_review":
            tags.add("urgent")
        if clean(row.get("biological_plausibility_status")) == "likely_CHIP" or "chip" in " ".join(
            [
                clean(row.get("mrd_relevance")),
                clean(row.get("biomarker_type")),
                clean(row.get("quoted_or_paraphrased_evidence")),
            ]
        ).lower():
            tags.add("chip")
        if clean(row.get("source_document_type")).startswith("supplementary"):
            tags.add("supplementary")

        genes = [g for g in split_terms(row.get("gene_symbol_normalized") or row.get("gene"), 10) if g in selected_genes]
        cancers = [clean(row.get("cancer_type"))] if clean(row.get("cancer_type")) in selected_cancers else []
        assays = [clean(row.get("assay_type"))] if clean(row.get("assay_type")) in selected_assays else []
        variants = [v for v in split_terms(row.get("variant"), 5) if v in selected_variants]
        mrd_terms = [clean(row.get("mrd_relevance"))] if clean(row.get("mrd_relevance")) in selected_mrd else []
        outcomes = [(label, outcome) for label, outcome in make_outcomes(row) if outcome in selected_outcomes]

        gene_ids = [add_node(node_stats, "Gene", gene, evidence_id, confidence, audit, tags) for gene in genes]
        cancer_ids = [add_node(node_stats, "Cancer", cancer, evidence_id, confidence, audit, tags) for cancer in cancers]
        assay_ids = [add_node(node_stats, "Assay", assay, evidence_id, confidence, audit, tags) for assay in assays]
        variant_ids = [add_node(node_stats, "Variant", variant, evidence_id, confidence, audit, tags) for variant in variants]
        mrd_ids = [add_node(node_stats, "MRDSignal", mrd, evidence_id, confidence, audit, tags) for mrd in mrd_terms]
        outcome_ids = [
            add_node(node_stats, label, outcome, evidence_id, confidence, audit, tags) for label, outcome in outcomes
        ]

        for cancer_id in cancer_ids:
            for gene_id in gene_ids:
                add_edge(edge_stats, cancer_id, gene_id, "ASSOCIATED_WITH", row, confidence, audit, tags)
        for gene_id in gene_ids:
            for variant_id in variant_ids:
                add_edge(edge_stats, gene_id, variant_id, "HAS_VARIANT", row, confidence, audit, tags)
        for variant_id in variant_ids:
            for cancer_id in cancer_ids:
                add_edge(edge_stats, variant_id, cancer_id, "OBSERVED_IN", row, confidence, audit, tags)
        for gene_id in gene_ids:
            for mrd_id in mrd_ids:
                add_edge(edge_stats, gene_id, mrd_id, "HAS_MRD_RELEVANCE", row, confidence, audit, tags)
        for variant_id in variant_ids:
            for mrd_id in mrd_ids:
                add_edge(edge_stats, variant_id, mrd_id, "HAS_MRD_RELEVANCE_VARIANT", row, confidence, audit, tags)
        for mrd_id in mrd_ids:
            for outcome_id in outcome_ids:
                outcome_name = node_stats[outcome_id]["name"]
                rel = {
                    "recurrence": "ASSOCIATED_WITH_RECURRENCE",
                    "relapse": "ASSOCIATED_WITH_RELAPSE",
                    "treatment_response": "USED_FOR_TREATMENT_MONITORING",
                }.get(outcome_name, "ASSOCIATED_WITH_PROGNOSIS")
                add_edge(edge_stats, mrd_id, outcome_id, rel, row, confidence, audit, tags)
        for assay_id in assay_ids:
            for gene_id in gene_ids:
                add_edge(edge_stats, assay_id, gene_id, "DETECTS", row, confidence, audit, tags)
            for variant_id in variant_ids:
                add_edge(edge_stats, assay_id, variant_id, "DETECTS_VARIANT", row, confidence, audit, tags)
        for mrd_id in mrd_ids:
            for assay_id in assay_ids:
                add_edge(edge_stats, mrd_id, assay_id, "DETECTED_BY", row, confidence, audit, tags)

    rule_meta = add_rule_layers(node_stats, edge_stats, baseline_csv, cooccurrence_csv, mutex_csv)

    selected_edges: list[dict[str, Any]] = []
    by_type: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for edge in edge_stats.values():
        edge["score"] = edge_score(edge)
        by_type[edge["type"]].append(edge)
    for edge_type, items in by_type.items():
        limit = DEFAULT_EDGE_LIMITS.get(edge_type, 120)
        selected_edges.extend(sorted(items, key=lambda item: item["score"], reverse=True)[:limit])

    endpoint_ids = {edge["source"] for edge in selected_edges} | {edge["target"] for edge in selected_edges}
    selected_nodes = [node for node_id, node in node_stats.items() if node_id in endpoint_ids]

    layout_nodes(selected_nodes)
    final_node_ids = {node["id"] for node in selected_nodes}
    selected_edges = [
        edge for edge in selected_edges if edge["source"] in final_node_ids and edge["target"] in final_node_ids
    ]

    nodes = [serialize_node(node) for node in selected_nodes]
    edges = [serialize_edge(edge, index) for index, edge in enumerate(selected_edges)]
    meta = {
        "selected_genes": len(selected_genes),
        "selected_cancers": len(selected_cancers),
        "selected_assays": len(selected_assays),
        "selected_variants": len(selected_variants),
        "selected_mrd_signals": len(selected_mrd),
        "selected_outcomes": len(selected_outcomes),
        "aggregated_edges_before_pruning": len(edge_stats),
        "visual_edges": len(edges),
        "visual_nodes": len(nodes),
        "visual_baseline_rules": rule_meta["baseline_rows"],
        "visual_cooccurrence_rules": rule_meta["cooccurrence_rows"],
        "visual_mutex_rules": rule_meta["mutex_rows"],
    }
    return nodes, edges, meta


def layout_nodes(nodes: list[dict[str, Any]]) -> None:
    theta_by_label = {
        "Cancer": 0.36 * math.pi,
        "Gene": 0.50 * math.pi,
        "MRDSignal": 0.57 * math.pi,
        "ClinicalOutcome": 0.42 * math.pi,
        "Assay": 0.32 * math.pi,
        "Variant": 0.68 * math.pi,
    }
    offset_by_label = {
        "Cancer": 0.15,
        "Gene": 0.85,
        "MRDSignal": 1.55,
        "ClinicalOutcome": 2.20,
        "Assay": 2.85,
        "Variant": 3.55,
    }
    golden_angle = math.pi * (3 - math.sqrt(5))
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for node in nodes:
        grouped[node["label"]].append(node)

    for label, items in grouped.items():
        items.sort(key=lambda item: (-item["count"], item["name"]))
        theta_center = theta_by_label.get(label, 0.55 * math.pi)
        offset = offset_by_label.get(label, 0.0)
        spread = 0.42 if len(items) > 18 else 0.24
        for index, node in enumerate(items):
            centered = 0 if len(items) <= 1 else (index / (len(items) - 1) - 0.5)
            theta = max(0.08, min(math.pi - 0.08, theta_center + centered * spread))
            phi = offset + index * golden_angle
            shell = 0.88 + 0.12 * ((index % 5) / 4)
            node["sx"] = shell * math.sin(theta) * math.cos(phi)
            node["sy"] = shell * math.cos(theta)
            node["sz"] = shell * math.sin(theta) * math.sin(phi)
            node["x"] = node["sx"]
            node["y"] = node["sy"]
            node["z"] = node["sz"]
            node["radius"] = 6 + min(18, math.log1p(node["count"]) * 2.2)


def serialize_counter(counter: Counter[str], limit: int | None = None) -> dict[str, int]:
    items = counter.most_common(limit) if limit else counter.items()
    return {key: value for key, value in items}


def serialize_node(node: dict[str, Any]) -> dict[str, Any]:
    confidence = serialize_counter(node["confidence_counts"])
    dominant = dominant_confidence(node["confidence_counts"])
    return {
        "id": node["id"],
        "label": node["label"],
        "name": node["name"],
        "count": node["count"],
        "confidence_counts": confidence,
        "audit_counts": serialize_counter(node["audit_counts"]),
        "dominant_confidence": dominant,
        "color": NODE_COLORS.get(node["label"], "#5c6672"),
        "x": round(node["x"], 4),
        "y": round(node["y"], 4),
        "z": round(node["z"], 4),
        "sx": round(node["sx"], 4),
        "sy": round(node["sy"], 4),
        "sz": round(node["sz"], 4),
        "radius": round(node["radius"], 2),
        "sample_evidence": node["sample_evidence"],
        "tags": sorted(node["tags"]),
    }


def serialize_edge(edge: dict[str, Any], index: int) -> dict[str, Any]:
    confidence = serialize_counter(edge["confidence_counts"])
    dominant = dominant_confidence(edge["confidence_counts"])
    return {
        "id": f"E{index:05d}",
        "source": edge["source"],
        "target": edge["target"],
        "type": edge["type"],
        "label": edge["label"],
        "count": edge["count"],
        "confidence_counts": confidence,
        "audit_counts": serialize_counter(edge["audit_counts"]),
        "dominant_confidence": dominant,
        "color": EDGE_TYPE_COLORS.get(edge["type"], EDGE_COLORS.get(dominant, "#7a808a")),
        "samples": edge["samples"],
        "tags": sorted(edge["tags"]),
        "score": round(edge["score"], 2),
    }


def serializable_summary(summary: dict[str, Any]) -> dict[str, Any]:
    return {
        "counts": summary["counts"],
        "confidence": serialize_counter(summary["confidence"]),
        "inclusion": serialize_counter(summary["inclusion"]),
        "audit": serialize_counter(summary["audit"]),
        "evidence_type": serialize_counter(summary["evidence_type"]),
        "clinical_relevance": serialize_counter(summary["clinical_relevance"]),
        "biological_plausibility": serialize_counter(summary["biological_plausibility"]),
        "source_document_type": serialize_counter(summary["source_document_type"], 16),
        "parser_status": serialize_counter(summary["parser_status"]),
        "chunk_status": serialize_counter(summary["chunk_status"]),
        "node_labels": serialize_counter(summary["node_labels"], 16),
        "relationship_types": serialize_counter(summary["relationship_types"], 20),
        "relationship_confidence": serialize_counter(summary["relationship_confidence"]),
        "review_priority": serialize_counter(summary["review_priority"]),
        "top_genes": serialize_counter(summary["top_genes"], 20),
        "top_cancers": serialize_counter(summary["top_cancers"], 16),
        "top_assays": serialize_counter(summary["top_assays"], 16),
    }


def build_visualization_data(args: argparse.Namespace) -> dict[str, Any]:
    root = artifact_root() if args.root == "." else Path(args.root)
    paths = {
        "raw_csv": artifact_path(root / args.raw_csv),
        "scored_csv": artifact_path(root / args.scored_csv),
        "sources_csv": artifact_path(root / args.sources_csv),
        "chunks_csv": artifact_path(root / args.chunks_csv),
        "nodes_csv": artifact_path(root / args.nodes_csv),
        "relationships_csv": artifact_path(root / args.relationships_csv),
        "review_csv": artifact_path(root / args.review_csv),
        "rag_jsonl": artifact_path(root / args.rag_jsonl),
        "baseline_csv": artifact_path(root / args.baseline_csv),
        "cooccurrence_csv": artifact_path(root / args.cooccurrence_csv),
        "mutex_csv": artifact_path(root / args.mutex_csv),
    }
    summary_paths = {key: paths[key] for key in ["scored_csv", "sources_csv", "chunks_csv", "nodes_csv", "relationships_csv", "review_csv", "rag_jsonl", "raw_csv"]}
    summary = collect_summary(**summary_paths)
    nodes, edges, meta = build_network(paths["scored_csv"], paths["baseline_csv"], paths["cooccurrence_csv"], paths["mutex_csv"])
    summary["counts"]["baseline_rules"] = csv_row_count(paths["baseline_csv"]) if paths["baseline_csv"].exists() else 0
    summary["counts"]["cooccurrence_rules"] = csv_row_count(paths["cooccurrence_csv"]) if paths["cooccurrence_csv"].exists() else 0
    summary["counts"]["mutex_rules"] = csv_row_count(paths["mutex_csv"]) if paths["mutex_csv"].exists() else 0
    return {
        "title": "MRD Full Corpus Support KG",
        "owner": "Phd.W",
        "generated_from": {name: str(path) for name, path in paths.items()},
        "summary": serializable_summary(summary),
        "network_meta": meta,
        "nodes": nodes,
        "edges": edges,
        "relation_types": sorted({edge["type"] for edge in edges}),
        "node_labels": sorted({node["label"] for node in nodes}),
    }


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def js_data(data: dict[str, Any]) -> str:
    return json.dumps(data, ensure_ascii=False).replace("</", "<\\/")


HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>MRD Full Corpus Support KG</title>
  <style>
    :root {
      --bg: #eef2f6;
      --ink: #18212f;
      --muted: #667085;
      --line: #d8dee8;
      --panel: #ffffff;
      --accent: #2b7c86;
      --focus: #174a54;
      --high: #25905a;
      --medium: #b8860b;
      --low: #b94a48;
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      background: var(--bg);
      color: var(--ink);
      font: 14px/1.45 "Segoe UI", Arial, sans-serif;
    }
    header {
      display: grid;
      grid-template-columns: 1fr auto;
      gap: 16px;
      align-items: center;
      padding: 14px 18px;
      border-bottom: 1px solid var(--line);
      background: #fff;
      position: sticky;
      top: 0;
      z-index: 10;
    }
    h1 {
      font-size: 18px;
      margin: 0;
      font-weight: 650;
      letter-spacing: 0;
    }
    .subhead {
      margin-top: 2px;
      color: var(--muted);
      font-size: 12px;
    }
    .top-actions {
      display: flex;
      gap: 8px;
      align-items: center;
      flex-wrap: wrap;
      justify-content: end;
    }
    button, select, input[type="search"] {
      border: 1px solid var(--line);
      background: #fff;
      color: var(--ink);
      border-radius: 6px;
      padding: 7px 9px;
      font: inherit;
      min-height: 34px;
    }
    button:hover, select:hover { border-color: #a7b0c0; }
    button:focus, select:focus, input:focus {
      outline: 2px solid rgba(45,108,223,.18);
      border-color: var(--accent);
    }
    input[type="search"] { width: min(28vw, 280px); }
    .workspace {
      display: grid;
      grid-template-columns: 300px minmax(520px, 1fr) 330px;
      gap: 12px;
      padding: 12px;
      min-height: calc(100vh - 66px);
    }
    .panel {
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      overflow: hidden;
      min-width: 0;
    }
    .panel h2 {
      margin: 0;
      padding: 10px 12px;
      font-size: 13px;
      letter-spacing: 0;
      border-bottom: 1px solid var(--line);
      background: #fafbfc;
    }
    .panel-body { padding: 10px 12px; }
    .stat-grid {
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 8px;
    }
    .stat {
      border: 1px solid #e5e9f0;
      border-radius: 7px;
      padding: 8px;
      min-width: 0;
      background: #fbfcfe;
    }
    .stat-value {
      font-size: 19px;
      font-weight: 700;
      white-space: nowrap;
    }
    .stat-label {
      color: var(--muted);
      font-size: 11px;
      margin-top: 2px;
    }
    .controls {
      display: grid;
      gap: 10px;
    }
    .control-row {
      display: grid;
      grid-template-columns: 1fr auto;
      gap: 8px;
      align-items: center;
    }
    .checks {
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 6px;
      color: #344054;
    }
    label {
      display: flex;
      align-items: center;
      gap: 6px;
      min-width: 0;
    }
    input[type="range"] { width: 100%; }
    .bar-list {
      display: grid;
      gap: 6px;
      margin-top: 8px;
    }
    .bar-row {
      display: grid;
      grid-template-columns: minmax(84px, 1fr) 58px;
      gap: 8px;
      align-items: center;
      min-height: 22px;
    }
    .bar-track {
      height: 8px;
      background: #edf1f6;
      border-radius: 999px;
      overflow: hidden;
      margin-top: 2px;
    }
    .bar-fill {
      height: 100%;
      background: var(--accent);
      width: 0;
    }
    .bar-label {
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
      font-size: 12px;
    }
    .bar-value {
      text-align: right;
      color: var(--muted);
      font-variant-numeric: tabular-nums;
      font-size: 12px;
    }
    .canvas-panel {
      display: grid;
      grid-template-rows: auto minmax(0, 1fr);
      align-self: start;
      height: calc(100vh - 92px);
      min-height: 760px;
      max-height: 920px;
    }
    .canvas-toolbar {
      display: flex;
      gap: 8px;
      align-items: center;
      padding: 10px 12px;
      border-bottom: 1px solid var(--line);
      background: #fff;
      flex-wrap: wrap;
    }
    .canvas-wrap {
      position: relative;
      background:
        radial-gradient(circle at 50% 48%, #ffffff 0%, #f7fafc 58%, #eef3f8 100%);
      min-height: 0;
    }
    canvas {
      display: block;
      width: 100%;
      height: 100%;
      min-height: 0;
      cursor: grab;
    }
    canvas.dragging { cursor: grabbing; }
    .sphere-hint {
      position: absolute;
      right: 12px;
      top: 12px;
      padding: 6px 8px;
      border: 1px solid rgba(91,111,132,.24);
      border-radius: 999px;
      color: rgba(38,51,68,.78);
      background: rgba(255,255,255,.78);
      font-size: 11px;
      pointer-events: none;
    }
    .legend {
      position: absolute;
      left: 12px;
      bottom: 12px;
      display: flex;
      gap: 8px;
      flex-wrap: wrap;
      max-width: calc(100% - 24px);
      pointer-events: none;
    }
    .legend-item {
      background: rgba(255,255,255,.86);
      color: #243244;
      border: 1px solid rgba(91,111,132,.20);
      border-radius: 999px;
      padding: 4px 8px;
      font-size: 11px;
      display: inline-flex;
      align-items: center;
      gap: 5px;
    }
    .dot {
      width: 9px;
      height: 9px;
      border-radius: 50%;
      display: inline-block;
    }
    .detail-title {
      font-weight: 700;
      font-size: 15px;
      margin-bottom: 6px;
      overflow-wrap: anywhere;
    }
    .kv {
      display: grid;
      grid-template-columns: 112px 1fr;
      gap: 4px 8px;
      font-size: 12px;
      margin: 8px 0;
    }
    .kv div:nth-child(odd) { color: var(--muted); }
    .samples {
      display: grid;
      gap: 8px;
      margin-top: 10px;
    }
    .sample {
      border-top: 1px solid #edf0f4;
      padding-top: 8px;
      font-size: 12px;
    }
    .sample strong { display: block; margin-bottom: 3px; }
    .table {
      width: 100%;
      border-collapse: collapse;
      font-size: 12px;
    }
    .table th, .table td {
      border-bottom: 1px solid #edf0f4;
      padding: 6px 4px;
      text-align: left;
      vertical-align: top;
    }
    .table th {
      color: var(--muted);
      font-weight: 600;
    }
    .pill {
      display: inline-block;
      border-radius: 999px;
      border: 1px solid var(--line);
      padding: 2px 6px;
      font-size: 11px;
      color: #344054;
      background: #fff;
    }
    .scroll { max-height: 340px; overflow: auto; }
    .muted { color: var(--muted); }
    .section-gap { margin-top: 12px; }
    @media (max-width: 1180px) {
      .workspace { grid-template-columns: 1fr; }
      input[type="search"] { width: min(70vw, 320px); }
      .canvas-panel { min-height: 720px; }
    }
  </style>
</head>
<body>
  <header>
    <div>
      <h1>MRD Full Corpus Support KG</h1>
      <div class="subhead">Owner Phd.W - candidate support graph - full-corpus build</div>
    </div>
    <div class="top-actions">
      <input id="searchBox" type="search" placeholder="Search node">
      <button id="resetView" type="button" title="Reset viewport">Reset</button>
    </div>
  </header>

  <main class="workspace">
    <aside class="panel">
      <h2>Build Summary</h2>
      <div class="panel-body">
        <div id="stats" class="stat-grid"></div>

        <div class="section-gap">
          <h2 style="margin:0 -12px 0 -12px;">Confidence</h2>
          <div id="confidenceBars" class="bar-list"></div>
        </div>

        <div class="section-gap">
          <h2 style="margin:0 -12px 0 -12px;">Evidence Source</h2>
          <div id="evidenceBars" class="bar-list"></div>
        </div>

        <div class="section-gap">
          <h2 style="margin:0 -12px 0 -12px;">Top Genes</h2>
          <div id="geneBars" class="bar-list"></div>
        </div>
      </div>
    </aside>

    <section class="panel canvas-panel">
      <div class="canvas-toolbar">
        <div class="controls" style="grid-template-columns: 170px 210px 1fr; flex:1;">
          <label>
            View
            <select id="viewSelect">
              <option value="core">Core</option>
              <option value="rules" selected>Gene Rules</option>
              <option value="cooccurrence">Co-occurrence</option>
              <option value="mutex">Mutual Exclusion</option>
              <option value="variants">Variants</option>
              <option value="review">Review</option>
              <option value="chip">CHIP</option>
              <option value="all">All</option>
            </select>
          </label>
          <label>
            Min edge evidence
            <input id="minCount" type="range" min="1" max="80" value="1">
            <span id="minCountLabel" class="pill">1</span>
          </label>
          <div class="checks" id="confidenceChecks">
            <label><input type="checkbox" value="High" checked> High</label>
            <label><input type="checkbox" value="Medium" checked> Medium</label>
            <label><input type="checkbox" value="Low" checked> Low</label>
            <label><input type="checkbox" value="Exclude"> Exclude</label>
          </div>
        </div>
      </div>
      <div class="canvas-wrap">
        <canvas id="graphCanvas"></canvas>
        <div class="sphere-hint">Drag rotate - Wheel zoom - Click node</div>
        <div id="legend" class="legend"></div>
      </div>
    </section>

    <aside class="panel">
      <h2>Selection</h2>
      <div id="detail" class="panel-body">
        <div class="muted">No selection</div>
      </div>
      <h2>Top Relationships</h2>
      <div class="panel-body scroll">
        <table class="table" id="relationshipTable">
          <thead><tr><th>Relation</th><th>Evidence</th><th>Confidence</th></tr></thead>
          <tbody></tbody>
        </table>
      </div>
      <h2>Top Cancer Types</h2>
      <div class="panel-body">
        <div id="cancerBars" class="bar-list"></div>
      </div>
      <h2>Top Assays</h2>
      <div class="panel-body">
        <div id="assayBars" class="bar-list"></div>
      </div>
    </aside>
  </main>

  <script id="kg-data" type="application/json">__DATA_JSON__</script>
  <script>
    const kg = JSON.parse(document.getElementById('kg-data').textContent);
    const nodesById = new Map(kg.nodes.map(node => [node.id, node]));
    const canvas = document.getElementById('graphCanvas');
    const ctx = canvas.getContext('2d');
    const state = {
      view: 'rules',
      minCount: 1,
      confidence: new Set(['High', 'Medium', 'Low']),
      search: '',
      zoom: 1,
      panX: 0,
      panY: 0,
      rotX: -0.22,
      rotY: 0.62,
      autoRotate: true,
      dragging: false,
      lastX: 0,
      lastY: 0,
      hover: null,
      selected: null
    };

    function fmt(value) {
      return Number(value || 0).toLocaleString('en-US');
    }

    function shortName(value, limit = 26) {
      value = String(value || '');
      return value.length <= limit ? value : value.slice(0, limit - 3) + '...';
    }

    function escapeHtml(value) {
      return String(value ?? '').replace(/[&<>"']/g, ch => ({
        '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'
      }[ch]));
    }

    function renderStats() {
      const counts = kg.summary.counts;
      const stats = [
        ['Evidence', counts.scored_evidence],
        ['Nodes', counts.graph_nodes],
        ['Relationships', counts.graph_relationships],
        ['Co-occur', counts.cooccurrence_rules],
        ['Mutex', counts.mutex_rules],
        ['Baselines', counts.baseline_rules],
        ['Review Tasks', counts.review_tasks],
        ['Source Files', counts.source_files],
        ['Supplementary', counts.supplementary_evidence]
      ];
      document.getElementById('stats').innerHTML = stats.map(([label, value]) => `
        <div class="stat">
          <div class="stat-value">${fmt(value)}</div>
          <div class="stat-label">${escapeHtml(label)}</div>
        </div>
      `).join('');
    }

    function renderBars(id, data, limit = 10, color = '#2d6cdf') {
      const entries = Object.entries(data || {}).sort((a, b) => b[1] - a[1]).slice(0, limit);
      const max = Math.max(1, ...entries.map(item => item[1]));
      document.getElementById(id).innerHTML = entries.map(([label, value]) => {
        const pct = Math.max(3, Math.round(value / max * 100));
        return `
          <div class="bar-row" title="${escapeHtml(label)}: ${fmt(value)}">
            <div>
              <div class="bar-label">${escapeHtml(label)}</div>
              <div class="bar-track"><div class="bar-fill" style="width:${pct}%;background:${color};"></div></div>
            </div>
            <div class="bar-value">${fmt(value)}</div>
          </div>
        `;
      }).join('');
    }

    function visibleByView(edge) {
      if (state.view === 'all') return true;
      if (state.view === 'rules') return edge.tags.includes('rules') || edge.type.includes('RULE') || edge.type.includes('GENE_') || edge.type === 'CANCER_HAS_BASELINE';
      if (state.view === 'cooccurrence') return edge.tags.includes('cooccurrence') || edge.type === 'GENE_COOCCURS_WITH';
      if (state.view === 'mutex') return edge.tags.includes('mutex') || edge.type === 'GENE_MUTUALLY_EXCLUSIVE_WITH';
      if (state.view === 'variants') return edge.type.includes('VARIANT') || edge.type === 'HAS_VARIANT' || edge.type === 'OBSERVED_IN';
      if (state.view === 'review') return Number(edge.audit_counts.urgent_review || 0) > 0 || Number(edge.audit_counts.high_review || 0) > 0;
      if (state.view === 'chip') return edge.tags.includes('chip');
      return !edge.type.includes('VARIANT') && edge.type !== 'HAS_VARIANT' && edge.type !== 'OBSERVED_IN';
    }

    function confidenceVisible(edge) {
      return Array.from(state.confidence).some(conf => Number(edge.confidence_counts[conf] || 0) > 0);
    }

    function visibleEdges() {
      return kg.edges.filter(edge =>
        edge.count >= state.minCount &&
        visibleByView(edge) &&
        confidenceVisible(edge)
      );
    }

    function visibleNodes(edges) {
      const ids = new Set();
      edges.forEach(edge => {
        ids.add(edge.source);
        ids.add(edge.target);
      });
      let nodes = kg.nodes.filter(node => ids.has(node.id));
      const query = state.search.trim().toLowerCase();
      if (query) {
        const matches = new Set(nodes.filter(node => node.name.toLowerCase().includes(query)).map(node => node.id));
        const neighbourIds = new Set(matches);
        edges.forEach(edge => {
          if (matches.has(edge.source)) neighbourIds.add(edge.target);
          if (matches.has(edge.target)) neighbourIds.add(edge.source);
        });
        nodes = nodes.filter(node => neighbourIds.has(node.id));
      }
      return nodes;
    }

    function rotatePoint(node) {
      const cy = Math.cos(state.rotY);
      const sy = Math.sin(state.rotY);
      const cx = Math.cos(state.rotX);
      const sx = Math.sin(state.rotX);
      const x0 = node.sx ?? node.x;
      const y0 = node.sy ?? node.y;
      const z0 = node.sz ?? node.z ?? 0;
      const x1 = x0 * cy + z0 * sy;
      const z1 = -x0 * sy + z0 * cy;
      const y1 = y0 * cx - z1 * sx;
      const z2 = y0 * sx + z1 * cx;
      return { x: x1, y: y1, z: z2 };
    }

    function sphereRadius() {
      const rect = canvas.getBoundingClientRect();
      return Math.max(190, Math.min(rect.width, rect.height) * 0.34 * state.zoom);
    }

    function screen(node) {
      const rect = canvas.getBoundingClientRect();
      const p3 = rotatePoint(node);
      const depth = (p3.z + 1) / 2;
      const perspective = 0.72 + depth * 0.34;
      const radius = sphereRadius();
      return {
        x: rect.width / 2 + state.panX + p3.x * radius * perspective,
        y: rect.height * 0.49 + state.panY + p3.y * radius * perspective,
        z: p3.z,
        depth,
        scale: perspective,
        r: Math.max(3.5, node.radius * 0.54 * perspective * state.zoom)
      };
    }

    function resizeCanvas() {
      const rect = canvas.parentElement.getBoundingClientRect();
      const dpr = window.devicePixelRatio || 1;
      canvas.width = Math.max(300, rect.width) * dpr;
      canvas.height = Math.max(520, rect.height) * dpr;
      canvas.style.width = `${Math.max(300, rect.width)}px`;
      canvas.style.height = `${Math.max(520, rect.height)}px`;
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      draw();
    }

    function hexToRgb(hex) {
      const value = hex.replace('#', '');
      return {
        r: parseInt(value.slice(0, 2), 16),
        g: parseInt(value.slice(2, 4), 16),
        b: parseInt(value.slice(4, 6), 16)
      };
    }

    function rgba(hex, alpha) {
      const c = hexToRgb(hex);
      return `rgba(${c.r},${c.g},${c.b},${alpha})`;
    }

    function drawSphereBase(rect) {
      const radius = sphereRadius();
      const cx = rect.width / 2 + state.panX;
      const cy = rect.height * 0.49 + state.panY;
      const glow = ctx.createRadialGradient(cx - radius * 0.28, cy - radius * 0.32, radius * 0.08, cx, cy, radius * 1.05);
      glow.addColorStop(0, 'rgba(255,255,255,0.95)');
      glow.addColorStop(0.42, 'rgba(224,239,246,0.58)');
      glow.addColorStop(1, 'rgba(214,224,235,0.26)');
      ctx.save();
      ctx.beginPath();
      ctx.arc(cx, cy, radius, 0, Math.PI * 2);
      ctx.fillStyle = glow;
      ctx.fill();
      ctx.lineWidth = 1.1;
      ctx.strokeStyle = 'rgba(75,91,110,0.22)';
      ctx.stroke();
      ctx.globalAlpha = 0.36;
      ctx.strokeStyle = 'rgba(77,126,148,0.18)';
      for (let i = -3; i <= 3; i++) {
        ctx.beginPath();
        ctx.ellipse(cx, cy, radius * Math.sqrt(Math.max(0, 1 - i * i / 16)), radius * 0.18, 0, 0, Math.PI * 2);
        ctx.stroke();
      }
      for (let i = 0; i < 6; i++) {
        ctx.beginPath();
        ctx.ellipse(cx, cy, radius, radius * 0.18, i * Math.PI / 6, 0, Math.PI * 2);
        ctx.stroke();
      }
      ctx.restore();
    }

    function drawEdge(edge) {
      const source = nodesById.get(edge.source);
      const target = nodesById.get(edge.target);
      if (!source || !target) return;
      const a = screen(source);
      const b = screen(target);
      const rect = canvas.getBoundingClientRect();
      const cx = rect.width / 2 + state.panX;
      const cy = rect.height / 2 + state.panY;
      const midx = (a.x + b.x) / 2;
      const midy = (a.y + b.y) / 2;
      const curve = 0.18 + Math.min(0.28, Math.log1p(edge.count) / 32);
      const c1x = midx + (midx - cx) * curve;
      const c1y = midy + (midy - cy) * curve;
      const avgDepth = (a.depth + b.depth) / 2;
      const width = Math.min(5.5, 0.55 + Math.log1p(edge.count) * 0.42) * (0.72 + avgDepth * 0.36);
      ctx.save();
      ctx.globalAlpha = state.selected === edge.id || state.hover === edge.id ? 0.9 : 0.13 + avgDepth * 0.28;
      ctx.strokeStyle = edge.color;
      ctx.lineWidth = width;
      ctx.beginPath();
      ctx.moveTo(a.x, a.y);
      ctx.quadraticCurveTo(c1x, c1y, b.x, b.y);
      ctx.stroke();
      ctx.restore();
    }

    function drawNode(node) {
      const p = screen(node);
      const query = state.search.trim().toLowerCase();
      const matched = query && node.name.toLowerCase().includes(query);
      const active = state.selected === node.id || state.hover === node.id || matched;
      ctx.save();
      ctx.beginPath();
      ctx.arc(p.x, p.y, p.r, 0, Math.PI * 2);
      ctx.fillStyle = node.color;
      ctx.globalAlpha = active ? 1 : 0.46 + p.depth * 0.48;
      ctx.fill();
      ctx.lineWidth = active ? 3 : 1.1;
      ctx.strokeStyle = active ? '#1f2937' : `rgba(255,255,255,${0.55 + p.depth * 0.35})`;
      ctx.stroke();
      if (active || (node.count > 1600 && p.depth > 0.68) || (p.depth > 0.76 && node.label !== 'Variant' && node.count > 650)) {
        ctx.font = `${active ? 13 : 11}px Segoe UI, Arial, sans-serif`;
        ctx.fillStyle = active ? '#111827' : 'rgba(31,41,55,0.82)';
        ctx.globalAlpha = active ? 1 : 0.52 + p.depth * 0.42;
        ctx.fillText(shortName(node.name, active ? 42 : 24), p.x + p.r + 4, p.y + 4);
      }
      ctx.restore();
    }

    function draw() {
      const rect = canvas.getBoundingClientRect();
      ctx.clearRect(0, 0, rect.width, rect.height);
      const edges = visibleEdges();
      const nodes = visibleNodes(edges);
      const nodeIds = new Set(nodes.map(node => node.id));
      const filteredEdges = edges.filter(edge => nodeIds.has(edge.source) && nodeIds.has(edge.target));
      drawSphereBase(rect);
      filteredEdges
        .slice()
        .sort((a, b) => {
          const az = (screen(nodesById.get(a.source)).z + screen(nodesById.get(a.target)).z) / 2;
          const bz = (screen(nodesById.get(b.source)).z + screen(nodesById.get(b.target)).z) / 2;
          return az - bz;
        })
        .forEach(drawEdge);
      nodes
        .slice()
        .sort((a, b) => screen(a).z - screen(b).z)
        .forEach(drawNode);
      updateRelationshipTable(filteredEdges.slice().sort((a, b) => b.count - a.count).slice(0, 80));
    }

    function updateRelationshipTable(edges) {
      const tbody = document.querySelector('#relationshipTable tbody');
      tbody.innerHTML = edges.slice(0, 26).map(edge => {
        const source = nodesById.get(edge.source);
        const target = nodesById.get(edge.target);
        const high = edge.confidence_counts.High || 0;
        const med = edge.confidence_counts.Medium || 0;
        const low = edge.confidence_counts.Low || 0;
        return `
          <tr data-edge="${escapeHtml(edge.id)}">
            <td>${escapeHtml(shortName(source?.name, 18))} -> ${escapeHtml(shortName(target?.name, 18))}<br><span class="muted">${escapeHtml(edge.label)}</span></td>
            <td>${fmt(edge.count)}</td>
            <td><span class="pill">H ${fmt(high)}</span> <span class="pill">M ${fmt(med)}</span> <span class="pill">L ${fmt(low)}</span></td>
          </tr>
        `;
      }).join('');
      tbody.querySelectorAll('tr[data-edge]').forEach(row => {
        row.addEventListener('click', () => {
          state.selected = row.getAttribute('data-edge');
          renderDetail();
          draw();
        });
      });
    }

    function findNearest(x, y) {
      const edges = visibleEdges();
      const nodes = visibleNodes(edges);
      let best = null;
      let bestDist = Infinity;
      nodes.forEach(node => {
        const p = screen(node);
        const d = Math.hypot(x - p.x, y - p.y);
        if (d < p.r + 8 && d < bestDist) {
          best = node.id;
          bestDist = d;
        }
      });
      return best;
    }

    function renderDetail() {
      const detail = document.getElementById('detail');
      const selected = state.selected;
      if (!selected) {
        detail.innerHTML = '<div class="muted">No selection</div>';
        return;
      }
      const node = nodesById.get(selected);
      if (node) {
        detail.innerHTML = `
          <div class="detail-title">${escapeHtml(node.name)}</div>
          <span class="pill">${escapeHtml(node.label)}</span>
          <div class="kv">
            <div>Evidence</div><div>${fmt(node.count)}</div>
            <div>Dominant</div><div>${escapeHtml(node.dominant_confidence)}</div>
            <div>High</div><div>${fmt(node.confidence_counts.High || 0)}</div>
            <div>Medium</div><div>${fmt(node.confidence_counts.Medium || 0)}</div>
            <div>Low</div><div>${fmt(node.confidence_counts.Low || 0)}</div>
            <div>Tags</div><div>${escapeHtml(node.tags.join(', ') || 'none')}</div>
          </div>
          <div class="muted">Sample evidence: ${escapeHtml(node.sample_evidence.join(', ') || 'none')}</div>
        `;
        return;
      }
      const edge = kg.edges.find(item => item.id === selected);
      if (edge) {
        const source = nodesById.get(edge.source);
        const target = nodesById.get(edge.target);
        detail.innerHTML = `
          <div class="detail-title">${escapeHtml(shortName(source?.name, 36))} -> ${escapeHtml(shortName(target?.name, 36))}</div>
          <span class="pill">${escapeHtml(edge.label)}</span>
          <div class="kv">
            <div>Evidence</div><div>${fmt(edge.count)}</div>
            <div>Dominant</div><div>${escapeHtml(edge.dominant_confidence)}</div>
            <div>High</div><div>${fmt(edge.confidence_counts.High || 0)}</div>
            <div>Medium</div><div>${fmt(edge.confidence_counts.Medium || 0)}</div>
            <div>Low</div><div>${fmt(edge.confidence_counts.Low || 0)}</div>
            <div>Urgent</div><div>${fmt(edge.audit_counts.urgent_review || 0)}</div>
          </div>
          <div class="samples">
            ${edge.samples.map(sample => `
              <div class="sample">
                <strong>${escapeHtml(sample.evidence_id)}</strong>
                <div>${escapeHtml(sample.title)}</div>
                <div class="muted">${escapeHtml(sample.confidence_score)} - ${escapeHtml(sample.audit_priority)}</div>
                <div class="muted">${escapeHtml(sample.source_file)}</div>
                <div class="muted">DOI ${escapeHtml(sample.doi || 'not_reported')} | PMID ${escapeHtml(sample.pmid || 'not_reported')}</div>
              </div>
            `).join('')}
          </div>
        `;
      }
    }

    function renderLegend() {
      const colors = {
        Cancer: '#2f6f73',
        Gene: '#6d5bd0',
        Variant: '#9a6a2f',
        Assay: '#2d6cdf',
        MRDSignal: '#b33f62',
        ClinicalOutcome: '#3f7d3f',
        GeneAssociationRule: '#c23b7a',
        CancerGeneBaseline: '#138a72',
        LocusFeature: '#8b5a2b'
      };
      document.getElementById('legend').innerHTML = Object.entries(colors).map(([label, color]) => `
        <span class="legend-item"><span class="dot" style="background:${color};"></span>${label}</span>
      `).join('');
    }

    function setInitialViewport() {
      state.zoom = 1;
      state.panX = 0;
      state.panY = 0;
      state.rotX = -0.22;
      state.rotY = 0.62;
      state.autoRotate = true;
    }

    document.getElementById('viewSelect').addEventListener('change', event => {
      state.view = event.target.value;
      draw();
    });
    document.getElementById('minCount').addEventListener('input', event => {
      state.minCount = Number(event.target.value);
      document.getElementById('minCountLabel').textContent = state.minCount;
      draw();
    });
    document.getElementById('searchBox').addEventListener('input', event => {
      state.search = event.target.value;
      draw();
    });
    document.querySelectorAll('#confidenceChecks input').forEach(input => {
      input.addEventListener('change', () => {
        state.confidence = new Set(Array.from(document.querySelectorAll('#confidenceChecks input:checked')).map(item => item.value));
        draw();
      });
    });
    document.getElementById('resetView').addEventListener('click', () => {
      setInitialViewport();
      state.selected = null;
      renderDetail();
      draw();
    });

    canvas.addEventListener('mousedown', event => {
      state.dragging = true;
      state.lastX = event.offsetX;
      state.lastY = event.offsetY;
      canvas.classList.add('dragging');
    });
    window.addEventListener('mouseup', () => {
      state.dragging = false;
      canvas.classList.remove('dragging');
    });
    canvas.addEventListener('mousemove', event => {
      if (state.dragging) {
        const dx = event.offsetX - state.lastX;
        const dy = event.offsetY - state.lastY;
        if (event.shiftKey) {
          state.panX += dx;
          state.panY += dy;
        } else {
          state.rotY += dx * 0.008;
          state.rotX += dy * 0.008;
          state.rotX = Math.max(-1.35, Math.min(1.35, state.rotX));
        }
        state.lastX = event.offsetX;
        state.lastY = event.offsetY;
        state.autoRotate = false;
        draw();
        return;
      }
      state.hover = findNearest(event.offsetX, event.offsetY);
      draw();
    });
    canvas.addEventListener('click', event => {
      const nearest = findNearest(event.offsetX, event.offsetY);
      if (nearest) {
        state.selected = nearest;
        renderDetail();
        draw();
      }
    });
    canvas.addEventListener('wheel', event => {
      event.preventDefault();
      const delta = event.deltaY < 0 ? 1.1 : 0.9;
      state.zoom = Math.max(0.25, Math.min(2.4, state.zoom * delta));
      draw();
    }, { passive: false });

    function animate() {
      if (state.autoRotate && !state.dragging && !state.search.trim()) {
        state.rotY += 0.0022;
        draw();
      }
      requestAnimationFrame(animate);
    }

    renderStats();
    renderBars('confidenceBars', kg.summary.confidence, 6, '#2d6cdf');
    renderBars('evidenceBars', kg.summary.evidence_type, 8, '#6d5bd0');
    renderBars('geneBars', kg.summary.top_genes, 12, '#6d5bd0');
    renderBars('cancerBars', kg.summary.top_cancers, 10, '#2f6f73');
    renderBars('assayBars', kg.summary.top_assays, 10, '#2d6cdf');
    renderLegend();
    setInitialViewport();
    resizeCanvas();
    animate();
    window.addEventListener('resize', resizeCanvas);
  </script>
</body>
</html>
"""


def write_html(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    html_text = HTML_TEMPLATE.replace("__DATA_JSON__", js_data(data))
    path.write_text(html_text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Build an MRD KG visualization HTML.")
    parser.add_argument("--root", default=".")
    parser.add_argument("--raw-csv", default="data/evidence/mrd_evidence_raw_full.csv")
    parser.add_argument("--scored-csv", default="data/evidence/mrd_evidence_scored_full.csv")
    parser.add_argument("--sources-csv", default="data/imported/full_corpus_sources.csv")
    parser.add_argument("--chunks-csv", default="data/imported/full_corpus_excel_chunks.csv")
    parser.add_argument("--nodes-csv", default="data/graph/nodes_full.csv")
    parser.add_argument("--relationships-csv", default="data/graph/relationships_full.csv")
    parser.add_argument("--review-csv", default="data/graph/review_tasks_full.csv")
    parser.add_argument("--rag-jsonl", default="data/rag/evidence_summaries_full.jsonl")
    parser.add_argument("--baseline-csv", default="data/qc_rules/literature_cancer_gene_baselines.csv")
    parser.add_argument("--cooccurrence-csv", default="data/qc_rules/literature_gene_cooccurrence_rules.csv")
    parser.add_argument("--mutex-csv", default="data/qc_rules/literature_gene_mutex_rules.csv")
    parser.add_argument("--out-json", default="data/graph/mrd_kg_visualization_data.json")
    parser.add_argument("--out-html", default="docs/graph/mrd_kg_visualization.html")
    args = parser.parse_args()

    root = artifact_root() if args.root == "." else Path(args.root)
    data = build_visualization_data(args)
    out_json = artifact_path(root / args.out_json)
    out_html = artifact_path(root / args.out_html)
    write_json(out_json, data)
    write_html(out_html, data)
    print(f"Wrote visualization data: {out_json.resolve()} ({len(data['nodes'])} nodes, {len(data['edges'])} edges)")
    print(f"Wrote visualization HTML: {out_html.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
