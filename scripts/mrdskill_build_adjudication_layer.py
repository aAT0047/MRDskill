#!/usr/bin/env python3
"""Add MRD call adjudication templates and multi-agent RAG protocol nodes."""

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


AGENTS = [
    {
        "id": "retriever",
        "name": "Retriever agent",
        "role": "Retrieve KG, SQLite rule, and RAG evidence with paper/source/row traceability.",
        "checks": ["supporting_evidence", "source_traceability", "rule_support"],
    },
    {
        "id": "variant_origin",
        "name": "Variant-origin agent",
        "role": "Assess tumor-derived, CHIP/hematopoietic, germline, and artifact explanations for variant-based calls.",
        "checks": ["chip_risk", "matched_wbc", "tumor_tissue_consistency", "vaf_context", "cancer_baseline"],
    },
    {
        "id": "assay_context",
        "name": "Assay-context agent",
        "role": "Evaluate tumor-informed/tumor-naive design, LOD, panel size, sample type, and post-treatment sampling window.",
        "checks": ["assay_lod", "panel_coverage", "sample_type", "timepoint", "technical_no_call"],
    },
    {
        "id": "clinical_context",
        "name": "Clinical-context agent",
        "role": "Evaluate consistency with cancer type, stage, treatment setting, recurrence pattern, and MRD-guided treatment evidence.",
        "checks": ["cancer_consistency", "treatment_setting", "recurrence_pattern", "actionability_evidence"],
    },
    {
        "id": "contradiction_uncertainty",
        "name": "Contradiction and uncertainty agent",
        "role": "Actively search for false-positive, false-negative, assay limitation, and conflicting-literature evidence.",
        "checks": ["contradictory_evidence", "uncertainty_flags", "false_positive_risk", "false_negative_risk"],
    },
    {
        "id": "adjudicator",
        "name": "Adjudicator agent",
        "role": "Integrate support and counter-evidence into an auditable case-level MRD interpretation category and expert review checklist.",
        "checks": ["decision_category", "review_priority", "missing_information", "expert_checklist"],
    },
]


DECISIONS = [
    {
        "id": "actionable_positive",
        "name": "actionable_positive",
        "meaning": "Positive call has coherent molecular, assay, clinical, and literature support for expert-reviewed MRD actionability.",
        "review_action": "expert sign-out review; consider management options supported by cancer/stage-specific evidence",
    },
    {
        "id": "likely_positive_but_repeat",
        "name": "likely_positive_but_repeat",
        "meaning": "Positive call is plausible but confirmation is recommended because signal is low, single timepoint, early window, or context-limited.",
        "review_action": "repeat plasma sampling and review longitudinal trend before escalation",
    },
    {
        "id": "indeterminate",
        "name": "indeterminate",
        "meaning": "Evidence is insufficient or conflicting; neither true MRD positivity nor false positivity can be assigned with confidence.",
        "review_action": "manual review; obtain missing clinical, WBC, tissue, timepoint, or assay information",
    },
    {
        "id": "likely_false_positive",
        "name": "likely_false_positive",
        "meaning": "Call is better explained by technical artifact, germline, non-specific background, or inconsistent cancer/assay context.",
        "review_action": "do not use as standalone MRD-positive evidence; confirm by orthogonal evidence",
    },
    {
        "id": "likely_chip_or_hematopoietic",
        "name": "likely_CHIP_or_hematopoietic",
        "meaning": "Variant pattern, gene context, WBC evidence, VAF behavior, or literature support suggests hematopoietic origin.",
        "review_action": "confirm with matched WBC/buffy coat when available; downgrade MRD actionability",
    },
    {
        "id": "technical_no_call",
        "name": "technical_no_call",
        "meaning": "Assay, sample, coverage, LOD, or QC limitations prevent reliable MRD adjudication.",
        "review_action": "repeat assay or obtain technically valid sample",
    },
    {
        "id": "negative_but_high_false_negative_risk",
        "name": "negative_but_high_false_negative_risk",
        "meaning": "Negative result should not be over-interpreted because of low shedding, early timepoint, inadequate panel/LOD, or high-risk clinical context.",
        "review_action": "continue surveillance; consider repeat sampling or alternative assay",
    },
]


FLAGS = [
    ("transient_signal_possible", "Single or declining positive signal may be transient and requires longitudinal confirmation."),
    ("early_postoperative_window", "Sampling near surgery may reflect cfDNA background perturbation rather than stable residual disease."),
    ("chip_false_positive_risk", "Gene/VAF/WBC context suggests clonal hematopoiesis or hematopoietic origin."),
    ("germline_possible", "Call may require germline filtering or matched normal review."),
    ("technical_artifact_possible", "Call may reflect sequencing, alignment, panel, or low-support artifact."),
    ("non_specific_cfDNA_background", "Inflammation, necrosis, immunotherapy, radiation, or tissue injury may create non-specific cfDNA signal."),
    ("cancer_origin_uncertain", "Variant or signal is not specific enough to localize cancer or lesion origin."),
    ("recurrence_pattern_not_assessable", "Local versus distant recurrence risk cannot be inferred without clinical/imaging evidence."),
    ("serial_confirmation_needed", "Actionability depends on repeat positivity or longitudinal increase."),
    ("deescalation_not_proven", "Negative MRD alone is insufficient for de-escalation without setting-specific evidence."),
    ("low_shedding_false_negative_risk", "Cancer type/stage/sample context may increase false-negative risk."),
    ("assay_lod_or_panel_limitation", "LOD, panel coverage, tumor-informed status, or sample input limits interpretation."),
]


QUESTION_TEMPLATES = [
    {
        "id": "transient_signal_query",
        "question": "Is this MRD-positive call a transient signal or a persistent molecular recurrence signal?",
        "required_inputs": ["sample_id", "timepoint", "variant_calls", "VAF", "serial_samples"],
        "agents": ["retriever", "variant_origin", "clinical_context", "contradiction_uncertainty", "adjudicator"],
        "flags": ["transient_signal_possible", "serial_confirmation_needed"],
        "outputs": ["longitudinal_interpretation", "repeat_sampling_recommended", "supporting_evidence", "contradictory_evidence"],
    },
    {
        "id": "early_postoperative_noise_query",
        "question": "Could this positive call reflect early postoperative cfDNA background perturbation?",
        "required_inputs": ["surgery_date", "blood_draw_date", "assay_type", "cfDNA_quantity", "variant_calls"],
        "agents": ["retriever", "assay_context", "clinical_context", "contradiction_uncertainty", "adjudicator"],
        "flags": ["early_postoperative_window", "non_specific_cfDNA_background"],
        "outputs": ["sampling_window_flag", "repeat_sampling_recommended", "missing_information", "evidence_level"],
    },
    {
        "id": "chip_false_positive_query",
        "question": "Is this variant-based MRD call more consistent with tumor-derived signal or CHIP/hematopoietic origin?",
        "required_inputs": ["gene", "variant", "VAF", "matched_WBC_status", "tumor_tissue_status", "patient_age"],
        "agents": ["retriever", "variant_origin", "contradiction_uncertainty", "adjudicator"],
        "flags": ["chip_false_positive_risk", "germline_possible", "technical_artifact_possible"],
        "outputs": ["variant_origin_assessment", "WBC_confirmation_needed", "MRD_actionability_downgrade", "supporting_evidence"],
    },
    {
        "id": "non_specific_cfDNA_background_query",
        "question": "Could inflammation, necrosis, radiotherapy, or immunotherapy create a non-specific cfDNA signal?",
        "required_inputs": ["recent_treatment", "infection_or_inflammation", "cfDNA_quantity", "fragmentomics_score", "variant_calls"],
        "agents": ["retriever", "assay_context", "clinical_context", "contradiction_uncertainty", "adjudicator"],
        "flags": ["non_specific_cfDNA_background", "technical_artifact_possible"],
        "outputs": ["background_noise_risk", "needed_clinical_context", "repeat_or_orthogonal_test_recommendation"],
    },
    {
        "id": "cancer_origin_consistency_query",
        "question": "Is the positive molecular signal consistent with the reported cancer type or another possible origin?",
        "required_inputs": ["cancer_type", "variant_calls", "gene_set", "assay_type", "known_primary_tumor"],
        "agents": ["retriever", "variant_origin", "clinical_context", "contradiction_uncertainty", "adjudicator"],
        "flags": ["cancer_origin_uncertain", "chip_false_positive_risk"],
        "outputs": ["cancer_consistency", "alternate_origin_flags", "supporting_baseline_rules", "contradictory_evidence"],
    },
    {
        "id": "recurrence_pattern_query",
        "question": "Does the MRD evidence support local recurrence risk, distant/systemic recurrence risk, or is site-specific inference not assessable?",
        "required_inputs": ["cancer_type", "stage", "treatment_setting", "MRD_status", "imaging", "recurrence_site_if_known"],
        "agents": ["retriever", "clinical_context", "contradiction_uncertainty", "adjudicator"],
        "flags": ["recurrence_pattern_not_assessable"],
        "outputs": ["recurrence_pattern_inference", "site_specific_evidence", "missing_information", "evidence_level"],
    },
    {
        "id": "single_vs_serial_positive_query",
        "question": "Is a single positive MRD result sufficient, or is serial confirmation required in this setting?",
        "required_inputs": ["cancer_type", "stage", "treatment_setting", "single_or_serial_result", "assay_type"],
        "agents": ["retriever", "assay_context", "clinical_context", "contradiction_uncertainty", "adjudicator"],
        "flags": ["serial_confirmation_needed", "transient_signal_possible"],
        "outputs": ["single_positive_actionability", "serial_confirmation_needed", "trial_or_guideline_support"],
    },
    {
        "id": "deescalation_safety_query",
        "question": "Does a negative MRD result safely support treatment de-escalation in this cancer/stage/assay context?",
        "required_inputs": ["cancer_type", "stage", "treatment_setting", "MRD_negative_status", "assay_LOD", "followup_duration"],
        "agents": ["retriever", "assay_context", "clinical_context", "contradiction_uncertainty", "adjudicator"],
        "flags": ["deescalation_not_proven", "low_shedding_false_negative_risk", "assay_lod_or_panel_limitation"],
        "outputs": ["deescalation_support_level", "false_negative_risk", "required_followup", "missing_information"],
    },
    {
        "id": "mrd_guided_treatment_evidence_query",
        "question": "Does this cancer type and stage have evidence that MRD-guided treatment escalation, de-escalation, or surveillance improves outcomes?",
        "required_inputs": ["cancer_type", "stage", "treatment_setting", "MRD_status", "candidate_management_options"],
        "agents": ["retriever", "clinical_context", "contradiction_uncertainty", "adjudicator"],
        "flags": ["deescalation_not_proven", "serial_confirmation_needed"],
        "outputs": ["management_options_supported_by_evidence", "evidence_level", "trial_support", "expert_review_recommendation"],
    },
]


WORKFLOW_ROUNDS = [
    ("support_retrieval", "Retrieve supporting literature, KG edges, baseline rules, and source-row evidence for the MRD call."),
    ("counter_evidence_retrieval", "Actively retrieve CHIP, germline, artifact, assay-limit, false-negative, and contradictory evidence."),
    ("context_reassessment", "Reassess the call using assay, timepoint, sample type, cancer context, and longitudinal data."),
    ("case_adjudication", "Map the case to an adjudication category with confidence, missing information, and review priority."),
    ("expert_checklist_generation", "Generate an expert review checklist with traceable support and counter-evidence."),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Add MRD-KGQC adjudication framework nodes and RAG chunks.")
    parser.add_argument("--base-nodes", default="data/graph/nodes_augmented.csv")
    parser.add_argument("--base-relationships", default="data/graph/relationships_augmented.csv")
    parser.add_argument("--out-nodes-extension", default="data/graph/adjudication_nodes.csv")
    parser.add_argument("--out-relationships-extension", default="data/graph/adjudication_relationships.csv")
    parser.add_argument("--out-nodes-augmented", default="data/graph/nodes_augmented.csv")
    parser.add_argument("--out-relationships-augmented", default="data/graph/relationships_augmented.csv")
    parser.add_argument("--rag-jsonl", default="data/rag/kg_rag_chunks.jsonl")
    parser.add_argument("--out-template-json", default="data/rag/mrd_adjudication_templates.json")
    parser.add_argument("--out-report-md", default="docs/rag/mrd_adjudication_layer_usage.md")
    return parser.parse_args()


def clean(value: Any) -> str:
    return "" if value is None else str(value).strip()


def hash_id(prefix: str, *parts: Any) -> str:
    text = "|".join(clean(part) for part in parts)
    return f"{prefix}:{hashlib.sha1(text.encode('utf-8', 'ignore')).hexdigest()[:16]}"


def props(**kwargs: Any) -> str:
    return json.dumps({k: v for k, v in kwargs.items() if clean(v) != ""}, ensure_ascii=False, sort_keys=True)


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


def rel(start: str, end: str, rtype: str, properties: str = "{}", confidence: str = "High", audit: str = "no_review_needed") -> dict[str, Any]:
    return {
        "relationship_id": hash_id("REL", start, end, rtype, properties),
        "start_node_id": start,
        "end_node_id": end,
        "relationship_type": rtype,
        "confidence_score": confidence,
        "evidence_strength_score": "framework",
        "evidence_ids": "",
        "source_pmids": "",
        "source_dois": "",
        "graph_inclusion_status": "include_as_adjudication_framework",
        "audit_priority": audit,
        "properties_json": properties,
    }


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def upsert_csv(base: Path, extension: Path, out: Path, key_field: str) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    temp = out.with_suffix(out.suffix + ".tmp")
    with extension.open("r", encoding="utf-8", newline="") as ext_handle:
        extension_reader = csv.DictReader(ext_handle)
        fieldnames = extension_reader.fieldnames or []
        extension_rows = list(extension_reader)
    extension_keys = {clean(row.get(key_field)) for row in extension_rows if clean(row.get(key_field))}

    base_rows: list[dict[str, Any]] = []
    if base.exists():
        with base.open("r", encoding="utf-8", newline="") as base_handle:
            base_reader = csv.DictReader(base_handle)
            fieldnames = list(base_reader.fieldnames or fieldnames)
            for row in base_reader:
                if clean(row.get(key_field)) not in extension_keys:
                    base_rows.append(row)

    with temp.open("w", encoding="utf-8", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in [*base_rows, *extension_rows]:
            writer.writerow({field: row.get(field, "") for field in fieldnames})
    shutil.move(str(temp), str(out))


def upsert_jsonl(path: Path, chunks: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    chunk_ids = {clean(chunk.get("chunk_id")) for chunk in chunks if clean(chunk.get("chunk_id"))}
    preserved_lines: list[str] = []
    if path.exists():
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped:
                    continue
                try:
                    payload = json.loads(stripped)
                except json.JSONDecodeError:
                    preserved_lines.append(line if line.endswith("\n") else line + "\n")
                    continue
                metadata = payload.get("metadata") if isinstance(payload, dict) else {}
                kind = clean(metadata.get("kind")) if isinstance(metadata, dict) else ""
                chunk_id = clean(payload.get("chunk_id")) if isinstance(payload, dict) else ""
                if chunk_id in chunk_ids or kind.startswith("adjudication_"):
                    continue
                preserved_lines.append(json.dumps(payload, ensure_ascii=False, sort_keys=True) + "\n")

    with temp.open("w", encoding="utf-8") as handle:
        for line in preserved_lines:
            handle.write(line)
        for chunk in chunks:
            handle.write(json.dumps(chunk, ensure_ascii=False, sort_keys=True) + "\n")
    shutil.move(str(temp), str(path))


def make_chunk(chunk_id: str, title: str, text: str, metadata: dict[str, Any]) -> dict[str, Any]:
    return {"chunk_id": chunk_id, "title": title, "text": text, "metadata": metadata}


def build_layer() -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    nodes: dict[str, dict[str, Any]] = {}
    rels: list[dict[str, Any]] = []
    chunks: list[dict[str, Any]] = []

    framework_id = "AdjudicationFramework:MRD-KGQC"
    nodes[framework_id] = node(
        framework_id,
        "AdjudicationFramework",
        "MRD-KGQC evidence adjudication framework",
        "mrd_call_adjudication_framework",
        "mrd_kgqc_adjudication_layer",
        props(
            purpose="Transform variant-based ctDNA MRD calls into auditable evidence-backed interpretation categories.",
            scope="expert decision support before MRD report sign-out; not autonomous treatment decision-making",
        ),
    )

    for agent in AGENTS:
        agent_id = f"AgentRole:{agent['id']}"
        nodes[agent_id] = node(agent_id, "AgentRole", agent["name"], "multi_agent_rag_role", "mrd_kgqc_adjudication_layer", props(**agent))
        rels.append(rel(framework_id, agent_id, "USES_AGENT_ROLE"))
        text = f"{agent['name']}: {agent['role']} Checks: {', '.join(agent['checks'])}."
        chunks.append(make_chunk(agent_id, agent["name"], text, {"kind": "adjudication_agent_role", **agent}))

    for decision in DECISIONS:
        decision_id = f"AdjudicationDecision:{decision['id']}"
        nodes[decision_id] = node(decision_id, "AdjudicationDecision", decision["name"], "mrd_call_decision_category", "mrd_kgqc_adjudication_layer", props(**decision))
        rels.append(rel(framework_id, decision_id, "CAN_OUTPUT_DECISION"))
        text = f"{decision['name']}: {decision['meaning']} Review action: {decision['review_action']}."
        chunks.append(make_chunk(decision_id, f"Decision category: {decision['name']}", text, {"kind": "adjudication_decision_category", **decision}))

    for flag_id, description in FLAGS:
        node_id = f"UncertaintyFlag:{flag_id}"
        nodes[node_id] = node(node_id, "UncertaintyFlag", flag_id, "mrd_adjudication_uncertainty_flag", "mrd_kgqc_adjudication_layer", props(description=description))
        rels.append(rel(framework_id, node_id, "USES_UNCERTAINTY_FLAG"))

    for round_index, (round_id, description) in enumerate(WORKFLOW_ROUNDS, start=1):
        node_id = f"RetrievalRound:{round_id}"
        nodes[node_id] = node(node_id, "RetrievalRound", round_id, "multi_round_rag_step", "mrd_kgqc_adjudication_layer", props(order=round_index, description=description))
        rels.append(rel(framework_id, node_id, "HAS_RETRIEVAL_ROUND", props(order=round_index)))

    for template in QUESTION_TEMPLATES:
        template_id = f"QCQuestionTemplate:{template['id']}"
        nodes[template_id] = node(template_id, "QCQuestionTemplate", template["question"], "mrd_qc_query_template", "mrd_kgqc_adjudication_layer", props(**template))
        rels.append(rel(framework_id, template_id, "HAS_QC_QUESTION_TEMPLATE"))
        for agent_id in template["agents"]:
            rels.append(rel(template_id, f"AgentRole:{agent_id}", "REQUIRES_AGENT_ROLE"))
        for flag_id in template["flags"]:
            rels.append(rel(template_id, f"UncertaintyFlag:{flag_id}", "CHECKS_UNCERTAINTY_FLAG"))
        text = (
            f"MRD-KGQC query template: {template['question']} "
            f"Required inputs: {', '.join(template['required_inputs'])}. "
            f"Agents: {', '.join(template['agents'])}. "
            f"Flags to check: {', '.join(template['flags'])}. "
            f"Expected outputs: {', '.join(template['outputs'])}. "
            "Answers must include supporting evidence, contradictory evidence, missing information, confidence, and recommended expert review action."
        )
        chunks.append(make_chunk(template_id, f"MRD-KGQC query: {template['id']}", text, {"kind": "adjudication_query_template", **template}))

    workflow_text = (
        "MRD-KGQC uses multi-agent, multi-round, contradiction-seeking RAG for variant-based ctDNA MRD call adjudication. "
        "Round 1 retrieves supporting evidence; round 2 retrieves CHIP, germline, artifact, assay-limit, and false-negative counter-evidence; "
        "round 3 reassesses assay, timepoint, sample type, cancer context, and longitudinal data; round 4 maps the case to an adjudication category; "
        "round 5 generates an expert checklist with source traceability. The workflow supports expert review before report sign-out and does not autonomously decide treatment."
    )
    chunks.append(
        make_chunk(
            "AdjudicationWorkflow:MRD-KGQC",
            "MRD-KGQC multi-agent adjudication workflow",
            workflow_text,
            {"kind": "adjudication_workflow_protocol", "rounds": WORKFLOW_ROUNDS, "decision_categories": [d["id"] for d in DECISIONS]},
        )
    )

    template_payload = {
        "framework_id": framework_id,
        "agents": AGENTS,
        "decisions": DECISIONS,
        "uncertainty_flags": [{"id": flag_id, "description": description} for flag_id, description in FLAGS],
        "question_templates": QUESTION_TEMPLATES,
        "workflow_rounds": [{"id": round_id, "description": description, "order": index} for index, (round_id, description) in enumerate(WORKFLOW_ROUNDS, start=1)],
    }
    return list(nodes.values()), rels, chunks, template_payload


def write_report(path: Path, stats: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = f"""# MRD-KGQC Adjudication Layer

## Purpose

This layer reframes the KG/RAG from passive literature retrieval to active MRD call adjudication support.

It supports variant-based ctDNA MRD review questions such as:

- Is the positive signal transient?
- Could it be early postoperative cfDNA background?
- Is it CHIP / hematopoietic / germline / artifact rather than tumor-derived?
- Is the call consistent with the reported cancer type and assay context?
- Is a single positive enough, or is serial confirmation needed?
- Is MRD-guided escalation/de-escalation supported in this cancer/stage?

## Generated Files

- Adjudication nodes: `data/graph/adjudication_nodes.csv`
- Adjudication relationships: `data/graph/adjudication_relationships.csv`
- Augmented graph: `data/graph/nodes_augmented.csv`, `data/graph/relationships_augmented.csv`
- Template JSON: `data/rag/mrd_adjudication_templates.json`
- RAG chunks upserted into: `data/rag/kg_rag_chunks.jsonl`

## Counts

- Agent roles: {stats['agents']}
- Query templates: {stats['question_templates']}
- Decision categories: {stats['decisions']}
- Uncertainty flags: {stats['flags']}
- Extension nodes: {stats['nodes']}
- Extension relationships: {stats['relationships']}
- RAG chunks upserted: {stats['chunks']}

## Use

Answers should include `answer`, `confidence`, `supporting_evidence`, `contradictory_evidence`, `missing_information`, and `recommended_review_action`.

Run a multi-round retrieval packet with:

```powershell
python scripts\\mrdskill_run_adjudication_rag.py `
  --query "CRC postoperative ctDNA TP53 low VAF positive, matched WBC unavailable, should this be actionable MRD?" `
  --cancer-type CRC `
  --gene TP53 `
  --vaf 0.05% `
  --matched-wbc-status unavailable `
  --out-md docs/rag/example_mrd_kgqc_adjudication_packet.md
```

This layer supports expert review before MRD report sign-out. It must not be used as autonomous treatment decision-making.
"""
    path.write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    nodes, rels, chunks, template_payload = build_layer()
    write_csv(artifact_path(args.out_nodes_extension), nodes, NODE_FIELDS)
    write_csv(artifact_path(args.out_relationships_extension), rels, REL_FIELDS)
    upsert_csv(artifact_path(args.base_nodes), artifact_path(args.out_nodes_extension), artifact_path(args.out_nodes_augmented), "node_id")
    upsert_csv(
        artifact_path(args.base_relationships),
        artifact_path(args.out_relationships_extension),
        artifact_path(args.out_relationships_augmented),
        "relationship_id",
    )
    upsert_jsonl(artifact_path(args.rag_jsonl), chunks)
    artifact_path(args.out_template_json).parent.mkdir(parents=True, exist_ok=True)
    artifact_path(args.out_template_json).write_text(json.dumps(template_payload, ensure_ascii=False, indent=2), encoding="utf-8")
    stats = {
        "agents": len(AGENTS),
        "question_templates": len(QUESTION_TEMPLATES),
        "decisions": len(DECISIONS),
        "flags": len(FLAGS),
        "nodes": len(nodes),
        "relationships": len(rels),
        "chunks": len(chunks),
    }
    write_report(artifact_path(args.out_report_md), stats)
    print(json.dumps(stats, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
