#!/usr/bin/env python3
"""Run multi-round MRD-KGQC source-tracing retrieval over the local KG/RAG index."""

from __future__ import annotations

import argparse
import heapq
import json
import math
import re
from datetime import datetime
from pathlib import Path
from typing import Any

from mrdskill_project_paths import artifact_path
from mrdskill_query_mrd_kg_rag import existing_chunks_path, safe_print, search_rule_db, tokenize
from collections import Counter


_CHUNK_TOKEN_CACHE: dict[tuple[str, int], list[tuple[dict[str, Any], Counter[str], str]]] = {}
_ROUND_CANDIDATE_CACHE: dict[
    tuple[str, int, str, bool, tuple[str, ...]],
    list[tuple[dict[str, Any], Counter[str], str]],
] = {}


CORE_SOURCE_TRACE_PHRASES = [
    "ctdna",
    "cfdna",
    "cell-free dna",
    "cell free dna",
    "circulating tumor dna",
    "circulating tumour dna",
    "molecular residual disease",
    "minimal residual disease",
    "measurable residual disease",
    "mrd",
    "chip",
    "clonal hematopo",
    "clonal haematopo",
]

SOURCE_CONTEXT_PHRASES = [
    "assay",
    "artifact",
    "artefact",
    "buffy",
    "clearance",
    "germline",
    "hematopoietic",
    "haematopoietic",
    "limit of detection",
    "lod",
    "longitudinal",
    "low vaf",
    "plasma",
    "recurrence",
    "relapse",
    "serial",
    "tumor informed",
    "tumor-informed",
    "tumor tissue",
    "tumour tissue",
    "variant",
    "vaf",
    "wbc",
]

IRRELEVANT_MRD_ACRONYM_PHRASES = [
    "gaussian process",
    "kronecker product",
    "latent dimensions",
    "multi-output gaussian",
    "rbf kernel",
    "single-cell multi-omics",
]

LOW_VALUE_REVIEW_PHRASES = [
    "autologous hematopoietic stem cell transplantation",
    "car t cell",
    "multiple myeloma",
]


ROUNDS = [
    {
        "round_id": "support_retrieval",
        "agent": "retriever",
        "purpose": "Find source-traceable evidence for caller-derived MRD molecular signals.",
        "query_terms": "ctDNA MRD molecular residual disease caller variant source tracing recurrence relapse persistence clearance prognosis",
        "required_terms": ["ctdna", "cfdna", "circulating tumor dna", "circulating tumour dna", "molecular residual disease", "minimal residual disease"],
        "boost_terms": ["tumor informed", "tumor-informed", "serial", "recurrence", "relapse", "variant", "vaf"],
    },
    {
        "round_id": "variant_origin",
        "agent": "variant_origin",
        "purpose": "Trace possible tumor-derived, CHIP, hematopoietic, germline, or artifact origin.",
        "query_terms": "CHIP clonal hematopoiesis matched WBC buffy coat germline artifact low VAF tumor tissue",
        "required_terms": ["chip", "clonal hematopo", "clonal haematopo", "wbc", "buffy", "germline", "artifact", "artefact", "low vaf", "ctdna", "cfdna"],
        "boost_terms": ["matched wbc", "buffy coat", "hematopoietic", "haematopoietic", "germline", "artifact", "low vaf"],
    },
    {
        "round_id": "assay_timepoint_context",
        "agent": "assay_context",
        "purpose": "Retrieve assay design, LOD, sample type, and post-treatment timing evidence that affects source tracing.",
        "query_terms": "tumor informed tumor naive assay LOD limit of detection plasma postoperative surgery treatment timepoint cfDNA background",
        "required_terms": ["assay", "limit of detection", "lod", "tumor informed", "tumor-informed", "tumor naive", "tumor-naive", "ctdna", "cfdna", "plasma"],
        "boost_terms": ["limit of detection", "lod", "tumor informed", "tumor-informed", "postoperative", "plasma"],
    },
    {
        "round_id": "longitudinal_clonal_context",
        "agent": "longitudinal_context",
        "purpose": "Find evidence for serial confirmation, tumor-informed tracking, clonal evolution, recurrence origin, and treatment-selected clones.",
        "query_terms": "serial sampling longitudinal ctDNA trajectory tumor informed tracking clonal evolution recurrence origin resistant clone treatment selected variant",
        "required_terms": ["serial", "longitudinal", "trajectory", "persistence", "clearance", "recurrence", "relapse", "ctdna", "cfdna"],
        "boost_terms": ["serial", "longitudinal", "trajectory", "persistence", "clearance", "tumor informed", "recurrence"],
    },
    {
        "round_id": "targeted_resistance_context",
        "agent": "therapy_context",
        "purpose": "Retrieve treatment-related MRD evidence where specific caller alterations affect resistance monitoring, targeted therapy review, or trial eligibility.",
        "query_terms": "targeted therapy resistance trial eligibility actionable alteration EGFR KRAS NRAS BRAF MET ALK ROS1 RET ERBB2 ESR1 BRCA reversion IDH ctDNA MRD",
        "required_terms": ["ctdna", "cfdna", "variant", "resistance", "targeted therapy", "trial eligibility", "actionable", "mrd"],
        "boost_terms": ["resistance", "targeted therapy", "trial eligibility", "actionable", "egfr", "kras", "braf", "met", "brca"],
    },
    {
        "round_id": "contradiction_uncertainty",
        "agent": "contradiction_uncertainty",
        "purpose": "Retrieve alternative source explanations, uncertainty, confounding, low-shedding, and no-call conditions.",
        "query_terms": "artifact confounder CHIP hematopoietic germline no-call uncertainty inflammation necrosis radiotherapy immunotherapy non-specific cfDNA low shedding",
        "required_terms": ["artifact", "artefact", "confound", "chip", "hematopoietic", "haematopoietic", "germline", "low shedding", "cfdna", "ctdna"],
        "boost_terms": ["artifact", "confound", "chip", "hematopoietic", "germline", "low shedding", "non-specific"],
    },
    {
        "round_id": "source_trace_protocol",
        "agent": "source_trace_protocol",
        "purpose": "Retrieve MRD-KGQC source-tracing templates, uncertainty flags, and expert checklist requirements.",
        "query_terms": "MRD-KGQC source tracing caller signal tumor origin CHIP hematopoietic artifact serial confirmation expert checklist",
        "required_terms": ["source", "trace", "tracing", "caller", "chip", "artifact", "serial", "confirmation", "ctdna", "mrd"],
        "boost_terms": ["source", "tracing", "caller", "confirmation", "expert", "checklist", "chip", "artifact"],
    },
]


DECISION_CATEGORIES = [
    "traceable_caller_evidence_candidate",
    "traceable_caller_signal_needs_serial_confirmation",
    "hematopoietic_or_CHIP_origin_suspected",
    "technical_artifact_origin_suspected",
    "technical_rescue_required_before_source_assignment",
    "cross_sectional_signal_needs_contextualization",
    "untraceable_missing_annotation_or_context",
    "discordant_signal_origin_needs_reconciliation",
    "low_shedding_or_no_caller_signal_context_needed",
]


UNCERTAINTY_FLAGS = [
    "transient_signal_possible",
    "early_postoperative_window",
    "chip_or_hematopoietic_origin_risk",
    "germline_possible",
    "technical_artifact_possible",
    "non_specific_cfDNA_background",
    "cancer_origin_uncertain",
    "recurrence_pattern_not_assessable",
    "serial_confirmation_needed",
    "deescalation_not_proven",
    "low_shedding_false_negative_risk",
    "assay_lod_or_panel_limitation",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run multi-agent, multi-round MRD-KGQC source-tracing retrieval.")
    parser.add_argument("--query", default="", help="Free-text MRD case question.")
    parser.add_argument("--case-json", default="", help="Optional JSON file with case fields.")
    parser.add_argument("--cancer-type", default="")
    parser.add_argument("--stage", default="")
    parser.add_argument("--treatment-setting", default="")
    parser.add_argument("--assay-type", default="")
    parser.add_argument("--sample-type", default="")
    parser.add_argument("--timepoint", default="")
    parser.add_argument("--gene", default="")
    parser.add_argument("--variant", default="")
    parser.add_argument("--vaf", default="")
    parser.add_argument("--mrd-status", default="")
    parser.add_argument("--matched-wbc-status", default="")
    parser.add_argument("--tumor-tissue-status", default="")
    parser.add_argument("--chunks-jsonl", default="data/rag/mrd_latest_unified_rag_chunks.jsonl")
    parser.add_argument("--rule-db", default="data/rag/mrd_latest_kg_rule_index.sqlite")
    parser.add_argument("--top-k-per-round", type=int, default=5)
    parser.add_argument("--rule-top-k-per-round", type=int, default=3)
    parser.add_argument("--max-chunks", type=int, default=0, help="Optional cap for quick smoke tests; 0 scans all chunks.")
    parser.add_argument("--min-score", type=float, default=0.0, help="Minimum denoised retrieval score to keep.")
    parser.add_argument("--no-denoise", action="store_true", help="Disable MRD source-tracing denoising filters.")
    parser.add_argument("--out-json", default="")
    parser.add_argument("--out-md", default="")
    return parser.parse_args()


def clean(value: Any) -> str:
    return "" if value is None else str(value).strip()


def is_missing_value(value: Any) -> bool:
    return clean(value).lower() in {
        "",
        ".",
        "na",
        "n/a",
        "none",
        "null",
        "unknown",
        "uncertain",
        "unavailable",
        "not_available",
        "needs_review",
        "not_provided",
        "not_reported",
        "not assessable",
        "not_assessable",
    }


def load_case_payload(args: argparse.Namespace) -> dict[str, str]:
    payload: dict[str, str] = {}
    if args.case_json:
        path = Path(args.case_json)
        if not path.exists():
            path = artifact_path(args.case_json)
        if path.exists():
            raw = json.loads(path.read_text(encoding="utf-8"))
            if isinstance(raw, dict):
                payload.update({str(key): clean(value) for key, value in raw.items()})
    for key in [
        "query",
        "cancer_type",
        "stage",
        "treatment_setting",
        "assay_type",
        "sample_type",
        "timepoint",
        "gene",
        "variant",
        "vaf",
        "mrd_status",
        "matched_wbc_status",
        "tumor_tissue_status",
    ]:
        value = clean(getattr(args, key, ""))
        if value:
            payload[key] = value
    return payload


def case_text(payload: dict[str, str]) -> str:
    labels = [
        "query",
        "cancer_type",
        "stage",
        "treatment_setting",
        "assay_type",
        "sample_type",
        "timepoint",
        "gene",
        "variant",
        "vaf",
        "mrd_status",
        "matched_wbc_status",
        "tumor_tissue_status",
    ]
    parts = []
    for label in labels:
        value = clean(payload.get(label))
        if value:
            parts.append(f"{label}: {value}")
    return " | ".join(parts)


def truncate(text: str, limit: int = 1200) -> str:
    text = clean(text)
    return text if len(text) <= limit else text[: limit - 3] + "..."


def text_blob(chunk: dict[str, Any]) -> str:
    return " ".join(
        [
            clean(chunk.get("title")),
            clean(chunk.get("text")),
            json.dumps(chunk.get("metadata", {}), ensure_ascii=False),
        ]
    )


def phrase_hits(text_lower: str, phrases: list[str]) -> list[str]:
    return [phrase for phrase in phrases if phrase.lower() in text_lower]


def metadata_score(metadata: dict[str, Any]) -> float:
    score = 0.0
    confidence = clean(metadata.get("confidence_score")).lower()
    inclusion = clean(metadata.get("graph_inclusion_status")).lower()
    audit = clean(metadata.get("audit_priority")).lower()
    if confidence == "high":
        score += 6
    elif confidence == "medium":
        score += 3
    elif confidence == "low":
        score -= 1
    if inclusion == "auto_include_high_confidence":
        score += 4
    elif inclusion == "include_as_medium_confidence":
        score += 2
    elif inclusion == "include_as_low_confidence_annotation":
        score += 0.5
    if audit in {"urgent_review", "high_review"}:
        score += 1
    return score


def hit_payload(score: float, chunk: dict[str, Any]) -> dict[str, Any]:
    metadata = chunk.get("metadata", {})
    if not isinstance(metadata, dict):
        metadata = {}
    return {
        "score": round(float(score), 4),
        "chunk_id": clean(chunk.get("chunk_id")),
        "title": clean(chunk.get("title")),
        "text": truncate(clean(chunk.get("text"))),
        "metadata": {
            "kind": clean(metadata.get("kind")),
            "evidence_id": clean(metadata.get("evidence_id")),
            "paper_id": clean(metadata.get("paper_id")),
            "doi": clean(metadata.get("doi") or metadata.get("supporting_dois")),
            "pmid": clean(metadata.get("pmid") or metadata.get("supporting_pmids")),
            "journal": clean(metadata.get("journal") or metadata.get("supporting_journals")),
            "confidence_score": clean(metadata.get("confidence_score")),
            "graph_inclusion_status": clean(metadata.get("graph_inclusion_status")),
            "audit_priority": clean(metadata.get("audit_priority")),
        },
    }


def iter_chunks(path: Path, max_chunks: int = 0):
    seen = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue
            try:
                chunk = json.loads(line)
            except json.JSONDecodeError:
                continue
            yield chunk
            seen += 1
            if max_chunks and seen >= max_chunks:
                break


def chunk_token_counts(chunk: dict[str, Any]) -> Counter[str]:
    text = text_blob(chunk)
    return Counter(tokenize(text))


def iter_chunk_search_records(path: Path, max_chunks: int = 0):
    cache_key = (str(path.resolve()), max_chunks)
    if cache_key not in _CHUNK_TOKEN_CACHE:
        records: list[tuple[dict[str, Any], Counter[str], str]] = []
        for chunk in iter_chunks(path, max_chunks):
            token_counts = chunk_token_counts(chunk)
            if token_counts:
                records.append((chunk, token_counts, text_blob(chunk).lower()))
        _CHUNK_TOKEN_CACHE[cache_key] = records
    yield from _CHUNK_TOKEN_CACHE[cache_key]


def passes_source_trace_filter(
    lower: str,
    spec: dict[str, Any],
    required_terms: list[str] | None = None,
    denoise: bool = True,
) -> bool:
    if not denoise:
        return True
    if phrase_hits(lower, IRRELEVANT_MRD_ACRONYM_PHRASES):
        return False
    required = list(required_terms or spec.get("required_terms") or CORE_SOURCE_TRACE_PHRASES)
    if required and not phrase_hits(lower, required):
        return False
    core_hits = phrase_hits(lower, CORE_SOURCE_TRACE_PHRASES)
    context_hits = phrase_hits(lower, SOURCE_CONTEXT_PHRASES)
    if not core_hits and len(context_hits) < 2:
        return False
    return True


def round_candidate_records(
    path: Path,
    max_chunks: int,
    spec: dict[str, Any],
    required_terms: list[str] | None,
    denoise: bool,
) -> list[tuple[dict[str, Any], Counter[str], str]]:
    required_key = tuple(sorted(term.lower() for term in (required_terms or spec.get("required_terms") or [])))
    cache_key = (str(path.resolve()), max_chunks, clean(spec.get("round_id")), denoise, required_key)
    if cache_key not in _ROUND_CANDIDATE_CACHE:
        _ROUND_CANDIDATE_CACHE[cache_key] = [
            record
            for record in iter_chunk_search_records(path, max_chunks)
            if passes_source_trace_filter(record[2], spec, required_terms, denoise)
        ]
    return _ROUND_CANDIDATE_CACHE[cache_key]


def score_counts(query_terms: Counter[str], token_counts: Counter[str]) -> float:
    score = 0.0
    for term, q_count in query_terms.items():
        tf = token_counts.get(term, 0)
        if tf:
            score += (1 + math.log(tf)) * q_count
    return score


def denoised_score(
    base_score: float,
    chunk: dict[str, Any],
    spec: dict[str, Any],
    required_terms: list[str] | None = None,
    boost_terms: list[str] | None = None,
    negative_terms: list[str] | None = None,
    denoise: bool = True,
    lower_blob: str | None = None,
) -> float:
    if base_score <= 0:
        return 0.0
    lower = lower_blob if lower_blob is not None else text_blob(chunk).lower()
    metadata = chunk.get("metadata", {})
    if not isinstance(metadata, dict):
        metadata = {}

    required = list(required_terms or spec.get("required_terms") or CORE_SOURCE_TRACE_PHRASES)
    boosts = list(boost_terms or []) + list(spec.get("boost_terms") or []) + SOURCE_CONTEXT_PHRASES
    negatives = list(negative_terms or []) + IRRELEVANT_MRD_ACRONYM_PHRASES

    if denoise:
        if phrase_hits(lower, IRRELEVANT_MRD_ACRONYM_PHRASES):
            return 0.0
        required_hits = phrase_hits(lower, required)
        core_hits = phrase_hits(lower, CORE_SOURCE_TRACE_PHRASES)
        context_hits = phrase_hits(lower, SOURCE_CONTEXT_PHRASES)
        if required and not required_hits:
            return 0.0
        if not core_hits and len(context_hits) < 2:
            return 0.0
        if phrase_hits(lower, LOW_VALUE_REVIEW_PHRASES) and not any(
            phrase in lower
            for phrase in ["ctdna", "cfdna", "variant", "vaf", "chip", "wbc", "germline", "artifact", "source"]
        ):
            base_score *= 0.25

    score = base_score
    score += len(phrase_hits(lower, boosts)) * 2.5
    score += metadata_score(metadata)
    score -= len(phrase_hits(lower, negatives)) * 8
    if re.search(r"\bmrd\b", lower, flags=re.I) and not any(
        phrase in lower
        for phrase in ["molecular residual disease", "minimal residual disease", "measurable residual disease", "ctdna", "cfdna"]
    ):
        score -= 12
    return max(0.0, score)


def retrieve_rounds_streaming(
    chunks_path: Path,
    base_text: str,
    top_k: int,
    max_chunks: int = 0,
    required_terms: list[str] | None = None,
    boost_terms: list[str] | None = None,
    negative_terms: list[str] | None = None,
    min_score: float = 0.0,
    denoise: bool = True,
) -> list[dict[str, Any]]:
    round_queries = [" ".join([base_text, spec["query_terms"]]).strip() for spec in ROUNDS]
    round_query_terms = [Counter(tokenize(query)) for query in round_queries]
    heaps: list[list[tuple[float, int, dict[str, Any]]]] = [[] for _ in ROUNDS]
    seq = 0
    all_records = list(iter_chunk_search_records(chunks_path, max_chunks))
    chunks_scanned = len(all_records)
    candidate_counts = [
        len(round_candidate_records(chunks_path, max_chunks, spec, required_terms, denoise))
        for spec in ROUNDS
    ]
    if top_k <= 0:
        return [
            {
                "round_id": spec["round_id"],
                "agent": spec["agent"],
                "purpose": spec["purpose"],
                "query": query,
                "chunks_scanned": 0,
                "hits": [],
                "retrieval_profile": {"denoise": denoise, "min_score": min_score},
            }
            for spec, query in zip(ROUNDS, round_queries)
        ]
    for index, (spec, query_terms) in enumerate(zip(ROUNDS, round_query_terms)):
        candidates = round_candidate_records(chunks_path, max_chunks, spec, required_terms, denoise)
        for chunk, token_counts, lower_blob in candidates:
            raw_score = score_counts(query_terms, token_counts)
            score = denoised_score(raw_score, chunk, spec, required_terms, boost_terms, negative_terms, denoise, lower_blob)
            if score <= 0:
                continue
            if min_score and score < min_score:
                continue
            seq += 1
            item = (float(score), seq, chunk)
            heap = heaps[index]
            if len(heap) < top_k:
                heapq.heappush(heap, item)
            elif score > heap[0][0]:
                heapq.heapreplace(heap, item)

    results: list[dict[str, Any]] = []
    for spec, query, heap in zip(ROUNDS, round_queries, heaps):
        ranked = sorted(heap, key=lambda item: item[0], reverse=True)
        results.append(
            {
                "round_id": spec["round_id"],
                "agent": spec["agent"],
                "purpose": spec["purpose"],
                "query": query,
                "chunks_scanned": chunks_scanned,
                "candidate_chunks_scored": candidate_counts[ROUNDS.index(spec)],
                "hits": [hit_payload(score, chunk) for score, _, chunk in ranked],
                "retrieval_profile": {"denoise": denoise, "min_score": min_score},
            }
        )
    return results


def checklist(payload: dict[str, str]) -> list[dict[str, str]]:
    fields = [
        ("matched_wbc_status", "Matched WBC/buffy coat status is needed to downgrade CHIP/hematopoietic origin."),
        ("tumor_tissue_status", "Tumor tissue confirmation helps separate tumor-derived MRD from plasma-only noise."),
        ("timepoint", "Sampling timepoint is needed to evaluate postoperative or post-treatment cfDNA perturbation."),
        ("vaf", "VAF or molecular burden is needed for LOD, low-level signal, and longitudinal trend review."),
        ("assay_type", "Assay type is needed for tumor-informed/tumor-naive and LOD interpretation."),
        ("stage", "Cancer stage is needed for MRD-guided management evidence and false-negative risk."),
        ("treatment_setting", "Treatment setting is needed before discussing escalation, de-escalation, or surveillance."),
    ]
    rows = []
    for field, reason in fields:
        rows.append(
            {
                "field": field,
                "status": "missing" if is_missing_value(payload.get(field)) else "provided",
                "reason": reason,
            }
        )
    return rows


def write_markdown(path: Path, result: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# MRD-KGQC Caller Signal Source-Tracing RAG Packet",
        "",
        f"- Generated: {result['generated_at']}",
        f"- Case: {result['case_text']}",
        "",
        "## Intended Use",
        "",
        "This packet supports caller-derived molecular evidence source tracing. It does not classify MRD positivity or make autonomous treatment decisions.",
        "",
        "## Candidate Decision Categories",
        "",
    ]
    lines.extend(f"- `{category}`" for category in DECISION_CATEGORIES)
    lines.extend(["", "## Uncertainty Flags To Check", ""])
    lines.extend(f"- `{flag}`" for flag in UNCERTAINTY_FLAGS)
    lines.extend(["", "## Missing-Information Checklist", ""])
    for row in result["missing_information_checklist"]:
        lines.append(f"- `{row['field']}`: {row['status']} - {row['reason']}")
    for round_result in result["rounds"]:
        lines.extend(["", f"## {round_result['round_id']} ({round_result['agent']})", "", round_result["purpose"], ""])
        lines.append(f"Query: `{round_result['query']}`")
        lines.append("")
        for index, hit in enumerate(round_result["hits"], start=1):
            meta = hit["metadata"]
            lines.append(f"{index}. score={hit['score']} `{hit['chunk_id']}`")
            lines.append(f"   Title: {hit['title']}")
            lines.append(f"   DOI: {meta.get('doi')} | PMID: {meta.get('pmid')} | confidence: {meta.get('confidence_score')}")
            lines.append(f"   Text: {truncate(hit['text'], 500)}")
        if round_result.get("rule_hits"):
            lines.extend(["", "Structured rule hits:", ""])
            for index, row in enumerate(round_result["rule_hits"], start=1):
                lines.append(
                    f"{index}. `{row.get('rule_id')}` {row.get('rule_type')} {row.get('cancer_type')} "
                    f"{row.get('feature_a_label') or row.get('gene_a') or row.get('feature_a_id')} / "
                    f"{row.get('feature_b_label') or row.get('gene_b') or row.get('feature_b_id')} "
                    f"pair={row.get('pair_count')}/{row.get('sample_count')}"
                )
    lines.extend(
        [
            "",
            "## Required Final Answer Shape",
            "",
            "- source_trace_answer",
            "- source_trace_category",
            "- confidence",
            "- supporting_evidence",
            "- contradictory_evidence",
            "- uncertainty_flags",
            "- missing_information",
            "- recommended_expert_review_action",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    payload = load_case_payload(args)
    base_text = case_text(payload)
    if not base_text:
        raise SystemExit("Provide --query, --case-json, or case fields such as --cancer-type/--gene/--variant.")
    chunks_path = existing_chunks_path(args.chunks_jsonl)
    rounds = retrieve_rounds_streaming(
        chunks_path,
        base_text,
        args.top_k_per_round,
        args.max_chunks,
        min_score=args.min_score,
        denoise=not args.no_denoise,
    )
    allowed_genes = [term for term in tokenize(clean(payload.get("gene"))) if re.fullmatch(r"[A-Z][A-Z0-9.-]{1,15}", term)]
    for round_result in rounds:
        round_result["rule_hits"] = search_rule_db(
            artifact_path(args.rule_db),
            round_result["query"],
            args.rule_top_k_per_round,
            allowed_genes=allowed_genes,
            require_allowed_gene=bool(allowed_genes),
        )
    result = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "chunks_jsonl": str(chunks_path),
        "rule_db": str(artifact_path(args.rule_db)),
        "case": payload,
        "case_text": base_text,
        "decision_categories": DECISION_CATEGORIES,
        "uncertainty_flags": UNCERTAINTY_FLAGS,
        "missing_information_checklist": checklist(payload),
        "rounds": rounds,
    }
    if args.out_json:
        out_json = artifact_path(args.out_json)
        out_json.parent.mkdir(parents=True, exist_ok=True)
        out_json.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    if args.out_md:
        write_markdown(artifact_path(args.out_md), result)
    safe_print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
