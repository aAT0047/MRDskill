---
name: mrd-graphrag
description: Run multi-round GraphRAG for MRD source tracing, retrieving support evidence, counter-evidence, assay constraints, missing context and review actions from evidence summaries, KG chunks and rule chunks. Use when Codex needs a source-aware evidence packet rather than a one-shot RAG answer.
---

# MRD GraphRAG

Use this skill after the evidence graph and rule chunks are available.

## Inputs

- Evidence summaries, KG/RAG chunks and rule chunks.
- Variant, gene, case query or clinical scenario.
- Optional case context, assay context and prior retrieval state.

## Workflow

1. Retrieve MRD support, recurrence, persistence, clearance and prognosis evidence.
2. Retrieve CHIP/WBC, germline, artifact and tumor-concordance counter-evidence.
3. Retrieve assay type, LOD, sample type and time-window constraints.
4. Retrieve serial-sampling and clonal-dynamics context.
5. Retrieve false-positive mechanisms and technical limitations.
6. Produce a structured evidence packet with support, counter-evidence, rules, missing context, uncertainty flags and review actions.

## Bundled Commands

- `scripts/build_index.py`: build the latest exact rule index and unified retrieval chunks.
- `scripts/query.py`: run a multi-round source-tracing retrieval packet.
- `scripts/build_unified_index.py`: build the unified KG and cohort-rule chunk index directly.
- `scripts/query_lexical.py`: run the lightweight local retriever without the full adjudication protocol.

## Outputs

- Source-aware evidence packet.
- Retrieval-round log.
- Literature hits, rule hits and missing-context checklist.

## QC Boundary

Retrieval output does not directly equal final source assignment. Contradictions and missing context must remain visible for review-aware adjudication.
