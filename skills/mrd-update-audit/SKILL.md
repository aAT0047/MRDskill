---
name: mrd-update-audit
description: Run incremental MRDskill evidence updates and archive-first rebuild audits for new papers, new supplements and existing KG/RAG layers. Use when Codex needs to update the evidence corpus while preserving provenance, old versions and coverage reports.
---

# MRD Update Audit

Use this skill when adding new literature or rebuilding an existing MRDskill evidence layer.

## Inputs

- New papers, new supplementary files or staged update manifests.
- Existing evidence tables, KG/RAG chunks, rule chunks and archive path.
- Update policy and coverage checks.

## Workflow

1. Archive the previous evidence tables, KG/RAG chunks, rule chunks and manifests before rebuilding.
2. Ingest new papers and supplementary files through the same extraction, confidence-rating and graph-readiness rules.
3. Preserve stable IDs for unchanged evidence, nodes and relationships when possible.
4. Record new, updated, rejected and failed records separately.
5. Run coverage audits for source files, evidence rows, chunks, graph IDs, rule chunks and retrieval readiness.
6. Output the latest KG/RAG layer and an audit report.

## Bundled Commands

- `scripts/update.py`: run an incremental ingest or a coordinated update.
- `scripts/rebuild.py`: run an archive-first KG/RAG rebuild.
- `scripts/audit_coverage.py`: write a compact coverage audit for the latest evidence and graph artifacts.
- `scripts/ingest_new_literature.py`: stage new literature through extraction and confidence scoring.
- `scripts/audit_coverage_full.py`: run the detailed coverage audit for a completed build.

## Outputs

- Latest KG/RAG layer.
- Update audit report.
- Archive of the previous version.

## QC Boundary

Keep the old version before rebuilding and record every update step. Do not silently overwrite evidence or remove failed records from the audit trail.
