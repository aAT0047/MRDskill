---
name: review-aware-adjudication
description: Integrate caller QC records, KG hits, rule hits, GraphRAG evidence packets, missing-context flags and patient history into variant-level, timepoint-level, case-level and longitudinal MRD QC conclusions, including molecular evidence verdicts, review actions, evidence lanes and HTML review reports.
---

# Review-Aware Adjudication

Use this skill as the final MRDskill output layer.

## Inputs

- Caller QC records.
- KG hits, rule hits and GraphRAG evidence packets.
- Missing-context flags: matched WBC, tumor concordance, serial confirmation, assay LOD and clinical timepoint.
- Optional patient history for longitudinal monitoring.

## Workflow

1. Combine caller support, graph support, rule support, RAG support, counter-evidence and missing-context penalties.
2. Assign a provisional source state: tumor-origin, WBC/CHIP-origin, germline-like, artifact-like, transient or unresolved.
3. Apply the review-aware gate to decide whether the evidence is resolved, review-needed, repeat-required, contradiction-review, technical no-call or not assessable.
4. Aggregate variant states into timepoint-level and case-level MRD QC conclusions.
5. For serial cases, update the patient trajectory as persistent molecular recurrence, clearance, new positivity, transient signal, WBC/CHIP-only persistence, mixed-source persistence or not assessable.
6. Render machine-readable outputs and an HTML review report.

## Bundled Commands

- `scripts/build_layer.py`: create or refresh the adjudication nodes, relationships, templates, and retrieval chunks.
- `scripts/run_adjudication.py`: integrate caller, graph, rule, retrieval, and missing-context signals into QC conclusions.
- `scripts/render_report.py`: render a provenance-linked review report from QC outputs.
- `scripts/build_decision_report.py`: render a compact decision report from variant QC outputs.
- `scripts/collect_case_outputs.py`: gather available QC and review files into a case-specific artifact folder.

## Outputs

- Per-variant QC table.
- Timepoint MRD state.
- Case-level QC conclusion.
- Molecular evidence verdict.
- Recommended review actions.
- Missing-for-sign-out checklist.
- Evidence lane and provenance path.
- HTML review report.

## QC Boundary

Do not force tumor-origin calls when matched WBC, tumor concordance, serial confirmation, assay LOD or clinical timing is missing. HTML output supports expert review and does not replace clinical sign-out.
