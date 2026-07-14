---
name: supplementary-rule
description: Convert public supplementary cohort tables into sample matrices and materialized QC-support rules such as co-occurrence, mutual-exclusion and review-trigger rules. Use when Codex needs cohort-table-derived rule chunks for MRDskill.
---

# Supplementary Rule

Use this skill for structured supplementary tables after they have been harvested and parsed.

## Inputs

- Supplementary tables, table manifests, parser logs and sample identifiers.
- Sample-variant, sample-gene or sample-feature matrices.
- Rule thresholds and provenance fields.

## Workflow

1. Normalize table fields into sample, gene, variant, feature, cancer type and endpoint columns.
2. Build sample matrices where row and column provenance can be retained.
3. Compute co-occurrence, mutual-exclusion and review-trigger rules when support is sufficient.
4. Attach rule weights, statistical fields, source table, row/column provenance and audit status.
5. Downgrade or mark rules with weak sample support, missing fields or unclear outcome definitions.

## Bundled Commands

- `scripts/build_rules.py`: materialize sample matrices, baseline tables, and rule candidates.
- `scripts/augment_graph.py`: add provenance-linked rule nodes and relationships to the graph layer.

## Outputs

- Co-occurrence rules.
- Mutual-exclusion rules.
- Rule chunks for KG/RAG use.

## QC Boundary

Rules are QC support only. They are not standalone clinical truth and must not override matched WBC, tumor concordance, assay LOD or serial evidence.
