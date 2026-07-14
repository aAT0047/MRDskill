---
name: mrd-vcf-qc
description: Normalize caller-derived VCF records and case context into MRDskill variant QC records, including gene, VAF, depth, caller filters, matched-control signal, review flags and case reports. Use when Codex needs caller-aware MRD QC before GraphRAG or adjudication.
---

# MRD VCF QC

Use this skill for VCF-derived input before source-aware retrieval or adjudication.

## Inputs

- Caller-derived VCF files.
- Case metadata, sample metadata and optional matched-control fields.
- Gene annotation and caller configuration.

## Workflow

1. Parse variant identity: chromosome, position, reference allele, alternate allele and caller source.
2. Extract gene, VAF, depth, alternate-read support, quality, filter status and annotation fields.
3. Detect matched-control signals and caller-specific risk flags when available.
4. Produce variant-level QC records and case-level summaries.
5. Mark missing WBC, tumor tissue, serial sample, LOD or clinical timepoint context.

## Bundled Command

Run `scripts/run_qc.py` to normalize VCF and case-manifest inputs, retrieve graph and rule context, and write variant-level QC and case outputs.

## Outputs

- Variant QC records.
- Caller burden summary.
- Case report or QC manifest.

## QC Boundary

VCF records represent caller signal, not biological truth. Do not infer final tumor origin without source-aware evidence and context.
