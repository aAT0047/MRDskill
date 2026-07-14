---
name: mrd-confidence-rater
description: Use when the user wants to assign evidence confidence, evidence source type, audit status, and graph-readiness labels to structured MRD, ctDNA, CHIP, cancer-gene, variant, assay, recurrence, relapse, prognosis, or false-positive evidence records extracted from biomedical papers. This skill assumes evidence records have already been extracted by mrd-evidence-miner. It does not perform large-scale paper retrieval, raw evidence extraction, or knowledge graph construction.
allowed-tools:
  - readPaper
  - searchArticles
---

# MRD Confidence Rater

Use this skill to rate evidence confidence for MRD biological QC knowledge graph construction.

This skill is the third step in the MRD biological QC pipeline.

## Python Helpers

Prefer bundled Python scripts for deterministic first-pass scoring before manual or LLM-assisted review.

- `scripts/mrdskill_score_mrd_evidence.py`: add final source type, confidence score, evidence strength score, graph inclusion status, audit priority, traceability status, downgrade/upgrade reasons, and biological plausibility labels to `mrd_evidence_raw.csv`.

Typical use:

```bash
python /path/to/mrd-confidence-rater/scripts/mrdskill_score_mrd_evidence.py --input-csv data/evidence/mrd_evidence_raw.csv --out-csv data/evidence/mrd_evidence_scored.csv --out-jsonl data/evidence/mrd_evidence_scored.jsonl
```

Treat script output as rule-scored evidence. Do not use it as final clinical truth without reviewing urgent/high priority rows, possible CHIP rows, weakly traceable rows, and rows that would affect automated MRD QC.

It should take structured evidence records from `mrd-evidence-miner` and assign:

- final evidence source type
- confidence score
- evidence strength
- graph inclusion status
- audit priority
- downgrade reasons
- limitations
- human review status

Do not perform broad paper retrieval.
Do not rebuild the paper corpus.
Do not perform primary evidence extraction from scratch unless a specific evidence record needs verification.
Do not build Neo4j schemas.
Do not generate final graph ontology.

Those tasks belong to other skills:

- paper-search: corpus construction
- mrd-evidence-miner: raw evidence extraction
- mrd-kg-builder: ontology, Neo4j schema, and graph construction

## Primary Goal

Convert raw extracted evidence into confidence-scored evidence records.

Input may include:

- mrd_evidence_raw.csv
- mrd_evidence_raw.jsonl
- extracted evidence table
- evidence records from mrd-evidence-miner
- paper metadata
- DOI / PMID / PMCID
- source section, table, figure, or supplementary reference

Output should be a scored evidence table suitable for downstream knowledge graph construction.

## Core Principle

Confidence must be assigned per evidence record or per graph edge.

Do not assign confidence only at the paper level.

A high-quality paper may contain weak claims.
A lower-impact paper may contain a strong table-supported result.
A review may mention many useful ideas but should not be treated as primary evidence.

Every confidence decision must be traceable.

## Input Record Requirements

Each input evidence record may include:

- evidence_id
- paper_id
- title
- year
- journal
- doi
- pmid
- pmcid
- cancer_type
- cancer_subtype
- disease_stage
- treatment_setting
- gene
- variant
- variant_type
- assay_type
- sample_type
- mrd_relevance
- relapse_association
- recurrence_association
- prognosis_association
- chip_association
- false_positive_risk
- source_section
- source_table
- source_figure
- supplementary_file
- supplementary_table
- quoted_or_paraphrased_evidence
- cohort_type
- sample_size
- hazard_ratio
- confidence_interval
- p_value
- validation_status
- evidence_source_type_preliminary
- limitation
- extraction_notes

If essential fields are missing, do not invent them.

Use:

- not_reported
- unavailable
- uncertain

## Output Fields

For each evidence record, add:

- evidence_source_type_final
- confidence_score
- confidence_reason
- evidence_strength_score
- graph_inclusion_status
- audit_priority
- downgrade_reason
- upgrade_reason
- contradiction_status
- reproducibility_status
- clinical_relevance_status
- biological_plausibility_status
- source_traceability_status
- final_review_note

Recommended output columns:

- evidence_id
- paper_id
- title
- year
- journal
- doi
- pmid
- pmcid
- cancer_type
- cancer_subtype
- gene
- variant
- variant_type
- assay_type
- sample_type
- mrd_relevance
- relapse_association
- recurrence_association
- prognosis_association
- chip_association
- false_positive_risk
- source_section
- source_table
- source_figure
- supplementary_file
- supplementary_table
- cohort_type
- sample_size
- hazard_ratio
- confidence_interval
- p_value
- validation_status
- evidence_source_type_preliminary
- evidence_source_type_final
- confidence_score
- confidence_reason
- evidence_strength_score
- graph_inclusion_status
- audit_priority
- downgrade_reason
- upgrade_reason
- contradiction_status
- reproducibility_status
- clinical_relevance_status
- biological_plausibility_status
- source_traceability_status
- final_review_note

## Evidence Source Type

Assign one `evidence_source_type_final` per evidence record.

Allowed values:

### clinical_validation

Use when evidence comes from:

- prospective clinical study
- independent validation cohort
- multicenter cohort
- clinical trial
- guideline
- consensus recommendation
- clinically validated assay

Usually strong evidence if directly relevant to MRD / recurrence / relapse / prognosis.

### primary_result

Use when evidence comes from:

- Results section
- main table
- main figure
- statistical analysis
- original cohort result

This can be High or Medium depending on sample size, validation, statistics, and direct MRD relevance.

### supplementary_result

Use when evidence comes from:

- supplementary table
- supplementary figure
- supplementary dataset
- appendix with original data

This can be High or Medium if the supplementary evidence is original, traceable, and statistically supported.

### meta_analysis

Use when evidence comes from:

- systematic review
- meta-analysis
- pooled analysis

This can be High if methods are rigorous and the relationship is directly supported.

### methods_assay_detail

Use when the record describes:

- assay design
- sequencing depth
- limit of detection
- VAF cutoff
- matched WBC filtering
- germline filtering
- CHIP filtering
- technical validation

Usually Medium unless clinically validated.

### abstract_only_claim

Use when only the abstract supports the claim.

Never assign High confidence to abstract-only claims.

### discussion_or_background_claim

Use when the claim appears only in:

- Introduction
- Discussion
- background text
- author speculation

Usually Low.

### review_statement

Use when the claim appears in a review but is not supported by directly reviewed original data in the current extraction.

Usually Low or Medium depending on citation support, but not High unless primary evidence is verified.

### inferred_candidate

Use when the relation is inferred indirectly from multiple statements.

Usually Low unless later verified by direct evidence.

### uncertain

Use when source support is unclear.

Usually Low or needs_review.

## Confidence Score

Assign one `confidence_score` per evidence record.

Allowed values:

- High
- Medium
- Low
- Exclude
- Uncertain

## High Confidence Criteria

Assign High only when most of the following are true:

- direct MRD / ctDNA / recurrence / relapse / prognosis relevance
- original clinical or experimental data
- Results, Table, Figure, or Supplementary Table support
- clear cancer context
- clear gene / variant / assay context
- cohort size is adequate for the claim
- statistical support is reported
- validation cohort or prospective design is present, or strong reproducibility exists
- DOI / PMID / PMCID traceability exists
- source table / figure / section is traceable

Examples:

- A prospective study shows postoperative ctDNA positivity predicts recurrence in colorectal cancer.
- A validation cohort shows MRD-positive patients have worse disease-free survival.
- A supplementary table lists recurrent ctDNA variants with cancer type, assay, recurrence outcome, and statistics.
- A guideline recommends MRD testing in a specific clinical setting.

## Medium Confidence Criteria

Assign Medium when evidence is useful but limited.

Examples:

- original data exists but sample size is small
- single-center retrospective cohort
- no independent validation
- statistical support is incomplete
- MRD relevance is indirect
- gene association is clear but clinical outcome is weak
- supplementary data is available but not fully described
- assay performance is reported but not linked to clinical outcome

## Low Confidence Criteria

Assign Low when evidence is weak or indirect.

Examples:

- abstract-only claim
- discussion-only claim
- review statement without primary verification
- inferred relation
- no clear cancer context
- no clear assay context
- no statistical support
- unclear cohort
- missing DOI / PMID
- low traceability
- claim is plausible but not demonstrated

## Exclude Criteria

Assign Exclude when the record should not enter the graph.

Use Exclude if:

- the record is not about MRD / ctDNA / relapse / recurrence / oncology evidence
- gene or cancer type is missing and cannot be resolved
- evidence is unsupported or hallucinated
- source is not traceable
- the claim contradicts the cited paper
- the record is a duplicate artifact with no added evidence
- the evidence is purely speculative and not useful even as Low confidence

## Evidence Strength Score

Also assign a numeric `evidence_strength_score` from 0 to 5.

Use:

### 5

Very strong evidence.

Usually:

- prospective or validation cohort
- direct MRD / recurrence endpoint
- strong statistics
- traceable table / figure / supplementary evidence
- clinically meaningful effect

### 4

Strong evidence.

Usually:

- original data
- clear MRD relevance
- good traceability
- some statistical support
- may lack independent validation

### 3

Moderate evidence.

Usually:

- original data but limited cohort
- indirect MRD relevance
- partial statistics
- useful but not definitive

### 2

Weak evidence.

Usually:

- abstract-only
- review statement
- discussion claim
- inferred relation
- unclear statistics

### 1

Very weak evidence.

Usually:

- speculative
- poor traceability
- missing key fields
- not directly MRD-related

### 0

Exclude.

Unsupported, irrelevant, contradictory, or unusable.

## Graph Inclusion Status

Assign `graph_inclusion_status`.

Allowed values:

### auto_include_high_confidence

Use when evidence can enter the high-confidence graph automatically.

Usually requires:

- confidence_score = High
- evidence_strength_score >= 4
- traceable source
- direct MRD / recurrence / relapse relevance

### include_as_medium_confidence

Use when evidence is useful but should be tagged Medium in the graph.

### include_as_low_confidence_annotation

Use when evidence may be useful as background or weak annotation but not for automated QC decisions.

### candidate_only

Use when evidence is interesting but not ready for graph inclusion.

### needs_human_review

Use when evidence may be important but has ambiguity, missing fields, or possible conflict.

### exclude

Use when record should not enter the graph.

## Audit Priority

Assign `audit_priority`.

Allowed values:

- urgent_review
- high_review
- medium_review
- low_review
- no_review_needed

Use urgent_review when:

- claim may strongly affect MRD QC decision
- CHIP / false positive interpretation is ambiguous
- high clinical impact but evidence is weak
- contradiction exists
- gene is often both tumor-derived and CHIP-related

Use high_review when:

- record may become high-confidence after verification
- supplementary table needs manual checking
- key field is uncertain

Use medium_review when:

- useful but not urgent

Use low_review when:

- weak background evidence

Use no_review_needed when:

- evidence is strong and traceable
- or excluded clearly

## MRD-Specific Rating Rules

### Direct MRD evidence

Upgrade confidence when:

- MRD positivity is linked to recurrence
- ctDNA persistence predicts relapse
- ctDNA clearance predicts response
- longitudinal monitoring shows molecular relapse before radiology
- MRD status is linked to survival endpoint
- validation cohort supports the finding

### Supplementary variant-table evidence

Treat rows from supplementary tables as graph-eligible evidence when they contain structured gene, variant, VAF, sample, ctDNA, tissue, timepoint, assay, response, recurrence, relapse, or survival fields with sheet/row/file traceability.

Do not automatically upgrade every supplementary row to High. Use Medium or candidate status for structured variant evidence that supports entity/assay/variant edges but lacks direct recurrence, relapse, survival, validation, or statistical support. Upgrade only when the row or linked table clearly provides clinical endpoint support, validation, or strong statistics.

### Indirect MRD evidence

Keep Medium or Low when:

- gene is a known cancer driver but not shown in MRD context
- assay detects variants but no recurrence/prognosis data
- paper discusses MRD conceptually without patient outcome data

### Tumor-informed assays

Upgrade when:

- tumor tissue was used to design assay
- patient-specific variants were tracked
- matched WBC or germline filtering was used
- recurrence outcomes were reported

### Tumor-naive assays

Be more cautious.

Upgrade only if:

- assay was clinically validated
- false positives were controlled
- CHIP filtering was addressed
- recurrence/prognosis data are strong

## CHIP / False Positive Rating Rules

CHIP-related evidence requires special caution.

Common CHIP-associated genes include:

- DNMT3A
- TET2
- ASXL1
- JAK2
- TP53
- PPM1D
- SF3B1
- SRSF2
- U2AF1
- IDH1
- IDH2
- ATM
- CHEK2
- GNAS
- KRAS

Do not automatically mark every mutation in these genes as CHIP.

Consider:

- cancer type
- variant context
- VAF pattern
- matched WBC evidence
- tumor tissue confirmation
- longitudinal dynamics
- treatment exposure
- age
- assay method
- plasma-only status

### Upgrade CHIP false-positive evidence when:

- matched WBC sequencing confirms CHIP origin
- plasma-only variant is absent from tumor tissue
- variant appears in known CHIP gene with supporting evidence
- paper directly reports CHIP-related false positives
- filtering method is described and validated

### Downgrade CHIP interpretation when:

- no matched WBC sequencing
- no tumor tissue comparison
- claim is only discussion-based
- gene is both tumor driver and CHIP-associated
- cancer context makes tumor origin plausible
- VAF dynamics support tumor-derived ctDNA

## Biological Plausibility Rating

Assign `biological_plausibility_status`.

Allowed values:

- biologically_plausible
- biologically_uncertain
- biologically_questionable
- likely_false_positive
- likely_CHIP
- likely_technical_noise
- not_assessable

Consider:

- cancer type
- gene
- variant type
- assay type
- tumor tissue confirmation
- matched WBC confirmation
- recurrence association
- VAF trend
- longitudinal dynamics
- known driver biology
- known CHIP biology

## Source Traceability Status

Assign `source_traceability_status`.

Allowed values:

- fully_traceable
- partially_traceable
- weakly_traceable
- not_traceable

### fully_traceable

Requires:

- DOI or PMID
- paper title
- source section
- table / figure / supplementary reference if applicable
- clear evidence text or paraphrase

### partially_traceable

Paper is identifiable but exact table / figure / section is incomplete.

### weakly_traceable

Only title or abstract is available.

### not_traceable

No reliable source.

## Downgrade Reasons

Use one or more `downgrade_reason` values:

- abstract_only
- discussion_only
- review_only
- no_full_text
- no_supplementary_access
- missing_statistics
- small_sample_size
- no_validation_cohort
- retrospective_only
- unclear_cancer_context
- unclear_gene_context
- unclear_assay_context
- no_matched_WBC
- no_tumor_tissue_confirmation
- possible_CHIP
- possible_germline
- possible_technical_noise
- weak_traceability
- conflicting_evidence
- not_directly_MRD_related

## Upgrade Reasons

Use one or more `upgrade_reason` values:

- prospective_study
- validation_cohort
- large_cohort
- guideline_or_consensus
- direct_recurrence_endpoint
- direct_relapse_endpoint
- survival_endpoint
- longitudinal_monitoring
- matched_WBC_filtering
- tumor_tissue_confirmation
- strong_statistics
- supplementary_table_support
- main_table_or_figure_support
- independent_replication
- clinically_validated_assay

## Contradiction Handling

Assign `contradiction_status`.

Allowed values:

- no_contradiction_detected
- possible_contradiction
- conflicting_literature
- internally_inconsistent
- not_assessable

If contradiction exists:

- do not assign auto_include_high_confidence
- mark needs_human_review
- explain the conflict briefly

## Output Format

Default output should be a structured table.

Recommended CSV filename:

- data/evidence/mrd_evidence_scored.csv

Recommended JSONL filename:

- data/evidence/mrd_evidence_scored.jsonl

Recommended audit notes file:

- docs/evidence/confidence_rating_notes.md

## JSON Output Schema

If the user asks for JSON, use:

{
  "evidence_id": "",
  "paper_id": "",
  "source": {
    "title": "",
    "year": "",
    "journal": "",
    "doi": "",
    "pmid": "",
    "pmcid": "",
    "source_section": "",
    "source_table": "",
    "source_figure": "",
    "supplementary_file": "",
    "supplementary_table": ""
  },
  "edge_candidate": {
    "cancer_type": "",
    "gene": "",
    "variant": "",
    "variant_type": "",
    "assay_type": "",
    "mrd_relevance": "",
    "clinical_meaning": ""
  },
  "rating": {
    "evidence_source_type_final": "",
    "confidence_score": "",
    "evidence_strength_score": "",
    "confidence_reason": "",
    "graph_inclusion_status": "",
    "audit_priority": "",
    "downgrade_reason": [],
    "upgrade_reason": [],
    "contradiction_status": "",
    "biological_plausibility_status": "",
    "source_traceability_status": ""
  },
  "final_review_note": ""
}

## Anti-Hallucination Rules

Do not fabricate:

- confidence justification
- statistics
- cohort size
- validation status
- PMID
- DOI
- table number
- figure number
- supplementary table number
- CHIP status
- biological interpretation

If evidence is missing, mark:

- uncertain
- not_reported
- unavailable
- not_assessable

Do not upgrade evidence because the claim sounds biologically plausible.

Do not upgrade evidence solely because the journal is prestigious.

Do not downgrade evidence solely because the journal is not a top journal.

Rate the actual evidence record, not the prestige of the paper.

## Boundary

This skill ends after evidence confidence rating.

Next step:

- Use mrd-kg-builder to build ontology, Neo4j schema, graph-ready nodes and relationships, and MRD biological QC rules.
