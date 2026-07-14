---
name: mrd-evidence-miner
description: Use when the user wants to extract structured MRD, ctDNA, CHIP, cancer-gene, variant, assay, relapse, recurrence, prognosis, or false-positive evidence from biomedical papers, OA full text, tables, figures, or supplementary materials. This skill assumes papers have already been retrieved by paper-search. It does not perform large-scale paper retrieval, final confidence scoring, or knowledge graph construction.
allowed-tools:
  - readPaper
  - searchArticles
---

# MRD Evidence Miner

Use this skill to extract structured evidence from biomedical literature for MRD biological QC.

This skill is the second step in the MRD biological QC pipeline.

## Python Helpers

Prefer bundled Python scripts for repeatable corpus processing before doing manual or LLM-assisted review.

- `scripts/mrdskill_inventory_kg_ready.py`: inventory a `kg_ready_no_pdf` corpus, summarize file types, manifests, years, journals, and topic tags.
- `scripts/mrdskill_extract_mrd_candidates.py`: extract conservative MRD evidence candidates from converted `.txt` paper files and candidate-table metadata.
- `scripts/mrdskill_extract_full_corpus_evidence.py`: extract conservative MRD evidence candidates from all parseable `kg_ready_no_pdf` materials, including paper TXT/HTML/XML, supplementary converted TXT/HTML/XML, native CSV/TSV/XLSX/DOCX/TXT/ZIP, and best-effort DOC/XLS strings.
- `scripts/mrdskill_split_scan_excel.py`: for very large supplementary `.xlsx` files, split sheets into row-chunk CSV files, scan chunks in parallel, show split/scan progress bars, and preserve sheet/row/column traceability.
- `scripts/convert_pdf_opendataloader.py`: optional PDF-to-HTML/JSON/text conversion through the external `opendataloader-pdf` CLI.
- `scripts/convert_pdf_pdfminer.py`: portable text-first fallback when the optional converter is unavailable or times out.

Typical use:

```bash
python /path/to/mrd-evidence-miner/scripts/mrdskill_inventory_kg_ready.py --corpus-dir /path/to/kg_ready_no_pdf --out-json data/imported/kg_ready_inventory.json
python /path/to/mrd-evidence-miner/scripts/mrdskill_extract_mrd_candidates.py --corpus-dir /path/to/kg_ready_no_pdf --out-csv data/evidence/mrd_evidence_raw.csv --out-jsonl data/evidence/mrd_evidence_raw.jsonl
python /path/to/mrd-evidence-miner/scripts/mrdskill_extract_full_corpus_evidence.py --corpus-dir /path/to/kg_ready_no_pdf --out-csv data/evidence/mrd_evidence_raw_full.csv --out-jsonl data/evidence/mrd_evidence_raw_full.jsonl --sources-csv data/imported/full_corpus_sources.csv --report-json data/imported/full_corpus_ingest_report.json
python /path/to/mrd-evidence-miner/scripts/mrdskill_split_scan_excel.py --xlsx /path/to/supplement.xlsx --candidate-table /path/to/kg_ready_no_pdf/manifests/keyword_research_candidate_table.csv --chunk-dir data/imported/excel_chunks --out-csv data/evidence/mrd_evidence_raw_excel.csv --chunk-log-csv data/imported/excel_chunk_log.csv --workers 4
```

Treat script output as candidate evidence. Preserve source traceability and review rows with weak source sections, no detected gene, possible CHIP, or missing statistics.

It should read and mine evidence from:

- OA full text
- Results sections
- Methods sections
- Tables
- Figures
- Supplementary tables
- Supplementary figures
- Supplementary datasets
- appendices

It should transform papers into structured evidence records.

Do not build the full paper corpus.
Do not perform broad paper discovery unless a paper is missing and needs metadata completion.
Do not assign final High / Medium / Low edge confidence.
Do not build Neo4j schemas.
Do not generate final MRD QC rules.

Those tasks belong to other skills:

- paper-search: corpus construction
- mrd-confidence-rater: final evidence confidence scoring
- mrd-kg-builder: ontology, Neo4j schema, and graph construction

## Primary Goal

Extract structured Cancer-Gene-Variant-Assay-MRD evidence from papers.

Input may include:

- a list of papers
- DOI
- PMID
- PMCID
- paper title
- OA full text
- supplementary tables
- supplementary files
- CSV from paper-search
- markdown paper list
- user-provided paper notes

Output should be a structured evidence table suitable for later confidence scoring and knowledge graph construction.

## Scope

Focus on MRD / ctDNA / oncology evidence, including:

- minimal residual disease
- measurable residual disease
- molecular residual disease
- ctDNA
- circulating tumor DNA
- liquid biopsy
- recurrence monitoring
- relapse prediction
- molecular relapse
- treatment monitoring
- tumor-informed MRD
- tumor-naive MRD
- CHIP
- clonal hematopoiesis
- false positive variants
- germline filtering
- matched white blood cell sequencing
- ultra-low VAF detection
- methylation MRD
- fragmentomics
- longitudinal ctDNA tracking

## Target Cancer Types

If the user does not specify cancer types, prioritize:

- colorectal cancer
- non-small cell lung cancer
- breast cancer
- AML
- MDS
- lymphoma
- pancreatic cancer

Also allow:

- gastric cancer
- ovarian cancer
- bladder cancer
- melanoma
- hepatocellular carcinoma
- esophageal cancer
- head and neck cancer

## Core Extraction Principle

Extract evidence only when the paper provides traceable support.

Prefer evidence from:

1. Results
2. Tables
3. Figures
4. Supplementary tables
5. Supplementary figures
6. Validation cohorts
7. Statistical analysis
8. Methods, only for assay or technical details

Treat structured supplementary variant tables as evidence candidates when they contain traceable gene, variant, VAF, sample, ctDNA, tissue, timepoint, assay, response, recurrence, relapse, or survival columns. These rows may support graph edges such as Gene-HAS_VARIANT, Variant-OBSERVED_IN, Assay-DETECTS_VARIANT, or MRDSignal-DETECTED_BY even when they do not by themselves prove recurrence or prognosis. Preserve sheet name, row number, column hits, table title, and file name.

Be cautious with:

- Abstract
- Introduction
- Discussion
- Review statements
- author speculation
- background statements

If a claim appears only in Discussion or Introduction, label it as discussion_or_background_claim.

If a claim appears only in Abstract and no full-text support is available, label it as abstract_only_claim.

## What To Extract

For each evidence item, extract the following fields.

### Paper Metadata

- paper_id
- title
- authors
- year
- journal
- doi
- pmid
- pmcid
- source_url
- source_database
- full_text_status
- supplementary_status

### Cancer Information

- cancer_type
- cancer_subtype
- disease_stage
- treatment_setting
- sample_context

Examples:

- colorectal cancer
- stage II colon cancer
- resected NSCLC
- metastatic breast cancer
- AML post-transplant
- MDS monitoring
- diffuse large B-cell lymphoma
- pancreatic ductal adenocarcinoma

### Gene / Biomarker Information

- gene
- gene_symbol_normalized
- transcript
- variant
- variant_type
- mutation_hotspot
- genomic_coordinate
- protein_change
- pathway
- biomarker_type

variant_type values may include:

- SNV
- indel
- fusion
- CNV
- SV
- methylation
- fragmentomics
- copy-number pattern
- mutational signature
- expression marker
- epigenetic marker

biomarker_type values may include:

- tumor-derived_ctDNA
- CHIP_related
- germline_risk
- resistance_marker
- recurrence_marker
- prognostic_marker
- treatment_response_marker
- technical_noise_candidate

### MRD / Clinical Meaning

Extract:

- mrd_relevance
- relapse_association
- recurrence_association
- prognosis_association
- treatment_monitoring_value
- resistance_association
- molecular_relapse_signal
- recurrence_lead_time
- survival_endpoint
- clinical_actionability
- biological_plausibility

mrd_relevance values may include:

- recurrence_marker
- relapse_driver
- prognostic_marker
- treatment_monitoring_marker
- resistance_marker
- CHIP_related
- age_related_background
- low_specificity_marker
- likely_false_positive
- tumor_specific_marker
- methylation_marker
- fragmentomics_marker
- uncertain

### Assay / Technical Information

Extract:

- assay_type
- assay_platform
- tumor_informed_or_tumor_naive
- sample_type
- sequencing_depth
- limit_of_detection
- VAF_cutoff
- assay_sensitivity
- assay_specificity
- positive_predictive_value
- negative_predictive_value
- matched_WBC_used
- germline_filtering_used
- CHIP_filtering_used

assay_type values may include:

- ddPCR
- targeted_NGS
- ultra_deep_sequencing
- WES
- WGS
- methylation_assay
- personalized_panel
- tumor_informed_assay
- tumor_naive_assay
- fragmentomics_assay
- hybrid_capture_panel
- amplicon_panel

sample_type values may include:

- plasma
- serum
- bone_marrow
- peripheral_blood
- tissue
- matched_WBC
- buffy_coat

### Evidence Location

Every extracted record must include traceability fields:

- source_section
- source_subsection
- source_table
- source_figure
- supplementary_file
- supplementary_table
- quoted_or_paraphrased_evidence
- page_or_line_reference_if_available

source_section values may include:

- Abstract
- Introduction
- Methods
- Results
- Discussion
- Table
- Figure
- Supplementary_Table
- Supplementary_Figure
- Supplementary_Dataset
- Appendix

## Evidence Type Labeling

This skill does not assign final confidence scores, but it must label the evidence type.

Use `evidence_source_type_preliminary`:

- primary_result
- clinical_validation
- meta_analysis
- supplementary_result
- methods_assay_detail
- abstract_only_claim
- discussion_or_background_claim
- review_statement
- inferred_candidate
- uncertain

Definitions:

### primary_result

Use when the evidence comes from Results, main tables, main figures, or statistical analysis.

### clinical_validation

Use when the evidence is from a validation cohort, prospective trial, independent cohort, guideline, or clinically validated setting.

### meta_analysis

Use when the evidence comes from a meta-analysis or systematic review.

### supplementary_result

Use when the evidence comes from supplementary tables, supplementary figures, or supplementary datasets.

### methods_assay_detail

Use for assay design, sample processing, sequencing depth, LOD, VAF cutoff, WBC filtering, or technical validation details from Methods.

### abstract_only_claim

Use when only the abstract supports the statement and no full-text/table/figure/supplementary evidence is available.

### discussion_or_background_claim

Use when the claim appears only in Introduction, Discussion, or background text.

### review_statement

Use when the claim is from a review without original data.

### inferred_candidate

Use when the relation is indirectly inferred from multiple statements, not directly stated.

### uncertain

Use when source support is unclear.

## CHIP and False Positive Mining

Pay special attention to CHIP and false-positive evidence.

Extract whether a gene or variant is associated with:

- CHIP
- clonal hematopoiesis
- age-related clonal hematopoiesis
- matched WBC signal
- germline signal
- sequencing artifact
- low VAF noise
- plasma-only false positive
- non-tumor origin
- technical error

For CHIP-related evidence, extract:

- chip_gene
- chip_variant
- chip_context
- matched_WBC_confirmation
- tumor_tissue_confirmation
- plasma_only_status
- false_positive_risk
- filtering_recommendation

Common CHIP-related genes to watch for include, but are not limited to:

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
Always preserve cancer context, tissue confirmation, VAF pattern, matched WBC status, and paper evidence.

## Longitudinal MRD Pattern Mining

Extract longitudinal patterns when available:

- preoperative_ctDNA_status
- postoperative_ctDNA_status
- post_treatment_ctDNA_status
- serial_monitoring_status
- ctDNA_clearance
- ctDNA_persistence
- ctDNA_reappearance
- time_to_recurrence
- recurrence_lead_time
- molecular_relapse_before_radiology
- association_with_survival

## Extraction Output Format

Default output should be a structured table or CSV-like schema.

Recommended columns:

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
- gene_symbol_normalized
- variant
- variant_type
- mutation_hotspot
- genomic_coordinate
- pathway
- biomarker_type
- assay_type
- assay_platform
- tumor_informed_or_tumor_naive
- sample_type
- VAF_cutoff
- limit_of_detection
- assay_sensitivity
- assay_specificity
- matched_WBC_used
- germline_filtering_used
- CHIP_filtering_used
- mrd_relevance
- relapse_association
- recurrence_association
- prognosis_association
- treatment_monitoring_value
- resistance_association
- recurrence_lead_time
- survival_endpoint
- hazard_ratio
- confidence_interval
- p_value
- cohort_type
- sample_size
- validation_status
- evidence_source_type_preliminary
- source_section
- source_table
- source_figure
- supplementary_file
- supplementary_table
- quoted_or_paraphrased_evidence
- limitation
- extraction_notes

## JSON Output Option

If the user asks for JSON, use this schema:

{
  "evidence_id": "",
  "paper_id": "",
  "paper_metadata": {
    "title": "",
    "authors": "",
    "year": "",
    "journal": "",
    "doi": "",
    "pmid": "",
    "pmcid": "",
    "source_url": "",
    "full_text_status": "",
    "supplementary_status": ""
  },
  "cancer": {
    "cancer_type": "",
    "cancer_subtype": "",
    "disease_stage": "",
    "treatment_setting": ""
  },
  "biomarker": {
    "gene": "",
    "gene_symbol_normalized": "",
    "variant": "",
    "variant_type": "",
    "mutation_hotspot": "",
    "genomic_coordinate": "",
    "pathway": "",
    "biomarker_type": ""
  },
  "assay": {
    "assay_type": "",
    "assay_platform": "",
    "tumor_informed_or_tumor_naive": "",
    "sample_type": "",
    "VAF_cutoff": "",
    "limit_of_detection": "",
    "assay_sensitivity": "",
    "assay_specificity": "",
    "matched_WBC_used": "",
    "germline_filtering_used": "",
    "CHIP_filtering_used": ""
  },
  "clinical_meaning": {
    "mrd_relevance": "",
    "relapse_association": "",
    "recurrence_association": "",
    "prognosis_association": "",
    "treatment_monitoring_value": "",
    "resistance_association": "",
    "recurrence_lead_time": "",
    "survival_endpoint": ""
  },
  "statistics": {
    "hazard_ratio": "",
    "confidence_interval": "",
    "p_value": "",
    "sample_size": "",
    "cohort_type": "",
    "validation_status": ""
  },
  "traceability": {
    "evidence_source_type_preliminary": "",
    "source_section": "",
    "source_table": "",
    "source_figure": "",
    "supplementary_file": "",
    "supplementary_table": "",
    "quoted_or_paraphrased_evidence": ""
  },
  "limitations": "",
  "extraction_notes": ""
}

## Normalization Rules

Normalize where possible:

- gene symbols to HGNC-style symbols
- cancer names to common oncology terminology
- assay names to controlled vocabulary
- variant types to controlled vocabulary
- MRD relevance to controlled vocabulary

If normalization is uncertain, preserve the original term and add:

- normalized_term_uncertain

Do not silently replace ambiguous terms.

## Deduplication Rules

If the same paper contains multiple cancer-gene-MRD relationships, create multiple evidence records.

If the same relationship appears in multiple papers, keep separate records.

Do not merge across papers during this skill.

Deduplication and evidence consolidation belong to later skills.

## Anti-Hallucination Rules

Do not fabricate:

- gene
- variant
- cancer type
- cohort size
- hazard ratio
- p value
- DOI
- PMID
- table number
- supplementary table number
- assay sensitivity
- assay specificity
- recurrence lead time

If missing, write:

- not_reported
- unavailable
- uncertain

Do not infer a gene-level relationship from general text unless clearly marked as inferred_candidate.

Do not treat a review statement as primary evidence.

Do not treat abstract-only evidence as table-supported evidence.

## Output Recommendations

When working in a project directory, recommend saving outputs as:

- data/evidence/mrd_evidence_raw.csv
- data/evidence/mrd_evidence_raw.jsonl
- docs/evidence/extraction_notes.md

Do not overwrite existing files without user instruction.

## Boundary

This skill ends after structured evidence extraction.

Next step:

- Use mrd-confidence-rater to assign final High / Medium / Low confidence.
- Use mrd-kg-builder to build ontology, Neo4j schema, and biological QC rules.
