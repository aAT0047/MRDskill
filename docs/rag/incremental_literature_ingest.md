# MRDskill Incremental Literature Ingest

This project keeps a stable interface for adding new papers and supplementary
materials into the MRD knowledge graph and RAG layer.

## Prepared Input

Place each new literature batch in a folder that follows the existing
`kg_ready_no_pdf` style as much as possible:

- `manifests/keyword_research_candidate_table.csv`
- `manifests/downloaded_pdfs_converted_manifest.csv`
- converted full text under `papers/`, `converted/`, or equivalent folders
- supplementary materials under `supplementary_materials/`

The ingest script can still scan parseable files when some manifest metadata is
missing, but paper/journal/DOI/PMID traceability is best when manifests are
present.

## One-Step Incremental Update

```powershell
python scripts\mrdskill_update_mrd_kg.py --mode incremental --input ".\artifacts\inputs\new_papers\batch_001"
```

This runs:

1. `mrdskill_ingest_new_literature.py`
2. evidence extraction from full text and supplementary materials
3. confidence scoring
4. de-duplicated merge into the full evidence/import layers
5. augmented KG/RAG/latest SQLite index rebuild while preserving existing cohort-rule CSVs
6. latest SQLite rule index and unified RAG JSONL rebuild

Use `--rebuild-cohort-rules` only when you intentionally want to re-scan all
registered supplementary sources and overwrite the current baseline,
co-occurrence, and mutual-exclusion rule CSVs.

## Ingest Only

```powershell
python scripts\mrdskill_ingest_new_literature.py --input ".\artifacts\inputs\new_papers\batch_001"
```

Outputs are staged under:

```text
artifacts\imports\literature_ingest\<run_id>\
```

The script then merges de-duplicated rows into:

- `data/evidence/mrd_evidence_raw_full.csv`
- `data/evidence/mrd_evidence_scored_full.csv`
- `data/evidence/mrd_evidence_scored_full.jsonl`
- `data/imported/full_corpus_sources.csv`
- `data/imported/full_corpus_excel_chunks.csv`

## Rebuild Only

```powershell
python scripts\mrdskill_update_mrd_kg.py --mode rebuild
```

This rebuilds KG/RAG/index from the current full evidence/import layers without
adding a new literature batch, preserving existing cohort-rule CSVs by default.

## Important Behavior

- Ingest is non-destructive and archive-first.
- Evidence IDs for incremental rows use stable content hashes:
  `EV-INCR-<hash>`.
- Re-ingesting the same batch should skip duplicate evidence/source/chunk rows.
- KG/RAG rebuild archives graph/RAG/cohort outputs, but does not move case QC
  folders or existing sample-level QC reports.
- Every ingest run writes a JSON report with command outputs and merge counts.

## Main Latest Outputs

- `data/graph/nodes_full.csv`
- `data/graph/relationships_full.csv`
- `data/graph/nodes_augmented.csv`
- `data/graph/relationships_augmented.csv`
- `data/rag/mrd_latest_unified_rag_chunks.jsonl`
- `data/rag/mrd_latest_kg_rule_index.sqlite`
- `docs/graph/mrd_kg_visualization.html`
