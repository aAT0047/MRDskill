---
name: literature-harvest
description: Download and organize legally accessible MRD, ctDNA, CHIP and oncology papers plus supplementary materials; generate paper manifests, supplement manifests and download logs. Use when Codex needs to run or maintain the literature-harvest workflow after candidate papers have been selected.
---

# Literature Harvest

Use this skill after candidate papers have been selected by `paper-search`.

## Inputs

- Candidate paper table, DOI/PMID/PMCID list or source manifest.
- Access policy, output root and download configuration.
- Optional retry list for failed downloads.

## Workflow

1. Resolve legally accessible full text, preferring publisher HTML/XML or PMC OA files when available and PDF as fallback.
2. Download supplementary materials exposed by API metadata, PMC packages, accessible article HTML/XML or user-provided links.
3. Normalize file names with a shared record identifier so papers and supplements remain grouped.
4. Write paper manifests, supplement manifests, package logs and failed-download logs.
5. Deduplicate paper and supplement files separately.

## Bundled Commands

- `scripts/harvest.py`: search and retrieve accessible papers and supplements.
- `scripts/continue_harvest.py`: retry failed records, discover supplements, and deduplicate.
- `scripts/clean_outputs.py`: quarantine invalid files and rebuild `papers_manifest.csv` and `supplements_manifest.csv`.
- `scripts/fetch_latest_oa.py`: fetch a small Europe PMC OA batch for an incremental ingest test.

## Outputs

- `papers_manifest.csv`.
- `supplements_manifest.csv`.
- Download logs, retry logs and deduplication manifests.

## QC Boundary

Do not bypass access restrictions. Do not scrape Google Scholar HTML directly. Do not interpret biological meaning or extract evidence records.
