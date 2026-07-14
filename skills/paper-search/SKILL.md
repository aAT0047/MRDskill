---
name: paper-search
description: Build auditable MRD, ctDNA, CHIP, liquid-biopsy, recurrence-monitoring, and oncology candidate literature sets. Use when Codex needs a database search strategy, journal/date filters, candidate paper table, or corpus-entry decision before full-text retrieval.
---

# Paper Search

Use this skill to create the candidate literature layer for MRDskill.

## Bundled Command

Run `scripts/search.py` with the same arguments as the shared harvest runner. The wrapper always appends `--search-only`, so it writes search outputs without downloading papers or supplements.

```powershell
python .\skills\paper-search\scripts\search.py --output-root ".\artifacts\runs" --config ".\config\literature_harvest_config.template.json" --run-name "mrd_2016_2026"
```

## Workflow

1. Define MRD/ctDNA/CHIP terms, cancer context, journal scope, date range, and exclusions.
2. Search PubMed, PMC, Europe PMC, Crossref, OpenAlex, and optionally a permitted Scholar export or SerpAPI result.
3. Deduplicate records by DOI, PMID, PMCID, normalized title, year, and first author.
4. Write a candidate table, high-priority table, medium-priority table, query configuration, and search summary.

## Outputs

- Candidate literature table.
- High-priority and medium-priority candidate tables.
- Search configuration and query log.

## QC Boundary

Do not download restricted full text, extract biological evidence, assign evidence confidence, or make source-origin decisions.
