# Data Availability

This release contains source code and configuration templates, plus an explicitly synthetic, CC0-1.0 software-test fixture under `examples/demo/` and its small precomputed outputs under `examples/expected/`. These fixtures contain no patient-level data, external paper text, real evidence claims, or benchmark observations. They exercise the KG/RAG code path only and are not for biological interpretation or clinical use.

The release does not include full text, supplementary files, patient-level data, VCFs, clinical records, benchmark cohorts, the complete graph export, or the complete retrieval index.

`full_kg_overview.html` is a static aggregate report from a complete local build. It contains row counts and distributions, not graph entities, evidence text, source documents, or retrievable RAG chunks. Local machine paths are omitted from the public copy.

| Input class | Expected local location | Release status | Access requirement |
| --- | --- | --- | --- |
| OA literature and supplements | `artifacts/runs/` or `artifacts/imports/literature_search/` | Not included | Retrieve only from lawful access points and preserve source URL and license metadata. |
| User VCF and case context | User-selected path or `artifacts/cases/` | Not included | Keep under the user's approved clinical or research governance. |
| Scored evidence and KG/RAG artifacts | `artifacts/data/` | Not included | Build locally from permitted inputs; do not commit generated outputs. |
| Synthetic software demo | `examples/demo/`, `examples/expected/` | Included (CC0-1.0) | Synthetic test data only; not research evidence, benchmark data, or for clinical use. |
| Razavi and other external benchmark inputs | User-selected path under `artifacts/benchmarks/` | Not included | Obtain directly under the original study, repository, and data-use terms. |

The reproducibility scripts intentionally fail when required external inputs are absent. They must not be used to publish or redistribute restricted data.
