# Data Availability

This release contains source code and configuration templates only. It does not include full text, supplementary files, patient-level data, VCFs, clinical records, benchmark cohorts, graph exports, or retrieval indexes.

| Input class | Expected local location | Release status | Access requirement |
| --- | --- | --- | --- |
| OA literature and supplements | `artifacts/runs/` or `artifacts/imports/literature_search/` | Not included | Retrieve only from lawful access points and preserve source URL and license metadata. |
| User VCF and case context | User-selected path or `artifacts/cases/` | Not included | Keep under the user's approved clinical or research governance. |
| Scored evidence and KG/RAG artifacts | `artifacts/data/` | Not included | Build locally from permitted inputs; do not commit generated outputs. |
| Razavi and other external benchmark inputs | User-selected path under `artifacts/benchmarks/` | Not included | Obtain directly under the original study, repository, and data-use terms. |

The reproducibility scripts intentionally fail when required external inputs are absent. They must not be used to publish or redistribute restricted data.
