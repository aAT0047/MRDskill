# MRDskill Source Package

This repository contains the source code for the MRDskill Agent-Skill workflow.

The operational contract is the ten skills in `skills/registry.csv`:

- paper search
- literature harvest
- evidence extraction
- confidence rating
- evidence graph construction
- supplementary-rule materialization
- multi-round GraphRAG
- VCF quality control
- review-aware adjudication
- incremental update and audit

Keep shared implementations in `scripts/` and expose them through the matching `skills/<skill>/scripts/` wrapper or bundled script. Keep the skill registry synchronized whenever an operational skill changes.

Do not commit downloaded papers, supplementary files, patient data, VCFs, evidence tables, graph exports, RAG indexes, browser profiles, API keys, local paths, or generated artifacts. Store those under `artifacts/`, which is intentionally ignored.

MRDskill outputs are research and review-support artifacts. They must preserve provenance and uncertainty, and must not be presented as autonomous clinical sign-out or treatment decisions.
