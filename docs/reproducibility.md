# Reproducibility

## Offline smoke demo (no external data)

Requirements: Python 3.10+ and Git; no third-party Python packages, network access, API key, downloaded paper, or patient/benchmark input is needed.

From a clean checkout:

```powershell
python .\scripts\run_demo.py
python .\scripts\run_demo.py --output-dir .\artifacts\demo --verify-expected
```

The first command builds a tiny KG from the synthetic scored-evidence fixture, turns its evidence summary plus a synthetic retrieval-vocabulary fixture into a JSONL RAG index, runs the real seven-round source-tracing retrieval implementation, and writes an interactive KG HTML. The second command also compares deterministic outputs to `examples/expected/` using SHA-256. `artifacts/` is ignored by Git, so generated files stay local. On Windows, `scripts/run_demo.ps1` forwards the same arguments.

Expected outputs include `nodes.csv`, `relationships.csv`, `review_tasks.csv`, `evidence_summaries.jsonl`, `rag_chunks.jsonl`, `rag_retrieval.json`, `kg_visualization_data.json`, `kg_demo.html`, a compact `summary.json`, and a SHA-256 `manifest.json`. The checked-in expected output reports seven graph nodes, eleven relations, two retrieval chunks and seven rounds with hits. These are software-path counts from synthetic text—not model accuracy, clinical performance, biology, or benchmark results.

Every demo input is labeled synthetic and has no real PMID/DOI, patient, source paper, or cohort. `examples/demo/manifest.json` records provenance and CC0-1.0 terms. Do not infer that the generic vocabulary snippets are factual evidence or medical rules. Demo outputs are for software verification only, never clinical use.

## Full research build

The full pipeline is intentionally separate from the offline demo. It requires a lawfully assembled literature/source corpus, extracted and reviewed evidence, and any external benchmark inputs under their original access and data-use terms. These data and generated full KG/RAG files are not in this repository. See `docs/data_availability.md` and the relevant skill instructions before running it.

The repository exposes the pipeline stages and entrypoints in `skills/registry.csv` and `README.md`. A full build follows this broad sequence:

1. Search/harvest only through permitted access routes; retain the source manifests and licenses.
2. Parse the locally available corpus with `skills/mrd-evidence-miner/scripts/mrdskill_extract_full_corpus_evidence.py`.
3. Review and score candidates using `skills/mrd-confidence-rater/scripts/mrdskill_score_mrd_evidence.py`.
4. Build provenance-linked graph rows with `skills/mrd-kg-builder/scripts/mrdskill_build_mrd_kg_csv.py`.
5. Build the full literature/RAG and any separately governed supplementary-rule layer using the `scripts/` and `skills/mrd-graphrag/` entrypoints; audit outputs and their source manifests.
6. For VCF QC or source-tracing benchmark reproduction, supply the relevant user-authorized VCF/cohort or original study supplementary input explicitly. The Razavi benchmark builder defaults to a local supplementary XLSX that is not included; obtain and cite it lawfully before invoking `skills/mrd-vcf-qc/scripts/reproduce_razavi_benchmark.py`.

Record exact input versions and checksums, command lines, Python/dependency versions, run date, parser errors, excluded sources and reviewer decisions with each research build. Never replace missing evidence with demo fixtures, and never describe the offline demo as reproducing full-corpus or manuscript benchmark metrics.
