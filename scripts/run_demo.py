#!/usr/bin/env python3
"""Run the bundled, offline synthetic MRDskill KG + multi-round RAG demo."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEMO = ROOT / "examples" / "demo"
EXPECTED = ROOT / "examples" / "expected"
KG_BUILDER = ROOT / "skills" / "mrd-kg-builder" / "scripts" / "mrdskill_build_mrd_kg_csv.py"
UNIFIED_INDEX = ROOT / "scripts" / "mrdskill_build_unified_rag_index.py"
ADJUDICATION = ROOT / "scripts" / "mrdskill_run_adjudication_rag.py"
VISUALIZATION = ROOT / "scripts" / "mrdskill_build_mrd_kg_visualization.py"
STABLE_FILES = (
    "nodes.csv",
    "relationships.csv",
    "review_tasks.csv",
    "evidence_summaries.jsonl",
    "rag_chunks.jsonl",
    "rag_retrieval.json",
    "summary.json",
    "kg_visualization_data.json",
    "kg_demo.html",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(command: list[str], env: dict[str, str]) -> None:
    completed = subprocess.run(command, cwd=ROOT, env=env, capture_output=True, text=True, encoding="utf-8", errors="replace")
    if completed.returncode:
        sys.stderr.write(completed.stdout)
        sys.stderr.write(completed.stderr)
        raise subprocess.CalledProcessError(completed.returncode, command)
    if "mrdskill_run_adjudication_rag.py" not in command[1]:
        sys.stdout.write(completed.stdout)


def csv_rows(path: Path) -> int:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return max(0, sum(1 for _ in csv.DictReader(handle)))


def canonical_retrieval(source: Path, target: Path) -> dict[str, Any]:
    result = json.loads(source.read_text(encoding="utf-8"))
    rounds = []
    for item in result.get("rounds", []):
        rounds.append(
            {
                key: item[key]
                for key in (
                    "round_id",
                    "agent",
                    "purpose",
                    "candidate_chunks_scored",
                    "hits",
                    "rule_hits",
                )
                if key in item
            }
        )
    canonical = {
        "scope": "synthetic software fixture; not biological evidence or benchmark",
        "case": result.get("case", {}),
        "missing_information_checklist": result.get("missing_information_checklist", []),
        "rounds": rounds,
    }
    target.write_text(json.dumps(canonical, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return canonical


def make_outputs(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="mrdskill-demo-") as temporary:
        work = Path(temporary)
        env = os.environ.copy()
        env["PYTHONUTF8"] = "1"

        nodes = output_dir / "nodes.csv"
        relationships = output_dir / "relationships.csv"
        review = output_dir / "review_tasks.csv"
        summaries = output_dir / "evidence_summaries.jsonl"
        rag_chunks = output_dir / "rag_chunks.jsonl"
        rag_report = work / "rag_index_report.md"
        raw_retrieval = work / "rag_retrieval.full.json"

        run(
            [
                sys.executable,
                str(KG_BUILDER),
                "--input-csv",
                str(DEMO / "scored_evidence.csv"),
                "--nodes-csv",
                str(nodes),
                "--relationships-csv",
                str(relationships),
                "--review-tasks-csv",
                str(review),
                "--evidence-summaries-jsonl",
                str(summaries),
            ],
            env,
        )

        run(
            [
                sys.executable,
                str(UNIFIED_INDEX),
                "--evidence-jsonl",
                str(summaries),
                "--rule-jsonl",
                str(DEMO / "rule_chunks.jsonl"),
                "--out-jsonl",
                str(rag_chunks),
                "--out-report-md",
                str(rag_report),
            ],
            env,
        )
        run(
            [
                sys.executable,
                str(ADJUDICATION),
                "--case-json",
                str(DEMO / "case.json"),
                "--chunks-jsonl",
                str(rag_chunks),
                "--rule-db",
                str(work / "not-used.sqlite"),
                "--out-json",
                str(raw_retrieval),
            ],
            env,
        )
        canonical = canonical_retrieval(raw_retrieval, output_dir / "rag_retrieval.json")

        graph_json = output_dir / "kg_visualization_data.json"
        graph_html = output_dir / "kg_demo.html"
        run(
            [
                sys.executable,
                str(VISUALIZATION),
                "--root",
                str(work),
                "--raw-csv",
                str(DEMO / "scored_evidence.csv"),
                "--scored-csv",
                str(DEMO / "scored_evidence.csv"),
                "--nodes-csv",
                str(nodes),
                "--relationships-csv",
                str(relationships),
                "--review-csv",
                str(review),
                "--rag-jsonl",
                str(rag_chunks),
                "--sources-csv",
                str(DEMO / "source_manifest.csv"),
                "--chunks-csv",
                str(DEMO / "excel_chunks_manifest.csv"),
                "--out-json",
                str(graph_json),
                "--out-html",
                str(graph_html),
            ],
            env,
        )

        sys.path.insert(0, str(ROOT / "scripts"))
        import mrdskill_build_mrd_kg_visualization as kg_visualizer

        graph_data = json.loads(graph_json.read_text(encoding="utf-8"))
        graph_data["title"] = "MRDskill Synthetic Demo KG"
        graph_data["generated_from"] = {
            "raw_csv": "examples/demo/scored_evidence.csv",
            "scored_csv": "examples/demo/scored_evidence.csv",
            "sources_csv": "examples/demo/source_manifest.csv",
            "chunks_csv": "examples/demo/excel_chunks_manifest.csv",
            "nodes_csv": "nodes.csv in selected output directory",
            "relationships_csv": "relationships.csv in selected output directory",
            "review_csv": "review_tasks.csv in selected output directory",
            "rag_jsonl": "rag_chunks.jsonl in selected output directory",
            "baseline_csv": "not supplied by the synthetic demo",
            "cooccurrence_csv": "not supplied by the synthetic demo",
            "mutex_csv": "not supplied by the synthetic demo",
        }
        kg_visualizer.HTML_TEMPLATE = (
            kg_visualizer.HTML_TEMPLATE
            .replace("MRD Full Corpus Support KG", "MRDskill Synthetic Demo KG")
            .replace("Owner Phd.W - candidate support graph - full-corpus build", "Synthetic software fixture - not evidence or clinical output")
        )
        kg_visualizer.write_json(graph_json, graph_data)
        kg_visualizer.write_html(graph_html, graph_data)

        summary = {
            "scope": "synthetic software fixture; not biological evidence or benchmark",
            "kg_nodes": csv_rows(nodes),
            "kg_relationships": csv_rows(relationships),
            "review_tasks": csv_rows(review),
            "rag_chunks": sum(1 for line in rag_chunks.read_text(encoding="utf-8").splitlines() if line.strip()),
            "rag_rounds": len(canonical["rounds"]),
            "rag_rounds_with_hits": sum(bool(item.get("hits")) for item in canonical["rounds"]),
            "patient_data_included": False,
            "external_literature_included": False,
            "clinical_use": False,
        }
        (output_dir / "summary.json").write_text(
            json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )

        inputs = [
            DEMO / name
            for name in (
                "README.md",
                "case.json",
                "scored_evidence.csv",
                "rule_chunks.jsonl",
                "source_manifest.csv",
                "excel_chunks_manifest.csv",
                "manifest.json",
            )
        ]
        output_names = list(STABLE_FILES)
        manifest = {
            "demo_version": "1.0.0",
            "base_commit_before_demo_changes": subprocess.run(
                ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, capture_output=True, text=True
            ).stdout.strip(),
            "working_tree_demo_changes_included": True,
            "python": sys.version.split()[0],
            "generated_by": "scripts/run_demo.py",
            "scope": "Synthetic offline demo only. Not the complete MRDskill corpus, KG, RAG index, or benchmark.",
            "demo_fixture_license": "CC0-1.0; see examples/demo/README.md",
            "implementation_sha256": {
                str(path.relative_to(ROOT)).replace("\\", "/"): sha256(path)
                for path in (Path(__file__).resolve(), KG_BUILDER, UNIFIED_INDEX, ADJUDICATION, VISUALIZATION)
            },
            "inputs_sha256": {str(path.relative_to(ROOT)).replace("\\", "/"): sha256(path) for path in inputs},
            "outputs_sha256": {name: sha256(output_dir / name) for name in output_names},
        }
        (output_dir / "manifest.json").write_text(
            json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )


def verify_expected(actual_dir: Path) -> None:
    failures = []
    actual_manifest = json.loads((actual_dir / "manifest.json").read_text(encoding="utf-8"))
    expected_manifest = json.loads((EXPECTED / "manifest.json").read_text(encoding="utf-8"))
    for name in STABLE_FILES:
        actual = actual_dir / name
        expected = EXPECTED / name
        if not expected.is_file():
            failures.append(f"missing committed expected output: {expected.relative_to(ROOT)}")
        elif sha256(actual) != sha256(expected):
            failures.append(f"output differs from precomputed fixture: {name}")
    for section in ("inputs_sha256", "implementation_sha256"):
        if actual_manifest.get(section) != expected_manifest.get(section):
            failures.append(f"precomputed {section} do not match current demo inputs/code")
        for name, digest in actual_manifest.get(section, {}).items():
            path = ROOT / name
            if not path.is_file() or sha256(path) != digest:
                failures.append(f"{section} checksum mismatch: {name}")
    for name, digest in actual_manifest.get("outputs_sha256", {}).items():
        path = actual_dir / name
        if not path.is_file() or sha256(path) != digest:
            failures.append(f"generated output checksum mismatch: {name}")
    if failures:
        raise SystemExit("Demo verification failed:\n- " + "\n- ".join(failures))
    print("PASS: deterministic KG/RAG demo outputs match examples/expected.")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "artifacts" / "demo",
        help="Output directory (default: artifacts/demo/).",
    )
    parser.add_argument(
        "--verify-expected",
        action="store_true",
        help="Also compare deterministic outputs with committed precomputed examples.",
    )
    args = parser.parse_args()
    output_dir = args.output_dir.expanduser().resolve()
    make_outputs(output_dir)
    print(f"Synthetic MRDskill demo outputs: {output_dir}")
    print(f"Interactive KG preview: {output_dir / 'kg_demo.html'}")
    if args.verify_expected:
        verify_expected(output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
