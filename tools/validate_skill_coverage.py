#!/usr/bin/env python3
"""Verify that release implementations are owned by the ten manuscript skills."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


EXPECTED_SKILLS = {
    "paper-search",
    "literature-harvest",
    "mrd-evidence-miner",
    "mrd-confidence-rater",
    "mrd-kg-builder",
    "supplementary-rule",
    "mrd-graphrag",
    "mrd-vcf-qc",
    "review-aware-adjudication",
    "mrd-update-audit",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate MRDskill implementation ownership.")
    parser.add_argument("--quiet", action="store_true", help="Print only validation errors.")
    return parser.parse_args()


def source_files(root: Path) -> set[str]:
    locations = [root / "scripts", root / "literature_harvest" / "scripts"]
    return {
        path.relative_to(root).as_posix()
        for location in locations
        for path in location.rglob("*.py")
        if "__pycache__" not in path.parts
    }


def main() -> int:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    registry_path = root / "skills" / "registry.csv"
    mapping_path = root / "skills" / "implementation_map.csv"
    errors: list[str] = []

    with registry_path.open(encoding="utf-8-sig", newline="") as handle:
        registry = list(csv.DictReader(handle))
    registry_skills = {row["skill_name"].strip() for row in registry if row.get("skill_name")}
    if registry_skills != EXPECTED_SKILLS:
        errors.append(
            "registry skill set differs from the ten supplement skills: "
            f"expected={sorted(EXPECTED_SKILLS)} actual={sorted(registry_skills)}"
        )
    if len(registry) != 10:
        errors.append(f"registry must contain exactly 10 skills, found {len(registry)}")
    for skill in EXPECTED_SKILLS:
        if not (root / "skills" / skill / "SKILL.md").is_file():
            errors.append(f"missing SKILL.md for {skill}")

    with mapping_path.open(encoding="utf-8-sig", newline="") as handle:
        mappings = list(csv.DictReader(handle))
    mapped_sources = {row["source_path"].strip() for row in mappings if row.get("source_path")}
    actual_sources = source_files(root)
    for path in sorted(actual_sources - mapped_sources):
        errors.append(f"unowned implementation source: {path}")
    for path in sorted(mapped_sources - actual_sources):
        errors.append(f"mapping points to missing implementation source: {path}")

    for row in mappings:
        source = row.get("source_path", "").strip()
        owners = [value.strip() for value in row.get("skill_owners", "").split(";") if value.strip()]
        if not owners:
            errors.append(f"missing owner for {source}")
        invalid_owners = [owner for owner in owners if owner not in EXPECTED_SKILLS and owner != "shared"]
        if invalid_owners:
            errors.append(f"invalid owner for {source}: {', '.join(invalid_owners)}")
        launchers = [value.strip() for value in row.get("launcher", "").split(";") if value.strip()]
        for launcher in launchers:
            launcher_path = root / launcher
            if not launcher_path.is_file():
                errors.append(f"missing launcher for {source}: {launcher}")
            elif not launcher.startswith("skills/"):
                errors.append(f"launcher must be packaged in a skill: {launcher}")

    if errors:
        print("MRDskill coverage validation failed:", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 1
    if not args.quiet:
        print(
            "MRDskill coverage validation passed: "
            f"{len(registry)} skills, {len(mapped_sources)} implementation sources."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
