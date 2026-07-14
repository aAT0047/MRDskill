from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts"))

from mrdskill_runner import run_target


if __name__ == "__main__":
    run_target("scripts/mrdskill_augment_kg_with_cohort_rules.py")
