from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts"))

from mrdskill_runner import run_target


if __name__ == "__main__":
    run_target("scripts/mrd_convert_pdfs_with_opendataloader.py")
