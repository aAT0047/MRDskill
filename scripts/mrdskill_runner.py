"""Run a project script from a bundled MRDskill wrapper."""

from __future__ import annotations

import runpy
import sys
from pathlib import Path
from typing import Sequence


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def run_target(relative_script: str, extra_args: Sequence[str] = ()) -> None:
    target = PROJECT_ROOT / relative_script
    if not target.is_file():
        raise FileNotFoundError(f"Bundled MRDskill script not found: {target}")
    sys.path.insert(0, str(target.parent))
    sys.argv = [str(target), *sys.argv[1:], *extra_args]
    runpy.run_path(str(target), run_name="__main__")
