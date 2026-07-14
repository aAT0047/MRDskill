#!/usr/bin/env python3
"""Portable paths for the MRDskill source package."""

from __future__ import annotations

import os
import re
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROJECT_CODE_ROOT = PROJECT_ROOT
SKILLS_ROOT = PROJECT_ROOT / "skills"
DEFAULT_ARTIFACT_ROOT = PROJECT_ROOT / "artifacts"


def artifact_root() -> Path:
    """Return the writable artifact root, overridable without editing source."""
    configured = os.environ.get("MRDSKILL_ARTIFACT_ROOT") or os.environ.get("MRDSKILL_DATA_ROOT")
    return Path(configured).expanduser() if configured else DEFAULT_ARTIFACT_ROOT


def code_backup_root() -> Path:
    """Compatibility alias for code that formerly used a separate backup root."""
    configured = os.environ.get("MRDSKILL_CODE_ROOT")
    return Path(configured).expanduser() if configured else PROJECT_CODE_ROOT


def artifact_path(value: str | Path) -> Path:
    """Resolve data/docs/logs paths against the MRDskill artifact root."""
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    if path.parts and path.parts[0].lower() in {"cases", "data", "docs", "imports", "logs"}:
        return artifact_root() / path
    return path


def safe_path_component(value: str | Path | None, fallback: str = "unknown_case") -> str:
    """Return a filesystem-safe folder/file component while preserving readable IDs."""
    text = "" if value is None else str(value).strip()
    if not text:
        text = fallback
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", text)
    text = re.sub(r"_+", "_", text).strip(" ._")
    return text or fallback


def qc_case_dir(sample_id: str | Path | None, case_id: str | Path | None = None) -> Path:
    """Folder that keeps all QC artifacts for one biological sample."""
    folder_name = safe_path_component(sample_id, safe_path_component(case_id, "unknown_case"))
    return artifact_root() / "cases" / folder_name


def script_path(script_name: str) -> Path:
    return PROJECT_CODE_ROOT / "scripts" / script_name


def skill_script_path(skill_name: str, script_name: str) -> Path:
    """Return a bundled skill script without relying on an external data directory."""
    return SKILLS_ROOT / skill_name / "scripts" / script_name


def display_path(path: str | Path) -> str:
    return str(Path(path))
