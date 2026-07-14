from __future__ import annotations

import argparse
import hashlib
import shutil
import sys
from pathlib import Path
from typing import Any

import pandas as pd

SKILL_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = SKILL_ROOT / "literature_harvest" / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import harvest_utils as utils  # noqa: E402


PAPER_DIR_NAMES = ["downloaded_pdfs", "downloaded_pdfs_deduplicated"]
SUPPLEMENTARY_DIR_NAMES = ["supplementary_materials", "supplementary_materials_deduplicated"]


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def write_csv(path: Path, df: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, encoding="utf-8-sig")


def record_url_from_log(log_df: pd.DataFrame, record_id: str) -> str:
    if log_df.empty or "record_id" not in log_df.columns:
        return ""
    rows = log_df.loc[log_df["record_id"].astype(str).eq(record_id)]
    if rows.empty:
        return ""
    return utils.clean_text(rows.iloc[-1].get("final_pdf_url"))


def doi_from_candidates(candidate_df: pd.DataFrame, record_id: str) -> str:
    if candidate_df.empty or "record_id" not in candidate_df.columns:
        return ""
    rows = candidate_df.loc[candidate_df["record_id"].astype(str).eq(record_id)]
    if rows.empty:
        return ""
    return utils.clean_text(rows.iloc[0].get("doi"))


def second_pass_url(second_df: pd.DataFrame, path: Path) -> str:
    if second_df.empty or "secondary_pdf_path" not in second_df.columns:
        return ""
    rows = second_df.loc[second_df["secondary_pdf_path"].fillna("").astype(str).eq(str(path))]
    if rows.empty:
        return ""
    return utils.clean_text(rows.iloc[-1].get("selected_pdf_url"))


def move_to_quarantine(path: Path, quarantine_dir: Path | None) -> Path:
    if not path.exists():
        return path
    if quarantine_dir is None:
        path.unlink(missing_ok=True)
        return path
    quarantine_dir.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha1(str(path).encode("utf-8", errors="ignore")).hexdigest()[:12]
    record_id = utils.safe_filename(utils.record_id_from_filename(path.name), max_len=24)
    destination = quarantine_dir / f"{record_id}__{digest}{path.suffix.lower()}"
    counter = 2
    while destination.exists():
        destination = quarantine_dir / f"{record_id}__{digest}__{counter}{path.suffix.lower()}"
        counter += 1
    shutil.move(str(path), str(destination))
    return destination


def clean_paper_file(
    path: Path,
    log_df: pd.DataFrame,
    second_df: pd.DataFrame,
    candidate_df: pd.DataFrame,
    quarantine_dir: Path | None,
) -> dict[str, Any] | None:
    suffix = path.suffix.lower()
    if not path.exists():
        return None
    record_id = utils.record_id_from_filename(path.name)
    reason = ""
    source_url = ""
    if suffix in {".html", ".htm"}:
        source_url = record_url_from_log(log_df, record_id)
        payload = path.read_bytes()
        if not utils.is_likely_article_html(payload, source_url):
            reason = "non_article_html"
    elif suffix == ".pdf":
        payload = path.read_bytes()[:5]
        if payload != b"%PDF-":
            reason = "invalid_pdf_signature"
        elif path.stem.endswith("__secondary"):
            source_url = second_pass_url(second_df, path)
            doi = doi_from_candidates(candidate_df, record_id)
            if not utils.is_likely_article_pdf_url(source_url, "", doi):
                reason = "secondary_pdf_not_article"
    elif suffix not in {".xml"}:
        reason = "unsupported_paper_extension"

    if not reason:
        return None
    destination = move_to_quarantine(path, quarantine_dir)
    return {
        "record_id": record_id,
        "file_type": "paper",
        "source_path": str(path),
        "quarantine_path": str(destination) if quarantine_dir else "",
        "reason": reason,
        "source_url": source_url,
    }


def clean_supplementary_file(path: Path, quarantine_dir: Path | None) -> dict[str, Any] | None:
    suffix = path.suffix.lower()
    if not path.exists():
        return None
    reason = ""
    if suffix in {".html", ".htm"}:
        reason = "supplement_html"
    elif suffix not in utils.SUPPLEMENTARY_TEXT_TABLE_EXTENSIONS:
        reason = "unsupported_supplement_extension"
    elif suffix == ".pdf" and path.read_bytes()[:5] != b"%PDF-":
        reason = "invalid_pdf_signature"
    elif not utils.is_supplementary_attachment(str(path), path.name):
        reason = "non_text_table_supplement"

    if not reason:
        return None
    destination = move_to_quarantine(path, quarantine_dir)
    return {
        "record_id": utils.record_id_from_filename(path.name),
        "file_type": "supplementary",
        "source_path": str(path),
        "quarantine_path": str(destination) if quarantine_dir else "",
        "reason": reason,
        "source_url": "",
    }


def update_logs(run_root: Path, removed_paths: set[str]) -> None:
    log_path = run_root / "download_logs" / "keyword_research_download_log.csv"
    log_df = read_csv(log_path)
    if not log_df.empty and "final_pdf_path" in log_df.columns:
        mask = log_df["final_pdf_path"].fillna("").astype(str).isin(removed_paths)
        if mask.any():
            log_df.loc[mask, "download_status"] = "metadata_only"
            log_df.loc[mask, "failure_reason"] = "cleaned_invalid_saved_file"
            log_df.loc[mask, "final_pdf_path"] = ""
            log_df.loc[mask, "content_format"] = ""
            write_csv(log_path, log_df)

    supp_log_path = run_root / "download_logs" / "keyword_research_supplementary_download_log.csv"
    supp_df = read_csv(supp_log_path)
    if not supp_df.empty and "local_file_path" in supp_df.columns:
        mask = supp_df["local_file_path"].fillna("").astype(str).isin(removed_paths)
        if mask.any():
            supp_df.loc[mask, "download_status"] = "metadata_only"
            supp_df.loc[mask, "failure_reason"] = "cleaned_invalid_supplement"
            supp_df.loc[mask, "local_file_path"] = ""
            supp_df.loc[mask, "content_format"] = ""
            write_csv(supp_log_path, supp_df)

    second_path = run_root / "download_logs" / "keyword_research_html_second_pass.csv"
    second_df = read_csv(second_path)
    if not second_df.empty and "secondary_pdf_path" in second_df.columns:
        valid_mask = (
            second_df["second_pass_status"].fillna("").eq("pdf_downloaded")
            & second_df["secondary_pdf_path"].fillna("").map(lambda value: bool(value) and Path(value).exists())
        )
        cleaned_df = second_df.loc[valid_mask].copy()
        if len(cleaned_df) != len(second_df):
            write_csv(second_path, cleaned_df)


def rebuild_manifest(directory: Path, output_path: Path) -> None:
    rows = []
    if directory.exists():
        for path in directory.iterdir():
            if path.is_file() and not path.name.startswith("~$"):
                rows.append(
                    {
                        "record_id": utils.record_id_from_filename(path.name),
                        "source_path": str(path.resolve().relative_to(output_path.parent.resolve())),
                        "filename": path.name,
                        "extension": path.suffix.lower(),
                        "sha256": file_sha256(path),
                    }
                )
    write_csv(output_path, pd.DataFrame(rows))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True)
    parser.add_argument("--delete", action="store_true", help="Delete invalid files instead of moving them to quarantine.")
    args = parser.parse_args()

    run_root = Path(args.run_root).resolve()
    candidate_df = read_csv(run_root / "keyword_research_candidate_table.csv")
    log_df = read_csv(run_root / "download_logs" / "keyword_research_download_log.csv")
    second_df = read_csv(run_root / "download_logs" / "keyword_research_html_second_pass.csv")
    quarantine_dir = None if args.delete else run_root / "quarantine_invalid_downloads"

    removed_rows: list[dict[str, Any]] = []
    for dirname in PAPER_DIR_NAMES:
        directory = run_root / dirname
        if not directory.exists():
            continue
        for path in list(directory.iterdir()):
            if path.is_file():
                row = clean_paper_file(path, log_df, second_df, candidate_df, quarantine_dir)
                if row:
                    row["directory"] = dirname
                    removed_rows.append(row)

    for dirname in SUPPLEMENTARY_DIR_NAMES:
        directory = run_root / dirname
        if not directory.exists():
            continue
        for path in list(directory.iterdir()):
            if path.is_file():
                if path.name.startswith("~$"):
                    continue
                row = clean_supplementary_file(path, quarantine_dir)
                if row:
                    row["directory"] = dirname
                    removed_rows.append(row)

    removed_paths = {row["source_path"] for row in removed_rows}
    update_logs(run_root, removed_paths)
    rebuild_manifest(run_root / "downloaded_pdfs", run_root / "papers_manifest.csv")
    rebuild_manifest(run_root / "supplementary_materials", run_root / "supplements_manifest.csv")

    report_path = run_root / "download_logs" / "cleaned_invalid_downloads.csv"
    if removed_rows:
        write_csv(report_path, pd.DataFrame(removed_rows))
    print(f"Cleaned invalid files: {len(removed_rows)}")
    if removed_rows:
        print(f"Report: {report_path}")


if __name__ == "__main__":
    main()
