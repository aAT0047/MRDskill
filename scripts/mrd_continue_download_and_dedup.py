from __future__ import annotations

import argparse
import hashlib
import html
import json
import re
import shutil
import ssl
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

import pandas as pd

SKILL_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = SKILL_ROOT / "literature_harvest" / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import harvest_utils as utils  # noqa: E402


USER_AGENT = "keyword-research-harvest/0.1"
PAYWALL_MARKERS = [
    "subscription required",
    "purchase this article",
    "access through your institution",
    "institutional access",
    "paywall",
    "buy article",
    "rent this article",
]


def clean_text(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and pd.isna(value):
        return ""
    text = html.unescape(str(value))
    text = re.sub(r"<[^>]+>", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def normalize_doi(value: Any) -> str:
    text = clean_text(value).lower()
    text = text.replace("https://doi.org/", "").replace("http://doi.org/", "")
    text = text.replace("doi:", "").strip().rstrip(".")
    return text


def normalize_title(value: Any) -> str:
    text = clean_text(value).lower()
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def content_extension(content_type: str, url: str) -> str:
    lowered_type = clean_text(content_type).lower()
    lowered_url = url.lower()
    if "pdf" in lowered_type or lowered_url.endswith(".pdf"):
        return ".pdf"
    if "xml" in lowered_type or lowered_url.endswith(".xml"):
        return ".xml"
    if "html" in lowered_type or lowered_type.startswith("text/"):
        return ".html"
    return ".bin"


def classify_payload(payload: bytes, content_type: str) -> str:
    sample = payload[:3000].decode("utf-8", errors="ignore").lower()
    if payload.startswith(b"%PDF") or "pdf" in content_type.lower():
        return "pdf"
    if "xml" in content_type.lower() or sample.lstrip().startswith("<?xml"):
        return "xml"
    if "html" in content_type.lower() or "<html" in sample:
        return "html"
    return "other"


def looks_paywalled(payload: bytes) -> bool:
    text = payload[:50000].decode("utf-8", errors="ignore").lower()
    return any(marker in text for marker in PAYWALL_MARKERS)


def request_url(url: str, timeout: int) -> dict[str, Any]:
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": USER_AGENT,
            "Accept": "text/html, application/xml, text/xml;q=0.95, application/pdf;q=0.6, */*;q=0.1",
        },
    )
    with urllib.request.urlopen(request, timeout=timeout, context=ssl.create_default_context()) as response:
        return {
            "payload": response.read(),
            "content_type": clean_text(response.headers.get("Content-Type")).lower(),
            "final_url": clean_text(response.geturl()) or url,
        }


def load_run_config(run_root: Path) -> dict[str, Any]:
    for candidate in [
        run_root / "config.json",
        run_root / "run_config.json",
        run_root.parent.parent / "references" / "demo_10_config.json",
        run_root.parent.parent / "references" / "config_template.json",
    ]:
        if candidate.exists():
            try:
                return json.loads(candidate.read_text(encoding="utf-8"))
            except Exception:  # noqa: BLE001
                continue
    return {}


def enrich_row(row: pd.Series, config: dict[str, Any]) -> pd.Series:
    return pd.Series(utils.append_external_links_to_row(row, config))


def candidate_urls(row: pd.Series) -> list[str]:
    urls: list[str] = []
    for key in ["fulltext_url", "landing_page_url", "pdf_url_candidate", "pdf_url"]:
        value = clean_text(row.get(key))
        if value.startswith("http://") or value.startswith("https://"):
            urls.append(value)
    raw_json = clean_text(row.get("candidate_urls_json"))
    if raw_json:
        try:
            parsed = json.loads(raw_json)
            if isinstance(parsed, list):
                for item in parsed:
                    if isinstance(item, str):
                        value = clean_text(item)
                    elif isinstance(item, dict):
                        value = clean_text(item.get("url"))
                    else:
                        value = ""
                    if value.startswith("http://") or value.startswith("https://"):
                        urls.append(value)
        except json.JSONDecodeError:
            pass
    doi = normalize_doi(row.get("doi"))
    if doi:
        urls.append(f"https://doi.org/{doi}")
    unique: list[str] = []
    seen: set[str] = set()
    for url in sorted(urls, key=utils.article_fulltext_preference_rank):
        if url not in seen:
            seen.add(url)
            unique.append(url)
    return unique


def supplementary_urls(row: pd.Series) -> list[dict[str, str]]:
    items = utils.parse_supplementary_urls(row.get("supplementary_urls_json"))
    raw_json = clean_text(row.get("candidate_urls_json"))
    if raw_json:
        try:
            parsed = json.loads(raw_json)
        except json.JSONDecodeError:
            parsed = []
        for item in parsed if isinstance(parsed, list) else []:
            if isinstance(item, str):
                url = clean_text(item)
                label = ""
            elif isinstance(item, dict):
                url = clean_text(item.get("url"))
                label = clean_text(item.get("label") or item.get("style") or item.get("type"))
            else:
                continue
            if url and utils.is_supplementary_attachment(url, label):
                items.append({"url": url, "label": label or "supplementary_material", "source": "candidate_urls"})
    return utils.unique_url_items(items)


def discover_supplementary_from_saved_file(path: Path, row: pd.Series) -> list[dict[str, str]]:
    if path.suffix.lower() not in {".html", ".htm", ".xml"}:
        return []
    text = path.read_text(encoding="utf-8", errors="ignore")
    return utils.extract_supplementary_links(text, candidate_urls(row))


def discover_supplementary_from_candidate_pages(
    row: pd.Series,
    timeout: int,
    max_pages: int,
    exclude_urls: set[str] | None = None,
) -> list[dict[str, str]]:
    if max_pages <= 0:
        return []
    exclude = {url.lower() for url in (exclude_urls or set()) if url}
    discovered: list[dict[str, str]] = []
    checked = 0
    for url in candidate_urls(row):
        if url.lower() in exclude:
            continue
        lowered = url.lower()
        if lowered.endswith(".pdf") or "/pdf" in lowered or "pdf=render" in lowered:
            continue
        try:
            response = request_url(url, timeout)
            payload = response["payload"]
            final_url = response["final_url"] or url
            kind = classify_payload(payload, response["content_type"])
            if kind not in {"html", "xml"}:
                continue
            if kind == "html":
                if looks_paywalled(payload):
                    continue
                if not utils.is_likely_article_html(payload, final_url):
                    continue
            discovered.extend(utils.extract_supplementary_links(payload.decode("utf-8", errors="ignore"), [final_url, url]))
            checked += 1
            if checked >= max_pages:
                break
        except Exception:  # noqa: BLE001
            continue
        finally:
            time.sleep(0.1)
    return utils.unique_url_items(discovered)


def download_supplementary_files(
    row: pd.Series,
    supp_dir: Path,
    timeout: int,
    max_items: int,
    discovered_items: list[dict[str, str]] | None = None,
    already_downloaded_urls: set[str] | None = None,
) -> list[dict[str, Any]]:
    supp_dir.mkdir(parents=True, exist_ok=True)
    items = supplementary_urls(row)
    if discovered_items:
        items.extend(discovered_items)
    items = utils.unique_url_items(items)
    items = sorted(items, key=utils.supplementary_download_preference_rank)
    done_urls = {url.lower() for url in (already_downloaded_urls or set()) if url}
    items = [item for item in items if item["url"].lower() not in done_urls]
    if max_items > 0:
        items = items[:max_items]

    rows: list[dict[str, Any]] = []
    for index, item in enumerate(items, start=1):
        url = item["url"]
        label = item.get("label") or "supplementary_material"
        status = "metadata_only"
        reason = ""
        final_url = url
        output_path = ""
        content_format = ""
        try:
            response = request_url(url, timeout)
            payload = response["payload"]
            final_url = response["final_url"] or url
            kind = classify_payload(payload, response["content_type"])
            if kind == "html" and looks_paywalled(payload):
                status = "inaccessible"
                reason = "paywall_detected"
            else:
                extension = utils.attachment_extension(response["content_type"], final_url, payload)
                if kind == "html" or extension in {".html", ".bin"}:
                    status = "metadata_only"
                    reason = "supplement_html_skipped" if extension == ".html" else "unknown_binary_skipped"
                    rows.append(
                        {
                            "record_id": row.get("record_id"),
                            "supplement_index": index,
                            "doi": row.get("doi"),
                            "title": row.get("title"),
                            "source_url": url,
                            "final_url": final_url,
                            "download_status": status,
                            "failure_reason": reason,
                            "content_format": kind,
                            "local_file_path": output_path,
                            "label": label,
                        }
                    )
                    continue
                output = supp_dir / f"{utils.supplementary_file_stem(row, index, final_url, label)}{extension}"
                output.write_bytes(payload)
                status = "success"
                content_format = kind
                output_path = str(output)
        except urllib.error.HTTPError as exc:
            code = int(getattr(exc, "code", 0))
            status = "inaccessible" if code in {401, 402, 403} else ("broken_link" if code in {404, 410} else "metadata_only")
            if code == 429:
                status = "rate_limited"
            reason = f"http_{code}"
        except Exception as exc:  # noqa: BLE001
            status = "metadata_only"
            reason = type(exc).__name__
        finally:
            time.sleep(0.1)
        rows.append(
            {
                "record_id": row.get("record_id"),
                "supplement_index": index,
                "doi": row.get("doi"),
                "title": row.get("title"),
                "source_url": url,
                "final_url": final_url,
                "download_status": status,
                "failure_reason": reason,
                "content_format": content_format,
                "local_file_path": output_path,
                "label": label,
            }
        )
    return rows


def download_pmc_oa_packages(
    row: pd.Series,
    package_dir: Path,
    extract_root: Path,
    pdf_dir: Path,
    supp_dir: Path,
    timeout: int,
) -> list[dict[str, Any]]:
    package_dir.mkdir(parents=True, exist_ok=True)
    extract_root.mkdir(parents=True, exist_ok=True)
    pdf_dir.mkdir(parents=True, exist_ok=True)
    supp_dir.mkdir(parents=True, exist_ok=True)
    record_extract_dir = extract_root / utils.safe_filename(clean_text(row.get("record_id")), 40)
    rows: list[dict[str, Any]] = []
    for index, item in enumerate(utils.pmc_oa_package_urls(row), start=1):
        url = item["url"]
        label = item.get("label") or "pmc_oa_package"
        status = "metadata_only"
        reason = ""
        final_url = url
        package_path = ""
        extracted_count = 0
        paper_files: list[str] = []
        supplementary_files: list[str] = []
        try:
            response = request_url(url, timeout)
            payload = response["payload"]
            final_url = response["final_url"] or url
            extension = utils.attachment_extension(response["content_type"], final_url, payload)
            if extension not in {".tar.gz", ".tgz", ".gz"}:
                reason = f"unexpected_package_extension:{extension}"
            else:
                package = package_dir / f"{utils.supplementary_file_stem(row, index, final_url, label)}{extension}"
                package.write_bytes(payload)
                package_path = str(package)
                extracted = utils.safe_extract_tar_members(package, record_extract_dir)
                extracted_count = len(extracted)
                for file_index, path in enumerate(extracted, start=1):
                    suffix = path.suffix.lower()
                    if suffix == ".pdf" and path.read_bytes()[:5] == b"%PDF-":
                        destination = pdf_dir / f"{utils.paper_file_stem(row)}__pmc_oa_package_{file_index:02d}.pdf"
                        shutil.copy2(path, destination)
                        paper_files.append(str(destination))
                    elif suffix in {".nxml", ".xml"}:
                        destination = pdf_dir / f"{utils.paper_file_stem(row)}__pmc_oa_package_{file_index:02d}.xml"
                        shutil.copy2(path, destination)
                        paper_files.append(str(destination))
                    elif suffix in utils.SUPPLEMENTARY_ATTACHMENT_EXTENSIONS and suffix not in {".tar", ".gz", ".tgz", ".tar.gz"}:
                        destination = supp_dir / f"{utils.supplementary_file_stem(row, file_index, final_url, path.name)}{suffix}"
                        shutil.copy2(path, destination)
                        supplementary_files.append(str(destination))
                status = "success"
        except Exception as exc:  # noqa: BLE001
            reason = type(exc).__name__
        finally:
            time.sleep(0.1)
        rows.append(
            {
                "record_id": row.get("record_id"),
                "doi": row.get("doi"),
                "pmcid": row.get("pmcid"),
                "title": row.get("title"),
                "source_url": url,
                "final_url": final_url,
                "download_status": status,
                "failure_reason": reason,
                "package_path": package_path,
                "extracted_file_count": extracted_count,
                "paper_files": json.dumps(paper_files, ensure_ascii=False),
                "supplementary_files": json.dumps(supplementary_files, ensure_ascii=False),
            }
        )
    return rows


def request_row_download(row: pd.Series, pdf_dir: Path, timeout: int, max_attempts: int) -> dict[str, Any]:
    status = "metadata_only"
    reason = "no_candidate_url"
    final_path = ""
    final_url = ""
    content_format = ""
    attempts = 0
    for url in candidate_urls(row)[:max_attempts]:
        attempts += 1
        try:
            response = request_url(url, timeout)
            payload = response["payload"]
            final_url = response["final_url"] or url
            kind = classify_payload(payload, response["content_type"])
            if kind == "html" and looks_paywalled(payload):
                status = "inaccessible"
                reason = "paywall_detected"
                continue
            if kind == "html" and not utils.is_likely_article_html(payload, final_url):
                status = "metadata_only"
                reason = "non_article_html"
                continue
            extension = ".pdf" if kind == "pdf" else content_extension(response["content_type"], final_url)
            if kind == "html":
                extension = ".html"
            elif kind == "xml":
                extension = ".xml"
            output = pdf_dir / f"{utils.paper_file_stem(row)}{extension}"
            output.write_bytes(payload)
            status = "success"
            reason = ""
            final_path = str(output)
            content_format = kind
            break
        except urllib.error.HTTPError as exc:
            code = int(getattr(exc, "code", 0))
            status = "inaccessible" if code in {401, 402, 403} else ("broken_link" if code in {404, 410} else "metadata_only")
            if code == 429:
                status = "rate_limited"
            reason = f"http_{code}"
            final_url = url
            if code == 429:
                break
            continue
        except Exception as exc:  # noqa: BLE001
            status = "metadata_only"
            reason = type(exc).__name__
            final_url = url
            continue
        finally:
            time.sleep(0.15)
    return {
        "record_id": row.get("record_id"),
        "doi": row.get("doi"),
        "title": row.get("title"),
        "final_pdf_path": final_path,
        "final_pdf_url": final_url,
        "download_status": status,
        "failure_reason": reason,
        "content_format": content_format,
        "attempt_count": attempts,
        "supplementary_url_count": len(supplementary_urls(row)),
    }


def extract_pdf_links(html_text: str, bases: list[str], doi: Any = "") -> list[str]:
    return utils.extract_article_pdf_links(html_text, bases, doi)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def relative_manifest_path(path: Path, run_root: Path) -> str:
    """Store run-relative paths in manifests so they remain portable."""
    try:
        return str(path.resolve().relative_to(run_root.resolve()))
    except ValueError:
        return str(path)


def prepare_download_targets(candidate_df: pd.DataFrame) -> pd.DataFrame:
    data = candidate_df.copy()
    title_text = data.get("title", pd.Series("", index=data.index)).fillna("").astype(str).str.lower()
    publication_text = data.get("publication_type", pd.Series("", index=data.index)).fillna("").astype(str).str.lower()
    abstract_like = (
        title_text.str.contains(r"\babstract\s+(?:lb|ct|p|po|o|e|eps|suppl)?\s*-?\s*\d+", regex=True)
        | title_text.str.contains(r"^\s*(?:lb|ct|p|po|o|e|eps)\s*-?\s*\d+", regex=True)
        | title_text.str.contains(r"^\s*\d{2,4}o\s+", regex=True)
        | publication_text.str.contains("conference abstract", regex=False)
    )
    data = data.loc[~abstract_like].copy()
    if "priority_journal_tier" not in data.columns:
        data["priority_journal_tier"] = data.apply(
            lambda row: utils.priority_journal_tier(row.get("matched_priority_journal") or row.get("journal")),
            axis=1,
        )
    if "dedup_download_key" not in data.columns:
        data["dedup_download_key"] = data.apply(
            lambda row: normalize_doi(row.get("doi"))
            or normalize_title(row.get("title"))
            or clean_text(row.get("record_id")),
            axis=1,
        )
    data["_downloadable_rank"] = (
        data.get("pdf_url", pd.Series("", index=data.index)).fillna("").ne("")
        | data.get("fulltext_url", pd.Series("", index=data.index)).fillna("").ne("")
        | data.get("landing_page_url", pd.Series("", index=data.index)).fillna("").ne("")
        | data.get("candidate_urls_json", pd.Series("", index=data.index)).fillna("").astype(str).str.len().gt(2)
        | data.get("pmc_oa_package_urls_json", pd.Series("", index=data.index)).fillna("").astype(str).str.len().gt(2)
    )
    data["_year_sort"] = pd.to_numeric(data.get("year", pd.Series(pd.NA, index=data.index)), errors="coerce")
    data = data.sort_values(
        [
            "priority_journal_tier",
            "keyword_relevance_score",
            "_downloadable_rank",
            "_year_sort",
        ],
        ascending=[True, False, False, False],
        na_position="last",
    )
    data = data.drop_duplicates(subset=["dedup_download_key"], keep="first")
    return data.drop(columns=["_downloadable_rank", "_year_sort"], errors="ignore").reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True)
    parser.add_argument("--retry-failed", action="store_true")
    parser.add_argument("--timeout", type=int, default=35)
    parser.add_argument("--max-attempts", type=int, default=3)
    parser.add_argument("--max-supplementary", type=int, default=0, help="Maximum supplementary files per record; 0 means no limit")
    parser.add_argument("--max-supplementary-discovery-pages", type=int, default=3)
    parser.add_argument("--max-records", type=int, default=0, help="Maximum records to process in this continuation run; 0 means no limit")
    args = parser.parse_args()

    run_root = Path(args.run_root).resolve()
    candidate_path = run_root / "keyword_research_candidate_table.csv"
    pdf_dir = run_root / "downloaded_pdfs"
    log_dir = run_root / "download_logs"
    log_path = log_dir / "keyword_research_download_log.csv"
    supp_log_path = log_dir / "keyword_research_supplementary_download_log.csv"
    package_log_path = log_dir / "keyword_research_pmc_oa_package_log.csv"
    second_pass_path = log_dir / "keyword_research_html_second_pass.csv"
    dedup_dir = run_root / "downloaded_pdfs_deduplicated"
    supp_dir = run_root / "supplementary_materials"
    supp_dedup_dir = run_root / "supplementary_materials_deduplicated"
    package_dir = run_root / "pmc_oa_packages"
    package_extract_dir = run_root / "pmc_oa_extracted"
    manifest_path = run_root / "papers_manifest.csv"
    supp_manifest_path = run_root / "supplements_manifest.csv"

    pdf_dir.mkdir(parents=True, exist_ok=True)
    supp_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    dedup_dir.mkdir(parents=True, exist_ok=True)
    supp_dedup_dir.mkdir(parents=True, exist_ok=True)

    config = load_run_config(run_root)
    candidate_df = pd.read_csv(candidate_path)
    log_df = pd.read_csv(log_path) if log_path.exists() else pd.DataFrame()
    success_ids = set(log_df.loc[log_df.get("download_status", pd.Series(dtype=str)).eq("success"), "record_id"].astype(str)) if not log_df.empty else set()
    processed_ids = set(log_df.get("record_id", pd.Series(dtype=str)).astype(str)) if not log_df.empty else set()

    eligible_df = prepare_download_targets(
        candidate_df.loc[candidate_df["exclusion_reason_if_any"].fillna("").eq("")].copy()
    )
    if args.retry_failed:
        target = eligible_df.loc[~eligible_df["record_id"].astype(str).isin(success_ids)].copy()
    else:
        target = eligible_df.loc[~eligible_df["record_id"].astype(str).isin(processed_ids)].copy()
    target = utils.apply_year_quota(target, config)
    if args.max_records > 0:
        target = target.head(args.max_records).copy()
    if not target.empty and config.get("external_fulltext", {}).get("enabled", True):
        enriched_target = []
        for _, item in target.iterrows():
            enriched = utils.append_external_links_to_row(item, config)
            enriched_target.append(enriched)
            record_id = clean_text(enriched.get("record_id"))
            if record_id:
                mask = candidate_df["record_id"].astype(str).eq(record_id)
                for column, value in enriched.items():
                    if column not in candidate_df.columns:
                        candidate_df[column] = ""
                    candidate_df.loc[mask, column] = value
        target = pd.DataFrame(enriched_target)
        candidate_df.to_csv(candidate_path, index=False, encoding="utf-8-sig")

    existing_supp_log = pd.read_csv(supp_log_path) if supp_log_path.exists() else pd.DataFrame()
    downloaded_supp_urls = set(
        existing_supp_log.loc[
            existing_supp_log.get("download_status", pd.Series(dtype=str)).eq("success"),
            "source_url",
        ].astype(str).str.lower()
    ) if not existing_supp_log.empty and "source_url" in existing_supp_log.columns else set()

    rows: list[dict[str, Any]] = []
    supp_rows: list[dict[str, Any]] = []
    package_rows: list[dict[str, Any]] = []
    for idx, (_, row) in enumerate(target.iterrows(), start=1):
        result = request_row_download(row, pdf_dir, args.timeout, args.max_attempts)
        rows.append(result)
        discovered = discover_supplementary_from_saved_file(Path(result["final_pdf_path"]), row) if result.get("final_pdf_path") else []
        discovered.extend(
            discover_supplementary_from_candidate_pages(
                row,
                args.timeout,
                args.max_supplementary_discovery_pages,
                {result.get("final_pdf_url", "")},
            )
        )
        discovered = utils.unique_url_items(discovered)
        new_supp_rows = download_supplementary_files(
            row,
            supp_dir,
            args.timeout,
            args.max_supplementary,
            discovered,
            downloaded_supp_urls,
        )
        downloaded_supp_urls.update(
            clean_text(item.get("source_url")).lower()
            for item in new_supp_rows
            if item.get("download_status") == "success"
        )
        supp_rows.extend(new_supp_rows)
        package_rows.extend(
            download_pmc_oa_packages(row, package_dir, package_extract_dir, pdf_dir, supp_dir, args.timeout)
        )
        if idx % 50 == 0:
            chunk = pd.DataFrame(rows)
            if log_path.exists():
                out = pd.concat([pd.read_csv(log_path), chunk], ignore_index=True)
            else:
                out = chunk
            out = out.drop_duplicates(subset=["record_id"], keep="last")
            out.to_csv(log_path, index=False, encoding="utf-8-sig")
            rows = []
            if supp_rows:
                supp_chunk = pd.DataFrame(supp_rows)
                if supp_log_path.exists():
                    supp_out = pd.concat([pd.read_csv(supp_log_path), supp_chunk], ignore_index=True)
                else:
                    supp_out = supp_chunk
                supp_out = supp_out.drop_duplicates(subset=["record_id", "source_url"], keep="last")
                supp_out.to_csv(supp_log_path, index=False, encoding="utf-8-sig")
                supp_rows = []
            if package_rows:
                package_chunk = pd.DataFrame(package_rows)
                if package_log_path.exists():
                    package_out = pd.concat([pd.read_csv(package_log_path), package_chunk], ignore_index=True)
                else:
                    package_out = package_chunk
                package_out = package_out.drop_duplicates(subset=["record_id", "source_url"], keep="last")
                package_out.to_csv(package_log_path, index=False, encoding="utf-8-sig")
                package_rows = []
    if rows:
        if log_path.exists():
            out = pd.concat([pd.read_csv(log_path), pd.DataFrame(rows)], ignore_index=True)
        else:
            out = pd.DataFrame(rows)
        out = out.drop_duplicates(subset=["record_id"], keep="last")
        out.to_csv(log_path, index=False, encoding="utf-8-sig")
    if supp_rows:
        if supp_log_path.exists():
            supp_out = pd.concat([pd.read_csv(supp_log_path), pd.DataFrame(supp_rows)], ignore_index=True)
        else:
            supp_out = pd.DataFrame(supp_rows)
        supp_out = supp_out.drop_duplicates(subset=["record_id", "source_url"], keep="last")
        supp_out.to_csv(supp_log_path, index=False, encoding="utf-8-sig")
    if package_rows:
        if package_log_path.exists():
            package_out = pd.concat([pd.read_csv(package_log_path), pd.DataFrame(package_rows)], ignore_index=True)
        else:
            package_out = pd.DataFrame(package_rows)
        package_out = package_out.drop_duplicates(subset=["record_id", "source_url"], keep="last")
        package_out.to_csv(package_log_path, index=False, encoding="utf-8-sig")

    # HTML second pass
    second_rows: list[dict[str, Any]] = []
    html_files = [p for p in pdf_dir.iterdir() if p.is_file() and p.suffix.lower() in {".html", ".htm"}]
    done_second = set(pd.read_csv(second_pass_path)["record_id"].astype(str)) if second_pass_path.exists() else set()
    meta = candidate_df.set_index("record_id", drop=False)
    for idx, path in enumerate(html_files, start=1):
        record_id = utils.record_id_from_filename(path.name)
        if record_id in done_second:
            continue
        row = meta.loc[record_id] if record_id in meta.index else pd.Series({"record_id": record_id})
        html_text = path.read_text(encoding="utf-8", errors="ignore")
        links = extract_pdf_links(html_text, candidate_urls(row), row.get("doi"))
        status = "no_pdf_link_found"
        out_path = ""
        selected_url = ""
        for link in links[: args.max_attempts]:
            try:
                response = request_url(link, args.timeout)
                if classify_payload(response["payload"], response["content_type"]) == "pdf":
                    dst = pdf_dir / f"{path.stem}__secondary.pdf"
                    dst.write_bytes(response["payload"])
                    status = "pdf_downloaded"
                    out_path = str(dst)
                    selected_url = response["final_url"] or link
                    break
            except Exception:  # noqa: BLE001
                continue
            finally:
                time.sleep(0.1)
        second_rows.append(
            {
                "record_id": record_id,
                "source_html_path": str(path),
                "second_pass_status": status,
                "secondary_pdf_path": out_path,
                "selected_pdf_url": selected_url,
            }
        )
        if idx % 50 == 0:
            chunk = pd.DataFrame(second_rows)
            if second_pass_path.exists():
                out = pd.concat([pd.read_csv(second_pass_path), chunk], ignore_index=True)
            else:
                out = chunk
            out = out.drop_duplicates(subset=["record_id"], keep="last")
            out.to_csv(second_pass_path, index=False, encoding="utf-8-sig")
            second_rows = []
    if second_rows:
        if second_pass_path.exists():
            out = pd.concat([pd.read_csv(second_pass_path), pd.DataFrame(second_rows)], ignore_index=True)
        else:
            out = pd.DataFrame(second_rows)
        out = out.drop_duplicates(subset=["record_id"], keep="last")
        out.to_csv(second_pass_path, index=False, encoding="utf-8-sig")

    # Dedup
    file_rows: list[dict[str, Any]] = []
    meta = candidate_df.set_index("record_id", drop=False)
    for path in pdf_dir.iterdir():
        if not path.is_file():
            continue
        record_id = utils.record_id_from_filename(path.name)
        row = meta.loc[record_id] if record_id in meta.index else pd.Series({})
        file_rows.append(
            {
                "record_id": record_id,
                "source_path": str(path),
                "filename": path.name,
                "extension": path.suffix.lower(),
                "sha256": sha256_file(path),
                "doi": normalize_doi(row.get("doi")),
                "title_normalized": normalize_title(row.get("title")),
            }
        )
    manifest = pd.DataFrame(file_rows)
    if not manifest.empty:
        keep_by_group: dict[str, str] = {}
        for _, row in manifest.iterrows():
            key = f"doi:{row['doi']}" if row["doi"] else (f"title:{row['title_normalized']}" if row["title_normalized"] else f"hash:{row['sha256']}")
            current = keep_by_group.get(key)
            if current is None:
                keep_by_group[key] = row["source_path"]
            else:
                current_path = Path(current)
                new_path = Path(row["source_path"])
                current_rank = (
                    utils.article_file_preference_rank(current_path),
                    -current_path.stat().st_size,
                )
                new_rank = (
                    utils.article_file_preference_rank(new_path),
                    -new_path.stat().st_size,
                )
                if new_rank < current_rank:
                    keep_by_group[key] = row["source_path"]
            manifest.loc[manifest.index == row.name, "dedup_group_key"] = key
        manifest["dedup_decision"] = manifest["source_path"].map(lambda p: "keep" if p in set(keep_by_group.values()) else "duplicate")
        for _, row in manifest.loc[manifest["dedup_decision"].eq("keep")].iterrows():
            src = Path(row["source_path"])
            dst = dedup_dir / src.name
            if not dst.exists():
                shutil.copy2(src, dst)
    manifest_for_output = manifest.copy()
    if not manifest_for_output.empty:
        manifest_for_output["source_path"] = manifest_for_output["source_path"].map(
            lambda value: relative_manifest_path(Path(value), run_root)
        )
    manifest_for_output.to_csv(manifest_path, index=False, encoding="utf-8-sig")

    supp_file_rows: list[dict[str, Any]] = []
    for path in supp_dir.iterdir():
        if not path.is_file():
            continue
        if path.name.startswith("~$"):
            continue
        record_id = utils.record_id_from_filename(path.name)
        row = meta.loc[record_id] if record_id in meta.index else pd.Series({})
        supp_file_rows.append(
            {
                "record_id": record_id,
                "source_path": str(path),
                "filename": path.name,
                "extension": path.suffix.lower(),
                "sha256": sha256_file(path),
                "doi": normalize_doi(row.get("doi")),
                "title_normalized": normalize_title(row.get("title")),
            }
        )
    supp_manifest = pd.DataFrame(supp_file_rows)
    if not supp_manifest.empty:
        supp_manifest["dedup_group_key"] = supp_manifest["sha256"].map(lambda value: f"hash:{value}")
        supp_manifest["dedup_decision"] = "duplicate"
        keep_hashes: set[str] = set()
        for idx, row in supp_manifest.iterrows():
            if row["sha256"] not in keep_hashes:
                keep_hashes.add(row["sha256"])
                supp_manifest.loc[idx, "dedup_decision"] = "keep"
        for _, row in supp_manifest.loc[supp_manifest["dedup_decision"].eq("keep")].iterrows():
            src = Path(row["source_path"])
            dst = supp_dedup_dir / src.name
            if not dst.exists():
                shutil.copy2(src, dst)
    supp_manifest_for_output = supp_manifest.copy()
    if not supp_manifest_for_output.empty:
        supp_manifest_for_output["source_path"] = supp_manifest_for_output["source_path"].map(
            lambda value: relative_manifest_path(Path(value), run_root)
        )
    supp_manifest_for_output.to_csv(supp_manifest_path, index=False, encoding="utf-8-sig")


if __name__ == "__main__":
    main()
