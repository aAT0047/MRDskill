from __future__ import annotations

import argparse
import json
import ssl
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

from harvest_utils import (
    LOG_DIR,
    MERGED_DIR,
    PDF_DIR,
    SUPPLEMENTARY_DIR,
    append_rows_csv,
    attachment_extension,
    clean_text,
    content_extension,
    extract_supplementary_links,
    ensure_directories,
    is_likely_article_html,
    article_fulltext_preference_rank,
    is_supplementary_attachment,
    load_config,
    paper_file_stem,
    parse_supplementary_urls,
    supplementary_file_stem,
    unique_url_items,
    write_csv,
    write_markdown,
)


PAYWALL_MARKERS = [
    "subscription required",
    "purchase this article",
    "access through your institution",
    "institutional access",
    "paywall",
    "buy article",
]


def parse_candidate_urls(row: pd.Series) -> list[str]:
    candidates = []
    for key in ["fulltext_url", "landing_page_url", "pdf_url_candidate", "pdf_url"]:
        url = clean_text(row.get(key))
        if url:
            candidates.append(url)
    raw_json = row.get("candidate_urls_json")
    if isinstance(raw_json, str) and raw_json.strip():
        try:
            parsed = json.loads(raw_json)
            if isinstance(parsed, list):
                for item in parsed:
                    if isinstance(item, str):
                        url = clean_text(item)
                    elif isinstance(item, dict):
                        url = clean_text(item.get("url"))
                    else:
                        url = ""
                    if url:
                        candidates.append(url)
        except json.JSONDecodeError:
            pass
    unique = []
    seen = set()
    for url in sorted(candidates, key=article_fulltext_preference_rank):
        if not (url.startswith("http://") or url.startswith("https://")):
            continue
        if url not in seen:
            seen.add(url)
            unique.append(url)
    return unique


def attempt_download(url: str, timeout: int) -> dict:
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "oncology-literature-harvest/0.1",
            "Accept": "text/html, application/xml, text/xml;q=0.95, application/pdf;q=0.6, */*;q=0.1",
        },
    )
    with urllib.request.urlopen(request, timeout=timeout, context=ssl.create_default_context()) as response:
        payload = response.read()
        content_type = clean_text(response.headers.get("Content-Type")).lower()
        final_url = clean_text(response.geturl())
        status_code = int(getattr(response, "status", 200))
    return {
        "status_code": status_code,
        "content_type": content_type,
        "final_url": final_url,
        "payload": payload,
    }


def classify_payload(payload: bytes, content_type: str) -> str:
    sample = payload[:2000].decode("utf-8", errors="ignore").lower()
    if "pdf" in content_type or payload.startswith(b"%PDF"):
        return "pdf"
    if "xml" in content_type or sample.lstrip().startswith("<?xml"):
        return "xml"
    if "html" in content_type or "<html" in sample:
        return "html"
    return "other"


def looks_paywalled(text: str) -> bool:
    lowered = text.lower()
    return any(marker in lowered for marker in PAYWALL_MARKERS)


def build_download_report(status_df: pd.DataFrame) -> str:
    success_count = int((status_df["status"] == "success").sum())
    lines = [
        "# Fulltext Download Report",
        "",
        f"- Total attempted records: `{len(status_df)}`",
        f"- Successful downloads: `{success_count}`",
        f"- Metadata-only outcomes: `{int((status_df['status'] == 'metadata_only').sum())}`",
        f"- Inaccessible outcomes: `{int((status_df['status'] == 'inaccessible').sum())}`",
        f"- Broken link outcomes: `{int((status_df['status'] == 'broken_link').sum())}`",
        f"- Rate-limited outcomes: `{int((status_df['status'] == 'rate_limited').sum())}`",
        "",
        "## Notes",
        "",
        "- Only direct, legal URLs from scholarly APIs were used.",
        "- No paywall bypassing was attempted.",
        "- XML/HTML full text was saved when PDF was not directly accessible.",
    ]
    return "\n".join(lines) + "\n"


def download_supplementary_for_row(row: pd.Series, discovered_items: list[dict], timeout: int, delay_seconds: float) -> list[dict]:
    SUPPLEMENTARY_DIR.mkdir(parents=True, exist_ok=True)
    items = parse_supplementary_urls(row.get("supplementary_urls_json"))
    items.extend(discovered_items)
    items = unique_url_items(items)
    items = [item for item in items if is_supplementary_attachment(item["url"], item.get("label", ""))]
    items = sorted(items, key=supplementary_download_preference_rank)
    rows = []
    for index, item in enumerate(items, start=1):
        url = item["url"]
        label = item.get("label") or "supplementary_material"
        status = "metadata_only"
        reason = ""
        final_url = url
        output_path = ""
        content_format = ""
        try:
            response = attempt_download(url, timeout)
            payload = response["payload"]
            final_url = response["final_url"] or url
            kind = classify_payload(payload, response["content_type"])
            if kind == "html" and looks_paywalled(payload[:25000].decode("utf-8", errors="ignore")):
                status = "inaccessible"
                reason = "paywall_detected"
            elif kind == "html":
                status = "metadata_only"
                reason = "supplement_html_skipped"
            else:
                extension = attachment_extension(response["content_type"], final_url, payload)
                if extension in {".html", ".bin"}:
                    status = "metadata_only"
                    reason = "supplement_html_skipped" if extension == ".html" else "unknown_binary_skipped"
                    continue
                output = SUPPLEMENTARY_DIR / f"{supplementary_file_stem(row, index, final_url, label)}{extension}"
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
        except Exception as exc:
            reason = type(exc).__name__
        finally:
            time.sleep(delay_seconds)
        rows.append(
            {
                "dedup_group_key": clean_text(row.get("dedup_group_key")),
                "title": clean_text(row.get("title")),
                "doi": clean_text(row.get("doi")),
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=None)
    parser.add_argument("--max-records", type=int, default=None)
    args = parser.parse_args()

    ensure_directories()
    config = load_config(args.config)
    if not config.get("download", {}).get("enabled", True):
        write_csv(pd.DataFrame(), LOG_DIR / "download_status.csv")
        write_markdown(
            "# Fulltext Download Report\n\nDownloading is disabled in config.\n",
            LOG_DIR.parent / "reports" / "fulltext_download_report.md",
        )
        return

    priority_path = MERGED_DIR / "oncology_candidate_priority_table.csv"
    if not priority_path.exists():
        priority_path = MERGED_DIR / "keyword_research_candidate_table.csv"
    if not priority_path.exists():
        raise FileNotFoundError(
            "No MRD candidate-priority table was found. Run the search and "
            "deduplication stages before requesting full-text downloads."
        )
    priority_df = pd.read_csv(priority_path)
    allowed_priorities = set(config.get("download", {}).get("priority_levels", ["high", "medium"]))
    max_attempts = int(config.get("download", {}).get("max_attempts_per_record", 5))
    timeout_seconds = int(config.get("download", {}).get("timeout_seconds", 45))
    delay_seconds = float(config.get("delay_seconds", {}).get("download", 0.25))
    log_path = LOG_DIR / "download_status.csv"

    existing_df = pd.read_csv(log_path) if log_path.exists() else pd.DataFrame()
    processed_keys = set(existing_df.get("dedup_group_key", pd.Series(dtype=str)).fillna("").tolist())
    existing_files = {
        path.stem.lower(): str(path)
        for path in PDF_DIR.glob("*")
        if path.is_file()
    }

    target_df = priority_df.loc[
        priority_df["likely_priority"].isin(allowed_priorities)
    ].copy()
    target_df = target_df.loc[
        ~target_df["dedup_group_key"].fillna("").isin(processed_keys)
    ].copy()
    if args.max_records is not None and args.max_records > 0:
        target_df = target_df.head(args.max_records).copy()
    status_rows = []

    for _, row in target_df.iterrows():
        dedup_key = clean_text(row.get("dedup_group_key"))
        stem = paper_file_stem(row).lower()
        if stem in existing_files:
            row_status = {
                "dedup_group_key": dedup_key,
                "priority_rank": row.get("priority_rank"),
                "likely_priority": clean_text(row.get("likely_priority")),
                "title": clean_text(row.get("title")),
                "doi": clean_text(row.get("doi")),
                "pmid": clean_text(row.get("pmid")),
                "pmcid": clean_text(row.get("pmcid")),
                "status": "success",
                "reason": "existing_file",
                "attempt_count": 0,
                "attempted_url_count": 0,
                "selected_url": "",
                "content_format": Path(existing_files[stem]).suffix.lstrip(".").lower(),
                "local_file_path_if_downloaded": existing_files[stem],
            }
            status_rows.append(row_status)
            append_rows_csv(log_path, [row_status])
            continue

        candidate_urls = parse_candidate_urls(row)[:max_attempts]
        record_status = "metadata_only"
        reason = "no_candidate_url"
        downloaded_path = ""
        chosen_url = ""
        content_format = ""
        attempt_count = 0

        for url in candidate_urls:
            attempt_count += 1
            try:
                response = attempt_download(url, timeout_seconds)
                payload = response["payload"]
                content_type = response["content_type"]
                final_url = response["final_url"] or url
                kind = classify_payload(payload, content_type)
                if kind == "pdf":
                    extension = ".pdf"
                elif kind == "xml":
                    extension = ".xml"
                elif kind == "html":
                    html_text = payload[:25000].decode("utf-8", errors="ignore")
                    if looks_paywalled(html_text):
                        record_status = "inaccessible"
                        reason = "paywall_detected"
                        chosen_url = final_url
                        content_format = "html"
                        continue
                    if not is_likely_article_html(payload, final_url):
                        record_status = "metadata_only"
                        reason = "non_article_html"
                        chosen_url = final_url
                        content_format = "html"
                        continue
                    extension = ".html"
                else:
                    extension = content_extension(content_type, final_url)

                filename = paper_file_stem(row) + extension
                output_path = PDF_DIR / filename
                output_path.write_bytes(payload)
                discovered_items = []
                if extension in {".html", ".xml"}:
                    text = payload.decode("utf-8", errors="ignore")
                    discovered_items = extract_supplementary_links(text, [final_url])
                supp_rows = download_supplementary_for_row(row, discovered_items, timeout_seconds, delay_seconds)
                if supp_rows:
                    append_rows_csv(LOG_DIR / "supplementary_download_status.csv", supp_rows)
                record_status = "success"
                reason = "downloaded"
                downloaded_path = str(output_path)
                chosen_url = final_url
                content_format = kind
                break
            except urllib.error.HTTPError as exc:
                code = int(getattr(exc, "code", 0))
                if code == 429:
                    record_status = "rate_limited"
                    reason = f"http_{code}"
                    chosen_url = url
                    break
                if code in {401, 402, 403}:
                    record_status = "inaccessible"
                    reason = f"http_{code}"
                    chosen_url = url
                    continue
                if code in {404, 410}:
                    record_status = "broken_link"
                    reason = f"http_{code}"
                    chosen_url = url
                    continue
                record_status = "metadata_only"
                reason = f"http_{code}"
                chosen_url = url
                continue
            except Exception as exc:
                record_status = "metadata_only"
                reason = type(exc).__name__
                chosen_url = url
                continue
            finally:
                time.sleep(delay_seconds)

        row_status = {
            "dedup_group_key": dedup_key,
            "priority_rank": row.get("priority_rank"),
            "likely_priority": clean_text(row.get("likely_priority")),
            "title": clean_text(row.get("title")),
            "doi": clean_text(row.get("doi")),
            "pmid": clean_text(row.get("pmid")),
            "pmcid": clean_text(row.get("pmcid")),
            "status": record_status,
            "reason": reason,
            "attempt_count": attempt_count,
            "attempted_url_count": len(candidate_urls),
            "selected_url": chosen_url,
            "content_format": content_format,
            "local_file_path_if_downloaded": downloaded_path,
        }
        status_rows.append(row_status)
        append_rows_csv(log_path, [row_status])

    status_df = pd.read_csv(log_path) if log_path.exists() else pd.DataFrame(status_rows)
    if not status_df.empty:
        # Keep the latest result per dedup key when resumed across runs.
        status_df = status_df.drop_duplicates(subset=["dedup_group_key"], keep="last")
    write_csv(status_df, LOG_DIR / "download_status.csv")
    write_markdown(
        build_download_report(status_df),
        LOG_DIR.parent / "reports" / "fulltext_download_report.md",
    )


if __name__ == "__main__":
    main()
