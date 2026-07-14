from __future__ import annotations

import argparse
import shutil
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
from pathlib import Path
from typing import Any

import pandas as pd


CURRENT_CONFIG: dict[str, Any] = {}


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def detect_terms(text: str, terms: list[str]) -> list[str]:
    lowered = text.lower()
    return [term for term in terms if term.lower() in lowered]


def article_type_guess(title: str, publication_type: str) -> str:
    text = f"{title} {publication_type}".lower()
    if (
        re.search(r"\babstract\s+(lb|ct|p|po|o|e|eps|suppl)?\s*-?\s*\d+", text)
        or re.search(r"^\s*(lb|ct|p|po|o|e|eps)\s*-?\s*\d+", text)
        or re.search(r"^\s*\d{2,4}o\s+", text)
        or "conference abstract" in text
    ):
        return "exclude_non_article"
    if any(token in text for token in ["review", "systematic review", "meta-analysis"]):
        return "review"
    if any(token in text for token in ["editorial", "commentary", "news", "conference abstract", "patent"]):
        return "exclude_non_article"
    if any(token in text for token in ["journal article", "research article", "article"]):
        return "research_article"
    return "research_article_like"


def route_label(url: str) -> str:
    lowered = (url or "").lower()
    if "unpaywall" in lowered:
        return "unpaywall"
    if "ftp.ncbi.nlm.nih.gov/pub/pmc" in lowered:
        return "pmc_oa"
    if "europepmc" in lowered:
        return "europe_pmc"
    if "pmc.ncbi.nlm.nih.gov" in lowered:
        return "pmc"
    if "doi.org" in lowered:
        return "doi_resolve"
    return "publisher_direct"


def candidate_urls(row: pd.Series, utils: Any) -> list[str]:
    urls: list[str] = []
    for key in ["fulltext_url", "landing_page_url", "pdf_url_candidate", "pdf_url"]:
        value = utils.clean_text(row.get(key))
        if value.startswith("http://") or value.startswith("https://"):
            urls.append(value)
    raw_json = utils.clean_text(row.get("candidate_urls_json"))
    if raw_json:
        try:
            parsed = json.loads(raw_json)
            if isinstance(parsed, list):
                for item in parsed:
                    if isinstance(item, str):
                        value = utils.clean_text(item)
                    elif isinstance(item, dict):
                        value = utils.clean_text(item.get("url"))
                    else:
                        value = ""
                    if value.startswith("http://") or value.startswith("https://"):
                        urls.append(value)
        except json.JSONDecodeError:
            pass
    doi = utils.normalize_doi(row.get("doi"))
    if doi:
        urls.append(f"https://doi.org/{doi}")
    unique: list[str] = []
    seen: set[str] = set()
    for url in sorted(urls, key=utils.article_fulltext_preference_rank):
        if url not in seen:
            seen.add(url)
            unique.append(url)
    return unique


def supplementary_urls(row: pd.Series, utils: Any) -> list[dict[str, str]]:
    items = utils.parse_supplementary_urls(row.get("supplementary_urls_json"))
    raw_json = utils.clean_text(row.get("candidate_urls_json"))
    if raw_json:
        try:
            parsed = json.loads(raw_json)
        except json.JSONDecodeError:
            parsed = []
        for item in parsed if isinstance(parsed, list) else []:
            if isinstance(item, str):
                url = utils.clean_text(item)
                label = ""
            elif isinstance(item, dict):
                url = utils.clean_text(item.get("url"))
                label = utils.clean_text(item.get("label") or item.get("style") or item.get("type"))
            else:
                continue
            if url and utils.is_supplementary_attachment(url, label):
                items.append({"url": url, "label": label or "supplementary_material", "source": "candidate_urls"})
    return utils.unique_url_items(items)


def discover_supplementary_from_saved_file(
    path: Path,
    row: pd.Series,
    utils: Any,
) -> list[dict[str, str]]:
    if path.suffix.lower() not in {".html", ".htm", ".xml"}:
        return []
    text = path.read_text(encoding="utf-8", errors="ignore")
    bases = [path.name]
    bases.extend(candidate_urls(row, utils))
    return utils.extract_supplementary_links(text, bases)


def discover_supplementary_from_candidate_pages(
    row: pd.Series,
    utils: Any,
    downloader: Any,
    timeout_seconds: int,
    delay_seconds: float,
    max_pages: int,
    exclude_urls: set[str] | None = None,
) -> list[dict[str, str]]:
    if max_pages <= 0:
        return []
    exclude = {url.lower() for url in (exclude_urls or set()) if url}
    discovered: list[dict[str, str]] = []
    checked = 0
    for url in candidate_urls(row, utils):
        if url.lower() in exclude:
            continue
        lowered = url.lower()
        if lowered.endswith(".pdf") or "/pdf" in lowered or "pdf=render" in lowered:
            continue
        try:
            response = downloader.attempt_download(url, timeout_seconds)
            payload = response["payload"]
            content_type = response["content_type"]
            final_url = response["final_url"] or url
            kind = downloader.classify_payload(payload, content_type)
            if kind not in {"html", "xml"}:
                continue
            if kind == "html":
                html_text = payload[:25000].decode("utf-8", errors="ignore")
                if downloader.looks_paywalled(html_text):
                    continue
                if not utils.is_likely_article_html(payload, final_url):
                    continue
            text = payload.decode("utf-8", errors="ignore")
            discovered.extend(utils.extract_supplementary_links(text, [final_url, url]))
            checked += 1
            if checked >= max_pages:
                break
        except Exception:  # noqa: BLE001
            continue
        finally:
            time.sleep(delay_seconds)
    return utils.unique_url_items(discovered)


def download_supplementary_files(
    row: pd.Series,
    run_root: Path,
    utils: Any,
    downloader: Any,
    timeout_seconds: int,
    delay_seconds: float,
    max_items: int | None = None,
    discovered_items: list[dict[str, str]] | None = None,
) -> list[dict[str, Any]]:
    supp_dir = run_root / "supplementary_materials"
    supp_dir.mkdir(parents=True, exist_ok=True)
    items = supplementary_urls(row, utils)
    if discovered_items:
        items.extend(discovered_items)
    items = utils.unique_url_items(items)
    items = sorted(items, key=utils.supplementary_download_preference_rank)
    if max_items and max_items > 0:
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
            response = downloader.attempt_download(url, timeout_seconds)
            payload = response["payload"]
            content_type = response["content_type"]
            final_url = response["final_url"] or url
            kind = downloader.classify_payload(payload, content_type)
            if kind == "html":
                html_text = payload[:25000].decode("utf-8", errors="ignore")
                if downloader.looks_paywalled(html_text):
                    status = "inaccessible"
                    reason = "paywall_detected"
                    rows.append(
                        {
                            "record_id": row["record_id"],
                            "supplement_index": index,
                            "title": row.get("title"),
                            "doi": row.get("doi"),
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
                if not utils.is_likely_article_html(payload, final_url):
                    status = "metadata_only"
                    reason = "non_attachment_html"
                    rows.append(
                        {
                            "record_id": row["record_id"],
                            "supplement_index": index,
                            "title": row.get("title"),
                            "doi": row.get("doi"),
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
            if kind == "html":
                status = "metadata_only"
                reason = "supplement_html_skipped"
                rows.append(
                    {
                        "record_id": row["record_id"],
                        "supplement_index": index,
                        "title": row.get("title"),
                        "doi": row.get("doi"),
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
            extension = utils.attachment_extension(content_type, final_url, payload)
            if extension in {".html", ".bin"}:
                status = "metadata_only"
                reason = "supplement_html_skipped" if extension == ".html" else "unknown_binary_skipped"
                rows.append(
                    {
                        "record_id": row["record_id"],
                        "supplement_index": index,
                        "title": row.get("title"),
                        "doi": row.get("doi"),
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
            time.sleep(delay_seconds)
        rows.append(
            {
                "record_id": row["record_id"],
                "supplement_index": index,
                "title": row.get("title"),
                "doi": row.get("doi"),
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
    run_root: Path,
    utils: Any,
    downloader: Any,
    timeout_seconds: int,
    delay_seconds: float,
) -> list[dict[str, Any]]:
    package_dir = run_root / "pmc_oa_packages"
    extracted_dir = run_root / "pmc_oa_extracted" / utils.safe_filename(utils.clean_text(row.get("record_id")), 40)
    paper_dir = run_root / "downloaded_pdfs"
    supp_dir = run_root / "supplementary_materials"
    package_dir.mkdir(parents=True, exist_ok=True)
    extracted_dir.mkdir(parents=True, exist_ok=True)
    paper_dir.mkdir(parents=True, exist_ok=True)
    supp_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    for index, item in enumerate(utils.pmc_oa_package_urls(row), start=1):
        url = item["url"]
        label = item.get("label") or "pmc_oa_package"
        status = "metadata_only"
        reason = ""
        package_path = ""
        extracted_count = 0
        paper_files: list[str] = []
        supplementary_files: list[str] = []
        final_url = url
        try:
            response = downloader.attempt_download(url, timeout_seconds)
            payload = response["payload"]
            final_url = response["final_url"] or url
            extension = utils.attachment_extension(response["content_type"], final_url, payload)
            if extension not in {".tar.gz", ".tgz", ".gz"}:
                status = "metadata_only"
                reason = f"unexpected_package_extension:{extension}"
            else:
                package_path_obj = package_dir / f"{utils.supplementary_file_stem(row, index, final_url, label)}{extension}"
                package_path_obj.write_bytes(payload)
                package_path = str(package_path_obj)
                extracted = utils.safe_extract_tar_members(package_path_obj, extracted_dir)
                extracted_count = len(extracted)
                for file_index, path in enumerate(extracted, start=1):
                    suffix = path.suffix.lower()
                    if suffix == ".pdf" and path.read_bytes()[:5] == b"%PDF-":
                        destination = paper_dir / f"{utils.paper_file_stem(row)}__pmc_oa_package_{file_index:02d}.pdf"
                        shutil.copy2(path, destination)
                        paper_files.append(str(destination))
                    elif suffix in {".nxml", ".xml"}:
                        destination = paper_dir / f"{utils.paper_file_stem(row)}__pmc_oa_package_{file_index:02d}.xml"
                        shutil.copy2(path, destination)
                        paper_files.append(str(destination))
                    elif suffix in utils.SUPPLEMENTARY_ATTACHMENT_EXTENSIONS and suffix not in {".tar", ".gz", ".tgz", ".tar.gz"}:
                        destination = supp_dir / f"{utils.supplementary_file_stem(row, file_index, final_url, path.name)}{suffix}"
                        shutil.copy2(path, destination)
                        supplementary_files.append(str(destination))
                status = "success"
        except Exception as exc:  # noqa: BLE001
            status = "metadata_only"
            reason = type(exc).__name__
        finally:
            time.sleep(delay_seconds)
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
                "paper_files": utils.json_compatible(paper_files),
                "supplementary_files": utils.json_compatible(supplementary_files),
            }
        )
    return rows


def build_candidate_rows(raw_df: pd.DataFrame, config: dict[str, Any], utils: Any) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    include_terms = config.get("include_terms", [])
    secondary_terms = config.get("secondary_terms", [])
    exclude_terms = config.get("exclude_terms", [])
    exclude_reviews = bool(config.get("exclude_reviews_by_default", False))

    for idx, row in raw_df.reset_index(drop=True).iterrows():
        title = utils.clean_text(row.get("title"))
        abstract = utils.clean_text(row.get("abstract"))
        publication_type = utils.clean_text(row.get("publication_type"))
        journal = utils.clean_text(row.get("journal"))
        doi = utils.normalize_doi(row.get("doi"))
        publisher = utils.clean_text(row.get("publisher") or row.get("source_publisher") or row.get("publisher_name"))
        blacklist_reason = utils.source_blacklist_reason(journal=journal, doi=doi, publisher=publisher, config=config)
        year_filter_reason = utils.publication_year_exclusion_reason(row.get("year"), config)
        matched_journal = "" if blacklist_reason else (
            utils.clean_text(row.get("matched_priority_journal")) or utils.match_allowed_journal(journal, config)
        )
        journal_tier = utils.priority_journal_tier(matched_journal or journal, config)
        query_name = utils.clean_text(row.get("query_name"))
        search_query = utils.clean_text(row.get("search_query"))
        evidence_text = " ".join(part for part in [title, abstract, publication_type] if part)
        query_text = " ".join(part for part in [query_name, search_query] if part)
        type_guess = article_type_guess(title, publication_type)
        include_hits = detect_terms(evidence_text, include_terms)
        secondary_hits = detect_terms(evidence_text, secondary_terms)
        exclude_hits = detect_terms(evidence_text, exclude_terms)
        query_hits = detect_terms(query_text, include_terms + secondary_terms)

        score = 0
        if matched_journal:
            score += 2
        if include_hits:
            score += 2 + min(len(include_hits), 3)
        elif query_hits:
            score += 1
        if secondary_hits:
            score += min(len(secondary_hits), 2)
        if bool(row.get("open_access_flag")):
            score += 1
        if type_guess == "review":
            score -= 1
        if type_guess == "exclude_non_article":
            score -= 3
        if exclude_hits and not include_hits:
            score -= 2

        exclusion_reason = ""
        if blacklist_reason:
            exclusion_reason = f"blacklisted_source:{blacklist_reason}"
        elif year_filter_reason:
            exclusion_reason = year_filter_reason
        elif not title:
            exclusion_reason = "missing_title"
        elif config.get("restrict_to_allowed_journals", True) and not matched_journal:
            exclusion_reason = "journal_not_in_priority_list"
        elif type_guess == "exclude_non_article":
            exclusion_reason = "non_article_type"
        elif exclude_reviews and type_guess == "review":
            exclusion_reason = "review_deprioritized"
        elif exclude_hits and not include_hits:
            exclusion_reason = "excluded_by_title_abstract_terms"
        elif not include_hits and not secondary_hits and not query_hits:
            exclusion_reason = "query_only_without_keyword_signal"
        elif score < 1:
            exclusion_reason = "low_relevance"

        priority = "high" if score >= 4 else ("medium" if score >= 2 else "low")
        record = {
                "record_id": f"KWRAW-{idx + 1:06d}",
                "title": title,
                "authors": utils.clean_text(row.get("authors")),
                "year": row.get("year"),
                "journal": journal,
                "matched_priority_journal": matched_journal,
                "priority_journal_tier": journal_tier,
                "doi": doi,
                "pmid": utils.normalize_pmid(row.get("pmid")),
                "pmcid": utils.normalize_pmcid(row.get("pmcid"))
                or utils.extract_pmcid_from_text(row.get("fulltext_url"))
                or utils.extract_pmcid_from_text(row.get("pdf_url"))
                or utils.extract_pmcid_from_text(row.get("candidate_urls_json")),
                "publisher": publisher,
                "source_blacklist_reason": blacklist_reason,
                "publication_year_filter_reason": year_filter_reason,
                "source_database": utils.clean_text(row.get("source_database")),
                "abstract_if_available": abstract,
                "article_type_guess": type_guess,
                "keyword_include_hits": "; ".join(include_hits),
                "keyword_secondary_hits": "; ".join(secondary_hits),
                "keyword_exclude_hits": "; ".join(exclude_hits),
                "keyword_relevance_score": score,
                "research_article_flag": "yes" if type_guess in {"research_article", "research_article_like"} else "no",
                "likely_topic_tags": "; ".join(sorted(set(include_hits + secondary_hits))),
                "open_access_status_if_detectable": "open_access" if bool(row.get("open_access_flag")) else "not_detected",
                "pdf_url_candidate": utils.clean_text(row.get("pdf_url")) or utils.clean_text(row.get("fulltext_url")) or utils.clean_text(row.get("landing_page_url")),
                "landing_page_url": utils.clean_text(row.get("landing_page_url")),
                "download_status": "excluded" if exclusion_reason else "pending",
                "exclusion_reason_if_any": exclusion_reason,
                "query_name": query_name,
                "search_query": search_query,
                "source_record_id": utils.clean_text(row.get("source_record_id")),
                "candidate_urls_json": row.get("candidate_urls_json", "[]"),
                "supplementary_urls_json": row.get("supplementary_urls_json", "[]"),
                "supplementary_url_count": len(utils.parse_supplementary_urls(row.get("supplementary_urls_json"))),
                "pdf_url": utils.clean_text(row.get("pdf_url")),
                "fulltext_url": utils.clean_text(row.get("fulltext_url")),
                "publication_type": publication_type,
                "pmc_oa_package_urls_json": row.get("pmc_oa_package_urls_json", "[]"),
                "external_fulltext_sources": "",
            }
        rows.append(record)
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["dedup_download_key"] = out.apply(
        lambda item: utils.normalize_doi(item.get("doi"))
        or utils.normalize_title(item.get("title"))
        or utils.clean_text(item.get("record_id")),
        axis=1,
    )
    return out


def rank_candidate_df(candidate_df: pd.DataFrame) -> pd.DataFrame:
    if candidate_df.empty:
        return candidate_df
    data = candidate_df.copy()
    data["_candidate_url_rank"] = data.get("candidate_urls_json", pd.Series("", index=data.index)).fillna("").astype(str).str.len().gt(2)
    data["_package_rank"] = data.get("pmc_oa_package_urls_json", pd.Series("", index=data.index)).fillna("").astype(str).str.len().gt(2)
    data["_downloadable_rank"] = (
        data["pdf_url"].fillna("").ne("")
        | data["fulltext_url"].fillna("").ne("")
        | data["landing_page_url"].fillna("").ne("")
        | data["_candidate_url_rank"]
        | data["_package_rank"]
    )
    data["_year_sort"] = pd.to_numeric(data["year"], errors="coerce")
    data = data.sort_values(
        [
            "exclusion_reason_if_any",
            "priority_journal_tier",
            "keyword_relevance_score",
            "_downloadable_rank",
            "_year_sort",
        ],
        ascending=[True, True, False, False, False],
        na_position="last",
    )
    data = data.drop_duplicates(subset=["dedup_download_key"], keep="first")
    return data.drop(columns=["_candidate_url_rank", "_package_rank", "_downloadable_rank", "_year_sort"]).reset_index(drop=True)


def run_downloads(candidate_df: pd.DataFrame, run_root: Path, utils: Any, downloader: Any, config: dict[str, Any]) -> pd.DataFrame:
    pdf_dir = run_root / "downloaded_pdfs"
    log_dir = run_root / "download_logs"
    supp_log_path = log_dir / "keyword_research_supplementary_download_log.csv"
    package_log_path = log_dir / "keyword_research_pmc_oa_package_log.csv"
    pdf_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / "keyword_research_download_log.csv"

    timeout_seconds = int(config.get("download", {}).get("timeout_seconds", 30))
    max_attempts = int(config.get("download", {}).get("max_attempts_per_record", 3))
    max_records = int(config.get("download", {}).get("max_records", 0) or 0)
    max_supplementary = int(config.get("download", {}).get("max_supplementary_per_record", 0) or 0)
    max_supplementary_discovery_pages = int(config.get("download", {}).get("max_supplementary_discovery_pages", 3) or 0)
    delay_seconds = float(config.get("delay_seconds", {}).get("download", 0.2))

    existing_log = pd.read_csv(log_path) if log_path.exists() else pd.DataFrame()
    processed = set(existing_log.get("record_id", pd.Series(dtype=str)).astype(str))
    target_df = candidate_df.loc[
        candidate_df["exclusion_reason_if_any"].fillna("").eq("")
        & ~candidate_df["record_id"].astype(str).isin(processed)
    ].copy()
    if "dedup_download_key" in target_df.columns:
        target_df = target_df.drop_duplicates(subset=["dedup_download_key"], keep="first")
    target_df = utils.apply_year_quota(target_df, config)
    if max_records > 0 and not config.get("download", {}).get("year_quota", {}).get("enabled", False):
        target_df = target_df.head(max_records).copy()
    elif max_records > 0:
        target_df = target_df.head(max_records).copy()
    if not target_df.empty and config.get("external_fulltext", {}).get("enabled", True):
        enriched_rows = []
        for _, row in target_df.iterrows():
            enriched = utils.append_external_links_to_row(row, config)
            enriched["supplementary_url_count"] = len(utils.parse_supplementary_urls(enriched.get("supplementary_urls_json")))
            enriched_rows.append(enriched)
        target_df = pd.DataFrame(enriched_rows)
        for _, enriched in target_df.iterrows():
            record_id = utils.clean_text(enriched.get("record_id"))
            if not record_id:
                continue
            mask = candidate_df["record_id"].astype(str).eq(record_id)
            for column, value in enriched.items():
                if column in candidate_df.columns:
                    candidate_df.loc[mask, column] = value
                else:
                    candidate_df[column] = ""
                    candidate_df.loc[mask, column] = value

    rows: list[dict[str, Any]] = []
    for idx, (_, row) in enumerate(target_df.iterrows(), start=1):
        status = "metadata_only"
        reason = "no_candidate_url"
        final_path = ""
        final_url = ""
        access_route = "failed"
        content_format = ""
        attempts = 0
        discovered_supplementary: list[dict[str, str]] = []

        for url in candidate_urls(row, utils)[:max_attempts]:
            attempts += 1
            try:
                response = downloader.attempt_download(url, timeout_seconds)
                payload = response["payload"]
                content_type = response["content_type"]
                final_url = response["final_url"] or url
                kind = downloader.classify_payload(payload, content_type)
                if kind == "html":
                    html_text = payload[:25000].decode("utf-8", errors="ignore")
                    if downloader.looks_paywalled(html_text):
                        status = "inaccessible"
                        reason = "paywall_detected"
                        access_route = route_label(final_url)
                        continue
                    if not utils.is_likely_article_html(payload, final_url):
                        status = "metadata_only"
                        reason = "non_article_html"
                        access_route = route_label(final_url)
                        continue
                extension = ".pdf" if kind == "pdf" else utils.content_extension(content_type, final_url)
                if kind == "html":
                    extension = ".html"
                elif kind == "xml":
                    extension = ".xml"
                output = pdf_dir / f"{utils.paper_file_stem(row)}{extension}"
                output.write_bytes(payload)
                discovered_supplementary = discover_supplementary_from_saved_file(output, row, utils)
                status = "success"
                reason = ""
                final_path = str(output)
                access_route = route_label(final_url)
                content_format = kind
                break
            except urllib.error.HTTPError as exc:
                code = int(getattr(exc, "code", 0))
                status = "inaccessible" if code in {401, 402, 403} else ("broken_link" if code in {404, 410} else "metadata_only")
                if code == 429:
                    status = "rate_limited"
                reason = f"http_{code}"
                final_url = url
                access_route = route_label(url)
                if code == 429:
                    break
                continue
            except Exception as exc:  # noqa: BLE001
                status = "metadata_only"
                reason = type(exc).__name__
                final_url = url
                access_route = route_label(url)
                continue
            finally:
                time.sleep(delay_seconds)

        discovered_supplementary.extend(
            discover_supplementary_from_candidate_pages(
                row,
                utils,
                downloader,
                timeout_seconds,
                delay_seconds,
                max_supplementary_discovery_pages,
                {final_url},
            )
        )
        discovered_supplementary = utils.unique_url_items(discovered_supplementary)

        rows.append(
            {
                "record_id": row["record_id"],
                "doi": row.get("doi"),
                "title": row.get("title"),
                "final_pdf_path": final_path,
                "final_pdf_url": final_url,
                "download_status": status,
                "failure_reason": reason,
                "access_route_used": access_route,
                "content_format": content_format,
                "attempt_count": attempts,
                "supplementary_url_count": len(supplementary_urls(row, utils)) + len(discovered_supplementary),
            }
        )
        supp_rows = download_supplementary_files(
            row,
            run_root,
            utils,
            downloader,
            timeout_seconds,
            delay_seconds,
            max_supplementary,
            discovered_supplementary,
        )
        package_rows = download_pmc_oa_packages(
            row,
            run_root,
            utils,
            downloader,
            timeout_seconds,
            delay_seconds,
        )
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
        if idx % 50 == 0:
            chunk = pd.DataFrame(rows)
            if log_path.exists():
                combined = pd.concat([pd.read_csv(log_path), chunk], ignore_index=True)
            else:
                combined = chunk
            combined = combined.drop_duplicates(subset=["record_id"], keep="last")
            combined.to_csv(log_path, index=False, encoding="utf-8-sig")
            rows = []

    if log_path.exists():
        combined = pd.concat([pd.read_csv(log_path), pd.DataFrame(rows)], ignore_index=True)
    else:
        combined = pd.DataFrame(rows)
    if not combined.empty:
        combined = combined.drop_duplicates(subset=["record_id"], keep="last")
    combined.to_csv(log_path, index=False, encoding="utf-8-sig")
    return combined


def write_summary(run_root: Path, candidate_df: pd.DataFrame, high_df: pd.DataFrame, medium_df: pd.DataFrame, log_df: pd.DataFrame) -> None:
    success = int(log_df["download_status"].eq("success").sum()) if not log_df.empty else 0
    pdf_count = int(log_df.loc[log_df["download_status"].eq("success"), "content_format"].fillna("").eq("pdf").sum()) if not log_df.empty else 0
    non_pdf = success - pdf_count
    pending = int(candidate_df["download_status"].eq("pending").sum())
    supp_log = run_root / "download_logs" / "keyword_research_supplementary_download_log.csv"
    supp_df = pd.read_csv(supp_log) if supp_log.exists() else pd.DataFrame()
    supp_success = int(supp_df.get("download_status", pd.Series(dtype=str)).eq("success").sum()) if not supp_df.empty else 0
    lines = [
        "# Keyword Research Harvest Summary",
        "",
        f"- Candidate rows: `{len(candidate_df)}`",
        f"- High-priority rows: `{len(high_df)}`",
        f"- Medium-priority rows: `{len(medium_df)}`",
        f"- Priority journal rows: `{int(candidate_df['matched_priority_journal'].fillna('').ne('').sum()) if 'matched_priority_journal' in candidate_df.columns else 0}`",
        f"- Successful downloads/full texts: `{success}`",
        f"- True PDFs: `{pdf_count}`",
        f"- HTML/XML full texts: `{non_pdf}`",
        f"- Successful supplementary downloads: `{supp_success}`",
        f"- Remaining pending: `{pending}`",
        f"- Download log: `{run_root / 'download_logs' / 'keyword_research_download_log.csv'}`",
        f"- Supplementary log: `{supp_log}`",
    ]
    (run_root / "keyword_research_harvest_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", required=True, help="Parent folder where the new run folder should be created")
    parser.add_argument("--config", required=True, help="JSON config path")
    parser.add_argument("--run-name", required=True, help="New run folder name under literature_harvest")
    parser.add_argument("--skip-search", action="store_true", help="Skip the search stage and build candidates from existing raw API output.")
    parser.add_argument("--search-only", action="store_true", help="Build candidate tables and a search summary without downloading full text or supplements.")
    parser.add_argument("--download-only", action="store_true", help="Use an existing candidate table and download MRD literature only.")
    parser.add_argument("--candidate-table", default=None, help="Path to an existing candidate table CSV for download-only mode.")
    args = parser.parse_args()

    if args.search_only and args.download_only:
        parser.error("--search-only cannot be combined with --download-only")

    output_root = Path(args.output_root).resolve()
    run_root = output_root / args.run_name
    run_root.mkdir(parents=True, exist_ok=True)
    os.environ["MRDSKILL_HARVEST_ROOT"] = str(run_root)
    os.environ["LITERATURE_HARVEST_ROOT"] = str(run_root)

    skill_root = Path(__file__).resolve().parents[1]
    scripts_dir = skill_root / "literature_harvest" / "scripts"
    sys.path.insert(0, str(scripts_dir))

    from download_fulltexts import attempt_download, classify_payload, looks_paywalled  # noqa: WPS433
    from harvest_utils import ensure_directories, load_config, write_csv  # noqa: WPS433
    from merge_and_deduplicate import load_sources, normalize_columns  # noqa: WPS433
    from search_crossref import search_crossref  # noqa: WPS433
    from search_europepmc import search_europepmc  # noqa: WPS433
    from search_google_scholar import search_google_scholar  # noqa: WPS433
    from search_openalex import search_openalex  # noqa: WPS433
    from search_pubmed import search_pubmed_and_pmc  # noqa: WPS433
    import harvest_utils as utils  # noqa: WPS433

    config = load_config(args.config)
    global CURRENT_CONFIG
    CURRENT_CONFIG = config
    ensure_directories()
    (run_root / "config.json").write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")

    downloader = type(
        "DownloaderNamespace",
        (),
        {
            "attempt_download": staticmethod(attempt_download),
            "classify_payload": staticmethod(classify_payload),
            "looks_paywalled": staticmethod(looks_paywalled),
        },
    )

    if args.download_only:
        candidate_path = Path(args.candidate_table) if args.candidate_table else run_root / "keyword_research_candidate_table.csv"
        if not candidate_path.exists():
            raise FileNotFoundError(f"Candidate table not found at {candidate_path}")
        candidate_df = pd.read_csv(candidate_path)
        log_df = run_downloads(candidate_df, run_root, utils, downloader, config)

        status = log_df[["record_id", "download_status"]].drop_duplicates("record_id", keep="last")
        candidate_df = candidate_df.drop(columns=["download_status"], errors="ignore").merge(status, on="record_id", how="left")
        candidate_df["download_status"] = candidate_df["download_status"].fillna(candidate_df["exclusion_reason_if_any"].map(lambda x: "excluded" if str(x).strip() else "pending"))
        write_csv(candidate_df, run_root / "keyword_research_candidate_table.csv")

        high_df = candidate_df.loc[candidate_df["exclusion_reason_if_any"].fillna("").eq("") & candidate_df["keyword_relevance_score"].ge(4)].copy()
        medium_df = candidate_df.loc[candidate_df["exclusion_reason_if_any"].fillna("").eq("") & candidate_df["keyword_relevance_score"].between(2, 3)].copy()
        write_csv(high_df, run_root / "keyword_research_high_priority.csv")
        write_csv(medium_df, run_root / "keyword_research_medium_priority.csv")
        write_summary(run_root, candidate_df, high_df, medium_df, log_df)
        return

    if not args.skip_search:
        search_pubmed_and_pmc(args.config)
        search_europepmc(args.config)
        search_crossref(args.config)
        search_openalex(args.config)
        search_google_scholar(args.config)

    raw_df, _ = load_sources()
    raw_df = normalize_columns(raw_df, config)
    candidate_df = build_candidate_rows(raw_df, config, utils)
    candidate_df = rank_candidate_df(candidate_df)
    write_csv(candidate_df, run_root / "keyword_research_candidate_table.csv")

    high_df = candidate_df.loc[candidate_df["exclusion_reason_if_any"].eq("") & candidate_df["keyword_relevance_score"].ge(4)].copy()
    medium_df = candidate_df.loc[candidate_df["exclusion_reason_if_any"].eq("") & candidate_df["keyword_relevance_score"].between(2, 3)].copy()
    write_csv(high_df, run_root / "keyword_research_high_priority.csv")
    write_csv(medium_df, run_root / "keyword_research_medium_priority.csv")

    if args.search_only or not config.get("download", {}).get("enabled", True):
        write_summary(run_root, candidate_df, high_df, medium_df, pd.DataFrame())
        return

    log_df = run_downloads(candidate_df, run_root, utils, downloader, config)

    status = log_df[["record_id", "download_status"]].drop_duplicates("record_id", keep="last")
    candidate_df = candidate_df.drop(columns=["download_status"], errors="ignore").merge(status, on="record_id", how="left")
    candidate_df["download_status"] = candidate_df["download_status"].fillna(candidate_df["exclusion_reason_if_any"].map(lambda x: "excluded" if x else "pending"))
    write_csv(candidate_df, run_root / "keyword_research_candidate_table.csv")
    write_summary(run_root, candidate_df, high_df, medium_df, log_df)


if __name__ == "__main__":
    main()
