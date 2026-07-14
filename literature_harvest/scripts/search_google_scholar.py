from __future__ import annotations

import argparse
import csv
import os
import re
from pathlib import Path
from typing import Any

import pandas as pd

from harvest_utils import (
    RAW_DIR,
    clean_text,
    ensure_directories,
    extract_year,
    flatten_authors,
    json_compatible,
    load_config,
    match_allowed_journal,
    normalize_doi,
    request_json,
    source_blacklist_reason,
    write_csv,
    write_json,
)


SERPAPI_URL = "https://serpapi.com/search.json"


def parse_publication_summary(summary: str) -> dict[str, str]:
    text = clean_text(summary)
    year = extract_year(text)
    journal = ""
    # Google Scholar summaries often look like:
    # "A Author, B Author - Journal Name, 2024 - publisher.com"
    parts = [clean_text(part) for part in text.split(" - ") if clean_text(part)]
    if len(parts) >= 2:
        candidate = parts[1]
        if year:
            candidate = re.sub(rf",?\s*{year}\b.*$", "", candidate).strip()
        journal = candidate
    return {"journal": journal, "year": str(year or "")}


def parse_serpapi_result(item: dict[str, Any], query_name: str, search_query: str, config: dict[str, Any]) -> dict[str, Any]:
    publication_info = item.get("publication_info") or {}
    summary = clean_text(publication_info.get("summary"))
    parsed_summary = parse_publication_summary(summary)
    resources = item.get("resources") or []
    pdf_url = ""
    resource_urls: list[str] = []
    for resource in resources:
        url = clean_text(resource.get("link"))
        if not url:
            continue
        resource_urls.append(url)
        file_format = clean_text(resource.get("file_format")).lower()
        if not pdf_url and ("pdf" in file_format or url.lower().endswith(".pdf")):
            pdf_url = url

    link = clean_text(item.get("link"))
    inline_links = item.get("inline_links") or {}
    doi = normalize_doi(item.get("doi") or inline_links.get("doi"))
    authors = flatten_authors(
        [author.get("name") for author in publication_info.get("authors", []) if isinstance(author, dict)]
    )
    journal = parsed_summary["journal"]
    if source_blacklist_reason(journal=journal, doi=doi, config=config):
        return {}
    matched_journal = match_allowed_journal(journal, config)
    year = parsed_summary["year"] or extract_year(summary)
    return {
        "query_name": query_name,
        "search_query": search_query,
        "source_database": "google_scholar",
        "source_record_id": clean_text(item.get("result_id")) or link or clean_text(item.get("title")),
        "title": clean_text(item.get("title")),
        "abstract": clean_text(item.get("snippet")),
        "authors": authors,
        "journal": journal,
        "matched_priority_journal": matched_journal,
        "year": year,
        "doi": doi,
        "pmid": "",
        "pmcid": "",
        "publisher": "",
        "publication_type": "scholar_result",
        "open_access_flag": bool(pdf_url),
        "fulltext_url": link,
        "pdf_url": pdf_url,
        "landing_page_url": link,
        "license": "",
        "licenses_json": json_compatible([]),
        "candidate_urls_json": json_compatible([url for url in [pdf_url, link, *resource_urls] if url]),
        "supplementary_urls_json": json_compatible([]),
    }


def read_manual_csv(path: Path, config: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for index, raw in enumerate(reader, start=1):
            normalized = {clean_text(key).lower().replace(" ", "_"): value for key, value in raw.items()}
            title = clean_text(normalized.get("title") or normalized.get("article_title"))
            journal = clean_text(
                normalized.get("journal")
                or normalized.get("publication")
                or normalized.get("source")
                or normalized.get("venue")
            )
            url = clean_text(normalized.get("url") or normalized.get("link"))
            pdf_url = clean_text(normalized.get("pdf_url") or normalized.get("pdf"))
            doi = normalize_doi(normalized.get("doi"))
            if source_blacklist_reason(journal=journal, doi=doi, config=config):
                continue
            matched_journal = match_allowed_journal(journal, config)
            rows.append(
                {
                    "query_name": f"manual_export:{path.stem}",
                    "search_query": str(path),
                    "source_database": "google_scholar_manual",
                    "source_record_id": doi or url or f"{path.stem}:{index}",
                    "title": title,
                    "abstract": clean_text(normalized.get("abstract") or normalized.get("snippet")),
                    "authors": clean_text(normalized.get("authors") or normalized.get("author")),
                    "journal": journal,
                    "matched_priority_journal": matched_journal,
                    "year": extract_year(normalized.get("year") or normalized.get("date")),
                    "doi": doi,
                    "pmid": "",
                    "pmcid": "",
                    "publisher": "",
                    "publication_type": clean_text(normalized.get("type")) or "scholar_manual_export",
                    "open_access_flag": bool(pdf_url),
                    "fulltext_url": url,
                    "pdf_url": pdf_url,
                    "landing_page_url": url,
                    "license": "",
                    "licenses_json": json_compatible([]),
                    "candidate_urls_json": json_compatible([value for value in [pdf_url, url] if value]),
                    "supplementary_urls_json": json_compatible([]),
                }
            )
    return rows


def parse_bibtex_entries(text: str) -> list[dict[str, str]]:
    entries: list[dict[str, str]] = []
    for match in re.finditer(r"@\w+\s*\{[^,]+,(.*?)(?=\n@\w+\s*\{|\Z)", text, flags=re.DOTALL):
        body = match.group(1)
        fields: dict[str, str] = {}
        for field_match in re.finditer(
            r"(\w+)\s*=\s*[\{\"](.+?)[\}\"]\s*,?\s*(?=\n\s*\w+\s*=|\Z)",
            body,
            flags=re.DOTALL,
        ):
            key = field_match.group(1).lower()
            value = re.sub(r"\s+", " ", field_match.group(2)).strip()
            fields[key] = value
        if fields:
            entries.append(fields)
    return entries


def read_manual_bibtex(path: Path, config: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    entries = parse_bibtex_entries(path.read_text(encoding="utf-8", errors="ignore"))
    for index, entry in enumerate(entries, start=1):
        title = clean_text(entry.get("title"))
        journal = clean_text(entry.get("journal") or entry.get("booktitle"))
        url = clean_text(entry.get("url"))
        doi = normalize_doi(entry.get("doi"))
        if source_blacklist_reason(journal=journal, doi=doi, config=config):
            continue
        matched_journal = match_allowed_journal(journal, config)
        rows.append(
            {
                "query_name": f"manual_bibtex:{path.stem}",
                "search_query": str(path),
                "source_database": "google_scholar_bibtex",
                "source_record_id": doi or url or f"{path.stem}:{index}",
                "title": title,
                "abstract": "",
                "authors": clean_text(entry.get("author")),
                "journal": journal,
                "matched_priority_journal": matched_journal,
                "year": extract_year(entry.get("year")),
                "doi": doi,
                "pmid": "",
                "pmcid": "",
                "publisher": "",
                "publication_type": "scholar_bibtex_export",
                "open_access_flag": False,
                "fulltext_url": url,
                "pdf_url": "",
                "landing_page_url": url,
                "license": "",
                "licenses_json": json_compatible([]),
                "candidate_urls_json": json_compatible([url] if url else []),
                "supplementary_urls_json": json_compatible([]),
            }
        )
    return rows


def read_manual_exports(config: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    scholar_config = config.get("google_scholar", {})
    for raw_path in scholar_config.get("manual_exports", []):
        path = Path(raw_path)
        if not path.is_absolute():
            path = Path.cwd() / path
        if not path.exists():
            continue
        if path.suffix.lower() == ".csv":
            rows.extend(read_manual_csv(path, config))
        elif path.suffix.lower() in {".bib", ".bibtex"}:
            rows.extend(read_manual_bibtex(path, config))
    return rows


def search_serpapi(config: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    scholar_config = config.get("google_scholar", {})
    api_key_env = clean_text(scholar_config.get("api_key_env")) or "SERPAPI_KEY"
    api_key = clean_text(scholar_config.get("api_key")) or clean_text(os.environ.get(api_key_env))
    if not api_key:
        return [], [{"db": "google_scholar_serpapi", "status": "skipped_missing_api_key", "api_key_env": api_key_env}]

    rows: list[dict[str, Any]] = []
    logs: list[dict[str, Any]] = []
    max_results = int(config.get("max_results_per_query", {}).get("google_scholar", 20))
    page_size = min(int(config.get("page_size", {}).get("google_scholar", 20)), 20)
    delay = float(config.get("delay_seconds", {}).get("google_scholar", 1.0))
    year_low = clean_text(scholar_config.get("as_ylo"))
    year_high = clean_text(scholar_config.get("as_yhi"))

    for query in config.get("queries", []):
        retrieved = 0
        start = 0
        first_payload: dict[str, Any] | None = None
        while retrieved < max_results:
            payload = request_json(
                SERPAPI_URL,
                params={
                    "engine": "google_scholar",
                    "q": query["query"],
                    "api_key": api_key,
                    "num": min(page_size, max_results - retrieved),
                    "start": start,
                    **({"as_ylo": year_low} if year_low else {}),
                    **({"as_yhi": year_high} if year_high else {}),
                },
                delay_seconds=delay,
            )
            if first_payload is None:
                first_payload = payload
                write_json(payload, RAW_DIR / f"google_scholar_serpapi_{query['name']}.json")
            organic = payload.get("organic_results", []) or []
            if not organic:
                break
            for item in organic:
                row = parse_serpapi_result(item, query["name"], query["query"], config)
                if row:
                    rows.append(row)
                retrieved += 1
                if retrieved >= max_results:
                    break
            if len(organic) < page_size:
                break
            start += page_size
        logs.append(
            {
                "query_name": query["name"],
                "search_query": query["query"],
                "db": "google_scholar_serpapi",
                "retrieved_row_count": retrieved,
                "status": "ok",
            }
        )
    return rows, logs


def search_google_scholar(config_path: str | None = None) -> pd.DataFrame:
    ensure_directories()
    config = load_config(config_path)
    if not config.get("sources", {}).get("google_scholar", False):
        return pd.DataFrame()

    rows: list[dict[str, Any]] = []
    logs: list[dict[str, Any]] = []
    rows.extend(read_manual_exports(config))
    serp_rows, serp_logs = search_serpapi(config)
    rows.extend(serp_rows)
    logs.extend(serp_logs)
    logs.append({"db": "google_scholar_manual", "retrieved_row_count": len(rows) - len(serp_rows)})

    df = pd.DataFrame(rows)
    write_csv(df, RAW_DIR / "google_scholar_records.csv")
    write_csv(pd.DataFrame(logs), RAW_DIR / "google_scholar_query_log.csv")
    return df


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=None)
    args = parser.parse_args()
    search_google_scholar(args.config)


if __name__ == "__main__":
    main()
