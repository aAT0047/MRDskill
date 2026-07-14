from __future__ import annotations

import argparse

import pandas as pd

from harvest_utils import (
    CROSSREF_URL,
    RAW_DIR,
    USER_AGENT,
    clean_text,
    crossref_publication_date_filter,
    ensure_directories,
    extract_year,
    is_supplementary_url,
    json_compatible,
    load_config,
    match_allowed_journal,
    normalize_doi,
    request_json,
    source_blacklist_reason,
    write_csv,
    write_json,
)


def parse_crossref_item(item: dict, query_name: str, search_query: str, config: dict | None = None) -> dict:
    authors = []
    for author in item.get("author", []) or []:
        full_name = " ".join(
            part
            for part in [clean_text(author.get("given")), clean_text(author.get("family"))]
            if part
        )
        if full_name:
            authors.append(full_name)
    links = item.get("link", []) or []
    license_urls = [
        clean_text(entry.get("URL"))
        for entry in item.get("license", []) or []
        if clean_text(entry.get("URL"))
    ]
    pdf_url = ""
    fulltext_url = ""
    supplementary_urls = []
    for link in links:
        link_url = clean_text(link.get("URL"))
        content_type = clean_text(link.get("content-type")).lower()
        intended = clean_text(link.get("intended-application"))
        if not pdf_url and "pdf" in content_type:
            pdf_url = link_url
        if not fulltext_url and link_url:
            fulltext_url = link_url
        if link_url and (
            "supplement" in intended.lower()
            or is_supplementary_url(link_url, f"{content_type} {intended}")
        ):
            supplementary_urls.append(
                {
                    "url": link_url,
                    "label": intended or content_type or "supplementary_material",
                    "source": "crossref_link",
                }
            )
    if not fulltext_url:
        fulltext_url = clean_text((item.get("resource") or {}).get("primary", {}).get("URL"))
    year = None
    for date_key in ["published-print", "published-online", "issued", "created"]:
        parts = ((item.get(date_key) or {}).get("date-parts") or [])
        if parts and parts[0]:
            year = parts[0][0]
            break
    if year is None:
        year = extract_year((item.get("created") or {}).get("date-time"))
    doi = normalize_doi(item.get("DOI"))
    publisher = clean_text(item.get("publisher"))
    journal = clean_text(" ".join(item.get("container-title", []) or []))
    if source_blacklist_reason(journal=journal, doi=doi, publisher=publisher, config=config):
        return {}
    matched_journal = match_allowed_journal(journal, config)
    if config and config.get("restrict_to_allowed_journals", True) and not matched_journal:
        return {}
    return {
        "query_name": query_name,
        "search_query": search_query,
        "source_database": "crossref",
        "source_record_id": clean_text(item.get("DOI")) or clean_text(item.get("URL")),
        "title": clean_text(" ".join(item.get("title", []) or [])),
        "abstract": clean_text(item.get("abstract")),
        "authors": "; ".join(authors),
        "journal": journal,
        "matched_priority_journal": matched_journal,
        "year": year,
        "doi": doi,
        "pmid": "",
        "pmcid": "",
        "publisher": publisher,
        "publication_type": clean_text(item.get("type")),
        "open_access_flag": bool(license_urls or links),
        "fulltext_url": fulltext_url,
        "pdf_url": pdf_url,
        "landing_page_url": clean_text(item.get("URL")),
        "license": "; ".join(license_urls),
        "licenses_json": json_compatible(license_urls),
        "candidate_urls_json": json_compatible(
            [url for url in [pdf_url, fulltext_url, clean_text(item.get("URL"))] if url]
        ),
        "supplementary_urls_json": json_compatible(supplementary_urls),
    }


def search_crossref(config_path: str | None = None) -> pd.DataFrame:
    ensure_directories()
    config = load_config(config_path)
    if not config["sources"].get("crossref", True):
        return pd.DataFrame()

    rows: list[dict] = []
    logs: list[dict] = []
    max_results = int(config["max_results_per_query"].get("crossref", 100))
    page_size = min(int(config["page_size"].get("crossref", 100)), 100)
    delay = float(config["delay_seconds"].get("crossref", 0.5))
    headers = {"User-Agent": f"{USER_AGENT} (mailto:{config.get('email') or 'none@example.com'})"}
    date_filter = crossref_publication_date_filter(config)

    for query in config["queries"]:
        retrieved = 0
        first_payload = None
        while retrieved < max_results:
            request_size = min(page_size, max_results - retrieved)
            params = {
                "query.bibliographic": query["query"],
                "rows": request_size,
                "offset": retrieved,
            }
            if date_filter:
                params["filter"] = date_filter
            payload = request_json(
                CROSSREF_URL,
                params=params,
                headers=headers,
                delay_seconds=delay,
            )
            if first_payload is None:
                first_payload = payload
                write_json(payload, RAW_DIR / f"crossref_search_{query['name']}.json")
            items = payload.get("message", {}).get("items", []) or []
            if not items:
                break
            for item in items:
                row = parse_crossref_item(item, query["name"], query["query"], config)
                if row:
                    rows.append(row)
            retrieved += len(items)
            if len(items) < request_size:
                break
        logs.append(
            {
                "query_name": query["name"],
                "search_query": f"{query['query']} | filter:{date_filter}" if date_filter else query["query"],
                "db": "crossref",
                "retrieved_row_count": retrieved,
                "reported_total_results": clean_text((first_payload or {}).get("message", {}).get("total-results")),
            }
        )

    df = pd.DataFrame(rows)
    write_csv(df, RAW_DIR / "crossref_records.csv")
    write_csv(pd.DataFrame(logs), RAW_DIR / "crossref_query_log.csv")
    return df


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=None)
    args = parser.parse_args()
    search_crossref(args.config)


if __name__ == "__main__":
    main()
