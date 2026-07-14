from __future__ import annotations

import argparse
import math

import pandas as pd

from harvest_utils import (
    EUROPEPMC_SEARCH_URL,
    RAW_DIR,
    clean_text,
    ensure_directories,
    extract_year,
    is_supplementary_url,
    json_compatible,
    load_config,
    match_allowed_journal,
    europepmc_publication_date_query,
    normalize_doi,
    normalize_pmcid,
    normalize_pmid,
    request_json,
    source_blacklist_reason,
    write_csv,
    write_json,
)


def parse_result(result: dict, query_name: str, search_query: str, config: dict | None = None) -> dict:
    journal_title = clean_text(result.get("journalTitle")) or clean_text(
        result.get("journalInfo", {}).get("journal", {}).get("title")
    )
    publisher = clean_text(result.get("journalInfo", {}).get("journal", {}).get("publisher"))
    doi = normalize_doi(result.get("doi"))
    if source_blacklist_reason(journal=journal_title, doi=doi, publisher=publisher, config=config):
        return {}
    fulltext_urls = []
    for item in (result.get("fullTextUrlList") or {}).get("fullTextUrl", []):
        url = clean_text(item.get("url"))
        if url:
            fulltext_urls.append(
                {
                    "url": url,
                    "style": clean_text(item.get("documentStyle")),
                    "site": clean_text(item.get("site")),
                    "availability": clean_text(item.get("availability")),
                }
            )
    supplementary_urls = [
        {
            "url": item["url"],
            "label": item.get("style") or "supplementary_material",
            "source": item.get("site") or "europepmc",
        }
        for item in fulltext_urls
        if is_supplementary_url(item.get("url"), f"{item.get('style')} {item.get('availability')}")
    ]
    pmcid = normalize_pmcid(result.get("pmcid"))
    pdf_url = ""
    if clean_text(result.get("hasPDF")).upper() == "Y" and pmcid:
        pdf_url = f"https://europepmc.org/articles/{pmcid}?pdf=render"
    fulltext_url = ""
    if pmcid:
        fulltext_url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
    elif fulltext_urls:
        fulltext_url = fulltext_urls[0]["url"]
    matched_journal = match_allowed_journal(journal_title, config)
    if config and config.get("restrict_to_allowed_journals", True) and not matched_journal:
        return {}
    return {
        "query_name": query_name,
        "search_query": search_query,
        "source_database": "europepmc",
        "source_record_id": clean_text(result.get("id")),
        "title": clean_text(result.get("title")),
        "abstract": clean_text(result.get("abstractText")),
        "authors": clean_text(result.get("authorString")),
        "journal": journal_title,
        "matched_priority_journal": matched_journal,
        "year": extract_year(result.get("pubYear")),
        "doi": doi,
        "pmid": normalize_pmid(result.get("pmid")),
        "pmcid": pmcid,
        "publisher": publisher,
        "publication_type": clean_text(result.get("pubType")),
        "open_access_flag": clean_text(result.get("isOpenAccess")).upper() == "Y",
        "fulltext_url": fulltext_url,
        "pdf_url": pdf_url,
        "landing_page_url": f"https://europepmc.org/article/MED/{normalize_pmid(result.get('pmid'))}" if normalize_pmid(result.get("pmid")) else clean_text(fulltext_url),
        "license": clean_text(result.get("license")),
        "licenses_json": json_compatible([clean_text(result.get("license"))] if clean_text(result.get("license")) else []),
        "candidate_urls_json": json_compatible([item["url"] for item in fulltext_urls] + [url for url in [pdf_url, fulltext_url] if url]),
        "supplementary_urls_json": json_compatible(supplementary_urls),
    }


def search_europepmc(config_path: str | None = None) -> pd.DataFrame:
    ensure_directories()
    config = load_config(config_path)
    rows: list[dict] = []
    logs: list[dict] = []
    if not config["sources"].get("europepmc", True):
        return pd.DataFrame()

    max_results = config["max_results_per_query"].get("europepmc", 100)
    page_size = min(config["page_size"].get("europepmc", 100), 1000)
    delay = config["delay_seconds"].get("europepmc", 0.34)
    date_filter = europepmc_publication_date_query(config)

    for query in config["queries"]:
        retrieved = 0
        page = 1
        first_payload = None
        search_query = query["query"]
        if date_filter:
            search_query = f"({search_query}) AND {date_filter}"
        while retrieved < max_results:
            payload = request_json(
                EUROPEPMC_SEARCH_URL,
                params={
                    "query": search_query,
                    "format": "json",
                    "resultType": "core",
                    "pageSize": page_size,
                    "page": page,
                },
                delay_seconds=delay,
            )
            if first_payload is None:
                first_payload = payload
            result_list = payload.get("resultList", {}).get("result", [])
            if not result_list:
                break
            for result in result_list[: max_results - retrieved]:
                row = parse_result(result, query["name"], search_query, config)
                if row:
                    rows.append(row)
                retrieved += 1
            if len(result_list) < page_size:
                break
            page += 1
        hit_count = int((first_payload or {}).get("hitCount", 0))
        if first_payload is not None:
            write_json(first_payload, RAW_DIR / f"europepmc_search_{query['name']}.json")
        logs.append(
            {
                "query_name": query["name"],
                "search_query": search_query,
                "db": "europepmc",
                "retrieved_row_count": retrieved,
                "reported_hit_count": hit_count,
                "page_count_estimate": math.ceil(min(hit_count, max_results) / page_size) if page_size else 0,
            }
        )

    df = pd.DataFrame(rows)
    write_csv(df, RAW_DIR / "europepmc_records.csv")
    write_csv(pd.DataFrame(logs), RAW_DIR / "europepmc_query_log.csv")
    return df


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=None)
    args = parser.parse_args()
    search_europepmc(args.config)


if __name__ == "__main__":
    main()
