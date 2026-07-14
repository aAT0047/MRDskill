from __future__ import annotations

import csv
import html
import json
import os
import re
import ssl
import tarfile
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from html.parser import HTMLParser
from pathlib import Path
from typing import Any, Iterable

import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
HARVEST_DIR = (
    Path(os.environ.get("MRDSKILL_HARVEST_ROOT") or os.environ.get("LITERATURE_HARVEST_ROOT")).resolve()
    if os.environ.get("MRDSKILL_HARVEST_ROOT") or os.environ.get("LITERATURE_HARVEST_ROOT")
    else SCRIPT_DIR.parent
)
RAW_DIR = HARVEST_DIR / "raw_api_results"
MERGED_DIR = HARVEST_DIR / "merged_tables"
PDF_DIR = HARVEST_DIR / "downloaded_pdfs"
SUPPLEMENTARY_DIR = HARVEST_DIR / "supplementary_materials"
LOG_DIR = HARVEST_DIR / "download_logs"
REPORT_DIR = HARVEST_DIR / "reports"
DEFAULT_CONFIG_PATH = SCRIPT_DIR / "config_search_terms.yaml"
USER_AGENT = "oncology-literature-harvest/0.1"
REQUEST_TIMEOUT = 45

NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EUROPEPMC_SEARCH_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
EUROPEPMC_FULLTEXT_XML = "https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
CROSSREF_URL = "https://api.crossref.org/works"
OPENALEX_URL = "https://api.openalex.org/works"
UNPAYWALL_URL = "https://api.unpaywall.org/v2/{doi}"
PMC_OA_FILE_LIST = "https://pmc.ncbi.nlm.nih.gov/utils/oa/oa.fcgi"
PMC_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/pub/pmc/"

PRIORITY_JOURNAL_GROUPS = {
    "highest_priority_oncology_clinical_translational": [
        "Nature",
        "Nature Medicine",
        "Nature Cancer",
        "Nature Communications",
        "Cell",
        "Cancer Cell",
        "NEJM",
        "Lancet",
        "Lancet Oncology",
        "Journal of Clinical Oncology",
        "JAMA Oncology",
        "Annals of Oncology",
        "Cancer Discovery",
        "Clinical Cancer Research",
        "Cancer Research",
        "Blood",
        "Leukemia",
    ],
    "genomics_precision_oncology_molecular_diagnostics": [
        "Genome Medicine",
        "Nature Genetics",
        "Nature Biotechnology",
        "Nature Biomedical Engineering",
        "npj Precision Oncology",
        "Molecular Cancer",
        "Clinical Chemistry",
        "Journal of Molecular Diagnostics",
        "Modern Pathology",
        "Journal of Hematology & Oncology",
    ],
    "additional_oa_or_high_quality": [
        "ESMO Open",
        "BMC Medicine",
        "Genome Biology",
        "Communications Medicine",
        "Scientific Reports",
        "British Journal of Cancer",
        "Molecular Oncology",
        "Oncogene",
    ],
}

DEFAULT_ALLOWED_JOURNALS = [
    journal
    for journals in PRIORITY_JOURNAL_GROUPS.values()
    for journal in journals
]

DEFAULT_JOURNAL_TIERS = {
    journal: tier
    for tier, journals in enumerate(PRIORITY_JOURNAL_GROUPS.values(), start=1)
    for journal in journals
}

DEFAULT_JOURNAL_ALIASES = {
    "Nature": ["Nature"],
    "Nature Medicine": ["Nature Medicine", "Nat Med"],
    "Nature Cancer": ["Nature Cancer", "Nat Cancer"],
    "Nature Communications": ["Nature Communications", "Nat Commun"],
    "Cell": ["Cell"],
    "Cancer Cell": ["Cancer Cell"],
    "NEJM": ["NEJM", "The New England Journal of Medicine", "New England Journal of Medicine", "N Engl J Med"],
    "Lancet": ["Lancet", "The Lancet"],
    "Lancet Oncology": ["Lancet Oncology", "The Lancet Oncology", "The Lancet. Oncology", "Lancet Oncol"],
    "Journal of Clinical Oncology": ["Journal of Clinical Oncology", "J Clin Oncol", "JCO"],
    "JAMA Oncology": ["JAMA Oncology", "JAMA Oncol"],
    "Annals of Oncology": ["Annals of Oncology", "Ann Oncol"],
    "Cancer Discovery": ["Cancer Discovery", "Cancer Discov"],
    "Clinical Cancer Research": ["Clinical Cancer Research", "Clin Cancer Res"],
    "Cancer Research": ["Cancer Research", "Cancer Res"],
    "Blood": ["Blood"],
    "Leukemia": ["Leukemia"],
    "Genome Medicine": ["Genome Medicine", "Genome Med"],
    "Nature Genetics": ["Nature Genetics", "Nat Genet"],
    "Nature Biotechnology": ["Nature Biotechnology", "Nat Biotechnol"],
    "Nature Biomedical Engineering": ["Nature Biomedical Engineering", "Nat Biomed Eng"],
    "npj Precision Oncology": ["npj Precision Oncology", "NPJ Precision Oncology", "NPJ Precis Oncol"],
    "Molecular Cancer": ["Molecular Cancer", "Mol Cancer"],
    "Clinical Chemistry": ["Clinical Chemistry", "Clin Chem"],
    "Journal of Molecular Diagnostics": ["Journal of Molecular Diagnostics", "J Mol Diagn"],
    "Modern Pathology": ["Modern Pathology", "Mod Pathol"],
    "Journal of Hematology & Oncology": [
        "Journal of Hematology & Oncology",
        "Journal of Hematology and Oncology",
        "J Hematol Oncol",
    ],
    "ESMO Open": ["ESMO Open"],
    "BMC Medicine": ["BMC Medicine", "BMC Med"],
    "Genome Biology": ["Genome Biology", "Genome Biol"],
    "Communications Medicine": ["Communications Medicine", "Commun Med"],
    "Scientific Reports": ["Scientific Reports", "Sci Rep"],
    "British Journal of Cancer": ["British Journal of Cancer", "Br J Cancer"],
    "Molecular Oncology": ["Molecular Oncology", "Mol Oncol"],
    "Oncogene": ["Oncogene"],
}

DEFAULT_JOURNAL_POLICY_PATH = SCRIPT_DIR.parent.parent / "config" / "mrd_high_impact_journal_policy.json"

DEFAULT_JOURNAL_BLACKLIST = [
    {
        "rule_id": "BL-PUBLISHER-FRONTIERS",
        "status": "active",
        "scope": "publisher",
        "match_type": "contains",
        "pattern": "Frontiers Media",
        "reason": "User blacklist: Frontiers publisher excluded from candidate literature.",
    },
    {
        "rule_id": "BL-JOURNAL-FRONTIERS",
        "status": "active",
        "scope": "journal",
        "match_type": "contains",
        "pattern": "Frontiers",
        "reason": "User blacklist: Frontiers journals excluded from candidate literature.",
    },
    {
        "rule_id": "BL-DOI-FRONTIERS",
        "status": "active",
        "scope": "doi",
        "match_type": "prefix",
        "pattern": "10.3389",
        "reason": "User blacklist: Frontiers DOI prefix.",
    },
    {
        "rule_id": "BL-PUBLISHER-MDPI",
        "status": "active",
        "scope": "publisher",
        "match_type": "contains",
        "pattern": "MDPI",
        "reason": "User blacklist: MDPI publisher excluded from candidate literature.",
    },
    {
        "rule_id": "BL-DOI-MDPI",
        "status": "active",
        "scope": "doi",
        "match_type": "prefix",
        "pattern": "10.3390",
        "reason": "User blacklist: MDPI DOI prefix.",
    },
    {
        "rule_id": "BL-JOURNAL-CANCERS",
        "status": "active",
        "scope": "journal",
        "match_type": "exact",
        "pattern": "Cancers",
        "reason": "User blacklist: MDPI journal.",
    },
    {
        "rule_id": "BL-JOURNAL-CANCERS-BASEL",
        "status": "active",
        "scope": "journal",
        "match_type": "contains",
        "pattern": "Cancers (Basel)",
        "reason": "User blacklist: MDPI journal.",
    },
    {
        "rule_id": "BL-JOURNAL-DIAGNOSTICS",
        "status": "active",
        "scope": "journal",
        "match_type": "exact",
        "pattern": "Diagnostics",
        "reason": "User blacklist: MDPI journal.",
    },
    {
        "rule_id": "BL-JOURNAL-GENES",
        "status": "active",
        "scope": "journal",
        "match_type": "exact",
        "pattern": "Genes",
        "reason": "User blacklist: MDPI journal.",
    },
    {
        "rule_id": "BL-JOURNAL-BIOMEDICINES",
        "status": "active",
        "scope": "journal",
        "match_type": "exact",
        "pattern": "Biomedicines",
        "reason": "User blacklist: MDPI journal.",
    },
    {
        "rule_id": "BL-JOURNAL-JCM",
        "status": "active",
        "scope": "journal",
        "match_type": "exact",
        "pattern": "Journal of Clinical Medicine",
        "reason": "User blacklist: MDPI journal.",
    },
    {
        "rule_id": "BL-JOURNAL-IJMS",
        "status": "active",
        "scope": "journal",
        "match_type": "contains",
        "pattern": "International Journal of Molecular Sciences",
        "reason": "User blacklist: MDPI journal.",
    },
    {
        "rule_id": "BL-JOURNAL-CELLS",
        "status": "active",
        "scope": "journal",
        "match_type": "exact",
        "pattern": "Cells",
        "reason": "User blacklist: MDPI journal.",
    },
]

SUPPLEMENTARY_LINK_TERMS = [
    "supplement",
    "supplementary",
    "supplemental",
    "supporting information",
    "supporting-information",
    "supporting_information",
    "data supplement",
    "appendix",
    "additional file",
    "additional-file",
    "additional_file",
    "supplementary file",
    "supplementary files",
    "supplementary table",
    "supplementary tables",
    "supplementary data",
    "supplementary information",
    "supplemental table",
    "supplemental tables",
    "supplemental data",
    "supplemental information",
    "supplementary appendix",
    "web appendix",
    "source data",
    "reporting summary",
    "suppl",
    "etable",
    "sdataset",
    "mmc",
    "moesm",
    "mediaobjects",
]

SUPPLEMENTARY_ATTACHMENT_EXTENSIONS = {
    ".pdf",
    ".xlsx",
    ".xls",
    ".csv",
    ".tsv",
    ".docx",
    ".doc",
    ".zip",
    ".gz",
    ".tar",
    ".tar.gz",
    ".tgz",
    ".xml",
    ".txt",
    ".xlsm",
    ".xlsb",
    ".ods",
    ".rtf",
}

NON_ARTICLE_HTML_MARKERS = [
    "preparing to download",
    "redirecting",
    "captcha",
    "access denied",
    "enable javascript",
    "just a moment",
    "checking your browser",
]

NON_ARTICLE_HTML_HOSTS = [
    "openalex.org",
    "scholar.google.",
    "pubmed.ncbi.nlm.nih.gov",
    "pure.au.dk",
]

LINK_BLOCKLIST_HOSTS = [
    "scholar.google.",
    "openalex.org",
]

LINK_BLOCKLIST_PATTERNS = [
    "/about/copyright",
    "creativecommons.org/licenses",
    "scholar_lookup",
    "/references",
]

SUPPLEMENTARY_ROUTE_HINTS = [
    "/esm/",
    "/supplement",
    "/supplementary",
    "/suppinfo",
    "/suppl",
    "/supp/",
    "/mmc",
    "moesm",
    "_esm.",
    "_suppl",
    "static-content.springer.com",
    "ftp.ncbi.nlm.nih.gov/pub/pmc/",
    "pmc/oa_package/",
    "oup.silverchair-cdn.com",
    "article_deployments",
    "supplemental-data",
    "supporting-information",
    "additional_file",
    "data_supplement",
]

NON_SUPPLEMENTARY_LABEL_HINTS = [
    "figures/",
    "figure/",
    "og:image",
    "twitter:image",
    "image_src",
    "src img",
    "meta property",
    "meta name",
    "thumbnail",
    "cover image",
    "altmetric",
]

SUPPLEMENTARY_TEXT_TABLE_HINTS = [
    "table",
    "tables",
    "data",
    "dataset",
    "source data",
    "appendix",
    "methods",
    "protocol",
    "reporting summary",
    "checklist",
    "doc",
    "document",
    "spreadsheet",
    "excel",
    "csv",
    "tsv",
    "xlsx",
    "xls",
    "xml",
    "txt",
]

SUPPLEMENTARY_TEXT_TABLE_EXTENSIONS = {
    ".pdf",
    ".xlsx",
    ".xls",
    ".csv",
    ".tsv",
    ".docx",
    ".doc",
    ".zip",
    ".gz",
    ".tar",
    ".tar.gz",
    ".tgz",
    ".xml",
    ".txt",
    ".xlsm",
    ".xlsb",
    ".ods",
    ".rtf",
}

FIGURE_ONLY_HINTS = [
    "figure",
    "fig.",
    "fig ",
    "figures/",
    "figure/",
    "extended data fig",
    "supplementary figure",
    "supplemental figure",
    "sfig",
    "full size image",
    "view figure",
]

NON_ATTACHMENT_ROUTES = [
    "/figures/",
    "/figure/",
    "/tables/",
    "/table/",
    "/articles/",
]

ARTICLE_PDF_LABEL_TERMS = [
    "citation_pdf_url",
    "article pdf",
    "download pdf",
    "download-pdf",
    "full text pdf",
    "full-text pdf",
    "data-article-pdf",
]

NON_ARTICLE_PDF_HINTS = [
    "static-content.springer.com/esm",
    "/mediaobjects/",
    "moesm",
    "_esm.",
    "supplement",
    "supplementary",
    "supplemental",
    "supporting-information",
    "reporting summary",
    "appendix",
    "additional_file",
    "mmc",
]

ARTICLE_HTML_MARKERS = [
    "citation_title",
    "citation_journal_title",
    "dc.title",
    "article-title",
    "abstract",
    "full text",
    "main-content",
    "article__body",
    "article-body",
    "pmc-article",
    "<article",
]

def ensure_directories() -> None:
    for path in [RAW_DIR, MERGED_DIR, PDF_DIR, SUPPLEMENTARY_DIR, LOG_DIR, REPORT_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def load_config(config_path: str | Path | None = None) -> dict[str, Any]:
    path = Path(config_path) if config_path else DEFAULT_CONFIG_PATH
    config = json.loads(path.read_text(encoding="utf-8"))
    return apply_journal_policy(config, path.parent)


def write_json(data: Any, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def write_csv(df: pd.DataFrame, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, encoding="utf-8-sig")


def write_markdown(text: str, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def json_compatible(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True)


def clean_text(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and pd.isna(value):
        return ""
    text = html.unescape(str(value))
    text = re.sub(r"<[^>]+>", " ", text)
    text = re.sub(r"\s+", " ", text).strip()
    return text


def normalize_text(value: Any) -> str:
    return clean_text(value).lower()


def normalize_title(value: Any) -> str:
    text = normalize_text(value)
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def normalize_doi(value: Any) -> str:
    text = normalize_text(value)
    text = text.replace("https://doi.org/", "").replace("http://doi.org/", "")
    text = text.replace("doi:", "").strip().rstrip(".")
    return text


def normalize_pmid(value: Any) -> str:
    return re.sub(r"\D+", "", clean_text(value))


def normalize_pmcid(value: Any) -> str:
    text = clean_text(value).upper().replace("PMC", "")
    digits = re.sub(r"\D+", "", text)
    return f"PMC{digits}" if digits else ""


def extract_pmcid_from_text(value: Any) -> str:
    text = clean_text(value)
    match = re.search(r"\bPMC\d+\b", text, flags=re.IGNORECASE)
    if match:
        return normalize_pmcid(match.group(0))
    match = re.search(r"/articles/(\d{5,})", text, flags=re.IGNORECASE)
    if match:
        return normalize_pmcid(match.group(1))
    return ""


def extract_year(value: Any) -> int | None:
    text = clean_text(value)
    match = re.search(r"(19|20)\d{2}", text)
    return int(match.group(0)) if match else None


def chunked(items: list[Any], size: int) -> Iterable[list[Any]]:
    for start in range(0, len(items), size):
        yield items[start : start + size]


def request_json(
    url: str,
    params: dict[str, Any] | None = None,
    headers: dict[str, str] | None = None,
    delay_seconds: float = 0.34,
    retries: int = 3,
) -> Any:
    if params:
        url = f"{url}?{urllib.parse.urlencode(params, doseq=True)}"
    request_headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
    if headers:
        request_headers.update(headers)
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(url, headers=request_headers)
            with urllib.request.urlopen(request, timeout=REQUEST_TIMEOUT, context=ssl.create_default_context()) as response:
                payload = json.load(response)
            time.sleep(delay_seconds)
            return payload
        except Exception as exc:
            last_error = exc
            time.sleep(delay_seconds * attempt)
    raise RuntimeError(f"Failed request_json for {url}: {last_error}")


def request_text(
    url: str,
    params: dict[str, Any] | None = None,
    headers: dict[str, str] | None = None,
    delay_seconds: float = 0.34,
    retries: int = 3,
) -> str:
    if params:
        url = f"{url}?{urllib.parse.urlencode(params, doseq=True)}"
    request_headers = {"User-Agent": USER_AGENT}
    if headers:
        request_headers.update(headers)
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(url, headers=request_headers)
            with urllib.request.urlopen(request, timeout=REQUEST_TIMEOUT, context=ssl.create_default_context()) as response:
                payload = response.read().decode("utf-8", errors="replace")
            time.sleep(delay_seconds)
            return payload
        except Exception as exc:
            last_error = exc
            time.sleep(delay_seconds * attempt)
    raise RuntimeError(f"Failed request_text for {url}: {last_error}")


def request_bytes(
    url: str,
    params: dict[str, Any] | None = None,
    headers: dict[str, str] | None = None,
    delay_seconds: float = 0.34,
    retries: int = 3,
) -> tuple[bytes, str, str]:
    if params:
        url = f"{url}?{urllib.parse.urlencode(params, doseq=True)}"
    request_headers = {"User-Agent": USER_AGENT, "Accept": "*/*"}
    if headers:
        request_headers.update(headers)
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(url, headers=request_headers)
            with urllib.request.urlopen(request, timeout=REQUEST_TIMEOUT, context=ssl.create_default_context()) as response:
                payload = response.read()
                content_type = clean_text(response.headers.get("Content-Type")).lower()
                final_url = clean_text(response.geturl()) or url
            time.sleep(delay_seconds)
            return payload, content_type, final_url
        except Exception as exc:
            last_error = exc
            time.sleep(delay_seconds * attempt)
    raise RuntimeError(f"Failed request_bytes for {url}: {last_error}")


def flatten_authors(author_items: Iterable[Any]) -> str:
    names: list[str] = []
    for author in author_items or []:
        if isinstance(author, str):
            author_name = clean_text(author)
        elif isinstance(author, dict):
            author_name = clean_text(
                author.get("name")
                or author.get("fullName")
                or author.get("display_name")
                or " ".join(
                    part
                    for part in [clean_text(author.get("given")), clean_text(author.get("family"))]
                    if part
                )
            )
        else:
            author_name = clean_text(author)
        if author_name:
            names.append(author_name)
    return "; ".join(dict.fromkeys(names))


def reconstruct_openalex_abstract(inverted_index: dict[str, list[int]] | None) -> str:
    if not inverted_index:
        return ""
    positions: dict[int, str] = {}
    for token, values in inverted_index.items():
        for position in values:
            positions[position] = token
    return clean_text(" ".join(positions[index] for index in sorted(positions)))


def source_priority(source_name: str) -> int:
    order = {"pubmed": 5, "pmc": 4, "europepmc": 3, "crossref": 2, "openalex": 1}
    return order.get(source_name, 0)


def record_richness_score(row: pd.Series | dict[str, Any]) -> tuple[int, int]:
    data = row if isinstance(row, dict) else row.to_dict()
    completeness = 0
    for column in [
        "title",
        "abstract",
        "authors",
        "journal",
        "year",
        "doi",
        "pmid",
        "pmcid",
        "license",
        "fulltext_url",
        "pdf_url",
        "supplementary_urls_json",
    ]:
        if clean_text(data.get(column)):
            completeness += 1
    return completeness, source_priority(clean_text(data.get("source_database")))


def combine_unique(values: Iterable[Any], separator: str = "; ") -> str:
    cleaned = [clean_text(value) for value in values if clean_text(value)]
    unique = list(dict.fromkeys(cleaned))
    return separator.join(unique)


def combine_json_lists(values: Iterable[Any]) -> str:
    merged: list[Any] = []
    seen: set[str] = set()
    for value in values:
        if value is None or (isinstance(value, float) and pd.isna(value)):
            continue
        if isinstance(value, str):
            text = value.strip()
            if not text:
                continue
            try:
                parsed = json.loads(text)
                items = parsed if isinstance(parsed, list) else [parsed]
            except json.JSONDecodeError:
                items = [text]
        elif isinstance(value, list):
            items = value
        else:
            items = [value]
        for item in items:
            marker = json.dumps(item, ensure_ascii=False, sort_keys=True)
            if marker not in seen:
                seen.add(marker)
                merged.append(item)
    return json.dumps(merged, ensure_ascii=False)


def parse_json_list(value: Any) -> list[Any]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    if isinstance(value, list):
        return value
    if isinstance(value, str):
        text = value.strip()
        if not text:
            return []
        try:
            parsed = json.loads(text)
        except json.JSONDecodeError:
            return []
        return parsed if isinstance(parsed, list) else [parsed]
    return [value]


def normalize_journal_name(value: Any) -> str:
    text = normalize_text(value)
    text = text.replace("&", " and ")
    text = re.sub(r"\bthe\b", " ", text)
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def _year_or_none(value: Any) -> int | None:
    text = clean_text(value)
    if not text:
        return None
    try:
        year = int(text)
    except Exception:
        year = extract_year(text)
    return year if year and 1800 <= year <= 3000 else None


def publication_year_bounds(config: dict[str, Any] | None = None) -> tuple[int | None, int | None]:
    year_config = (config or {}).get("publication_year") or (config or {}).get("publication_year_range") or {}
    if not isinstance(year_config, dict):
        year_config = {}
    quota = year_quota_config(config)
    start = (
        year_config.get("start_year")
        or year_config.get("from_year")
        or year_config.get("min_year")
        or year_config.get("start")
        or quota.get("start_year")
    )
    end = (
        year_config.get("end_year")
        or year_config.get("to_year")
        or year_config.get("max_year")
        or year_config.get("end")
        or quota.get("end_year")
    )
    return _year_or_none(start), _year_or_none(end)


def publication_year_exclusion_reason(year: Any, config: dict[str, Any] | None = None) -> str:
    start_year, end_year = publication_year_bounds(config)
    if start_year is None and end_year is None:
        return ""
    parsed_year = _year_or_none(year)
    year_config = (config or {}).get("publication_year") or (config or {}).get("publication_year_range") or {}
    exclude_unknown = True
    if isinstance(year_config, dict) and "exclude_unknown_year" in year_config:
        exclude_unknown = bool(year_config.get("exclude_unknown_year"))
    if parsed_year is None:
        return "publication_year_missing" if exclude_unknown else ""
    if start_year is not None and parsed_year < start_year:
        return f"publication_year_before_{start_year}"
    if end_year is not None and parsed_year > end_year:
        return f"publication_year_after_{end_year}"
    return ""


def pubmed_publication_date_query(config: dict[str, Any] | None = None) -> str:
    start_year, end_year = publication_year_bounds(config)
    if start_year is None and end_year is None:
        return ""
    return f"{start_year or 1800}:{end_year or 3000}[dp]"


def europepmc_publication_date_query(config: dict[str, Any] | None = None) -> str:
    start_year, end_year = publication_year_bounds(config)
    if start_year is None and end_year is None:
        return ""
    return f"FIRST_PDATE:[{start_year or 1800}-01-01 TO {end_year or 3000}-12-31]"


def crossref_publication_date_filter(config: dict[str, Any] | None = None) -> str:
    start_year, end_year = publication_year_bounds(config)
    filters = []
    if start_year is not None:
        filters.append(f"from-pub-date:{start_year}-01-01")
    if end_year is not None:
        filters.append(f"until-pub-date:{end_year}-12-31")
    return ",".join(filters)


def openalex_publication_date_filter(config: dict[str, Any] | None = None) -> str:
    start_year, end_year = publication_year_bounds(config)
    filters = []
    if start_year is not None:
        filters.append(f"from_publication_date:{start_year}-01-01")
    if end_year is not None:
        filters.append(f"to_publication_date:{end_year}-12-31")
    return ",".join(filters)


def _policy_path(value: Any, base_dir: Path | None = None) -> Path:
    text = clean_text(value)
    if not text:
        return DEFAULT_JOURNAL_POLICY_PATH
    path = Path(text)
    return path if path.is_absolute() else (base_dir or SCRIPT_DIR) / path


def _load_journal_policy(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _flatten_journal_groups(groups: Any) -> tuple[list[str], dict[str, int]]:
    allowed: list[str] = []
    tiers: dict[str, int] = {}
    if not isinstance(groups, dict):
        return allowed, tiers
    for tier, journals in enumerate(groups.values(), start=1):
        if not isinstance(journals, list):
            continue
        for journal in journals:
            name = clean_text(journal)
            if not name or name in allowed:
                continue
            allowed.append(name)
            tiers[name] = tier
    return allowed, tiers


def apply_journal_policy(config: dict[str, Any], base_dir: Path | None = None) -> dict[str, Any]:
    policy_file = config.get("journal_policy_file") or config.get("journal_policy_path")
    policy = _load_journal_policy(_policy_path(policy_file, base_dir))
    if not policy:
        policy = _load_journal_policy(DEFAULT_JOURNAL_POLICY_PATH)
    if not policy:
        return config

    policy_allowed, policy_tiers = _flatten_journal_groups(policy.get("allowed_journal_groups"))
    if isinstance(policy.get("allowed_journals"), list):
        for journal in policy["allowed_journals"]:
            name = clean_text(journal)
            if name and name not in policy_allowed:
                policy_allowed.append(name)
                policy_tiers[name] = policy_tiers.get(name, 99)

    local_additions = [clean_text(item) for item in config.get("additional_allowed_journals", []) if clean_text(item)]
    for journal in local_additions:
        if journal not in policy_allowed:
            policy_allowed.append(journal)
            policy_tiers[journal] = 99

    if policy_allowed:
        config["allowed_journals"] = policy_allowed
    config["journal_tiers"] = {**policy_tiers, **(config.get("journal_tiers") or {})}

    policy_aliases = policy.get("journal_aliases") if isinstance(policy.get("journal_aliases"), dict) else {}
    config_aliases = config.get("journal_aliases") if isinstance(config.get("journal_aliases"), dict) else {}
    merged_aliases = dict(policy_aliases)
    for journal, values in config_aliases.items():
        existing = merged_aliases.get(journal, [])
        if not isinstance(existing, list):
            existing = [existing]
        incoming = values if isinstance(values, list) else [values]
        merged_aliases[journal] = list(dict.fromkeys([clean_text(v) for v in existing + incoming if clean_text(v)]))
    if merged_aliases:
        config["journal_aliases"] = merged_aliases

    policy_blacklist = policy.get("journal_blacklist") if isinstance(policy.get("journal_blacklist"), list) else []
    config_blacklist = config.get("journal_blacklist") if isinstance(config.get("journal_blacklist"), list) else []
    config["journal_blacklist"] = policy_blacklist + config_blacklist
    return config


def configured_allowed_journals(config: dict[str, Any] | None = None) -> list[str]:
    if config and config.get("allowed_journals"):
        return [clean_text(journal) for journal in config.get("allowed_journals", []) if clean_text(journal)]
    return DEFAULT_ALLOWED_JOURNALS.copy()


def configured_journal_tiers(config: dict[str, Any] | None = None) -> dict[str, int]:
    tiers = DEFAULT_JOURNAL_TIERS.copy()
    if config and isinstance(config.get("journal_tiers"), dict):
        for journal, tier in config["journal_tiers"].items():
            try:
                tiers[clean_text(journal)] = int(tier)
            except Exception:
                tiers[clean_text(journal)] = 99
    return tiers


def configured_journal_aliases(config: dict[str, Any] | None = None) -> dict[str, list[str]]:
    allowed = configured_allowed_journals(config)
    aliases: dict[str, list[str]] = {}
    for journal in allowed:
        values = DEFAULT_JOURNAL_ALIASES.get(journal, [journal])
        aliases[journal] = list(dict.fromkeys([journal] + values))
    extra_aliases = (config or {}).get("journal_aliases", {})
    if isinstance(extra_aliases, dict):
        for journal, values in extra_aliases.items():
            canonical = clean_text(journal)
            if not canonical:
                continue
            if isinstance(values, list):
                alias_values = [clean_text(value) for value in values if clean_text(value)]
            else:
                alias_values = [clean_text(values)]
            aliases[canonical] = list(dict.fromkeys(aliases.get(canonical, [canonical]) + alias_values))
    return aliases


def configured_journal_blacklist(config: dict[str, Any] | None = None) -> list[dict[str, Any]]:
    rules = DEFAULT_JOURNAL_BLACKLIST.copy()
    if config and isinstance(config.get("journal_blacklist"), list):
        rules.extend(rule for rule in config["journal_blacklist"] if isinstance(rule, dict))
    return rules


def _rule_value_matches(value: str, pattern: str, match_type: str, normalize_as_journal: bool = False) -> bool:
    if normalize_as_journal:
        value_norm = normalize_journal_name(value)
        pattern_norm = normalize_journal_name(pattern)
    else:
        value_norm = clean_text(value).lower()
        pattern_norm = clean_text(pattern).lower()
    if not value_norm or not pattern_norm:
        return False
    if match_type == "exact":
        return value_norm == pattern_norm
    if match_type == "prefix":
        return value_norm.startswith(pattern_norm)
    return pattern_norm in value_norm


def source_blacklist_reason(
    journal: Any = "",
    doi: Any = "",
    publisher: Any = "",
    config: dict[str, Any] | None = None,
) -> str:
    values = {
        "journal": clean_text(journal),
        "doi": normalize_doi(doi),
        "publisher": clean_text(publisher),
    }
    for rule in configured_journal_blacklist(config):
        if clean_text(rule.get("status") or "active").lower() not in {"", "active"}:
            continue
        scope = clean_text(rule.get("scope")).lower()
        pattern = clean_text(rule.get("pattern"))
        match_type = clean_text(rule.get("match_type") or "contains").lower()
        if scope not in values:
            continue
        if _rule_value_matches(values[scope], pattern, match_type, normalize_as_journal=(scope == "journal")):
            return f"{clean_text(rule.get('rule_id')) or 'blacklist'}:{clean_text(rule.get('reason'))}"
    return ""


def match_allowed_journal(journal: Any, config: dict[str, Any] | None = None) -> str:
    if source_blacklist_reason(journal=journal, config=config):
        return ""
    normalized = normalize_journal_name(journal)
    if not normalized:
        return ""
    for canonical, aliases in configured_journal_aliases(config).items():
        for alias in aliases:
            if normalized == normalize_journal_name(alias):
                return canonical
    return ""


def priority_journal_tier(journal: Any, config: dict[str, Any] | None = None) -> int:
    canonical = clean_text(journal)
    tiers = configured_journal_tiers(config)
    if canonical not in tiers:
        canonical = match_allowed_journal(canonical, config)
    if canonical in tiers:
        return tiers[canonical]
    if canonical and canonical in configured_allowed_journals(config):
        return 3
    return 99


def is_allowed_journal(journal: Any, config: dict[str, Any] | None = None) -> bool:
    if source_blacklist_reason(journal=journal, config=config):
        return False
    if config and not config.get("restrict_to_allowed_journals", True):
        return True
    return bool(match_allowed_journal(journal, config))


def priority_journal_query_filter(config: dict[str, Any] | None = None, field: str = "journal") -> str:
    query_filter_config = (config or {}).get("journal_query_filter", {})
    if isinstance(query_filter_config, dict) and query_filter_config.get("enabled") is False:
        return ""
    terms: list[str] = []
    for aliases in configured_journal_aliases(config).values():
        terms.extend(aliases)
    unique = list(dict.fromkeys(clean_text(term) for term in terms if clean_text(term)))
    max_terms = 0
    if isinstance(query_filter_config, dict):
        try:
            max_terms = int(query_filter_config.get("max_terms") or 0)
        except Exception:
            max_terms = 0
    if max_terms > 0 and len(unique) > max_terms:
        return ""
    if not unique:
        return ""
    if field.lower() in {"pubmed_journal", "pubmed"}:
        return " OR ".join(f'"{term}"[Journal]' for term in unique)
    return " OR ".join(f'{field}:"{term}"' for term in unique)


def is_supplementary_url(url: Any, label: Any = "") -> bool:
    text = normalize_text(f"{url} {label}")
    return any(term in text for term in SUPPLEMENTARY_LINK_TERMS)


def _url_host(url: str) -> str:
    return urllib.parse.urlparse(clean_text(url)).netloc.lower()


def _url_path_query(url: str) -> str:
    parsed = urllib.parse.urlparse(clean_text(url))
    return f"{parsed.netloc}{parsed.path}?{parsed.query}".lower()


def _has_host(url: str, host_fragments: Iterable[str]) -> bool:
    host = _url_host(url)
    return any(fragment in host for fragment in host_fragments)


def is_blocked_link_url(url: Any) -> bool:
    text = clean_text(url)
    lowered = text.lower()
    if not text or lowered.startswith(("#", "mailto:", "javascript:", "tel:")):
        return True
    if _has_host(text, LINK_BLOCKLIST_HOSTS):
        return True
    return any(pattern in lowered for pattern in LINK_BLOCKLIST_PATTERNS)


def has_supplementary_hint(url: Any, label: Any = "") -> bool:
    text = normalize_text(f"{url} {label}")
    return any(term in text for term in SUPPLEMENTARY_LINK_TERMS) or any(
        hint in text for hint in SUPPLEMENTARY_ROUTE_HINTS
    )


def looks_like_attachment_route(url: Any) -> bool:
    text = _url_path_query(clean_text(url))
    if any(route in text for route in NON_ATTACHMENT_ROUTES) and not any(
        hint in text for hint in SUPPLEMENTARY_ROUTE_HINTS
    ):
        return False
    return any(hint in text for hint in SUPPLEMENTARY_ROUTE_HINTS) or "download" in text


def is_likely_article_html(payload: bytes | str, url: str = "") -> bool:
    text = payload.decode("utf-8", errors="ignore") if isinstance(payload, bytes) else payload
    if _has_host(url, NON_ARTICLE_HTML_HOSTS):
        return False
    lowered = normalize_text(text[:60000])
    if any(marker in lowered for marker in NON_ARTICLE_HTML_MARKERS):
        return False
    return any(marker in lowered for marker in ARTICLE_HTML_MARKERS)


def is_probably_redirect_or_challenge_html(payload: bytes | str, url: str = "") -> bool:
    text = payload.decode("utf-8", errors="ignore") if isinstance(payload, bytes) else payload
    lowered = normalize_text(f"{url} {text[:20000]}")
    return any(marker in lowered for marker in NON_ARTICLE_HTML_MARKERS)


def supplementary_extension_from_url(url: str, default: str = ".bin") -> str:
    lowered = urllib.parse.urlparse(url).path.lower()
    if lowered.endswith(".tar.gz"):
        return ".tar.gz"
    for extension in [
        ".pdf",
        ".xlsx",
        ".xls",
        ".csv",
        ".tsv",
        ".docx",
        ".doc",
        ".zip",
        ".gz",
        ".tar",
        ".tgz",
        ".xml",
        ".txt",
        ".ppt",
        ".pptx",
        ".xlsm",
        ".xlsb",
        ".ods",
        ".rtf",
    ]:
        if lowered.endswith(extension):
            return extension
    return default


def attachment_extension(content_type: str, url: str, payload: bytes | None = None) -> str:
    lowered_type = clean_text(content_type).lower()
    lowered_url = clean_text(url).lower()
    if lowered_url.endswith(".tar.gz"):
        return ".tar.gz"
    if payload and payload.startswith(b"%PDF"):
        return ".pdf"
    if "pdf" in lowered_type:
        return ".pdf"
    if "spreadsheetml" in lowered_type or "excel" in lowered_type:
        return ".xlsx"
    if "msword" in lowered_type:
        return ".doc"
    if "wordprocessingml" in lowered_type:
        return ".docx"
    if "gzip" in lowered_type or "x-gzip" in lowered_type:
        return supplementary_extension_from_url(url, default=".tar.gz")
    if "zip" in lowered_type:
        return ".zip"
    if "csv" in lowered_type:
        return ".csv"
    if "tab-separated" in lowered_type or "tsv" in lowered_type:
        return ".tsv"
    if "xml" in lowered_type:
        return ".xml"
    if "plain" in lowered_type:
        return supplementary_extension_from_url(url, default=".txt")
    if "html" in lowered_type or lowered_type.startswith("text/"):
        return ".html"
    return supplementary_extension_from_url(url, default=content_extension(content_type, url))


def supplementary_download_preference_rank(item: dict[str, Any]) -> tuple[int, int, int, int, str]:
    url = clean_text(item.get("url"))
    label = clean_text(item.get("label"))
    combined = normalize_text(f"{url} {label}")
    extension = supplementary_extension_from_url(url, default="")
    table_extensions = {".xlsx", ".xls", ".xlsm", ".xlsb", ".ods", ".csv", ".tsv"}
    table_hints = ["table", "tables", "data", "dataset", "source data", "spreadsheet", "excel", "csv", "tsv"]
    figure_only = any(hint in combined for hint in FIGURE_ONLY_HINTS) and not any(hint in combined for hint in table_hints)
    table_rank = 0 if extension in table_extensions or any(hint in combined for hint in table_hints) else 1
    direct_attachment_rank = 0 if extension in SUPPLEMENTARY_TEXT_TABLE_EXTENSIONS else 1
    figure_rank = 1 if figure_only else 0
    html_route_rank = 1 if "pmc.ncbi.nlm.nih.gov/articles/instance/" in url.lower() and "/bin/" in url.lower() else 0
    return (table_rank, direct_attachment_rank, figure_rank, html_route_rank, url.lower())


def is_supplementary_attachment(url: str, label: str = "", content_type: str = "") -> bool:
    if is_blocked_link_url(url):
        return False
    normalized_label = normalize_text(label)
    normalized_url = normalize_text(url)
    combined = f"{normalized_url} {normalized_label}"
    if any(hint in combined for hint in NON_SUPPLEMENTARY_LABEL_HINTS):
        return False
    extension = supplementary_extension_from_url(url, default="")
    has_text_table_hint = any(hint in combined for hint in SUPPLEMENTARY_TEXT_TABLE_HINTS)
    if any(hint in combined for hint in FIGURE_ONLY_HINTS) and not has_text_table_hint:
        return False
    if extension and extension not in SUPPLEMENTARY_TEXT_TABLE_EXTENSIONS:
        return False
    lowered_type = clean_text(content_type).lower()
    type_is_attachment = any(
        token in lowered_type
        for token in ["spreadsheet", "excel", "msword", "wordprocessingml", "zip", "pdf", "csv", "xml", "text/plain"]
    )
    if extension in SUPPLEMENTARY_ATTACHMENT_EXTENSIONS:
        return has_supplementary_hint(url, label) and has_text_table_hint
    if type_is_attachment:
        return (has_supplementary_hint(url, label) or looks_like_attachment_route(url)) and has_text_table_hint
    if has_supplementary_hint(url, label) and looks_like_attachment_route(url) and has_text_table_hint:
        return True
    return False


class _HtmlLinkExtractor(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.links: list[dict[str, str]] = []
        self._anchor_stack: list[dict[str, Any]] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attr_map = {name.lower(): value or "" for name, value in attrs}
        label_parts = [
            attr_map.get(name, "")
            for name in [
                "title",
                "aria-label",
                "download",
                "data-track-label",
                "data-track-action",
                "data-test",
                "name",
                "property",
                "rel",
            ]
            if attr_map.get(name)
        ]
        for attr_name in ["href", "src", "data-href", "data-url", "data-download-url", "content"]:
            raw_url = html.unescape(attr_map.get(attr_name, "").strip())
            if not raw_url:
                continue
            item: dict[str, Any] = {
                "url": raw_url,
                "label_parts": label_parts + [attr_name, tag],
                "text_parts": [],
            }
            self.links.append(item)
            if tag == "a":
                self._anchor_stack.append(item)

    def handle_data(self, data: str) -> None:
        if self._anchor_stack and data.strip():
            self._anchor_stack[-1]["text_parts"].append(data)

    def handle_endtag(self, tag: str) -> None:
        if tag == "a" and self._anchor_stack:
            self._anchor_stack.pop()

    def normalized_links(self) -> list[dict[str, str]]:
        output: list[dict[str, str]] = []
        for item in self.links:
            label = clean_text(" ".join(item.get("label_parts", []) + item.get("text_parts", [])))
            output.append({"url": clean_text(item.get("url")), "label": label})
        return output


def html_links(text: str) -> list[dict[str, str]]:
    parser = _HtmlLinkExtractor()
    try:
        parser.feed(text or "")
    except Exception:  # noqa: BLE001
        return []
    return parser.normalized_links()


def extract_supplementary_links(text: str, bases: list[str] | None = None) -> list[dict[str, str]]:
    values: list[dict[str, str]] = []
    source_text = text or ""
    for item in html_links(source_text):
        raw_url = item["url"]
        label = item["label"]
        if is_blocked_link_url(raw_url):
            continue
        for base in bases or [""]:
            resolved = urllib.parse.urljoin(base, raw_url) if base else raw_url
            if resolved.startswith("http://") or resolved.startswith("https://"):
                if not is_supplementary_attachment(resolved, label):
                    continue
                values.append({"url": resolved, "label": clean_text(label)[:140], "source": "article_page"})
                break
    return unique_url_items(values)


def is_likely_article_pdf_url(url: str, label: str = "", doi: Any = "") -> bool:
    lowered_url = clean_text(url).lower()
    lowered_label = normalize_text(label)
    if not lowered_url or is_blocked_link_url(lowered_url):
        return False
    if not (
        lowered_url.endswith(".pdf")
        or ".pdf?" in lowered_url
        or "/pdf/" in lowered_url
        or "pdfdirect" in lowered_url
        or "pdf=render" in lowered_url
    ):
        return False
    if any(hint in lowered_url or hint in lowered_label for hint in NON_ARTICLE_PDF_HINTS):
        return False
    if any(term in lowered_label for term in ARTICLE_PDF_LABEL_TERMS):
        return True
    parsed = urllib.parse.urlparse(lowered_url)
    if parsed.netloc.endswith("pmc.ncbi.nlm.nih.gov") and "/pdf/" in parsed.path:
        return True
    doi_suffix = normalize_doi(doi).split("/", 1)[-1]
    if doi_suffix and doi_suffix.replace("/", "-").lower() in lowered_url.replace("%2f", "-"):
        return True
    return False


def extract_article_pdf_links(text: str, bases: list[str] | None = None, doi: Any = "") -> list[str]:
    links: list[str] = []
    seen: set[str] = set()
    for item in html_links(text or ""):
        raw_url = item["url"]
        label = item["label"]
        if is_blocked_link_url(raw_url):
            continue
        for base in bases or [""]:
            resolved = urllib.parse.urljoin(base, raw_url) if base else raw_url
            if not (resolved.startswith("http://") or resolved.startswith("https://")):
                continue
            if not is_likely_article_pdf_url(resolved, label, doi):
                continue
            marker = resolved.lower()
            if marker not in seen:
                seen.add(marker)
                links.append(resolved)
            break
    return links


def normalize_pmc_oa_href(href: Any) -> str:
    url = clean_text(href)
    if not url:
        return ""
    if url.startswith("ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/"):
        return PMC_FTP_BASE + url.removeprefix("ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/")
    if url.startswith("ftp://"):
        parsed = urllib.parse.urlparse(url)
        if parsed.netloc == "ftp.ncbi.nlm.nih.gov" and parsed.path.startswith("/pub/pmc/"):
            return f"https://ftp.ncbi.nlm.nih.gov{parsed.path}"
    if url.startswith("/"):
        return urllib.parse.urljoin("https://pmc.ncbi.nlm.nih.gov/", url)
    return url


def unpaywall_email(config: dict[str, Any] | None = None) -> str:
    external = (config or {}).get("external_fulltext", {})
    return (
        clean_text(external.get("unpaywall_email"))
        or clean_text((config or {}).get("email"))
        or clean_text(os.environ.get("UNPAYWALL_EMAIL"))
    )


def external_fulltext_enabled(config: dict[str, Any] | None, key: str, default: bool = True) -> bool:
    external = (config or {}).get("external_fulltext", {})
    if not external.get("enabled", True):
        return False
    return bool(external.get(key, default))


def unpaywall_locations_for_doi(doi: Any, config: dict[str, Any] | None = None) -> dict[str, Any]:
    normalized = normalize_doi(doi)
    email = unpaywall_email(config)
    if not normalized or not email or not external_fulltext_enabled(config, "unpaywall", True):
        return {"candidate_urls": [], "supplementary_urls": [], "metadata": {}}
    try:
        payload = request_json(
            UNPAYWALL_URL.format(doi=urllib.parse.quote(normalized, safe="/")),
            params={"email": email},
            delay_seconds=float((config or {}).get("delay_seconds", {}).get("unpaywall", 0.34)),
        )
    except Exception:  # noqa: BLE001
        return {"candidate_urls": [], "supplementary_urls": [], "metadata": {}}

    candidate_items: list[dict[str, str]] = []
    supplementary_items: list[dict[str, str]] = []
    locations = []
    best = payload.get("best_oa_location")
    if isinstance(best, dict):
        locations.append(best)
    locations.extend(item for item in payload.get("oa_locations", []) if isinstance(item, dict))
    for index, location in enumerate(locations, start=1):
        source = "unpaywall_best_oa_location" if index == 1 and location is best else "unpaywall_oa_location"
        for key, label in [
            ("url_for_pdf", "unpaywall pdf"),
            ("url", "unpaywall url"),
            ("url_for_landing_page", "unpaywall landing page"),
        ]:
            url = clean_text(location.get(key))
            if not url.startswith(("http://", "https://")):
                continue
            item = {
                "url": url,
                "label": f"{label}; host_type={clean_text(location.get('host_type'))}; license={clean_text(location.get('license'))}",
                "source": source,
            }
            if is_supplementary_attachment(url, item["label"]):
                supplementary_items.append(item)
            else:
                candidate_items.append(item)

    return {
        "candidate_urls": unique_url_items(candidate_items),
        "supplementary_urls": unique_url_items(supplementary_items),
        "metadata": {
            "is_oa": bool(payload.get("is_oa")),
            "oa_status": clean_text(payload.get("oa_status")),
            "has_repository_copy": bool(payload.get("has_repository_copy")),
        },
    }


def pmc_oa_links_for_pmcid(pmcid: Any, config: dict[str, Any] | None = None) -> dict[str, Any]:
    normalized = normalize_pmcid(pmcid)
    if not normalized or not external_fulltext_enabled(config, "pmc_oa", True):
        return {"candidate_urls": [], "supplementary_urls": [], "package_urls": [], "metadata": {}}
    try:
        xml_text = request_text(
            PMC_OA_FILE_LIST,
            params={"id": normalized},
            delay_seconds=float((config or {}).get("delay_seconds", {}).get("pmc_oa", 0.34)),
        )
    except Exception:  # noqa: BLE001
        return {"candidate_urls": [], "supplementary_urls": [], "package_urls": [], "metadata": {}}

    candidate_items: list[dict[str, str]] = []
    supplementary_items: list[dict[str, str]] = []
    package_items: list[dict[str, str]] = []
    license_text = ""
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return {"candidate_urls": [], "supplementary_urls": [], "package_urls": [], "metadata": {}}

    record = root.find(".//record")
    if record is not None:
        license_text = clean_text(record.attrib.get("license"))
    for link in root.findall(".//link"):
        href = normalize_pmc_oa_href(link.attrib.get("href"))
        if not href.startswith(("http://", "https://")):
            continue
        fmt = clean_text(link.attrib.get("format")).lower()
        updated = clean_text(link.attrib.get("updated"))
        label = f"pmc oa {fmt}; license={license_text}; updated={updated}"
        item = {"url": href, "label": label, "source": "pmc_oa"}
        if fmt in {"tgz", "tar.gz"} or href.lower().endswith((".tar.gz", ".tgz")):
            package_items.append(item)
        elif fmt == "pdf" or href.lower().endswith(".pdf"):
            candidate_items.append(item)
        elif is_supplementary_attachment(href, label):
            supplementary_items.append(item)
        else:
            candidate_items.append(item)

    return {
        "candidate_urls": unique_url_items(candidate_items),
        "supplementary_urls": unique_url_items(supplementary_items),
        "package_urls": unique_url_items(package_items),
        "metadata": {"license": license_text, "pmcid": normalized},
    }


def enrich_candidate_urls_with_external_sources(
    row: pd.Series | dict[str, Any],
    config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    data = row if isinstance(row, dict) else row.to_dict()
    candidate_items: list[dict[str, str]] = []
    supplementary_items: list[dict[str, str]] = []
    package_items: list[dict[str, str]] = []
    sources: list[str] = []

    unpaywall = unpaywall_locations_for_doi(data.get("doi"), config)
    if unpaywall["candidate_urls"] or unpaywall["supplementary_urls"]:
        sources.append("unpaywall")
    candidate_items.extend(unpaywall["candidate_urls"])
    supplementary_items.extend(unpaywall["supplementary_urls"])

    pmcid = (
        normalize_pmcid(data.get("pmcid"))
        or extract_pmcid_from_text(data.get("fulltext_url"))
        or extract_pmcid_from_text(data.get("pdf_url"))
        or extract_pmcid_from_text(data.get("landing_page_url"))
        or extract_pmcid_from_text(data.get("candidate_urls_json"))
    )
    pmc_oa = pmc_oa_links_for_pmcid(pmcid, config)
    if pmc_oa["candidate_urls"] or pmc_oa["supplementary_urls"] or pmc_oa["package_urls"]:
        sources.append("pmc_oa")
    candidate_items.extend(pmc_oa["candidate_urls"])
    supplementary_items.extend(pmc_oa["supplementary_urls"])
    package_items.extend(pmc_oa["package_urls"])

    return {
        "candidate_urls": unique_url_items(candidate_items),
        "supplementary_urls": unique_url_items(supplementary_items),
        "package_urls": unique_url_items(package_items),
        "sources": "; ".join(dict.fromkeys(sources)),
    }


def append_external_links_to_row(
    row: pd.Series | dict[str, Any],
    config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    data = row.to_dict() if isinstance(row, pd.Series) else dict(row)
    enrichment = enrich_candidate_urls_with_external_sources(data, config)
    current_candidates = unique_url_items(parse_json_list(data.get("candidate_urls_json")))
    current_supplementary = parse_supplementary_urls(data.get("supplementary_urls_json"))
    current_packages = unique_url_items(parse_json_list(data.get("pmc_oa_package_urls_json")))
    data["candidate_urls_json"] = json_compatible(current_candidates + enrichment["candidate_urls"])
    data["supplementary_urls_json"] = json_compatible(current_supplementary + enrichment["supplementary_urls"])
    data["pmc_oa_package_urls_json"] = json_compatible(current_packages + enrichment["package_urls"])
    data["external_fulltext_sources"] = enrichment["sources"]
    data["pmcid"] = (
        normalize_pmcid(data.get("pmcid"))
        or extract_pmcid_from_text(data.get("fulltext_url"))
        or extract_pmcid_from_text(data.get("pdf_url"))
        or extract_pmcid_from_text(data.get("landing_page_url"))
        or extract_pmcid_from_text(data.get("candidate_urls_json"))
    )
    if enrichment["candidate_urls"] and not clean_text(data.get("pdf_url")):
        first_pdf = next(
            (item["url"] for item in enrichment["candidate_urls"] if is_likely_article_pdf_url(item["url"], item.get("label", ""), data.get("doi"))),
            "",
        )
        if first_pdf:
            data["pdf_url"] = first_pdf
    if enrichment["candidate_urls"] and not clean_text(data.get("fulltext_url")):
        data["fulltext_url"] = enrichment["candidate_urls"][0]["url"]
    return data


def pmc_oa_package_urls(row: pd.Series | dict[str, Any]) -> list[dict[str, str]]:
    data = row if isinstance(row, dict) else row.to_dict()
    return unique_url_items(parse_json_list(data.get("pmc_oa_package_urls_json")))


def safe_extract_tar_members(package_path: Path, output_dir: Path) -> list[Path]:
    extracted: list[Path] = []
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(package_path, "r:*") as archive:
        for member in archive.getmembers():
            if not member.isfile():
                continue
            member_name = member.name.replace("\\", "/").lstrip("/")
            if ".." in Path(member_name).parts:
                continue
            target = (output_dir / member_name).resolve()
            if output_dir not in target.parents and target != output_dir:
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            file_obj = archive.extractfile(member)
            if file_obj is None:
                continue
            target.write_bytes(file_obj.read())
            extracted.append(target)
    return extracted


def record_id_from_filename(filename: str) -> str:
    name = Path(filename).name
    if "__" in name:
        return name.split("__", 1)[0]
    return name.split("_", 1)[0]


def paper_file_stem(row: pd.Series | dict[str, Any], record_id: str | None = None) -> str:
    data = row if isinstance(row, dict) else row.to_dict()
    rid = safe_filename(clean_text(record_id or data.get("record_id")) or "record", max_len=40)
    journal = clean_text(data.get("matched_priority_journal")) or clean_text(data.get("journal")) or "unknown_journal"
    year = clean_text(data.get("year")) or "unknown_year"
    author = first_author(data.get("authors"))
    title = clean_text(data.get("title")) or clean_text(data.get("doi")) or "untitled"
    parts = [
        rid,
        "paper",
        safe_filename(journal.lower(), 36),
        safe_filename(year, 12),
        author,
        safe_filename(title.lower(), 70),
    ]
    return "__".join(parts)[:180].strip("._")


def supplementary_file_stem(
    row: pd.Series | dict[str, Any],
    index: int,
    url: str,
    label: str = "",
    record_id: str | None = None,
) -> str:
    data = row if isinstance(row, dict) else row.to_dict()
    rid = safe_filename(clean_text(record_id or data.get("record_id")) or "record", max_len=40)
    label_text = clean_text(label) or Path(urllib.parse.urlparse(url).path).stem or "supplementary_material"
    return "__".join([rid, f"supp{index:02d}", safe_filename(label_text.lower(), 80)])[:160].strip("._")


def unique_url_items(items: Iterable[Any]) -> list[dict[str, str]]:
    output: list[dict[str, str]] = []
    seen: set[str] = set()
    for item in items:
        if isinstance(item, str):
            url = clean_text(item)
            label = ""
            source = ""
        elif isinstance(item, dict):
            url = clean_text(item.get("url") or item.get("URL"))
            label = clean_text(item.get("label") or item.get("type") or item.get("name") or item.get("title"))
            source = clean_text(item.get("source") or item.get("site") or item.get("availability"))
        else:
            continue
        if not (url.startswith("http://") or url.startswith("https://")):
            continue
        marker = url.lower()
        if marker in seen:
            continue
        seen.add(marker)
        output.append({"url": url, "label": label, "source": source})
    return output


def parse_supplementary_urls(value: Any) -> list[dict[str, str]]:
    return unique_url_items(parse_json_list(value))


def year_quota_config(config: dict[str, Any] | None = None) -> dict[str, Any]:
    download_cfg = (config or {}).get("download", {})
    quota = download_cfg.get("year_quota", {})
    if not quota:
        quota = {
            "enabled": bool(download_cfg.get("year_quota_enabled", False)),
            "per_year": download_cfg.get("max_records_per_year", 0),
            "years": download_cfg.get("year_quota_years", 0),
            "start_year": download_cfg.get("year_quota_start_year", ""),
            "end_year": download_cfg.get("year_quota_end_year", ""),
        }
    return quota if isinstance(quota, dict) else {}


def apply_year_quota(df: pd.DataFrame, config: dict[str, Any] | None = None) -> pd.DataFrame:
    quota = year_quota_config(config)
    if df.empty or not quota.get("enabled", False):
        return df
    per_year = int(quota.get("per_year", 0) or 0)
    if per_year <= 0:
        return df

    data = df.copy()
    data["_quota_year"] = pd.to_numeric(data.get("year", pd.Series(pd.NA, index=data.index)), errors="coerce").astype("Int64")
    data = data.loc[data["_quota_year"].notna()].copy()
    if data.empty:
        return data.drop(columns=["_quota_year"], errors="ignore")

    start_year = quota.get("start_year", "")
    end_year = quota.get("end_year", "")
    if clean_text(start_year):
        data = data.loc[data["_quota_year"].astype(int).ge(int(start_year))].copy()
    if clean_text(end_year):
        data = data.loc[data["_quota_year"].astype(int).le(int(end_year))].copy()
    years = int(quota.get("years", 0) or 0)
    if years > 0 and not clean_text(start_year):
        max_year = int(data["_quota_year"].max()) if not data.empty else 0
        min_year = max_year - years + 1
        data = data.loc[data["_quota_year"].astype(int).between(min_year, max_year)].copy()

    if data.empty:
        return data.drop(columns=["_quota_year"], errors="ignore")
    selected = data.groupby("_quota_year", group_keys=False, sort=False).head(per_year)
    return selected.drop(columns=["_quota_year"], errors="ignore").reset_index(drop=True)


def pick_first_url(values: Iterable[Any]) -> str:
    for value in values:
        text = clean_text(value)
        if text.startswith("http://") or text.startswith("https://"):
            return text
    return ""


def merge_source_rows(group: pd.DataFrame) -> dict[str, Any]:
    ranked = group.copy()
    ranked["_completeness"] = ranked.apply(lambda row: record_richness_score(row)[0], axis=1)
    ranked["_source_priority"] = ranked["source_database"].map(source_priority).fillna(0)
    ranked = ranked.sort_values(
        ["_completeness", "_source_priority", "open_access_flag"],
        ascending=[False, False, False],
        na_position="last",
    )
    primary = ranked.iloc[0].to_dict()
    merged = {key: value for key, value in primary.items() if not str(key).startswith("_")}
    for column in ["title", "abstract", "authors", "journal", "publication_type", "license", "matched_priority_journal"]:
        merged[column] = combine_unique(group[column].tolist(), separator=" || ")
    for column in ["doi", "pmid", "pmcid"]:
        cleaned_values = [clean_text(value) for value in group[column].tolist() if clean_text(value)]
        merged[column] = cleaned_values[0] if cleaned_values else ""
    merged["source_databases"] = combine_unique(group["source_database"].tolist())
    merged["query_names"] = combine_unique(group["query_name"].tolist())
    merged["search_queries"] = combine_unique(group["search_query"].tolist(), separator=" || ")
    merged["source_record_ids"] = combine_unique(group["source_record_id"].tolist())
    merged["open_access_flag"] = bool(group["open_access_flag"].fillna(False).any())
    merged["pdf_url"] = pick_first_url(group["pdf_url"].tolist())
    merged["fulltext_url"] = pick_first_url(group["fulltext_url"].tolist())
    merged["landing_page_url"] = pick_first_url(group["landing_page_url"].tolist())
    merged["candidate_urls_json"] = combine_json_lists(group["candidate_urls_json"].tolist())
    merged["supplementary_urls_json"] = combine_json_lists(group["supplementary_urls_json"].tolist())
    merged["licenses_json"] = combine_json_lists(group["licenses_json"].tolist())
    merged["matched_sources_count"] = int(group["source_database"].nunique())
    merged["matched_rows_count"] = int(len(group))
    return merged


def likely_review(title: Any, publication_type: Any) -> bool:
    text = normalize_text(title) + " " + normalize_text(publication_type)
    return any(token in text for token in ["review", "systematic review", "meta-analysis"])


def likely_original_research(title: Any, publication_type: Any) -> bool:
    if likely_review(title, publication_type):
        return False
    text = normalize_text(publication_type)
    if not text:
        return True
    if any(token in text for token in ["journal article", "research article", "article"]):
        return True
    return not any(token in text for token in ["editorial", "comment", "news"])


def classify_priority(score: int) -> str:
    if score >= 5:
        return "high"
    if score >= 2:
        return "medium"
    return "low"


def safe_filename(text: str, max_len: int = 120) -> str:
    cleaned = re.sub(r"[^0-9A-Za-z._-]+", "_", clean_text(text))
    cleaned = re.sub(r"_+", "_", cleaned).strip("._")
    return cleaned[:max_len] or "record"


def first_author(authors_text: Any) -> str:
    text = clean_text(authors_text)
    if not text:
        return "unknown_author"
    first = text.split(";")[0].strip()
    surname = re.split(r"[\s,]+", first)[0]
    return safe_filename(surname.lower(), max_len=40)


def stable_file_stem(row: pd.Series | dict[str, Any]) -> str:
    data = row if isinstance(row, dict) else row.to_dict()
    author = first_author(data.get("authors"))
    year = clean_text(data.get("year")) or "unknown_year"
    title = clean_text(data.get("title")) or clean_text(data.get("doi")) or "untitled"
    title_slug = safe_filename(title.lower(), max_len=60)
    return safe_filename(f"{author}_{year}_{title_slug}", max_len=110)


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


def article_fulltext_preference_rank(url: Any) -> tuple[int, int, str]:
    text = clean_text(url)
    lowered = text.lower()
    if not text:
        return (9, 0, "")
    if lowered.endswith((".xml", ".nxml")) or "fulltextxml" in lowered:
        return (0, len(text), text)
    if any(token in lowered for token in ["/full", "fulltext", "/article/", "/articles/"]):
        if not is_likely_article_pdf_url(text):
            return (1, len(text), text)
    if lowered.endswith((".html", ".htm")):
        return (1, len(text), text)
    if "doi.org/" in lowered:
        return (3, len(text), text)
    if lowered.endswith(".pdf") or "/pdf" in lowered or "pdf=render" in lowered:
        return (4, len(text), text)
    return (2, len(text), text)


def article_file_preference_rank(path: Path) -> int:
    suffix = path.suffix.lower()
    if suffix in {".html", ".htm", ".xml", ".nxml"}:
        return 0
    if suffix == ".pdf":
        return 1
    return 2


def append_rows_csv(path: str | Path, rows: list[dict[str, Any]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    fieldnames = list(pd.DataFrame(rows).columns)
    exists = path.exists()
    with path.open("a", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        if not exists:
            writer.writeheader()
        writer.writerows(rows)
