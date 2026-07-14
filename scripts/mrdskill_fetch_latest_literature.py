#!/usr/bin/env python3
"""Fetch a small latest OA literature batch for MRDskill ingest testing.

The fetcher uses Europe PMC search to discover recent open-access MRD/ctDNA
papers, downloads PMC full-text XML where available, and attempts to mirror
open supplementary links listed in the article XML.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import html
import json
import re
import tarfile
import urllib.parse
import urllib.request
from datetime import datetime
from io import BytesIO
from pathlib import Path
from typing import Any
from xml.etree import ElementTree as ET

from mrdskill_project_paths import artifact_path, safe_path_component


EUROPE_PMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
NCBI_OA_XML = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fetch latest OA MRD/ctDNA literature into a kg_ready-style batch.")
    parser.add_argument("--limit", type=int, default=10)
    parser.add_argument("--query", default='((ctDNA OR "circulating tumor DNA" OR "circulating tumour DNA") AND ("minimal residual disease" OR "molecular residual disease" OR "measurable residual disease" OR MRD) AND OPEN_ACCESS:y AND HAS_FT:y)')
    parser.add_argument("--out-dir", default="", help="Defaults to imports/literature_search/latest_mrd_ctdna_<timestamp>.")
    parser.add_argument("--prefer-supplement", action="store_true", help="Fetch HAS_SUPPL:y papers first, then fill remaining recent papers.")
    parser.add_argument("--timeout", type=int, default=40)
    parser.add_argument("--user-agent", default="MRDskill-KGRAT/1.0 (research QC; contact: local)")
    return parser.parse_args()


def clean(value: Any) -> str:
    return "" if value is None else str(value).strip()


def request_json(url: str, timeout: int, user_agent: str) -> dict[str, Any]:
    req = urllib.request.Request(url, headers={"User-Agent": user_agent})
    with urllib.request.urlopen(req, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8", errors="replace"))


def request_text(url: str, timeout: int, user_agent: str) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": user_agent})
    with urllib.request.urlopen(req, timeout=timeout) as response:
        data = response.read()
        encoding = response.headers.get_content_charset() or "utf-8"
        return data.decode(encoding, errors="replace")


def request_bytes(url: str, timeout: int, user_agent: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": user_agent})
    with urllib.request.urlopen(req, timeout=timeout) as response:
        return response.read()


def europe_pmc_search(query: str, limit: int, timeout: int, user_agent: str) -> list[dict[str, Any]]:
    url = EUROPE_PMC_SEARCH + "?" + urllib.parse.urlencode(
        {
            "query": query,
            "format": "json",
            "pageSize": str(max(limit, 1)),
            "resultType": "core",
            "sort": "FIRST_PDATE_D desc",
        }
    )
    data = request_json(url, timeout, user_agent)
    return list(data.get("resultList", {}).get("result", []))


def dedupe_results(results: list[dict[str, Any]]) -> list[dict[str, Any]]:
    seen: set[str] = set()
    unique: list[dict[str, Any]] = []
    for row in results:
        key = clean(row.get("pmcid")) or clean(row.get("doi")) or clean(row.get("pmid")) or clean(row.get("title"))
        if not key or key in seen:
            continue
        seen.add(key)
        unique.append(row)
    return unique


def pmc_oa_xml_url(pmcid: str, timeout: int, user_agent: str) -> str:
    url = NCBI_OA_XML + "?" + urllib.parse.urlencode({"id": pmcid})
    text = request_text(url, timeout, user_agent)
    root = ET.fromstring(text)
    for link in root.findall(".//link"):
        href = link.attrib.get("href", "")
        fmt = link.attrib.get("format", "")
        if fmt.lower() == "xml" and href:
            return href
    return ""


def pmc_oa_package_url(pmcid: str, timeout: int, user_agent: str) -> str:
    url = NCBI_OA_XML + "?" + urllib.parse.urlencode({"id": pmcid})
    text = request_text(url, timeout, user_agent)
    root = ET.fromstring(text)
    for link in root.findall(".//link"):
        href = link.attrib.get("href", "")
        fmt = link.attrib.get("format", "")
        if fmt.lower() in {"tgz", "tar.gz"} and href:
            return href
    return ""


def europe_pmc_fulltext_xml_url(pmcid: str) -> str:
    return f"https://www.ebi.ac.uk/europepmc/webservices/rest/{urllib.parse.quote(pmcid)}/fullTextXML"


def fetch_fulltext_xml(pmcid: str, timeout: int, user_agent: str) -> tuple[str, str, str]:
    """Return XML text, source URL, and status."""
    url = europe_pmc_fulltext_xml_url(pmcid)
    try:
        xml_text = request_text(url, timeout, user_agent)
        if "<article" in xml_text[:1000] or "<pmc-articleset" in xml_text[:1000]:
            return xml_text, url, "downloaded"
    except Exception:
        pass
    try:
        oa_url = pmc_oa_xml_url(pmcid, timeout, user_agent)
        if oa_url and oa_url.endswith(".xml"):
            xml_text = request_text(oa_url, timeout, user_agent)
            return xml_text, oa_url, "downloaded"
        if oa_url:
            return "", oa_url, "oa_package_not_xml"
    except Exception as exc:
        return "", "", f"failed:{type(exc).__name__}"
    return "", "", "oa_xml_not_found"


def strip_xml_text(xml_text: str) -> str:
    try:
        root = ET.fromstring(xml_text.encode("utf-8"))
    except ET.ParseError:
        return re.sub(r"<[^>]+>", " ", xml_text)
    text = " ".join(part.strip() for part in root.itertext() if part and part.strip())
    return re.sub(r"\s+", " ", html.unescape(text)).strip()


def suffix_from_url(url: str, fallback: str = ".dat") -> str:
    path = urllib.parse.urlparse(url).path
    suffix = Path(path).suffix
    if suffix:
        return suffix[:12]
    return fallback


def supplementary_links_from_xml(xml_text: str) -> list[dict[str, str]]:
    links: list[dict[str, str]] = []
    try:
        root = ET.fromstring(xml_text.encode("utf-8"))
    except ET.ParseError:
        return links
    for elem in root.iter():
        tag = elem.tag.split("}")[-1]
        if tag not in {"supplementary-material", "media"}:
            continue
        href = ""
        for key, value in elem.attrib.items():
            if key.endswith("href"):
                href = value
                break
        if not href:
            continue
        label = clean(elem.attrib.get("id")) or clean(elem.attrib.get("content-type")) or "supplement"
        links.append({"label": safe_path_component(label, "supplement"), "href": href})
    return links


def safe_archive_member_name(name: str) -> str:
    path = Path(name)
    safe_parts = [safe_path_component(part, "part") for part in path.parts if part not in {"", ".", ".."}]
    return "_".join(safe_parts) or "member"


def extract_oa_package(
    pmcid: str,
    record_id: str,
    package_url: str,
    fulltext_dir: Path,
    converted_dir: Path,
    supplement_dir: Path,
    timeout: int,
    user_agent: str,
) -> dict[str, Any]:
    package_url_http = package_url.replace("ftp://", "https://", 1)
    data = request_bytes(package_url_http, timeout, user_agent)
    xml_text = ""
    xml_path = ""
    text_path = ""
    supplement_rows: list[dict[str, str]] = []
    with tarfile.open(fileobj=BytesIO(data), mode="r:gz") as archive:
        for member in archive.getmembers():
            if not member.isfile():
                continue
            member_name = safe_archive_member_name(member.name)
            suffix = Path(member.name).suffix.lower()
            extracted = archive.extractfile(member)
            if extracted is None:
                continue
            payload = extracted.read()
            if suffix in {".nxml", ".xml"} and not xml_text:
                xml_text = payload.decode("utf-8", errors="replace")
                xml_file = fulltext_dir / f"{record_id}_{pmcid}.xml"
                txt_file = converted_dir / f"{record_id}_{pmcid}.txt"
                xml_file.write_text(xml_text, encoding="utf-8")
                txt_file.write_text(strip_xml_text(xml_text), encoding="utf-8")
                xml_path = str(xml_file)
                text_path = str(txt_file)
                continue
            if suffix in {".pdf", ".xlsx", ".xls", ".csv", ".tsv", ".txt", ".docx", ".doc", ".zip"}:
                target = supplement_dir / f"{record_id}_{pmcid}_{member_name}"
                target.write_bytes(payload)
                supplement_rows.append(
                    {
                        "record_id": record_id,
                        "pmcid": pmcid,
                        "source_url": package_url,
                        "source_filename": member.name,
                        "local_path": str(target),
                        "status": "downloaded_from_oa_package",
                    }
                )
    return {
        "xml_text": xml_text,
        "xml_path": xml_path,
        "text_path": text_path,
        "supplement_rows": supplement_rows,
        "package_url": package_url,
        "package_bytes": len(data),
    }


def absolute_pmc_url(pmcid: str, href: str) -> str:
    if href.startswith("http://") or href.startswith("https://"):
        return href
    pmc_number = pmcid[3:] if pmcid.upper().startswith("PMC") else pmcid
    return f"https://pmc.ncbi.nlm.nih.gov/articles/instance/{pmc_number}/bin/{href.lstrip('/')}"


def looks_like_blocked_download(data: bytes, content_type: str, suffix: str) -> bool:
    head = data[:500].lower()
    if b"preparing to download" in head or b"pow_challenge" in head:
        return True
    if "text/html" in content_type.lower() and suffix.lower() not in {".html", ".htm"}:
        return True
    return False


def validate_supplement_payload(data: bytes, content_type: str, suffix: str) -> None:
    suffix = suffix.lower()
    if looks_like_blocked_download(data, content_type, suffix):
        raise ValueError("blocked_or_html_download_page")
    if suffix == ".pdf" and not data.startswith(b"%PDF"):
        raise ValueError("pdf_magic_missing")
    if suffix in {".xlsx", ".docx", ".pptx", ".zip"} and not data.startswith(b"PK"):
        raise ValueError("zip_magic_missing")


def write_csv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = artifact_path(args.out_dir) if args.out_dir else artifact_path(f"imports/literature_search/latest_mrd_ctdna_{timestamp}")
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_dir = out_dir / "manifests"
    fulltext_dir = out_dir / "downloaded_pdfs" / "publisher_fulltext_html_xml"
    converted_dir = out_dir / "downloaded_pdfs" / "converted_from_pdf"
    supplement_dir = out_dir / "supplementary_materials" / "native_non_pdf"
    for directory in [manifest_dir, fulltext_dir, converted_dir, supplement_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    queries = []
    if args.prefer_supplement:
        queries.append(f"({args.query}) AND HAS_SUPPL:y")
    queries.append(args.query)

    results: list[dict[str, Any]] = []
    for query in queries:
        results.extend(europe_pmc_search(query, args.limit * 3, args.timeout, args.user_agent))
    selected = [row for row in dedupe_results(results) if clean(row.get("pmcid"))][: args.limit]

    candidate_rows: list[dict[str, str]] = []
    converted_manifest: list[dict[str, str]] = []
    fetch_rows: list[dict[str, str]] = []
    supplement_rows: list[dict[str, str]] = []

    for idx, result in enumerate(selected, start=1):
        record_id = f"KWRAW-INCR-{idx:04d}"
        pmcid = clean(result.get("pmcid"))
        title = clean(result.get("title"))
        year = clean(result.get("pubYear")) or clean(result.get("firstPublicationDate"))[:4]
        journal = clean(result.get("journalTitle"))
        doi = clean(result.get("doi"))
        pmid = clean(result.get("pmid"))
        xml_status = "not_attempted"
        xml_url = ""
        xml_path = ""
        text_path = ""
        supplement_count = 0
        error = ""
        try:
            package_url = pmc_oa_package_url(pmcid, args.timeout, args.user_agent)
            package_result: dict[str, Any] = {}
            if package_url:
                try:
                    package_result = extract_oa_package(
                        pmcid,
                        record_id,
                        package_url,
                        fulltext_dir,
                        converted_dir,
                        supplement_dir,
                        args.timeout,
                        args.user_agent,
                    )
                    xml_text = clean(package_result.get("xml_text"))
                    xml_url = clean(package_result.get("package_url"))
                    xml_path = clean(package_result.get("xml_path"))
                    text_path = clean(package_result.get("text_path"))
                    supplement_rows.extend(package_result.get("supplement_rows", []))
                    supplement_count = len(package_result.get("supplement_rows", []))
                    xml_status = "downloaded_oa_package" if xml_text else "oa_package_no_xml"
                except Exception as exc:
                    xml_status = f"oa_package_failed:{type(exc).__name__}"

            if xml_status != "downloaded_oa_package":
                xml_text, xml_url, xml_status = fetch_fulltext_xml(pmcid, args.timeout, args.user_agent)
            if xml_status == "downloaded":
                xml_name = f"{record_id}_{pmcid}.xml"
                txt_name = f"{record_id}_{pmcid}.txt"
                xml_file = fulltext_dir / xml_name
                txt_file = converted_dir / txt_name
                xml_file.write_text(xml_text, encoding="utf-8")
                txt_file.write_text(strip_xml_text(xml_text), encoding="utf-8")
                xml_path = str(xml_file)
                text_path = str(txt_file)
                for link_idx, link in enumerate(supplementary_links_from_xml(xml_text), start=1):
                    href = absolute_pmc_url(pmcid, link["href"])
                    suffix = suffix_from_url(href)
                    supp_name = f"{record_id}_{pmcid}_{link_idx:02d}_{link['label']}{suffix}"
                    supp_path = supplement_dir / supp_name
                    try:
                        req = urllib.request.Request(href, headers={"User-Agent": args.user_agent})
                        with urllib.request.urlopen(req, timeout=args.timeout) as response:
                            data = response.read()
                            content_type = response.headers.get("content-type", "")
                        validate_supplement_payload(data, content_type, suffix)
                        supp_path.write_bytes(data)
                        supplement_count += 1
                        supplement_rows.append(
                            {
                                "record_id": record_id,
                                "pmcid": pmcid,
                                "source_url": href,
                                "source_filename": supp_name,
                                "local_path": str(supp_path),
                                "status": "downloaded",
                            }
                        )
                    except Exception as exc:
                        supplement_rows.append(
                            {
                                "record_id": record_id,
                                "pmcid": pmcid,
                                "source_url": href,
                                "source_filename": supp_name,
                                "local_path": str(supp_path),
                                "status": f"failed:{type(exc).__name__}",
                            }
                        )
            elif xml_status == "downloaded_oa_package":
                pass
        except Exception as exc:
            xml_status = f"failed:{type(exc).__name__}"
            error = str(exc)[:300]

        candidate_rows.append(
            {
                "record_id": record_id,
                "title": title,
                "year": year,
                "journal": journal,
                "doi": doi,
                "pmid": pmid,
                "pmcid": pmcid,
                "source_filename": f"{record_id}_{pmcid}.xml",
                "has_supplement": clean(result.get("hasSuppl")),
                "first_publication_date": clean(result.get("firstPublicationDate")),
                "source": "Europe PMC latest OA search",
            }
        )
        converted_manifest.append(
            {
                "source_filename": f"{record_id}_{pmcid}.xml",
                "output_stem": f"{record_id}_{pmcid}",
                "title": title,
                "year": year,
                "journal": journal,
                "doi": doi,
                "pmid": pmid,
                "pmcid": pmcid,
                "conversion_status": xml_status,
            }
        )
        fetch_rows.append(
            {
                "record_id": record_id,
                "pmcid": pmcid,
                "title": title,
                "journal": journal,
                "doi": doi,
                "pmid": pmid,
                "xml_url": xml_url,
                "xml_path": xml_path,
                "text_path": text_path,
                "supplement_count": str(supplement_count),
                "status": xml_status,
                "error": error,
            }
        )
        print(f"[{idx}/{len(selected)}] {record_id} {pmcid} {xml_status} supplements={supplement_count} {title[:80]}")

    write_csv(
        manifest_dir / "keyword_research_candidate_table.csv",
        candidate_rows,
        ["record_id", "title", "year", "journal", "doi", "pmid", "pmcid", "source_filename", "has_supplement", "first_publication_date", "source"],
    )
    write_csv(
        manifest_dir / "downloaded_pdfs_converted_manifest.csv",
        converted_manifest,
        ["source_filename", "output_stem", "title", "year", "journal", "doi", "pmid", "pmcid", "conversion_status"],
    )
    write_csv(
        manifest_dir / "supplementary_materials_manifest.csv",
        supplement_rows,
        ["record_id", "pmcid", "source_url", "source_filename", "local_path", "status"],
    )
    report = {
        "query": args.query,
        "limit": args.limit,
        "out_dir": str(out_dir),
        "selected_count": len(selected),
        "fulltext_downloaded": sum(1 for row in fetch_rows if row["status"] == "downloaded"),
        "supplements_downloaded": sum(int(row["supplement_count"]) for row in fetch_rows),
        "records": fetch_rows,
    }
    (out_dir / "fetch_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"out_dir": str(out_dir), "selected_count": len(selected), "fulltext_downloaded": report["fulltext_downloaded"], "supplements_downloaded": report["supplements_downloaded"]}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
