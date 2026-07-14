from __future__ import annotations

import argparse
import csv
import hashlib
import re
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Any


FORMATS = ("html", "json", "text")
FORMAT_EXTENSIONS = {
    "html": "html",
    "json": "json",
    "text": "txt",
}


def short_stem(pdf_path: Path) -> str:
    digest = hashlib.sha1(pdf_path.name.encode("utf-8")).hexdigest()[:12]
    match = re.match(r"^(KWRAW-\d+)", pdf_path.stem)
    if match:
        prefix = match.group(1)
    else:
        prefix = re.sub(r"[^A-Za-z0-9]+", "_", pdf_path.stem).strip("_")[:32] or "pdf"
    return f"{prefix}__{digest}"


def output_paths(output_dir: Path, stem: str) -> dict[str, Path]:
    return {
        "html": output_dir / f"{stem}.html",
        "json": output_dir / f"{stem}.json",
        "text": output_dir / f"{stem}.txt",
    }


def append_manifest(path: Path, row: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    exists = path.exists()
    with path.open("a", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row.keys()))
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def converted(output_dir: Path, stem: str) -> bool:
    paths = output_paths(output_dir, stem)
    return all(path.exists() and path.stat().st_size > 0 for path in paths.values())


def run_command(command: list[str], timeout: int) -> subprocess.CompletedProcess[str]:
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        stdout, stderr = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        process.kill()
        stdout, stderr = process.communicate()
        raise subprocess.TimeoutExpired(command, timeout, output=stdout, stderr=stderr) from None
    return subprocess.CompletedProcess(command, process.returncode, stdout, stderr)


def convert_pdf(pdf_path: Path, output_dir: Path, stem: str, threads: int, timeout: int) -> tuple[str, str]:
    outputs = output_paths(output_dir, stem)
    if converted(output_dir, stem):
        return "success", "already_converted"

    with tempfile.TemporaryDirectory(prefix="mrd_odl_pdf_", ignore_cleanup_errors=True) as tmp:
        tmp_dir = Path(tmp)
        tmp_pdf = tmp_dir / f"{stem}.pdf"
        tmp_out = tmp_dir / "out"
        tmp_out.mkdir(parents=True, exist_ok=True)
        shutil.copy2(pdf_path, tmp_pdf)

        command = [
            "opendataloader-pdf",
            str(tmp_pdf),
            "-o",
            str(tmp_out),
            "-f",
            ",".join(FORMATS),
            "--threads",
            str(threads),
            "--image-output",
            "off",
            "--quiet",
        ]
        if shutil.which(command[0]) is None:
            return "failed", "opendataloader_pdf_cli_not_found"
        try:
            result = run_command(command, timeout)
        except subprocess.TimeoutExpired:
            return "failed", f"timeout_after_{timeout}_seconds"

        if result.returncode != 0:
            reason = (result.stderr or result.stdout or f"exit_{result.returncode}").strip()[:1000]
            return "failed", reason

        for fmt in FORMATS:
            src = tmp_out / f"{stem}.{FORMAT_EXTENSIONS[fmt]}"
            if not src.exists() or src.stat().st_size <= 0:
                return "failed", f"missing_{fmt}_output"
            shutil.copy2(src, outputs[fmt])

    return "success", ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--max-files", type=int, default=0)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--timeout", type=int, default=600)
    parser.add_argument("--progress-every", type=int, default=25)
    args = parser.parse_args()

    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    manifest = Path(args.manifest).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    pdfs = sorted(input_dir.glob("*.pdf"))
    if args.max_files > 0:
        pdfs = pdfs[: args.max_files]

    for index, pdf_path in enumerate(pdfs, start=1):
        stem = short_stem(pdf_path)
        outputs = output_paths(output_dir, stem)
        started = time.time()
        status, reason = convert_pdf(pdf_path, output_dir, stem, args.threads, args.timeout)
        elapsed = round(time.time() - started, 3)
        if index == 1 or index == len(pdfs) or index % args.progress_every == 0 or status != "success":
            print(f"[{index}/{len(pdfs)}] {status}: {pdf_path.name} ({elapsed}s)", flush=True)
        append_manifest(
            manifest,
            {
                "source_pdf": str(pdf_path),
                "source_filename": pdf_path.name,
                "source_dir": input_dir.name,
                "output_dir": str(output_dir),
                "output_stem": stem,
                "html_path": str(outputs["html"]) if outputs["html"].exists() else "",
                "json_path": str(outputs["json"]) if outputs["json"].exists() else "",
                "text_path": str(outputs["text"]) if outputs["text"].exists() else "",
                "status": status,
                "reason": reason,
                "elapsed_seconds": elapsed,
                "sequence_index": index,
            },
        )


if __name__ == "__main__":
    main()
