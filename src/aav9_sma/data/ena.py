"""Resolve and download public ENA FASTQ files with integrity checks."""

from __future__ import annotations

import hashlib
from collections.abc import Iterable
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import pandas as pd
import requests

ENA_FILE_REPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_FIELDS = "run_accession,fastq_ftp,fastq_md5,fastq_bytes"


def _https_url(ena_path: str) -> str:
    if ena_path.startswith("ftp://"):
        return "https://" + ena_path.removeprefix("ftp://")
    if ena_path.startswith("http://") or ena_path.startswith("https://"):
        return ena_path
    return "https://" + ena_path


def fetch_ena_fastq_manifest(accession: str = "PRJNA1131359") -> pd.DataFrame:
    """Fetch ENA FASTQ locations and checksums for a study or run accession."""
    response = requests.get(
        ENA_FILE_REPORT,
        params={
            "accession": accession,
            "result": "read_run",
            "fields": ENA_FIELDS,
            "format": "tsv",
        },
        timeout=120,
    )
    response.raise_for_status()
    rows = pd.read_csv(pd.io.common.StringIO(response.text), sep="\t", dtype=str)
    if rows.empty:
        raise ValueError(f"ENA returned no FASTQ records for {accession}")
    if rows["run_accession"].duplicated().any():
        raise ValueError("ENA returned duplicate run accessions")
    for column in ("fastq_ftp", "fastq_md5", "fastq_bytes"):
        if rows[column].str.contains(";", na=False).any():
            raise ValueError("This workflow currently expects one single-end FASTQ per run")
    rows = rows.rename(
        columns={
            "fastq_ftp": "fastq_url",
            "fastq_md5": "fastq_md5",
            "fastq_bytes": "fastq_bytes",
        }
    )
    rows["fastq_url"] = rows["fastq_url"].map(_https_url)
    rows["fastq_bytes"] = pd.to_numeric(rows["fastq_bytes"], errors="raise").astype("int64")
    return rows[["run_accession", "fastq_url", "fastq_md5", "fastq_bytes"]]


def resolve_ena_fastqs(target_manifest: pd.DataFrame, accession: str) -> pd.DataFrame:
    """Attach ENA FASTQ URLs, MD5 hashes, and compressed byte sizes to target runs."""
    if "run_accession" not in target_manifest:
        raise ValueError("Target manifest must contain run_accession")
    if target_manifest["run_accession"].duplicated().any():
        raise ValueError("Target manifest contains duplicate run accessions")
    ena = fetch_ena_fastq_manifest(accession)
    resolved = target_manifest.merge(ena, on="run_accession", how="left", validate="one_to_one")
    missing = resolved.loc[resolved["fastq_url"].isna(), "run_accession"].tolist()
    if missing:
        raise ValueError(f"Missing ENA FASTQ records for: {missing}")
    resolved["fastq_bytes"] = resolved["fastq_bytes"].astype("int64")
    return resolved


def _md5_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.md5()  # noqa: S324 - ENA publishes MD5 for transfer integrity.
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_fastq(
    row: dict[str, Any],
    output_dir: str | Path,
    *,
    chunk_size: int = 8 * 1024 * 1024,
) -> dict[str, object]:
    """Download one FASTQ with HTTP range resume, size check, and ENA MD5 check."""
    run_accession = str(row["run_accession"])
    expected_bytes = int(row["fastq_bytes"])
    expected_md5 = str(row["fastq_md5"]).lower()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir / f"{run_accession}.fastq.gz"
    partial = output_dir / f"{run_accession}.fastq.gz.part"

    if destination.exists():
        if destination.stat().st_size != expected_bytes:
            raise ValueError(f"Existing {destination} has the wrong byte size")
        observed_md5 = _md5_file(destination, chunk_size)
        if observed_md5 != expected_md5:
            raise ValueError(f"Existing {destination} failed its ENA MD5 check")
        return {
            "run_accession": run_accession,
            "status": "cached",
            "path": str(destination),
            "bytes": expected_bytes,
            "md5": observed_md5,
        }

    offset = partial.stat().st_size if partial.exists() else 0
    if offset > expected_bytes:
        raise ValueError(f"Partial {partial} is larger than the ENA file")
    if offset < expected_bytes:
        headers = {"Range": f"bytes={offset}-"} if offset else {}
        response = requests.get(
            str(row["fastq_url"]), headers=headers, stream=True, timeout=(30, 180)
        )
        response.raise_for_status()
        if offset and response.status_code != 206:
            offset = 0
        mode = "ab" if offset and response.status_code == 206 else "wb"
        with partial.open(mode) as stream:
            for chunk in response.iter_content(chunk_size=chunk_size):
                if chunk:
                    stream.write(chunk)

    observed_bytes = partial.stat().st_size
    if observed_bytes != expected_bytes:
        raise ValueError(
            f"Incomplete download for {run_accession}: {observed_bytes} != {expected_bytes} bytes"
        )
    observed_md5 = _md5_file(partial, chunk_size)
    if observed_md5 != expected_md5:
        raise ValueError(f"Downloaded {run_accession} failed its ENA MD5 check")
    partial.replace(destination)
    return {
        "run_accession": run_accession,
        "status": "downloaded",
        "path": str(destination),
        "bytes": observed_bytes,
        "md5": observed_md5,
    }


def download_fastq_manifest(
    rows: Iterable[dict[str, Any]], output_dir: str | Path, workers: int = 4
) -> list[dict[str, object]]:
    """Download a resolved manifest concurrently and return results in manifest order."""
    rows = list(rows)
    if workers < 1:
        raise ValueError("workers must be at least 1")
    results: dict[str, dict[str, object]] = {}
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(download_fastq, row, output_dir): str(row["run_accession"])
            for row in rows
        }
        for future in as_completed(futures):
            run_accession = futures[future]
            results[run_accession] = future.result()
    return [results[str(row["run_accession"])] for row in rows]
