#!/usr/bin/env python3
"""Stream public NCBI SRA runs through FASTQ conversion and Bowtie2 counting.

This is a recovery/acceleration path for environments where ENA FASTQ downloads
are throttled. Each SRA archive is validated before conversion, and temporary
SRA/FASTQ files are removed only after the count table and QC summary are written.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from aav9_sma.data.alignment import count_bowtie2_fastq

S3_HOST = "sra-pub-run-odp.s3.amazonaws.com"


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--roles", nargs="+", required=True)
    parser.add_argument("--runs", nargs="+")
    parser.add_argument("--sra-bin", type=Path, required=True)
    parser.add_argument("--s3-ip", required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--index-prefix", type=Path, required=True)
    parser.add_argument("--whitelist-csv", type=Path, required=True)
    parser.add_argument("--whitelist-column", default="AA")
    parser.add_argument("--counts-dir", type=Path, required=True)
    parser.add_argument("--summaries-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--threads-per-worker", type=int, default=2)
    parser.add_argument("--connections-per-download", type=int, default=4)
    parser.add_argument("--keep-temporary", action="store_true")
    return parser.parse_args()


def _run_checked(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, check=True, capture_output=True, text=True)


def _process_run(
    row: dict[str, object],
    *,
    sra_bin: Path,
    s3_ip: str,
    work_dir: Path,
    index_prefix: Path,
    whitelist: set[str],
    counts_dir: Path,
    summaries_dir: Path,
    threads: int,
    connections: int,
    keep_temporary: bool,
) -> dict[str, object]:
    run = str(row["run_accession"])
    counts_path = counts_dir / f"{run}.counts.csv.gz"
    summary_path = summaries_dir / f"{run}.json"
    if counts_path.exists() and summary_path.exists():
        return {"run_accession": run, "status": "cached"}

    run_dir = work_dir / run
    run_dir.mkdir(parents=True, exist_ok=True)
    sra_path = run_dir / f"{run}.sra"
    fastq_path = run_dir / f"{run}.fastq"
    url = f"https://{s3_ip}/sra/{run}/{run}"
    _run_checked(
        [
            "aria2c",
            url,
            f"--header=Host: {S3_HOST}",
            "--check-certificate=false",
            f"--dir={run_dir}",
            f"--out={sra_path.name}",
            "--file-allocation=none",
            "--continue=true",
            "--allow-overwrite=true",
            "--auto-file-renaming=false",
            f"--max-connection-per-server={connections}",
            f"--split={connections}",
            "--min-split-size=1M",
            "--max-tries=20",
            "--retry-wait=3",
            "--summary-interval=0",
            "--console-log-level=warn",
        ]
    )
    validation = _run_checked([str(sra_bin / "vdb-validate"), str(sra_path)])
    _run_checked(
        [
            str(sra_bin / "fasterq-dump"),
            str(sra_path),
            "--threads",
            str(threads),
            "--outdir",
            str(run_dir),
            "--temp",
            str(run_dir),
            "--force",
        ]
    )
    result = count_bowtie2_fastq(
        fastq_path,
        index_prefix,
        peptide_whitelist=whitelist,
        counts_output=counts_path,
        threads=threads,
    )
    result.update(
        {
            "path": f"ncbi-sra://{run}",
            "raw_source": "NCBI SRA Open Data (AWS S3)",
            "raw_source_url": f"https://{S3_HOST}/sra/{run}/{run}",
            "sra_validation": "consistent",
            "sra_validation_log": validation.stderr.strip(),
        }
    )
    summary_path.write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    if not keep_temporary:
        fastq_path.unlink(missing_ok=True)
        sra_path.unlink(missing_ok=True)
        (run_dir / f"{run}.sra.aria2").unlink(missing_ok=True)
        try:
            run_dir.rmdir()
        except OSError:
            pass
    return {
        "run_accession": run,
        "status": "processed",
        "reads_examined": result["reads_examined"],
        "whitelist_matched_reads": result["whitelist_matched_reads"],
    }


def main() -> None:
    args = _parse_args()
    if args.workers < 1 or args.threads_per_worker < 1:
        raise ValueError("workers and threads-per-worker must be positive")
    manifest = pd.read_csv(args.manifest)
    manifest = manifest.loc[manifest["project_role"].isin(args.roles)].copy()
    if args.runs:
        requested = set(args.runs)
        manifest = manifest.loc[manifest["run_accession"].isin(requested)].copy()
        missing = requested.difference(manifest["run_accession"])
        if missing:
            raise ValueError(f"Requested runs not found for selected roles: {sorted(missing)}")
    whitelist = set(
        pd.read_csv(args.whitelist_csv)[args.whitelist_column].dropna().astype(str)
    )
    args.counts_dir.mkdir(parents=True, exist_ok=True)
    args.summaries_dir.mkdir(parents=True, exist_ok=True)
    args.work_dir.mkdir(parents=True, exist_ok=True)

    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(
                _process_run,
                row,
                sra_bin=args.sra_bin,
                s3_ip=args.s3_ip,
                work_dir=args.work_dir,
                index_prefix=args.index_prefix,
                whitelist=whitelist,
                counts_dir=args.counts_dir,
                summaries_dir=args.summaries_dir,
                threads=args.threads_per_worker,
                connections=args.connections_per_download,
                keep_temporary=args.keep_temporary,
            ): str(row["run_accession"])
            for row in manifest.to_dict(orient="records")
        }
        for future in as_completed(futures):
            result = future.result()
            print(json.dumps(result, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
