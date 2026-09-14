"""Command-line entry points for data auditing and candidate ranking."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from aav9_sma.data.audit import audit_csv
from aav9_sma.data.fastq import count_fastq
from aav9_sma.data.fit4function import audit_official_release
from aav9_sma.data.sra import fetch_sra_manifest, summarize_sra_manifest
from aav9_sma.models.evaluate import (
    SCREEN_TASKS,
    benchmark_production_generalization,
    benchmark_screen_models,
)
from aav9_sma.screening.score import rank_candidates


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="aav9-sma")
    subparsers = parser.add_subparsers(dest="command", required=True)

    audit_parser = subparsers.add_parser("audit-data", help="Audit a canonical CSV file")
    audit_parser.add_argument("input", type=Path)
    audit_parser.add_argument("--output", type=Path)

    release_parser = subparsers.add_parser(
        "audit-fit4function", help="Audit an official Fit4Function checkout"
    )
    release_parser.add_argument("repository_root", type=Path)
    release_parser.add_argument("--output", type=Path)

    sra_parser = subparsers.add_parser(
        "fetch-sra-manifest", help="Fetch public SRA metadata for a BioProject"
    )
    sra_parser.add_argument("--bioproject", default="PRJNA1131359")
    sra_parser.add_argument("--output-csv", type=Path, required=True)
    sra_parser.add_argument("--output-summary", type=Path, required=True)

    fastq_parser = subparsers.add_parser(
        "count-fastq", help="Pilot extraction of Fit4Function insertions from FASTQ"
    )
    fastq_parser.add_argument("input", type=Path)
    fastq_parser.add_argument("--max-reads", type=int)
    fastq_parser.add_argument("--whitelist-csv", type=Path)
    fastq_parser.add_argument("--whitelist-column", default="AA")
    fastq_parser.add_argument("--counts-output", type=Path)
    fastq_parser.add_argument("--output", type=Path)

    benchmark_parser = subparsers.add_parser(
        "benchmark-fit4function", help="Benchmark models on the sequence-linked 100K screen"
    )
    benchmark_parser.add_argument("input", type=Path)
    benchmark_parser.add_argument(
        "--models", nargs="+", choices=("ridge", "random_forest"), default=["ridge"]
    )
    benchmark_parser.add_argument("--tasks", nargs="+", default=list(SCREEN_TASKS))
    benchmark_parser.add_argument("--output", type=Path, required=True)

    production_parser = subparsers.add_parser(
        "benchmark-production", help="Benchmark production on the independent assessment library"
    )
    production_parser.add_argument("modeling_csv", type=Path)
    production_parser.add_argument("assessment_csv", type=Path)
    production_parser.add_argument(
        "--models", nargs="+", choices=("ridge", "random_forest"), default=["ridge"]
    )
    production_parser.add_argument("--output", type=Path, required=True)

    rank_parser = subparsers.add_parser("rank-candidates", help="Rank model predictions")
    rank_parser.add_argument("input", type=Path)
    rank_parser.add_argument("--packaging-threshold", type=float, required=True)
    rank_parser.add_argument("--output", type=Path, required=True)
    return parser


def _write_json(payload: dict[str, object], output: Path | None) -> None:
    rendered = json.dumps(payload, indent=2, ensure_ascii=False)
    if output is None:
        print(rendered)
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(rendered + "\n", encoding="utf-8")


def main() -> None:
    args = _build_parser().parse_args()
    if args.command == "audit-data":
        _write_json(audit_csv(args.input).to_dict(), args.output)
        return
    if args.command == "audit-fit4function":
        _write_json(audit_official_release(args.repository_root), args.output)
        return
    if args.command == "fetch-sra-manifest":
        rows = fetch_sra_manifest(args.bioproject)
        args.output_csv.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(args.output_csv, index=False)
        _write_json(summarize_sra_manifest(rows), args.output_summary)
        return
    if args.command == "count-fastq":
        whitelist = None
        if args.whitelist_csv is not None:
            whitelist = set(pd.read_csv(args.whitelist_csv)[args.whitelist_column].dropna())
        _write_json(
            count_fastq(
                args.input,
                max_reads=args.max_reads,
                peptide_whitelist=whitelist,
                counts_output=args.counts_output,
            ),
            args.output,
        )
        return
    if args.command == "benchmark-fit4function":
        rows = benchmark_screen_models(
            args.input,
            model_names=tuple(args.models),
            tasks=tuple(args.tasks),
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(args.output, index=False)
        return
    if args.command == "benchmark-production":
        rows = benchmark_production_generalization(
            args.modeling_csv,
            args.assessment_csv,
            model_names=tuple(args.models),
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(args.output, index=False)
        return
    if args.command == "rank-candidates":
        predictions = pd.read_csv(args.input)
        ranked = rank_candidates(predictions, packaging_threshold=args.packaging_threshold)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        ranked.to_csv(args.output, index=False)
        return
    raise RuntimeError(f"Unhandled command: {args.command}")


if __name__ == "__main__":
    main()
