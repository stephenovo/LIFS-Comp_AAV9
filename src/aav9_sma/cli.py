"""Command-line entry points for data auditing and candidate ranking."""

from __future__ import annotations

import argparse
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from aav9_sma.blind import freeze_blind_study, verify_blind_freeze
from aav9_sma.data.alignment import count_bowtie2_fastq, write_bowtie2_reference
from aav9_sma.data.audit import audit_csv
from aav9_sma.data.ena import download_fastq_manifest, resolve_ena_fastqs
from aav9_sma.data.fastq import count_fastq
from aav9_sma.data.fit4function import audit_official_release
from aav9_sma.data.reconstruct import (
    load_count_matrix,
    reconstruct_liver,
    reconstruct_multiorgan,
)
from aav9_sma.data.sra import fetch_sra_manifest, summarize_sra_manifest
from aav9_sma.demo import run_demo
from aav9_sma.models.evaluate import (
    MULTIORGAN_ENDPOINTS,
    SCREEN_TASKS,
    benchmark_masked_multitask_animal_holdout,
    benchmark_multiorgan_animal_holdout,
    benchmark_multitask_animal_holdout,
    benchmark_multitask_ensemble_animal_holdout,
    benchmark_multitask_ensemble_leave_one_animal_out,
    benchmark_production_generalization,
    benchmark_screen_models,
)
from aav9_sma.repro import verify_manifest
from aav9_sma.screening.audit import (
    build_composition_audit,
    build_gate_sensitivity_audit,
    select_composition_challenge_panel,
    summarize_funnel_audit,
)
from aav9_sma.screening.score import rank_candidates
from aav9_sma.screening.virtual import (
    fit_calibrated_packaging_model,
    run_control_audit,
    run_virtual_screen,
)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="aav9-sma")
    subparsers = parser.add_subparsers(dest="command", required=True)

    demo_parser = subparsers.add_parser(
        "demo",
        help="Run a self-contained synthetic smoke test without external research data",
    )
    demo_parser.add_argument("--output-dir", type=Path, default=Path("artifacts/demo"))
    demo_parser.add_argument("--random-state", type=int, default=42)

    manifest_parser = subparsers.add_parser(
        "verify-manifest", help="Verify SHA256 checksums in a reproducibility manifest"
    )
    manifest_parser.add_argument("manifest", type=Path)
    manifest_parser.add_argument("--root", type=Path, default=Path("."))

    blind_freeze_parser = subparsers.add_parser(
        "freeze-blind-study",
        help="Freeze a one-shot external blind-study plan before outcome access",
    )
    blind_freeze_parser.add_argument("config", type=Path)
    blind_freeze_parser.add_argument("--root", type=Path, default=Path("."))
    blind_freeze_parser.add_argument("--output", type=Path, required=True)

    blind_verify_parser = subparsers.add_parser(
        "verify-blind-freeze",
        help="Verify the frozen files and Git state before locked blind analysis",
    )
    blind_verify_parser.add_argument("manifest", type=Path)
    blind_verify_parser.add_argument("--root", type=Path, default=Path("."))
    blind_verify_parser.add_argument("--output", type=Path)

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

    ena_parser = subparsers.add_parser(
        "resolve-ena-fastq", help="Attach ENA FASTQ URLs and checksums to an SRA manifest"
    )
    ena_parser.add_argument("input", type=Path)
    ena_parser.add_argument("--accession", default="PRJNA1131359")
    ena_parser.add_argument("--output", type=Path, required=True)

    download_parser = subparsers.add_parser(
        "download-ena-fastq", help="Download resolved ENA FASTQs with resume and MD5 checks"
    )
    download_parser.add_argument("input", type=Path)
    download_parser.add_argument("--output-dir", type=Path, required=True)
    download_parser.add_argument("--workers", type=int, default=4)
    download_parser.add_argument("--project-roles", nargs="+")
    download_parser.add_argument("--exclude-alias-regex")
    download_parser.add_argument("--output", type=Path)

    fastq_parser = subparsers.add_parser(
        "count-fastq", help="Pilot extraction of Fit4Function insertions from FASTQ"
    )
    fastq_parser.add_argument("input", type=Path)
    fastq_parser.add_argument("--max-reads", type=int)
    fastq_parser.add_argument("--whitelist-csv", type=Path)
    fastq_parser.add_argument("--whitelist-column", default="AA")
    fastq_parser.add_argument("--counts-output", type=Path)
    fastq_parser.add_argument("--output", type=Path)

    reference_parser = subparsers.add_parser(
        "prepare-bowtie2-reference", help="Build the published short-reference Bowtie2 index"
    )
    reference_parser.add_argument("--output-prefix", type=Path, required=True)
    reference_parser.add_argument("--output", type=Path)

    bowtie_parser = subparsers.add_parser(
        "count-bowtie2", help="Count Q20 7-mers after the published Bowtie2 alignment"
    )
    bowtie_parser.add_argument("input", type=Path)
    bowtie_parser.add_argument("--index-prefix", type=Path, required=True)
    bowtie_parser.add_argument("--whitelist-csv", type=Path, required=True)
    bowtie_parser.add_argument("--whitelist-column", default="AA")
    bowtie_parser.add_argument("--counts-output", type=Path, required=True)
    bowtie_parser.add_argument("--threads", type=int, default=2)
    bowtie_parser.add_argument("--output", type=Path)

    batch_parser = subparsers.add_parser(
        "count-fastq-batch", help="Count whitelist 7-mers across downloaded ENA FASTQs"
    )
    batch_parser.add_argument("manifest", type=Path)
    batch_parser.add_argument("--fastq-dir", type=Path, required=True)
    batch_parser.add_argument("--whitelist-csv", type=Path, required=True)
    batch_parser.add_argument("--whitelist-column", default="AA")
    batch_parser.add_argument("--counts-dir", type=Path, required=True)
    batch_parser.add_argument("--summaries-dir", type=Path, required=True)
    batch_parser.add_argument("--workers", type=int, default=4)
    batch_parser.add_argument("--project-roles", nargs="+")
    batch_parser.add_argument("--exclude-alias-regex")
    batch_parser.add_argument("--overwrite", action="store_true")

    bowtie_batch_parser = subparsers.add_parser(
        "count-bowtie2-batch", help="Run published Bowtie2 counting across an ENA manifest"
    )
    bowtie_batch_parser.add_argument("manifest", type=Path)
    bowtie_batch_parser.add_argument("--fastq-dir", type=Path, required=True)
    bowtie_batch_parser.add_argument("--index-prefix", type=Path, required=True)
    bowtie_batch_parser.add_argument("--whitelist-csv", type=Path, required=True)
    bowtie_batch_parser.add_argument("--whitelist-column", default="AA")
    bowtie_batch_parser.add_argument("--counts-dir", type=Path, required=True)
    bowtie_batch_parser.add_argument("--summaries-dir", type=Path, required=True)
    bowtie_batch_parser.add_argument("--workers", type=int, default=2)
    bowtie_batch_parser.add_argument("--threads-per-worker", type=int, default=2)
    bowtie_batch_parser.add_argument("--project-roles", nargs="+")
    bowtie_batch_parser.add_argument("--exclude-alias-regex")
    bowtie_batch_parser.add_argument("--overwrite", action="store_true")

    reconstruct_parser = subparsers.add_parser(
        "reconstruct-liver", help="Reconstruct liver enrichments and validate public labels"
    )
    reconstruct_parser.add_argument("manifest", type=Path)
    reconstruct_parser.add_argument("--counts-dir", type=Path, required=True)
    reconstruct_parser.add_argument("--summaries-dir", type=Path, required=True)
    reconstruct_parser.add_argument("--public-screens", type=Path, required=True)
    reconstruct_parser.add_argument("--output-reconstruction", type=Path, required=True)
    reconstruct_parser.add_argument("--output-metrics", type=Path, required=True)
    reconstruct_parser.add_argument("--output-qc", type=Path, required=True)
    reconstruct_parser.add_argument("--exclude-alias-regex")

    multiorgan_parser = subparsers.add_parser(
        "reconstruct-multiorgan",
        help="Reconstruct animal-level multi-organ enrichments from raw counts",
    )
    multiorgan_parser.add_argument("manifest", type=Path)
    multiorgan_parser.add_argument("--counts-dir", type=Path, required=True)
    multiorgan_parser.add_argument("--summaries-dir", type=Path, required=True)
    multiorgan_parser.add_argument("--virus-round", type=int, default=2)
    multiorgan_parser.add_argument("--output-reconstruction", type=Path, required=True)
    multiorgan_parser.add_argument("--output-metrics", type=Path, required=True)
    multiorgan_parser.add_argument("--output-qc", type=Path, required=True)

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

    multiorgan_benchmark_parser = subparsers.add_parser(
        "benchmark-multiorgan",
        help="Evaluate unseen sequences against held-out animal 4 labels",
    )
    multiorgan_benchmark_parser.add_argument("input", type=Path)
    multiorgan_benchmark_parser.add_argument(
        "--models",
        nargs="+",
        choices=("ridge", "random_forest", "lightgbm"),
        default=["ridge"],
    )
    multiorgan_benchmark_parser.add_argument(
        "--endpoints", nargs="+", default=list(MULTIORGAN_ENDPOINTS)
    )
    multiorgan_benchmark_parser.add_argument(
        "--feature-set",
        choices=("one_hot", "one_hot_physchem"),
        default="one_hot",
    )
    multiorgan_benchmark_parser.add_argument("--bootstrap-resamples", type=int, default=0)
    multiorgan_benchmark_parser.add_argument("--output", type=Path, required=True)

    multitask_benchmark_parser = subparsers.add_parser(
        "benchmark-multitask",
        help="Evaluate a shared multi-output MLP against held-out animal 4 labels",
    )
    multitask_benchmark_parser.add_argument("input", type=Path)
    multitask_benchmark_parser.add_argument(
        "--endpoints", nargs="+", default=list(MULTIORGAN_ENDPOINTS)
    )
    multitask_benchmark_parser.add_argument("--max-iter", type=int, default=80)
    multitask_benchmark_parser.add_argument("--output", type=Path, required=True)

    ensemble_benchmark_parser = subparsers.add_parser(
        "benchmark-multitask-ensemble",
        help="Evaluate the shared MLP ensemble used by virtual screening",
    )
    ensemble_benchmark_parser.add_argument("input", type=Path)
    ensemble_benchmark_parser.add_argument(
        "--endpoints", nargs="+", default=list(MULTIORGAN_ENDPOINTS)
    )
    ensemble_benchmark_parser.add_argument("--ensemble-size", type=int, default=5)
    ensemble_benchmark_parser.add_argument("--max-iter", type=int, default=80)
    ensemble_benchmark_parser.add_argument("--bootstrap-resamples", type=int, default=0)
    ensemble_benchmark_parser.add_argument("--output", type=Path, required=True)

    cross_animal_parser = subparsers.add_parser(
        "benchmark-cross-animal",
        help="Retrospectively leave out each animal with the frozen shared-MLP ensemble",
    )
    cross_animal_parser.add_argument("input", type=Path)
    cross_animal_parser.add_argument(
        "--endpoints", nargs="+", default=list(MULTIORGAN_ENDPOINTS)
    )
    cross_animal_parser.add_argument("--animals", nargs="+", type=int, default=[1, 2, 3, 4])
    cross_animal_parser.add_argument("--ensemble-size", type=int, default=5)
    cross_animal_parser.add_argument("--max-iter", type=int, default=80)
    cross_animal_parser.add_argument("--bootstrap-resamples", type=int, default=500)
    cross_animal_parser.add_argument("--output", type=Path, required=True)

    masked_benchmark_parser = subparsers.add_parser(
        "benchmark-masked-multitask",
        help="Evaluate the optional PyTorch multi-task MLP with missing-label masks",
    )
    masked_benchmark_parser.add_argument("input", type=Path)
    masked_benchmark_parser.add_argument(
        "--endpoints", nargs="+", default=list(MULTIORGAN_ENDPOINTS)
    )
    masked_benchmark_parser.add_argument("--max-epochs", type=int, default=80)
    masked_benchmark_parser.add_argument("--bootstrap-resamples", type=int, default=500)
    masked_benchmark_parser.add_argument("--device", choices=("cpu", "mps"), default="cpu")
    masked_benchmark_parser.add_argument("--output", type=Path, required=True)

    rank_parser = subparsers.add_parser("rank-candidates", help="Rank model predictions")
    rank_parser.add_argument("input", type=Path)
    rank_parser.add_argument("--packaging-threshold", type=float, required=True)
    rank_parser.add_argument("--output", type=Path, required=True)

    virtual_parser = subparsers.add_parser(
        "screen-virtual",
        help="Generate, predict, rank, and diversify virtual 7-mer candidates",
    )
    virtual_parser.add_argument("screen_csv", type=Path)
    virtual_parser.add_argument("reconstructed_csv", type=Path)
    virtual_parser.add_argument("--pool-size", type=int, default=200_000)
    virtual_parser.add_argument("--ensemble-size", type=int, default=5)
    virtual_parser.add_argument("--max-iter", type=int, default=80)
    virtual_parser.add_argument("--output-ranked", type=Path, required=True)
    virtual_parser.add_argument("--output-pareto", type=Path)
    virtual_parser.add_argument("--output-shortlist", type=Path, required=True)
    virtual_parser.add_argument("--output-summary", type=Path, required=True)

    control_parser = subparsers.add_parser(
        "audit-controls",
        help="Run empirical positive and negative controls through the fitted funnel",
    )
    control_parser.add_argument("screen_csv", type=Path)
    control_parser.add_argument("reconstructed_csv", type=Path)
    control_parser.add_argument("--per-group", type=int, default=25)
    control_parser.add_argument("--ensemble-size", type=int, default=5)
    control_parser.add_argument("--max-iter", type=int, default=80)
    control_parser.add_argument("--output-controls", type=Path, required=True)
    control_parser.add_argument("--output-summary", type=Path, required=True)

    funnel_audit_parser = subparsers.add_parser(
        "audit-screen-funnel",
        help="Audit packaging-gate sensitivity and residue composition without reranking",
    )
    funnel_audit_parser.add_argument("ranked_csv", type=Path)
    funnel_audit_parser.add_argument("shortlist_csv", type=Path)
    funnel_audit_parser.add_argument("screen_csv", type=Path)
    funnel_audit_parser.add_argument("--random-state", type=int, default=42)
    funnel_audit_parser.add_argument("--output-gates", type=Path, required=True)
    funnel_audit_parser.add_argument("--output-composition", type=Path, required=True)
    funnel_audit_parser.add_argument("--output-summary", type=Path, required=True)

    composition_challenge_parser = subparsers.add_parser(
        "select-composition-challenge",
        help="Select a separate packaging-only panel for residues missing from the shortlist",
    )
    composition_challenge_parser.add_argument("ranked_csv", type=Path)
    composition_challenge_parser.add_argument("shortlist_csv", type=Path)
    composition_challenge_parser.add_argument("screen_csv", type=Path)
    composition_challenge_parser.add_argument(
        "--target-residues", nargs="+", default=list("CFIMWY")
    )
    composition_challenge_parser.add_argument("--minimum-pairwise-distance", type=int, default=3)
    composition_challenge_parser.add_argument("--random-state", type=int, default=42)
    composition_challenge_parser.add_argument("--output-panel", type=Path, required=True)
    composition_challenge_parser.add_argument("--output-summary", type=Path, required=True)
    return parser


def _write_json(payload: dict[str, object], output: Path | None) -> None:
    rendered = json.dumps(payload, indent=2, ensure_ascii=False)
    if output is None:
        print(rendered)
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(rendered + "\n", encoding="utf-8")


def _count_one_fastq(
    fastq_path: Path,
    whitelist: set[str],
    counts_path: Path,
    summary_path: Path,
) -> dict[str, object]:
    result = count_fastq(
        fastq_path,
        peptide_whitelist=whitelist,
        counts_output=counts_path,
    )
    _write_json(result, summary_path)
    return result


def _count_one_bowtie2_fastq(
    fastq_path: Path,
    index_prefix: Path,
    whitelist: set[str],
    counts_path: Path,
    summary_path: Path,
    threads: int,
) -> dict[str, object]:
    result = count_bowtie2_fastq(
        fastq_path,
        index_prefix,
        peptide_whitelist=whitelist,
        counts_output=counts_path,
        threads=threads,
    )
    _write_json(result, summary_path)
    return result


def _filter_manifest(
    frame: pd.DataFrame, roles: list[str] | None, excluded_alias_pattern: str | None
) -> pd.DataFrame:
    filtered = frame
    if roles is not None:
        if "project_role" not in filtered:
            raise ValueError("Manifest must contain project_role for role filtering")
        filtered = filtered.loc[filtered["project_role"].isin(roles)]
    if excluded_alias_pattern is not None:
        if "experiment_alias" not in filtered:
            raise ValueError("Manifest must contain experiment_alias for alias filtering")
        filtered = filtered.loc[
            ~filtered["experiment_alias"].str.contains(excluded_alias_pattern, regex=True, na=False)
        ]
    return filtered.copy()


def main() -> None:
    args = _build_parser().parse_args()
    if args.command == "demo":
        summary = run_demo(args.output_dir, random_state=args.random_state)
        print(json.dumps(summary, indent=2, ensure_ascii=False))
        return
    if args.command == "verify-manifest":
        result = verify_manifest(args.manifest, args.root)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        if not result["ok"]:
            raise SystemExit(1)
        return
    if args.command == "freeze-blind-study":
        payload = freeze_blind_study(args.config, args.root, args.output)
        print(json.dumps(payload, indent=2, ensure_ascii=False))
        return
    if args.command == "verify-blind-freeze":
        result = verify_blind_freeze(args.manifest, args.root)
        _write_json(result, args.output)
        if not result["ok"]:
            raise SystemExit(1)
        return
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
    if args.command == "resolve-ena-fastq":
        target = pd.read_csv(args.input)
        resolved = resolve_ena_fastqs(target, args.accession)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        resolved.to_csv(args.output, index=False)
        return
    if args.command == "download-ena-fastq":
        manifest = _filter_manifest(
            pd.read_csv(args.input), args.project_roles, args.exclude_alias_regex
        )
        results = download_fastq_manifest(
            manifest.to_dict(orient="records"), args.output_dir, workers=args.workers
        )
        _write_json({"downloads": results}, args.output)
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
    if args.command == "prepare-bowtie2-reference":
        _write_json(write_bowtie2_reference(args.output_prefix), args.output)
        return
    if args.command == "count-bowtie2":
        whitelist = set(pd.read_csv(args.whitelist_csv)[args.whitelist_column].dropna())
        _write_json(
            count_bowtie2_fastq(
                args.input,
                args.index_prefix,
                peptide_whitelist=whitelist,
                counts_output=args.counts_output,
                threads=args.threads,
            ),
            args.output,
        )
        return
    if args.command == "count-fastq-batch":
        manifest = _filter_manifest(
            pd.read_csv(args.manifest), args.project_roles, args.exclude_alias_regex
        )
        whitelist = set(pd.read_csv(args.whitelist_csv)[args.whitelist_column].dropna())
        args.counts_dir.mkdir(parents=True, exist_ok=True)
        args.summaries_dir.mkdir(parents=True, exist_ok=True)
        futures = {}
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            for row in manifest.to_dict(orient="records"):
                run = str(row["run_accession"])
                fastq_path = args.fastq_dir / f"{run}.fastq.gz"
                if not fastq_path.exists():
                    raise FileNotFoundError(f"Missing FASTQ: {fastq_path}")
                counts_path = args.counts_dir / f"{run}.counts.csv.gz"
                summary_path = args.summaries_dir / f"{run}.json"
                if counts_path.exists() and summary_path.exists() and not args.overwrite:
                    print(f"{run}: cached")
                    continue
                futures[
                    executor.submit(
                        _count_one_fastq,
                        fastq_path,
                        whitelist,
                        counts_path,
                        summary_path,
                    )
                ] = run
            for future in as_completed(futures):
                run = futures[future]
                result = future.result()
                print(f"{run}: {result['reads_examined']} reads")
        return
    if args.command == "count-bowtie2-batch":
        manifest = _filter_manifest(
            pd.read_csv(args.manifest), args.project_roles, args.exclude_alias_regex
        )
        whitelist = set(pd.read_csv(args.whitelist_csv)[args.whitelist_column].dropna())
        args.counts_dir.mkdir(parents=True, exist_ok=True)
        args.summaries_dir.mkdir(parents=True, exist_ok=True)
        futures = {}
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            for row in manifest.to_dict(orient="records"):
                run = str(row["run_accession"])
                fastq_path = args.fastq_dir / f"{run}.fastq.gz"
                if not fastq_path.exists():
                    raise FileNotFoundError(f"Missing FASTQ: {fastq_path}")
                counts_path = args.counts_dir / f"{run}.counts.csv.gz"
                summary_path = args.summaries_dir / f"{run}.json"
                if counts_path.exists() and summary_path.exists() and not args.overwrite:
                    print(f"{run}: cached")
                    continue
                futures[
                    executor.submit(
                        _count_one_bowtie2_fastq,
                        fastq_path,
                        args.index_prefix,
                        whitelist,
                        counts_path,
                        summary_path,
                        args.threads_per_worker,
                    )
                ] = run
            for future in as_completed(futures):
                run = futures[future]
                result = future.result()
                print(f"{run}: {result['reads_examined']} aligned records")
        return
    if args.command == "reconstruct-liver":
        manifest = _filter_manifest(
            pd.read_csv(args.manifest),
            ["Liver", "Virus reference"],
            args.exclude_alias_regex,
        )
        counts, qc = load_count_matrix(manifest, args.counts_dir, args.summaries_dir)
        public = pd.read_csv(args.public_screens, usecols=["AA", "Liver"])
        reconstructed, metrics = reconstruct_liver(counts, qc, public)
        for path in (args.output_reconstruction, args.output_metrics, args.output_qc):
            path.parent.mkdir(parents=True, exist_ok=True)
        reconstructed.to_csv(args.output_reconstruction, index=False)
        metrics.to_csv(args.output_metrics, index=False)
        qc.to_csv(args.output_qc, index=False)
        return
    if args.command == "reconstruct-multiorgan":
        manifest = pd.read_csv(args.manifest)
        roles = ["Brain", "Spinal cord", "Liver", "Heart", "Kidney"]
        organ = manifest.loc[manifest["project_role"].isin(roles)]
        virus = manifest.loc[
            manifest["experiment_alias"].str.contains(
                f"virus_prod{args.virus_round}", regex=False, na=False
            )
        ]
        selected = pd.concat([organ, virus], ignore_index=True)
        counts, qc = load_count_matrix(selected, args.counts_dir, args.summaries_dir)
        reconstructed, metrics = reconstruct_multiorgan(counts, qc, virus_round=args.virus_round)
        for path in (args.output_reconstruction, args.output_metrics, args.output_qc):
            path.parent.mkdir(parents=True, exist_ok=True)
        reconstructed.to_csv(args.output_reconstruction, index=False)
        metrics.to_csv(args.output_metrics, index=False)
        qc.to_csv(args.output_qc, index=False)
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
    if args.command == "benchmark-multiorgan":
        rows = benchmark_multiorgan_animal_holdout(
            args.input,
            model_names=tuple(args.models),
            endpoints=tuple(args.endpoints),
            feature_set=args.feature_set,
            bootstrap_resamples=args.bootstrap_resamples,
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(args.output, index=False)
        return
    if args.command == "benchmark-multitask":
        rows = benchmark_multitask_animal_holdout(
            args.input,
            endpoints=tuple(args.endpoints),
            max_iter=args.max_iter,
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(args.output, index=False)
        return
    if args.command == "benchmark-multitask-ensemble":
        rows = benchmark_multitask_ensemble_animal_holdout(
            args.input,
            endpoints=tuple(args.endpoints),
            ensemble_size=args.ensemble_size,
            max_iter=args.max_iter,
            bootstrap_resamples=args.bootstrap_resamples,
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(args.output, index=False)
        return
    if args.command == "benchmark-cross-animal":
        rows = benchmark_multitask_ensemble_leave_one_animal_out(
            args.input,
            endpoints=tuple(args.endpoints),
            animals=tuple(args.animals),
            ensemble_size=args.ensemble_size,
            max_iter=args.max_iter,
            bootstrap_resamples=args.bootstrap_resamples,
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(args.output, index=False)
        return
    if args.command == "benchmark-masked-multitask":
        rows = benchmark_masked_multitask_animal_holdout(
            args.input,
            endpoints=tuple(args.endpoints),
            max_epochs=args.max_epochs,
            bootstrap_resamples=args.bootstrap_resamples,
            device=args.device,
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
    if args.command == "screen-virtual":
        screen = pd.read_csv(args.screen_csv)
        reconstructed = pd.read_csv(args.reconstructed_csv)
        ranked, shortlist, summary = run_virtual_screen(
            screen,
            reconstructed,
            pool_size=args.pool_size,
            ensemble_size=args.ensemble_size,
            max_iter=args.max_iter,
        )
        output_paths = [args.output_ranked, args.output_shortlist, args.output_summary]
        if args.output_pareto is not None:
            output_paths.append(args.output_pareto)
        for path in output_paths:
            path.parent.mkdir(parents=True, exist_ok=True)
        ranked.to_csv(args.output_ranked, index=False)
        if args.output_pareto is not None:
            ranked.loc[ranked["is_pareto"]].to_csv(args.output_pareto, index=False)
        shortlist.to_csv(args.output_shortlist, index=False)
        _write_json(summary, args.output_summary)
        return
    if args.command == "audit-controls":
        screen = pd.read_csv(args.screen_csv)
        reconstructed = pd.read_csv(args.reconstructed_csv)
        controls, summary = run_control_audit(
            screen,
            reconstructed,
            per_group=args.per_group,
            ensemble_size=args.ensemble_size,
            max_iter=args.max_iter,
        )
        for path in (args.output_controls, args.output_summary):
            path.parent.mkdir(parents=True, exist_ok=True)
        controls.to_csv(args.output_controls, index=False)
        _write_json(summary, args.output_summary)
        return
    if args.command == "audit-screen-funnel":
        ranked = pd.read_csv(args.ranked_csv)
        shortlist = pd.read_csv(args.shortlist_csv)
        screen = pd.read_csv(args.screen_csv)
        offsets: dict[str, float] = {}
        packaging_threshold = None
        for gate, coverage in (("sensitivity_90_lcb", 0.90), ("strict_95_lcb", 0.95)):
            _, threshold, offset, _ = fit_calibrated_packaging_model(
                screen,
                random_state=args.random_state,
                coverage=coverage,
            )
            if packaging_threshold is not None and threshold != packaging_threshold:
                raise RuntimeError("Packaging threshold changed between calibration runs")
            packaging_threshold = threshold
            offsets[gate] = offset
        offsets["sensitivity_point_prediction"] = 0.0
        gate_audit = build_gate_sensitivity_audit(
            ranked,
            packaging_threshold=float(packaging_threshold),
            lower_bound_offsets=offsets,
        )
        composition_audit = build_composition_audit(ranked, shortlist, screen)
        for path in (args.output_gates, args.output_composition, args.output_summary):
            path.parent.mkdir(parents=True, exist_ok=True)
        gate_audit.to_csv(args.output_gates, index=False)
        composition_audit.to_csv(args.output_composition, index=False)
        _write_json(summarize_funnel_audit(gate_audit, composition_audit), args.output_summary)
        return
    if args.command == "select-composition-challenge":
        ranked = pd.read_csv(args.ranked_csv)
        shortlist = pd.read_csv(args.shortlist_csv)
        screen = pd.read_csv(args.screen_csv)
        _, packaging_threshold, strict_95_offset, _ = fit_calibrated_packaging_model(
            screen,
            random_state=args.random_state,
            coverage=0.95,
        )
        _, threshold_90, sensitivity_90_offset, _ = fit_calibrated_packaging_model(
            screen,
            random_state=args.random_state,
            coverage=0.90,
        )
        if threshold_90 != packaging_threshold:
            raise RuntimeError("Packaging threshold changed between calibration runs")
        panel, summary = select_composition_challenge_panel(
            ranked,
            shortlist,
            packaging_threshold=packaging_threshold,
            strict_95_offset=strict_95_offset,
            sensitivity_90_offset=sensitivity_90_offset,
            target_residues=args.target_residues,
            minimum_pairwise_distance=args.minimum_pairwise_distance,
        )
        for path in (args.output_panel, args.output_summary):
            path.parent.mkdir(parents=True, exist_ok=True)
        panel.to_csv(args.output_panel, index=False)
        _write_json(summary, args.output_summary)
        return
    raise RuntimeError(f"Unhandled command: {args.command}")


if __name__ == "__main__":
    main()
