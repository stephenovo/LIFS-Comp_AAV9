"""Reconstruct sequence-linked Fit4Function organ enrichment from raw counts."""

from __future__ import annotations

import json
import re
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd

ORGAN_ALIAS = re.compile(
    r"^hammerhead_(?P<organ>[A-Za-z_]+)_a(?P<animal>\d+)_r(?P<technical>\d+)$"
)
VIRUS_ALIAS = re.compile(
    r"^hammerhead_start(?:er)?_virus_prod(?P<production>\d+)_(?:r)?(?P<technical>\d+)$"
)


@dataclass(frozen=True)
class RunDesign:
    kind: str
    endpoint: str
    biological_replicate: int | None
    technical_replicate: int
    production_round: int | None


def parse_hammerhead_alias(alias: str) -> RunDesign:
    """Parse the biological and technical replicate encoded in an SRA alias."""
    organ = ORGAN_ALIAS.fullmatch(alias)
    if organ:
        endpoint = organ.group("organ").replace("_", " ")
        return RunDesign(
            kind="organ",
            endpoint=endpoint,
            biological_replicate=int(organ.group("animal")),
            technical_replicate=int(organ.group("technical")),
            production_round=None,
        )
    virus = VIRUS_ALIAS.fullmatch(alias)
    if virus:
        production = int(virus.group("production"))
        return RunDesign(
            kind="virus",
            endpoint="Virus DNA",
            biological_replicate=None,
            technical_replicate=int(virus.group("technical")),
            production_round=production,
        )
    raise ValueError(f"Unrecognized Hammerhead experiment alias: {alias}")


def _safe_log2_ratio(numerator: pd.Series, denominator: pd.Series) -> pd.Series:
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log2(numerator / denominator)


def _load_run_summary(summary_path: Path) -> dict[str, object]:
    return json.loads(summary_path.read_text(encoding="utf-8"))


def load_count_matrix(
    manifest: pd.DataFrame,
    counts_dir: str | Path,
    summaries_dir: str | Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load dense whitelist counts and extraction QC for every run in a manifest."""
    counts_dir = Path(counts_dir)
    summaries_dir = Path(summaries_dir)
    count_columns: list[pd.Series] = []
    qc_rows: list[dict[str, object]] = []
    expected_index: pd.Index | None = None

    for row in manifest.to_dict(orient="records"):
        run = str(row["run_accession"])
        count_candidates = [counts_dir / f"{run}.counts.csv.gz", counts_dir / f"{run}.counts.csv"]
        count_path = next((path for path in count_candidates if path.exists()), None)
        summary_path = summaries_dir / f"{run}.json"
        if count_path is None or not summary_path.exists():
            raise FileNotFoundError(f"Missing count table or summary for {run}")
        count_frame = pd.read_csv(count_path)
        if count_frame["peptide_7mer"].duplicated().any():
            raise ValueError(f"Duplicate peptide rows in {count_path}")
        series = count_frame.set_index("peptide_7mer")["raw_read_count"].astype("int64")
        if expected_index is None:
            expected_index = series.index
        elif not series.index.equals(expected_index):
            raise ValueError(f"Whitelist order differs in {count_path}")
        series.name = run
        count_columns.append(series)

        summary = _load_run_summary(summary_path)
        design = parse_hammerhead_alias(str(row["experiment_alias"]))
        valid_reads = int(summary["status_counts"].get("valid", 0))
        whitelist_reads = int(summary["whitelist_matched_reads"])
        if whitelist_reads != int(series.sum()):
            raise ValueError(f"Count total does not match extraction summary for {run}")
        qc_rows.append(
            {
                **{key: row.get(key) for key in manifest.columns},
                **asdict(design),
                "reads_examined": int(summary["reads_examined"]),
                "valid_reads_all_7mer": valid_reads,
                "whitelist_matched_reads": whitelist_reads,
                "valid_fraction": float(summary["valid_fraction"]),
                "whitelist_share_of_valid": (
                    whitelist_reads / valid_reads if valid_reads else np.nan
                ),
            }
        )

    if not count_columns:
        raise ValueError("Manifest contains no runs")
    return pd.concat(count_columns, axis=1), pd.DataFrame(qc_rows)


def _aggregate_rpm(
    count_matrix: pd.DataFrame, qc: pd.DataFrame, denominator_mode: str
) -> dict[str, pd.Series]:
    if denominator_mode == "all_valid":
        denominators = qc.set_index("run_accession")["valid_reads_all_7mer"]
    elif denominator_mode == "whitelist":
        denominators = qc.set_index("run_accession")["whitelist_matched_reads"]
    else:
        raise ValueError("denominator_mode must be all_valid or whitelist")
    rpm = count_matrix.div(denominators.reindex(count_matrix.columns), axis=1) * 1_000_000
    aggregated: dict[str, pd.Series] = {}
    for (kind, replicate), group in qc.groupby(["kind", "biological_replicate"], dropna=False):
        if kind != "organ" or pd.isna(replicate):
            continue
        runs = group.sort_values("technical_replicate")["run_accession"].tolist()
        endpoint = str(group["endpoint"].iloc[0]).lower().replace(" ", "_")
        aggregated[f"{endpoint}_a{int(replicate)}"] = rpm[runs].mean(axis=1)
    for production, group in qc.loc[qc["kind"] == "virus"].groupby("production_round"):
        runs = group.sort_values("technical_replicate")["run_accession"].tolist()
        aggregated[f"virus_prod{int(production)}"] = rpm[runs].mean(axis=1)
    return aggregated


def _validation_row(public: pd.Series, reconstructed: pd.Series, label: str) -> dict[str, object]:
    finite = np.isfinite(public) & np.isfinite(reconstructed)
    observed = public.loc[finite].astype(float)
    predicted = reconstructed.loc[finite].astype(float)
    if len(observed) < 3:
        raise ValueError(f"Not enough finite observations for {label}")
    slope, intercept = np.polyfit(predicted, observed, 1)
    calibrated = slope * predicted + intercept
    return {
        "reconstruction": label,
        "finite_sequences": int(finite.sum()),
        "pearson_r": float(observed.corr(predicted, method="pearson")),
        "spearman_r": float(observed.corr(predicted, method="spearman")),
        "calibration_slope": float(slope),
        "calibration_intercept": float(intercept),
        "calibrated_rmse": float(np.sqrt(np.mean((calibrated - observed) ** 2))),
    }


def reconstruct_liver(
    count_matrix: pd.DataFrame,
    qc: pd.DataFrame,
    public_screens: pd.DataFrame,
    denominator_modes: Iterable[str] = ("all_valid", "whitelist"),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build liver/virus RPM scores and compare plausible aggregations to public Liver labels."""
    required_public = {"AA", "Liver"}
    if not required_public.issubset(public_screens.columns):
        raise ValueError("Public screen must contain AA and Liver")
    public = public_screens.set_index("AA")["Liver"].reindex(count_matrix.index)
    output = pd.DataFrame(index=count_matrix.index)
    output.index.name = "AA"
    output["public_Liver"] = public
    metrics: list[dict[str, object]] = []

    for denominator_mode in denominator_modes:
        aggregated = _aggregate_rpm(count_matrix, qc, denominator_mode)
        liver_keys = sorted(key for key in aggregated if key.startswith("liver_a"))
        virus_keys = sorted(key for key in aggregated if key.startswith("virus_prod"))
        if not liver_keys or not virus_keys:
            raise ValueError("Liver and virus runs are both required")
        prefix = f"rpm_{denominator_mode}"
        for key, values in aggregated.items():
            output[f"{prefix}__{key}"] = values

        liver_groups = {"liver_all_animals": liver_keys}
        if len(liver_keys) >= 3:
            liver_groups["liver_animals_1_3"] = liver_keys[:3]
        for key in liver_keys:
            liver_groups[key] = [key]
        virus_groups = {key: [key] for key in virus_keys}
        if len(virus_keys) > 1:
            virus_groups["virus_available_rounds"] = virus_keys

        for liver_label, selected_liver in liver_groups.items():
            liver_rpm = pd.concat([aggregated[key] for key in selected_liver], axis=1).mean(axis=1)
            output[f"{prefix}__{liver_label}"] = liver_rpm
            for virus_label, selected_virus in virus_groups.items():
                virus_rpm = pd.concat(
                    [aggregated[key] for key in selected_virus], axis=1
                ).mean(axis=1)
                enrichment = _safe_log2_ratio(liver_rpm, virus_rpm)
                label = f"{prefix}__{liver_label}__over__{virus_label}"
                output[label] = enrichment
                metrics.append(_validation_row(public, enrichment, label))

    metrics_frame = pd.DataFrame(metrics).sort_values(
        ["pearson_r", "spearman_r"], ascending=False
    )
    return output.reset_index(), metrics_frame.reset_index(drop=True)


def reconstruct_multiorgan(
    count_matrix: pd.DataFrame,
    qc: pd.DataFrame,
    denominator_modes: Iterable[str] = ("all_valid", "whitelist"),
    virus_round: int = 2,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Reconstruct animal-level organ enrichments using a validated virus round."""
    output = pd.DataFrame(index=count_matrix.index)
    output.index.name = "AA"
    metrics: list[dict[str, object]] = []

    for denominator_mode in denominator_modes:
        aggregated = _aggregate_rpm(count_matrix, qc, denominator_mode)
        virus_key = f"virus_prod{virus_round}"
        if virus_key not in aggregated:
            raise ValueError(f"Missing virus reference: {virus_key}")
        virus_rpm = aggregated[virus_key]
        prefix = f"rpm_{denominator_mode}"
        output[f"{prefix}__{virus_key}"] = virus_rpm
        endpoints = sorted(
            {
                key.rsplit("_a", maxsplit=1)[0]
                for key in aggregated
                if not key.startswith("virus_prod")
            }
        )
        for endpoint in endpoints:
            animal_keys = sorted(
                key for key in aggregated if key.startswith(f"{endpoint}_a")
            )
            animal_enrichments: dict[str, pd.Series] = {}
            for key in animal_keys:
                rpm = aggregated[key]
                enrichment = _safe_log2_ratio(rpm, virus_rpm)
                output[f"{prefix}__{key}"] = rpm
                enrichment_key = (
                    f"log2enr_{denominator_mode}__{key}__over__{virus_key}"
                )
                output[enrichment_key] = enrichment
                animal_enrichments[key] = enrichment

            groups = {"all_animals": animal_keys}
            if len(animal_keys) >= 3:
                groups["animals_1_3"] = animal_keys[:3]
            for group_name, selected in groups.items():
                mean_rpm = pd.concat([aggregated[key] for key in selected], axis=1).mean(
                    axis=1
                )
                output[f"{prefix}__{endpoint}_{group_name}"] = mean_rpm
                output[
                    f"log2enr_{denominator_mode}__{endpoint}_{group_name}__over__{virus_key}"
                ] = _safe_log2_ratio(mean_rpm, virus_rpm)

            animal_frame = pd.DataFrame(animal_enrichments).replace(
                [np.inf, -np.inf], np.nan
            )
            correlations = animal_frame.corr(method="pearson")
            upper = correlations.where(
                np.triu(np.ones(correlations.shape), k=1).astype(bool)
            ).stack()
            metrics.append(
                {
                    "denominator_mode": denominator_mode,
                    "endpoint": endpoint,
                    "virus_round": virus_round,
                    "animals": len(animal_keys),
                    "median_pairwise_animal_pearson": float(upper.median()),
                    "minimum_pairwise_animal_pearson": float(upper.min()),
                    "maximum_pairwise_animal_pearson": float(upper.max()),
                    "mean_finite_sequences_per_animal": float(
                        np.isfinite(animal_frame).sum().mean()
                    ),
                }
            )

    return output.reset_index(), pd.DataFrame(metrics)
