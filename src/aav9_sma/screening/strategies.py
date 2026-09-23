"""Parallel joint and sequential organ-screening strategy utilities.

The existing model predicts all organs together. This module changes only the
decision logic applied to those predictions: a configurable sequential funnel
prioritizes spinal cord first, then brain, then off-target organs. Keeping this
logic separate makes it possible to compare both strategies on exactly the
same predictions without retraining or regenerating the virtual pool.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd

from aav9_sma.constants import AMINO_ACIDS
from aav9_sma.screening.virtual import hamming_distance


@dataclass(frozen=True)
class SequentialStage:
    """One configurable organ filter in a sequential screening funnel."""

    name: str
    column: str
    direction: Literal["maximize", "minimize"]
    retain_fraction: float


def default_sma_sequential_stages(
    *,
    spinal_retain_fraction: float = 0.25,
    brain_retain_fraction: float = 0.50,
    liver_retain_fraction: float = 0.50,
    heart_retain_fraction: float = 0.75,
    kidney_retain_fraction: float = 0.75,
) -> tuple[SequentialStage, ...]:
    """Return the default SMA organ order with independently tunable stringency."""
    return (
        SequentialStage(
            "spinal_cord", "pred_spinal_cord_mouse", "maximize", spinal_retain_fraction
        ),
        SequentialStage("brain", "pred_brain_mouse", "maximize", brain_retain_fraction),
        SequentialStage("liver", "pred_liver_mouse", "minimize", liver_retain_fraction),
        SequentialStage("heart", "pred_heart_mouse", "minimize", heart_retain_fraction),
        SequentialStage("kidney", "pred_kidney_mouse", "minimize", kidney_retain_fraction),
    )


def _validate_stages(frame: pd.DataFrame, stages: tuple[SequentialStage, ...]) -> None:
    if not stages:
        raise ValueError("at least one sequential stage is required")
    names = [stage.name for stage in stages]
    if len(names) != len(set(names)):
        raise ValueError("sequential stage names must be unique")
    for stage in stages:
        if stage.column not in frame:
            raise ValueError(f"Missing sequential-stage column: {stage.column}")
        if not 0 < stage.retain_fraction <= 1:
            raise ValueError(
                f"retain_fraction for {stage.name} must be greater than 0 and at most 1"
            )


def apply_sequential_organ_funnel(
    ranked: pd.DataFrame,
    stages: tuple[SequentialStage, ...] | None = None,
    *,
    minimum_training_distance: int = 2,
) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    """Apply ordered quantile filters to one shared set of organ predictions.

    Packaging and training-distance eligibility are evaluated before organ
    stages. At each stage, the quantile cutoff is estimated only among the
    candidates that survived all previous stages. Ties at the cutoff are kept,
    so the observed retained fraction can be slightly larger than requested.
    """
    stages = stages or default_sma_sequential_stages()
    _validate_stages(ranked, stages)
    required = {"passes_packaging_gate", "training_distance_lower_bound"}
    missing = required - set(ranked.columns)
    if missing:
        raise ValueError(f"Missing sequential-funnel columns: {sorted(missing)}")

    output = ranked.copy()
    base_mask = output["passes_packaging_gate"].astype(bool) & output[
        "training_distance_lower_bound"
    ].ge(minimum_training_distance)
    output["passes_sequential_base"] = base_mask
    output["sequential_stages_passed"] = 0
    active = base_mask.copy()
    audit: list[dict[str, object]] = [
        {
            "stage": "base",
            "column": None,
            "direction": None,
            "requested_retain_fraction": None,
            "input_rows": len(output),
            "survivor_rows": int(active.sum()),
            "cutoff": None,
        }
    ]

    for stage_number, stage in enumerate(stages, start=1):
        stage_column = pd.to_numeric(output[stage.column], errors="coerce")
        finite_active = active & np.isfinite(stage_column)
        input_rows = int(active.sum())
        finite_rows = int(finite_active.sum())
        pass_column = f"passes_sequential_{stage_number}_{stage.name}"
        output[pass_column] = False
        if finite_rows == 0:
            active[:] = False
            cutoff = float("nan")
        else:
            if stage.direction == "maximize":
                cutoff = float(
                    stage_column.loc[finite_active].quantile(
                        1.0 - stage.retain_fraction, interpolation="higher"
                    )
                )
                stage_pass = finite_active & stage_column.ge(cutoff)
            else:
                cutoff = float(
                    stage_column.loc[finite_active].quantile(
                        stage.retain_fraction, interpolation="lower"
                    )
                )
                stage_pass = finite_active & stage_column.le(cutoff)
            active = stage_pass
            output.loc[active, "sequential_stages_passed"] = stage_number
        output.loc[active, pass_column] = True
        audit.append(
            {
                "stage": stage.name,
                "column": stage.column,
                "direction": stage.direction,
                "requested_retain_fraction": stage.retain_fraction,
                "input_rows": input_rows,
                "finite_input_rows": finite_rows,
                "survivor_rows": int(active.sum()),
                "observed_retain_fraction": (
                    float(active.sum() / input_rows) if input_rows else 0.0
                ),
                "cutoff": cutoff,
            }
        )

    output["passes_sequential_funnel"] = active
    sort_columns = [stage.column for stage in stages]
    ascending = [stage.direction == "minimize" for stage in stages]
    ordered_indices = output.loc[active].sort_values(
        sort_columns, ascending=ascending, kind="stable"
    ).index
    output["sequential_rank"] = pd.Series(pd.NA, index=output.index, dtype="Int64")
    output.loc[ordered_indices, "sequential_rank"] = np.arange(1, len(ordered_indices) + 1)
    return output, audit


def select_sequential_shortlist(
    annotated: pd.DataFrame,
    *,
    candidate_count: int = 30,
    minimum_pairwise_distance: int = 3,
) -> pd.DataFrame:
    """Greedily diversify the lexicographically ordered sequential survivors."""
    if candidate_count < 1:
        raise ValueError("candidate_count must be positive")
    required = {"AA", "passes_sequential_funnel", "sequential_rank"}
    missing = required - set(annotated.columns)
    if missing:
        raise ValueError(f"Missing sequential-shortlist columns: {sorted(missing)}")

    eligible = annotated.loc[annotated["passes_sequential_funnel"]].sort_values(
        "sequential_rank", kind="stable"
    )
    selected_indices: list[int] = []
    selected_peptides: list[str] = []
    for index, row in eligible.iterrows():
        peptide = str(row["AA"])
        if any(
            hamming_distance(peptide, selected) < minimum_pairwise_distance
            for selected in selected_peptides
        ):
            continue
        selected_indices.append(index)
        selected_peptides.append(peptide)
        if len(selected_indices) == candidate_count:
            break

    shortlist = annotated.loc[selected_indices].copy()
    shortlist["selection_strategy"] = "sequential_spinal_first"
    shortlist["strategy_selection_rank"] = np.arange(1, len(shortlist) + 1)
    shortlist["minimum_pairwise_hamming_required"] = minimum_pairwise_distance
    return shortlist.reset_index(drop=True)


def _residue_frequencies(peptides: pd.Series) -> np.ndarray:
    counts = Counter("".join(peptides.dropna().astype(str)))
    total = sum(counts.values())
    if total == 0:
        return np.zeros(len(AMINO_ACIDS), dtype=float)
    return np.array([counts[residue] / total for residue in AMINO_ACIDS], dtype=float)


def _jensen_shannon_divergence(left: np.ndarray, right: np.ndarray) -> float:
    midpoint = 0.5 * (left + right)

    def kl_divergence(values: np.ndarray, reference: np.ndarray) -> float:
        observed = values > 0
        return float(np.sum(values[observed] * np.log2(values[observed] / reference[observed])))

    return 0.5 * kl_divergence(left, midpoint) + 0.5 * kl_divergence(right, midpoint)


def _minimum_pairwise_hamming(peptides: pd.Series) -> int | None:
    values = peptides.dropna().astype(str).tolist()
    if len(values) < 2:
        return None
    return min(
        hamming_distance(left, right)
        for index, left in enumerate(values)
        for right in values[index + 1 :]
    )


def compare_strategy_shortlists(
    joint_shortlist: pd.DataFrame,
    sequential_shortlist: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Quantify overlap, predicted phenotypes, ranks, and sequence composition."""
    for name, frame in (
        ("joint_shortlist", joint_shortlist),
        ("sequential_shortlist", sequential_shortlist),
    ):
        if "AA" not in frame:
            raise ValueError(f"{name} must contain an AA column")
        if frame["AA"].duplicated().any():
            raise ValueError(f"{name} contains duplicate AA sequences")

    joint = joint_shortlist.copy()
    sequential = sequential_shortlist.copy()
    joint["joint_selected"] = True
    sequential["sequential_selected"] = True
    joint = joint.sort_values("display_score", ascending=False, kind="stable")
    joint["joint_strategy_rank"] = np.arange(1, len(joint) + 1)
    if "strategy_selection_rank" not in sequential:
        sequential["strategy_selection_rank"] = np.arange(1, len(sequential) + 1)

    joint_columns = [
        "AA",
        "joint_selected",
        "joint_strategy_rank",
        *[
            column
            for column in (
                "selection_group",
                "display_score",
                "f_sma_target",
                "pred_spinal_cord_mouse",
                "pred_brain_mouse",
                "pred_liver_mouse",
                "pred_heart_mouse",
                "pred_kidney_mouse",
                "organ_uncertainty_mean",
            )
            if column in joint
        ],
    ]
    sequential_columns = [
        "AA",
        "sequential_selected",
        "strategy_selection_rank",
        *[
            column
            for column in (
                "sequential_rank",
                "display_score",
                "f_sma_target",
                "pred_spinal_cord_mouse",
                "pred_brain_mouse",
                "pred_liver_mouse",
                "pred_heart_mouse",
                "pred_kidney_mouse",
                "organ_uncertainty_mean",
            )
            if column in sequential
        ],
    ]
    comparison = joint[joint_columns].merge(
        sequential[sequential_columns],
        on="AA",
        how="outer",
        suffixes=("_joint", "_sequential"),
        validate="one_to_one",
    )
    comparison["joint_selected"] = comparison["joint_selected"].eq(True)
    comparison["sequential_selected"] = comparison["sequential_selected"].eq(True)
    comparison["selection_relationship"] = np.select(
        [
            comparison["joint_selected"] & comparison["sequential_selected"],
            comparison["joint_selected"],
            comparison["sequential_selected"],
        ],
        ["shared", "joint_only", "sequential_only"],
        default="unselected",
    )

    joint_set = set(joint["AA"].astype(str))
    sequential_set = set(sequential["AA"].astype(str))
    intersection = joint_set & sequential_set
    union = joint_set | sequential_set
    common = comparison.loc[comparison["selection_relationship"].eq("shared")]
    rank_correlation = None
    if len(common) >= 2:
        rank_correlation = float(
            common["joint_strategy_rank"].corr(
                common["strategy_selection_rank"], method="spearman"
            )
        )

    metric_columns = [
        column
        for column in (
            "pred_spinal_cord_mouse",
            "pred_brain_mouse",
            "pred_liver_mouse",
            "pred_heart_mouse",
            "pred_kidney_mouse",
            "f_sma_target",
            "display_score",
            "organ_uncertainty_mean",
        )
        if column in joint and column in sequential
    ]
    phenotype_summary = {
        column: {
            "joint_median": float(pd.to_numeric(joint[column], errors="coerce").median()),
            "sequential_median": float(
                pd.to_numeric(sequential[column], errors="coerce").median()
            ),
            "sequential_minus_joint": float(
                pd.to_numeric(sequential[column], errors="coerce").median()
                - pd.to_numeric(joint[column], errors="coerce").median()
            ),
        }
        for column in metric_columns
    }
    joint_frequencies = _residue_frequencies(joint["AA"])
    sequential_frequencies = _residue_frequencies(sequential["AA"])
    summary: dict[str, object] = {
        "joint_rows": len(joint),
        "sequential_rows": len(sequential),
        "shared_rows": len(intersection),
        "joint_only_rows": len(joint_set - sequential_set),
        "sequential_only_rows": len(sequential_set - joint_set),
        "jaccard_similarity": len(intersection) / len(union) if union else 1.0,
        "overlap_coefficient": (
            len(intersection) / min(len(joint_set), len(sequential_set))
            if joint_set and sequential_set
            else 0.0
        ),
        "common_sequence_rank_spearman": rank_correlation,
        "predicted_phenotype_medians": phenotype_summary,
        "residue_frequency_jensen_shannon_divergence": _jensen_shannon_divergence(
            joint_frequencies, sequential_frequencies
        ),
        "minimum_pairwise_hamming": {
            "joint": _minimum_pairwise_hamming(joint["AA"]),
            "sequential": _minimum_pairwise_hamming(sequential["AA"]),
        },
        "interpretation": (
            "Overlap measures agreement; phenotype medians show each strategy's selection "
            "pressure; Jensen-Shannon divergence quantifies amino-acid composition shift."
        ),
    }
    return comparison.sort_values(
        ["selection_relationship", "joint_strategy_rank", "strategy_selection_rank"],
        kind="stable",
    ).reset_index(drop=True), summary


def run_sequential_strategy(
    ranked: pd.DataFrame,
    stages: tuple[SequentialStage, ...] | None = None,
    *,
    candidate_count: int = 30,
    minimum_training_distance: int = 2,
    minimum_pairwise_distance: int = 3,
) -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, object]]]:
    """Apply the sequential funnel and return annotations, shortlist, and audit."""
    annotated, audit = apply_sequential_organ_funnel(
        ranked,
        stages,
        minimum_training_distance=minimum_training_distance,
    )
    shortlist = select_sequential_shortlist(
        annotated,
        candidate_count=candidate_count,
        minimum_pairwise_distance=minimum_pairwise_distance,
    )
    return annotated, shortlist, audit
