"""Reproducible virtual generation, prediction, and diverse shortlisting."""

from __future__ import annotations

import math
import re
from collections.abc import Collection

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.preprocessing import StandardScaler

from aav9_sma.constants import AMINO_ACIDS, PEPTIDE_LENGTH
from aav9_sma.features.encode import one_hot_7mer
from aav9_sma.models.baseline import build_shared_mlp
from aav9_sma.models.evaluate import MULTIORGAN_ENDPOINTS, sequence_distance_split
from aav9_sma.screening.score import rank_candidates

HUMAN_LIVER_TASKS = ("HepG2_bind", "HepG2_tr", "THLE_bind", "THLE_tr")
HYDROPHOBIC_RESIDUES = frozenset("AVILMFWY")


def generate_candidate_peptides(
    count: int,
    excluded: Collection[str] = (),
    random_state: int = 42,
) -> list[str]:
    """Uniformly sample unique legal 7-mers outside the observed library."""
    if count < 1:
        raise ValueError("count must be positive")
    maximum = len(AMINO_ACIDS) ** PEPTIDE_LENGTH - len(set(excluded))
    if count > maximum:
        raise ValueError("count exceeds the available 7-mer sequence space")
    rng = np.random.default_rng(random_state)
    amino_acids = np.array(list(AMINO_ACIDS))
    seen = set(excluded)
    candidates: list[str] = []
    while len(candidates) < count:
        batch_size = max(1024, int((count - len(candidates)) * 1.1))
        encoded = rng.integers(0, len(amino_acids), size=(batch_size, PEPTIDE_LENGTH))
        for row in encoded:
            peptide = "".join(amino_acids[row])
            if peptide in seen:
                continue
            seen.add(peptide)
            candidates.append(peptide)
            if len(candidates) == count:
                break
    return candidates


def training_distance_lower_bound(
    candidates: Collection[str], references: Collection[str]
) -> np.ndarray:
    """Return 0, 1, or 2, where 2 means distance two or greater."""
    reference_set = set(references)
    one_wildcard_signatures = {
        peptide[:position] + "*" + peptide[position + 1 :]
        for peptide in reference_set
        for position in range(PEPTIDE_LENGTH)
    }
    distances = np.full(len(candidates), 2, dtype=np.int8)
    for row, peptide in enumerate(candidates):
        if peptide in reference_set:
            distances[row] = 0
            continue
        if any(
            peptide[:position] + "*" + peptide[position + 1 :]
            in one_wildcard_signatures
            for position in range(PEPTIDE_LENGTH)
        ):
            distances[row] = 1
    return distances


def annotate_sequence_liabilities(peptides: Collection[str]) -> pd.DataFrame:
    """Calculate transparent sequence-composition annotations without filtering."""

    def longest_run(peptide: str) -> int:
        return max(len(match.group(0)) for match in re.finditer(r"(.)\1*", peptide))

    return pd.DataFrame(
        {
            "max_homopolymer_run": [longest_run(peptide) for peptide in peptides],
            "hydrophobic_fraction": [
                sum(residue in HYDROPHOBIC_RESIDUES for residue in peptide)
                / PEPTIDE_LENGTH
                for peptide in peptides
            ],
            "net_charge_proxy": [
                peptide.count("K")
                + peptide.count("R")
                - peptide.count("D")
                - peptide.count("E")
                for peptide in peptides
            ],
            "cysteine_count": [peptide.count("C") for peptide in peptides],
            "proline_count": [peptide.count("P") for peptide in peptides],
            "has_n_linked_motif": [
                bool(re.search(r"N[^P][ST]", peptide)) for peptide in peptides
            ],
        }
    )


def fit_calibrated_packaging_model(
    screen: pd.DataFrame,
    random_state: int = 42,
    coverage: float = 0.95,
) -> tuple[Ridge, float, float, dict[str, float | int | str]]:
    """Fit Ridge and estimate a one-sided lower-bound offset on held-out sequences."""
    production = pd.to_numeric(screen["Production2"], errors="coerce")
    observed = screen.loc[np.isfinite(production)].copy()
    peptides = observed["AA"].tolist()
    targets = pd.to_numeric(observed["Production2"], errors="coerce").to_numpy(float)
    features = one_hot_7mer(peptides)
    train_rows, calibration_rows = sequence_distance_split(
        peptides, random_state=random_state, minimum_hamming_distance=2
    )
    calibration_model = Ridge(alpha=1.0)
    calibration_model.fit(features[train_rows], targets[train_rows])
    calibration_predictions = calibration_model.predict(features[calibration_rows])
    overprediction = calibration_predictions - targets[calibration_rows]
    quantile_offset = max(
        0.0,
        float(np.quantile(overprediction, coverage, method="higher")),
    )
    lower_bounds = calibration_predictions - quantile_offset
    final_model = Ridge(alpha=1.0).fit(features, targets)
    threshold = float(np.median(targets))
    metrics: dict[str, float | int | str] = {
        "model": "ridge",
        "split": "distance-2-sequence-calibration",
        "train_rows": int(train_rows.sum()),
        "calibration_rows": int(calibration_rows.sum()),
        "pearson_r": float(np.corrcoef(targets[calibration_rows], calibration_predictions)[0, 1]),
        "r2": float(r2_score(targets[calibration_rows], calibration_predictions)),
        "mae": float(mean_absolute_error(targets[calibration_rows], calibration_predictions)),
        "lower_bound_target_coverage": coverage,
        "lower_bound_empirical_coverage": float(
            np.mean(targets[calibration_rows] >= lower_bounds)
        ),
        "lower_bound_offset": quantile_offset,
        "packaging_threshold": threshold,
        "threshold_definition": "median observed Production2",
    }
    return final_model, threshold, quantile_offset, metrics


def predict_organ_ensemble(
    reconstructed: pd.DataFrame,
    candidate_features: np.ndarray,
    ensemble_size: int = 5,
    random_state: int = 42,
    max_iter: int = 80,
) -> tuple[np.ndarray, np.ndarray, dict[str, int | list[int]]]:
    """Fit shared multi-output MLPs on animals 1-3 and predict candidates."""
    target_columns = [
        f"log2enr_whitelist__{endpoint}_animals_1_3__over__virus_prod2"
        for endpoint in MULTIORGAN_ENDPOINTS
    ]
    targets = reconstructed[target_columns].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    complete = np.isfinite(targets).all(axis=1)
    train_features = one_hot_7mer(reconstructed.loc[complete, "AA"].tolist())
    scaler = StandardScaler().fit(targets[complete])
    scaled_targets = scaler.transform(targets[complete])
    predictions = []
    iterations = []
    for member in range(ensemble_size):
        model = build_shared_mlp(
            random_state=random_state + member,
            max_iter=max_iter,
        )
        model.fit(train_features, scaled_targets)
        predictions.append(scaler.inverse_transform(model.predict(candidate_features)))
        iterations.append(int(model.n_iter_))
    stacked = np.stack(predictions)
    return (
        stacked.mean(axis=0),
        stacked.std(axis=0),
        {"train_rows": int(complete.sum()), "iterations": iterations},
    )


def predict_human_liver_annotations(
    screen: pd.DataFrame, candidate_features: np.ndarray
) -> tuple[pd.DataFrame, float]:
    """Predict human liver-cell labels as secondary annotations only."""
    output: dict[str, np.ndarray] = {}
    numeric = screen[list(HUMAN_LIVER_TASKS)].apply(pd.to_numeric, errors="coerce")
    numeric = numeric.replace([np.inf, -np.inf], np.nan)
    observed_means = numeric.mean(axis=1, skipna=True)
    for task in HUMAN_LIVER_TASKS:
        targets = numeric[task].to_numpy(float)
        observed = np.isfinite(targets)
        model = Ridge(alpha=1.0).fit(
            one_hot_7mer(screen.loc[observed, "AA"].tolist()), targets[observed]
        )
        output[f"pred_{task.lower()}"] = model.predict(candidate_features)
    frame = pd.DataFrame(output)
    frame["pred_human_liver_mean"] = frame.mean(axis=1)
    warning_threshold = float(observed_means.quantile(0.75))
    frame["human_liver_warning"] = frame["pred_human_liver_mean"] >= warning_threshold
    return frame, warning_threshold


def add_weight_sensitivity(
    ranked: pd.DataFrame,
    top_fraction: float = 0.05,
) -> pd.DataFrame:
    """Measure rank stability over a predeclared ±0.10 score-weight grid."""
    output = ranked.copy()
    eligible_indices = np.flatnonzero(output["passes_packaging_gate"].to_numpy(bool))
    counts = np.zeros(len(output), dtype=int)
    percentile_sum = np.zeros(len(output), dtype=float)
    scenarios = 0
    if not len(eligible_indices):
        output["weight_stability_top_fraction"] = 0.0
        output["weight_mean_percentile"] = 0.0
        return output
    top_count = max(1, math.ceil(len(eligible_indices) * top_fraction))
    for cns_weight in (0.35, 0.45, 0.55):
        for liver_weight in (0.25, 0.35, 0.45):
            for off_weight in (0.10, 0.20, 0.30):
                total = cns_weight + liver_weight + off_weight
                scores = (
                    cns_weight / total * output.loc[eligible_indices, "f_cns"].to_numpy()
                    - liver_weight / total * output.loc[eligible_indices, "f_liv"].to_numpy()
                    - off_weight / total * output.loc[eligible_indices, "f_off"].to_numpy()
                )
                local_order = np.argsort(-scores, kind="stable")
                counts[eligible_indices[local_order[:top_count]]] += 1
                percentiles = np.empty(len(local_order), dtype=float)
                percentiles[local_order] = 1.0 - np.arange(len(local_order)) / max(
                    1, len(local_order) - 1
                )
                percentile_sum[eligible_indices] += percentiles
                scenarios += 1
    output["weight_stability_top_fraction"] = counts / scenarios
    output["weight_mean_percentile"] = percentile_sum / scenarios
    return output


def hamming_distance(left: str, right: str) -> int:
    return sum(a != b for a, b in zip(left, right, strict=True))


def select_diverse_shortlist(
    ranked: pd.DataFrame,
    per_group: int = 10,
    minimum_pairwise_distance: int = 3,
) -> pd.DataFrame:
    """Select CNS, low-liver, and balanced groups with greedy diversity."""
    required = {"AA", "passes_packaging_gate", "is_pareto", "training_distance_lower_bound"}
    missing = required - set(ranked.columns)
    if missing:
        raise ValueError(f"Missing shortlist columns: {sorted(missing)}")
    eligible = ranked.loc[
        ranked["passes_packaging_gate"]
        & ranked["training_distance_lower_bound"].ge(2)
    ].copy()
    if eligible.empty:
        return eligible.assign(
            selection_group=pd.Series(dtype=str),
            selection_rank=pd.Series(dtype=int),
            selected_from_pareto=pd.Series(dtype=bool),
        )
    quality_cutoff = float(eligible["display_score"].quantile(0.95))
    quality_pool = eligible.loc[eligible["display_score"].ge(quality_cutoff)].copy()
    selected_indices: list[int] = []
    selected_peptides: list[str] = []
    assignments: list[tuple[int, str, int, bool, float]] = []
    group_orders = {
        "cns_favoring": ["f_cns", "f_liv", "organ_uncertainty_mean"],
        "low_liver": ["f_liv", "f_cns", "organ_uncertainty_mean"],
        "balanced": [
            "display_score",
            "weight_stability_top_fraction",
            "organ_uncertainty_mean",
        ],
    }
    group_ascending = {
        "cns_favoring": [False, True, True],
        "low_liver": [True, False, True],
        "balanced": [False, False, True],
    }
    ordered_pools = {
        group_name: quality_pool.sort_values(
            columns, ascending=group_ascending[group_name], kind="stable"
        )
        for group_name, columns in group_orders.items()
    }
    positions = dict.fromkeys(group_orders, 0)
    for group_rank in range(1, per_group + 1):
        for group_name in group_orders:
            ordered = ordered_pools[group_name]
            chosen_index = None
            while positions[group_name] < len(ordered):
                index = ordered.index[positions[group_name]]
                positions[group_name] += 1
                if index in selected_indices:
                    continue
                peptide = str(ranked.loc[index, "AA"])
                if any(
                    hamming_distance(peptide, chosen) < minimum_pairwise_distance
                    for chosen in selected_peptides
                ):
                    continue
                chosen_index = index
                break
            if chosen_index is None:
                continue
            peptide = str(ranked.loc[chosen_index, "AA"])
            selected_indices.append(chosen_index)
            selected_peptides.append(peptide)
            assignments.append(
                (
                    chosen_index,
                    group_name,
                    group_rank,
                    bool(ranked.loc[chosen_index, "is_pareto"]),
                    quality_cutoff,
                )
            )
    shortlist = ranked.loc[selected_indices].copy()
    assignment_frame = pd.DataFrame(
        assignments,
        columns=[
            "source_index",
            "selection_group",
            "selection_rank",
            "selected_from_pareto",
            "selection_quality_cutoff",
        ],
    ).set_index("source_index")
    shortlist = shortlist.join(assignment_frame)
    shortlist["minimum_pairwise_hamming_required"] = minimum_pairwise_distance
    return shortlist.sort_values(["selection_group", "selection_rank"]).reset_index(drop=True)


def run_virtual_screen(
    screen: pd.DataFrame,
    reconstructed: pd.DataFrame,
    pool_size: int = 200_000,
    ensemble_size: int = 5,
    random_state: int = 42,
    max_iter: int = 80,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    """Run the complete in-silico candidate funnel."""
    observed_sequences = set(screen["AA"]) | set(reconstructed["AA"])
    peptides = generate_candidate_peptides(
        pool_size, excluded=observed_sequences, random_state=random_state
    )
    candidate_features = one_hot_7mer(peptides)
    packaging_model, packaging_threshold, lower_bound_offset, packaging_metrics = (
        fit_calibrated_packaging_model(screen, random_state=random_state)
    )
    packaging_predictions = packaging_model.predict(candidate_features)
    organ_mean, organ_std, organ_metadata = predict_organ_ensemble(
        reconstructed,
        candidate_features,
        ensemble_size=ensemble_size,
        random_state=random_state,
        max_iter=max_iter,
    )
    candidates = pd.DataFrame(
        {
            "variant_id": [f"VIRT_{index:07d}" for index in range(1, pool_size + 1)],
            "AA": peptides,
            "pred_pack": packaging_predictions,
            "pred_pack_lcb": packaging_predictions - lower_bound_offset,
        }
    )
    for endpoint_index, endpoint in enumerate(MULTIORGAN_ENDPOINTS):
        candidates[f"pred_{endpoint}_mouse"] = organ_mean[:, endpoint_index]
        candidates[f"uncertainty_{endpoint}_mouse"] = organ_std[:, endpoint_index]
    candidates["organ_uncertainty_mean"] = organ_std.mean(axis=1)
    candidates["training_distance_lower_bound"] = training_distance_lower_bound(
        peptides, observed_sequences
    )
    human_liver, human_liver_threshold = predict_human_liver_annotations(
        screen, candidate_features
    )
    candidates = pd.concat(
        [
            candidates,
            annotate_sequence_liabilities(peptides),
            human_liver,
        ],
        axis=1,
    )
    ranked = rank_candidates(candidates, packaging_threshold=packaging_threshold)
    ranked = add_weight_sensitivity(ranked)
    shortlist = select_diverse_shortlist(ranked)
    summary: dict[str, object] = {
        "random_state": random_state,
        "pool_size": pool_size,
        "sequence_space_size": len(AMINO_ACIDS) ** PEPTIDE_LENGTH,
        "sequence_space_fraction": pool_size / (len(AMINO_ACIDS) ** PEPTIDE_LENGTH),
        "observed_sequences_excluded": len(observed_sequences),
        "packaging": packaging_metrics,
        "organ_ensemble": {
            "architecture": "shared_mlp_64_32",
            "members": ensemble_size,
            **organ_metadata,
        },
        "human_liver_warning_threshold": human_liver_threshold,
        "weight_scenarios": 27,
        "packaging_gate_passed": int(ranked["passes_packaging_gate"].sum()),
        "pareto_candidates": int(ranked["is_pareto"].sum()),
        "distance_two_or_more": int(ranked["training_distance_lower_bound"].ge(2).sum()),
        "shortlist_rows": len(shortlist),
        "shortlist_groups": shortlist["selection_group"].value_counts().to_dict(),
        "shortlist_pareto_rows": int(shortlist["selected_from_pareto"].sum()),
        "shortlist_human_liver_warnings": int(shortlist["human_liver_warning"].sum()),
        "shortlist_quality_cutoff": float(shortlist["selection_quality_cutoff"].iloc[0]),
        "shortlist_minimum_pairwise_hamming": int(
            min(
                hamming_distance(left, right)
                for index, left in enumerate(shortlist["AA"])
                for right in shortlist["AA"].iloc[index + 1 :]
            )
            if len(shortlist) > 1
            else 0
        ),
        "immune_annotation": "not_scored_no_validated_7mer_neutralization_table",
        "claim_boundary": (
            "Computational mouse-organ biodistribution hypotheses; not human motor-neuron "
            "specificity, reduced hepatotoxicity, or therapeutic efficacy."
        ),
    }
    return ranked, shortlist, summary
