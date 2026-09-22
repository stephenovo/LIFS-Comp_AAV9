"""Deterministic baseline evaluation for sequence-linked Fit4Function tables."""

from __future__ import annotations

import hashlib
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

from aav9_sma.features.encode import combined_7mer_features, one_hot_7mer
from aav9_sma.models.baseline import build_regressor, build_shared_mlp
from aav9_sma.models.metrics import (
    bootstrap_confidence_intervals,
    pearson_r,
    regression_metrics,
)

MULTIORGAN_ENDPOINTS = ("brain", "spinal_cord", "liver", "heart", "kidney")

SCREEN_TASKS = (
    "Production2",
    "Liver",
    "HepG2_bind",
    "HepG2_tr",
    "THLE_bind",
    "THLE_tr",
)

# Animal 4 was inspected during model comparison and screening-policy selection.
# Keep this role explicit so downstream reports cannot silently relabel it as a
# final blind test.
DEVELOPMENT_TEST_ROLE = "development_holdout_animal4"
CROSS_ANIMAL_TEST_ROLE = "development_cross_animal_robustness"
FINAL_BLIND_TEST_ROLE = "final_blind_external"


def validate_test_role(test_role: str, *, includes_animal4: bool = False) -> None:
    """Reject an invalid claim that includes Animal 4 as a final blind test."""
    valid_roles = {DEVELOPMENT_TEST_ROLE, CROSS_ANIMAL_TEST_ROLE, FINAL_BLIND_TEST_ROLE}
    if test_role not in valid_roles:
        raise ValueError(f"Unknown test role: {test_role}")
    if test_role == FINAL_BLIND_TEST_ROLE and includes_animal4:
        raise ValueError(
            "Animal 4 was inspected during development and cannot be labeled final_blind_external"
        )


def _pearson(targets: np.ndarray, predictions: np.ndarray) -> float:
    return pearson_r(targets, predictions)


def _encode(peptides: list[str], feature_set: str) -> np.ndarray:
    if feature_set == "one_hot":
        return one_hot_7mer(peptides)
    if feature_set == "one_hot_physchem":
        return combined_7mer_features(peptides)
    raise ValueError(f"Unknown feature set: {feature_set}")


def _bootstrap_columns(
    targets: np.ndarray,
    predictions: np.ndarray,
    *,
    n_resamples: int,
    random_state: int,
    prefix: str = "model_vs_animal4",
) -> dict[str, float | int]:
    metrics = (
        bootstrap_confidence_intervals(
            targets,
            predictions,
            n_resamples=n_resamples,
            random_state=random_state,
        )
        if n_resamples
        else regression_metrics(targets, predictions)
    )
    return {f"{prefix}_{name}": value for name, value in metrics.items()}


def _other_animal_mean_enrichment(
    frame: pd.DataFrame,
    endpoint: str,
    training_animals: tuple[int, ...],
    *,
    denominator_mode: str,
    virus_round: int,
) -> np.ndarray:
    """Rebuild an organ enrichment from the mean RPM of selected animals."""
    organ_columns = [
        f"rpm_{denominator_mode}__{endpoint}_a{animal}" for animal in training_animals
    ]
    virus_column = f"rpm_{denominator_mode}__virus_prod{virus_round}"
    missing = [column for column in organ_columns + [virus_column] if column not in frame]
    if missing:
        raise ValueError(f"Missing reconstructed columns: {missing}")
    organ_rpm = frame[organ_columns].apply(pd.to_numeric, errors="coerce").mean(axis=1)
    virus_rpm = pd.to_numeric(frame[virus_column], errors="coerce")
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log2(organ_rpm.to_numpy(float) / virus_rpm.to_numpy(float))


def benchmark_screen_models(
    screen_csv: str | Path,
    model_names: tuple[str, ...] = ("ridge",),
    tasks: tuple[str, ...] = SCREEN_TASKS,
    random_state: int = 42,
    test_fraction: float = 0.2,
) -> list[dict[str, object]]:
    """Benchmark baselines on the public 100K sequence-linked screen sample.

    This deliberately uses a random holdout to mirror an early sanity check. It
    is not a substitute for the sequence-distance and animal holdouts required
    for final project claims.
    """
    frame = pd.read_csv(screen_csv)
    features = one_hot_7mer(frame["AA"].tolist())
    rows: list[dict[str, object]] = []
    for task in tasks:
        targets = pd.to_numeric(frame[task], errors="coerce").to_numpy(dtype=float)
        observed = np.isfinite(targets)
        indices = np.flatnonzero(observed)
        train_indices, test_indices = train_test_split(
            indices,
            test_size=test_fraction,
            random_state=random_state,
        )
        for model_name in model_names:
            model = build_regressor(model_name, random_state=random_state)
            model.fit(features[train_indices], targets[train_indices])
            predictions = model.predict(features[test_indices])
            rows.append(
                {
                    "task": task,
                    "model": model_name,
                    "split": "random-holdout-sanity-check",
                    "random_state": random_state,
                    "train_rows": len(train_indices),
                    "test_rows": len(test_indices),
                    "pearson_r": _pearson(targets[test_indices], predictions),
                    "r2": float(r2_score(targets[test_indices], predictions)),
                    "mae": float(mean_absolute_error(targets[test_indices], predictions)),
                }
            )
    return rows


def benchmark_production_generalization(
    modeling_csv: str | Path,
    assessment_csv: str | Path,
    model_names: tuple[str, ...] = ("ridge",),
    random_state: int = 42,
    train_rows: int = 24_000,
) -> list[dict[str, object]]:
    """Train on the modeling library and test on unique assessment variants."""
    modeling = pd.read_csv(modeling_csv, usecols=["AA", "Label", "Production"])
    assessment = pd.read_csv(assessment_csv, usecols=["AA", "Label", "Production"])
    modeling = modeling[modeling["Label"].eq("Designed")].copy()
    assessment = assessment[assessment["Label"].eq("Designed")].copy()
    modeling["target"] = np.log2(modeling["Production"].where(modeling["Production"] > 0))
    assessment["target"] = np.log2(assessment["Production"].where(assessment["Production"] > 0))
    modeling = modeling[np.isfinite(modeling["target"])].reset_index(drop=True)
    assessment = assessment[
        np.isfinite(assessment["target"]) & ~assessment["AA"].isin(set(modeling["AA"]))
    ].reset_index(drop=True)
    train = modeling.sample(n=train_rows, random_state=random_state)
    train_features = one_hot_7mer(train["AA"].tolist())
    test_features = one_hot_7mer(assessment["AA"].tolist())
    rows: list[dict[str, object]] = []
    for model_name in model_names:
        model = build_regressor(model_name, random_state=random_state)
        model.fit(train_features, train["target"].to_numpy())
        predictions = model.predict(test_features)
        targets = assessment["target"].to_numpy()
        rows.append(
            {
                "task": "Production",
                "model": model_name,
                "split": "independent-assessment-library",
                "random_state": random_state,
                "train_rows": len(train),
                "test_rows": len(assessment),
                "pearson_r": _pearson(targets, predictions),
                "r2": float(r2_score(targets, predictions)),
                "mae": float(mean_absolute_error(targets, predictions)),
            }
        )
    return rows


def sequence_distance_split(
    peptides: list[str],
    test_fraction: float = 0.2,
    random_state: int = 42,
    minimum_hamming_distance: int = 2,
) -> tuple[np.ndarray, np.ndarray]:
    """Create a deterministic holdout and remove distance-1 neighbors from training."""
    if not 0 < test_fraction < 1:
        raise ValueError("test_fraction must be between 0 and 1")
    if minimum_hamming_distance not in {1, 2}:
        raise ValueError("minimum_hamming_distance currently supports 1 or 2")
    threshold = int(test_fraction * (2**64))
    test_mask = np.array(
        [
            int.from_bytes(
                hashlib.sha256(f"{random_state}:{peptide}".encode()).digest()[:8],
                "big",
            )
            < threshold
            for peptide in peptides
        ]
    )
    train_sequences = {
        peptide for peptide, is_test in zip(peptides, test_mask, strict=True) if not is_test
    }
    if minimum_hamming_distance == 2:
        alphabet = "ACDEFGHIKLMNPQRSTVWY"
        for peptide, is_test in zip(peptides, test_mask, strict=True):
            if not is_test:
                continue
            for position, replacement in product(range(len(peptide)), alphabet):
                if replacement == peptide[position]:
                    continue
                neighbor = peptide[:position] + replacement + peptide[position + 1 :]
                train_sequences.discard(neighbor)
    train_mask = np.array([peptide in train_sequences for peptide in peptides])
    return train_mask, test_mask


def benchmark_multiorgan_animal_holdout(
    reconstructed_csv: str | Path,
    model_names: tuple[str, ...] = ("ridge",),
    endpoints: tuple[str, ...] = MULTIORGAN_ENDPOINTS,
    denominator_mode: str = "whitelist",
    virus_round: int = 2,
    random_state: int = 42,
    test_fraction: float = 0.2,
    feature_set: str = "one_hot",
    bootstrap_resamples: int = 0,
) -> list[dict[str, object]]:
    """Train on animals 1-3 and unseen sequences; evaluate against animal 4."""
    frame = pd.read_csv(reconstructed_csv)
    peptides = frame["AA"].tolist()
    train_split, test_split = sequence_distance_split(
        peptides,
        test_fraction=test_fraction,
        random_state=random_state,
        minimum_hamming_distance=2,
    )
    features = _encode(peptides, feature_set)
    rows: list[dict[str, object]] = []
    for endpoint in endpoints:
        train_column = (
            f"log2enr_{denominator_mode}__{endpoint}_animals_1_3__over__virus_prod{virus_round}"
        )
        animal4_column = f"log2enr_{denominator_mode}__{endpoint}_a4__over__virus_prod{virus_round}"
        if train_column not in frame or animal4_column not in frame:
            raise ValueError(f"Missing reconstructed columns for {endpoint}")
        train_targets = pd.to_numeric(frame[train_column], errors="coerce").to_numpy(float)
        animal4_targets = pd.to_numeric(frame[animal4_column], errors="coerce").to_numpy(float)
        train_rows = train_split & np.isfinite(train_targets)
        test_rows = test_split & np.isfinite(train_targets) & np.isfinite(animal4_targets)
        for model_name in model_names:
            model = build_regressor(model_name, random_state=random_state)
            model.fit(features[train_rows], train_targets[train_rows])
            predictions = model.predict(features[test_rows])
            held_out_targets = animal4_targets[test_rows]
            row: dict[str, object] = {
                "task": endpoint,
                "model": model_name,
                "feature_set": feature_set,
                "split": "distance-2-sequence-holdout__train-a1-a3__test-a4",
                "test_role": DEVELOPMENT_TEST_ROLE,
                "random_state": random_state,
                "train_rows": int(train_rows.sum()),
                "test_rows": int(test_rows.sum()),
                "model_vs_animal4_pearson_r": _pearson(held_out_targets, predictions),
                "model_vs_animal4_r2": float(r2_score(held_out_targets, predictions)),
                "model_vs_animal4_mae": float(mean_absolute_error(held_out_targets, predictions)),
                "animal_mean_vs_animal4_pearson_r": _pearson(
                    held_out_targets, train_targets[test_rows]
                ),
                "model_vs_animal_mean_pearson_r": _pearson(train_targets[test_rows], predictions),
            }
            if bootstrap_resamples:
                row.update(
                    _bootstrap_columns(
                        held_out_targets,
                        predictions,
                        n_resamples=bootstrap_resamples,
                        random_state=random_state + len(rows),
                    )
                )
            rows.append(row)
    return rows


def benchmark_multitask_animal_holdout(
    reconstructed_csv: str | Path,
    endpoints: tuple[str, ...] = MULTIORGAN_ENDPOINTS,
    denominator_mode: str = "whitelist",
    virus_round: int = 2,
    random_state: int = 42,
    test_fraction: float = 0.2,
    max_iter: int = 80,
) -> list[dict[str, object]]:
    """Evaluate one shared multi-output MLP on the strict animal-4 holdout.

    The network shares both hidden layers across organs and has one output per
    endpoint. Targets are standardized using training rows only so high-range
    organs do not dominate optimization. Since scikit-learn's MLP requires a
    complete target matrix, training uses rows observed for every endpoint.
    """
    frame = pd.read_csv(reconstructed_csv)
    peptides = frame["AA"].tolist()
    train_split, test_split = sequence_distance_split(
        peptides,
        test_fraction=test_fraction,
        random_state=random_state,
        minimum_hamming_distance=2,
    )
    features = one_hot_7mer(peptides)
    train_columns = [
        f"log2enr_{denominator_mode}__{endpoint}_animals_1_3__over__virus_prod{virus_round}"
        for endpoint in endpoints
    ]
    animal4_columns = [
        f"log2enr_{denominator_mode}__{endpoint}_a4__over__virus_prod{virus_round}"
        for endpoint in endpoints
    ]
    missing = [column for column in train_columns + animal4_columns if column not in frame]
    if missing:
        raise ValueError(f"Missing reconstructed columns: {missing}")

    train_targets = frame[train_columns].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    animal4_targets = frame[animal4_columns].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    complete_train_rows = train_split & np.isfinite(train_targets).all(axis=1)
    if complete_train_rows.sum() < 2:
        raise ValueError("Fewer than two complete multi-organ training rows")

    target_scaler = StandardScaler()
    scaled_targets = target_scaler.fit_transform(train_targets[complete_train_rows])
    model = build_shared_mlp(random_state=random_state, max_iter=max_iter)
    model.fit(features[complete_train_rows], scaled_targets)
    predictions = target_scaler.inverse_transform(model.predict(features))

    rows: list[dict[str, object]] = []
    for endpoint_index, endpoint in enumerate(endpoints):
        test_rows = (
            test_split
            & np.isfinite(train_targets[:, endpoint_index])
            & np.isfinite(animal4_targets[:, endpoint_index])
        )
        held_out_targets = animal4_targets[test_rows, endpoint_index]
        endpoint_predictions = predictions[test_rows, endpoint_index]
        rows.append(
            {
                "task": endpoint,
                "model": "shared_mlp_64_32",
                "split": "distance-2-sequence-holdout__train-a1-a3__test-a4",
                "test_role": DEVELOPMENT_TEST_ROLE,
                "random_state": random_state,
                "train_rows": int(complete_train_rows.sum()),
                "test_rows": int(test_rows.sum()),
                "iterations": int(model.n_iter_),
                "model_vs_animal4_pearson_r": _pearson(held_out_targets, endpoint_predictions),
                "model_vs_animal4_r2": float(r2_score(held_out_targets, endpoint_predictions)),
                "model_vs_animal4_mae": float(
                    mean_absolute_error(held_out_targets, endpoint_predictions)
                ),
                "animal_mean_vs_animal4_pearson_r": _pearson(
                    held_out_targets, train_targets[test_rows, endpoint_index]
                ),
                "model_vs_animal_mean_pearson_r": _pearson(
                    train_targets[test_rows, endpoint_index], endpoint_predictions
                ),
            }
        )
    return rows


def benchmark_multitask_ensemble_animal_holdout(
    reconstructed_csv: str | Path,
    endpoints: tuple[str, ...] = MULTIORGAN_ENDPOINTS,
    ensemble_size: int = 5,
    denominator_mode: str = "whitelist",
    virus_round: int = 2,
    random_state: int = 42,
    test_fraction: float = 0.2,
    max_iter: int = 80,
    bootstrap_resamples: int = 0,
) -> list[dict[str, object]]:
    """Evaluate the exact shared-MLP ensemble used by virtual screening."""
    if ensemble_size < 2:
        raise ValueError("ensemble_size must be at least two")
    frame = pd.read_csv(reconstructed_csv)
    peptides = frame["AA"].tolist()
    train_split, test_split = sequence_distance_split(
        peptides,
        test_fraction=test_fraction,
        random_state=random_state,
        minimum_hamming_distance=2,
    )
    features = one_hot_7mer(peptides)
    train_columns = [
        f"log2enr_{denominator_mode}__{endpoint}_animals_1_3__over__virus_prod{virus_round}"
        for endpoint in endpoints
    ]
    animal4_columns = [
        f"log2enr_{denominator_mode}__{endpoint}_a4__over__virus_prod{virus_round}"
        for endpoint in endpoints
    ]
    missing = [column for column in train_columns + animal4_columns if column not in frame]
    if missing:
        raise ValueError(f"Missing reconstructed columns: {missing}")
    train_targets = frame[train_columns].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    animal4_targets = frame[animal4_columns].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    complete_train_rows = train_split & np.isfinite(train_targets).all(axis=1)
    scaler = StandardScaler().fit(train_targets[complete_train_rows])
    scaled_targets = scaler.transform(train_targets[complete_train_rows])
    member_predictions = []
    iterations = []
    for member in range(ensemble_size):
        model = build_shared_mlp(
            random_state=random_state + member,
            max_iter=max_iter,
        )
        model.fit(features[complete_train_rows], scaled_targets)
        member_predictions.append(scaler.inverse_transform(model.predict(features)))
        iterations.append(int(model.n_iter_))
    stacked = np.stack(member_predictions)
    prediction_mean = stacked.mean(axis=0)
    prediction_std = stacked.std(axis=0)

    rows: list[dict[str, object]] = []
    for endpoint_index, endpoint in enumerate(endpoints):
        test_rows = (
            test_split
            & np.isfinite(train_targets[:, endpoint_index])
            & np.isfinite(animal4_targets[:, endpoint_index])
        )
        targets = animal4_targets[test_rows, endpoint_index]
        predictions = prediction_mean[test_rows, endpoint_index]
        disagreement = prediction_std[test_rows, endpoint_index]
        absolute_error = np.abs(targets - predictions)
        row: dict[str, object] = {
            "task": endpoint,
            "model": f"shared_mlp_64_32_ensemble_{ensemble_size}",
            "split": "distance-2-sequence-holdout__train-a1-a3__test-a4",
            "test_role": DEVELOPMENT_TEST_ROLE,
            "random_state": random_state,
            "train_rows": int(complete_train_rows.sum()),
            "test_rows": int(test_rows.sum()),
            "ensemble_size": ensemble_size,
            "member_iterations": "|".join(map(str, iterations)),
            "model_vs_animal4_pearson_r": _pearson(targets, predictions),
            "model_vs_animal4_r2": float(r2_score(targets, predictions)),
            "model_vs_animal4_mae": float(mean_absolute_error(targets, predictions)),
            "animal_mean_vs_animal4_pearson_r": _pearson(
                targets, train_targets[test_rows, endpoint_index]
            ),
            "model_vs_animal_mean_pearson_r": _pearson(
                train_targets[test_rows, endpoint_index], predictions
            ),
            "disagreement_vs_absolute_error_pearson_r": _pearson(absolute_error, disagreement),
        }
        if bootstrap_resamples:
            row.update(
                _bootstrap_columns(
                    targets,
                    predictions,
                    n_resamples=bootstrap_resamples,
                    random_state=random_state + endpoint_index,
                )
            )
        rows.append(row)
    return rows


def benchmark_multitask_ensemble_leave_one_animal_out(
    reconstructed_csv: str | Path,
    endpoints: tuple[str, ...] = MULTIORGAN_ENDPOINTS,
    animals: tuple[int, ...] = (1, 2, 3, 4),
    ensemble_size: int = 5,
    denominator_mode: str = "whitelist",
    virus_round: int = 2,
    random_state: int = 42,
    test_fraction: float = 0.2,
    max_iter: int = 80,
    bootstrap_resamples: int = 0,
) -> list[dict[str, object]]:
    """Run a retrospective leave-one-animal-out robustness audit.

    Every fold also keeps the distance-2 sequence holdout used by the formal
    Animal 4 development benchmark. This is stronger than rerunning Animal 4,
    but it remains retrospective development evidence rather than a new blind
    animal because all four animals' labels were already available to the team.
    """
    if ensemble_size < 2:
        raise ValueError("ensemble_size must be at least two")
    if len(animals) < 3 or len(set(animals)) != len(animals):
        raise ValueError("animals must contain at least three distinct identifiers")

    frame = pd.read_csv(reconstructed_csv)
    peptides = frame["AA"].tolist()
    train_split, test_split = sequence_distance_split(
        peptides,
        test_fraction=test_fraction,
        random_state=random_state,
        minimum_hamming_distance=2,
    )
    features = one_hot_7mer(peptides)
    rows: list[dict[str, object]] = []

    for held_out_animal in animals:
        training_animals = tuple(animal for animal in animals if animal != held_out_animal)
        train_targets = np.column_stack(
            [
                _other_animal_mean_enrichment(
                    frame,
                    endpoint,
                    training_animals,
                    denominator_mode=denominator_mode,
                    virus_round=virus_round,
                )
                for endpoint in endpoints
            ]
        )
        held_out_columns = [
            (
                f"log2enr_{denominator_mode}__{endpoint}_a{held_out_animal}"
                f"__over__virus_prod{virus_round}"
            )
            for endpoint in endpoints
        ]
        missing = [column for column in held_out_columns if column not in frame]
        if missing:
            raise ValueError(f"Missing reconstructed columns: {missing}")
        held_out_targets = (
            frame[held_out_columns].apply(pd.to_numeric, errors="coerce").to_numpy(float)
        )
        complete_train_rows = train_split & np.isfinite(train_targets).all(axis=1)
        if complete_train_rows.sum() < 2:
            raise ValueError(
                f"Fewer than two complete training rows for held-out animal {held_out_animal}"
            )

        scaler = StandardScaler().fit(train_targets[complete_train_rows])
        scaled_targets = scaler.transform(train_targets[complete_train_rows])
        member_predictions = []
        iterations = []
        for member in range(ensemble_size):
            model = build_shared_mlp(
                random_state=random_state + member,
                max_iter=max_iter,
            )
            model.fit(features[complete_train_rows], scaled_targets)
            member_predictions.append(scaler.inverse_transform(model.predict(features)))
            iterations.append(int(model.n_iter_))
        stacked = np.stack(member_predictions)
        prediction_mean = stacked.mean(axis=0)
        prediction_std = stacked.std(axis=0)

        for endpoint_index, endpoint in enumerate(endpoints):
            test_rows = (
                test_split
                & np.isfinite(train_targets[:, endpoint_index])
                & np.isfinite(held_out_targets[:, endpoint_index])
            )
            targets = held_out_targets[test_rows, endpoint_index]
            predictions = prediction_mean[test_rows, endpoint_index]
            disagreement = prediction_std[test_rows, endpoint_index]
            absolute_error = np.abs(targets - predictions)
            point_metrics = regression_metrics(targets, predictions)
            row: dict[str, object] = {
                "task": endpoint,
                "model": f"shared_mlp_64_32_ensemble_{ensemble_size}",
                "split": (
                    "distance-2-sequence-holdout__leave-one-animal-out__"
                    f"test-a{held_out_animal}"
                ),
                "test_role": CROSS_ANIMAL_TEST_ROLE,
                "held_out_animal": held_out_animal,
                "training_animals": "|".join(map(str, training_animals)),
                "random_state": random_state,
                "train_rows": int(complete_train_rows.sum()),
                "test_rows": int(test_rows.sum()),
                "ensemble_size": ensemble_size,
                "member_iterations": "|".join(map(str, iterations)),
                **{
                    f"model_vs_heldout_{metric}": value
                    for metric, value in point_metrics.items()
                },
                "training_mean_vs_heldout_pearson_r": _pearson(
                    targets, train_targets[test_rows, endpoint_index]
                ),
                "model_vs_training_mean_pearson_r": _pearson(
                    train_targets[test_rows, endpoint_index], predictions
                ),
                "disagreement_vs_absolute_error_pearson_r": _pearson(
                    absolute_error, disagreement
                ),
            }
            if bootstrap_resamples:
                row.update(
                    _bootstrap_columns(
                        targets,
                        predictions,
                        n_resamples=bootstrap_resamples,
                        random_state=(
                            random_state + held_out_animal * 100 + endpoint_index
                        ),
                        prefix="model_vs_heldout",
                    )
                )
            rows.append(row)
    return rows


def benchmark_masked_multitask_animal_holdout(
    reconstructed_csv: str | Path,
    endpoints: tuple[str, ...] = MULTIORGAN_ENDPOINTS,
    denominator_mode: str = "whitelist",
    virus_round: int = 2,
    random_state: int = 42,
    test_fraction: float = 0.2,
    validation_fraction: float = 0.15,
    max_epochs: int = 80,
    bootstrap_resamples: int = 500,
    device: str = "cpu",
) -> list[dict[str, object]]:
    """Evaluate a shared PyTorch MLP that masks missing labels per endpoint.

    Epoch count is chosen using a distance-separated validation subset drawn
    only from the outer training partition. The final model is refit on the
    union of inner-train and validation rows before Animal 4 evaluation.
    """
    from aav9_sma.models.multitask_torch import fit_masked_multitask_mlp

    frame = pd.read_csv(reconstructed_csv)
    peptides = frame["AA"].tolist()
    outer_train, outer_test = sequence_distance_split(
        peptides,
        test_fraction=test_fraction,
        random_state=random_state,
        minimum_hamming_distance=2,
    )
    outer_indices = np.flatnonzero(outer_train)
    inner_train_local, validation_local = sequence_distance_split(
        [peptides[index] for index in outer_indices],
        test_fraction=validation_fraction,
        random_state=random_state + 1009,
        minimum_hamming_distance=2,
    )
    inner_train = np.zeros(len(frame), dtype=bool)
    validation = np.zeros(len(frame), dtype=bool)
    inner_train[outer_indices] = inner_train_local
    validation[outer_indices] = validation_local

    train_columns = [
        f"log2enr_{denominator_mode}__{endpoint}_animals_1_3__over__virus_prod{virus_round}"
        for endpoint in endpoints
    ]
    animal4_columns = [
        f"log2enr_{denominator_mode}__{endpoint}_a4__over__virus_prod{virus_round}"
        for endpoint in endpoints
    ]
    missing = [column for column in train_columns + animal4_columns if column not in frame]
    if missing:
        raise ValueError(f"Missing reconstructed columns: {missing}")
    train_targets = frame[train_columns].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    animal4_targets = frame[animal4_columns].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    features = one_hot_7mer(peptides)
    result = fit_masked_multitask_mlp(
        features,
        train_targets,
        inner_train,
        validation,
        max_epochs=max_epochs,
        random_state=random_state,
        device=device,
    )

    rows: list[dict[str, object]] = []
    for endpoint_index, endpoint in enumerate(endpoints):
        test_rows = (
            outer_test
            & np.isfinite(train_targets[:, endpoint_index])
            & np.isfinite(animal4_targets[:, endpoint_index])
        )
        targets = animal4_targets[test_rows, endpoint_index]
        predictions = result.predictions[test_rows, endpoint_index]
        row: dict[str, object] = {
            "task": endpoint,
            "model": "torch_masked_shared_mlp_64_32",
            "feature_set": "one_hot",
            "split": "distance-2-sequence-holdout__train-a1-a3__test-a4",
            "test_role": DEVELOPMENT_TEST_ROLE,
            "random_state": random_state,
            "train_rows": result.training_rows,
            "task_train_observations": result.task_observations[endpoint_index],
            "internal_validation_rows": result.validation_rows,
            "test_rows": int(test_rows.sum()),
            "best_epoch": result.best_epoch,
            "internal_validation_loss": result.validation_loss,
            "device": result.device,
            "animal_mean_vs_animal4_pearson_r": _pearson(
                targets, train_targets[test_rows, endpoint_index]
            ),
            "model_vs_animal_mean_pearson_r": _pearson(
                train_targets[test_rows, endpoint_index], predictions
            ),
        }
        row.update(
            _bootstrap_columns(
                targets,
                predictions,
                n_resamples=bootstrap_resamples,
                random_state=random_state + endpoint_index,
            )
        )
        rows.append(row)
    return rows
