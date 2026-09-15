"""Deterministic baseline evaluation for sequence-linked Fit4Function tables."""

from __future__ import annotations

import hashlib
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import train_test_split

from aav9_sma.features.encode import one_hot_7mer
from aav9_sma.models.baseline import build_regressor

MULTIORGAN_ENDPOINTS = ("brain", "spinal_cord", "liver", "heart", "kidney")

SCREEN_TASKS = (
    "Production2",
    "Liver",
    "HepG2_bind",
    "HepG2_tr",
    "THLE_bind",
    "THLE_tr",
)


def _pearson(targets: np.ndarray, predictions: np.ndarray) -> float:
    if np.std(targets) == 0 or np.std(predictions) == 0:
        return float("nan")
    return float(np.corrcoef(targets, predictions)[0, 1])


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
    features = one_hot_7mer(peptides)
    rows: list[dict[str, object]] = []
    for endpoint in endpoints:
        train_column = (
            f"log2enr_{denominator_mode}__{endpoint}_animals_1_3"
            f"__over__virus_prod{virus_round}"
        )
        animal4_column = (
            f"log2enr_{denominator_mode}__{endpoint}_a4__over__virus_prod{virus_round}"
        )
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
            rows.append(
                {
                    "task": endpoint,
                    "model": model_name,
                    "split": "distance-2-sequence-holdout__train-a1-a3__test-a4",
                    "random_state": random_state,
                    "train_rows": int(train_rows.sum()),
                    "test_rows": int(test_rows.sum()),
                    "model_vs_animal4_pearson_r": _pearson(held_out_targets, predictions),
                    "model_vs_animal4_r2": float(r2_score(held_out_targets, predictions)),
                    "model_vs_animal4_mae": float(
                        mean_absolute_error(held_out_targets, predictions)
                    ),
                    "animal_mean_vs_animal4_pearson_r": _pearson(
                        held_out_targets, train_targets[test_rows]
                    ),
                    "model_vs_animal_mean_pearson_r": _pearson(
                        train_targets[test_rows], predictions
                    ),
                }
            )
    return rows
