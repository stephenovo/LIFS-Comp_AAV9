"""Regression metrics and paired bootstrap confidence intervals."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score


def pearson_r(targets: np.ndarray, predictions: np.ndarray) -> float:
    """Return Pearson correlation, or NaN for a constant vector."""
    if np.std(targets) == 0 or np.std(predictions) == 0:
        return float("nan")
    return float(np.corrcoef(targets, predictions)[0, 1])


def regression_metrics(targets: np.ndarray, predictions: np.ndarray) -> dict[str, float]:
    """Return the common point metrics used across every model comparison."""
    targets = np.asarray(targets, dtype=float)
    predictions = np.asarray(predictions, dtype=float)
    if targets.shape != predictions.shape or targets.ndim != 1:
        raise ValueError("targets and predictions must be equal-length one-dimensional arrays")
    if not np.isfinite(targets).all() or not np.isfinite(predictions).all():
        raise ValueError("metrics require finite targets and predictions")
    return {
        "pearson_r": pearson_r(targets, predictions),
        "spearman_r": float(spearmanr(targets, predictions).statistic),
        "r2": float(r2_score(targets, predictions)),
        "mae": float(mean_absolute_error(targets, predictions)),
        "rmse": float(mean_squared_error(targets, predictions) ** 0.5),
    }


def bootstrap_confidence_intervals(
    targets: np.ndarray,
    predictions: np.ndarray,
    n_resamples: int = 500,
    confidence: float = 0.95,
    random_state: int = 42,
) -> dict[str, float | int]:
    """Estimate paired row-bootstrap intervals for all regression metrics."""
    if n_resamples < 1:
        raise ValueError("n_resamples must be positive")
    if not 0 < confidence < 1:
        raise ValueError("confidence must be between zero and one")
    targets = np.asarray(targets, dtype=float)
    predictions = np.asarray(predictions, dtype=float)
    point = regression_metrics(targets, predictions)
    metric_functions: dict[str, Callable[[np.ndarray, np.ndarray], float]] = {
        "pearson_r": pearson_r,
        "spearman_r": lambda y, p: float(spearmanr(y, p).statistic),
        "r2": lambda y, p: float(r2_score(y, p)),
        "mae": lambda y, p: float(mean_absolute_error(y, p)),
        "rmse": lambda y, p: float(mean_squared_error(y, p) ** 0.5),
    }
    rng = np.random.default_rng(random_state)
    samples = {name: np.empty(n_resamples, dtype=float) for name in metric_functions}
    for bootstrap_index in range(n_resamples):
        indices = rng.integers(0, len(targets), size=len(targets))
        sample_targets = targets[indices]
        sample_predictions = predictions[indices]
        for name, function in metric_functions.items():
            samples[name][bootstrap_index] = function(sample_targets, sample_predictions)
    alpha = (1.0 - confidence) / 2.0
    output: dict[str, float | int] = {
        **point,
        "bootstrap_resamples": n_resamples,
        "bootstrap_confidence": confidence,
    }
    for name, values in samples.items():
        finite = values[np.isfinite(values)]
        output[f"{name}_ci_low"] = float(np.quantile(finite, alpha))
        output[f"{name}_ci_high"] = float(np.quantile(finite, 1.0 - alpha))
    return output
