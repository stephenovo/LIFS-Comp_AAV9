"""Small, interpretable baselines used before neural-network experiments."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
from sklearn.base import RegressorMixin
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge


def build_regressor(name: str, random_state: int = 42) -> RegressorMixin:
    """Create a supported single-task regression baseline."""
    if name == "ridge":
        return Ridge(alpha=1.0)
    if name == "random_forest":
        return RandomForestRegressor(
            n_estimators=300,
            min_samples_leaf=2,
            n_jobs=-1,
            random_state=random_state,
        )
    raise ValueError(f"Unknown baseline model: {name}")


def fit_task_models(
    features: np.ndarray,
    targets: Mapping[str, np.ndarray],
    model_name: str = "ridge",
    random_state: int = 42,
) -> dict[str, RegressorMixin]:
    """Fit one model per task while masking missing labels independently."""
    models: dict[str, RegressorMixin] = {}
    for task_name, values in targets.items():
        values = np.asarray(values, dtype=float)
        observed = np.isfinite(values)
        if observed.sum() < 2:
            raise ValueError(f"Task {task_name!r} has fewer than two observed labels")
        model = build_regressor(model_name, random_state=random_state)
        model.fit(features[observed], values[observed])
        models[task_name] = model
    return models


def predict_task_models(
    models: Mapping[str, RegressorMixin], features: np.ndarray
) -> dict[str, np.ndarray]:
    """Predict all independently fitted task models."""
    return {task_name: model.predict(features) for task_name, model in models.items()}

