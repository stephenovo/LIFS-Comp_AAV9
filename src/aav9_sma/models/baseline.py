"""Small, interpretable baselines used before neural-network experiments."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
from sklearn.base import RegressorMixin
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.neural_network import MLPRegressor


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
    if name == "lightgbm":
        try:
            from lightgbm import LGBMRegressor
        except ImportError as error:
            raise ImportError(
                "LightGBM is optional; install the project with the 'models' extra"
            ) from error
        return LGBMRegressor(
            objective="regression",
            n_estimators=500,
            learning_rate=0.04,
            num_leaves=31,
            max_depth=-1,
            min_child_samples=40,
            subsample=0.9,
            colsample_bytree=0.9,
            reg_lambda=1.0,
            random_state=random_state,
            n_jobs=-1,
            verbosity=-1,
        )
    raise ValueError(f"Unknown baseline model: {name}")


def build_shared_mlp(
    random_state: int = 42,
    hidden_layer_sizes: tuple[int, ...] = (64, 32),
    max_iter: int = 80,
) -> MLPRegressor:
    """Create a shared-encoder, multi-output MLP for the five organ tasks."""
    return MLPRegressor(
        hidden_layer_sizes=hidden_layer_sizes,
        activation="relu",
        solver="adam",
        batch_size=512,
        learning_rate_init=1.0e-3,
        early_stopping=True,
        validation_fraction=0.1,
        n_iter_no_change=8,
        max_iter=max_iter,
        random_state=random_state,
    )


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
