"""Sequence-linked motor-neuron transduction head.

This module intentionally requires cell-type-resolved evidence. It does not
convert bulk spinal-cord labels into motor-neuron labels.
"""

from __future__ import annotations

import pickle
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.preprocessing import OneHotEncoder

from aav9_sma.data.audit import is_valid_peptide
from aav9_sma.data.motor_neuron import (
    EVIDENCE_RANK,
    normalize_motor_neuron_evidence,
    study_level_split,
)
from aav9_sma.features.encode import one_hot_7mer
from aav9_sma.models.metrics import regression_metrics

_CATEGORICAL_CONTEXT = ("species", "strain", "route", "payload", "promoter_or_enhancer")
_NUMERIC_CONTEXT = ("dose", "timepoint")


def _validate_comparable_readouts(frame: pd.DataFrame) -> None:
    """Reject raw assay scales that the current head cannot compare safely."""
    readout_types = sorted(frame["readout_type"].str.lower().dropna().unique().tolist())
    if len(readout_types) != 1:
        raise ValueError(
            "Motor-neuron regression requires one comparable readout_type per fit; "
            f"found {readout_types}. Normalize to a common target or fit separate heads."
        )
    if "readout_unit" in frame:
        units = sorted(frame["readout_unit"].astype("string").dropna().unique().tolist())
        if len(units) > 1:
            raise ValueError(
                "Motor-neuron regression requires one readout_unit per fit; "
                f"found {units}"
            )


def _model_for_name(name: str, random_state: int):
    if name == "ridge":
        return Ridge(alpha=1.0)
    if name == "random_forest":
        return RandomForestRegressor(
            n_estimators=300,
            min_samples_leaf=2,
            n_jobs=-1,
            random_state=random_state,
        )
    raise ValueError(f"Unknown motor-neuron model: {name}")


@dataclass(frozen=True)
class MotorNeuronFitResult:
    """Fitted head and study-level validation metadata."""

    model: MotorNeuronHead
    validation_metrics: dict[str, float]
    training_rows: int
    validation_rows: int
    held_out_studies: tuple[str, ...]


class MotorNeuronHead:
    """Regress cell-resolved motor-neuron readouts from capsid and context."""

    def __init__(self, model_name: str = "ridge", random_state: int = 42) -> None:
        self.model_name = model_name
        self.random_state = random_state
        self._model = None
        self._context_encoder: OneHotEncoder | None = None
        self._numeric_means: np.ndarray | None = None
        self._numeric_scales: np.ndarray | None = None

    def _context_values(self, frame: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
        categorical = pd.DataFrame(index=frame.index)
        for column in _CATEGORICAL_CONTEXT:
            if column in frame:
                values = frame[column]
            else:
                values = pd.Series("<unknown>", index=frame.index)
            categorical[column] = values.astype("string").fillna("<unknown>").astype(str)

        numeric = np.zeros((len(frame), len(_NUMERIC_CONTEXT)), dtype=np.float32)
        for index, column in enumerate(_NUMERIC_CONTEXT):
            if column in frame:
                numeric[:, index] = pd.to_numeric(frame[column], errors="coerce").to_numpy(
                    dtype=float
                )
        if self._numeric_means is None or self._numeric_scales is None:
            raise RuntimeError("Context transformer has not been fitted")
        for index in range(numeric.shape[1]):
            missing = ~np.isfinite(numeric[:, index])
            numeric[missing, index] = self._numeric_means[index]
            numeric[:, index] = (
                numeric[:, index] - self._numeric_means[index]
            ) / self._numeric_scales[index]
        return categorical, numeric

    def _fit_features(self, frame: pd.DataFrame) -> np.ndarray:
        categorical = pd.DataFrame(index=frame.index)
        for column in _CATEGORICAL_CONTEXT:
            if column in frame:
                values = frame[column]
            else:
                values = pd.Series("<unknown>", index=frame.index)
            categorical[column] = values.astype("string").fillna("<unknown>").astype(str)
        self._context_encoder = OneHotEncoder(handle_unknown="ignore", sparse_output=False)
        context_categorical = self._context_encoder.fit_transform(categorical)

        numeric = np.zeros((len(frame), len(_NUMERIC_CONTEXT)), dtype=np.float32)
        for index, column in enumerate(_NUMERIC_CONTEXT):
            if column in frame:
                numeric[:, index] = pd.to_numeric(frame[column], errors="coerce").to_numpy(
                    dtype=float
                )
        means = np.nanmedian(np.where(np.isfinite(numeric), numeric, np.nan), axis=0)
        means = np.where(np.isfinite(means), means, 0.0).astype(np.float32)
        scales = np.nanstd(np.where(np.isfinite(numeric), numeric, np.nan), axis=0)
        scales = np.where(np.isfinite(scales) & (scales > 0), scales, 1.0).astype(np.float32)
        self._numeric_means = means
        self._numeric_scales = scales
        for index in range(numeric.shape[1]):
            missing = ~np.isfinite(numeric[:, index])
            numeric[missing, index] = means[index]
            numeric[:, index] = (numeric[:, index] - means[index]) / scales[index]
        sequences = one_hot_7mer(frame["peptide_7mer"].tolist())
        return np.concatenate([sequences, context_categorical, numeric], axis=1)

    def _transform_features(self, frame: pd.DataFrame) -> np.ndarray:
        if self._context_encoder is None:
            raise RuntimeError("Model has not been fitted")
        categorical, numeric = self._context_values(frame)
        context_categorical = self._context_encoder.transform(categorical)
        sequences = one_hot_7mer(frame["peptide_7mer"].tolist())
        return np.concatenate([sequences, context_categorical, numeric], axis=1)

    def fit(
        self,
        frame: pd.DataFrame,
        *,
        target_column: str = "readout_value",
    ) -> MotorNeuronHead:
        """Fit on rows already selected as motor-neuron evidence."""
        normalized = normalize_motor_neuron_evidence(
            frame,
            require_sequence=True,
            require_motor_labels=True,
        )
        normalized = normalized[normalized["is_motor_neuron"]].copy()
        _validate_comparable_readouts(normalized)
        if target_column not in normalized:
            raise ValueError(f"Missing target column: {target_column}")
        if len(normalized) < 2 or normalized["variant_id"].nunique() < 2:
            raise ValueError("At least two sequence-linked variants are required")
        features = self._fit_features(normalized)
        targets = pd.to_numeric(normalized[target_column], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(targets).all():
            raise ValueError(f"Target column {target_column!r} must be finite")
        self._model = _model_for_name(self.model_name, self.random_state)
        self._model.fit(features, targets)
        return self

    def predict(self, frame: pd.DataFrame) -> np.ndarray:
        """Predict a motor-neuron readout for sequence-linked contexts."""
        if self._model is None:
            raise RuntimeError("Model has not been fitted")
        if "peptide_7mer" not in frame:
            raise ValueError("Prediction requires peptide_7mer")
        prediction_frame = frame.copy()
        prediction_frame["peptide_7mer"] = (
            prediction_frame["peptide_7mer"].astype("string").str.strip().str.upper()
        )
        invalid = ~prediction_frame["peptide_7mer"].map(is_valid_peptide)
        if invalid.any():
            raise ValueError(
                f"Invalid peptide_7mer values at rows {prediction_frame.index[invalid].tolist()}"
            )
        return np.asarray(
            self._model.predict(self._transform_features(prediction_frame)),
            dtype=float,
        )

    def save(self, path: str | Path) -> None:
        """Persist the fitted head without adding a new dependency."""
        if self._model is None:
            raise RuntimeError("Cannot save an unfitted model")
        with Path(path).open("wb") as handle:
            pickle.dump(self, handle)

    @classmethod
    def load(cls, path: str | Path) -> MotorNeuronHead:
        with Path(path).open("rb") as handle:
            model = pickle.load(handle)
        if not isinstance(model, cls):
            raise TypeError("Serialized object is not a MotorNeuronHead")
        return model


def fit_motor_neuron_head(
    frame: pd.DataFrame,
    *,
    model_name: str = "ridge",
    target_column: str = "readout_value",
    min_evidence_level: str = "L3",
    validation_fraction: float = 0.2,
    random_state: int = 42,
) -> MotorNeuronFitResult:
    """Fit a head and evaluate it on studies withheld as complete groups.

    The returned model is refit on all eligible motor-neuron rows after
    validation. This makes it suitable for later screening while retaining the
    held-out metrics in the result object.
    """
    if min_evidence_level not in EVIDENCE_RANK:
        raise ValueError(f"Unknown minimum evidence level: {min_evidence_level}")
    normalized = normalize_motor_neuron_evidence(
        frame,
        require_sequence=True,
        require_motor_labels=True,
    )
    eligible = normalized[
        normalized["is_motor_neuron"]
        & normalized["evidence_level"].map(EVIDENCE_RANK.__getitem__).ge(
            EVIDENCE_RANK[min_evidence_level]
        )
    ].copy()
    if len(eligible) < 4 or eligible["variant_id"].nunique() < 2:
        raise ValueError("Too few direct motor-neuron rows for a sequence-linked head")
    _validate_comparable_readouts(eligible)
    train_rows, validation_rows = study_level_split(
        eligible,
        validation_fraction=validation_fraction,
        random_state=random_state,
    )
    validation_model = MotorNeuronHead(model_name=model_name, random_state=random_state)
    validation_model.fit(eligible.loc[train_rows], target_column=target_column)
    validation_predictions = validation_model.predict(eligible.loc[validation_rows])
    validation_targets = eligible.loc[validation_rows, target_column].to_numpy(dtype=float)
    validation_metrics = (
        regression_metrics(validation_targets, validation_predictions)
        if len(validation_targets) >= 2
        else {}
    )

    final_model = MotorNeuronHead(model_name=model_name, random_state=random_state)
    final_model.fit(eligible, target_column=target_column)
    held_out_studies = tuple(
        sorted(eligible.loc[validation_rows, "study_id"].astype(str).unique().tolist())
    )
    return MotorNeuronFitResult(
        model=final_model,
        validation_metrics=validation_metrics,
        training_rows=int(train_rows.sum()),
        validation_rows=int(validation_rows.sum()),
        held_out_studies=held_out_studies,
    )
