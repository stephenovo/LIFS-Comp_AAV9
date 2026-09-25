"""Canonical data handling for cell-type-resolved motor-neuron evidence."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass

import numpy as np
import pandas as pd

from aav9_sma.data.audit import is_valid_peptide

REQUIRED_EVIDENCE_COLUMNS = (
    "study_id",
    "species",
    "route",
    "tissue",
    "cell_type",
    "readout_type",
    "readout_value",
    "evidence_level",
)
SEQUENCE_LINKAGE_COLUMNS = ("variant_id", "peptide_7mer")
EVIDENCE_LEVELS = ("L1", "L2", "L3", "L4")
EVIDENCE_RANK = {level: rank for rank, level in enumerate(EVIDENCE_LEVELS, start=1)}

_COLUMN_ALIASES = {
    "study": "study_id",
    "animal": "animal_id",
    "variant": "variant_id",
    "AA": "peptide_7mer",
    "aa": "peptide_7mer",
    "capsid_7mer": "peptide_7mer",
    "cell_type_label": "cell_type",
    "readout": "readout_value",
    "readout_type_label": "readout_type",
    "evidence": "evidence_level",
}

_MOTOR_NEURON_LABELS = frozenset(
    {
        "motorneuron",
        "motoneuron",
        "spinalmotorneuron",
        "spinalmotoneuron",
        "smn",
        "alphamotorneuron",
        "betamotorneuron",
        "gammamotorneuron",
    }
)

_EVIDENCE_SCORE = {"L1": 0.25, "L2": 0.50, "L3": 0.80, "L4": 1.00}
_READOUT_SCORE = {
    "dna": 0.40,
    "vector_genome": 0.40,
    "rna": 0.90,
    "protein": 0.90,
    "reporter": 0.90,
    "vector_barcode": 1.00,
    "function": 1.00,
}


@dataclass(frozen=True)
class MotorNeuronEvidenceAudit:
    """Summary of a normalized cell-type evidence table."""

    row_count: int
    motor_neuron_rows: int
    sequence_linked_rows: int
    direct_cell_level_rows: int
    studies: tuple[str, ...]
    invalid_evidence_rows: int


def _normalize_cell_label(value: object) -> str:
    text = str(value).strip().lower()
    return re.sub(r"[^a-z0-9]+", "", text)


def is_motor_neuron_label(value: object) -> bool:
    """Return whether a cell-type label identifies a motor neuron."""
    normalized = _normalize_cell_label(value)
    return normalized in _MOTOR_NEURON_LABELS or (
        "motorneuron" in normalized or "motoneuron" in normalized
    )


def _canonicalize_columns(frame: pd.DataFrame) -> pd.DataFrame:
    output = frame.copy()
    renames: dict[str, str] = {}
    for source, target in _COLUMN_ALIASES.items():
        if source in output.columns and target not in output.columns:
            renames[source] = target
    return output.rename(columns=renames)


def normalize_motor_neuron_evidence(
    frame: pd.DataFrame,
    *,
    require_sequence: bool = False,
    require_motor_labels: bool = False,
) -> pd.DataFrame:
    """Validate and normalize the cell-type evidence contract.

    The input is never modified. L1/L2 rows are retained for provenance, but
    direct cell-type modeling should explicitly filter to L3/L4 rows.
    """
    if not isinstance(frame, pd.DataFrame):
        raise TypeError("frame must be a pandas DataFrame")
    output = _canonicalize_columns(frame)
    missing = [column for column in REQUIRED_EVIDENCE_COLUMNS if column not in output]
    if missing:
        raise ValueError(f"Missing motor-neuron evidence columns: {missing}")
    if require_sequence:
        missing_sequence = [column for column in SEQUENCE_LINKAGE_COLUMNS if column not in output]
        if missing_sequence:
            raise ValueError(f"Sequence-linked modeling requires columns: {missing_sequence}")

    for column in REQUIRED_EVIDENCE_COLUMNS:
        output[column] = output[column].astype("string").str.strip()
    output["readout_value"] = pd.to_numeric(output["readout_value"], errors="coerce")
    output["evidence_level"] = output["evidence_level"].str.upper()
    invalid_levels = ~output["evidence_level"].isin(EVIDENCE_LEVELS)
    if invalid_levels.any():
        bad = sorted(output.loc[invalid_levels, "evidence_level"].dropna().unique().tolist())
        raise ValueError(f"Unknown evidence levels: {bad}; expected {EVIDENCE_LEVELS}")
    if output["readout_value"].isna().any():
        raise ValueError("readout_value must be numeric and finite")
    if not np.isfinite(output["readout_value"].to_numpy(dtype=float)).all():
        raise ValueError("readout_value must be numeric and finite")

    output["cell_type_normalized"] = output["cell_type"].map(_normalize_cell_label)
    output["is_motor_neuron"] = output["cell_type"].map(is_motor_neuron_label)
    output["direct_cell_level_evidence"] = output["evidence_level"].isin(("L3", "L4"))
    if require_motor_labels and not output["is_motor_neuron"].any():
        raise ValueError("No motor-neuron cell labels were found")
    if require_sequence:
        output["variant_id"] = output["variant_id"].astype("string").str.strip()
        output["peptide_7mer"] = output["peptide_7mer"].astype("string").str.strip().str.upper()
        if output["variant_id"].isna().any() or (output["variant_id"] == "").any():
            raise ValueError("variant_id must be present for sequence-linked modeling")
        invalid_peptides = ~output["peptide_7mer"].map(is_valid_peptide)
        if invalid_peptides.any():
            raise ValueError(
                f"Invalid peptide_7mer values at rows {output.index[invalid_peptides].tolist()}"
            )
        sequence_counts = output.groupby("variant_id")["peptide_7mer"].nunique()
        conflicting_variants = sequence_counts[sequence_counts > 1].index.tolist()
        if conflicting_variants:
            raise ValueError(
                "Each variant_id must map to exactly one peptide_7mer; conflicts: "
                f"{conflicting_variants}"
            )
        output["sequence_linked"] = True
    else:
        output["sequence_linked"] = (
            output.get("variant_id", pd.Series(index=output.index, dtype="string")).notna()
            & output.get("peptide_7mer", pd.Series(index=output.index, dtype="string")).notna()
        )
    return output.reset_index(drop=True)


def audit_motor_neuron_evidence(frame: pd.DataFrame) -> MotorNeuronEvidenceAudit:
    """Return a provenance-focused audit without changing the source table."""
    normalized = normalize_motor_neuron_evidence(frame)
    studies = tuple(sorted(normalized["study_id"].dropna().unique().tolist()))
    return MotorNeuronEvidenceAudit(
        row_count=len(normalized),
        motor_neuron_rows=int(normalized["is_motor_neuron"].sum()),
        sequence_linked_rows=int(normalized["sequence_linked"].sum()),
        direct_cell_level_rows=int(normalized["direct_cell_level_evidence"].sum()),
        studies=studies,
        invalid_evidence_rows=0,
    )


def aggregate_motor_neuron_evidence(frame: pd.DataFrame) -> pd.DataFrame:
    """Aggregate traceable motor-neuron evidence by variant.

    This is an evidence summary, not a biological prediction. A variant with
    no L3/L4 motor-neuron row is never assigned a direct cell-level score.
    """
    normalized = normalize_motor_neuron_evidence(frame, require_sequence=True)
    motor = normalized[normalized["is_motor_neuron"]].copy()
    if motor.empty:
        raise ValueError("Cannot aggregate evidence without motor-neuron rows")
    motor["evidence_score"] = motor["evidence_level"].map(_EVIDENCE_SCORE)
    readout_key = motor["readout_type"].str.lower().str.replace("-", "_", regex=False)
    motor["evidence_score"] *= readout_key.map(_READOUT_SCORE).fillna(0.5)
    grouped = (
        motor.groupby("variant_id", as_index=False)
        .agg(
            peptide_7mer=("peptide_7mer", "first"),
            motor_neuron_evidence_score=("evidence_score", "max"),
            evidence_rows=("variant_id", "size"),
            direct_cell_level_evidence=("direct_cell_level_evidence", "any"),
            max_evidence_level=(
                "evidence_level",
                lambda values: max(values, key=EVIDENCE_RANK.get),
            ),
            study_count=("study_id", "nunique"),
        )
        .reset_index(drop=True)
    )
    grouped["motor_neuron_evidence_score"] = grouped[
        "motor_neuron_evidence_score"
    ].clip(lower=0.0, upper=1.0)
    return grouped


def study_level_split(
    frame: pd.DataFrame,
    *,
    group_column: str = "study_id",
    validation_fraction: float = 0.2,
    random_state: int = 42,
) -> tuple[np.ndarray, np.ndarray]:
    """Create disjoint row masks whose groups never cross the split."""
    if not 0 < validation_fraction < 1:
        raise ValueError("validation_fraction must be between zero and one")
    if group_column not in frame:
        raise ValueError(f"Missing split group column: {group_column}")
    groups = frame[group_column].astype("string")
    unique_groups = np.asarray(sorted(groups.dropna().unique().tolist()), dtype=object)
    if len(unique_groups) < 2:
        raise ValueError("At least two study groups are required for external validation")
    validation_count = max(1, int(math.ceil(len(unique_groups) * validation_fraction)))
    if validation_count >= len(unique_groups):
        validation_count = len(unique_groups) - 1
    rng = np.random.default_rng(random_state)
    validation_groups = set(rng.permutation(unique_groups)[:validation_count].tolist())
    validation_rows = groups.isin(validation_groups).to_numpy(dtype=bool)
    return ~validation_rows, validation_rows
