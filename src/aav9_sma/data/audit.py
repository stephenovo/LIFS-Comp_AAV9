"""Audit data after source fields have been mapped to the canonical schema."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import pandas as pd

from aav9_sma.constants import (
    AMINO_ACIDS,
    CORE_ID_COLUMNS,
    PEPTIDE_LENGTH,
    PRIMARY_LABEL_COLUMNS,
)


@dataclass(frozen=True)
class AuditResult:
    row_count: int
    duplicate_variant_ids: int
    duplicate_peptides: int
    invalid_peptide_rows: int
    missing_columns: tuple[str, ...]
    non_null_counts: dict[str, int]
    complete_primary_rows: int

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def is_valid_peptide(value: object) -> bool:
    """Return whether a value is an uppercase 7-mer over standard amino acids."""
    if not isinstance(value, str) or len(value) != PEPTIDE_LENGTH:
        return False
    return all(residue in AMINO_ACIDS for residue in value)


def audit_dataframe(frame: pd.DataFrame) -> AuditResult:
    """Summarize schema coverage and obvious quality issues."""
    required = (*CORE_ID_COLUMNS, *PRIMARY_LABEL_COLUMNS)
    missing_columns = tuple(column for column in required if column not in frame.columns)

    duplicate_variant_ids = (
        int(frame["variant_id"].duplicated().sum()) if "variant_id" in frame else 0
    )
    duplicate_peptides = (
        int(frame["peptide_7mer"].duplicated().sum()) if "peptide_7mer" in frame else 0
    )
    invalid_peptide_rows = (
        int((~frame["peptide_7mer"].map(is_valid_peptide)).sum())
        if "peptide_7mer" in frame
        else len(frame)
    )

    observed_labels = [column for column in PRIMARY_LABEL_COLUMNS if column in frame]
    non_null_counts = {column: int(frame[column].notna().sum()) for column in observed_labels}
    complete_primary_rows = (
        int(frame[list(PRIMARY_LABEL_COLUMNS)].notna().all(axis=1).sum())
        if not missing_columns and len(frame)
        else 0
    )

    return AuditResult(
        row_count=len(frame),
        duplicate_variant_ids=duplicate_variant_ids,
        duplicate_peptides=duplicate_peptides,
        invalid_peptide_rows=invalid_peptide_rows,
        missing_columns=missing_columns,
        non_null_counts=non_null_counts,
        complete_primary_rows=complete_primary_rows,
    )


def audit_csv(path: str | Path) -> AuditResult:
    """Load a CSV file and return its canonical-schema audit."""
    return audit_dataframe(pd.read_csv(path))

