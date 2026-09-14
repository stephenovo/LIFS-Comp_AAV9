"""Deterministic baseline encodings for seven-amino-acid peptides."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from aav9_sma.constants import AMINO_ACIDS, PEPTIDE_LENGTH
from aav9_sma.data.audit import is_valid_peptide


def one_hot_7mer(peptides: Sequence[str]) -> np.ndarray:
    """Encode peptides as a flattened position-by-amino-acid one-hot matrix."""
    residue_to_index = {residue: index for index, residue in enumerate(AMINO_ACIDS)}
    output = np.zeros(
        (len(peptides), PEPTIDE_LENGTH * len(AMINO_ACIDS)),
        dtype=np.float32,
    )

    for row, peptide in enumerate(peptides):
        if not is_valid_peptide(peptide):
            raise ValueError(f"Invalid 7-mer at row {row}: {peptide!r}")
        for position, residue in enumerate(peptide):
            column = position * len(AMINO_ACIDS) + residue_to_index[residue]
            output[row, column] = 1.0
    return output

