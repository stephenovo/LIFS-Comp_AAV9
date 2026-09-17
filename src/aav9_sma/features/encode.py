"""Deterministic encodings for seven-amino-acid peptides."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from aav9_sma.constants import AMINO_ACIDS, PEPTIDE_LENGTH
from aav9_sma.data.audit import is_valid_peptide

# Transparent residue descriptors. Hydropathy uses the Kyte-Doolittle scale;
# the remaining values are deliberately simple chemical class indicators.
# Every continuous column is standardized across the 20-residue alphabet so
# the feature scale is fixed before any train/test split is observed.
_HYDROPATHY = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3,
}
_RESIDUE_VOLUME = {
    "A": 88.6,
    "C": 108.5,
    "D": 111.1,
    "E": 138.4,
    "F": 189.9,
    "G": 60.1,
    "H": 153.2,
    "I": 166.7,
    "K": 168.6,
    "L": 166.7,
    "M": 162.9,
    "N": 114.1,
    "P": 112.7,
    "Q": 143.8,
    "R": 173.4,
    "S": 89.0,
    "T": 116.1,
    "V": 140.0,
    "W": 227.8,
    "Y": 193.6,
}
_POSITIVE = frozenset("KRH")
_NEGATIVE = frozenset("DE")
_AROMATIC = frozenset("FWY")
_POLAR = frozenset("CDEHKNQRSTY")
_SPECIAL = frozenset("CGP")


def _standardized_lookup(values: dict[str, float]) -> dict[str, float]:
    ordered = np.asarray([values[residue] for residue in AMINO_ACIDS], dtype=np.float32)
    return {
        residue: float((values[residue] - ordered.mean()) / ordered.std())
        for residue in AMINO_ACIDS
    }


_HYDROPATHY_Z = _standardized_lookup(_HYDROPATHY)
_VOLUME_Z = _standardized_lookup(_RESIDUE_VOLUME)


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


def physicochemical_7mer(peptides: Sequence[str]) -> np.ndarray:
    """Encode positional and whole-peptide physicochemical descriptors.

    Each residue contributes seven fixed values: hydropathy, volume, charge,
    and aromatic/polar/special-class indicators. Seven whole-peptide summaries
    are appended. No statistic is fitted to the modeling dataset, preventing
    information leakage when this encoding is used on a held-out split.
    """
    per_residue_features = 7
    summary_features = 7
    output = np.zeros(
        (len(peptides), PEPTIDE_LENGTH * per_residue_features + summary_features),
        dtype=np.float32,
    )
    for row, peptide in enumerate(peptides):
        if not is_valid_peptide(peptide):
            raise ValueError(f"Invalid 7-mer at row {row}: {peptide!r}")
        hydropathy = []
        volume = []
        charge = []
        for position, residue in enumerate(peptide):
            residue_charge = 1.0 if residue in _POSITIVE else -1.0 if residue in _NEGATIVE else 0.0
            values = (
                _HYDROPATHY_Z[residue],
                _VOLUME_Z[residue],
                residue_charge,
                float(residue in _AROMATIC),
                float(residue in _POLAR),
                float(residue in _SPECIAL),
                float(residue == "C"),
            )
            start = position * per_residue_features
            output[row, start : start + per_residue_features] = values
            hydropathy.append(_HYDROPATHY_Z[residue])
            volume.append(_VOLUME_Z[residue])
            charge.append(residue_charge)
        output[row, -summary_features:] = (
            np.mean(hydropathy),
            np.std(hydropathy),
            np.mean(volume),
            np.std(volume),
            np.sum(charge),
            sum(residue in _AROMATIC for residue in peptide) / PEPTIDE_LENGTH,
            sum(residue in _POLAR for residue in peptide) / PEPTIDE_LENGTH,
        )
    return output


def combined_7mer_features(peptides: Sequence[str]) -> np.ndarray:
    """Concatenate one-hot and fixed physicochemical descriptors."""
    return np.concatenate([one_hot_7mer(peptides), physicochemical_7mer(peptides)], axis=1)
