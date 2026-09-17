import numpy as np
import pytest

from aav9_sma.features.encode import (
    combined_7mer_features,
    one_hot_7mer,
    physicochemical_7mer,
)


def test_one_hot_shape_and_position_sums() -> None:
    encoded = one_hot_7mer(["ACDEFGH", "YYYYYYY"])

    assert encoded.shape == (2, 140)
    np.testing.assert_allclose(encoded.sum(axis=1), [7.0, 7.0])


def test_one_hot_rejects_invalid_peptide() -> None:
    with pytest.raises(ValueError, match="Invalid 7-mer"):
        one_hot_7mer(["TOO-LONG"])


def test_physicochemical_encoding_is_fixed_and_composable() -> None:
    peptides = ["ACDEFGH", "KKKKKKK"]
    physicochemical = physicochemical_7mer(peptides)
    combined = combined_7mer_features(peptides)

    assert physicochemical.shape == (2, 56)
    assert combined.shape == (2, 196)
    np.testing.assert_allclose(combined[:, :140], one_hot_7mer(peptides))
    assert combined[1, -3] == 7.0
