import numpy as np
import pytest

from aav9_sma.features.encode import one_hot_7mer


def test_one_hot_shape_and_position_sums() -> None:
    encoded = one_hot_7mer(["ACDEFGH", "YYYYYYY"])

    assert encoded.shape == (2, 140)
    np.testing.assert_allclose(encoded.sum(axis=1), [7.0, 7.0])


def test_one_hot_rejects_invalid_peptide() -> None:
    with pytest.raises(ValueError, match="Invalid 7-mer"):
        one_hot_7mer(["TOO-LONG"])

