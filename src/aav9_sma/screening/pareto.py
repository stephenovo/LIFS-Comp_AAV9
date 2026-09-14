"""Small-dataset Pareto-front utilities."""

from __future__ import annotations

import numpy as np


def pareto_mask(objectives: np.ndarray) -> np.ndarray:
    """Return non-dominated rows when every objective is maximized."""
    values = np.asarray(objectives, dtype=float)
    if values.ndim != 2:
        raise ValueError("objectives must be a two-dimensional array")
    if not np.isfinite(values).all():
        raise ValueError("objectives contain missing or non-finite values")

    keep = np.ones(values.shape[0], dtype=bool)
    for row in range(values.shape[0]):
        dominates_row = np.all(values >= values[row], axis=1) & np.any(
            values > values[row], axis=1
        )
        if dominates_row.any():
            keep[row] = False
    return keep

