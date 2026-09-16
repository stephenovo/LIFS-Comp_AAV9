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

    if values.shape[1] == 3:
        return _pareto_mask_3d(values)

    keep = np.ones(values.shape[0], dtype=bool)
    for row in range(values.shape[0]):
        dominates_row = np.all(values >= values[row], axis=1) & np.any(
            values > values[row], axis=1
        )
        if dominates_row.any():
            keep[row] = False
    return keep


def _pareto_mask_3d(values: np.ndarray) -> np.ndarray:
    """Return a three-objective Pareto mask in O(n log n) time.

    Equal points are all retained. The implementation processes equal first
    objectives together, then uses a Fenwick tree to query the best third
    objective among points with an equal-or-better second objective.
    """
    row_count = len(values)
    if row_count == 0:
        return np.ones(0, dtype=bool)
    order = np.lexsort((-values[:, 2], -values[:, 1], -values[:, 0]))
    sorted_values = values[order]
    unique_y = np.unique(values[:, 1])[::-1]
    y_rank = np.searchsorted(-unique_y, -sorted_values[:, 1]) + 1
    tree = np.full(len(unique_y) + 1, -np.inf)
    keep_sorted = np.ones(row_count, dtype=bool)

    def query(index: int) -> float:
        best = -np.inf
        while index > 0:
            best = max(best, tree[index])
            index -= index & -index
        return best

    def update(index: int, value: float) -> None:
        while index < len(tree):
            tree[index] = max(tree[index], value)
            index += index & -index

    group_start = 0
    while group_start < row_count:
        group_end = group_start + 1
        x_value = sorted_values[group_start, 0]
        while group_end < row_count and sorted_values[group_end, 0] == x_value:
            group_end += 1

        best_z_at_higher_y = -np.inf
        y_start = group_start
        while y_start < group_end:
            y_end = y_start + 1
            y_value = sorted_values[y_start, 1]
            while y_end < group_end and sorted_values[y_end, 1] == y_value:
                y_end += 1
            same_y_max_z = sorted_values[y_start:y_end, 2].max()
            for row in range(y_start, y_end):
                z_value = sorted_values[row, 2]
                dominated_by_higher_x = query(int(y_rank[row])) >= z_value
                dominated_within_x = (
                    best_z_at_higher_y >= z_value or same_y_max_z > z_value
                )
                keep_sorted[row] = not (dominated_by_higher_x or dominated_within_x)
            best_z_at_higher_y = max(best_z_at_higher_y, same_y_max_z)
            y_start = y_end

        for row in range(group_start, group_end):
            update(int(y_rank[row]), float(sorted_values[row, 2]))
        group_start = group_end

    keep = np.empty(row_count, dtype=bool)
    keep[order] = keep_sorted
    return keep
