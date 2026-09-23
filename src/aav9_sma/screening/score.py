"""Hard packaging gate and transparent post-model candidate scoring."""

from __future__ import annotations

import numpy as np
import pandas as pd

from aav9_sma.screening.pareto import pareto_mask

PREDICTION_COLUMNS = (
    "pred_pack",
    "pred_brain_mouse",
    "pred_spinal_cord_mouse",
    "pred_liver_mouse",
    "pred_heart_mouse",
    "pred_kidney_mouse",
)


def rank_candidates(
    frame: pd.DataFrame,
    packaging_threshold: float,
    cns_weight: float = 0.45,
    liver_weight: float = 0.35,
    off_target_weight: float = 0.20,
    brain_target_weight: float = 0.30,
    spinal_target_weight: float = 0.70,
    epsilon: float = 1.0e-6,
) -> pd.DataFrame:
    """Apply the packaging gate and compute SMA-prioritized ranking fields.

    ``f_cns`` is retained as the unweighted brain/spinal mean for backwards-
    compatible reporting. Selection uses ``f_sma_target`` instead, with spinal
    cord prioritized by default. The two target weights are normalized so
    callers can supply any non-negative ratio.
    """
    missing = [column for column in PREDICTION_COLUMNS if column not in frame]
    if missing:
        raise ValueError(f"Missing prediction columns: {missing}")
    if brain_target_weight < 0 or spinal_target_weight < 0:
        raise ValueError("brain and spinal target weights must be non-negative")
    target_weight_total = brain_target_weight + spinal_target_weight
    if target_weight_total <= epsilon:
        raise ValueError("at least one brain/spinal target weight must be positive")
    normalized_brain_weight = brain_target_weight / target_weight_total
    normalized_spinal_weight = spinal_target_weight / target_weight_total

    ranked = frame.copy()
    packaging_value = (
        ranked["pred_pack_lcb"] if "pred_pack_lcb" in ranked else ranked["pred_pack"]
    )
    ranked["passes_packaging_gate"] = packaging_value >= packaging_threshold
    ranked["f_cns"] = ranked[["pred_brain_mouse", "pred_spinal_cord_mouse"]].mean(axis=1)
    ranked["f_sma_target"] = (
        normalized_brain_weight * ranked["pred_brain_mouse"]
        + normalized_spinal_weight * ranked["pred_spinal_cord_mouse"]
    )
    ranked["brain_target_weight"] = normalized_brain_weight
    ranked["spinal_target_weight"] = normalized_spinal_weight
    ranked["f_liv"] = ranked["pred_liver_mouse"]
    ranked["f_off"] = ranked[["pred_heart_mouse", "pred_kidney_mouse"]].mean(axis=1)
    ranked["display_score"] = (
        cns_weight * ranked["f_sma_target"]
        - liver_weight * ranked["f_liv"]
        - off_target_weight * ranked["f_off"]
    )
    ranked["log2_specificity"] = ranked["f_sma_target"] - ranked["f_liv"]
    ranked["specificity_index"] = np.exp2(
        ranked["log2_specificity"].clip(lower=-30, upper=30)
    )
    ranked["is_pareto"] = False

    eligible = ranked["passes_packaging_gate"]
    if eligible.any():
        objectives = ranked.loc[eligible, ["f_sma_target", "f_liv", "f_off"]].copy()
        objectives[["f_liv", "f_off"]] *= -1
        ranked.loc[eligible, "is_pareto"] = pareto_mask(objectives.to_numpy())

    return ranked.sort_values(
        ["passes_packaging_gate", "is_pareto", "display_score"],
        ascending=[False, False, False],
        kind="stable",
    ).reset_index(drop=True)
