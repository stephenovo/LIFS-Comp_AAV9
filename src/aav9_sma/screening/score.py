"""Hard packaging gate and transparent post-model candidate scoring."""

from __future__ import annotations

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
    epsilon: float = 1.0e-6,
) -> pd.DataFrame:
    """Apply the packaging gate and compute display/Pareto ranking fields."""
    missing = [column for column in PREDICTION_COLUMNS if column not in frame]
    if missing:
        raise ValueError(f"Missing prediction columns: {missing}")

    ranked = frame.copy()
    packaging_value = (
        ranked["pred_pack_lcb"] if "pred_pack_lcb" in ranked else ranked["pred_pack"]
    )
    ranked["passes_packaging_gate"] = packaging_value >= packaging_threshold
    ranked["f_cns"] = ranked[["pred_brain_mouse", "pred_spinal_cord_mouse"]].mean(axis=1)
    ranked["f_liv"] = ranked["pred_liver_mouse"]
    ranked["f_off"] = ranked[["pred_heart_mouse", "pred_kidney_mouse"]].mean(axis=1)
    ranked["display_score"] = (
        cns_weight * ranked["f_cns"]
        - liver_weight * ranked["f_liv"]
        - off_target_weight * ranked["f_off"]
    )
    ranked["specificity_index"] = ranked["f_cns"] / (ranked["f_liv"] + epsilon)
    ranked["is_pareto"] = False

    eligible = ranked["passes_packaging_gate"]
    if eligible.any():
        objectives = ranked.loc[eligible, ["f_cns", "f_liv", "f_off"]].copy()
        objectives[["f_liv", "f_off"]] *= -1
        ranked.loc[eligible, "is_pareto"] = pareto_mask(objectives.to_numpy())

    return ranked.sort_values(
        ["passes_packaging_gate", "is_pareto", "display_score"],
        ascending=[False, False, False],
        kind="stable",
    ).reset_index(drop=True)

