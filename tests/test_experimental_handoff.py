import pandas as pd

from aav9_sma.experimental.handoff import build_validation_panel, select_empirical_controls


def _controls() -> pd.DataFrame:
    rows = []
    for group in (
        "high_packaging",
        "low_packaging",
        "cns_high_liver_low",
        "liver_high_cns_low",
    ):
        for value in (1.0, 3.0):
            rows.append(
                {
                    "control_group": group,
                    "AA": f"SEQ{len(rows):04d}"[-7:],
                    "variant_id": f"CTRL_{len(rows):03d}",
                    "Production2": value,
                    "control_axis": value,
                    "empirical_cns": value,
                    "empirical_liver": -value,
                }
            )
    return pd.DataFrame(rows)


def _shortlist() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "variant_id": ["V1", "V2", "V3"],
            "AA": ["AAAAAAA", "CCCCCCC", "DDDDDDD"],
            "selection_group": ["balanced", "cns_favoring", "low_liver"],
            "selection_rank": [1, 1, 1],
            "strict_conservative": [True, False, False],
            "is_pareto": [True, True, False],
            "pred_pack_lcb": [1.0, 1.0, 1.0],
            "pred_brain_mouse": [1.0, 1.0, 1.0],
            "pred_spinal_cord_mouse": [1.0, 1.0, 1.0],
            "pred_liver_mouse": [-1.0, -1.0, -1.0],
            "pred_heart_mouse": [0.0, 0.0, 0.0],
            "pred_kidney_mouse": [0.0, 0.0, 0.0],
            "log2_specificity": [2.0, 2.0, 2.0],
            "organ_uncertainty_mean": [0.1, 0.1, 0.1],
            "human_liver_warning": [False, False, True],
        }
    )


def test_control_selection_uses_declared_extremes() -> None:
    selected = select_empirical_controls(_controls()).set_index("control_group")
    assert selected.loc["high_packaging", "Production2"] == 3.0
    assert selected.loc["low_packaging", "Production2"] == 1.0
    assert selected.loc["cns_high_liver_low", "control_axis"] == 3.0
    assert selected.loc["liver_high_cns_low", "control_axis"] == 3.0


def test_panel_contains_candidates_and_required_controls() -> None:
    panel = build_validation_panel(_shortlist(), _controls())
    assert len(panel) == 10
    assert panel["sample_id"].is_unique
    assert set(panel.loc[panel["panel_role"] == "computational_candidate", "priority_tier"]) == {
        "Tier 1",
        "Tier 2",
        "Tier 3",
    }
    assert {"PARENT_AAV9_K449R", "WILD_TYPE_AAV9", "NO_VECTOR_MOCK"} <= set(
        panel["variant_id"]
    )
