import numpy as np
import pandas as pd

from aav9_sma.screening.pareto import pareto_mask
from aav9_sma.screening.score import rank_candidates
from aav9_sma.screening.virtual import (
    add_weight_sensitivity,
    generate_candidate_peptides,
    select_diverse_shortlist,
    training_distance_lower_bound,
)


def test_packaging_gate_and_ranking_fields() -> None:
    frame = pd.DataFrame(
        {
            "variant_id": ["good", "low_pack", "low_liver"],
            "pred_pack": [0.8, 0.2, 0.8],
            "pred_brain_mouse": [0.8, 1.0, 0.6],
            "pred_spinal_cord_mouse": [0.8, 1.0, 0.6],
            "pred_liver_mouse": [0.2, 0.1, 0.05],
            "pred_heart_mouse": [0.2, 0.1, 0.1],
            "pred_kidney_mouse": [0.2, 0.1, 0.1],
        }
    )

    ranked = rank_candidates(frame, packaging_threshold=0.5)

    assert ranked.iloc[-1]["variant_id"] == "low_pack"
    assert not bool(ranked.iloc[-1]["passes_packaging_gate"])
    assert {"f_cns", "f_liv", "f_off", "display_score", "specificity_index"} <= set(
        ranked.columns
    )
    assert ranked.loc[ranked["variant_id"] == "good", "is_pareto"].item()
    assert ranked.loc[ranked["variant_id"] == "low_liver", "is_pareto"].item()
    good = ranked.loc[ranked["variant_id"] == "good"].iloc[0]
    assert np.isclose(good["log2_specificity"], 0.6)
    assert np.isclose(good["specificity_index"], 2**0.6)


def test_fast_three_dimensional_pareto_matches_brute_force() -> None:
    rng = np.random.default_rng(7)
    values = rng.integers(0, 8, size=(150, 3)).astype(float)
    expected = np.ones(len(values), dtype=bool)
    for row in range(len(values)):
        expected[row] = not np.any(
            np.all(values >= values[row], axis=1)
            & np.any(values > values[row], axis=1)
        )
    np.testing.assert_array_equal(pareto_mask(values), expected)


def test_generation_distance_and_diverse_shortlist() -> None:
    generated = generate_candidate_peptides(
        100, excluded={"AAAAAAA", "CAAAAAA"}, random_state=9
    )
    assert len(generated) == len(set(generated)) == 100
    assert not {"AAAAAAA", "CAAAAAA"}.intersection(generated)
    np.testing.assert_array_equal(
        training_distance_lower_bound(
            ["AAAAAAA", "CAAAAAA", "CCAAAAA", "CCCCCCC"],
            ["AAAAAAA"],
        ),
        [0, 1, 2, 2],
    )

    frame = pd.DataFrame(
        {
            "AA": ["AAAAAAA", "CCCCCCC", "DDDDDDD", "EEEEEEE", "FFFFFFF", "GGGGGGG"],
            "passes_packaging_gate": True,
            "is_pareto": [True, True, True, False, False, False],
            "training_distance_lower_bound": 2,
            "f_cns": [6, 5, 4, 3, 2, 1],
            "f_liv": [6, 5, 4, 3, 2, 1],
            "display_score": [6, 6, 6, 6, 6, 6],
            "weight_stability_top_fraction": [1, 0.8, 0.6, 0.4, 0.2, 0],
            "organ_uncertainty_mean": 0.1,
        }
    )
    shortlist = select_diverse_shortlist(frame, per_group=1)
    assert len(shortlist) == 3
    assert set(shortlist["selection_group"]) == {
        "cns_favoring",
        "low_liver",
        "balanced",
    }


def test_weight_sensitivity_is_bounded() -> None:
    frame = pd.DataFrame(
        {
            "passes_packaging_gate": [True, True, False],
            "f_cns": [2.0, 1.0, 5.0],
            "f_liv": [0.0, 0.5, 0.0],
            "f_off": [0.0, 0.5, 0.0],
        }
    )
    output = add_weight_sensitivity(frame, top_fraction=0.5)
    assert output["weight_stability_top_fraction"].between(0, 1).all()
    assert output["weight_mean_percentile"].between(0, 1).all()
    assert output.loc[2, "weight_stability_top_fraction"] == 0
