import pandas as pd

from aav9_sma.screening.score import rank_candidates


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

