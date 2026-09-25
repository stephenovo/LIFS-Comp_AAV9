import numpy as np
import pandas as pd

from aav9_sma.screening.audit import (
    build_composition_audit,
    build_gate_sensitivity_audit,
    select_composition_challenge_panel,
    summarize_funnel_audit,
)
from aav9_sma.screening.pareto import pareto_mask
from aav9_sma.screening.score import rank_candidates
from aav9_sma.screening.strategies import (
    SequentialStage,
    compare_strategy_shortlists,
    run_sequential_strategy,
)
from aav9_sma.screening.virtual import (
    add_conservative_annotations,
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
    assert {
        "f_cns",
        "f_sma_target",
        "f_liv",
        "f_off",
        "display_score",
        "specificity_index",
    } <= set(ranked.columns)
    assert ranked.loc[ranked["variant_id"] == "good", "is_pareto"].item()
    assert ranked.loc[ranked["variant_id"] == "low_liver", "is_pareto"].item()
    good = ranked.loc[ranked["variant_id"] == "good"].iloc[0]
    assert np.isclose(good["log2_specificity"], 0.6)
    assert np.isclose(good["specificity_index"], 2**0.6)


def test_sma_target_prioritizes_spinal_cord_over_brain() -> None:
    frame = pd.DataFrame(
        {
            "variant_id": ["brain_high", "spinal_high"],
            "pred_pack": [1.0, 1.0],
            "pred_brain_mouse": [1.0, 0.0],
            "pred_spinal_cord_mouse": [0.0, 1.0],
            "pred_liver_mouse": [0.0, 0.0],
            "pred_heart_mouse": [0.0, 0.0],
            "pred_kidney_mouse": [0.0, 0.0],
        }
    )

    baseline = rank_candidates(frame, packaging_threshold=0.5)
    assert baseline["f_sma_target"].eq(0.5).all()
    assert baseline["brain_target_weight"].eq(0.5).all()
    assert baseline["spinal_target_weight"].eq(0.5).all()

    ranked = rank_candidates(
        frame,
        packaging_threshold=0.5,
        brain_target_weight=0.30,
        spinal_target_weight=0.70,
    )

    assert ranked.iloc[0]["variant_id"] == "spinal_high"
    assert np.isclose(ranked.iloc[0]["f_sma_target"], 0.7)
    assert np.isclose(ranked.iloc[1]["f_sma_target"], 0.3)


def test_fast_three_dimensional_pareto_matches_brute_force() -> None:
    rng = np.random.default_rng(7)
    values = rng.integers(0, 8, size=(150, 3)).astype(float)
    expected = np.ones(len(values), dtype=bool)
    for row in range(len(values)):
        expected[row] = not np.any(
            np.all(values >= values[row], axis=1) & np.any(values > values[row], axis=1)
        )
    np.testing.assert_array_equal(pareto_mask(values), expected)


def test_generation_distance_and_diverse_shortlist() -> None:
    generated = generate_candidate_peptides(100, excluded={"AAAAAAA", "CAAAAAA"}, random_state=9)
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
            "pred_spinal_cord_mouse": [6, 5, 4, 3, 2, 1],
            "pred_brain_mouse": [1, 2, 3, 4, 5, 6],
            "f_cns": [6, 5, 4, 3, 2, 1],
            "f_sma_target": [6, 5, 4, 3, 2, 1],
            "f_liv": [6, 5, 4, 3, 2, 1],
            "display_score": [6, 6, 6, 6, 6, 6],
            "weight_stability_top_fraction": [1, 0.8, 0.6, 0.4, 0.2, 0],
            "organ_uncertainty_mean": 0.1,
        }
    )
    shortlist = select_diverse_shortlist(frame, per_group=1)
    assert len(shortlist) == 3
    assert set(shortlist["selection_group"]) == {
        "spinal_favoring",
        "low_liver",
        "balanced",
    }


def test_weight_sensitivity_is_bounded() -> None:
    frame = pd.DataFrame(
        {
            "passes_packaging_gate": [True, True, False],
            "f_cns": [2.0, 1.0, 5.0],
            "f_sma_target": [2.0, 1.0, 5.0],
            "f_liv": [0.0, 0.5, 0.0],
            "f_off": [0.0, 0.5, 0.0],
        }
    )
    output = add_weight_sensitivity(frame, top_fraction=0.5)
    assert output["weight_stability_top_fraction"].between(0, 1).all()
    assert output["weight_mean_percentile"].between(0, 1).all()
    assert output.loc[2, "weight_stability_top_fraction"] == 0


def test_sequential_funnel_applies_spinal_then_brain_then_liver() -> None:
    ranked = pd.DataFrame(
        {
            "AA": ["AAAAAAA", "CCCCCCC", "DDDDDDD", "EEEEEEE"],
            "passes_packaging_gate": [True] * 4,
            "training_distance_lower_bound": [2] * 4,
            "pred_spinal_cord_mouse": [4.0, 3.0, 2.0, 1.0],
            "pred_brain_mouse": [1.0, 4.0, 3.0, 2.0],
            "pred_liver_mouse": [4.0, 1.0, 2.0, 3.0],
            "pred_heart_mouse": [0.0] * 4,
            "pred_kidney_mouse": [0.0] * 4,
            "display_score": [1.0] * 4,
        }
    )
    stages = (
        SequentialStage("spinal", "pred_spinal_cord_mouse", "maximize", 0.50),
        SequentialStage("brain", "pred_brain_mouse", "maximize", 0.50),
        SequentialStage("liver", "pred_liver_mouse", "minimize", 1.00),
    )

    annotated, shortlist, audit = run_sequential_strategy(
        ranked,
        stages,
        candidate_count=2,
        minimum_pairwise_distance=1,
    )

    assert audit[1]["stage"] == "spinal"
    assert audit[2]["stage"] == "brain"
    assert audit[1]["survivor_rows"] == 2
    assert audit[2]["survivor_rows"] == 1
    assert annotated.loc[annotated["passes_sequential_funnel"], "AA"].tolist() == ["CCCCCCC"]
    assert shortlist["AA"].tolist() == ["CCCCCCC"]


def test_strategy_comparison_quantifies_overlap_and_composition() -> None:
    joint = pd.DataFrame(
        {
            "AA": ["AAAAAAA", "CCCCCCC"],
            "display_score": [2.0, 1.0],
            "pred_spinal_cord_mouse": [1.0, 2.0],
            "pred_brain_mouse": [2.0, 1.0],
        }
    )
    sequential = pd.DataFrame(
        {
            "AA": ["CCCCCCC", "DDDDDDD"],
            "display_score": [1.0, 0.5],
            "pred_spinal_cord_mouse": [2.0, 3.0],
            "pred_brain_mouse": [1.0, 0.0],
            "strategy_selection_rank": [1, 2],
        }
    )

    comparison, summary = compare_strategy_shortlists(joint, sequential)

    assert set(comparison["selection_relationship"]) == {
        "shared",
        "joint_only",
        "sequential_only",
    }
    assert summary["shared_rows"] == 1
    assert np.isclose(summary["jaccard_similarity"], 1 / 3)
    assert summary["predicted_phenotype_medians"]["pred_spinal_cord_mouse"][
        "sequential_minus_joint"
    ] == 1.0


def test_strict_conservative_annotation_uses_training_medians() -> None:
    ranked = pd.DataFrame(
        {
            "passes_packaging_gate": [True, True, False],
            "pred_brain_mouse": [1.0, -1.0, 1.0],
            "pred_spinal_cord_mouse": [1.0, -1.0, 1.0],
            "pred_liver_mouse": [-1.0, 1.0, -1.0],
            "organ_uncertainty_mean": [0.1, 0.2, 0.05],
        }
    )
    reconstructed = pd.DataFrame(
        {
            "log2enr_whitelist__brain_animals_1_3__over__virus_prod2": [-0.5, 0.5],
            "log2enr_whitelist__spinal_cord_animals_1_3__over__virus_prod2": [-0.5, 0.5],
            "log2enr_whitelist__liver_animals_1_3__over__virus_prod2": [-0.5, 0.5],
        }
    )

    output, thresholds = add_conservative_annotations(ranked, reconstructed)

    assert output["strict_conservative"].tolist() == [True, False, False]
    assert output["passes_cns_median"].tolist() == [True, False, True]
    assert thresholds["brain_training_median"] == 0.0
    assert thresholds["spinal_cord_training_median"] == 0.0
    assert thresholds["liver_training_median"] == 0.0


def test_funnel_audits_trace_gate_and_composition_changes() -> None:
    ranked = pd.DataFrame(
        {
            "AA": ["AAAAAAA", "CCCCCCC", "DDDDDDD", "EEEEEEE"],
            "pred_pack": [3.0, 2.0, 1.0, 0.0],
            "passes_packaging_gate": [True, False, False, False],
            "is_pareto": [True, False, False, False],
            "display_score": [4.0, 3.0, 2.0, 1.0],
            "training_distance_lower_bound": [2, 2, 2, 2],
        }
    )
    shortlist = ranked.iloc[[0]][["AA"]]
    observed = pd.DataFrame({"AA": ["ACDEFGH", "IKLMNPQ"]})

    gates = build_gate_sensitivity_audit(
        ranked,
        packaging_threshold=0.5,
        lower_bound_offsets={
            "strict_95_lcb": 2.0,
            "sensitivity_90_lcb": 1.0,
            "sensitivity_point_prediction": 0.0,
        },
    )
    composition = build_composition_audit(ranked, shortlist, observed)
    summary = summarize_funnel_audit(gates, composition)

    assert gates.set_index("gate").loc["strict_95_lcb", "passed"] == 1
    assert gates.set_index("gate").loc["sensitivity_point_prediction", "passed"] == 3
    assert set(composition["stage"]) == {
        "observed_100k_library",
        "virtual_pool",
        "packaging_95_lcb",
        "pareto_after_95_lcb",
        "quality_top5_after_95_lcb",
        "final_shortlist",
    }
    assert summary["interpretation"]["primary_shortlist_changed"] is False


def test_composition_challenge_keeps_primary_shortlist_separate() -> None:
    ranked = pd.DataFrame(
        {
            "variant_id": ["c95", "c90", "f95", "f90", "main"],
            "AA": ["CAAAAAD", "ACDDDDD", "FAAAAAG", "AFEEEEE", "AAAAAAA"],
            "pred_pack": [3.0, 2.4, 3.1, 2.3, 4.0],
            "pred_pack_lcb": [1.0, 0.4, 1.1, 0.3, 2.0],
            "display_score": [1.0, 0.0, 0.8, -0.1, 2.0],
            "f_cns": [0.0] * 5,
            "f_liv": [0.0] * 5,
            "f_off": [0.0] * 5,
            "organ_uncertainty_mean": [0.1] * 5,
            "training_distance_lower_bound": [2] * 5,
            "human_liver_warning": [False] * 5,
        }
    )
    shortlist = pd.DataFrame({"AA": ["AAAAAAA"]})

    panel, summary = select_composition_challenge_panel(
        ranked,
        shortlist,
        packaging_threshold=0.5,
        strict_95_offset=2.0,
        sensitivity_90_offset=1.0,
        target_residues=("C", "F"),
        minimum_pairwise_distance=1,
    )

    assert len(panel) == 4
    assert set(panel["gate_stratum"]) == {"strict_95_survivor", "sensitivity_90_only"}
    assert not panel["replaces_primary_shortlist"].any()
    assert summary["primary_shortlist_changed"] is False
