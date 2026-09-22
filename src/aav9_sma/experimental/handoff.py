"""Build a traceable wet-lab validation panel from screening outputs.

The generated panel is a planning artifact. It does not replace a vector core's
validated production procedures, biosafety review, or animal ethics approval.
"""

from __future__ import annotations

from collections.abc import Mapping

import pandas as pd

CONTROL_RULES: Mapping[str, tuple[str, str]] = {
    "high_packaging": ("Production2", "max"),
    "low_packaging": ("Production2", "min"),
    "cns_high_liver_low": ("control_axis", "max"),
    "liver_high_cns_low": ("control_axis", "max"),
}


def _as_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.lower().eq("true")


def select_empirical_controls(control_audit: pd.DataFrame) -> pd.DataFrame:
    """Select one deterministic empirical control for each audit direction."""
    required = {"control_group", "AA", "variant_id", "Production2", "control_axis"}
    missing = required - set(control_audit.columns)
    if missing:
        raise ValueError(f"control audit is missing columns: {sorted(missing)}")

    selected: list[pd.Series] = []
    for group, (column, direction) in CONTROL_RULES.items():
        group_rows = control_audit.loc[control_audit["control_group"] == group].copy()
        if group_rows.empty:
            raise ValueError(f"control audit has no rows for {group!r}")
        values = pd.to_numeric(group_rows[column], errors="coerce")
        if not values.notna().any():
            raise ValueError(f"control group {group!r} has no finite {column!r} values")
        index = values.idxmax() if direction == "max" else values.idxmin()
        selected.append(group_rows.loc[index])
    return pd.DataFrame(selected).reset_index(drop=True)


def build_validation_panel(shortlist: pd.DataFrame, control_audit: pd.DataFrame) -> pd.DataFrame:
    """Create the candidate and control manifest used by the validation SOP."""
    candidate_required = {
        "variant_id",
        "AA",
        "selection_group",
        "selection_rank",
        "strict_conservative",
        "is_pareto",
        "passes_cns_median",
        "pred_pack_lcb",
        "pred_brain_mouse",
        "pred_spinal_cord_mouse",
        "pred_liver_mouse",
        "pred_heart_mouse",
        "pred_kidney_mouse",
        "log2_specificity",
        "organ_uncertainty_mean",
        "human_liver_warning",
    }
    missing = candidate_required - set(shortlist.columns)
    if missing:
        raise ValueError(f"shortlist is missing columns: {sorted(missing)}")

    candidates = shortlist.copy()
    candidates["strict_conservative"] = _as_bool(candidates["strict_conservative"])
    candidates["is_pareto"] = _as_bool(candidates["is_pareto"])
    candidates["passes_cns_median"] = _as_bool(candidates["passes_cns_median"])
    candidates["human_liver_warning"] = _as_bool(candidates["human_liver_warning"])
    candidates = candidates.sort_values(
        ["strict_conservative", "is_pareto", "selection_group", "selection_rank"],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)
    candidates["sample_id"] = [f"CAND-{index:03d}" for index in range(1, len(candidates) + 1)]
    candidates["panel_role"] = "computational_candidate"
    candidates["priority_tier"] = "Tier 3"
    tier_2 = (
        candidates["is_pareto"]
        & candidates["passes_cns_median"]
        & candidates["log2_specificity"].gt(0)
    )
    candidates.loc[tier_2, "priority_tier"] = "Tier 2"
    candidates.loc[candidates["strict_conservative"], "priority_tier"] = "Tier 1"
    candidates["stage_1_plan"] = "package_and_qc"
    candidates["stage_2_plan"] = (
        "hold_after_stage_1; advance only by preregistered diversity challenge or measured QC"
    )
    candidates.loc[
        candidates["priority_tier"].isin(["Tier 1", "Tier 2"]), "stage_2_plan"
    ] = "advance_if_stage_1_passes; prioritize by measured QC"
    candidates["stage_3_plan"] = "rank_after_stage_2; test top 4-6"
    candidates["selection_basis"] = candidates.apply(
        lambda row: (
            "strict_conservative_shortlist"
            if row["strict_conservative"]
            else "target_plausible_pareto_shortlist"
            if row["priority_tier"] == "Tier 2"
            else "exploratory_pareto_boundary_case"
            if row["is_pareto"]
            else "diverse_high_scoring_near_front"
        ),
        axis=1,
    )
    candidates["source_status"] = "model_prediction_unvalidated"

    controls = select_empirical_controls(control_audit)
    controls = controls.rename(
        columns={
            "empirical_cns": "observed_cns_reference",
            "empirical_liver": "observed_liver_reference",
        }
    )
    controls["sample_id"] = [f"CTRL-EMP-{index:02d}" for index in range(1, len(controls) + 1)]
    controls["panel_role"] = "empirical_directionality_control"
    controls["priority_tier"] = "Control"
    controls["selection_group"] = controls["control_group"]
    controls["selection_rank"] = 1
    controls["strict_conservative"] = False
    controls["is_pareto"] = False
    controls["human_liver_warning"] = False
    controls["stage_1_plan"] = "package_and_qc"
    controls["stage_2_plan"] = "assay_as_directionality_control_if_stage_1_passes"
    controls["stage_3_plan"] = "optional; use only if production is adequate"
    controls["selection_basis"] = "extreme_observed_fit4function_control"
    controls["source_status"] = "in_sample_empirical_control_not_independent_validation"

    fixed_controls = pd.DataFrame(
        [
            {
                "sample_id": "CTRL-PARENT-01",
                "variant_id": "PARENT_AAV9_K449R",
                "AA": "",
                "panel_role": "parental_backbone_control",
                "priority_tier": "Control",
                "selection_group": "parental_control",
                "selection_rank": 1,
                "strict_conservative": False,
                "is_pareto": False,
                "human_liver_warning": False,
                "stage_1_plan": "package_and_qc",
                "stage_2_plan": "required_reference",
                "stage_3_plan": "required_reference",
                "selection_basis": "same_AAV9_K449R_backbone_without_7mer_insertion",
                "source_status": "required_experimental_reference",
            },
            {
                "sample_id": "CTRL-AAV9WT-01",
                "variant_id": "WILD_TYPE_AAV9",
                "AA": "",
                "panel_role": "external_clinical_context_control",
                "priority_tier": "Control",
                "selection_group": "external_control",
                "selection_rank": 1,
                "strict_conservative": False,
                "is_pareto": False,
                "human_liver_warning": False,
                "stage_1_plan": "optional_package_and_qc",
                "stage_2_plan": "optional_context_reference",
                "stage_3_plan": "optional_context_reference",
                "selection_basis": "wild_type_AAV9_context; not interchangeable_with_K449R_parent",
                "source_status": "optional_experimental_reference",
            },
            {
                "sample_id": "CTRL-MOCK-01",
                "variant_id": "NO_VECTOR_MOCK",
                "AA": "",
                "panel_role": "negative_control",
                "priority_tier": "Control",
                "selection_group": "negative_control",
                "selection_rank": 1,
                "strict_conservative": False,
                "is_pareto": False,
                "human_liver_warning": False,
                "stage_1_plan": "process_blank_where_applicable",
                "stage_2_plan": "required_negative_control",
                "stage_3_plan": "required_negative_control",
                "selection_basis": "background_and_assay_specificity_control",
                "source_status": "required_experimental_reference",
            },
        ]
    )

    output_columns = [
        "sample_id",
        "variant_id",
        "AA",
        "panel_role",
        "priority_tier",
        "selection_group",
        "selection_rank",
        "selection_basis",
        "strict_conservative",
        "is_pareto",
        "passes_cns_median",
        "pred_pack_lcb",
        "pred_brain_mouse",
        "pred_spinal_cord_mouse",
        "pred_liver_mouse",
        "pred_heart_mouse",
        "pred_kidney_mouse",
        "log2_specificity",
        "organ_uncertainty_mean",
        "human_liver_warning",
        "Production2",
        "observed_cns_reference",
        "observed_liver_reference",
        "stage_1_plan",
        "stage_2_plan",
        "stage_3_plan",
        "source_status",
    ]
    panel = pd.concat([candidates, controls, fixed_controls], ignore_index=True, sort=False)
    for column in output_columns:
        if column not in panel:
            panel[column] = pd.NA
    return panel.loc[:, output_columns]
