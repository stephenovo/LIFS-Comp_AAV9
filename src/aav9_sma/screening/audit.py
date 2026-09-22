"""Post-screen audits for packaging gates and amino-acid composition."""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from aav9_sma.constants import AMINO_ACIDS

DEFAULT_COMPOSITION_CHALLENGE_RESIDUES = tuple("CFIMWY")


def _hamming_distance(left: str, right: str) -> int:
    return sum(a != b for a, b in zip(left, right, strict=True))


def _residue_frequencies(peptides: Sequence[str]) -> tuple[Counter[str], dict[str, float]]:
    counts = Counter("".join(peptides))
    total = sum(counts.values())
    frequencies = {
        residue: counts[residue] / total if total else 0.0 for residue in AMINO_ACIDS
    }
    return counts, frequencies


def _jensen_shannon_divergence(
    frequencies: Mapping[str, float], reference: Mapping[str, float]
) -> float:
    """Return Jensen-Shannon divergence in bits for two residue distributions."""
    left = np.array([frequencies[residue] for residue in AMINO_ACIDS], dtype=float)
    right = np.array([reference[residue] for residue in AMINO_ACIDS], dtype=float)
    midpoint = 0.5 * (left + right)

    def divergence(values: np.ndarray) -> float:
        observed = values > 0
        return float(np.sum(values[observed] * np.log2(values[observed] / midpoint[observed])))

    return 0.5 * (divergence(left) + divergence(right))


def build_composition_audit(
    ranked: pd.DataFrame,
    shortlist: pd.DataFrame,
    observed_library: pd.DataFrame,
    *,
    quality_quantile: float = 0.95,
) -> pd.DataFrame:
    """Trace residue composition through each stage of the strict screening funnel."""
    required_ranked = {
        "AA",
        "passes_packaging_gate",
        "is_pareto",
        "display_score",
        "training_distance_lower_bound",
    }
    missing = required_ranked - set(ranked.columns)
    if missing:
        raise ValueError(f"Missing ranked columns: {sorted(missing)}")
    if "AA" not in shortlist or "AA" not in observed_library:
        raise ValueError("shortlist and observed_library must contain an AA column")

    eligible = ranked.loc[
        ranked["passes_packaging_gate"] & ranked["training_distance_lower_bound"].ge(2)
    ]
    cutoff = float(eligible["display_score"].quantile(quality_quantile))
    stages = {
        "observed_100k_library": observed_library["AA"].astype(str).tolist(),
        "virtual_pool": ranked["AA"].astype(str).tolist(),
        "packaging_95_lcb": ranked.loc[ranked["passes_packaging_gate"], "AA"]
        .astype(str)
        .tolist(),
        "pareto_after_95_lcb": ranked.loc[ranked["is_pareto"], "AA"].astype(str).tolist(),
        "quality_top5_after_95_lcb": eligible.loc[eligible["display_score"].ge(cutoff), "AA"]
        .astype(str)
        .tolist(),
        "final_shortlist": shortlist["AA"].astype(str).tolist(),
    }
    _, reference = _residue_frequencies(stages["virtual_pool"])
    rows: list[dict[str, object]] = []
    for stage, peptides in stages.items():
        counts, frequencies = _residue_frequencies(peptides)
        divergence = _jensen_shannon_divergence(frequencies, reference)
        for residue in AMINO_ACIDS:
            reference_frequency = reference[residue]
            rows.append(
                {
                    "stage": stage,
                    "sequence_count": len(peptides),
                    "residue": residue,
                    "residue_count": counts[residue],
                    "frequency": frequencies[residue],
                    "fold_vs_virtual_pool": (
                        frequencies[residue] / reference_frequency
                        if reference_frequency
                        else np.nan
                    ),
                    "js_divergence_bits_vs_virtual_pool": divergence,
                    "missing_from_stage": counts[residue] == 0,
                }
            )
    return pd.DataFrame(rows)


def build_gate_sensitivity_audit(
    ranked: pd.DataFrame,
    *,
    packaging_threshold: float,
    lower_bound_offsets: Mapping[str, float],
    strict_gate: str = "strict_95_lcb",
) -> pd.DataFrame:
    """Compare predeclared packaging gates without changing the primary shortlist."""
    if not {"AA", "pred_pack"} <= set(ranked.columns):
        raise ValueError("ranked must contain AA and pred_pack columns")
    if strict_gate not in lower_bound_offsets:
        raise ValueError(f"strict_gate {strict_gate!r} is not present in lower_bound_offsets")

    pool_peptides = ranked["AA"].astype(str).tolist()
    _, pool_frequencies = _residue_frequencies(pool_peptides)
    masks = {
        name: ranked["pred_pack"].sub(offset).ge(packaging_threshold).to_numpy()
        for name, offset in lower_bound_offsets.items()
    }
    strict_mask = masks[strict_gate]
    strict_count = int(strict_mask.sum())
    rows: list[dict[str, object]] = []
    for name, offset in lower_bound_offsets.items():
        mask = masks[name]
        passed = int(mask.sum())
        peptides = ranked.loc[mask, "AA"].astype(str).tolist()
        _, frequencies = _residue_frequencies(peptides)
        intersection = int(np.logical_and(mask, strict_mask).sum())
        union = int(np.logical_or(mask, strict_mask).sum())
        rows.append(
            {
                "gate": name,
                "lower_bound_offset": float(offset),
                "packaging_threshold": packaging_threshold,
                "passed": passed,
                "pass_rate": passed / len(ranked),
                "fold_count_vs_strict": passed / strict_count if strict_count else np.nan,
                "strict_candidates_retained": (
                    intersection / strict_count if strict_count else np.nan
                ),
                "jaccard_vs_strict": intersection / union if union else np.nan,
                "js_divergence_bits_vs_virtual_pool": _jensen_shannon_divergence(
                    frequencies, pool_frequencies
                ),
                "missing_residues": "".join(
                    residue for residue in AMINO_ACIDS if frequencies[residue] == 0
                ),
                "role": "primary_pre_registered" if name == strict_gate else "sensitivity_only",
            }
        )
    return pd.DataFrame(rows)


def summarize_funnel_audit(
    gate_audit: pd.DataFrame, composition_audit: pd.DataFrame
) -> dict[str, object]:
    """Build a compact, machine-readable interpretation of both audits."""
    gates = {
        str(row["gate"]): {
            "passed": int(row["passed"]),
            "pass_rate": float(row["pass_rate"]),
            "fold_count_vs_strict": float(row["fold_count_vs_strict"]),
            "role": str(row["role"]),
        }
        for row in gate_audit.to_dict(orient="records")
    }
    stage_summary: dict[str, object] = {}
    for stage, group in composition_audit.groupby("stage", sort=False):
        stage_summary[str(stage)] = {
            "sequence_count": int(group["sequence_count"].iloc[0]),
            "js_divergence_bits_vs_virtual_pool": float(
                group["js_divergence_bits_vs_virtual_pool"].iloc[0]
            ),
            "missing_residues": group.loc[group["missing_from_stage"], "residue"].tolist(),
        }
    strict = gates.get("strict_95_lcb", {})
    point = gates.get("sensitivity_point_prediction", {})
    return {
        "packaging_gates": gates,
        "composition_stages": stage_summary,
        "interpretation": {
            "strict_gate_is_high_false_negative_risk": bool(
                strict and point and point["passed"] >= 5 * strict["passed"]
            ),
            "composition_bias_first_visible_at": "packaging_95_lcb",
            "primary_shortlist_changed": False,
            "reason": (
                "Sensitivity gates diagnose robustness; they do not replace the pre-registered "
                "95% lower-bound shortlist without new experimental evidence."
            ),
        },
        "claim_boundary": (
            "This audit measures computational funnel behavior, not packaging, motor-neuron "
            "transduction, biodistribution, safety, or efficacy in wet-lab experiments."
        ),
    }


def select_composition_challenge_panel(
    ranked: pd.DataFrame,
    shortlist: pd.DataFrame,
    *,
    packaging_threshold: float,
    strict_95_offset: float,
    sensitivity_90_offset: float,
    target_residues: Sequence[str] = DEFAULT_COMPOSITION_CHALLENGE_RESIDUES,
    minimum_pairwise_distance: int = 3,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Select a separate packaging-only panel that challenges composition bias.

    For every target residue, one candidate is taken from the strict 95% lower-bound
    survivors and one from the 90%-only near-miss band. Candidates contain exactly
    one residue from the target set, which makes the composition challenge easier to
    interpret. This panel never changes or replaces the primary therapeutic shortlist.
    """
    required = {
        "variant_id",
        "AA",
        "pred_pack",
        "pred_pack_lcb",
        "display_score",
        "organ_uncertainty_mean",
        "training_distance_lower_bound",
    }
    missing = required - set(ranked.columns)
    if missing:
        raise ValueError(f"Missing ranked columns: {sorted(missing)}")
    if "AA" not in shortlist:
        raise ValueError("shortlist must contain an AA column")
    residues = tuple(dict.fromkeys(str(residue) for residue in target_residues))
    if not residues or any(residue not in AMINO_ACIDS for residue in residues):
        raise ValueError("target_residues must be non-empty standard amino-acid codes")
    if minimum_pairwise_distance < 1:
        raise ValueError("minimum_pairwise_distance must be positive")

    primary_sequences = set(shortlist["AA"].dropna().astype(str))
    candidates = ranked.loc[
        ranked["training_distance_lower_bound"].ge(2)
        & ~ranked["AA"].astype(str).isin(primary_sequences)
    ].copy()
    candidates["pred_pack_lcb_95"] = candidates["pred_pack"] - strict_95_offset
    candidates["pred_pack_lcb_90"] = candidates["pred_pack"] - sensitivity_90_offset
    candidates["packaging_margin_95"] = (
        candidates["pred_pack_lcb_95"] - packaging_threshold
    )
    candidates["packaging_threshold"] = packaging_threshold
    target_set = set(residues)
    candidates["challenge_residue_count"] = candidates["AA"].astype(str).map(
        lambda peptide: sum(residue in target_set for residue in peptide)
    )
    candidates = candidates.loc[candidates["challenge_residue_count"].eq(1)].copy()

    selected_rows: list[pd.Series] = []
    selected_sequences: list[str] = []
    strata = ("strict_95_survivor", "sensitivity_90_only")
    for residue in residues:
        contains_residue = candidates["AA"].astype(str).str.contains(residue, regex=False)
        for stratum in strata:
            if stratum == "strict_95_survivor":
                pool = candidates.loc[
                    contains_residue
                    & candidates["pred_pack_lcb_95"].ge(packaging_threshold)
                ].sort_values(
                    ["display_score", "organ_uncertainty_mean", "AA"],
                    ascending=[False, True, True],
                    kind="stable",
                )
            else:
                pool = candidates.loc[
                    contains_residue
                    & candidates["pred_pack_lcb_90"].ge(packaging_threshold)
                    & candidates["pred_pack_lcb_95"].lt(packaging_threshold)
                ].sort_values(
                    ["packaging_margin_95", "display_score", "organ_uncertainty_mean", "AA"],
                    ascending=[False, False, True, True],
                    kind="stable",
                )
            chosen = None
            for _, row in pool.iterrows():
                peptide = str(row["AA"])
                if all(
                    _hamming_distance(peptide, previous) >= minimum_pairwise_distance
                    for previous in selected_sequences
                ):
                    chosen = row.copy()
                    break
            if chosen is None:
                raise ValueError(
                    f"No composition challenge candidate available for {residue}/{stratum} "
                    f"at pairwise distance {minimum_pairwise_distance}"
                )
            chosen["target_residue"] = residue
            chosen["gate_stratum"] = stratum
            chosen["challenge_id"] = (
                f"COMP-{residue}-95" if stratum == "strict_95_survivor" else f"COMP-{residue}-90"
            )
            chosen["experimental_scope"] = "stage_1_packaging_and_qc_only"
            chosen["therapeutic_priority"] = "not_assigned_composition_challenge"
            chosen["replaces_primary_shortlist"] = False
            selected_rows.append(chosen)
            selected_sequences.append(str(chosen["AA"]))

    panel = pd.DataFrame(selected_rows)
    leading_columns = [
        "challenge_id",
        "target_residue",
        "gate_stratum",
        "variant_id",
        "AA",
        "pred_pack",
        "pred_pack_lcb_90",
        "pred_pack_lcb_95",
        "packaging_threshold",
        "packaging_margin_95",
        "display_score",
        "f_cns",
        "f_liv",
        "f_off",
        "organ_uncertainty_mean",
        "training_distance_lower_bound",
        "human_liver_warning",
        "experimental_scope",
        "therapeutic_priority",
        "replaces_primary_shortlist",
    ]
    output_columns = [column for column in leading_columns if column in panel]
    panel = panel.loc[:, output_columns]
    minimum_observed_distance = min(
        (
            _hamming_distance(left, right)
            for index, left in enumerate(selected_sequences)
            for right in selected_sequences[index + 1 :]
        ),
        default=0,
    )
    summary: dict[str, object] = {
        "panel_role": "separate_packaging_composition_challenge",
        "rows": len(panel),
        "target_residues": list(residues),
        "rows_per_residue": panel["target_residue"].value_counts().sort_index().to_dict(),
        "rows_per_gate_stratum": panel["gate_stratum"].value_counts().to_dict(),
        "minimum_pairwise_hamming_required": minimum_pairwise_distance,
        "minimum_pairwise_hamming_observed": minimum_observed_distance,
        "packaging_threshold": packaging_threshold,
        "strict_95_lower_bound_offset": strict_95_offset,
        "sensitivity_90_lower_bound_offset": sensitivity_90_offset,
        "primary_shortlist_changed": False,
        "selection_rule": (
            "One strict-95 survivor and one 90-only near miss per residue; exactly one target "
            "residue per sequence; distance >=2 from training and pairwise challenge distance "
            f">={minimum_pairwise_distance}."
        ),
        "claim_boundary": (
            "Packaging/QC challenge only. These rows are not therapeutic candidates and must "
            "not be pooled with the primary 30-candidate hit-rate analysis."
        ),
    }
    return panel, summary
