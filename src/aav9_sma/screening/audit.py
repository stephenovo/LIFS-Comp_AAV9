"""Post-screen audits for packaging gates and amino-acid composition."""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from aav9_sma.constants import AMINO_ACIDS


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
