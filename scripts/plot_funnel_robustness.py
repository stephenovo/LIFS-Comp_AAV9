#!/usr/bin/env python3
"""Plot packaging-gate sensitivity and sequence-composition drift."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

GATE_ORDER = (
    "strict_95_lcb",
    "sensitivity_90_lcb",
    "sensitivity_point_prediction",
)
GATE_LABELS = ("95% lower bound\n(primary)", "90% lower bound", "Point prediction")
STAGE_ORDER = (
    "virtual_pool",
    "packaging_95_lcb",
    "pareto_after_95_lcb",
    "quality_top5_after_95_lcb",
    "final_shortlist",
)
STAGE_LABELS = ("Virtual\npool", "Packaging\ngate", "Pareto", "Top 5%\nquality", "Final\n30")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gate-audit", type=Path, required=True)
    parser.add_argument("--composition-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    gates = pd.read_csv(args.gate_audit).set_index("gate").reindex(GATE_ORDER)
    composition = pd.read_csv(args.composition_audit)
    stages = (
        composition.groupby("stage", sort=False)
        .agg(
            sequence_count=("sequence_count", "first"),
            js_divergence=("js_divergence_bits_vs_virtual_pool", "first"),
            missing_residues=("missing_from_stage", "sum"),
        )
        .reindex(STAGE_ORDER)
    )

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
        }
    )
    figure, axes = plt.subplots(1, 2, figsize=(11.2, 4.8), constrained_layout=True)

    counts = gates["passed"].to_numpy()
    colors = ["#2855A6", "#7AA6C2", "#CBD5E1"]
    bars = axes[0].bar(np.arange(len(gates)), counts, color=colors)
    axes[0].set_yscale("log")
    axes[0].set_xticks(np.arange(len(gates)), GATE_LABELS)
    axes[0].set_ylabel("Candidates passing gate (log scale)")
    axes[0].set_ylim(5_000, counts.max() * 1.7)
    axes[0].set_title("A. Packaging uncertainty changes pool size")
    axes[0].grid(axis="y", alpha=0.2)
    for bar, count, pass_rate in zip(bars, counts, gates["pass_rate"], strict=True):
        axes[0].text(
            bar.get_x() + bar.get_width() / 2,
            count * 1.08,
            f"{count:,}\n({pass_rate:.2%})",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    x = np.arange(len(stages))
    divergence = stages["js_divergence"].to_numpy()
    sizes = 55 + 22 * np.log10(stages["sequence_count"].to_numpy())
    axes[1].plot(x, divergence, color="#7AA6C2", linewidth=1.8, zorder=1)
    axes[1].scatter(x, divergence, s=sizes, color="#D78B22", zorder=2)
    axes[1].set_xticks(x, STAGE_LABELS)
    axes[1].set_ylabel("Jensen–Shannon divergence (bits)")
    axes[1].set_ylim(-0.015, max(divergence) + 0.075)
    axes[1].set_title("B. Composition drift through the funnel")
    axes[1].grid(axis="y", alpha=0.2)
    for index, row in enumerate(stages.itertuples()):
        missing = int(row.missing_residues)
        note = f"n={int(row.sequence_count):,}"
        if missing:
            noun = "residue" if missing == 1 else "residues"
            note += f"\n{missing} {noun} absent"
        axes[1].annotate(
            note,
            (index, row.js_divergence),
            xytext=(0, 12),
            textcoords="offset points",
            ha="center",
            fontsize=8,
        )

    figure.suptitle(
        "Virtual-screen funnel robustness and composition pressure",
        fontsize=14,
        fontweight="bold",
    )
    figure.text(
        0.5,
        -0.02,
        "Sensitivity analyses diagnose the frozen funnel; "
        "they do not replace its primary shortlist.",
        ha="center",
        fontsize=9,
        color="#475569",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(figure)


if __name__ == "__main__":
    main()
