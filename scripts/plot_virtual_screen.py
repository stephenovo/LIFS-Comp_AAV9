#!/usr/bin/env python3
"""Create the tracked summary figure for the virtual-screen checkpoint."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ranked", type=Path, required=True)
    parser.add_argument("--shortlist", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--multitask", type=Path, required=True)
    parser.add_argument("--ensemble", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    ranked = pd.read_csv(
        args.ranked,
        usecols=["AA", "passes_packaging_gate", "is_pareto", "f_cns", "f_liv"],
    )
    eligible = ranked.loc[ranked["passes_packaging_gate"]]
    pareto = eligible.loc[eligible["is_pareto"]]
    shortlist = pd.read_csv(args.shortlist)
    baseline = pd.read_csv(args.baseline)
    multitask = pd.read_csv(args.multitask)
    ensemble = pd.read_csv(args.ensemble)
    summary = json.loads(args.summary.read_text(encoding="utf-8"))

    tasks = ["brain", "spinal_cord", "liver", "heart", "kidney"]
    labels = ["Brain", "Spinal cord", "Liver", "Heart", "Kidney"]
    baseline_values = [
        baseline.loc[baseline["task"].eq(task), "model_vs_animal4_pearson_r"].max()
        for task in tasks
    ]
    multitask_values = [
        multitask.loc[multitask["task"].eq(task), "model_vs_animal4_pearson_r"].iloc[0]
        for task in tasks
    ]
    ensemble_values = [
        ensemble.loc[ensemble["task"].eq(task), "model_vs_animal4_pearson_r"].iloc[0]
        for task in tasks
    ]

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
        }
    )
    figure, axes = plt.subplots(1, 3, figsize=(16, 5.2), constrained_layout=True)
    figure.suptitle("LIFS-Comp_AAV9 computational screening checkpoint", fontsize=16)

    x = np.arange(len(tasks))
    width = 0.24
    axes[0].bar(x - width, baseline_values, width, label="Best single-task", color="#A8B3C7")
    axes[0].bar(x, multitask_values, width, label="Shared MLP", color="#6C83B5")
    axes[0].bar(x + width, ensemble_values, width, label="5-model ensemble", color="#2855A6")
    axes[0].set_xticks(x, labels, rotation=20, ha="right")
    axes[0].set_ylim(0, 0.9)
    axes[0].set_ylabel("Pearson r against held-out Animal 4")
    axes[0].set_title("A. Model comparison")
    axes[0].grid(axis="y", alpha=0.2)
    axes[0].legend(frameon=False, fontsize=8)

    funnel_labels = ["Generated", "Packaging gate", "Strict Pareto", "Shortlist"]
    funnel_values = [
        summary["pool_size"],
        summary["packaging_gate_passed"],
        summary["pareto_candidates"],
        summary["shortlist_rows"],
    ]
    y = np.arange(len(funnel_labels))
    axes[1].barh(y, funnel_values, color=["#CBD5E1", "#7AA6C2", "#E0A458", "#2E7D5B"])
    axes[1].set_yticks(y, funnel_labels)
    axes[1].invert_yaxis()
    axes[1].set_xscale("log")
    axes[1].set_xlabel("Candidate count (log scale)")
    axes[1].set_title("B. Screening funnel")
    axes[1].grid(axis="x", alpha=0.2)
    for row, value in enumerate(funnel_values):
        axes[1].text(value * 1.12, row, f"{value:,}", va="center", fontsize=9)

    axes[2].scatter(
        eligible["f_liv"],
        eligible["f_cns"],
        s=8,
        alpha=0.16,
        color="#64748B",
        linewidths=0,
        label=f"Packaging eligible ({len(eligible):,})",
    )
    axes[2].scatter(
        pareto["f_liv"],
        pareto["f_cns"],
        s=14,
        alpha=0.75,
        color="#D78B22",
        linewidths=0,
        label=f"Strict Pareto ({len(pareto):,})",
    )
    group_styles = {
        "balanced": ("#2E7D5B", "o", "Balanced"),
        "cns_favoring": ("#7048A8", "^", "CNS-favoring"),
        "low_liver": ("#C43D4D", "s", "Low-liver"),
    }
    for group, (color, marker, label) in group_styles.items():
        selected = shortlist.loc[shortlist["selection_group"].eq(group)]
        axes[2].scatter(
            selected["f_liv"],
            selected["f_cns"],
            s=58,
            color=color,
            marker=marker,
            edgecolors="white",
            linewidths=0.7,
            label=label,
            zorder=5,
        )
    axes[2].axhline(0, color="#94A3B8", linewidth=0.8, linestyle="--")
    axes[2].axvline(0, color="#94A3B8", linewidth=0.8, linestyle="--")
    axes[2].set_xlabel("Predicted mouse liver enrichment (lower is preferred)")
    axes[2].set_ylabel("Predicted CNS enrichment (higher is preferred)")
    axes[2].set_title("C. Eligible prediction landscape")
    axes[2].grid(alpha=0.15)
    axes[2].legend(frameon=False, fontsize=7, loc="lower left")

    figure.text(
        0.5,
        -0.01,
        "Mouse-organ enrichment proxies and computational packaging bounds; "
        "not experimental validation.",
        ha="center",
        fontsize=9,
        color="#475569",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(figure)


if __name__ == "__main__":
    main()
