#!/usr/bin/env python3
"""Plot cross-animal model robustness against empirical replicate agreement."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ENDPOINTS = ("brain", "spinal_cord", "liver", "heart", "kidney")
ENDPOINT_LABELS = ("Brain", "Spinal cord", "Liver", "Heart", "Kidney")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cross-animal", type=Path, required=True)
    parser.add_argument("--replicate-metrics", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    cross = pd.read_csv(args.cross_animal)
    replicate = pd.read_csv(args.replicate_metrics)
    metric = "model_vs_heldout_pearson_r"
    matrix = (
        cross.pivot(index="task", columns="held_out_animal", values=metric)
        .reindex(ENDPOINTS)
        .reindex(columns=(1, 2, 3, 4))
    )
    repeatability = (
        replicate.loc[replicate["denominator_mode"].eq("whitelist")]
        .set_index("endpoint")
        .reindex(ENDPOINTS)["median_pairwise_animal_pearson"]
    )

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
        }
    )
    figure, axes = plt.subplots(1, 2, figsize=(12.2, 4.8), constrained_layout=True)

    image = axes[0].imshow(matrix.to_numpy(), cmap="cividis", vmin=0.0, vmax=0.85)
    axes[0].set_xticks(np.arange(4), ["Animal 1", "Animal 2", "Animal 3", "Animal 4"])
    axes[0].set_yticks(np.arange(len(ENDPOINT_LABELS)), ENDPOINT_LABELS)
    axes[0].set_xlabel("Completely held-out animal")
    axes[0].set_title("A. Leave-one-animal-out Pearson r")
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            value = matrix.iloc[row, column]
            text_color = "white" if value < 0.28 or value > 0.68 else "black"
            axes[0].text(
                column,
                row,
                f"{value:.2f}",
                ha="center",
                va="center",
                color=text_color,
                fontweight="bold",
            )
    colorbar = figure.colorbar(image, ax=axes[0], fraction=0.046, pad=0.04)
    colorbar.set_label("Pearson r")

    values = matrix.to_numpy()
    means = np.nanmean(values, axis=1)
    lower = means - np.nanmin(values, axis=1)
    upper = np.nanmax(values, axis=1) - means
    y = np.arange(len(ENDPOINTS))
    axes[1].errorbar(
        means,
        y,
        xerr=np.vstack([lower, upper]),
        fmt="o",
        color="#2855A6",
        ecolor="#7AA6C2",
        capsize=4,
        label="Model mean + animal range",
    )
    axes[1].scatter(
        repeatability.to_numpy(),
        y,
        marker="D",
        s=42,
        color="#D78B22",
        label="Median animal-pair agreement",
        zorder=3,
    )
    axes[1].set_yticks(y, ENDPOINT_LABELS)
    axes[1].invert_yaxis()
    axes[1].set_xlim(0.2, 0.9)
    axes[1].set_xlabel("Pearson r")
    axes[1].set_title("B. Robustness relative to replicate agreement")
    axes[1].grid(axis="x", alpha=0.2)
    axes[1].legend(
        frameon=False,
        fontsize=9,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
    )

    figure.suptitle(
        "Cross-animal robustness of the five-model shared-MLP ensemble",
        fontsize=14,
        fontweight="bold",
    )
    figure.text(
        0.5,
        -0.02,
        "Development audit only; agreement is not a wet-lab efficacy estimate.",
        ha="center",
        fontsize=9,
        color="#475569",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(figure)


if __name__ == "__main__":
    main()
