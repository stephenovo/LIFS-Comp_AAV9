"""Small deterministic smoke test that runs without external research data."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from aav9_sma.constants import AMINO_ACIDS, PRIMARY_LABEL_COLUMNS
from aav9_sma.data.audit import audit_dataframe
from aav9_sma.screening.score import rank_candidates


def _demo_peptides(count: int, random_state: int) -> list[str]:
    rng = np.random.default_rng(random_state)
    peptides: set[str] = set()
    while len(peptides) < count:
        peptides.add("".join(rng.choice(list(AMINO_ACIDS), size=7)))
    return sorted(peptides)


def _demo_frame(count: int = 64, random_state: int = 42) -> pd.DataFrame:
    """Create a synthetic canonical table for installation and CLI smoke tests."""
    rng = np.random.default_rng(random_state)
    peptides = _demo_peptides(count, random_state)
    charges = np.array(
        [sum(residue in "KRH" for residue in peptide) - sum(residue in "DE" for residue in peptide)
         for peptide in peptides],
        dtype=float,
    )
    hydrophobic = np.array(
        [sum(residue in "AVILMFWY" for residue in peptide) for peptide in peptides],
        dtype=float,
    )
    latent = charges + 0.25 * hydrophobic
    noise = rng.normal(0.0, 0.08, size=(count, len(PRIMARY_LABEL_COLUMNS)))
    labels = np.column_stack(
        [
            0.8 + 0.05 * latent,
            0.4 + 0.12 * latent,
            0.35 + 0.10 * latent,
            -0.25 - 0.14 * latent,
            0.05 + 0.05 * latent,
            0.02 + 0.04 * latent,
        ]
    ) + noise
    frame = pd.DataFrame(
        {
            "variant_id": [f"demo_{index:04d}" for index in range(count)],
            "peptide_7mer": peptides,
        }
    )
    frame[list(PRIMARY_LABEL_COLUMNS)] = labels
    return frame


def run_demo(output_dir: str | Path, *, random_state: int = 42) -> dict[str, object]:
    """Run the self-contained demo and write auditable CSV/JSON outputs."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    frame = _demo_frame(random_state=random_state)
    input_path = output / "demo_canonical.csv"
    audit_path = output / "demo_audit.json"
    predictions_path = output / "demo_predictions.csv"
    ranked_path = output / "demo_ranked_candidates.csv"
    summary_path = output / "demo_summary.json"

    frame.to_csv(input_path, index=False)
    audit = audit_dataframe(frame).to_dict()
    audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    predictions = pd.DataFrame({"variant_id": frame["variant_id"], "AA": frame["peptide_7mer"]})
    for label, prediction in zip(
        PRIMARY_LABEL_COLUMNS,
        (
            "pred_pack",
            "pred_brain_mouse",
            "pred_spinal_cord_mouse",
            "pred_liver_mouse",
            "pred_heart_mouse",
            "pred_kidney_mouse",
        ),
        strict=True,
    ):
        predictions[prediction] = frame[label]
    predictions["pred_pack_lcb"] = predictions["pred_pack"] - 0.20
    predictions.to_csv(predictions_path, index=False)

    threshold = float(frame["f_pack"].median() - 0.20)
    ranked = rank_candidates(predictions, packaging_threshold=threshold)
    ranked.to_csv(ranked_path, index=False)
    summary = {
        "mode": "synthetic_smoke_test",
        "random_state": random_state,
        "rows": len(frame),
        "packaging_threshold": threshold,
        "passed_packaging_gate": int(ranked["passes_packaging_gate"].sum()),
        "outputs": {
            "input": str(input_path),
            "audit": str(audit_path),
            "predictions": str(predictions_path),
            "ranked": str(ranked_path),
        },
        "warning": "Synthetic data only; never use demo outputs for biological claims.",
    }
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return summary
