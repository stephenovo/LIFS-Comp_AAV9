#!/usr/bin/env python3
"""Prepare existing project data for the motor-neuron module.

The current repository contains organ-level virtual-screen predictions, not
cell-type-resolved motor-neuron labels. This script imports those candidates
as explicitly marked L1 spinal-cord proxy rows so they can be carried through
the future pipeline without being used as motor-neuron training labels.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def prepare_existing_candidates(input_path: str | Path) -> pd.DataFrame:
    """Map an existing shortlist/panel into the motor-neuron context schema."""
    project_root = Path(__file__).resolve().parents[1]
    if str(project_root / "src") not in sys.path:
        sys.path.insert(0, str(project_root / "src"))
    from aav9_sma.data.motor_neuron import normalize_motor_neuron_evidence

    source = pd.read_csv(input_path)
    peptide_column = "peptide_7mer" if "peptide_7mer" in source else "AA"
    if peptide_column not in source:
        raise ValueError("Input must contain AA or peptide_7mer")
    variant_column = "variant_id" if "variant_id" in source else None
    if variant_column is None:
        source["variant_id"] = [f"imported_{index:06d}" for index in range(len(source))]
        variant_column = "variant_id"
    if "pred_spinal_cord_mouse" not in source:
        raise ValueError("Input must contain pred_spinal_cord_mouse")
    source = source[pd.to_numeric(source["pred_spinal_cord_mouse"], errors="coerce").notna()].copy()
    if source.empty:
        raise ValueError("Input has no finite pred_spinal_cord_mouse rows")

    output = pd.DataFrame(
        {
            "study_id": "fit4function_existing_virtual_screen",
            "animal_id": pd.NA,
            "variant_id": source[variant_column].astype("string"),
            "peptide_7mer": source[peptide_column].astype("string").str.strip().str.upper(),
            "species": "mouse",
            "strain": pd.NA,
            "route": "unknown",
            "dose": pd.NA,
            "timepoint": pd.NA,
            "payload": "unknown",
            "promoter_or_enhancer": "unknown",
            "tissue": "spinal_cord",
            "cell_type": "spinal_cord_bulk_proxy",
            "cell_subtype": pd.NA,
            "readout_type": "predicted_spinal_cord_proxy",
            "readout_value": pd.to_numeric(
                source["pred_spinal_cord_mouse"], errors="coerce"
            ),
            "evidence_level": "L1",
            "source_file": str(Path(input_path)),
            "source_status": "organ_proxy_not_motor_neuron",
        }
    )
    normalized = normalize_motor_neuron_evidence(output, require_sequence=True)
    normalized["eligible_for_motor_neuron_training"] = False
    normalized["import_note"] = (
        "Existing spinal-cord proxy; not a motor-neuron transduction label."
    )
    return normalized


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prepared = prepare_existing_candidates(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    prepared.to_csv(args.output, index=False)
    print(
        f"wrote {len(prepared)} rows; "
        f"motor_neuron_rows={int(prepared['is_motor_neuron'].sum())}; "
        f"training_eligible={int(prepared['eligible_for_motor_neuron_training'].sum())}"
    )


if __name__ == "__main__":
    main()
