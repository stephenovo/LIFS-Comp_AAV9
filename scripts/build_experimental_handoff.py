"""Build wet-lab candidate, result-entry, and decision-log CSV files."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from aav9_sma.experimental.handoff import build_validation_panel

RESULT_COLUMNS = [
    "study_id",
    "stage",
    "sample_id",
    "blind_id",
    "production_batch",
    "biological_replicate",
    "technical_replicate",
    "assay",
    "cell_type_or_tissue",
    "timepoint",
    "dose_value",
    "dose_unit",
    "raw_value",
    "raw_unit",
    "normalized_value",
    "normalized_unit",
    "normalization_reference",
    "qc_status",
    "included_in_primary_analysis",
    "exclusion_reason",
    "operator_blinded",
    "notes",
]

DECISION_COLUMNS = [
    "decision_date",
    "stage",
    "sample_id",
    "decision",
    "primary_endpoint_result",
    "packaging_gate_result",
    "target_signal_result",
    "liver_signal_result",
    "safety_result",
    "evidence_location",
    "decision_maker",
    "deviation_from_preregistered_rule",
    "deviation_reason",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shortlist", type=Path, required=True)
    parser.add_argument("--controls", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    shortlist = pd.read_csv(args.shortlist)
    controls = pd.read_csv(args.controls)
    panel = build_validation_panel(shortlist, controls)
    panel.to_csv(args.output_dir / "experimental_validation_panel.csv", index=False)
    pd.DataFrame(columns=RESULT_COLUMNS).to_csv(
        args.output_dir / "wet_lab_results_template.csv", index=False
    )
    pd.DataFrame(columns=DECISION_COLUMNS).to_csv(
        args.output_dir / "wet_lab_decision_log.csv", index=False
    )
    manifest = {
        "schema_version": "1.0",
        "shortlist_source": str(args.shortlist),
        "shortlist_sha256": hashlib.sha256(args.shortlist.read_bytes()).hexdigest(),
        "control_source": str(args.controls),
        "control_sha256": hashlib.sha256(args.controls.read_bytes()).hexdigest(),
        "panel_rows": int(len(panel)),
        "computational_candidates": int((panel["panel_role"] == "computational_candidate").sum()),
        "controls": int((panel["priority_tier"] == "Control").sum()),
        "tier_counts": {
            tier: int(count)
            for tier, count in panel.loc[
                panel["panel_role"] == "computational_candidate", "priority_tier"
            ]
            .value_counts()
            .sort_index()
            .items()
        },
    }
    (args.output_dir / "experimental_validation_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
