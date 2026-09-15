import json
from pathlib import Path

import numpy as np
import pandas as pd

from aav9_sma.data.reconstruct import (
    load_count_matrix,
    parse_hammerhead_alias,
    reconstruct_liver,
    reconstruct_multiorgan,
)


def test_parses_organ_and_virus_aliases() -> None:
    organ = parse_hammerhead_alias("hammerhead_Spinal_cord_a4_r3")
    compact_spinal = parse_hammerhead_alias("hammerhead_SpinalCord_a4_r1")
    virus = parse_hammerhead_alias("hammerhead_starter_virus_prod2_1")

    assert organ.endpoint == "Spinal cord"
    assert compact_spinal.endpoint == "Spinal cord"
    assert organ.biological_replicate == 4
    assert organ.technical_replicate == 3
    assert virus.kind == "virus"
    assert virus.production_round == 2


def _write_run(tmp_path: Path, run: str, counts: list[int]) -> None:
    counts_dir = tmp_path / "counts"
    summaries_dir = tmp_path / "summaries"
    counts_dir.mkdir(exist_ok=True)
    summaries_dir.mkdir(exist_ok=True)
    pd.DataFrame(
        {
            "peptide_7mer": ["AAAAAAA", "CCCCCCC", "DDDDDDD"],
            "raw_read_count": counts,
        }
    ).to_csv(counts_dir / f"{run}.counts.csv.gz", index=False)
    summary = {
        "reads_examined": 100,
        "status_counts": {"valid": 80},
        "valid_fraction": 0.8,
        "whitelist_matched_reads": sum(counts),
    }
    (summaries_dir / f"{run}.json").write_text(json.dumps(summary), encoding="utf-8")


def test_reconstructs_liver_enrichment_and_validation(tmp_path: Path) -> None:
    rows = []
    for technical in (1, 2, 3):
        liver_run = f"L{technical}"
        brain_run = f"B{technical}"
        virus_run = f"V{technical}"
        _write_run(tmp_path, liver_run, [40, 10, 20])
        _write_run(tmp_path, brain_run, [20, 20, 20])
        _write_run(tmp_path, virus_run, [10, 40, 20])
        rows.extend(
            [
                {
                    "run_accession": liver_run,
                    "experiment_alias": f"hammerhead_Liver_a1_r{technical}",
                },
                {
                    "run_accession": brain_run,
                    "experiment_alias": f"hammerhead_Brain_a1_r{technical}",
                },
                {
                    "run_accession": virus_run,
                    "experiment_alias": f"hammerhead_start_virus_prod3_r{technical}",
                },
            ]
        )
    manifest = pd.DataFrame(rows)
    counts, qc = load_count_matrix(manifest, tmp_path / "counts", tmp_path / "summaries")
    public = pd.DataFrame(
        {"AA": ["AAAAAAA", "CCCCCCC", "DDDDDDD"], "Liver": [2.0, -2.0, 0.0]}
    )

    reconstructed, metrics = reconstruct_liver(counts, qc, public)

    column = "rpm_whitelist__liver_all_animals__over__virus_prod3"
    assert np.allclose(reconstructed[column], [2.0, -2.0, 0.0])
    assert not metrics.empty

    multiorgan, animal_metrics = reconstruct_multiorgan(counts, qc, virus_round=3)
    multiorgan_column = "log2enr_whitelist__liver_a1__over__virus_prod3"
    assert np.allclose(multiorgan[multiorgan_column], [2.0, -2.0, 0.0])
    assert "log2enr_whitelist__brain_a1__over__virus_prod3" in multiorgan
    assert set(animal_metrics["endpoint"]) == {"brain", "liver"}
