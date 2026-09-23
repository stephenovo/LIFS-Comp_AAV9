from pathlib import Path

import pandas as pd
import pytest

from scripts.prepare_motor_neuron_data import prepare_existing_candidates


def test_imports_existing_shortlist_as_non_motor_neuron_proxy(tmp_path: Path) -> None:
    source = tmp_path / "shortlist.csv"
    pd.DataFrame(
        {
            "variant_id": ["V1", "V2"],
            "AA": ["ACDEFGH", "YYYYYYY"],
            "pred_spinal_cord_mouse": [0.2, -0.1],
        }
    ).to_csv(source, index=False)

    output = prepare_existing_candidates(source)

    assert output["evidence_level"].tolist() == ["L1", "L1"]
    assert not output["is_motor_neuron"].any()
    assert not output["eligible_for_motor_neuron_training"].any()
    assert output["readout_type"].eq("predicted_spinal_cord_proxy").all()


def test_import_rejects_data_without_spinal_proxy(tmp_path: Path) -> None:
    source = tmp_path / "not_a_screen.csv"
    pd.DataFrame({"AA": ["ACDEFGH"]}).to_csv(source, index=False)

    with pytest.raises(ValueError, match="pred_spinal_cord_mouse"):
        prepare_existing_candidates(source)
