import numpy as np
import pandas as pd
import pytest

from aav9_sma.data.motor_neuron import (
    aggregate_motor_neuron_evidence,
    audit_motor_neuron_evidence,
    normalize_motor_neuron_evidence,
    study_level_split,
)
from aav9_sma.models.motor_neuron import MotorNeuronHead, fit_motor_neuron_head


def _evidence_frame() -> pd.DataFrame:
    rows = []
    peptides = ("AAAAAAA", "CCCCCCC", "DDDDDDD", "EEEEEEE", "FFFFFFF", "GGGGGGG")
    for index, peptide in enumerate(peptides):
        study = f"study_{index // 2 + 1}"
        rows.append(
            {
                "study_id": study,
                "animal_id": f"animal_{index}",
                "variant_id": f"variant_{index}",
                "peptide_7mer": peptide,
                "species": "mouse",
                "strain": "C57BL/6J",
                "route": "intravenous",
                "dose": 1.0e13,
                "timepoint": 14,
                "payload": "GFP",
                "promoter_or_enhancer": "CAG",
                "tissue": "spinal_cord",
                "cell_type": "spinal motor neuron",
                "readout_type": "reporter",
                "readout_value": float(index),
                "evidence_level": "L3",
            }
        )
        rows.append(
            {
                **rows[-1],
                "cell_type": "astrocyte",
                "readout_value": float(index) / 10,
            }
        )
    return pd.DataFrame(rows)


def test_normalizes_and_audits_cell_type_evidence() -> None:
    frame = _evidence_frame()
    normalized = normalize_motor_neuron_evidence(frame, require_sequence=True)
    audit = audit_motor_neuron_evidence(frame)

    assert normalized["is_motor_neuron"].sum() == 6
    assert normalized["direct_cell_level_evidence"].all()
    assert audit.row_count == 12
    assert audit.motor_neuron_rows == 6
    assert audit.sequence_linked_rows == 12
    assert audit.studies == ("study_1", "study_2", "study_3")


def test_rejects_bulk_spinal_cord_as_motor_neuron_training_data() -> None:
    frame = _evidence_frame().iloc[[0]].copy()
    frame["cell_type"] = "spinal cord bulk"
    frame["evidence_level"] = "L1"

    with pytest.raises(ValueError, match="No motor-neuron"):
        normalize_motor_neuron_evidence(
            frame,
            require_sequence=True,
            require_motor_labels=True,
        )


def test_rejects_conflicting_sequences_for_one_variant() -> None:
    frame = _evidence_frame().iloc[:2].copy()
    frame.loc[frame.index[1], "peptide_7mer"] = "CCCCCCC"

    with pytest.raises(ValueError, match="exactly one peptide_7mer"):
        normalize_motor_neuron_evidence(frame, require_sequence=True)


def test_rejects_mixed_raw_readout_scales() -> None:
    frame = _evidence_frame()
    frame.loc[frame["study_id"].eq("study_3"), "readout_type"] = "function"

    with pytest.raises(ValueError, match="one comparable readout_type"):
        fit_motor_neuron_head(frame, validation_fraction=0.34, random_state=3)


def test_aggregates_traceable_evidence_without_fabricating_labels() -> None:
    frame = _evidence_frame()
    frame.loc[0, "evidence_level"] = "L4"
    frame.loc[0, "readout_type"] = "function"

    aggregate = aggregate_motor_neuron_evidence(frame)

    first = aggregate.loc[aggregate["variant_id"] == "variant_0"].iloc[0]
    assert first["max_evidence_level"] == "L4"
    assert first["motor_neuron_evidence_score"] == 1.0
    assert aggregate["direct_cell_level_evidence"].all()


def test_study_split_never_leaks_a_study() -> None:
    frame = _evidence_frame()
    train_rows, validation_rows = study_level_split(
        frame,
        validation_fraction=0.34,
        random_state=7,
    )

    train_studies = set(frame.loc[train_rows, "study_id"])
    validation_studies = set(frame.loc[validation_rows, "study_id"])
    assert train_studies
    assert validation_studies
    assert train_studies.isdisjoint(validation_studies)


def test_fits_motor_neuron_head_and_predicts_new_variants(tmp_path) -> None:
    frame = _evidence_frame()
    result = fit_motor_neuron_head(
        frame,
        model_name="ridge",
        validation_fraction=0.34,
        random_state=3,
    )
    predictions = result.model.predict(
        pd.DataFrame(
            {
                "peptide_7mer": ["ACDEFGH", "YYYYYYY"],
                "species": ["mouse", "macaque"],
                "route": ["intravenous", "intrathecal"],
                "dose": [1.0e13, 5.0e12],
                "timepoint": [14, 28],
                "payload": ["GFP", "SMN1"],
                "promoter_or_enhancer": ["CAG", "Syn1"],
            }
        )
    )
    path = tmp_path / "motor_neuron_head.pkl"
    result.model.save(path)
    loaded = MotorNeuronHead.load(path)

    assert predictions.shape == (2,)
    assert np.isfinite(predictions).all()
    np.testing.assert_allclose(
        predictions,
        loaded.predict(
            pd.DataFrame(
                {
                    "peptide_7mer": ["ACDEFGH", "YYYYYYY"],
                    "species": ["mouse", "macaque"],
                    "route": ["intravenous", "intrathecal"],
                    "dose": [1.0e13, 5.0e12],
                    "timepoint": [14, 28],
                    "payload": ["GFP", "SMN1"],
                    "promoter_or_enhancer": ["CAG", "Syn1"],
                }
            )
        ),
    )
    assert result.training_rows > 0
    assert result.validation_rows > 0
    assert result.held_out_studies
