import pandas as pd

from aav9_sma.data.audit import audit_dataframe


def test_audit_reports_schema_and_peptide_quality() -> None:
    frame = pd.DataFrame(
        {
            "variant_id": ["v1", "v2"],
            "peptide_7mer": ["ACDEFGH", "ACDEXGH"],
            "f_pack": [0.8, 0.7],
            "brain_mouse": [0.6, None],
            "spinal_cord_mouse": [0.5, 0.4],
            "liver_mouse": [0.2, 0.3],
            "heart_mouse": [0.1, 0.2],
            "kidney_mouse": [0.1, 0.2],
        }
    )

    result = audit_dataframe(frame)

    assert result.row_count == 2
    assert result.invalid_peptide_rows == 1
    assert result.complete_primary_rows == 1
    assert result.missing_columns == ()
