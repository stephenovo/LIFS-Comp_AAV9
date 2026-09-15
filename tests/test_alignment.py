from aav9_sma.data.alignment import (
    INSERTION_REFERENCE_START,
    REFERENCE_PREFIX,
    REFERENCE_SUFFIX,
    extract_aligned_insertion,
)


def test_extracts_insert_from_full_length_alignment() -> None:
    insertion = "GCTTGTGATGAATTTGGTCAT"
    sequence = REFERENCE_PREFIX + insertion + REFERENCE_SUFFIX

    status, nucleotide, peptide = extract_aligned_insertion(
        sequence=sequence,
        quality="I" * len(sequence),
        reference_start_one_based=1,
        cigar=f"{len(sequence)}M",
    )

    assert status == "valid"
    assert nucleotide == insertion
    assert peptide == "ACDEFGH"


def test_rejects_deletion_overlapping_insert() -> None:
    sequence = "A" * (INSERTION_REFERENCE_START + 20)
    status, _, _ = extract_aligned_insertion(
        sequence=sequence,
        quality="I" * len(sequence),
        reference_start_one_based=1,
        cigar=f"{INSERTION_REFERENCE_START + 10}M1D10M",
    )

    assert status == "indel_or_incomplete"
