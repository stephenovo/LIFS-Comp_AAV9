from pathlib import Path

from aav9_sma.data.fastq import LEFT_ANCHOR, RIGHT_ANCHOR, count_fastq, extract_insertion


def test_extracts_and_translates_valid_insertion() -> None:
    insertion = "GCTTGTGATGAATTTGGTCAT"  # ACDEF GH
    sequence = f"PREFIX{LEFT_ANCHOR}{insertion}{RIGHT_ANCHOR}SUFFIX"
    quality = "I" * len(sequence)

    result = extract_insertion(sequence, quality)

    assert result.status == "valid"
    assert result.nucleotide_21mer == insertion
    assert result.peptide_7mer == "ACDEFGH"


def test_rejects_low_quality_insertion() -> None:
    insertion = "GCTTGTGATGAATTTGGTCAT"
    sequence = f"{LEFT_ANCHOR}{insertion}{RIGHT_ANCHOR}"
    quality = "I" * len(LEFT_ANCHOR) + "!" + "I" * (len(sequence) - len(LEFT_ANCHOR) - 1)

    assert extract_insertion(sequence, quality).status == "low_quality"


def test_counts_whitelist_matches(tmp_path: Path) -> None:
    insertion = "GCTTGTGATGAATTTGGTCAT"
    sequence = f"{LEFT_ANCHOR}{insertion}{RIGHT_ANCHOR}"
    fastq = tmp_path / "sample.fastq"
    fastq.write_text(f"@read\n{sequence}\n+\n{'I' * len(sequence)}\n", encoding="ascii")

    counts_path = tmp_path / "counts.csv"
    result = count_fastq(
        fastq,
        peptide_whitelist={"ACDEFGH", "AAAAAAA"},
        counts_output=counts_path,
    )

    assert result["whitelist_matched_reads"] == 1
    assert result["whitelist_detected_peptides"] == 1
    assert "ACDEFGH,1" in counts_path.read_text()
