"""Conservative pilot extraction of Fit4Function 21-nt insertions from FASTQ."""

from __future__ import annotations

import csv
import gzip
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import TextIO

# Short exact anchors adjacent to the insertion in the AAV9 reference reported
# in the Fit4Function Methods. The published pipeline used Bowtie2; this lighter
# extractor is intentionally conservative and is for pilot QC, not parity claims.
LEFT_ANCHOR = "AAACCACCAGAGTGCCCAA"
RIGHT_ANCHOR = "GCACAGGCGCAGACCGGT"
INSERTION_NT_LENGTH = 21

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


@dataclass(frozen=True)
class ExtractedInsertion:
    status: str
    nucleotide_21mer: str | None = None
    peptide_7mer: str | None = None


def translate_21mer(sequence: str) -> str:
    if len(sequence) != INSERTION_NT_LENGTH:
        raise ValueError("Expected a 21-nt sequence")
    return "".join(CODON_TABLE[sequence[index : index + 3]] for index in range(0, 21, 3))


def extract_insertion(sequence: str, quality: str, minimum_phred: int = 20) -> ExtractedInsertion:
    """Extract one insertion using exact adjacent anchors and Q20 filtering."""
    left_index = sequence.find(LEFT_ANCHOR)
    if left_index < 0:
        return ExtractedInsertion("failed_anchor")
    start = left_index + len(LEFT_ANCHOR)
    stop = start + INSERTION_NT_LENGTH
    if sequence[stop : stop + len(RIGHT_ANCHOR)] != RIGHT_ANCHOR:
        return ExtractedInsertion("failed_anchor")
    insertion = sequence[start:stop]
    insertion_quality = quality[start:stop]
    if len(insertion_quality) != INSERTION_NT_LENGTH:
        return ExtractedInsertion("failed_anchor")
    if "N" in insertion:
        return ExtractedInsertion("undetermined", insertion)
    if any((ord(character) - 33) < minimum_phred for character in insertion_quality):
        return ExtractedInsertion("low_quality", insertion)
    if any(base not in "ACGT" for base in insertion):
        return ExtractedInsertion("invalid_base", insertion)
    peptide = translate_21mer(insertion)
    if "*" in peptide:
        return ExtractedInsertion("stop_codon", insertion, peptide)
    return ExtractedInsertion("valid", insertion, peptide)


def _open_fastq(path: Path) -> TextIO:
    if path.suffix in {".gz", ".gzip"}:
        return gzip.open(path, "rt", encoding="ascii")
    return path.open("r", encoding="ascii")


def count_fastq(
    path: str | Path,
    max_reads: int | None = None,
    peptide_whitelist: set[str] | None = None,
    counts_output: str | Path | None = None,
) -> dict[str, object]:
    """Count extraction outcomes and unique variants in a FASTQ file."""
    path = Path(path)
    statuses: Counter[str] = Counter()
    nucleotides: Counter[str] = Counter()
    peptides: Counter[str] = Counter()
    whitelist_peptides: Counter[str] = Counter()
    reads = 0
    with _open_fastq(path) as stream:
        while max_reads is None or reads < max_reads:
            header = stream.readline()
            if not header:
                break
            sequence = stream.readline().rstrip("\n")
            separator = stream.readline()
            quality = stream.readline().rstrip("\n")
            if not separator or not quality:
                raise ValueError("Truncated FASTQ record")
            result = extract_insertion(sequence, quality)
            statuses[result.status] += 1
            if result.nucleotide_21mer is not None:
                nucleotides[result.nucleotide_21mer] += 1
            if result.peptide_7mer is not None and result.status == "valid":
                peptides[result.peptide_7mer] += 1
                if peptide_whitelist is not None and result.peptide_7mer in peptide_whitelist:
                    whitelist_peptides[result.peptide_7mer] += 1
            reads += 1

    valid_reads = statuses["valid"]
    output: dict[str, object] = {
        "path": str(path),
        "mode": "conservative-exact-anchor-pilot",
        "reads_examined": reads,
        "status_counts": dict(statuses),
        "valid_fraction": valid_reads / reads if reads else 0.0,
        "unique_valid_nucleotide_21mers": len(nucleotides),
        "unique_valid_peptide_7mers": len(peptides),
        "top_peptide_7mers": peptides.most_common(20),
        "method_note": (
            "The publication used Bowtie2 alignment and a synthetic-library whitelist. "
            "This pilot uses exact adjacent anchors and does not claim exact parity."
        ),
        "anchors": asdict(
            _AnchorMetadata(left=LEFT_ANCHOR, right=RIGHT_ANCHOR, insertion_nt=21)
        ),
    }
    if peptide_whitelist is not None:
        output["whitelist_size"] = len(peptide_whitelist)
        output["whitelist_matched_reads"] = sum(whitelist_peptides.values())
        output["whitelist_detected_peptides"] = len(whitelist_peptides)
        output["whitelist_detection_fraction"] = (
            len(whitelist_peptides) / len(peptide_whitelist) if peptide_whitelist else 0.0
        )
        output["top_whitelist_peptide_7mers"] = whitelist_peptides.most_common(20)
        if counts_output is not None:
            counts_path = Path(counts_output)
            counts_path.parent.mkdir(parents=True, exist_ok=True)
            with counts_path.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.writer(stream)
                writer.writerow(["peptide_7mer", "raw_read_count"])
                writer.writerows(
                    (peptide, whitelist_peptides.get(peptide, 0))
                    for peptide in sorted(peptide_whitelist)
                )
            output["whitelist_counts_output"] = str(counts_path)
    elif counts_output is not None:
        raise ValueError("counts_output requires a peptide whitelist")
    return output


@dataclass(frozen=True)
class _AnchorMetadata:
    left: str
    right: str
    insertion_nt: int
