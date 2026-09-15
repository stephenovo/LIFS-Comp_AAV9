"""Bowtie2-aligned extraction matching the published Fit4Function workflow."""

from __future__ import annotations

import csv
import gzip
import re
import subprocess
import tempfile
from collections import Counter
from pathlib import Path
from typing import TextIO

from aav9_sma.data.fastq import INSERTION_NT_LENGTH, translate_21mer

REFERENCE_NAME = "fit4function_aav9_k449r_insert"
REFERENCE_PREFIX = (
    "CCAACGAAGAAGAAATTAAAACTACTAACCCGGTAGCAACGGAGTCCTATGGACAAGTGGCCACAAACCACCAGAGTGCCCAA"
)
REFERENCE_SUFFIX = "GCACAGGCGCAGACCGGTTGGGTTCAAAACCAAGGAATACTTCCG"
REFERENCE_SEQUENCE = REFERENCE_PREFIX + ("N" * INSERTION_NT_LENGTH) + REFERENCE_SUFFIX
INSERTION_REFERENCE_START = len(REFERENCE_PREFIX)
INSERTION_REFERENCE_STOP = INSERTION_REFERENCE_START + INSERTION_NT_LENGTH

BOWTIE2_ARGUMENTS = (
    "--end-to-end",
    "--very-sensitive",
    "--np",
    "0",
    "--n-ceil",
    "L,21,0.5",
    "--xeq",
    "-N",
    "1",
    "--reorder",
    "--score-min",
    "L,-0.6,-0.6",
    "-5",
    "8",
    "-3",
    "8",
)

CIGAR_TOKEN = re.compile(r"(\d+)([MIDNSHP=X])")


def write_bowtie2_reference(output_prefix: str | Path) -> dict[str, object]:
    """Write the published short reference and build a Bowtie2 index."""
    output_prefix = Path(output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fasta_path = output_prefix.with_suffix(".fasta")
    fasta_path.write_text(f">{REFERENCE_NAME}\n{REFERENCE_SEQUENCE}\n", encoding="ascii")
    completed = subprocess.run(
        ["bowtie2-build", "--quiet", str(fasta_path), str(output_prefix)],
        check=True,
        capture_output=True,
        text=True,
    )
    return {
        "reference_name": REFERENCE_NAME,
        "reference_length": len(REFERENCE_SEQUENCE),
        "insertion_reference_start_zero_based": INSERTION_REFERENCE_START,
        "insertion_reference_stop_zero_based_exclusive": INSERTION_REFERENCE_STOP,
        "fasta_path": str(fasta_path),
        "index_prefix": str(output_prefix),
        "bowtie2_build_stderr": completed.stderr.strip(),
    }


def _parse_cigar(cigar: str) -> list[tuple[int, str]]:
    tokens = [(int(length), operation) for length, operation in CIGAR_TOKEN.findall(cigar)]
    if not tokens or "".join(f"{length}{operation}" for length, operation in tokens) != cigar:
        raise ValueError(f"Unsupported CIGAR: {cigar}")
    return tokens


def extract_aligned_insertion(
    sequence: str,
    quality: str,
    reference_start_one_based: int,
    cigar: str,
    minimum_phred: int = 20,
) -> tuple[str, str | None, str | None]:
    """Extract the query bases aligned to the 21-N reference interval."""
    reference_position = reference_start_one_based - 1
    query_position = 0
    pieces: list[str] = []
    quality_pieces: list[str] = []
    indel_in_insertion = False

    for length, operation in _parse_cigar(cigar):
        if operation in {"M", "=", "X"}:
            overlap_start = max(reference_position, INSERTION_REFERENCE_START)
            overlap_stop = min(reference_position + length, INSERTION_REFERENCE_STOP)
            if overlap_start < overlap_stop:
                query_start = query_position + overlap_start - reference_position
                query_stop = query_start + overlap_stop - overlap_start
                pieces.append(sequence[query_start:query_stop])
                quality_pieces.append(quality[query_start:query_stop])
            reference_position += length
            query_position += length
        elif operation in {"D", "N"}:
            if max(reference_position, INSERTION_REFERENCE_START) < min(
                reference_position + length, INSERTION_REFERENCE_STOP
            ):
                indel_in_insertion = True
            reference_position += length
        elif operation == "I":
            if INSERTION_REFERENCE_START <= reference_position <= INSERTION_REFERENCE_STOP:
                indel_in_insertion = True
            query_position += length
        elif operation == "S":
            query_position += length
        elif operation in {"H", "P"}:
            continue

    insertion = "".join(pieces)
    insertion_quality = "".join(quality_pieces)
    if indel_in_insertion or len(insertion) != INSERTION_NT_LENGTH:
        return "indel_or_incomplete", insertion or None, None
    if "N" in insertion:
        return "undetermined", insertion, None
    if any(base not in "ACGT" for base in insertion):
        return "invalid_base", insertion, None
    if any((ord(character) - 33) < minimum_phred for character in insertion_quality):
        return "low_quality", insertion, None
    peptide = translate_21mer(insertion)
    if "*" in peptide:
        return "stop_codon", insertion, peptide
    return "valid", insertion, peptide


def _open_counts_output(path: Path) -> TextIO:
    if path.suffix in {".gz", ".gzip"}:
        return gzip.open(path, "wt", newline="", encoding="utf-8")
    return path.open("w", newline="", encoding="utf-8")


def count_bowtie2_fastq(
    fastq_path: str | Path,
    index_prefix: str | Path,
    peptide_whitelist: set[str],
    counts_output: str | Path,
    threads: int = 2,
) -> dict[str, object]:
    """Stream Bowtie2 SAM output and count Q20, in-frame, whitelist-matched 7-mers."""
    fastq_path = Path(fastq_path)
    counts_output = Path(counts_output)
    counts_output.parent.mkdir(parents=True, exist_ok=True)
    statuses: Counter[str] = Counter()
    peptides: Counter[str] = Counter()
    whitelist_peptides: Counter[str] = Counter()

    command = [
        "bowtie2",
        *BOWTIE2_ARGUMENTS,
        "--threads",
        str(threads),
        "-x",
        str(index_prefix),
        "-U",
        str(fastq_path),
    ]
    with tempfile.TemporaryFile(mode="w+", encoding="utf-8") as stderr:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=stderr,
            text=True,
            encoding="ascii",
        )
        if process.stdout is None:
            raise RuntimeError("Bowtie2 stdout was not captured")
        for line in process.stdout:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                raise ValueError("Malformed SAM alignment")
            flag = int(fields[1])
            if flag & 0x4:
                statuses["unaligned"] += 1
                continue
            status, _insertion, peptide = extract_aligned_insertion(
                sequence=fields[9],
                quality=fields[10],
                reference_start_one_based=int(fields[3]),
                cigar=fields[5],
            )
            statuses[status] += 1
            if status == "valid" and peptide is not None:
                peptides[peptide] += 1
                if peptide in peptide_whitelist:
                    whitelist_peptides[peptide] += 1
        return_code = process.wait()
        stderr.seek(0)
        bowtie2_stderr = stderr.read().strip()
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, command, stderr=bowtie2_stderr)

    with _open_counts_output(counts_output) as stream:
        writer = csv.writer(stream)
        writer.writerow(["peptide_7mer", "raw_read_count"])
        writer.writerows(
            (peptide, whitelist_peptides.get(peptide, 0))
            for peptide in sorted(peptide_whitelist)
        )
    reads_examined = sum(statuses.values())
    valid_reads = statuses["valid"]
    return {
        "path": str(fastq_path),
        "mode": "bowtie2-paper-aligned",
        "reads_examined": reads_examined,
        "status_counts": dict(statuses),
        "valid_fraction": valid_reads / reads_examined if reads_examined else 0.0,
        "unique_valid_peptide_7mers": len(peptides),
        "whitelist_size": len(peptide_whitelist),
        "whitelist_matched_reads": sum(whitelist_peptides.values()),
        "whitelist_detected_peptides": len(whitelist_peptides),
        "whitelist_detection_fraction": (
            len(whitelist_peptides) / len(peptide_whitelist) if peptide_whitelist else 0.0
        ),
        "whitelist_counts_output": str(counts_output),
        "bowtie2_arguments": list(BOWTIE2_ARGUMENTS),
        "bowtie2_stderr": bowtie2_stderr,
        "method_note": (
            "Alignment and Q20/indel filtering follow the published short-reference method. "
            "The public 100K amino-acid list is a subset of the original synthetic library."
        ),
    }
