from pathlib import Path

import pandas as pd

from aav9_sma.data.ena import download_fastq


def test_download_fastq_accepts_valid_cached_file(tmp_path: Path) -> None:
    path = tmp_path / "SRR1.fastq.gz"
    path.write_bytes(b"test")
    row = {
        "run_accession": "SRR1",
        "fastq_url": "https://example.invalid/SRR1.fastq.gz",
        "fastq_md5": "098f6bcd4621d373cade4e832627b4f6",
        "fastq_bytes": 4,
    }

    result = download_fastq(row, tmp_path)

    assert result["status"] == "cached"
    assert result["bytes"] == 4


def test_resolved_manifest_columns_are_machine_readable() -> None:
    frame = pd.DataFrame(
        [{"run_accession": "SRR1", "fastq_bytes": 4, "fastq_md5": "abc"}]
    )
    assert frame["fastq_bytes"].dtype.kind in "iu"
