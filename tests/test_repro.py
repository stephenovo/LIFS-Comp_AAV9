import hashlib
import json
from pathlib import Path

import pytest

from aav9_sma.repro import verify_manifest


def test_verify_manifest_reports_checksum_and_missing_file(tmp_path: Path) -> None:
    good = tmp_path / "good.txt"
    good.write_text("reproducible\n", encoding="utf-8")
    digest = hashlib.sha256(good.read_bytes()).hexdigest()
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "files": [
                    {"path": "good.txt", "sha256": digest},
                    {"path": "missing.dat", "sha256": "0" * 64},
                ]
            }
        ),
        encoding="utf-8",
    )

    result = verify_manifest(manifest, tmp_path)

    assert result["verified"] == ["good.txt"]
    assert result["missing"] == ["missing.dat"]
    assert result["ok"] is False


def test_verify_manifest_requires_files_list(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.json"
    manifest.write_text("{}", encoding="utf-8")
    with pytest.raises(ValueError, match="files"):
        verify_manifest(manifest, tmp_path)
