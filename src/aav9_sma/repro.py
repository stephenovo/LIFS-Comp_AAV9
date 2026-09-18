"""Utilities for checking the integrity of external reproducibility inputs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


def sha256_file(path: str | Path, *, chunk_size: int = 1024 * 1024) -> str:
    """Return the SHA256 digest of a file without loading it all into memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_manifest(manifest_path: str | Path, root: str | Path = ".") -> dict[str, object]:
    """Verify files listed by a JSON manifest.

    The manifest must contain ``{"files": [{"path": ..., "sha256": ...}]}``.
    Paths are resolved relative to ``root`` and are never downloaded or modified.
    """
    manifest = Path(manifest_path)
    root_path = Path(root).resolve()
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    entries = payload.get("files") if isinstance(payload, dict) else None
    if not isinstance(entries, list):
        raise ValueError("Manifest must contain a top-level 'files' list")

    verified: list[str] = []
    missing: list[str] = []
    mismatched: list[dict[str, str]] = []
    invalid: list[str] = []
    for entry in entries:
        if not isinstance(entry, dict) or not isinstance(entry.get("path"), str):
            invalid.append(str(entry))
            continue
        relative = Path(entry["path"])
        expected = entry.get("sha256")
        label = relative.as_posix()
        if not isinstance(expected, str) or len(expected) != 64:
            invalid.append(label)
            continue
        candidate = (root_path / relative).resolve()
        try:
            candidate.relative_to(root_path)
        except ValueError:
            invalid.append(label)
            continue
        if not candidate.is_file():
            missing.append(label)
            continue
        actual = sha256_file(candidate)
        if actual != expected.lower():
            mismatched.append({"path": label, "expected": expected.lower(), "actual": actual})
            continue
        verified.append(label)

    return {
        "manifest": str(manifest),
        "root": str(root_path),
        "files_listed": len(entries),
        "verified": verified,
        "missing": missing,
        "mismatched": mismatched,
        "invalid": invalid,
        "ok": not (missing or mismatched or invalid) and len(verified) == len(entries),
    }
