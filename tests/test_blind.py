from pathlib import Path

import pytest

from aav9_sma.blind import build_blind_freeze_manifest


def _config(tmp_path: Path) -> dict[str, object]:
    files = []
    for role in ("analysis_plan", "candidate_panel", "prediction_table", "evaluation_code"):
        path = tmp_path / f"{role}.txt"
        path.write_text(role + "\n", encoding="utf-8")
        files.append({"role": role, "path": path.name})
    return {
        "study_id": "BLIND-001",
        "final_blind_dataset_id": "independent_mouse_batch_2027_01",
        "outcome_data_accessed": False,
        "roles": {
            "model_team": "team-model",
            "custodian": "team-key",
            "wet_lab": "team-lab",
            "analysis": "team-stats",
        },
        "locked_decisions": {
            "model_and_encoder": "shared MLP ensemble and one-hot encoder",
            "training_data_version": "frozen manifest",
            "primary_endpoint": "spinal Spearman correlation",
            "analysis_method": "one-shot locked script",
            "success_rule": "point and lower confidence limits",
            "exclusion_rule": "QC failures defined before outcomes",
            "unblinding_condition": "QC lock and signed report hash",
        },
        "files": files,
    }


def test_blind_freeze_hashes_required_files(tmp_path: Path) -> None:
    payload = build_blind_freeze_manifest(
        _config(tmp_path),
        tmp_path,
        repository_commit="a" * 40,
        repository_clean=True,
    )

    assert payload["freeze_status"] == "frozen_pre_outcome"
    assert payload["test_role"] == "final_blind_external"
    assert len(payload["files"]) == 4
    assert all(len(entry["sha256"]) == 64 for entry in payload["files"])


def test_blind_freeze_rejects_animal4_or_prior_outcome_access(tmp_path: Path) -> None:
    animal4 = _config(tmp_path)
    animal4["final_blind_dataset_id"] = "Animal 4 rerun"
    with pytest.raises(ValueError, match="Animal 4"):
        build_blind_freeze_manifest(
            animal4,
            tmp_path,
            repository_commit="a" * 40,
            repository_clean=True,
        )

    seen = _config(tmp_path)
    seen["outcome_data_accessed"] = True
    with pytest.raises(ValueError, match="outcome_data_accessed"):
        build_blind_freeze_manifest(
            seen,
            tmp_path,
            repository_commit="a" * 40,
            repository_clean=True,
        )


def test_blind_freeze_requires_distinct_roles(tmp_path: Path) -> None:
    config = _config(tmp_path)
    config["roles"]["analysis"] = "team-model"
    with pytest.raises(ValueError, match="distinct"):
        build_blind_freeze_manifest(
            config,
            tmp_path,
            repository_commit="a" * 40,
            repository_clean=True,
        )
