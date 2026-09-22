"""Guards for freezing a genuinely external, one-shot blind evaluation."""

from __future__ import annotations

import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path

from aav9_sma.repro import sha256_file, verify_manifest

BLIND_FREEZE_STATUS = "frozen_pre_outcome"
FINAL_BLIND_TEST_ROLE = "final_blind_external"
REQUIRED_FILE_ROLES = {
    "analysis_plan",
    "candidate_panel",
    "prediction_table",
    "evaluation_code",
}
REQUIRED_DECISIONS = {
    "model_and_encoder",
    "training_data_version",
    "primary_endpoint",
    "analysis_method",
    "success_rule",
    "exclusion_rule",
    "unblinding_condition",
}


def _git_state(root: Path, *, include_untracked: bool = True) -> tuple[str, bool, list[str]]:
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    lines = [
        line
        for line in subprocess.run(
            [
                "git",
                "status",
                "--porcelain",
                f"--untracked-files={'all' if include_untracked else 'no'}",
            ],
            cwd=root,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.splitlines()
        if line
    ]
    return commit, not lines, lines


def build_blind_freeze_manifest(
    config: dict[str, object],
    root: str | Path,
    *,
    repository_commit: str,
    repository_clean: bool,
) -> dict[str, object]:
    """Validate a preregistration config and hash every frozen analysis input."""
    root_path = Path(root).resolve()
    study_id = str(config.get("study_id", "")).strip()
    dataset_id = str(config.get("final_blind_dataset_id", "")).strip()
    if not study_id or not dataset_id:
        raise ValueError("study_id and final_blind_dataset_id are required")
    normalized_dataset = "".join(
        character for character in dataset_id.lower() if character.isalnum()
    )
    if "animal4" in normalized_dataset:
        raise ValueError("Animal 4 cannot be registered as a final blind dataset")
    if config.get("outcome_data_accessed") is not False:
        raise ValueError("outcome_data_accessed must be explicitly false at freeze time")
    if not repository_clean:
        raise ValueError("Repository must be clean before a blind-study freeze")
    if len(repository_commit) != 40 or any(
        character not in "0123456789abcdefABCDEF" for character in repository_commit
    ):
        raise ValueError("repository_commit must be a full 40-character hexadecimal Git commit")

    roles = config.get("roles")
    if not isinstance(roles, dict):
        raise ValueError("roles must define model_team, custodian, wet_lab, and analysis")
    required_roles = ("model_team", "custodian", "wet_lab", "analysis")
    role_values = [str(roles.get(role, "")).strip() for role in required_roles]
    if any(not value for value in role_values):
        raise ValueError("roles must define model_team, custodian, wet_lab, and analysis")
    if len(set(role_values)) != len(role_values):
        raise ValueError("Blind-study roles must be assigned to distinct people or teams")

    decisions = config.get("locked_decisions")
    if not isinstance(decisions, dict):
        raise ValueError("locked_decisions must be an object")
    missing_decisions = sorted(
        decision for decision in REQUIRED_DECISIONS if not str(decisions.get(decision, "")).strip()
    )
    if missing_decisions:
        raise ValueError(f"Missing locked blind-study decisions: {missing_decisions}")

    entries = config.get("files")
    if not isinstance(entries, list):
        raise ValueError("files must be a list")
    file_roles = {str(entry.get("role", "")) for entry in entries if isinstance(entry, dict)}
    missing_file_roles = sorted(REQUIRED_FILE_ROLES - file_roles)
    if missing_file_roles:
        raise ValueError(f"Missing frozen file roles: {missing_file_roles}")
    frozen_files: list[dict[str, str]] = []
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValueError("Every files entry must be an object")
        role = str(entry.get("role", "")).strip()
        raw_path = str(entry.get("path", "")).strip()
        if not role or not raw_path:
            raise ValueError("Every files entry requires role and path")
        if "outcome" in role.lower() or "unblind" in role.lower():
            raise ValueError("Frozen inputs must not include outcomes or an unblinding key")
        relative = Path(raw_path)
        candidate = (root_path / relative).resolve()
        try:
            candidate.relative_to(root_path)
        except ValueError as error:
            raise ValueError(f"Frozen file escapes root: {relative}") from error
        if not candidate.is_file():
            raise FileNotFoundError(candidate)
        frozen_files.append(
            {
                "role": role,
                "path": relative.as_posix(),
                "sha256": sha256_file(candidate),
            }
        )

    return {
        "schema_version": "1.0",
        "freeze_status": BLIND_FREEZE_STATUS,
        "test_role": FINAL_BLIND_TEST_ROLE,
        "study_id": study_id,
        "final_blind_dataset_id": dataset_id,
        "frozen_at_utc": datetime.now(UTC).isoformat(),
        "repository_commit": repository_commit,
        "outcome_data_accessed_before_freeze": False,
        "roles": {role: value for role, value in zip(required_roles, role_values, strict=True)},
        "locked_decisions": decisions,
        "files": frozen_files,
        "policy": {
            "interim_outcome_review": "forbidden",
            "model_or_threshold_changes_after_freeze": "forbidden",
            "unblinding_before_locked_report": "forbidden",
            "failure_consequence": (
                "The cohort becomes development data and a new independent cohort is required "
                "for any final-blind claim."
            ),
        },
    }


def freeze_blind_study(
    config_path: str | Path,
    root: str | Path,
    output_path: str | Path,
) -> dict[str, object]:
    """Create a freeze manifest only from a clean Git checkout."""
    root_path = Path(root).resolve()
    commit, clean, dirty_entries = _git_state(root_path)
    if not clean:
        preview = dirty_entries[:10]
        raise ValueError(f"Repository must be clean before freeze; found: {preview}")
    config = json.loads(Path(config_path).read_text(encoding="utf-8"))
    payload = build_blind_freeze_manifest(
        config,
        root_path,
        repository_commit=commit,
        repository_clean=clean,
    )
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return payload


def verify_blind_freeze(manifest_path: str | Path, root: str | Path) -> dict[str, object]:
    """Verify hashes and confirm that the checked-out code still matches the freeze."""
    manifest = Path(manifest_path)
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    if payload.get("freeze_status") != BLIND_FREEZE_STATUS:
        raise ValueError("Manifest is not a frozen pre-outcome blind-study record")
    if payload.get("test_role") != FINAL_BLIND_TEST_ROLE:
        raise ValueError("Manifest is not labeled final_blind_external")
    commit, clean, dirty_entries = _git_state(Path(root).resolve(), include_untracked=False)
    file_result = verify_manifest(manifest, root)
    return {
        **file_result,
        "expected_repository_commit": payload.get("repository_commit"),
        "actual_repository_commit": commit,
        "repository_commit_matches": commit == payload.get("repository_commit"),
        "tracked_repository_clean": clean,
        "tracked_repository_changes": dirty_entries,
        "ok": bool(
            file_result["ok"] and commit == payload.get("repository_commit") and clean
        ),
    }
