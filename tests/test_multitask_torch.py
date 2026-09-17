import numpy as np
import pytest

pytest.importorskip("torch")

from aav9_sma.models.multitask_torch import fit_masked_multitask_mlp


def test_masked_multitask_model_uses_partially_observed_rows() -> None:
    rng = np.random.default_rng(4)
    features = rng.normal(size=(80, 12)).astype(np.float32)
    targets = np.column_stack([features[:, 0] + features[:, 1], features[:, 2] - features[:, 3]])
    targets[::3, 0] = np.nan
    targets[1::4, 1] = np.nan
    train_rows = np.zeros(80, dtype=bool)
    validation_rows = np.zeros(80, dtype=bool)
    train_rows[:55] = True
    validation_rows[55:70] = True

    result = fit_masked_multitask_mlp(
        features,
        targets,
        train_rows,
        validation_rows,
        hidden_sizes=(8, 4),
        max_epochs=4,
        batch_size=16,
        patience=2,
        random_state=11,
    )

    assert result.predictions.shape == targets.shape
    assert np.isfinite(result.predictions).all()
    assert result.training_rows == 64
    assert result.task_observations[0] < result.training_rows
    assert result.task_observations[1] < result.training_rows
