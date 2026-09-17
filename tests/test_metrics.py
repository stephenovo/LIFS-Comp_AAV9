import numpy as np

from aav9_sma.models.metrics import bootstrap_confidence_intervals, regression_metrics


def test_regression_metrics_and_bootstrap_are_deterministic() -> None:
    targets = np.linspace(-2, 2, 100)
    predictions = targets + np.sin(targets) * 0.1

    point = regression_metrics(targets, predictions)
    first = bootstrap_confidence_intervals(targets, predictions, n_resamples=50, random_state=7)
    second = bootstrap_confidence_intervals(targets, predictions, n_resamples=50, random_state=7)

    assert point["pearson_r"] > 0.99
    assert first == second
    assert first["pearson_r_ci_low"] <= first["pearson_r"] <= first["pearson_r_ci_high"]
    assert first["bootstrap_resamples"] == 50
