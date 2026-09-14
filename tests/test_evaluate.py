from pathlib import Path

import pandas as pd

from aav9_sma.models.evaluate import benchmark_screen_models


def test_benchmark_screen_models(tmp_path: Path) -> None:
    frame = pd.DataFrame(
        {
            "AA": ["AAAAAAA", "CAAAAAA", "DAAAAAA", "EAAAAAA", "FAAAAAA"],
            "Liver": [0.0, 1.0, 2.0, 3.0, 4.0],
        }
    )
    path = tmp_path / "screen.csv"
    frame.to_csv(path, index=False)

    rows = benchmark_screen_models(
        path,
        model_names=("ridge",),
        tasks=("Liver",),
        test_fraction=0.4,
    )

    assert len(rows) == 1
    assert rows[0]["train_rows"] == 3
    assert rows[0]["test_rows"] == 2
