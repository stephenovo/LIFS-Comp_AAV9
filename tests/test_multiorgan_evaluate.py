import numpy as np
import pandas as pd

from aav9_sma.models.evaluate import (
    benchmark_multitask_animal_holdout,
    sequence_distance_split,
)


def test_distance_split_removes_single_mutation_neighbors_from_training() -> None:
    peptides = ["AAAAAAA", "CAAAAAA", "CCCCCCC", "DDDDDDD"]
    found = False
    for random_state in range(1000):
        train, test = sequence_distance_split(
            peptides, test_fraction=0.5, random_state=random_state
        )
        if test[0] and not test[1]:
            assert not train[1]
            found = True
            break
    assert found


def test_multitask_benchmark_returns_one_row_per_endpoint(tmp_path) -> None:
    rng = np.random.default_rng(3)
    alphabet = np.array(list("ACDEFGHIKLMNPQRSTVWY"))
    peptides = ["".join(row) for row in rng.choice(alphabet, size=(240, 7))]
    frame = pd.DataFrame({"AA": peptides})
    endpoints = ("brain", "spinal_cord")
    for index, endpoint in enumerate(endpoints):
        signal = np.array([peptide.count("A") + index for peptide in peptides], dtype=float)
        frame[f"log2enr_whitelist__{endpoint}_animals_1_3__over__virus_prod2"] = signal
        frame[f"log2enr_whitelist__{endpoint}_a4__over__virus_prod2"] = signal + rng.normal(
            0, 0.1, len(frame)
        )
    source = tmp_path / "reconstructed.csv"
    frame.to_csv(source, index=False)

    rows = benchmark_multitask_animal_holdout(
        source,
        endpoints=endpoints,
        test_fraction=0.25,
        max_iter=8,
    )

    assert [row["task"] for row in rows] == list(endpoints)
    assert all(row["model"] == "shared_mlp_64_32" for row in rows)
    assert all(row["train_rows"] > 0 for row in rows)
