import numpy as np
import pandas as pd

from aav9_sma.models.evaluate import (
    CROSS_ANIMAL_TEST_ROLE,
    DEVELOPMENT_TEST_ROLE,
    benchmark_multitask_animal_holdout,
    benchmark_multitask_ensemble_leave_one_animal_out,
    sequence_distance_split,
    validate_test_role,
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
    assert all(row["test_role"] == DEVELOPMENT_TEST_ROLE for row in rows)


def test_animal4_cannot_be_labeled_final_blind() -> None:
    validate_test_role(DEVELOPMENT_TEST_ROLE, includes_animal4=True)
    try:
        validate_test_role("final_blind_external", includes_animal4=True)
    except ValueError as error:
        assert "Animal 4" in str(error)
    else:
        raise AssertionError("Animal 4 must not be labeled as a final blind test")


def test_leave_one_animal_out_audits_every_animal(tmp_path) -> None:
    rng = np.random.default_rng(11)
    alphabet = np.array(list("ACDEFGHIKLMNPQRSTVWY"))
    peptides = ["".join(row) for row in rng.choice(alphabet, size=(240, 7))]
    frame = pd.DataFrame({"AA": peptides, "rpm_whitelist__virus_prod2": 1.0})
    endpoints = ("brain", "spinal_cord")
    for endpoint_index, endpoint in enumerate(endpoints):
        base = np.array([peptide.count("A") + endpoint_index for peptide in peptides])
        for animal in (1, 2, 3, 4):
            enrichment = base + rng.normal(0, 0.1, len(frame)) + animal * 0.01
            frame[f"rpm_whitelist__{endpoint}_a{animal}"] = 2**enrichment
            frame[
                f"log2enr_whitelist__{endpoint}_a{animal}__over__virus_prod2"
            ] = enrichment
    source = tmp_path / "reconstructed.csv"
    frame.to_csv(source, index=False)

    rows = benchmark_multitask_ensemble_leave_one_animal_out(
        source,
        endpoints=endpoints,
        ensemble_size=2,
        test_fraction=0.25,
        max_iter=8,
    )

    assert len(rows) == 8
    assert {row["held_out_animal"] for row in rows} == {1, 2, 3, 4}
    assert all(row["test_role"] == CROSS_ANIMAL_TEST_ROLE for row in rows)
    assert all(row["train_rows"] > 0 and row["test_rows"] > 0 for row in rows)
