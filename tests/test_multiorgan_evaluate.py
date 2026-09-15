from aav9_sma.models.evaluate import sequence_distance_split


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
