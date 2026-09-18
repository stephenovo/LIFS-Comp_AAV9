# Animal 4 evaluation policy

## Current status

Animal 4 is a **development held-out test**. It was excluded from fitting and
protected by a distance-2 sequence split, but its metrics were inspected while
comparing Ridge, Random Forest, LightGBM, the shared MLP, the ensemble, and
screening rules. That means Animal 4 is useful for model debugging and
comparison, but it is not a final blind test.

The code writes `test_role=development_holdout_animal4` into every Animal 4
benchmark row. The `validate_test_role` guard rejects any attempt to label a
result that includes Animal 4 as `final_blind_external`.

This is a disclosure of evaluation history, not a cosmetic rename. Looking at a
test result can influence model class, hyperparameters, score weights,
packaging thresholds, shortlist rules, or the claims made from them. Re-running
the same split cannot undo that information flow.

## What qualifies as a final blind test

A final blind test must use an independent animal, batch, or NHP dataset that
was held back until all of the following were frozen:

- model family and encoder;
- hyperparameters and random seeds;
- missing-label policy and training data version;
- packaging threshold and score weights;
- shortlist and diversity rules;
- output claims and acceptance criteria.

The external file should be checksummed in a data manifest and evaluated once
with the frozen pipeline. Its result should be reported separately from the
Animal 4 development table.

## Required reporting language

Use:

> Animal 4 was a distance-2, development held-out evaluation set and was
> inspected during model selection. It is not a final blind validation.

Do not use “Animal 4 blind test”, “unseen final test”, or “validated in animals”
without an independent experimental dataset and the corresponding wet-lab
evidence.
