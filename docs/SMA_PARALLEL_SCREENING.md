# SMA-prioritized parallel screening strategies

## Purpose

The organ-prediction model continues to predict brain, spinal cord, liver,
heart, and kidney together. Two independent decision strategies are then run
on the same predictions:

1. **Joint multi-organ strategy** - configurable brain/spinal target score,
   liver/off-target penalties, Pareto status, weight stability, and diversity.
2. **Sequential strategy** - packaging/distance eligibility followed by spinal
   cord, brain, liver, heart, and kidney filters in that exact order.

Keeping the prediction layer shared makes the comparison interpretable: any
difference in selected sequences comes from selection logic rather than from
different fitted models or different virtual pools.

## Joint strategy change

The historical unweighted field is retained for reporting:

```text
f_cns = 0.5 * brain + 0.5 * spinal_cord
```

Selection now uses an explicit target score. The default preserves the
previously evaluated 50/50 baseline:

```text
f_sma_target = 0.50 * brain + 0.50 * spinal_cord
```

An explicit SMA sensitivity run can use:

```text
f_sma_target = 0.30 * brain + 0.70 * spinal_cord
```

`display_score`, specificity, Pareto analysis, weight sensitivity, and the
balanced shortlist use `f_sma_target`. The target-oriented shortlist group is
`spinal_favoring`, ordered first by spinal prediction and then by brain.
Non-equal weights are a screening hypothesis and are never enabled silently.

Whole-brain and whole-spinal-cord DNA enrichment remain organ-level mouse
proxies. They are not evidence of motor-neuron-specific functional
transduction.

## Sequential strategy

The currently configured comparison stages are:

| Order | Organ | Direction | Retain among current survivors |
|---:|---|---|---:|
| 1 | spinal cord | maximize | 25% |
| 2 | brain | maximize | 50% |
| 3 | liver | minimize | 50% |
| 4 | heart | minimize | 75% |
| 5 | kidney | minimize | 75% |

These fractions have not been promoted as a superior selector. They are
deliberately CLI-configurable. Each cutoff is recomputed
among candidates that survived the preceding stage. Ties at a cutoff are kept,
so the realized fraction can differ slightly from the requested fraction.

The sequential route is lexicographic: after all gates, spinal cord remains
the first sorting key, followed by brain, liver, heart, and kidney. A greedy
Hamming-distance step then produces a diverse final panel.

## Run both strategies on the same predictions

First create the joint results with `screen-virtual`. To run the optional
30/70 sensitivity, pass `--brain-target-weight 0.30` and
`--spinal-target-weight 0.70`. Then run the sequential route and comparison
without retraining:

```bash
aav9-sma compare-screening-strategies \
  artifacts/virtual_screen_ranked.csv.gz \
  artifacts/virtual_screen_shortlist.csv \
  --spinal-retain-fraction 0.25 \
  --brain-retain-fraction 0.50 \
  --liver-retain-fraction 0.50 \
  --heart-retain-fraction 0.75 \
  --kidney-retain-fraction 0.75 \
  --candidate-count 30 \
  --output-annotated artifacts/sequential_ranked.csv.gz \
  --output-sequential-shortlist artifacts/sequential_shortlist.csv \
  --output-comparison artifacts/strategy_comparison.csv \
  --output-summary artifacts/strategy_comparison_summary.json
```

This separation is intentional: selection stringency can be tuned repeatedly
without rerunning the expensive organ ensemble.

## Quantitative comparison

The comparison summary reports:

- exact shared, joint-only, and sequential-only sequence counts;
- Jaccard similarity and overlap coefficient;
- rank correlation among shared sequences, when at least two are shared;
- median predicted spinal, brain, liver, heart, kidney, SMA-target score,
  display score, and model disagreement for each shortlist;
- minimum within-shortlist pairwise Hamming distance;
- Jensen-Shannon divergence between amino-acid frequency distributions.

These quantities answer different questions:

| Measure | Interpretation |
|---|---|
| Overlap/Jaccard | Do the strategies choose the same sequences? |
| Organ medians | Which prediction pressure distinguishes each strategy? |
| Shared-sequence rank correlation | Do they agree on priority within their overlap? |
| Jensen-Shannon divergence | Do they select different amino-acid compositions? |
| Hamming diversity | Is one panel more sequence-concentrated? |

Future experimental results should add the decisive comparison: packaging,
motor-neuron functional transduction, liver-cell transduction, and in-vivo
spinal/liver outcomes by strategy. Until then, the comparison describes model
behavior rather than biological superiority.

## Promotion status

The sequential route is **comparison only**. On the tracked one-million-candidate
screen, its default funnel shared 8 of 30 sequences with the 30/70 joint list,
while its median predicted spinal-cord score was 0.36 lower and its median
predicted liver score was 0.86 higher. Those results do not support replacing
the joint shortlist. Likewise, changing the joint score from 50/50 to 30/70
produced a rank Spearman correlation of 0.99985 among packaging-eligible
candidates and changed only one of 30 shortlisted sequences; this is useful as
a sensitivity check, not evidence of a biological model improvement.
