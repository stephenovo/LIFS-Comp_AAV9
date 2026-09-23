# SMA-prioritized parallel screening strategies

## Purpose

The organ-prediction model continues to predict brain, spinal cord, liver,
heart, and kidney together. Two independent decision strategies are then run
on the same predictions:

1. **Joint multi-organ strategy** - spinal-prioritized SMA target score,
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

Selection now uses an explicit SMA target score:

```text
f_sma_target = 0.30 * brain + 0.70 * spinal_cord
```

`display_score`, specificity, Pareto analysis, weight sensitivity, and the
balanced shortlist use `f_sma_target`. The target-oriented shortlist group is
now `spinal_favoring`, ordered first by spinal prediction and then by brain.

Whole-brain and whole-spinal-cord DNA enrichment remain organ-level mouse
proxies. They are not evidence of motor-neuron-specific functional
transduction.

## Sequential strategy

Default stages are:

| Order | Organ | Direction | Retain among current survivors |
|---:|---|---|---:|
| 1 | spinal cord | maximize | 25% |
| 2 | brain | maximize | 50% |
| 3 | liver | minimize | 50% |
| 4 | heart | minimize | 75% |
| 5 | kidney | minimize | 75% |

The fractions are deliberately CLI-configurable. Each cutoff is recomputed
among candidates that survived the preceding stage. Ties at a cutoff are kept,
so the realized fraction can differ slightly from the requested fraction.

The sequential route is lexicographic: after all gates, spinal cord remains
the first sorting key, followed by brain, liver, heart, and kidney. A greedy
Hamming-distance step then produces a diverse final panel.

## Run both strategies on the same predictions

First create the spinal-prioritized joint results with `screen-virtual` as
usual. Then run the sequential route and comparison without retraining:

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
