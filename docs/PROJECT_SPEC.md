# Project specification

## Research question

Can a constrained virtual-screening workflow identify AAV9 7-mer insertion candidates that
retain predicted packaging fitness while shifting the available mouse-organ proxy profile toward
brain/spinal cord and away from liver and other off-target organs?

## Primary application context

The motivating application is future SMN1 delivery for spinal muscular atrophy (SMA). The project
engineers the delivery capsid; it does not redesign SMN1, dosage, route of administration, or a
clinical product.

## Primary endpoints

1. Packaging fitness (`F_pack`) — hard eligibility gate.
2. Mouse CNS enrichment proxy (`F_CNS`) — weighted brain and spinal-cord predictions.
3. Mouse liver burden proxy (`F_liv`) — primary negative endpoint.
4. Other off-target proxy (`F_off`) — initially heart and kidney.

Human liver-cell predictions and immune-footprint proximity are annotations unless the data audit
provides evidence for a validated role in the primary objective.

## Candidate selection

Eligible candidates first pass a conservative packaging threshold. Remaining candidates are
reported using both a presentation score and a Pareto analysis. The final set should contain
24–30 diverse candidates across CNS-favoring, liver-minimizing, and balanced groups.

## Claims we will not make

- Mouse brain/spinal-cord enrichment is not human motor-neuron specificity.
- Reduced predicted liver enrichment is not demonstrated reduction of hepatotoxicity.
- A computational score is not evidence of clinical efficacy or superiority to Zolgensma.
- No candidate is experimentally validated until packaging, transduction, biodistribution, and
  safety experiments have been performed.

## Decision gates

1. **Data gate:** required Fit4Function labels and sequence metadata are present and interpretable.
2. **Baseline gate:** held-out predictions outperform naive baselines without sequence leakage.
3. **Screening gate:** candidates pass packaging, distance, uncertainty, and diversity checks.
4. **Communication gate:** every figure and claim is consistent with the available evidence.

