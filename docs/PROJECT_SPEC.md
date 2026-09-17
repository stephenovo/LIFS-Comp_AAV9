# Project specification

## Research question

Can a constrained virtual-screening workflow identify AAV9 7-mer insertion candidates that
retain predicted packaging fitness while shifting the available mouse-organ proxy profile toward
brain/spinal cord and away from liver and other off-target organs?

## Primary application context

The motivating application is future SMN1 delivery for spinal muscular atrophy (SMA). The project
engineers the delivery capsid; it does not redesign SMN1, dosage, route of administration, or a
clinical product.

The Fit4Function labels were measured in an AAV9 (K449R) backbone with a 7-mer inserted between
VP1 residues 588 and 589. The project does not assume that learned effects transfer unchanged to a
different capsid backbone.

## Primary endpoints

1. Packaging fitness (`F_pack`) — hard eligibility gate.
2. Mouse CNS enrichment proxy (`F_CNS`) — weighted brain and spinal-cord predictions.
3. Mouse liver burden proxy (`F_liv`) — primary negative endpoint.
4. Other off-target proxy (`F_off`) — initially heart and kidney.

Brain and spinal-cord heads are trained separately and combined only after endpoint-specific
validation. Their labels are early (two-hour) mouse-organ vector-genome biodistribution proxies,
not motor-neuron transduction measurements.

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

1. **Data gate:** SRA-derived organ labels reproduce the public liver labels and expected replicate
   structure; the processed multi-organ workbooks alone do not pass this gate.
2. **Baseline gate:** held-out predictions outperform naive baselines without sequence leakage.
3. **Screening gate:** candidates pass packaging, distance, uncertainty, and diversity checks.
4. **Communication gate:** every figure and claim is consistent with the available evidence.

Current checkpoint: the data gate passed through public-Liver reconstruction (`r = 0.9778`), and
the baseline gate passed for all five organ heads under a distance-separated sequence holdout with
animal 4 held out from training. A shared multi-task MLP then improved Animal 4 Pearson `r` over the best
single-task baseline by `0.061–0.082` for all five endpoints; a five-model ensemble improved it by
`0.071–0.109`. A masked-loss PyTorch model and a LightGBM/physicochemical challenger were then
tested under the same split and were not promoted because neither improved the five-endpoint mean.
A one-million-sequence screen produced 6,016 packaging-eligible variants, 169 strict Pareto
candidates, eight strict conservative candidates, and a 30-row diverse computational shortlist
containing seven of the conservative candidates. The computational screening gate is complete;
experimental validation remains open.
