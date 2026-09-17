# Experimental validation handoff SOP

**Project:** LIFS-Comp_AAV9  
**Version:** 1.0  
**Purpose:** Test whether computationally selected AAV9 7-mer variants can be produced and whether they improve a CNS-versus-liver proxy profile relative to the matched parental capsid.

> This is a study-design and decision SOP for transfer to a qualified AAV vector core, cell laboratory, and animal facility. The receiving facility must execute cell culture, transfection, vector purification, dosing, sample collection, and waste handling under its own validated SOPs, risk assessment, biosafety approval, and animal ethics approval. This document does not authorize laboratory work and is not a substitute for facility-specific training.

## 1. What this validation can establish

The computational screen produced hypotheses, not a therapy. Each claim has a corresponding experiment:

| Computational claim | Required observation | Claim after a positive result |
|---|---|---|
| Packaging lower bound passes | Reproducible vector yield and acceptable product quality | The variant is producible in the tested process |
| Mouse spinal-cord/brain prediction is high | Vector genome, RNA, and reporter signal in CNS tissue | The variant has improved CNS delivery or expression in the tested model |
| Mouse liver prediction is low | Lower liver vector genome and expression than the parental control at matched dose | The variant reduces liver exposure in the tested model |
| High predicted CNS:liver specificity | Pre-registered CNS:liver ratio improves while CNS signal is retained | The observed distribution is more favorable under the tested conditions |
| Site-level immune annotation | A validated neutralization assay with individual sera | Only the measured neutralization phenotype; no general immune-escape claim |

Even a successful study does not establish human motor-neuron specificity, clinical efficacy, reduced human hepatotoxicity, or superiority to Zolgensma.

## 2. Required approvals and named owners

Do not start experimental work until the receiving institution records:

- biosafety review for recombinant AAV work and an approved facility risk assessment;
- animal ethics approval for any in vivo study, including humane endpoints and monitoring;
- approved procedures for human-derived cells or sera, if used;
- trained owners for vector production, analytical QC, cell assays, animal work, pathology, statistics, and data custody;
- a frozen protocol version, sample manifest, randomization seed, exclusion rules, and primary endpoint.

The project team supplies the candidate rationale and analysis plan. The qualified facility owns procedural parameters such as cell source, passage limits, media, transfection method, purification train, injection technique, anesthesia, euthanasia, and local acceptance specifications.

## 3. Test article and control panel

Use [`experimental_validation_panel.csv`](audit_data/experimental_validation_panel.csv) as the authoritative manifest.

### 3.1 Candidates

- **Stage 1:** test all 30 computational candidates for production and product quality.
- **Tier 1:** `GKEKGGE`, `KKDGQEG`, `QGSGKEG`, `VNKGKDE`, `KGDKNEG`, `SKDNKEG`, and `GAKGGEG`. These seven passed the strict conservative definition and enter the front of the queue.
- **Tier 2:** non-Tier-1 members of the strict Pareto front.
- **Tier 3:** diverse, high-scoring near-front candidates.
- If resources allow, represent each peptide with two synonymous DNA encodings. Treat these as construct replicates, not biological replicates. This helps detect nucleotide-level or barcode-related effects.

### 3.2 Required controls

1. **Matched parental control:** AAV9 `(K449R)` with no 7-mer insertion. This is the primary comparator because the Fit4Function library used that background.
2. **No-vector/mock control:** carries the assay background through every applicable stage.
3. **Wild-type AAV9:** optional external context control. Do not substitute it for the matched `(K449R)` parent.
4. **Empirical directionality controls:** one high-packaging, one low-packaging, one high-CNS/low-liver, and one high-liver/low-CNS sequence selected deterministically from the Fit4Function audit. These controls check assay direction, but they are in-sample references and are not independent validation.

All capsids compared within an experiment must carry the same vector genome cassette and be analyzed at matched genome-containing particle input. Use a reporter cassette for ranking. Introduce an SMN1 cassette only after a capsid has passed the distribution study.

## 4. Blinding, randomization, and batch design

1. A person outside sample processing assigns blind IDs. Keep the sequence-to-blind-ID key read-only until the primary analysis is locked.
2. Randomize construct order across cloning, production, assay plates, tissue processing, and analytical runs.
3. Produce candidates and the parental comparator in the same campaign when possible. Do not confound capsid with production day or operator.
4. Target three independent production batches for candidates that may advance. A second synonymous construct is useful, but it does not replace an independent production batch.
5. For animal studies, randomize animals after eligibility assessment and balance cage, sex, litter, and dosing order when applicable. Blind dosing labels, tissue processing, imaging, and primary analysis where feasible.
6. Record every deviation before unblinding.

## 5. Stage 0 — protocol lock and construct verification

### Inputs

- frozen candidate/control manifest;
- one matched reporter genome design;
- facility-approved plasmid system and production SOP;
- preregistration completed from [`WET_LAB_PREREGISTRATION_TEMPLATE.md`](WET_LAB_PREREGISTRATION_TEMPLATE.md).

### Procedure

1. Confirm the parental sequence, `(K449R)` background, and intended insertion between VP1 positions 588/589 using the facility's sequence convention.
2. Generate a construct map for every sample. Record peptide, DNA encoding, plasmid identifier, lot, and custody chain.
3. Verify the complete engineered region by sequencing before vector production. Investigate mixed traces, unintended substitutions, indels, and rearrangements.
4. Freeze the reporter cassette, promoter, enhancer, polyadenylation element, genome size, and purification approach across capsids.
5. Create blind IDs and release only blind IDs to operators.

### Gate 0

Advance only constructs with confirmed identity and an approved production record. A failed construct is rebuilt; it is not silently replaced or relabeled.

## 6. Stage 1 — production and packaging qualification

The receiving core performs production, harvest, purification, formulation, and storage under its validated AAV SOP. Addgene's HEK293 AAV production protocol is a useful public reference, but local validated procedures take precedence.

### Required measurements

For each production batch, record:

- genome-containing particle titer by a validated qPCR or ddPCR method with nuclease treatment and documented standard/reference material;
- capsid particle concentration by a validated capsid assay when available;
- genome-to-capsid ratio and a full/empty estimate by the core's validated orthogonal method when available;
- capsid and genome identity;
- visible aggregation, recovery, formulation, freeze/thaw history, and analytical run QC;
- sterility/bioburden, endotoxin, and mycoplasma status appropriate to the planned downstream use;
- yield normalized to production scale and to the matched parental batch.

### Pre-registered project gate

The following are project defaults, not universal release specifications. The receiving core may replace them before unblinding and must document the reason.

- identity is correct;
- all downstream-use QC specifications pass;
- normalized genome-containing yield is at least `0.5×` the matched parent in at least two of three independent batches;
- no major aggregation, integrity, or full/empty defect relative to the parent;
- the result is reproduced in an independent production batch.

Report exact values and confidence intervals. Do not convert a continuous yield into only pass/fail.

### Stage 1 output

Rank survivors by reproducibility and product quality. Advance Tier 1 first, then use measured production evidence to fill a practical Stage 2 panel. Keep the parental and mock controls in every Stage 2 experiment.

## 7. Stage 2 — cell-based target and liver off-target study

### Cell systems

Use at least:

- a qualified human iPSC-derived motor-neuron system or another well-characterized motor-neuron model for target relevance;
- HepG2 and/or THLE-2 for continuity with the available human liver-cell labels;
- a robust permissive cell line as an assay-performance control, not as evidence of motor-neuron tropism.

Confirm motor-neuron identity and culture quality with the facility's validated marker panel, such as ChAT, ISL1/2, and HB9/MNX1, plus viability and morphology criteria. Record donor, differentiation batch, maturity window, passage or day of differentiation, and cell density.

### Design

1. Use the same reporter genome across capsids.
2. Test a pre-registered multi-level vector-genome input range that avoids ceiling effects and excessive toxicity in each cell system. Report actual `vg/cell`; do not assume a range validated in HEK293 applies to iPSC-derived neurons.
3. Include at least three independent biological experiments. Technical wells quantify assay precision but do not increase biological sample size.
4. Balance blind IDs across plates. Include mock and parental controls on every plate.
5. Pre-register one uptake endpoint and one later expression endpoint using the facility's validated timing.

### Readouts

- intracellular vector genomes per cell or per diploid genome;
- reporter RNA normalized to validated reference genes;
- reporter protein or fluorescence with a predefined segmentation method;
- cell viability and morphology;
- in motor-neuron cultures, reporter colocalization with motor-neuron markers;
- expression normalized to vector-genome uptake, so entry and post-entry expression are not conflated.

### Stage 2 decision rule

Advance variants that:

- passed Stage 1;
- retain motor-neuron signal at least comparable to the matched parent;
- show liver-cell uptake or expression no higher than the parent, with preference for a reproducible reduction;
- improve the pre-registered motor-neuron:liver ratio;
- show no material viability penalty;
- reproduce across an independent vector batch and biological experiment.

Select four to six capsids for Stage 3 using the measured target:liver trade-off, not the original model score alone.

## 8. Stage 3 — systemic mouse biodistribution and transduction

Animal work requires an approved protocol and facility-specific procedures. The study below defines comparisons and measurements, not injection or surgical technique.

### Design

- Test four to six Stage 2 survivors, the matched `(K449R)` parent, mock, and one informative empirical comparator if product quality permits.
- Use the same reporter genome, matched administered genome-containing dose, formulation, route, and handling.
- Use the intravenous route when the question is systemic CNS-versus-liver distribution, matching the project hypothesis.
- Determine group size by an a priori power analysis on the primary endpoint. Fit4Function used groups of five for several individual-variant follow-ups; this is contextual evidence, not an automatic sample-size justification.
- Balance sex and pre-specified covariates or justify a restricted population.
- Record exclusions and humane endpoints before unblinding.

### Timepoints

Use two conceptually distinct endpoints if resources permit:

1. **Early biodistribution:** an early vector-genome endpoint aligned as closely as feasible with the Fit4Function two-hour organ label, to test the model's proxy directly.
2. **Established transduction:** a later endpoint around three weeks, aligned with the source paper's individual-variant follow-up, to measure vector genomes, RNA, reporter protein, histology, and clinical chemistry.

If only one timepoint is possible, choose it before unblinding and state which claim can no longer be tested.

### Tissues and measurements

Collect predefined brain regions, cervical/thoracic/lumbar spinal cord, liver, heart, kidney, spleen, lung, dorsal-root ganglia where approved, blood/serum, and any protocol-mandated safety tissues.

Measure:

- vector genomes per diploid genome with tissue-specific assay QC;
- transgene RNA and reporter protein;
- spinal-cord reporter colocalization with motor-neuron markers, including anterior-horn analysis where anatomically feasible;
- liver histology and serum ALT, AST, bilirubin, and other facility-approved safety chemistry;
- blinded pathology and image quantification;
- animal-level values. Do not treat multiple fields of view or qPCR wells as independent animals.

### Primary endpoint

Pre-register one primary endpoint, recommended as:

```text
log2 specificity = mean(log2 spinal-cord signal) - mean(log2 liver signal)
```

Brain may be a key secondary endpoint or combined with spinal cord only if that combination is defined before data review. Analyze vector genomes, RNA, and protein separately; do not merge them into one undocumented score.

### Stage 3 project success rule

The default rule is:

- at least a `2×` improvement in the spinal-cord:liver ratio versus the matched parent;
- spinal-cord signal is non-inferior to the parent under a pre-specified margin;
- liver signal is no higher than the parent;
- no worse predefined safety finding;
- the result is supported by more than one measurement layer and is not driven by one outlier or one production batch.

The effect-size threshold and non-inferiority margin must be finalized using pilot variance before unblinding. A ratio can improve because both organs decrease, so spinal-cord retention is mandatory.

## 9. Stage 4 — SMA relevance and immune characterization

These studies are separate from capsid ranking and should follow only after a capsid passes the reporter study.

### SMA relevance

- Package an SMN1 expression cassette in the parental and selected capsids using matched genome design and QC.
- In patient-derived SMA iPSC motor neurons, measure SMN protein, nuclear-gem restoration, cell health, and a validated disease-relevant phenotype.
- If an SMA animal model is justified, predefine molecular, histological, motor, survival, and safety endpoints. Do not infer therapeutic benefit from reporter biodistribution alone.

### Neutralizing-antibody study

- Use a validated cell-based neutralization assay with individual serum samples and matched AAV9 controls.
- Report the assay definition and neutralization titer or IC50 distribution, detection limits, sample provenance, and uncertainty.
- The current `site_context_only` annotation does not predict immune escape. Only measured data may support a sequence-specific neutralization statement.

## 10. Statistical analysis

1. Define the experimental unit before data collection: production batch for process reproducibility, cell differentiation/donor batch for cell studies, and animal for in vivo studies.
2. Use technical replicates to estimate assay precision, not to inflate `n`.
3. Base sample size on the primary endpoint, a scientifically meaningful effect, pilot or literature variance, power, alpha, and expected attrition.
4. Report individual observations, effect sizes, 95% confidence intervals, and exact sample sizes.
5. Use a model that reflects batch and repeated structure. Include production batch and cell differentiation batch as blocking or random effects when appropriate.
6. Control multiplicity for many capsids and secondary endpoints. Preserve the single pre-registered primary comparison.
7. Define below-detection handling, missing data, exclusions, transformations, and outlier policy before unblinding.
8. Keep raw, processed, and analysis-ready data separate. Never overwrite raw measurements.

Enter observations in [`wet_lab_results_template.csv`](audit_data/wet_lab_results_template.csv) and decisions in [`wet_lab_decision_log.csv`](audit_data/wet_lab_decision_log.csv). The Excel handoff workbook contains the same panel and entry schemas.

## 11. Failure interpretation

| Observation | Likely interpretation | Next action |
|---|---|---|
| Low yield in every batch | Packaging prediction failed or construct/process interaction | Stop; verify identity and repeat only with a justified process check |
| Good yield, poor motor-neuron signal | Packaging and target transduction are separable | Do not advance for the SMA objective |
| CNS:liver ratio improves because both fall | Apparent specificity without retained target delivery | Fail unless spinal-cord non-inferiority is met |
| Liver falls, heart/kidney rises | Off-target redistribution | Review full-organ profile; do not call it safer |
| Cell result succeeds, mouse result fails | In vitro model does not predict systemic distribution | Prioritize in vivo evidence and investigate mechanism |
| One batch succeeds, another fails | Manufacturing instability | Hold; investigate batch/process interaction |
| Reporter succeeds, SMN1 fails | Cargo-dependent performance | Limit the claim to reporter cargo and redesign the translational study |
| Neutralization changes without tropism benefit | Immune phenotype is separate from delivery objective | Report separately; do not rescue a failed delivery candidate |

## 12. Completion package from the receiving laboratory

The handoff is complete only when the project receives:

- final construct maps and identity records;
- unblinded sample key after analysis lock;
- production batch records and complete QC values;
- raw and processed assay data with units and detection limits;
- randomization, exclusion, and protocol-deviation logs;
- analysis code and software versions;
- animal protocol identifiers and ARRIVE Essential 10 reporting items;
- a signed stage decision for every candidate, including failed candidates.

## 13. References

- Bryant et al. [Fit4Function: data-driven AAV capsid engineering](https://www.nature.com/articles/s41467-024-50555-y), including the supplementary methods for individual-variant follow-up.
- Addgene. [AAV production in HEK293 cells](https://www.addgene.org/protocols/aav-production-hek293-cells/). Public reference only; use the receiving core's validated SOP.
- Percie du Sert et al. [ARRIVE Guidelines 2.0](https://arriveguidelines.org/arrive-guidelines).
- NIH. [Guidelines for Research Involving Recombinant or Synthetic Nucleic Acid Molecules](https://osp.od.nih.gov/wp-content/uploads/NIH_Guidelines.pdf), or the applicable local equivalent.
- FDA. [Zolgensma prescribing information](https://www.fda.gov/media/126109/download).

