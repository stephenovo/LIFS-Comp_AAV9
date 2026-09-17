# Four-person collaboration workflow

## Suggested ownership

| Track | Primary responsibility | Typical deliverables |
|---|---|---|
| A — Data and evidence | Source data, metadata, label definitions, literature | Data card, mapping table, audit report |
| B — Packaging | `F_pack`, thresholds, production-related validation | Baseline models, gate analysis |
| C — Tropism | `F_CNS`, `F_liv`, `F_off`, evaluation | Multi-task models, metrics, uncertainty |
| D — Screening and delivery | Generation, Pareto, diversity, reporting | Candidate table, figures, final report |

## Experimental-validation ownership

The wet-lab handoff is designed so that a qualified external laboratory can
execute it without changing the scientific question. Before contacting a core,
the team should name one owner for each row below. One person may hold more than
one role, but every role must have a named backup.

| Role | Owns | Must not change alone |
|---|---|---|
| Study lead | Scope, approvals, facility coordination | Primary endpoint or claim boundary |
| Construct lead | Sample manifest, maps, identity records | Sequence-to-blind-ID key after lock |
| Assay lead | Production and analytical QC transfer | Acceptance rules after unblinding |
| Biology lead | Cell/animal model suitability and tissue plan | Target population or timepoint after data review |
| Data lead | Randomization, data receipt, code, audit trail | Exclusion status without documented reason |
| Statistics lead | Power, model, multiplicity, confidence intervals | Primary analysis after outcome review |

Before the next group discussion, decide only the items that affect budget and
scientific interpretation:

1. whether Stage 1 starts with all 30 candidates or a seven-candidate Tier 1 pilot;
2. which qualified AAV core or collaborator can own production and QC;
3. which motor-neuron model is realistically available;
4. whether the intended endpoint is early Fit4Function-comparable biodistribution,
   later transduction, or both;
5. the maximum number of candidates affordable in Stage 2 and Stage 3;
6. who holds the blind key and who signs stage decisions.

Do not tune thresholds or remove failed candidates during the meeting. Record
changes in the preregistration and decision log before outcome data are reviewed.

## Git workflow

1. Create a short branch from `main`, for example `data/audit-fit4function`.
2. Keep raw data outside Git unless redistribution is explicitly allowed.
3. Put reusable logic in `src/aav9_sma`; notebooks are for exploration and communication.
4. Add or update tests for reusable code.
5. Open a pull request with the question, method, result, and limitations.
6. Require at least one teammate review before merging into `main`.

## Definition of done

A task is done only when:

- code runs from a clean environment;
- inputs and outputs are documented;
- tests pass;
- random seeds and data versions are recorded;
- biological claims are no stronger than the labels support;
- generated figures and tables can be recreated by a documented command.
