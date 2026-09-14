# Four-person collaboration workflow

## Suggested ownership

| Track | Primary responsibility | Typical deliverables |
|---|---|---|
| A — Data and evidence | Source data, metadata, label definitions, literature | Data card, mapping table, audit report |
| B — Packaging | `F_pack`, thresholds, production-related validation | Baseline models, gate analysis |
| C — Tropism | `F_CNS`, `F_liv`, `F_off`, evaluation | Multi-task models, metrics, uncertainty |
| D — Screening and delivery | Generation, Pareto, diversity, reporting | Candidate table, figures, final report |

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

