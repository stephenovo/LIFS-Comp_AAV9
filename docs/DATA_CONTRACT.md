# Canonical data contract

Raw Fit4Function column names must be mapped to this internal schema only after their definitions
are verified against the paper and supplementary metadata.

## Identity and sequence columns

| Column | Type | Meaning |
|---|---|---|
| `variant_id` | string | Stable source or project identifier |
| `peptide_7mer` | string | Seven standard amino acids, uppercase |
| `source_split` | string, optional | Source-defined train/test/library role |
| `experiment_batch` | string, optional | Experimental batch identifier |

## Primary label columns

| Column | Type | Intended use |
|---|---|---|
| `f_pack` | float | Packaging/production target |
| `brain_mouse` | float | Mouse brain proxy target |
| `spinal_cord_mouse` | float | Mouse spinal-cord proxy target |
| `liver_mouse` | float | Mouse liver negative target |
| `heart_mouse` | float | Mouse heart off-target |
| `kidney_mouse` | float | Mouse kidney off-target |

## Optional annotation columns

| Column | Type | Intended use |
|---|---|---|
| `liver_human` | float | Human liver-cell annotation/secondary penalty |
| `lung_mouse` | float | Secondary off-target record |
| `spleen_mouse` | float | Secondary off-target record |
| `immune_footprint_risk` | float | Annotation only |
| `species` | string | Species or experimental system |

## Required audit questions

- What physical quantity does each label measure: abundance, DNA, RNA, or functional expression?
- Are labels already normalized, and if so, against what reference?
- Can scores from different organs be compared directly?
- Which labels occur on the same variants?
- Are biological replicates available?
- Does the source define a held-out or non-human-primate set?
- Which residue-numbering convention defines the insertion site?
- What license or attribution requirements apply?

Do not rename an unclear source column into a canonical biological claim. Record unresolved fields
in the audit report until their meanings are confirmed.

