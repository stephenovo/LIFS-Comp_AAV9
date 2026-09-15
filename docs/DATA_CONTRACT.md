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
| `animal_id` | string, optional | Mouse identifier; required for organ-label reconstruction |
| `technical_replicate` | string, optional | FASTQ technical-replicate identifier |
| `sra_run` | string, optional | NCBI SRA run accession for provenance |

## Primary label columns

| Column | Type | Intended use |
|---|---|---|
| `f_pack` | float | Packaging/production target |
| `brain_mouse` | float | Mouse brain proxy target |
| `spinal_cord_mouse` | float | Mouse spinal-cord proxy target |
| `liver_mouse` | float | Mouse liver negative target |
| `heart_mouse` | float | Mouse heart off-target |
| `kidney_mouse` | float | Mouse kidney off-target |

For Fit4Function, the mouse-organ labels refer to two-hour post-injection vector-genome DNA
biodistribution, not cell-type-resolved transduction. Processed multi-organ workbooks do not contain
amino-acid sequences or a public ID-to-sequence lookup. These labels must therefore be rebuilt from
BioProject `PRJNA1131359` before being mapped into the canonical columns.

## Raw-count reconstruction columns

| Column | Type | Meaning |
|---|---|---|
| `raw_read_count` | integer | Valid reads assigned to a 7-mer in one SRA run |
| `rpm` | float | Reads per million within one technical replicate |
| `replicate_mean_rpm` | float | Mean RPM across technical replicates for one biological sample |
| `virus_reference_rpm` | float | Production-virus abundance used as enrichment denominator |
| `log2_enrichment` | float | `log2(replicate_mean_rpm / virus_reference_rpm)` after documented zero handling |
| `detected` | boolean | Whether the variant passed the analysis detection rule |

The virus-reference rule is fixed to production round 2 after reconstructed liver values reproduced
the public sequence-linked `Liver` column (`Pearson r = 0.9778`). Production round 3 is retained for
denominator sensitivity analysis, not used in the final five-organ table.

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
