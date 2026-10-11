# Captured sources

This folder contains the captured sources for this edition. The files keep their original paths relative to the project folder. The scripts find their dependencies through these paths. Do not move the files.

| Folder | Content |
| --- | --- |
| `notes/2026-10-10-*/` | Analyses of direction, composition, cross-species responses, marker references, naive states and potential mechanisms. |
| `notes/2026-10-10-divergence-main-figures/` | Scripts for Figures 1–4. |
| `notes/2026-10-10-supplement-v11/` | Supplementary Figures S1–S21, legends and the S1–S10 workbook. |
| `analysis_code/pas_rerun/` | The maintained PAS pipeline for input preparation and statistics, with its regression tests. |
| `reviews/2026-10-08-plan-a-redesign/` | The final protocols and code for the conditional-reference, cycle and cross-species analyses. |
| `reviews/2026-10-08-plan-a-analysis/` | The initial Plan A code. The redesign uses its extraction and bridge scripts to make its inputs. This folder also contains the profile-entry comparison that the manuscript uses. |
| `reviews/2026-10-08-narrative-logic-chain/netshift-local/` | The local classification of interval and whole-gene shifts. |
| `notes/2026-10-08-narrative-figures/` | Earlier figure compositors that prepare panels from source PDFs and recorded values. |
| `github/` and `external/` | Captured numerical dependencies. |
| `vendor/apadetect_annotation/` | The legacy APADETECT 0.3.0 annotation source. Only the comparison with the legacy annotation uses it. |

Some protocols and analysis notes in `reviews/` and `notes/` are in Chinese.

For analyses A1, A2, A3 and A8, the Plan A redesign gives the final definitions. The initial Plan A definitions of these analyses do not apply.

The current common annotation uses APADETECT 0.4.1. This package is external. It is not the 0.3.0 copy in `vendor/`.

For more information, see the [analysis map](../docs/current-analysis.md), the [execution guide](../docs/reproduction.md) and the [analysis environments](../docs/scientific-environments.md).
