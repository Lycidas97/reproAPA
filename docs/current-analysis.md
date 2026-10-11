# Analysis map

## Edition

This edition accompanies the manuscript version prepared for submission on 11 October 2026.

The manuscript contains four main figures, 21 supplementary figures and Tables S1–S10. The workbook has 57 data sheets, one Contents sheet and one Sources sheet (200 source blocks). These counts come from the workbook XML and its build record.

The [manuscript record](../provenance/release.json) identifies the version used for this edition. The [output map](manuscript-map.tsv) lists its figures, tables and source scripts. The [source-data index](../provenance/source-data-release.json) identifies the input archive.

## Analysis chain

1. **Expression framework.** Modules 01–04 retain the atlas, RNA QC, local states, metacells and integration workflows. The framework contains seven studies, 161 atlas libraries, 551,201 cells and 11,104 metacells. Module 04 includes the Moran metacell-scVI configuration, the 4,000-HVG list and the SEACells sources.
2. **PAS statistics.** `code/analysis_code/pas_rerun/` retains the input preparation, common annotation adapters and downstream statistics. Deployment entry points are `run_sources.py` and `prepare_pipeline.py`.
3. **Recurrence and composition.** The new `code/notes/2026-10-10-directional-recurrence/`, `recurrence-utr-ortholog/`, `recurrence-component-split/` and `route-step-decomposition/` sources analyze response direction, strict tandem compatibility and route-level composition.
4. **Cross-species responses.** `code/notes/2026-10-10-ortholog-threshold-scan/`, `utr-event-stratification/` and `conserved-response/` test direction pairing, event classes and corresponding PAS switches. The marker-reference and naive-state workflows provide the sensitivity analyses.
5. **Potential mechanisms.** The proliferation, CPA-factor, UGUA-signature, cycle-source-attribution and shared-proximal-mechanism workflows retain their protocols and input identities. RNA-factor expression and sequence signatures do not measure CPA activity. No analysis was preregistered.
6. **Local intervals.** `code/reviews/2026-10-08-narrative-logic-chain/netshift-local/compute.py` and `code/notes/2026-10-08-narrative-figures/opposite-sign/derive_opposite.py` retain the interval versus whole-gene comparison used by Figure 4.
7. **Figures and tables.** `code/notes/2026-10-10-divergence-main-figures/` makes Figures 1–4. `code/notes/2026-10-10-supplement-v11/` assembles the supplementary set, legends and S1–S10 workbook, using the S9f base workbook. The final workbook includes S5d.6–S5d.8. The scripts under `code/notes/2026-10-11-caption-figure-fixes/` apply the final main-figure label, panel-order and colour edits; `code/notes/2026-10-11-supplement-caption-revision/render.py` rebuilds S1 and S2 with native state IDs.

## Source dependencies

The supporting source PDFs and tables are recorded in the [input archive inventory](../provenance/source-data-dependencies.json).

The Plan A redesign protocols define the conditional-reference, cycle and cross-species analyses. The original extraction and bridge scripts provide their inputs. Source files preserve the project paths that their dependencies use.

The primary recurrence sets remain 483 human and 604 mouse genes, including 268 and 420 longer responses, respectively.
