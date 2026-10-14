# Inputs and execution

Full scientific reproduction requires the inputs and runtimes below. Repository validation does not rerun these analyses.

## Scientific inputs

To run the scientific analyses, you must have these inputs:

- the prepared PAS matrices
- the frozen GEX and metacell inputs
- the reference resources

Each input must agree with its original hash in the deployment configuration.

The current v11f deposit is prepared under reserved DOI `10.5281/zenodo.23347745`. It was not yet public when this edition was made. The [package record](../provenance/v11f-data-package.json) binds the 81 archives and six companion files by size and SHA-256. The [restored-file inventory](../provenance/v11f-FILES.tsv) binds the recovered files.

After publication, download all `reproAPA-v11f-part*.tar.gz` files and the companion files from the record. Extract the independent archives into one folder, then restore and verify:

```sh
for archive in reproAPA-v11f-part*.tar.gz; do tar -xzf "$archive"; done
cd reproAPA-v11f
python3 restore.py
python3 verify.py
```

Use the restored `project/` tree as the working project. It contains the final figures, workbook, source tables, input PDFs and final production scripts. The `code/` tree inside the deposit is the original public repository commit `ffa0cd20288699795702c5eebad77d1741298dab`. This updated GitHub edition additionally captures the final production scripts from `project/`. Copy its `code/notes/`, `code/reviews/` and `code/analysis_code/` into the corresponding project paths if using this repository for execution. The archived data package and its original code snapshot have not been rewritten.

The [historical input archive inventory](../provenance/source-data-dependencies.json) describes the earlier companion packages, which v11f already merges under `project/`:

- **Supporting source archive.** It keeps the project-relative paths. It contains supporting S1–S9 workbook sources and figure assets.
- **Generator inputs and current implementations.** `Supplementary_Source_Code_v11b.zip` supplies the frozen PDFs and legend sources under `revision-backups/20261011-before-v11-apply/`. Restore these files at their recorded paths before running the compositors.
- **S9 companion.** It contains the S9 sources that the manuscript uses, and the table filters that the companion declares.

The raw sequencing data and the complete upstream training matrices are not in these archives.

## Figures and tables for submission

The [current output map](current-manuscript-map.tsv) identifies all 25 final figures and the workbook. It is bound to the v11f restored-file inventory. The [11 October source-data index](../provenance/source-data-release.json) is a historical record, with the earlier 57-sheet workbook; it does not describe the current S4j edition.

Use the `revision-1011` [locked environment](../environments/reconstruction/revision-1011/requirements.lock) for the rendering and assembly sources. Use the separate [`revision-1011-statistics` lock](../environments/reconstruction/revision-1011-statistics/requirements.lock) for the statistical workflows with recorded NumPy 2.3.5 and pandas 2.2.3. Figure scripts also require Poppler (`pdftoppm`), `qpdf`, the recorded fonts and frozen production PDFs.

In the restored project layout, run the relevant `analyze.py` or extraction script under `notes/2026-10-10-*`. The main compositors are `compose.py` followed by `compose_v11.py` in `notes/2026-10-10-divergence-main-figures/`. For supplementary figures, run the individual scripts under `notes/2026-10-10-supplement-v11/figures/`, then `assemble.py`. The final S1/S2 renderer is `notes/2026-10-11-supplement-caption-revision/render.py`. The workbook entry point is `notes/2026-10-10-supplement-v11/tables/assemble_workbook.py`. The `verify_tables.py` script checks it against its frozen source tables and base workbook.

For the final main figures, first rebuild the v11 base figures in an isolated project layout. Then run `fix_fig1.py` through `fix_fig4.py` under `notes/2026-10-11-caption-figure-fixes/`; they read the unedited base PDFs from `Fig/` and write the final PDFs under `out/`. The source-data archive contains the final PDFs and edit receipts.

The 12–13 October production stages follow those original builders. Restore the pre-edit PDFs from the corresponding `revision-backups/` locations into an isolated project layout before each edit stage. Do not apply the edits again to the final PDFs in `project/Fig/`. The order is the 12 October main and supplementary cleanup, the 13 October supplementary text fixes and figure-density fixes, and the final S19 tightening. The final S19 stage uses `notes/2026-10-13-s19-layout/plan.py` to select its pre-density input. The current output map names each final producer, and archived receipts identify its inputs.

For the workbook, the original `assemble_workbook.py` produces the 57-sheet base. The final insertion uses `notes/2026-10-13-table-s4j/build_s4j.py`. In an isolated project, put the restored `notes/2026-10-10-supplement-v11/out/Supplementary_Tables_S1-S10.v11.xlsx` at `table/Supplementary_Tables_S1-S10.xlsx`, retain the two frozen source tables at their recorded paths, set `REPRO_ROOT` to that project and `S4J_OUT` to a new output directory, and run the insertion. `verify_s4j.py` uses the same environment variables and checks values and preservation against that base. It does not recalculate statistics. A rebuild may differ in ZIP timestamps; compare workbook part contents and values as well as file hashes.

These scripts use project-relative paths and write their recorded output directories. Full scientific recomputation still needs the frozen matrices and external resources described in their protocols. This code refresh checks inventories, syntax and frozen file identities; it does not rerun scientific analyses or regenerate the published figures.

## Software environments

[Analysis environments](scientific-environments.md) gives the installation recipes and the dependency locks with hashes.

[`verification.json`](../provenance/verification.json) records the runtime of the local tests separately.

For the interface runtime only, use `code/analysis_code/pas_rerun/requirements-interface.txt`. This file has fewer packages.

## PAS statistics

The current common annotation must have APADETECT 0.4.1. The lock identities of its code and cache are in `environments/evidence/pas-runtime-and-structural-cache.json`. The 0.3.0 source in `code/vendor/` is only for the legacy numerical comparison tests. It is not the installed 0.4.1 package.

To run the PAS statistics:

1. Use the captured configuration `code/notes/2026-10-05-apadetect-analysis-intake/pipeline.environment.json` with the contracts in the same folder.
2. Change `project_root`, `frozen_gex`, `python` and the resource paths to the paths on your system.
3. Do not change the pinned hashes of the scientific inputs.
4. Replace the `source_overrides` value with the content of [`config/source-overrides.json`](../config/source-overrides.json).
5. Add the repository root path to the start of each path from `source-overrides.json`.
6. Run the pipeline with a new output directory:

```sh
python3 code/analysis_code/pas_rerun/run_sources.py \
  --environment /path/to/deployment/pipeline.environment.json \
  --prepared-native /path/to/prepared-native \
  --out /path/to/new-statistics-run --workers 2
```

The output directory must not exist. If it exists, the script stops.

| Option | Function |
| --- | --- |
| `--prepare-only` | Examines and records the inputs. It does not calculate statistics. |
| `--inputs` | Adds independent prepared sources. These sources do not automatically become manuscript sources. |

## Final Plan A and local interval analyses

To run these analyses:

1. Follow each frozen protocol in `code/reviews/2026-10-08-plan-a-redesign/`.
2. Put the frozen input tables back at their recorded project paths.

The initial Plan A extraction scripts and the A4 and A6 bridge scripts make the inputs for these final definitions.

The final cycle results have zero accepted contexts and 26 planned statistical test entries with NA results. Thus, a successful run of the code does not mean that the balance passed or that the effects can be estimated.

## Supporting figure and workbook sources

These instructions cover the supporting figure assets, earlier compositors and S1–S9 workbook. The figure scripts use the original project paths. There is no portable command that makes all figures from the matrices.

To make the figures:

1. Put `code/analysis_code/`, `code/notes/` and `code/reviews/` back into the original project folder structure.
2. Supply the `main.tex` root marker and the `Fig/` folder.
3. Supply the pinned master PDFs, the fonts and the registered bundles of source tables.

The main standalone builders accept `--records` and `--out`. The supporting supplementary compositor applies its source patches in memory. The narrative compositors use the source PDFs to prepare panels. Use the file identities in the input archive inventory.

The supporting S1–S9 workbook has two parts:

- **Tables S1–S8.** pandas and xlsxwriter assemble these tables.
- **Table S9.** The original Python scripts and the `@oai/artifact-tool` Node scripts prepare the table and add it to the workbook. Then they make sure that the result is correct. These scripts also need the original payload and the reference workbook.

To assemble the workbook again, you must supply all the original inputs and runtimes.

## Inventory maintenance

To copy the pinned sources from the project workspace into this repository, run this command:

```sh
python3 scripts/sync_current_analysis.py --workspace /path/to/reproAPA
```

Before the script copies files, it compares each pinned source and the manuscript with the recorded hashes. If a source or the manuscript changed, the script stops. You must then update the edition records manually. The script does not find a new edition automatically.

The following packager is retained for the historical 11 October selection. It does not package current v11f:

```sh
python3 scripts/package_source_data.py --workspace /path/to/reproAPA --out /path/to/new-release.zip
```

The script puts the historical pinned source-data files into the archive and records the identities of the earlier companion archives. It stops if those inputs changed or the output exists. Use the current v11f record and its restoration instructions for the current edition.

