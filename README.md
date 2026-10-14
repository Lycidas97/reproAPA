# reproAPA analysis

This repository contains the analysis code for the manuscript **Recurrent Gene-Level Responses beneath Divergent Alternative Polyadenylation in Human and Mouse Reprogramming**.

## Edition

This edition accompanies the manuscript production files of 13 October 2026. It covers Figures 1–4, Supplementary Figures S1–S21 and Tables S1–S10 (58 data sheets, including S4j, and 202 source blocks).

The [current edition record](provenance/current-edition.json) binds the final figure and workbook hashes to their production scripts. The v11f data deposit uses 81 small archives and six companion files under reserved DOI `10.5281/zenodo.23347745`; at this edition's preparation it is awaiting publication. See the [execution guide](docs/reproduction.md) for restoration and version details.

The analysis uses polyadenylation site (PAS) counts from polyApipe, the integrated v2 PAS reference and the frozen expression-state framework. Scientific input matrices and source-data archives are distributed separately.

## Contents

| Folder | Content |
| --- | --- |
| [`modules/`](modules/01_atlas_membership/README.md) | Ten analysis entry guides and the upstream expression workflows. |
| [`code/`](code/README.md) | Captured analysis, figure and workbook sources, retaining their original project paths. |
| `docs/` | The [analysis map](docs/current-analysis.md), [execution guide](docs/reproduction.md) and [software environments](docs/scientific-environments.md). |
| `provenance/` | The [source inventory](provenance/source_manifest.tsv), [current manuscript record](provenance/current-edition.json), data manifest and verification records. |

## Main figures

| Figure | Analysis |
| --- | --- |
| 1 | Expression states and route-level PAS-position shifts |
| 2 | Directional recurrence, proximal–distal composition and strict tandem 3′ UTR changes |
| 3 | Cross-species direction pairing, threshold sensitivity, event classes and PDLIM7 |
| 4 | Interval inclusion, whole-gene shifts and predicted-site retention |

The [current output map](docs/current-manuscript-map.tsv) lists the manuscript figures and tables, their final production scripts and their hashes. Some final figures apply PDF edits after their original compositors; follow the ordered production stages in the execution guide.

## Validation

```sh
python3 scripts/validate_repository.py
```

Validation checks the source hashes, Python syntax, current v11f output/workbook identities and the retained historical environment records without loading scientific matrices. The [current verification record](provenance/current-verification.json) reports the inventory checks, rejection cases and isolated S4j rebuild. The [11 October verification record](provenance/verification.json) reports tests of the unchanged numerical sources; the [current edition](provenance/current-edition.json) records the later production updates. Packaging this update did not rerun scientific analyses. [Third-party sources](docs/third-party.md) gives the license status.

The original [11 October edition record](provenance/release.json), [output map](docs/manuscript-map.tsv) and source-data selections are historical. Its public commit `ffa0cd20288699795702c5eebad77d1741298dab` remains available. The v11f data package contains that original repository snapshot under `code/` and the final updated production sources under `project/`; this GitHub edition combines the unchanged analysis sources with those final production updates. No data archive was changed for this GitHub refresh.

The authors' code is released under the [MIT License](LICENSE). Third-party components retain their own terms.
