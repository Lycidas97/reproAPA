# reproAPA analysis

This repository contains the analysis code for the manuscript **Recurrent Gene-Level Responses beneath Divergent Alternative Polyadenylation in Human and Mouse Reprogramming**.

## Edition

This edition accompanies the manuscript version prepared for submission on 11 October 2026. It covers Figures 1–4, Supplementary Figures S1–S21 and Tables S1–S10 (57 data sheets).

The analysis uses polyadenylation site (PAS) counts from polyApipe, the integrated v2 PAS reference and the frozen expression-state framework. Scientific input matrices and source-data archives are distributed separately.

## Contents

| Folder | Content |
| --- | --- |
| [`modules/`](modules/01_atlas_membership/README.md) | Ten analysis entry guides and the upstream expression workflows. |
| [`code/`](code/README.md) | Captured analysis, figure and workbook sources, retaining their original project paths. |
| `docs/` | The [analysis map](docs/current-analysis.md), [execution guide](docs/reproduction.md) and [software environments](docs/scientific-environments.md). |
| `provenance/` | The [source inventory](provenance/source_manifest.tsv), [manuscript record](provenance/release.json) and verification records. |

## Main figures

| Figure | Analysis |
| --- | --- |
| 1 | Expression states and route-level PAS-position shifts |
| 2 | Directional recurrence, proximal–distal composition and strict tandem 3′ UTR changes |
| 3 | Cross-species direction pairing, threshold sensitivity, event classes and PDLIM7 |
| 4 | Interval inclusion, whole-gene shifts and predicted-site retention |

The [output map](docs/manuscript-map.tsv) lists the manuscript figures and tables, their source scripts and their hashes.

## Validation

```sh
python3 scripts/validate_repository.py
```

Validation checks the source hashes, Python syntax and output records without loading the scientific input matrices. [Verification records](provenance/verification.json) and [third-party sources](docs/third-party.md) give the checks and license status.

The authors' code is released under the [MIT License](LICENSE). Third-party components retain their own terms.
