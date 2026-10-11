# Analysis environments

## Version policy

- The environments keep the known execution versions and the original declared pins.
- If a version is not known, the reconstruction uses the latest compatible release identified in the 9 October 2026 metadata snapshot.
- The records show reconstructed versions separately from execution evidence.

The [environment index](../environments/reconstruction.json) gives the installation recipes. For each reconstruction lock, it gives all dependency pins and download hashes. The [registry](../environments/registry.json) maps each analysis component to its installation specification.

## Installation

The reconstructed environments are for Linux x86_64.

To install a reconstructed environment:

1. Find your analysis in the table headed Analysis below.
2. Create the conda environment from its `environment.yml` file.
3. Install the locked Python packages with hash checks.

This example installs the `plan-a` environment:

```sh
conda env create -f environments/reconstruction/plan-a/environment.yml
conda run -n reproapa-plan-a python -m pip install --require-hashes \
  -r environments/reconstruction/plan-a/requirements.lock
```

The environments use these Python versions:

| Environment | Python | Source of the version |
| --- | --- | --- |
| Other reconstructed environments | 3.12.14 | Selected reconstruction target |
| Common PAS and cross-null | 3.10.14 | Kept from the original records |
| Plan A | 3.12.12 | Kept from the original records |
| Legacy PAS | 3.10.12 | Kept from the original declaration |
| RNA QC | 3.10.20 | Original environment and full Linux conda lock |

To install the original RNA QC environment, run this command:

```sh
conda create -n reproapa-rna-qc \
  --file environments/locks/rna-qc/environment-linux-64.lock
```

| Analysis | Installation environment | Version choices |
| --- | --- | --- |
| Atlas membership | [`atlas`](../environments/reconstruction/atlas/environment.yml) | Python 3.12.14. Standard library only. |
| RNA QC | [Original release](../environments/locks/rna-qc/environment.yml) | Keeps the original conda and pip records. |
| Local states | [`local-states`](../environments/reconstruction/local-states/requirements.lock) | pyclustree 0.4.1, Scanpy 1.12.4. |
| Human SEACells | [`human-seacells`](../environments/reconstruction/human-seacells/requirements.lock) | Keeps the recorded SEACells 0.3.3, NumPy 2.5.2 and anndata 0.13.3.post0. |
| Mouse SEACells | [`mouse-seacells`](../environments/reconstruction/mouse-seacells/requirements.lock) | SEACells 0.3.3. The execution version is not known; 0.3.3 was the latest release in the metadata snapshot. |
| Moran metacell scVI | [`moran-scvi`](../environments/reconstruction/moran-scvi/requirements.lock) | Keeps scvi-tools 1.3.3 from the checkpoint. Adds the missing dependencies. |
| Species and temporal integration | [`species-temporal`](../environments/reconstruction/species-temporal/requirements.lock) | Keeps the original declared pins for Scanpy, scVI and moscot. |
| genes2genes (optional) | [`genes2genes`](../environments/reconstruction/genes2genes/requirements.lock) | A separate environment, because genes2genes needs NumPy <2. |
| Common PAS statistics and cross-null | [`common-pas`](../environments/reconstruction/common-pas/requirements.lock) | Keeps NumPy 1.26.4, pandas 2.3.3 and SciPy 1.15.3. Uses anndata 0.11.4, selected for compatibility with Python 3.10.14 and NumPy 1.26.4. |
| Final Plan A, bridges and cycle | [`plan-a`](../environments/reconstruction/plan-a/requirements.lock) | NumPy 2.5.3, pandas 3.0.6, SciPy 1.18.1. Adds the missing dependencies for the cycle and bridge analyses. |
| Common PAS annotation | [Current APADETECT snapshot](../environments/snapshots/apadetect-0.4.1/requirements.txt) | Python 3.10.12. You must get the original APADETECT 0.4.1 wheel separately. |
| Legacy PAS comparisons | [`legacy-pas`](../environments/reconstruction/legacy-pas/requirements.lock) | Keeps the original core pins and the captured 0.3.0 comparison source. [GSE221739](../environments/reconstruction/legacy-gse221739/requirements.lock) has a separate environment that keeps pysam 0.21.0. |
| Cis preparation | [`cis`](../environments/reconstruction/cis/requirements.lock) | Keeps the original declared pins. Adds the missing sequence dependencies. |
| Figures and workbooks | `plan-a` and the [S9 Node runtime](../environments/evidence/current-workbook-node.json) | Node 24.19.0, @oai/artifact-tool 2.8.89 (supplied externally), Arial font. |

For manuscript rendering and assembly, use the [`revision-1011` environment](../environments/reconstruction/revision-1011/environment.yml) and its [`requirements.lock`](../environments/reconstruction/revision-1011/requirements.lock). This reconstructed environment uses Python 3.12.14 and has 37 locked dependencies for Linux x86_64. It keeps all 33 Plan A dependency pins and adds pyBigWig 0.3.22, [pikepdf 10.16.0](https://pypi.org/project/pikepdf/10.16.0/) and [xlrd 2.0.2](https://pypi.org/project/xlrd/2.0.2/). The [additional release metadata](../environments/evidence/revision-1011-pypi-releases.json) records the new downloadable versions.

For the UTR/ortholog, naive-comparison, direct-specificity and RNA-module analyses, the [runtime records](../environments/evidence/revision-1011-runtimes.json) report Python 3.12.14, NumPy 2.3.5 and pandas 2.2.3. The separate [`revision-1011-statistics` lock](../environments/reconstruction/revision-1011-statistics/requirements.lock) keeps these versions and has 39 locked dependencies. Both reconstruction recipes leave the original locks unchanged. Neither environment was fully installed or used to rerun the analyses during this update.

The original mouse requirements file declares two packages that conflict:

- genes2genes 0.2.0 needs NumPy <2.
- decoupler 1.9.2 needs NumPy >=2.

Thus, you cannot install the two packages in one environment. The original declaration does not change. genes2genes has a separate reconstructed environment. The index records all other version changes.

## Current runtime records

The [server inventories](../environments/evidence/current-server-environments.json) record five Python installations on the current server. They record the package versions, operating system and architecture for each installation, plus conda builds and hashes where available. They also record FIMO, bedtools and the GPUs. `environments/snapshots/` contains the related pip and conda snapshots.

The [CUDA metadata](../environments/evidence/current-cuda.json) records Torch 2.8.0 in `/root/miniconda3/envs/scanpy/bin/python` as a CPU-only build, with no CUDA or cuDNN available in that Python environment. The server has GPUs.

The executable hash of FIMO 5.5.2 agrees with the pinned deployment resource. The server has bedtools 2.31.1 and samtools 1.13.

The inspected `scanpy` environment has scvi-tools 1.4.0; the original Moran checkpoint records 1.3.3. The model lock keeps 1.3.3. The current inventories give more information. They are not records of the historical runs.

Limits of this update:

- The dependency resolution and the hashes were examined.
- The reconstructed environments were not fully installed.
- The reconstructed environments were not used to fit the models again.

To fit a model again on a GPU, install a PyTorch CUDA build that is compatible with your driver. Refer to the [upstream installation instructions](https://pytorch.org/get-started/locally/).

## Frozen scientific references

The [Moran model record](../environments/evidence/moran-scvi-model-identity.json) keeps these items of the fit:

- 11,104 metacells and 4,000 HVGs
- a 30-dimensional latent space
- the seed and the training configuration
- the hashes of the model and the inputs

The record does not contain the optimizer, the learning rate, the stopping patience or the original GPU. These items are not package versions. Thus, the reconstructed environment does not supply them.

The common structural caches use GENCODE v44 (human) and GENCODE M33 (mouse), both from Ensembl 110. The cache signatures confirm these versions.

Other resources use different versions:

- GENCODE v40 and M25 are for the integrated external PAS support.
- Ensembl 116 is for the other transcript and ortholog resources.

The [resource versions](../environments/resource-versions.json) keep the MirGeneDB-curated inputs, CisBP-RNA build 2.00 and the UniProtKB entry identities.
