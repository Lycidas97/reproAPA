#!/usr/bin/env python3
"""Build frozen six-study QC manifests from reviewed local source registries."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
from pathlib import Path


MANIFEST_COLUMNS = [
    "study_id", "species", "dataset_id", "library_id", "sample_id",
    "sample_accession", "sample_title", "experimental_time", "condition",
    "author_state", "metadata_derivation", "source_url", "input_adapter",
    "input_path", "expected_called_cells", "basic_qc_evidence_adapter",
    "prior_basic_qc_path", "diagnostic_evidence_adapter",
    "prior_diagnostic_path",
]

ALIGNMENT_ROOT = Path(
    "/root/nfsdata/apadetect/data/reprogramming_differentiation_apa/analyses/"
    "raw_realignment_apa_2026-06-07/cellranger"
)
GSE118258_ROOT = Path(
    "/root/nfsdata/apadetect/data/reprogramming_differentiation_apa/processed/"
    "gse118258_human_pluripotency_reprogramming/cellranger_grch38_2024a_v1"
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_table(path: Path, rows: list[dict[str, object]], columns: list[str], delimiter: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter=delimiter, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sample_sources(project: Path, directory: str, series: str) -> list[dict[str, str]]:
    rows = read_csv(project / "study-state-mapping" / directory / "sample_manifest.csv")
    selected = [
        row for row in rows
        if row["series_accession"] == series and row["include_in_3prime_scrna"] == "yes"
    ]
    if not selected:
        raise ValueError(f"No included samples found for {series}")
    return selected


def snapshot_counts(project: Path) -> dict[str, int]:
    path = project / "data/multica-doc27-final-qc/cellranger_qc_snapshot.tsv"
    result: dict[str, int] = {}
    for row in read_tsv(path):
        shape = row["matrix_shape_features_barcodes_entries"].split("x")
        if len(shape) != 3:
            raise ValueError(f"Malformed matrix shape for {row['library_id']}")
        result[row["library_id"]] = int(shape[1])
    return result


def doc523_counts(project: Path) -> dict[str, int]:
    root = project / "dispatch/2026-08-21-complete-study-qc-consensus-rerun/common"
    result: dict[str, int] = {}
    for name in ("DOC523_existing_output_validation.tsv", "DOC523_final_validation.tsv"):
        for row in read_tsv(root / name):
            if row["status"] != "pass":
                raise ValueError(f"DOC-523 validation did not pass: {row['library_id']}")
            result[row["library_id"]] = int(row["called_cells"])
    return result


def gse103221_metadata(title: str) -> tuple[str, str, str]:
    chemical = re.fullmatch(r"10x_cip_d([0-9]+)", title)
    if chemical:
        return f"D{chemical.group(1)}", "chemical_CIP", "GEO title encodes CIP route and collection day"
    if title == "10x_osk_mef":
        return "MEF", "OSK_source", "GEO title identifies the uninduced MEF source"
    if title == "10x_osk_esc":
        return "ESC", "OSK_reference", "GEO title identifies the ESC reference"
    osk = re.fullmatch(r"10x_osk_d([0-9]+)", title)
    if osk:
        return f"D{osk.group(1)}", "OSK_reprogramming", "GEO title encodes OSK route and collection day"
    raise ValueError(f"Unrecognised GSE103221 title: {title}")


def wot_metadata(title: str) -> tuple[str, str, str]:
    endpoint = re.fullmatch(r"DiPSC_(2i|serum)_C[12]", title)
    if endpoint:
        return "DiPSC", endpoint.group(1), "GEO title identifies the DiPSC endpoint and culture condition"
    temporal = re.fullmatch(r"(D[0-9]+(?:\.[0-9]+)?)_(Dox|2i|serum)_C[12]", title)
    if temporal:
        return temporal.group(1), temporal.group(2), "GEO title encodes collection day, culture condition, and well"
    raise ValueError(f"Unrecognised GSE115943 title: {title}")


def gse118258_metadata(title: str) -> tuple[str, str, str]:
    mapping = {
        "D0_10X": ("D0", "OSKM_source"),
        "D2_10X": ("D2", "OSKM_reprogramming"),
        "D8_10X": ("D8", "OSKM_reprogramming"),
        "D12_10X": ("D12", "OSKM_reprogramming"),
        "D16-_10X": ("D16", "TRA-1-60_negative"),
        "D16+_10X": ("D16", "TRA-1-60_positive"),
    }
    if title not in mapping:
        raise ValueError(f"Unrecognised GSE118258 title: {title}")
    time, condition = mapping[title]
    return time, condition, "GEO title and characteristics encode collection day and TRA-1-60 selection"


def gse150311_metadata(title: str) -> tuple[str, str, str]:
    mapping = {
        "SC_RNA-seq_D0_D7": ("D0_D7_pool", "FM_reprogramming_intermediates"),
        "SC_RNA-seq_Naive": ("endpoint", "naive"),
        "SC_RNA-seq_Primed": ("endpoint", "primed"),
        "SC_RNA-seq_RSet": ("endpoint", "RSeT"),
    }
    if title not in mapping:
        raise ValueError(f"Unrecognised GSE150311 title: {title}")
    time, condition = mapping[title]
    return time, condition, "GEO title identifies pooled reprogramming intermediates or endpoint culture state"


def gse171127_metadata(accession: str) -> tuple[str, str, str]:
    mapping = {
        "GSM5219399": ("D12", "pooled_OSKM_GETM"),
        "GSM5219400": ("D6", "pooled_OSKM_GETM"),
        "GSM6071749": ("D3", "GETM"),
        "GSM6071750": ("D3", "OSKM"),
    }
    if accession not in mapping:
        raise ValueError(f"Unrecognised GSE171127 sample: {accession}")
    time, condition = mapping[accession]
    return time, condition, "GEO characteristics encode treatment day and reprogramming factors"


def raw_path(dataset_id: str, library_id: str) -> str:
    return str(ALIGNMENT_ROOT / dataset_id / library_id / "outs/filtered_feature_bc_matrix.h5")


def build_rows(
    source_rows: list[dict[str, str]], study_id: str, species: str,
    dataset_id: str, counts: dict[str, int], metadata,
    library_name=None, path_builder=None,
) -> list[dict[str, object]]:
    output = []
    for source in source_rows:
        accession = source["sample_accession"]
        library_id = library_name(accession) if library_name else f"{dataset_id}_{accession}"
        if library_id not in counts:
            raise ValueError(f"Missing called-cell count for {library_id}")
        args = (accession,) if study_id == "GSE171127" else (source["sample_title"],)
        time, condition, derivation = metadata(*args)
        input_path = path_builder(accession, library_id) if path_builder else raw_path(dataset_id, library_id)
        output.append({
            "study_id": study_id,
            "species": species,
            "dataset_id": dataset_id,
            "library_id": library_id,
            "sample_id": accession,
            "sample_accession": accession,
            "sample_title": source["sample_title"],
            "experimental_time": time,
            "condition": condition,
            "author_state": "NOT_USED",
            "metadata_derivation": derivation,
            "source_url": source["source_url"],
            "input_adapter": "cellranger_filtered_h5",
            "input_path": input_path,
            "expected_called_cells": counts[library_id],
            "basic_qc_evidence_adapter": "policy_document_recompute",
            "prior_basic_qc_path": "contract_inputs/historical_qc/report.md",
            "diagnostic_evidence_adapter": "scrublet_policy_recompute",
            "prior_diagnostic_path": "contract_inputs/historical_qc/production_policy_recommendation.tsv",
        })
    return output


def write_config(root: Path, name: str, study_id: str, species: str, label: str, rows: list[dict[str, object]], pattern: str, output_dir: str) -> str:
    config_path = root / "configs" / f"{name}.yaml"
    report_sha = sha256_file(root / "contract_inputs/historical_qc/report.md")
    policy_sha = sha256_file(root / "contract_inputs/historical_qc/production_policy_recommendation.tsv")
    text = f"""workflow_version: 3.2.0
python: /root/miniconda3/envs/single-cell-qc/bin/python
study_id: {study_id}
species: {species}
report_label: {label}
library_manifest: manifests/{name}.tsv
metadata_registry: metadata_sources/{name}.csv
metadata_contract: metadata_contracts/{name}.tsv
output_root: runs
output_dir: {output_dir}
output_mode: fresh
expected_library_count: {len(rows)}
expected_called_cell_count: {sum(int(row['expected_called_cells']) for row in rows)}
library_id_pattern: '{pattern}'
global_cell_id_template: '{{study_id}}|{{library_id}}|{{barcode}}'
adapters:
  raw_input: [cellranger_filtered_h5, cellranger_filtered_mtx]
  basic_qc_evidence: [validated_prior_run, policy_document_recompute]
  diagnostic_evidence: [validated_prior_run, scrublet_policy_recompute]
basic_qc:
  strict: {{hard_min_umi: 1500, hard_min_genes: 750, mad_n: 3, mito_cap_pct: 10, complexity_min: 0.75}}
  candidate: {{hard_min_umi: 1500, hard_min_genes: 750, mad_n: 3, mito_cap_pct: 20, complexity_min: 0.75}}
  retention_sensitivity: {{hard_min_umi: 500, hard_min_genes: 300, mad_n: 5, mito_cap_pct: 30, complexity_min: 0.75}}
sample_gates:
  primary_min_candidate_basic_qc: 100
  sensitivity_min_candidate_basic_qc: 50
scrublet:
  expected_doublet_rate: 0.06
  sim_doublet_ratio: 2
  synthetic_doublet_umi_subsampling: 1.0
  min_counts: 2
  min_cells: 3
  min_gene_variability_pctl: 85
  n_prin_comps: 30
  primary_seed: 0
  diagnostic_seeds: [1, 2]
  action: diagnostic_only
evidence_sha256:
  report.md: {report_sha}
  production_policy_recommendation.tsv: {policy_sha}
evidence_required_text:
  - 'ReproAPA RNA QC v1'
  - 'Scrublet'
"""
    config_path.write_text(text, encoding="utf-8")
    return str(config_path.relative_to(root))


def write_study(root: Path, name: str, rows: list[dict[str, object]], source_rows: list[dict[str, str]]) -> None:
    write_table(root / "manifests" / f"{name}.tsv", rows, MANIFEST_COLUMNS, "\t")
    registry = [{
        "series_accession": row["study_id"],
        "sample_accession": row["sample_accession"],
        "sample_title": row["sample_title"],
        "organism": row["species"],
        "source_url": row["source_url"],
    } for row in rows]
    write_table(
        root / "metadata_sources" / f"{name}.csv", registry,
        ["series_accession", "sample_accession", "sample_title", "organism", "source_url"], ",",
    )
    contract = [{key: row[key] for key in (
        "sample_accession", "experimental_time", "condition", "metadata_derivation"
    )} for row in rows]
    write_table(
        root / "metadata_contracts" / f"{name}.tsv", contract,
        ["sample_accession", "experimental_time", "condition", "metadata_derivation"], "\t",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    project = root.parent.parent
    snapshot = snapshot_counts(project)
    doc523 = doc523_counts(project)

    specifications = []

    wot_source = sample_sources(project, "02_GSE115943_GSE122662", "GSE115943")
    wot_dataset = "gse115943_pluripotency_reprogramming_core"
    wot_rows = build_rows(wot_source, "GSE115943", "Mus musculus", wot_dataset, snapshot, wot_metadata)
    specifications.append(("WOT", "GSE115943_production", "GSE115943", "Mus musculus", "WOT / GSE115943 complete 126-library QC", wot_rows, r"^gse115943_pluripotency_reprogramming_core_GSM[0-9]+$", "gse115943_wot_production", wot_source))

    gse103_source = sample_sources(project, "01_GSE103221", "GSE103221")
    gse103_dataset = "gse103221_pluripotency_reprogramming_core"
    gse103_counts = {**snapshot, **doc523}
    gse103_rows = build_rows(gse103_source, "GSE103221", "Mus musculus", gse103_dataset, gse103_counts, gse103221_metadata)
    specifications.append(("GSE103221", "GSE103221_production", "GSE103221", "Mus musculus", "GSE103221 complete 21-library QC", gse103_rows, r"^gse103221_pluripotency_reprogramming_core_GSM[0-9]+$", "gse103221_production", gse103_source))

    gse118_source = sample_sources(project, "03_GSE118258", "GSE118258")
    count_ids = {
        "GSM3323825": ("gse118258_D0", 7334),
        "GSM3323826": ("gse118258_D2", 5238),
        "GSM3323827": ("gse118258_D8", 4057),
        "GSM3323828": ("gse118258_D12", 7639),
        "GSM3323829": ("gse118258_D16_negative", 7029),
        "GSM3323830": ("gse118258_D16_positive", 3083),
    }
    gse118_counts = {library: count for library, count in count_ids.values()}
    gse118_rows = build_rows(
        gse118_source, "GSE118258", "Homo sapiens", "gse118258_human_pluripotency_reprogramming",
        gse118_counts, gse118258_metadata,
        library_name=lambda accession: count_ids[accession][0],
        path_builder=lambda accession, library: str(GSE118258_ROOT / library / "outs/filtered_feature_bc_matrix.h5"),
    )
    specifications.append(("GSE118258", "GSE118258_production", "GSE118258", "Homo sapiens", "GSE118258 complete 6-library QC", gse118_rows, r"^gse118258_D(?:0|2|8|12|16_negative|16_positive)$", "gse118258_production", gse118_source))

    gse150_source = sample_sources(project, "04_GSE150311", "GSE150311")
    gse150_dataset = "gse150311_human_pluripotency_reprogramming"
    gse150_rows = build_rows(gse150_source, "GSE150311", "Homo sapiens", gse150_dataset, snapshot, gse150311_metadata)
    specifications.append(("GSE150311", "GSE150311_production", "GSE150311", "Homo sapiens", "GSE150311 complete 4-library QC", gse150_rows, r"^gse150311_human_pluripotency_reprogramming_GSM[0-9]+$", "gse150311_production", gse150_source))

    gse171_source = sample_sources(project, "05_GSE171127", "GSE171127")
    gse171_dataset = "gse171127_pluripotency_vs_trophoblast_bifurcation"
    gse171_rows = build_rows(gse171_source, "GSE171127", "Mus musculus", gse171_dataset, snapshot, gse171127_metadata)
    specifications.append(("GSE171127", "GSE171127_production", "GSE171127", "Mus musculus", "GSE171127 complete 4-library QC", gse171_rows, r"^gse171127_pluripotency_vs_trophoblast_bifurcation_GSM[0-9]+$", "gse171127_production", gse171_source))

    registry_rows = []
    for scope, name, study_id, species, label, rows, pattern, output_dir, sources in specifications:
        write_study(root, name, rows, sources)
        config_path = write_config(root, name, study_id, species, label, rows, pattern, output_dir)
        registry_rows.append({
            "scope": scope, "config_path": config_path, "output_dir": output_dir,
            "expected_libraries": len(rows),
            "expected_called_cells": sum(int(row["expected_called_cells"]) for row in rows),
            "production_status": "pending",
        })

    registry_rows.append({
        "scope": "GSE178325", "config_path": "configs/GSE178325_regression.yaml",
        "output_dir": "gse178325_regression", "expected_libraries": 32,
        "expected_called_cells": 222726, "production_status": "recompute_with_exact_regression_gate",
    })
    write_table(
        root / "production/studies.tsv", registry_rows,
        ["scope", "config_path", "output_dir", "expected_libraries", "expected_called_cells", "production_status"], "\t",
    )
    total_libraries = sum(int(row["expected_libraries"]) for row in registry_rows)
    total_cells = sum(int(row["expected_called_cells"]) for row in registry_rows)
    if total_libraries != 193 or total_cells != 778207:
        raise ValueError(f"Six-study scope drift: {total_libraries} libraries, {total_cells} cells")
    print(f"Built six-study production inputs: {total_libraries} libraries, {total_cells} called cells")


if __name__ == "__main__":
    main()
