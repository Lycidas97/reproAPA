import hashlib
import json
import os
import pickle
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from ..utils.gtf import build_reference_intervals, read_gtf
from ..utils.genome import normalize_chromosome_column, validate_chrom_mapping
from ..utils.external_pas import compute_external_pas_overlap, summarize_external_status
from ..utils.sequence_qc import compute_sequence_qc
from .pas_qc import write_pas_annotation_qc
from .neighbor_policy_contract import (
    CORE_TES_EXONIC_AUDIT_STATE,
    CORE_TES_POLICY_CLASS,
    CORE_TES_STATE,
    INCOMPLETE_TES_STATE,
    INELIGIBLE_TES_POLICY_CLASS,
    OUT_OF_RANGE_TES_STATE,
    QUALIFIED_TES2KB_POLICY_CLASS,
    classify_selected_tes,
)


PAS_BED_COLUMNS = ["chrom", "start", "end", "name", "score", "strand"]
REGION_TABLES = [
    "utr3",
    "terminal_exon",
    "tez_250bp",
    "tez_1kb",
    "tez_2kb",
    "terminal_intron",
    "internal_intron",
    "internal_exon",
    "cds",
    "gene_downstream_2kb",
]

CANONICAL_TERMINAL_REGION_CLASSES = {"annotated_3utr_pas", "terminal_exon_pas"}
INTRON_REGION_CLASSES = {"terminal_intron_pas", "internal_intron_pas"}
EXONIC_CONFLICT_REGION_CLASSES = {"coding_exon_pas", "internal_exon_pas"}
TES_REGION_CLASS = "transcript_end_zone_pas"
TES_POLICY_MIN_DISTANCE_BP = 250
TES_POLICY_MAX_DISTANCE_BP = 2000
GENE_BIOTYPE_PRIORITY = {"protein_coding": 0, "lncRNA": 1, "lncRNA_gene": 1}

NEIGHBOR_DECISION_COLUMNS = [
    "same_strand_other_gene_intron_overlap",
    "same_strand_exonic_gene_conflict",
    "opposite_strand_gene_overlap",
    "neighbor_same_strand_gene_ids",
    "neighbor_same_strand_region_classes",
    "neighbor_opposite_strand_gene_ids",
    "neighbor_opposite_strand_region_classes",
    "selected_tes_candidate_id",
    "selected_tes_policy_class",
    "selected_tes_policy_reason",
    "neighbor_gene_decision_rule",
    "neighbor_gene_terminal_state",
    "neighbor_gene_exclusion_reason",
    "displaced_tes_candidate_ids",
    "displaced_tes_gene_ids",
    "displaced_tes_transcript_ids",
    "structural_priority_feature_family",
    "feature_first_candidate_id",
    "feature_first_gene_id",
    "biotype_priority_candidate_id",
    "biotype_priority_gene_id",
    "biotype_assignment_conflict",
    "biotype_tie_break_status",
]


@dataclass
class StructuralAnnotationResult:
    output_dir: str
    quantification_bed_df: pd.DataFrame
    pas_feature_map: pd.DataFrame
    pas_feature_annotation: pd.DataFrame


def annotate_pas_structural(
    filtered_pas_bed: str,
    gtf_path: str,
    fasta_path: str,
    chrom_map_path: str,
    output_dir: str,
    annotation_config: Optional[dict] = None,
) -> StructuralAnnotationResult:
    """Run the first-pass structural PAS annotation and write TSV artifacts."""
    config = validate_annotation_config(annotation_config or {})
    os.makedirs(output_dir, exist_ok=True)

    species = config["species"]
    genome_build = config["genome_build"]
    cut_site_window_bp = int(config.get("cut_site_window_bp", 25))
    transcript_end_zones_bp = config.get("transcript_end_zones_bp", [250, 1000, 2000])
    gene_downstream_bp = int(config.get("gene_downstream_bp", 2000))

    peaks = read_pas_bed(filtered_pas_bed)
    pas_map = normalize_pas_features(peaks, species=species, genome_build=genome_build)
    quant_bed = peaks.copy()
    quant_bed["name"] = pas_map["pas_id"]

    intervals, gtf_row_count, cache_status, cache_key = _load_reference_intervals(
        gtf_path,
        chrom_map_path,
        transcript_end_zones_bp=transcript_end_zones_bp,
        gene_downstream_bp=gene_downstream_bp,
        cache_dir=config.get("reference_cache_dir"),
    )
    if config.get("write_reference_intervals", True):
        _write_reference_intervals(intervals, os.path.join(output_dir, "reference_intervals"))

    candidates = generate_annotation_candidates(pas_map, intervals, cut_site_window_bp=cut_site_window_bp)
    feature_annotation = choose_primary_assignments(pas_map, candidates)
    decision = feature_annotation[["pas_id", *NEIGHBOR_DECISION_COLUMNS]].copy()
    candidates = candidates.drop(columns=NEIGHBOR_DECISION_COLUMNS, errors="ignore").merge(
        decision,
        on="pas_id",
        how="left",
        validate="many_to_one",
    )
    external = compute_external_pas_overlap(
        pas_map,
        config.get("external_pas_resources", []),
        annotation_species=species,
        default_window_bp=int(config.get("external_pas_window_bp", 50)),
    )
    candidates = _attach_external_to_candidates(candidates, external)
    feature_annotation = _attach_external_support(feature_annotation, external)
    if config.get("sequence_qc_enabled", True):
        fasta_index_path = str(
            config.get("sequence_qc_fasta_index_path")
            or os.path.join(output_dir, "sequence_qc_reference.fai")
        )
        sequence = compute_sequence_qc(
            feature_annotation[["pas_id", "chrom", "cut_site", "strand"]],
            fasta_path=fasta_path,
            fasta_index_path=fasta_index_path,
            motif_window_upstream_bp=tuple(config.get("motif_window_upstream_bp", [10, 50])),
            motif_search_downstream_bp=int(config.get("motif_search_downstream_bp", 200)),
            downstream_bp=int(config.get("internal_priming_downstream_bp", 20)),
            max_a_fraction=float(config.get("internal_priming_max_a_fraction", 0.8)),
            max_a_run=int(config.get("internal_priming_max_a_run", 8)),
        )
    else:
        sequence = compute_sequence_qc(feature_annotation[["pas_id"]], fasta_path=None)
    feature_annotation = _attach_sequence_qc(feature_annotation, sequence)
    confidence = build_confidence_tiers(feature_annotation)
    gene_pas_map = build_gene_pas_map(feature_annotation, confidence)
    gene_testability = build_gene_testability(feature_annotation, confidence)

    _write_outputs(
        output_dir=output_dir,
        pas_map=pas_map,
        candidates=candidates,
        feature_annotation=feature_annotation,
        external=external,
        sequence=sequence,
        confidence=confidence,
        gene_pas_map=gene_pas_map,
        gene_testability=gene_testability,
        manifest=_run_manifest(
            config,
            gtf_row_count,
            intervals,
            gtf_path,
            fasta_path,
            chrom_map_path,
            cache_status,
            cache_key,
        ),
    )
    if config.get("write_qc_outputs", True):
        write_pas_annotation_qc(
            output_dir=output_dir,
            annotation_dir=output_dir,
            write_plots=bool(config.get("write_qc_plots", True)),
        )
    return StructuralAnnotationResult(
        output_dir=output_dir,
        quantification_bed_df=quant_bed,
        pas_feature_map=pas_map,
        pas_feature_annotation=feature_annotation,
    )


def validate_annotation_config(config: dict) -> dict:
    """Validate the standalone annotation contract and return a shallow copy."""
    config = dict(config or {})
    for key in ["species", "genome_build"]:
        value = str(config.get(key, "")).strip()
        if not value:
            raise ValueError(f"Missing required annotation parameter: annotation.{key}")
        config[key] = value

    cut_window = int(config.get("cut_site_window_bp", 25))
    if cut_window < 0:
        raise ValueError("annotation.cut_site_window_bp must be >= 0")

    tez = [int(value) for value in config.get("transcript_end_zones_bp", [250, 1000, 2000])]
    if len(tez) != 3 or any(value <= 0 for value in tez) or tez != sorted(set(tez)):
        raise ValueError("annotation.transcript_end_zones_bp must contain 3 strictly increasing positive values")

    motif_window = [int(value) for value in config.get("motif_window_upstream_bp", [10, 50])]
    if len(motif_window) != 2 or min(motif_window) < 0 or motif_window[0] == motif_window[1]:
        raise ValueError("annotation.motif_window_upstream_bp must contain 2 distinct non-negative values")
    if int(config.get("motif_search_downstream_bp", 200)) < 0:
        raise ValueError("annotation.motif_search_downstream_bp must be >= 0")
    if int(config.get("internal_priming_downstream_bp", 20)) <= 0:
        raise ValueError("annotation.internal_priming_downstream_bp must be > 0")
    max_fraction = float(config.get("internal_priming_max_a_fraction", 0.8))
    if not 0 <= max_fraction <= 1:
        raise ValueError("annotation.internal_priming_max_a_fraction must be between 0 and 1")
    if int(config.get("internal_priming_max_a_run", 8)) < 0:
        raise ValueError("annotation.internal_priming_max_a_run must be >= 0")

    for resource in config.get("external_pas_resources", []):
        if isinstance(resource, str):
            path = resource
            status = "available"
        elif isinstance(resource, dict):
            status = str(resource.get("status", "available")).lower()
            if status not in {"available", "sparse", "unavailable"}:
                raise ValueError(f"Unsupported external PAS resource status: {status}")
            path = str(resource.get("path", ""))
        else:
            raise ValueError("annotation.external_pas_resources entries must be paths or mappings")
        if status in {"available", "sparse"} and (not path or not os.path.exists(path)):
            raise FileNotFoundError(f"External PAS resource not found: {path}")
    return config


def _file_signature(path: str) -> dict:
    resolved = Path(path).resolve()
    stat = resolved.stat()
    return {
        "path": str(resolved),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def _load_reference_intervals(
    gtf_path: str,
    chrom_map_path: str,
    transcript_end_zones_bp,
    gene_downstream_bp: int,
    cache_dir: Optional[str] = None,
) -> tuple[dict[str, pd.DataFrame], int, str, str]:
    payload = {
        "schema": 3,
        "gtf": _file_signature(gtf_path),
        "chromosome_map": _file_signature(chrom_map_path),
        "transcript_end_zones_bp": [int(value) for value in transcript_end_zones_bp],
        "gene_downstream_bp": int(gene_downstream_bp),
    }
    cache_key = hashlib.sha256(
        json.dumps(payload, sort_keys=True).encode("utf-8")
    ).hexdigest()[:20]
    cache_path = None
    if cache_dir:
        root = Path(cache_dir)
        root.mkdir(parents=True, exist_ok=True)
        cache_path = root / f"reference_intervals_{cache_key}.pickle"
        if cache_path.exists():
            try:
                with cache_path.open("rb") as handle:
                    cached = pickle.load(handle)
                if cached.get("payload") == payload:
                    return (
                        cached["intervals"],
                        int(cached["gtf_row_count"]),
                        "hit",
                        cache_key,
                    )
            except (EOFError, OSError, pickle.PickleError, AttributeError, KeyError):
                pass

    gtf = read_gtf(
        gtf_path,
        features=[
            "gene",
            "transcript",
            "exon",
            "CDS",
            "three_prime_utr",
            "five_prime_utr",
            "3UTR",
        ],
    )
    gtf = _apply_chrom_mapping_to_gtf(gtf, chrom_map_path)
    intervals = build_reference_intervals(
        gtf,
        transcript_end_zones_bp=transcript_end_zones_bp,
        gene_downstream_bp=gene_downstream_bp,
    )
    if cache_path is not None:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=cache_path.parent,
            prefix=f".{cache_path.name}.",
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            pickle.dump(
                {
                    "payload": payload,
                    "gtf_row_count": len(gtf),
                    "intervals": intervals,
                },
                handle,
                protocol=pickle.HIGHEST_PROTOCOL,
            )
        os.replace(temporary, cache_path)
    return intervals, len(gtf), "miss" if cache_path is not None else "disabled", cache_key


def finalize_pas_annotation(
    structural_dir: str,
    quantification_status_path: str,
    output_dir: Optional[str] = None,
    annotation_config: Optional[dict] = None,
) -> None:
    """Join quantification support into first-pass PAS annotation outputs."""
    output_dir = output_dir or structural_dir
    feature_path = os.path.join(structural_dir, "pas_feature_annotation.tsv")

    feature = pd.read_csv(feature_path, sep="\t", dtype=str).fillna("")
    if os.path.exists(quantification_status_path):
        status = pd.read_csv(quantification_status_path, sep="\t", dtype=str).fillna("")
        if "name" in status.columns:
            status = status.rename(columns={"name": "quantification_feature_id"})
            status["pas_id"] = _map_quantification_feature_ids(
                status["quantification_feature_id"],
                structural_dir,
            )
        keep_cols = [
            col
            for col in [
                "pas_id",
                "quantification_feature_id",
                "umi_count_pre_cell_filter",
                "cell_count_pre_cell_filter",
                "umi_count_post_cell_filter",
                "cell_count_post_cell_filter",
                "in_matrix",
                "drop_reason",
            ]
            if col in status.columns
        ]
        if "pas_id" in status.columns:
            status = status[keep_cols].drop_duplicates("pas_id")
            feature = feature.merge(status, on="pas_id", how="left")
            feature["total_pas_umi"] = _string_column(feature, "umi_count_pre_cell_filter", "", "")
            feature["n_cells_detected"] = _string_column(feature, "cell_count_pre_cell_filter", "", "")
            feature["quantification_status"] = _string_column(
                feature,
                "drop_reason",
                "quantification_status_unavailable",
                "not_observed",
            )
        else:
            feature["total_pas_umi"] = ""
            feature["n_cells_detected"] = ""
            feature["quantification_status"] = "quantification_status_unmatched"
    else:
        feature["quantification_status"] = "quantification_status_unavailable"

    confidence = build_confidence_tiers(feature)
    gene_pas_map = build_gene_pas_map(feature, confidence)
    gene_testability = build_gene_testability(feature, confidence)

    os.makedirs(output_dir, exist_ok=True)
    feature.to_csv(os.path.join(output_dir, "pas_feature_annotation.tsv"), sep="\t", index=False)
    confidence.to_csv(os.path.join(output_dir, "pas_confidence_tiers.tsv"), sep="\t", index=False)
    gene_pas_map.to_csv(os.path.join(output_dir, "gene_pas_map.tsv"), sep="\t", index=False)
    gene_testability.to_csv(os.path.join(output_dir, "gene_testability.tsv"), sep="\t", index=False)
    config = annotation_config or {}
    if config.get("write_qc_outputs", True):
        write_pas_annotation_qc(
            annotation_dir=output_dir,
            output_dir=os.path.join(output_dir, "qc"),
            write_plots=bool(config.get("write_qc_plots", True)),
        )


def read_pas_bed(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", header=None, names=PAS_BED_COLUMNS)


def normalize_pas_features(peaks: pd.DataFrame, species: str, genome_build: str) -> pd.DataFrame:
    rows = []
    for row in peaks.itertuples(index=False):
        cut_site = _cut_site(row)
        raw_peak_id = str(row.name)
        stable_key = "|".join(
            [
                str(species),
                str(genome_build),
                str(row.chrom),
                str(int(row.start)),
                str(int(row.end)),
                str(cut_site),
                str(row.strand),
                raw_peak_id,
            ]
        )
        cluster_id = f"pas_{hashlib.sha1(stable_key.encode('utf-8')).hexdigest()[:12]}"
        pas_id = f"{species}:{genome_build}:{row.chrom}:{cut_site}:{row.strand}:{cluster_id}"
        rows.append(
            {
                "raw_peak_id": raw_peak_id,
                "pas_id": pas_id,
                "species": species,
                "genome_build": genome_build,
                "chrom": row.chrom,
                "peak_start": int(row.start),
                "peak_end": int(row.end),
                "cut_site": cut_site,
                "strand": row.strand,
                "peak_width": int(row.end) - int(row.start),
                "score_field_original": row.score,
                "cut_site_source": "bed_score",
                "cluster_id": cluster_id,
            }
        )
    result = pd.DataFrame(rows)
    if not result.empty and result["pas_id"].duplicated().any():
        duplicates = result.loc[result["pas_id"].duplicated(keep=False), "raw_peak_id"].tolist()
        raise ValueError(f"Duplicate deterministic PAS IDs for input peaks: {duplicates[:5]}")
    return result


def generate_annotation_candidates(
    pas_map: pd.DataFrame,
    intervals: dict[str, pd.DataFrame],
    cut_site_window_bp: int = 25,
) -> pd.DataFrame:
    candidate_rows = []
    pas_ids_with_candidates = set()
    interval_index = _build_interval_index(intervals)
    for pas in pas_map.itertuples(index=False):
        pas_start = max(0, int(pas.cut_site) - cut_site_window_bp)
        pas_end = int(pas.cut_site) + cut_site_window_bp + 1
        for table_name in REGION_TABLES:
            hits = _query_interval_index(
                interval_index.get(table_name, {}),
                pas.chrom,
                pas_start,
                pas_end,
                pas.strand,
                same_strand=True,
            )
            for _, hit in hits.iterrows():
                candidate_rows.append(_candidate_row(pas, hit, "same_strand", pas_start, pas_end))
                pas_ids_with_candidates.add(pas.pas_id)

        for table_name in REGION_TABLES:
            opposite_hits = _query_interval_index(
                interval_index.get(table_name, {}),
                pas.chrom,
                pas_start,
                pas_end,
                pas.strand,
                same_strand=False,
            )
            for _, hit in opposite_hits.iterrows():
                candidate_rows.append(
                    _candidate_row(pas, hit, "opposite_strand", pas_start, pas_end)
                )
                pas_ids_with_candidates.add(pas.pas_id)

        antisense_hits = _query_interval_index(
            interval_index.get("gene_body", {}),
            pas.chrom,
            pas_start,
            pas_end,
            pas.strand,
            same_strand=False,
        )
        for _, hit in antisense_hits.iterrows():
            hit = hit.copy()
            hit["candidate_region_class"] = "opposite_strand_gene_body_pas"
            hit["region_priority"] = 12
            candidate_rows.append(
                _candidate_row(pas, hit, "opposite_strand", pas_start, pas_end)
            )
            pas_ids_with_candidates.add(pas.pas_id)

        if pas.pas_id not in pas_ids_with_candidates:
            candidate_rows.append(_unassigned_candidate(pas))

    return annotate_neighbor_gene_conflicts(pd.DataFrame(candidate_rows))


def annotate_neighbor_gene_conflicts(candidates: pd.DataFrame) -> pd.DataFrame:
    """Attach deterministic PAS-level candidate sets and eligible TES competition provenance."""
    result = candidates.copy()
    aggregate_columns = [
        "candidate_gene_ids",
        "candidate_region_classes",
        "same_strand_candidate_gene_ids",
        "same_strand_candidate_ids",
        "same_strand_candidate_transcript_ids",
        "same_strand_candidate_gene_biotypes",
        "same_strand_candidate_region_classes",
        "opposite_strand_candidate_gene_ids",
        "opposite_strand_candidate_region_classes",
        "original_tes_candidate_ids",
        "original_tes_gene_ids",
        "original_tes_transcript_ids",
        "valid_same_strand_core_tes_candidate_ids",
        "valid_same_strand_core_tes_gene_ids",
        "valid_same_strand_core_tes_transcript_ids",
        "qualified_same_strand_tes_candidate_ids",
        "qualified_same_strand_tes_gene_ids",
        "qualified_same_strand_tes_transcript_ids",
    ]
    count_columns = [
        "valid_same_strand_core_tes_candidate_count",
        "valid_same_strand_core_tes_gene_count",
        "qualified_same_strand_tes_candidate_count",
        "qualified_same_strand_tes_gene_count",
    ]
    if result.empty:
        for column in aggregate_columns:
            result[column] = pd.Series(dtype=object)
        for column in count_columns:
            result[column] = pd.Series(dtype=int)
        result["qualified_same_strand_tes"] = pd.Series(dtype=object)
        return result

    required = {
        "pas_id",
        "candidate_id",
        "candidate_gene_id",
        "candidate_transcript_id",
        "candidate_transcript_end",
        "candidate_region_class",
        "candidate_strand",
        "strand_relation",
        "downstream_distance_to_transcript_end",
    }
    missing = sorted(required - set(result.columns))
    if missing:
        raise ValueError(f"Candidate table missing neighbor-policy columns: {missing}")
    for column in [
        "candidate_gene_biotype",
        "candidate_transcript_biotype",
        "candidate_model_transcript_end",
        "candidate_model_3p_distance",
    ]:
        if column not in result.columns:
            result[column] = ""

    same = result["strand_relation"].eq("same_strand")
    opposite = result["strand_relation"].isin({"antisense", "opposite_strand"})
    tes = same & result["candidate_region_class"].eq(TES_REGION_CLASS)
    directed = pd.to_numeric(
        result["downstream_distance_to_transcript_end"], errors="coerce"
    )
    complete = (
        result["candidate_id"].astype(str).str.len().gt(0)
        & result["candidate_gene_id"].astype(str).str.len().gt(0)
        & result["candidate_transcript_id"].astype(str).str.len().gt(0)
        & result["candidate_transcript_end"].astype(str).str.len().gt(0)
        & result["candidate_strand"].isin({"+", "-"})
    )
    qualified = (
        tes
        & complete
        & directed.gt(TES_POLICY_MIN_DISTANCE_BP)
        & directed.le(TES_POLICY_MAX_DISTANCE_BP)
    )
    core = (
        tes
        & complete
        & directed.ge(0)
        & directed.le(TES_POLICY_MIN_DISTANCE_BP)
    )
    result["qualified_same_strand_tes"] = qualified.map(
        {True: "true", False: "false"}
    )

    result["candidate_gene_ids"] = result["pas_id"].map(
        _group_join(result, pd.Series(True, index=result.index), "candidate_gene_id")
    ).fillna("")
    result["candidate_region_classes"] = result["pas_id"].map(
        _group_join(result, pd.Series(True, index=result.index), "candidate_region_class")
    ).fillna("")
    result["same_strand_candidate_gene_ids"] = result["pas_id"].map(
        _group_join(result, same, "candidate_gene_id")
    ).fillna("")
    result["same_strand_candidate_ids"] = result["pas_id"].map(
        _group_join(result, same, "candidate_id")
    ).fillna("")
    result["same_strand_candidate_transcript_ids"] = result["pas_id"].map(
        _group_join(result, same, "candidate_transcript_id")
    ).fillna("")
    result["same_strand_candidate_gene_biotypes"] = result["pas_id"].map(
        _group_join(result, same, "candidate_gene_biotype")
    ).fillna("")
    result["same_strand_candidate_region_classes"] = result["pas_id"].map(
        _group_join(result, same, "candidate_region_class")
    ).fillna("")
    result["opposite_strand_candidate_gene_ids"] = result["pas_id"].map(
        _group_join(result, opposite, "candidate_gene_id")
    ).fillna("")
    result["opposite_strand_candidate_region_classes"] = result["pas_id"].map(
        _group_join(result, opposite, "candidate_region_class")
    ).fillna("")
    result["original_tes_candidate_ids"] = result["pas_id"].map(
        _group_join(result, tes, "candidate_id")
    ).fillna("")
    result["original_tes_gene_ids"] = result["pas_id"].map(
        _group_join(result, tes, "candidate_gene_id")
    ).fillna("")
    result["original_tes_transcript_ids"] = result["pas_id"].map(
        _group_join(result, tes, "candidate_transcript_id")
    ).fillna("")
    result["valid_same_strand_core_tes_candidate_ids"] = result["pas_id"].map(
        _group_join(result, core, "candidate_id")
    ).fillna("")
    result["valid_same_strand_core_tes_gene_ids"] = result["pas_id"].map(
        _group_join(result, core, "candidate_gene_id")
    ).fillna("")
    result["valid_same_strand_core_tes_transcript_ids"] = result["pas_id"].map(
        _group_join(result, core, "candidate_transcript_id")
    ).fillna("")
    result["valid_same_strand_core_tes_candidate_count"] = result["pas_id"].map(
        result.loc[core].groupby("pas_id")["candidate_id"].nunique()
    ).fillna(0).astype(int)
    result["valid_same_strand_core_tes_gene_count"] = result["pas_id"].map(
        result.loc[core].groupby("pas_id")["candidate_gene_id"].nunique()
    ).fillna(0).astype(int)
    result["qualified_same_strand_tes_candidate_ids"] = result["pas_id"].map(
        _group_join(result, qualified, "candidate_id")
    ).fillna("")
    result["qualified_same_strand_tes_gene_ids"] = result["pas_id"].map(
        _group_join(result, qualified, "candidate_gene_id")
    ).fillna("")
    result["qualified_same_strand_tes_transcript_ids"] = result["pas_id"].map(
        _group_join(result, qualified, "candidate_transcript_id")
    ).fillna("")
    result["qualified_same_strand_tes_candidate_count"] = result["pas_id"].map(
        result.loc[qualified].groupby("pas_id")["candidate_id"].nunique()
    ).fillna(0).astype(int)
    result["qualified_same_strand_tes_gene_count"] = result["pas_id"].map(
        result.loc[qualified].groupby("pas_id")["candidate_gene_id"].nunique()
    ).fillna(0).astype(int)
    return result


def _group_join(frame: pd.DataFrame, mask: pd.Series, column: str) -> pd.Series:
    subset = frame.loc[mask & frame[column].astype(str).str.len().gt(0), ["pas_id", column]]
    if subset.empty:
        return pd.Series(dtype=object)
    return subset.groupby("pas_id", sort=False)[column].agg(
        lambda values: "|".join(sorted(set(values.astype(str))))
    )


def choose_primary_assignments(pas_map: pd.DataFrame, candidates: pd.DataFrame) -> pd.DataFrame:
    if "candidate_gene_ids" not in candidates.columns:
        candidates = annotate_neighbor_gene_conflicts(candidates)
    rows = []
    empty_candidates = candidates.iloc[0:0].copy()
    candidate_groups = {pas_id: group.copy() for pas_id, group in candidates.groupby("pas_id", sort=False)}
    for pas in pas_map.itertuples(index=False):
        pas_candidates = candidate_groups.get(pas.pas_id, empty_candidates).copy()
        informative = pas_candidates[
            pas_candidates["strand_relation"].eq("same_strand")
            & ~pas_candidates["candidate_region_class"].isin(
                ["unassigned_pas", "intergenic_distal_pas"]
            )
        ].copy()
        if informative.empty:
            antisense = pas_candidates[
                pas_candidates["strand_relation"].isin({"antisense", "opposite_strand"})
            ]
            if not antisense.empty:
                primary = antisense.sort_values("region_priority").iloc[0]
                primary = primary.copy()
                primary["candidate_region_class"] = "antisense_pas"
                status = "excluded_antisense"
                confidence = "none"
                reason = "antisense_only"
            else:
                primary = pas_candidates.iloc[0]
                status = "unassigned"
                confidence = "none"
                reason = "no_reliable_gene_candidate"
        else:
            informative["region_priority"] = pd.to_numeric(informative["region_priority"], errors="coerce")
            canonical = _valid_canonical_terminal_candidates(informative)
            core_tes = _valid_same_strand_core_tes_candidates(informative)
            qualified_tes = _eligible_same_strand_tes_candidates(informative)
            other_exons = _same_strand_feature_family_candidates(
                informative, EXONIC_CONFLICT_REGION_CLASSES
            )
            introns = _same_strand_feature_family_candidates(
                informative, INTRON_REGION_CLASSES
            )
            if not canonical.empty:
                selection_pool = canonical
                top_priority = selection_pool["region_priority"].min()
                top = selection_pool[
                    selection_pool["region_priority"] == top_priority
                ]
                if top["candidate_gene_id"].nunique() > 1:
                    primary = top.sort_values("candidate_id").iloc[0]
                    status = "ambiguous_multi_gene"
                    confidence = "none"
                    reason = "multiple_same_priority_gene_candidates"
                    primary = primary.copy()
                    primary["candidate_region_class"] = "multi_gene_ambiguous_pas"
                else:
                    primary = top.sort_values(
                        ["distance_to_transcript_end", "overlap_bases"],
                        ascending=[True, False],
                    ).iloc[0]
                    if _terminal_exon_cds_overlap_without_utr(primary, pas_candidates):
                        primary = primary.copy()
                        primary["candidate_region_class"] = "terminal_exon_cds_overlap_pas"
                        status = "assigned_sensitivity"
                        confidence = "low"
                        reason = "terminal_exon_cds_overlap_retained_for_sensitivity"
                    else:
                        status, confidence, reason = _assignment_status(primary)
            elif not core_tes.empty:
                primary, status, confidence, reason = _select_tes_tier(
                    core_tes, tier="core"
                )
            elif not other_exons.empty:
                primary, status, confidence, reason = _select_structural_tier(
                    other_exons, family="other_exon"
                )
                if not qualified_tes.empty and status != "ambiguous_multi_gene":
                    reason = "other_exon_over_qualified_tes2kb"
            elif not introns.empty:
                terminal_introns = introns[
                    introns["candidate_region_class"].eq("terminal_intron_pas")
                ]
                structural_pool = terminal_introns if not terminal_introns.empty else introns
                primary, status, confidence, reason = _select_structural_tier(
                    structural_pool, family="intron"
                )
                if not qualified_tes.empty and status != "ambiguous_multi_gene":
                    reason = "intron_over_qualified_tes2kb"
            elif not qualified_tes.empty:
                primary, status, confidence, reason = _select_tes_tier(
                    qualified_tes, tier="qualified_rescue"
                )
            else:
                selection_pool = informative
                top_priority = selection_pool["region_priority"].min()
                top = selection_pool[
                    selection_pool["region_priority"] == top_priority
                ]
                if top["candidate_gene_id"].nunique() > 1:
                    primary = top.sort_values("candidate_id").iloc[0]
                    status = "ambiguous_multi_gene"
                    confidence = "none"
                    reason = "multiple_same_priority_gene_candidates"
                    primary = primary.copy()
                    primary["candidate_region_class"] = "multi_gene_ambiguous_pas"
                else:
                    primary = top.sort_values(
                        ["distance_to_transcript_end", "overlap_bases"],
                        ascending=[True, False],
                    ).iloc[0]
                    if _terminal_exon_cds_overlap_without_utr(primary, pas_candidates):
                        primary = primary.copy()
                        primary["candidate_region_class"] = "terminal_exon_cds_overlap_pas"
                        status = "assigned_sensitivity"
                        confidence = "low"
                        reason = "terminal_exon_cds_overlap_retained_for_sensitivity"
                    else:
                        status, confidence, reason = _assignment_status(primary)

        selected_tes_policy_class, selected_tes_policy_reason = classify_selected_tes(
            region_class=primary.get("candidate_region_class", ""),
            candidate_id=primary.get("candidate_id", ""),
            gene_id=primary.get("candidate_gene_id", ""),
            transcript_id=primary.get("candidate_transcript_id", ""),
            transcript_end=primary.get("candidate_transcript_end", ""),
            transcript_strand=primary.get("candidate_strand", ""),
            strand_relation=primary.get("strand_relation", ""),
            pas_strand=pas.strand,
            downstream_distance=primary.get(
                "downstream_distance_to_transcript_end", ""
            ),
        )
        if (
            primary.get("candidate_region_class", "") == TES_REGION_CLASS
            and selected_tes_policy_class == CORE_TES_POLICY_CLASS
            and status != "ambiguous_multi_gene_tes"
        ):
            status = "assigned_primary"
            confidence = "medium"
            reason = "core_tes_at_or_below_250bp"
        elif (
            primary.get("candidate_region_class", "") == TES_REGION_CLASS
            and selected_tes_policy_class == INELIGIBLE_TES_POLICY_CLASS
        ):
            status = "assigned_sensitivity"
            confidence = "low"
            reason = f"policy_ineligible_tes:{selected_tes_policy_reason}"

        neighbor = _neighbor_assignment_provenance(
            primary,
            pas_candidates,
            status=status,
            selected_tes_policy_class=selected_tes_policy_class,
            selected_tes_policy_reason=selected_tes_policy_reason,
        )

        rows.append(
            {
                **pas._asdict(),
                "primary_candidate_id": primary.get("candidate_id", ""),
                "primary_gene_id": primary.get("candidate_gene_id", ""),
                "primary_gene_symbol": primary.get("candidate_gene_symbol", ""),
                "primary_gene_biotype": primary.get("candidate_gene_biotype", ""),
                "primary_transcript_id": primary.get("candidate_transcript_id", ""),
                "primary_transcript_biotype": primary.get(
                    "candidate_transcript_biotype", ""
                ),
                "primary_model_transcript_end": primary.get(
                    "candidate_model_transcript_end", ""
                ),
                "primary_model_3p_distance": primary.get(
                    "candidate_model_3p_distance", ""
                ),
                "primary_transcript_strand": primary.get("candidate_strand", ""),
                "primary_strand_relation": primary.get("strand_relation", ""),
                "primary_transcript_end": primary.get("candidate_transcript_end", ""),
                "distance_to_primary_transcript_end": primary.get("distance_to_transcript_end", ""),
                "signed_distance_to_primary_transcript_end": primary.get(
                    "signed_distance_to_transcript_end",
                    "",
                ),
                "downstream_distance_to_primary_transcript_end": primary.get(
                    "downstream_distance_to_transcript_end",
                    "",
                ),
                "primary_region_class": primary.get("candidate_region_class", ""),
                "primary_assignment_status": status,
                "primary_assignment_confidence": confidence,
                "primary_assignment_reason": reason,
                "transcript_model_conflict": str(_has_transcript_conflict(pas_candidates)).lower(),
                "transcript_context_class": _transcript_context_class(
                    pas_candidates
                ),
                "multi_gene_ambiguous": str(
                    status in {"ambiguous_multi_gene", "ambiguous_multi_gene_tes"}
                ).lower(),
                "candidate_gene_ids": primary.get("candidate_gene_ids", ""),
                "candidate_region_classes": primary.get("candidate_region_classes", ""),
                "same_strand_candidate_gene_ids": primary.get(
                    "same_strand_candidate_gene_ids", ""
                ),
                "same_strand_candidate_ids": primary.get(
                    "same_strand_candidate_ids", ""
                ),
                "same_strand_candidate_transcript_ids": primary.get(
                    "same_strand_candidate_transcript_ids", ""
                ),
                "same_strand_candidate_gene_biotypes": primary.get(
                    "same_strand_candidate_gene_biotypes", ""
                ),
                "same_strand_candidate_region_classes": primary.get(
                    "same_strand_candidate_region_classes", ""
                ),
                "opposite_strand_candidate_gene_ids": primary.get(
                    "opposite_strand_candidate_gene_ids", ""
                ),
                "opposite_strand_candidate_region_classes": primary.get(
                    "opposite_strand_candidate_region_classes", ""
                ),
                "original_tes_candidate_ids": primary.get("original_tes_candidate_ids", ""),
                "original_tes_gene_ids": primary.get("original_tes_gene_ids", ""),
                "original_tes_transcript_ids": primary.get(
                    "original_tes_transcript_ids", ""
                ),
                "valid_same_strand_core_tes_candidate_ids": primary.get(
                    "valid_same_strand_core_tes_candidate_ids", ""
                ),
                "valid_same_strand_core_tes_gene_ids": primary.get(
                    "valid_same_strand_core_tes_gene_ids", ""
                ),
                "valid_same_strand_core_tes_transcript_ids": primary.get(
                    "valid_same_strand_core_tes_transcript_ids", ""
                ),
                "valid_same_strand_core_tes_candidate_count": int(
                    primary.get("valid_same_strand_core_tes_candidate_count", 0)
                    or 0
                ),
                "valid_same_strand_core_tes_gene_count": int(
                    primary.get("valid_same_strand_core_tes_gene_count", 0) or 0
                ),
                "qualified_same_strand_tes_candidate_ids": primary.get(
                    "qualified_same_strand_tes_candidate_ids", ""
                ),
                "qualified_same_strand_tes_gene_ids": primary.get(
                    "qualified_same_strand_tes_gene_ids", ""
                ),
                "qualified_same_strand_tes_transcript_ids": primary.get(
                    "qualified_same_strand_tes_transcript_ids", ""
                ),
                "qualified_same_strand_tes_candidate_count": int(
                    primary.get("qualified_same_strand_tes_candidate_count", 0) or 0
                ),
                "qualified_same_strand_tes_gene_count": int(
                    primary.get("qualified_same_strand_tes_gene_count", 0) or 0
                ),
                **neighbor,
                "internal_priming_flag": "false",
                "external_support_status": "unavailable",
                "total_pas_umi": "",
                "n_cells_detected": "",
                "quantification_status": "not_joined",
            }
        )
    return pd.DataFrame(rows)


def _eligible_same_strand_tes_candidates(candidates: pd.DataFrame) -> pd.DataFrame:
    directed = pd.to_numeric(
        candidates["downstream_distance_to_transcript_end"], errors="coerce"
    )
    complete = (
        candidates["candidate_id"].astype(str).str.len().gt(0)
        & candidates["candidate_gene_id"].astype(str).str.len().gt(0)
        & candidates["candidate_transcript_id"].astype(str).str.len().gt(0)
        & candidates["candidate_transcript_end"].astype(str).str.len().gt(0)
        & candidates["candidate_strand"].isin({"+", "-"})
    )
    return candidates[
        candidates["strand_relation"].eq("same_strand")
        & candidates["candidate_region_class"].eq(TES_REGION_CLASS)
        & complete
        & directed.gt(TES_POLICY_MIN_DISTANCE_BP)
        & directed.le(TES_POLICY_MAX_DISTANCE_BP)
    ].copy()


def _valid_canonical_terminal_candidates(candidates: pd.DataFrame) -> pd.DataFrame:
    complete = (
        candidates["candidate_id"].astype(str).str.len().gt(0)
        & candidates["candidate_gene_id"].astype(str).str.len().gt(0)
        & candidates["candidate_transcript_id"].astype(str).str.len().gt(0)
        & candidates["candidate_strand"].isin({"+", "-"})
    )
    return candidates[
        candidates["strand_relation"].eq("same_strand")
        & candidates["candidate_region_class"].isin(
            CANONICAL_TERMINAL_REGION_CLASSES
        )
        & complete
    ].copy()


def _valid_same_strand_core_tes_candidates(candidates: pd.DataFrame) -> pd.DataFrame:
    directed = pd.to_numeric(
        candidates["downstream_distance_to_transcript_end"], errors="coerce"
    )
    complete = (
        candidates["candidate_id"].astype(str).str.len().gt(0)
        & candidates["candidate_gene_id"].astype(str).str.len().gt(0)
        & candidates["candidate_transcript_id"].astype(str).str.len().gt(0)
        & candidates["candidate_transcript_end"].astype(str).str.len().gt(0)
        & candidates["candidate_strand"].isin({"+", "-"})
    )
    return candidates[
        candidates["strand_relation"].eq("same_strand")
        & candidates["candidate_region_class"].eq(TES_REGION_CLASS)
        & complete
        & directed.ge(0)
        & directed.le(TES_POLICY_MIN_DISTANCE_BP)
    ].copy()


def _select_tes_tier(
    candidates: pd.DataFrame, *, tier: str
) -> tuple[pd.Series, str, str, str]:
    """Select the nearest candidate inside one already-prioritized TES tier."""
    directed = pd.to_numeric(
        candidates["downstream_distance_to_transcript_end"], errors="raise"
    )
    minimum = directed.min()
    nearest = candidates.loc[directed.eq(minimum)]
    if nearest["candidate_gene_id"].nunique() > 1:
        primary = nearest.sort_values("candidate_id").iloc[0].copy()
        return (
            primary,
            "ambiguous_multi_gene_tes",
            "none",
            f"exact_nearest_{tier}_tes_distance_tie_across_genes",
        )
    primary = nearest.sort_values(
        ["overlap_bases", "candidate_id"], ascending=[False, True]
    ).iloc[0]
    status, confidence, reason = _assignment_status(primary)
    return primary, status, confidence, reason


def _same_strand_feature_family_candidates(
    candidates: pd.DataFrame, region_classes: set[str]
) -> pd.DataFrame:
    return candidates[
        candidates["strand_relation"].eq("same_strand")
        & candidates["candidate_region_class"].isin(region_classes)
        & candidates["candidate_gene_id"].astype(str).str.len().gt(0)
    ].copy()


def _gene_biotype_rank(value: object) -> int:
    biotype = str(value).strip()
    if biotype == "protein_coding":
        return 0
    if biotype in {"lncRNA", "lncRNA_gene", "lincRNA"}:
        return 1
    return 2


def _model_3p_distance(candidates: pd.DataFrame) -> pd.Series:
    model_distance = pd.to_numeric(
        candidates.get("candidate_model_3p_distance", ""), errors="coerce"
    )
    legacy_distance = pd.to_numeric(
        candidates["distance_to_transcript_end"], errors="coerce"
    )
    return model_distance.fillna(legacy_distance).fillna(float("inf"))


def _select_structural_tier(
    candidates: pd.DataFrame, *, family: str
) -> tuple[pd.Series, str, str, str]:
    """Select within one structural family by biotype, then 3' distance."""
    ranked = candidates.copy()
    ranked["_biotype_rank"] = ranked["candidate_gene_biotype"].map(
        _gene_biotype_rank
    )
    ranked["_model_3p_distance"] = _model_3p_distance(ranked)
    best_rank = ranked["_biotype_rank"].min()
    best_biotype = ranked.loc[ranked["_biotype_rank"].eq(best_rank)]
    minimum_distance = best_biotype["_model_3p_distance"].min()
    nearest = best_biotype.loc[
        best_biotype["_model_3p_distance"].eq(minimum_distance)
    ]
    if nearest["candidate_gene_id"].nunique() > 1:
        primary = nearest.sort_values("candidate_id").iloc[0].copy()
        return (
            primary,
            "ambiguous_multi_gene",
            "none",
            f"exact_{family}_biotype_and_3p_distance_tie_across_genes",
        )
    primary = nearest.sort_values(
        ["overlap_bases", "candidate_id"], ascending=[False, True]
    ).iloc[0]
    status, confidence, reason = _assignment_status(primary)
    return primary, status, confidence, reason


def _structural_selection_views(
    primary: pd.Series, candidates: pd.DataFrame, selection_status: str
) -> dict[str, str]:
    region = str(primary.get("candidate_region_class", ""))
    if region in EXONIC_CONFLICT_REGION_CLASSES:
        family = "other_exon"
        pool = _same_strand_feature_family_candidates(
            candidates, EXONIC_CONFLICT_REGION_CLASSES
        )
    elif region in INTRON_REGION_CLASSES:
        family = "intron"
        pool = _same_strand_feature_family_candidates(candidates, INTRON_REGION_CLASSES)
        terminal = pool[pool["candidate_region_class"].eq("terminal_intron_pas")]
        if not terminal.empty:
            pool = terminal
    else:
        return {
            "structural_priority_feature_family": "not_applicable",
            "feature_first_candidate_id": str(primary.get("candidate_id", "")),
            "feature_first_gene_id": str(primary.get("candidate_gene_id", "")),
            "biotype_priority_candidate_id": str(primary.get("candidate_id", "")),
            "biotype_priority_gene_id": str(primary.get("candidate_gene_id", "")),
            "biotype_assignment_conflict": "false",
            "biotype_tie_break_status": "guarded_canonical_core_or_nonstructural",
        }

    feature_pool = pool.copy()
    feature_pool["_model_3p_distance"] = _model_3p_distance(feature_pool)
    minimum = feature_pool["_model_3p_distance"].min()
    feature_nearest = feature_pool.loc[
        feature_pool["_model_3p_distance"].eq(minimum)
    ]
    if feature_nearest["candidate_gene_id"].nunique() == 1:
        feature_row = feature_nearest.sort_values(
            ["overlap_bases", "candidate_id"], ascending=[False, True]
        ).iloc[0]
        feature_candidate = str(feature_row["candidate_id"])
        feature_gene = str(feature_row["candidate_gene_id"])
    else:
        feature_candidate = ""
        feature_gene = ""

    biotype_candidate = str(primary.get("candidate_id", ""))
    biotype_gene = str(primary.get("candidate_gene_id", ""))
    conflict = bool(feature_gene and biotype_gene and feature_gene != biotype_gene)
    ranks = pool["candidate_gene_biotype"].map(_gene_biotype_rank)
    tie_status = (
        "biotype_selected"
        if ranks.nunique() > 1
        else "distance_selected"
    )
    if selection_status == "ambiguous_multi_gene":
        biotype_candidate = ""
        biotype_gene = ""
        tie_status = "exact_cross_gene_tie"
    return {
        "structural_priority_feature_family": family,
        "feature_first_candidate_id": feature_candidate,
        "feature_first_gene_id": feature_gene,
        "biotype_priority_candidate_id": biotype_candidate,
        "biotype_priority_gene_id": biotype_gene,
        "biotype_assignment_conflict": str(conflict).lower(),
        "biotype_tie_break_status": tie_status,
    }


def _has_other_gene_exonic_conflict(primary: pd.Series, candidates: pd.DataFrame) -> bool:
    gene_id = str(primary.get("candidate_gene_id", ""))
    conflict = candidates[
        candidates["strand_relation"].eq("same_strand")
        & candidates["candidate_gene_id"].astype(str).str.len().gt(0)
        & candidates["candidate_gene_id"].astype(str).ne(gene_id)
        & candidates["candidate_region_class"].isin(EXONIC_CONFLICT_REGION_CLASSES)
    ]
    return not conflict.empty


def _neighbor_assignment_provenance(
    primary: pd.Series,
    candidates: pd.DataFrame,
    status: str,
    selected_tes_policy_class: str,
    selected_tes_policy_reason: str,
) -> dict:
    selected_gene = str(primary.get("candidate_gene_id", ""))
    same_other = candidates[
        candidates["strand_relation"].eq("same_strand")
        & candidates["candidate_gene_id"].astype(str).str.len().gt(0)
        & candidates["candidate_gene_id"].astype(str).ne(selected_gene)
    ]
    opposite = candidates[
        candidates["strand_relation"].isin({"antisense", "opposite_strand"})
        & candidates["candidate_gene_id"].astype(str).str.len().gt(0)
    ]
    intron = same_other[
        same_other["candidate_region_class"].isin(INTRON_REGION_CLASSES)
    ]
    exonic = same_other[
        same_other["candidate_region_class"].isin(EXONIC_CONFLICT_REGION_CLASSES)
    ]
    qualified_tes = _eligible_same_strand_tes_candidates(candidates)
    structural_primary = primary.get("candidate_region_class", "") in (
        EXONIC_CONFLICT_REGION_CLASSES | INTRON_REGION_CLASSES
    )
    if structural_primary:
        displaced_tes = qualified_tes[
            qualified_tes["candidate_id"].astype(str).ne(
                str(primary.get("candidate_id", ""))
            )
        ]
    else:
        displaced_tes = qualified_tes.iloc[0:0]
    selection_views = _structural_selection_views(primary, candidates, status)
    competing_tes = qualified_tes["candidate_gene_id"].nunique() > 1
    canonical_over_tes = (
        primary.get("candidate_region_class", "") in CANONICAL_TERMINAL_REGION_CLASSES
        and not qualified_tes[
            qualified_tes["candidate_gene_id"].astype(str).ne(selected_gene)
        ].empty
    )

    if (
        status == "ambiguous_multi_gene"
        and selection_views["structural_priority_feature_family"]
        in {"other_exon", "intron"}
    ):
        rule = "select_structural_feature_biotype_distance_exact_tie"
        terminal = "ambiguous_multi_gene_structural"
        exclusion = "exact_feature_biotype_3p_distance_tie_across_genes"
    elif status == "ambiguous_multi_gene_tes":
        rule = "select_nearest_same_strand_tes_exact_multi_gene_tie"
        terminal = "ambiguous_multi_gene_tes"
        exclusion = "exact_minimum_tes_distance_tied_across_genes"
    elif canonical_over_tes:
        rule = "canonical_terminal_other_gene_over_tes2kb"
        terminal = "canonical_terminal_other_gene_over_tes2kb"
        exclusion = ""
    elif (
        selected_tes_policy_class == CORE_TES_POLICY_CLASS
        and not exonic.empty
    ):
        rule = "core_tes_exonic_overlap_audit_only"
        terminal = CORE_TES_EXONIC_AUDIT_STATE
        exclusion = ""
    elif selected_tes_policy_class == CORE_TES_POLICY_CLASS:
        rule = "core_tes_assignment_preserved"
        terminal = CORE_TES_STATE
        exclusion = ""
    elif (
        selected_tes_policy_class == INELIGIBLE_TES_POLICY_CLASS
        and selected_tes_policy_reason == "distance_above_2000bp"
    ):
        rule = "tes_candidate_outside_rescue_distance"
        terminal = OUT_OF_RANGE_TES_STATE
        exclusion = "tes_candidate_outside_rescue_distance"
    elif selected_tes_policy_class == INELIGIBLE_TES_POLICY_CLASS:
        rule = "tes_candidate_incomplete_provenance"
        terminal = INCOMPLETE_TES_STATE
        exclusion = "tes_candidate_incomplete_provenance"
    elif (
        primary.get("candidate_region_class", "") in EXONIC_CONFLICT_REGION_CLASSES
        and not qualified_tes.empty
    ):
        rule = "same_strand_other_exon_over_qualified_tes2kb"
        terminal = "other_exon_over_tes2kb"
        exclusion = "other_exon_over_tes2kb"
    elif (
        primary.get("candidate_region_class", "") in INTRON_REGION_CLASSES
        and not qualified_tes.empty
    ):
        rule = "same_strand_intron_over_qualified_tes2kb"
        terminal = "intron_over_tes2kb"
        exclusion = "intron_over_tes2kb"
    elif (
        selected_tes_policy_class == QUALIFIED_TES2KB_POLICY_CLASS
        and not exonic.empty
    ):
        rule = "same_strand_exonic_gene_conflict_to_sensitivity"
        terminal = "sensitivity_same_strand_exonic_gene_conflict"
        exclusion = "same_strand_exonic_gene_conflict"
    elif selected_tes_policy_class == QUALIFIED_TES2KB_POLICY_CLASS and competing_tes:
        rule = "select_unique_nearest_same_strand_tes"
        terminal = "tes2kb_rescue_nearest_competing_tes"
        exclusion = ""
    elif selected_tes_policy_class == QUALIFIED_TES2KB_POLICY_CLASS and not opposite.empty:
        rule = "tes2kb_by_pas_transcript_strand_over_opposite_gene"
        terminal = "tes2kb_rescue_opposite_strand_gene_overlap"
        exclusion = ""
    elif selected_tes_policy_class == QUALIFIED_TES2KB_POLICY_CLASS:
        rule = "tes2kb_no_disqualifying_neighbor_conflict"
        terminal = "tes2kb_rescue_no_neighbor_conflict"
        exclusion = ""
    else:
        rule = "region_priority_assignment"
        terminal = status
        exclusion = "" if status in {"assigned_primary", "assigned_sensitivity"} else status

    return {
        "same_strand_other_gene_intron_overlap": str(not intron.empty).lower(),
        "same_strand_exonic_gene_conflict": str(not exonic.empty).lower(),
        "opposite_strand_gene_overlap": str(not opposite.empty).lower(),
        "neighbor_same_strand_gene_ids": _sorted_join(same_other["candidate_gene_id"]),
        "neighbor_same_strand_region_classes": _sorted_join(
            same_other["candidate_region_class"]
        ),
        "neighbor_opposite_strand_gene_ids": _sorted_join(opposite["candidate_gene_id"]),
        "neighbor_opposite_strand_region_classes": _sorted_join(
            opposite["candidate_region_class"]
        ),
        "selected_tes_candidate_id": (
            str(primary.get("candidate_id", ""))
            if primary.get("candidate_region_class", "") == TES_REGION_CLASS
            else ""
        ),
        "selected_tes_policy_class": selected_tes_policy_class,
        "selected_tes_policy_reason": selected_tes_policy_reason,
        "neighbor_gene_decision_rule": rule,
        "neighbor_gene_terminal_state": terminal,
        "neighbor_gene_exclusion_reason": exclusion,
        "displaced_tes_candidate_ids": _sorted_join(
            displaced_tes["candidate_id"]
        ),
        "displaced_tes_gene_ids": _sorted_join(
            displaced_tes["candidate_gene_id"]
        ),
        "displaced_tes_transcript_ids": _sorted_join(
            displaced_tes["candidate_transcript_id"]
        ),
        **selection_views,
    }


def _sorted_join(values: pd.Series) -> str:
    return "|".join(sorted(set(value for value in values.astype(str) if value)))


def build_confidence_tiers(feature_annotation: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in feature_annotation.itertuples(index=False):
        status = row.primary_assignment_status
        region = row.primary_region_class
        internal_priming = str(getattr(row, "internal_priming_flag", "false")).lower() == "true"
        if internal_priming and status in {"assigned_primary", "assigned_sensitivity"}:
            tier = "relaxed"
            reason = "internal_priming_risk_retained_for_sensitivity"
        elif status == "assigned_primary":
            tier = "primary"
            reason = f"{region};assignment_{row.primary_assignment_confidence}"
        elif status == "assigned_sensitivity":
            tier = "relaxed"
            reason = f"{region};sensitivity_assignment"
        elif status == "excluded_artifact":
            tier = "exclude"
            reason = getattr(row, "primary_assignment_reason", "") or status
        else:
            tier = "exclude"
            reason = status
        rows.append(
            {
                "pas_id": row.pas_id,
                "confidence_tier": tier,
                "keep_for_pas_quant": str(tier in {"strict", "primary", "relaxed"}).lower(),
                "keep_for_apa_primary": str(tier in {"strict", "primary"}).lower(),
                "keep_for_sensitivity": str(tier in {"strict", "primary", "relaxed"}).lower(),
                "tier_reason": reason,
                "failure_reason": "" if tier != "exclude" else reason,
            }
        )
    return pd.DataFrame(rows)


def build_gene_pas_map(feature_annotation: pd.DataFrame, confidence: pd.DataFrame) -> pd.DataFrame:
    data = feature_annotation.merge(confidence, on="pas_id", how="left")
    data = data[data["primary_assignment_status"].isin(["assigned_primary", "assigned_sensitivity"])].copy()
    if data.empty:
        return pd.DataFrame(
            columns=[
                "species",
                "gene_id",
                "gene_symbol",
                "pas_id",
                "primary_region_class",
                "confidence_tier",
                "pas_rank_proximal_to_distal",
                "pas_rank_distal_to_proximal",
                "relative_position_within_gene",
                "distance_to_stop_codon",
                "distance_to_transcript_end",
                "is_proximal_pas",
                "is_distal_pas",
                "is_primary_usage_pas",
                "is_sensitivity_pas",
            ]
        )

    rows = []
    for _, group in data.groupby("primary_gene_id", dropna=False):
        strand = group["strand"].iloc[0]
        ascending = strand == "+"
        ordered = group.sort_values("cut_site", ascending=ascending).reset_index(drop=True)
        n = len(ordered)
        for idx, row in ordered.iterrows():
            rank = idx + 1
            rows.append(
                {
                    "species": row["species"],
                    "gene_id": row["primary_gene_id"],
                    "gene_symbol": row["primary_gene_symbol"],
                    "pas_id": row["pas_id"],
                    "primary_region_class": row["primary_region_class"],
                    "confidence_tier": row["confidence_tier"],
                    "pas_rank_proximal_to_distal": rank,
                    "pas_rank_distal_to_proximal": n - idx,
                    "relative_position_within_gene": "",
                    "distance_to_stop_codon": "",
                    "distance_to_transcript_end": "",
                    "is_proximal_pas": str(rank == 1).lower(),
                    "is_distal_pas": str(rank == n).lower(),
                    "is_primary_usage_pas": str(row["confidence_tier"] in {"strict", "primary"}).lower(),
                    "is_sensitivity_pas": str(row["confidence_tier"] == "relaxed").lower(),
                }
            )
    return pd.DataFrame(rows)


def build_gene_testability(feature_annotation: pd.DataFrame, confidence: pd.DataFrame) -> pd.DataFrame:
    data = feature_annotation.merge(confidence, on="pas_id", how="left")
    assigned = data[data["primary_gene_id"].astype(str).str.len() > 0].copy()
    rows = []
    for gene_id, group in assigned.groupby("primary_gene_id", dropna=False):
        n_primary = int(group["confidence_tier"].isin(["strict", "primary"]).sum())
        n_relaxed = int((group["confidence_tier"] == "relaxed").sum())
        rows.append(
            {
                "species": group["species"].iloc[0],
                "gene_id": gene_id,
                "gene_symbol": group["primary_gene_symbol"].iloc[0],
                "n_pas_total": len(group),
                "n_pas_strict": int((group["confidence_tier"] == "strict").sum()),
                "n_pas_primary": n_primary,
                "n_pas_relaxed": n_relaxed,
                "n_pas_excluded": int((group["confidence_tier"] == "exclude").sum()),
                "has_multi_pas_primary": str(n_primary >= 2).lower(),
                "is_apa_testable_primary": str(n_primary >= 2).lower(),
                "is_apa_testable_relaxed": str((n_primary + n_relaxed) >= 2).lower(),
                "testability_failure_reason": "" if n_primary >= 2 else "fewer_than_two_primary_pas",
            }
        )
    return pd.DataFrame(rows)


def _apply_chrom_mapping_to_gtf(gtf: pd.DataFrame, chrom_map_path: str) -> pd.DataFrame:
    if not chrom_map_path or not os.path.exists(chrom_map_path):
        return gtf
    mapping = validate_chrom_mapping(pd.read_csv(chrom_map_path, sep="\t", dtype=str))
    normalized, _, _ = normalize_chromosome_column(
        gtf,
        mapping,
        source_column="gtf_chrom",
    )
    return normalized


def _cut_site(row) -> int:
    score = pd.to_numeric(row.score, errors="coerce")
    if pd.notna(score):
        return int(score)
    return int((int(row.start) + int(row.end)) // 2)


def _overlapping_rows(table: pd.DataFrame, chrom: str, start: int, end: int, strand: str, same_strand: bool) -> pd.DataFrame:
    if same_strand:
        strand_mask = table["strand"] == strand
    else:
        strand_mask = table["strand"] != strand
    return table[
        (table["chrom"] == chrom)
        & strand_mask
        & (pd.to_numeric(table["start"], errors="coerce") < end)
        & (pd.to_numeric(table["end"], errors="coerce") > start)
    ].copy()


def _build_interval_index(intervals: dict[str, pd.DataFrame]) -> dict[str, dict[tuple[str, str], tuple[pd.DataFrame, np.ndarray]]]:
    index = {}
    for table_name, table in intervals.items():
        if table.empty or not {"chrom", "start", "end", "strand"}.issubset(table.columns):
            index[table_name] = {}
            continue
        table = table.copy()
        table["start"] = pd.to_numeric(table["start"], errors="coerce")
        table["end"] = pd.to_numeric(table["end"], errors="coerce")
        table = table.dropna(subset=["start", "end"])
        table["start"] = table["start"].astype(int)
        table["end"] = table["end"].astype(int)
        lookup = {}
        for key, group in table.sort_values(["chrom", "strand", "start"]).groupby(["chrom", "strand"], dropna=False):
            group = group.sort_values("start").reset_index(drop=True)
            lookup[(str(key[0]), str(key[1]))] = (group, group["start"].to_numpy())
        index[table_name] = lookup
    return index


def _query_interval_index(
    lookup: dict[tuple[str, str], tuple[pd.DataFrame, np.ndarray]],
    chrom: str,
    start: int,
    end: int,
    strand: str,
    same_strand: bool,
) -> pd.DataFrame:
    if not lookup:
        return pd.DataFrame()
    if same_strand:
        groups = [lookup.get((str(chrom), str(strand)))]
    else:
        groups = [value for (group_chrom, group_strand), value in lookup.items() if group_chrom == str(chrom) and group_strand != str(strand)]

    hits = []
    for item in groups:
        if item is None:
            continue
        group, starts = item
        limit = int(np.searchsorted(starts, end, side="left"))
        if limit == 0:
            continue
        subset = group.iloc[:limit]
        subset = subset[subset["end"] > start]
        if not subset.empty:
            hits.append(subset)
    if not hits:
        return pd.DataFrame()
    return pd.concat(hits, ignore_index=True)


def _candidate_row(pas, hit: pd.Series, strand_relation: str, pas_start: int, pas_end: int) -> dict:
    overlap = max(0, min(pas_end, int(hit["end"])) - max(pas_start, int(hit["start"])))
    distance = _distance_to_interval_point(int(pas.cut_site), hit)
    signed_distance = _signed_distance_to_interval_point(int(pas.cut_site), hit)
    downstream_distance = _downstream_distance_to_interval_point(int(pas.cut_site), hit)
    model_transcript_end = pd.to_numeric(
        hit.get("model_transcript_end", ""), errors="coerce"
    )
    model_3p_distance = (
        ""
        if pd.isna(model_transcript_end)
        else abs(int(pas.cut_site) - int(model_transcript_end))
    )
    return {
        "pas_id": pas.pas_id,
        "candidate_id": _stable_candidate_id(pas.pas_id, hit),
        "candidate_gene_id": hit.get("gene_id", ""),
        "candidate_gene_symbol": hit.get("gene_symbol", ""),
        "candidate_gene_biotype": hit.get("gene_biotype", ""),
        "candidate_transcript_id": hit.get("transcript_id", ""),
        "candidate_transcript_biotype": hit.get("transcript_biotype", ""),
        "candidate_model_transcript_end": (
            "" if pd.isna(model_transcript_end) else int(model_transcript_end)
        ),
        "candidate_model_3p_distance": model_3p_distance,
        "candidate_strand": hit.get("strand", ""),
        "strand_relation": strand_relation,
        "candidate_transcript_end": hit.get("transcript_end", ""),
        "candidate_region_class": hit.get("candidate_region_class", ""),
        "region_priority": hit.get("region_priority", ""),
        "overlap_bases": overlap,
        "distance_to_transcript_end": distance,
        "signed_distance_to_transcript_end": signed_distance,
        "downstream_distance_to_transcript_end": downstream_distance,
        "distance_to_gene_end": "",
        "distance_to_external_pas": "",
        "external_pas_source": "",
        "candidate_assignment_confidence": "medium" if hit.get("region_priority", 99) in [2, 3, 4] else "low",
        "candidate_assignment_reason": f"{strand_relation}_{hit.get('candidate_region_class', '')}",
    }


def _unassigned_candidate(pas) -> dict:
    return {
        "pas_id": pas.pas_id,
        "candidate_id": f"{pas.pas_id}:unassigned",
        "candidate_gene_id": "",
        "candidate_gene_symbol": "",
        "candidate_gene_biotype": "",
        "candidate_transcript_id": "",
        "candidate_transcript_biotype": "",
        "candidate_model_transcript_end": "",
        "candidate_model_3p_distance": "",
        "candidate_strand": "",
        "strand_relation": "none",
        "candidate_transcript_end": "",
        "candidate_region_class": "unassigned_pas",
        "region_priority": 13,
        "overlap_bases": 0,
        "distance_to_transcript_end": "",
        "signed_distance_to_transcript_end": "",
        "downstream_distance_to_transcript_end": "",
        "distance_to_gene_end": "",
        "distance_to_external_pas": "",
        "external_pas_source": "",
        "candidate_assignment_confidence": "none",
        "candidate_assignment_reason": "no_overlap",
    }


def _stable_candidate_id(pas_id: str, hit: pd.Series) -> str:
    raw = "|".join(
        [
            pas_id,
            str(hit.get("gene_id", "")),
            str(hit.get("transcript_id", "")),
            str(hit.get("candidate_region_class", "")),
            str(hit.get("start", "")),
            str(hit.get("end", "")),
        ]
    )
    return hashlib.sha1(raw.encode("utf-8")).hexdigest()[:16]


def _distance_to_interval_point(cut_site: int, hit: pd.Series):
    if "transcript_end" in hit and pd.notna(hit.get("transcript_end")):
        return abs(cut_site - int(hit["transcript_end"]))
    start = int(hit.get("start", cut_site))
    end = int(hit.get("end", cut_site))
    if start <= cut_site <= end:
        return 0
    return min(abs(cut_site - start), abs(cut_site - end))


def _signed_distance_to_interval_point(cut_site: int, hit: pd.Series):
    """Return genomic-coordinate signed distance (cut site minus transcript TES)."""
    transcript_end = pd.to_numeric(hit.get("transcript_end", ""), errors="coerce")
    if pd.isna(transcript_end):
        return ""
    return int(cut_site - int(transcript_end))


def _downstream_distance_to_interval_point(cut_site: int, hit: pd.Series):
    """Return strand-directed TES distance, positive in the downstream direction."""
    signed_distance = _signed_distance_to_interval_point(cut_site, hit)
    strand = str(hit.get("strand", ""))
    if signed_distance == "" or strand not in {"+", "-"}:
        return ""
    return int(signed_distance if strand == "+" else -signed_distance)


def _assignment_status(primary: pd.Series) -> tuple[str, str, str]:
    priority = int(primary.get("region_priority", 99))
    region = primary.get("candidate_region_class", "")
    if priority <= 4:
        return "assigned_primary", "medium", f"same_strand_{region}"
    if priority in {5, 6, 7, 8, 9, 11}:
        return "assigned_sensitivity", "low", f"sensitivity_{region}"
    if priority == 10:
        return "assigned_sensitivity", "low", "coding_exon_candidate_retained_for_sensitivity"
    return "unassigned", "none", "no_primary_candidate"


def _terminal_exon_cds_overlap_without_utr(primary: pd.Series, candidates: pd.DataFrame) -> bool:
    if primary.get("candidate_region_class", "") != "terminal_exon_pas":
        return False
    if candidates.empty:
        return False

    transcript_id = primary.get("candidate_transcript_id", "")
    gene_id = primary.get("candidate_gene_id", "")
    if not transcript_id and not gene_id:
        return False

    if transcript_id:
        same_model = candidates[candidates["candidate_transcript_id"] == transcript_id]
    else:
        same_model = candidates[candidates["candidate_gene_id"] == gene_id]

    has_cds_overlap = (same_model["candidate_region_class"] == "coding_exon_pas").any()
    has_utr_support = (same_model["candidate_region_class"] == "annotated_3utr_pas").any()
    return bool(has_cds_overlap and not has_utr_support)


def _has_transcript_conflict(candidates: pd.DataFrame) -> bool:
    if candidates.empty:
        return False

    same_strand = candidates[
        candidates["strand_relation"].eq("same_strand")
        & candidates["candidate_gene_id"].astype(str).str.len().gt(0)
        & candidates["candidate_transcript_id"].astype(str).str.len().gt(0)
    ].copy()
    if same_strand.empty:
        return False

    same_strand["region_priority_num"] = pd.to_numeric(
        same_strand["region_priority"], errors="coerce"
    ).fillna(99)
    for _, gene_candidates in same_strand.groupby("candidate_gene_id", sort=False):
        transcript_families = set()
        for _, transcript_candidates in gene_candidates.groupby("candidate_transcript_id", sort=False):
            best_priority = transcript_candidates["region_priority_num"].min()
            best_classes = transcript_candidates.loc[
                transcript_candidates["region_priority_num"].eq(best_priority),
                "candidate_region_class",
            ]
            transcript_families.update(
                family
                for family in best_classes.map(_transcript_region_family)
                if family
            )
        if len(transcript_families) > 1:
            return True
    return False


def _transcript_context_class(candidates: pd.DataFrame) -> str:
    if candidates.empty:
        return "none"
    same_strand = candidates[
        candidates["strand_relation"].eq("same_strand")
        & candidates["candidate_gene_id"].astype(str).str.len().gt(0)
        & candidates["candidate_transcript_id"].astype(str).str.len().gt(0)
    ].copy()
    if same_strand.empty:
        return "none"
    same_strand["region_priority_num"] = pd.to_numeric(
        same_strand["region_priority"], errors="coerce"
    ).fillna(99)
    gene_families: dict[str, set[str]] = {}
    for gene_id, gene_candidates in same_strand.groupby(
        "candidate_gene_id", sort=False
    ):
        transcript_families: set[str] = set()
        for _, transcript_candidates in gene_candidates.groupby(
            "candidate_transcript_id", sort=False
        ):
            best_priority = transcript_candidates["region_priority_num"].min()
            best_classes = transcript_candidates.loc[
                transcript_candidates["region_priority_num"].eq(best_priority),
                "candidate_region_class",
            ]
            transcript_families.update(
                family
                for family in best_classes.map(_transcript_region_family)
                if family
            )
        gene_families[str(gene_id)] = transcript_families
        if {"terminal_3end", "intronic"}.issubset(transcript_families):
            return "alternative_terminal_intron"

    terminal_genes = {
        gene_id
        for gene_id, families in gene_families.items()
        if "terminal_3end" in families
    }
    intron_genes = {
        gene_id for gene_id, families in gene_families.items() if "intronic" in families
    }
    if any(left != right for left in terminal_genes for right in intron_genes):
        return "cross_gene_terminal_intron"
    if _has_transcript_conflict(candidates):
        return "alternative_transcript_context"
    return "none"


def _transcript_region_family(region: str) -> str:
    if region in {"annotated_3utr_pas", "terminal_exon_pas", "transcript_end_zone_pas"}:
        return "terminal_3end"
    if region in {"terminal_intron_pas", "internal_intron_pas"}:
        return "intronic"
    if region == "coding_exon_pas":
        return "coding"
    return ""


def _attach_sequence_qc(feature: pd.DataFrame, sequence: pd.DataFrame) -> pd.DataFrame:
    feature = feature.drop(columns=["internal_priming_flag"], errors="ignore")
    keep = ["pas_id", "internal_priming_flag"]
    if "sequence_qc_status" in sequence.columns:
        keep.append("sequence_qc_status")
    return feature.merge(sequence[keep], on="pas_id", how="left")


def _attach_external_support(feature: pd.DataFrame, external: pd.DataFrame) -> pd.DataFrame:
    feature = feature.copy()
    if external.empty:
        feature["external_support_status"] = "unavailable"
        return feature
    summary = (
        external.groupby("pas_id", dropna=False)["overlap_status"]
        .agg(summarize_external_status)
        .reset_index()
        .rename(columns={"overlap_status": "external_support_status_new"})
    )
    feature = feature.merge(summary, on="pas_id", how="left")
    feature["external_support_status"] = feature["external_support_status_new"].fillna("unavailable")
    return feature.drop(columns=["external_support_status_new"])


def _attach_external_to_candidates(candidates: pd.DataFrame, external: pd.DataFrame) -> pd.DataFrame:
    if candidates.empty or external.empty:
        return candidates
    evidence = external.copy()
    evidence["distance_num"] = pd.to_numeric(evidence["distance_to_source_pas"], errors="coerce")
    evidence["support_rank"] = evidence["overlap_status"].map(
        {"supported": 0, "absent": 1, "unavailable_or_sparse": 2, "unavailable": 3}
    ).fillna(4)
    evidence = evidence.sort_values(["pas_id", "support_rank", "distance_num", "source_name"], na_position="last")
    evidence = evidence.drop_duplicates("pas_id")
    lookup = evidence.set_index("pas_id")[["distance_to_source_pas", "source_name"]]
    result = candidates.copy()
    result["distance_to_external_pas"] = result["pas_id"].map(lookup["distance_to_source_pas"]).fillna("")
    result["external_pas_source"] = result["pas_id"].map(lookup["source_name"]).fillna("")
    return result


def _string_column(df: pd.DataFrame, column: str, missing_column: str, missing_value: str) -> pd.Series:
    if column not in df.columns:
        return pd.Series([missing_column] * len(df), index=df.index, dtype=object)
    return df[column].fillna(missing_value)


def _map_quantification_feature_ids(feature_ids: pd.Series, structural_dir: str) -> pd.Series:
    feature_map_path = os.path.join(structural_dir, "pas_feature_map.tsv")
    if not os.path.exists(feature_map_path):
        return feature_ids
    feature_map = pd.read_csv(feature_map_path, sep="\t", dtype=str).fillna("")
    if not {"raw_peak_id", "pas_id"}.issubset(feature_map.columns):
        return feature_ids
    mapper = dict(zip(feature_map["raw_peak_id"], feature_map["pas_id"]))
    return feature_ids.map(mapper).fillna(feature_ids)


def _run_manifest(
    config: dict,
    gtf_row_count: int,
    intervals: dict[str, pd.DataFrame],
    gtf_path: str,
    fasta_path: str,
    chrom_map_path: str,
    cache_status: str,
    cache_key: str,
) -> pd.DataFrame:
    gene_span_source = ""
    gene_body = intervals.get("gene_body", pd.DataFrame())
    if not gene_body.empty and "gene_span_source" in gene_body.columns:
        gene_span_source = ";".join(sorted(set(gene_body["gene_span_source"].astype(str))))
    rows = [
        ("species", config.get("species", "")),
        ("genome_build", config.get("genome_build", "")),
        ("gtf_path", gtf_path),
        ("fasta_path", fasta_path),
        ("chrom_map_path", chrom_map_path),
        ("gtf_input_mode", "gtf_gzip_infer" if str(gtf_path).endswith(".gz") else "gtf_plain"),
        ("gene_span_source", gene_span_source),
        ("pas_id_version", "coordinate_peak_hash_v1"),
        ("external_pas_resource_count", len(config.get("external_pas_resources", []))),
        ("external_pas_window_bp", config.get("external_pas_window_bp", 50)),
        ("motif_search_downstream_bp", config.get("motif_search_downstream_bp", 200)),
        ("n_gtf_rows_loaded", gtf_row_count),
        ("n_gene_intervals", len(intervals.get("gene_body", []))),
        ("n_transcript_intervals", len(intervals.get("transcript_body", []))),
        ("reference_cache_status", cache_status),
        ("reference_cache_key", cache_key),
    ]
    return pd.DataFrame(rows, columns=["key", "value"])


def _write_reference_intervals(intervals: dict[str, pd.DataFrame], output_dir: str) -> None:
    os.makedirs(output_dir, exist_ok=True)
    for name, table in intervals.items():
        table.to_csv(os.path.join(output_dir, f"{name}.tsv"), sep="\t", index=False)


def _write_outputs(
    output_dir: str,
    pas_map: pd.DataFrame,
    candidates: pd.DataFrame,
    feature_annotation: pd.DataFrame,
    external: pd.DataFrame,
    sequence: pd.DataFrame,
    confidence: pd.DataFrame,
    gene_pas_map: pd.DataFrame,
    gene_testability: pd.DataFrame,
    manifest: pd.DataFrame,
) -> None:
    tables = {
        "pas_annotation_run_manifest.tsv": manifest,
        "pas_feature_map.tsv": pas_map,
        "pas_feature_annotation.tsv": feature_annotation,
        "pas_annotation_candidates.tsv": candidates,
        "pas_external_overlap.tsv": external,
        "pas_sequence_qc.tsv": sequence,
        "pas_confidence_tiers.tsv": confidence,
        "gene_pas_map.tsv": gene_pas_map,
        "gene_testability.tsv": gene_testability,
    }
    for filename, table in tables.items():
        table.to_csv(os.path.join(output_dir, filename), sep="\t", index=False)
