import hashlib
import os
from dataclasses import dataclass
from numbers import Integral, Real
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from .quantification import quantify_pas
from .neighbor_policy_contract import (
    NEIGHBOR_POLICY_REQUIRED_COLUMNS,
    validate_neighbor_policy_row,
)


STRICT_REGION_CLASSES = {"annotated_3utr_pas", "terminal_exon_pas"}
RESCUE_REGION_CLASSES = {"transcript_end_zone_pas"}
DEFAULT_TEZ_RESCUE_MIN_DISTANCE_BP = 250
DEFAULT_TEZ_RESCUE_MAX_DISTANCE_BP = 2000
EXCLUDED_REGION_CLASSES = {
    "antisense_pas",
    "multi_gene_ambiguous_pas",
    "unassigned_pas",
}

NEIGHBOR_PROVENANCE_COLUMNS = [
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
    "valid_same_strand_core_tes_candidate_count",
    "valid_same_strand_core_tes_gene_count",
    "qualified_same_strand_tes_candidate_ids",
    "qualified_same_strand_tes_gene_ids",
    "qualified_same_strand_tes_transcript_ids",
    "qualified_same_strand_tes_candidate_count",
    "qualified_same_strand_tes_gene_count",
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
    "transcript_model_conflict",
    "transcript_context_class",
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

SAMPLE_EVIDENCE_COLUMNS = [
    "sample_id",
    "species",
    "genome_build",
    "sample_pas_id",
    "raw_peak_id",
    "sample_strand",
    "sample_primary_gene_id",
    "chrom",
    "start",
    "end",
    "cut_site",
    "strand",
    "primary_gene_id",
    "primary_gene_symbol",
    "primary_gene_biotype",
    "primary_transcript_id",
    "primary_transcript_biotype",
    "primary_model_transcript_end",
    "primary_model_3p_distance",
    "primary_candidate_id",
    "primary_transcript_strand",
    "primary_strand_relation",
    "primary_transcript_end",
    "distance_to_primary_transcript_end",
    "signed_distance_to_primary_transcript_end",
    "downstream_distance_to_primary_transcript_end",
    "primary_assignment_reason",
    "primary_region_class",
    "primary_assignment_status",
    "confidence_tier",
    "evidence_tier_for_atlas",
    "evidence_group",
    "evidence_policy_status",
    "internal_priming_flag",
    "sequence_qc_status",
    "best_pas_motif",
    "best_pas_motif_search_class",
    "external_overlap_status",
    "total_pas_umi",
    "n_cells_detected",
    "in_matrix",
    "quantification_status",
    "candidate_status",
    "main_atlas_eligible",
    "main_atlas_exclusion_reason",
    *NEIGHBOR_PROVENANCE_COLUMNS,
]

PAS_ATLAS_COLUMNS = [
    "atlas_pas_id",
    "species",
    "genome_build",
    "chrom",
    "strand",
    "representative_cut_site",
    "cut_site_min",
    "cut_site_max",
    "cut_site_mad",
    "atlas_start",
    "atlas_end",
    "gene_id",
    "gene_symbol",
    "atlas_region_class",
    "atlas_confidence_tier",
    "support_library_count",
    "support_sample_count",
    "member_peak_count",
    "total_member_umi",
    "total_member_cells",
    "has_motif_support",
    "has_expected_motif_support",
    "has_expanded_motif_support",
    "has_external_support",
    "broad_or_unstable_cut_site",
    "atlas_status",
    "has_same_strand_other_gene_intron_overlap",
    "has_same_strand_exonic_gene_conflict",
    "has_opposite_strand_gene_overlap",
    "neighbor_gene_decision_rules",
    "neighbor_gene_terminal_states",
    "has_transcript_model_conflict",
    "transcript_context_classes",
    "has_biotype_assignment_conflict",
    "biotype_tie_break_statuses",
]

SAMPLE_TO_ATLAS_COLUMNS = [
    "sample_id",
    "sample_pas_id",
    "raw_peak_id",
    "sample_strand",
    "sample_primary_gene_id",
    "atlas_pas_id",
    "distance_to_atlas_cut_site",
    "mapping_status",
    "sample_main_atlas_eligible",
    "sample_main_atlas_exclusion_reason",
    "sample_total_umi",
    "sample_n_cells",
    "sample_confidence_tier",
    "sample_region_class",
    "sample_evidence_tier",
    "sample_evidence_policy_status",
    "sample_primary_transcript_id",
    "sample_primary_gene_biotype",
    "sample_primary_transcript_biotype",
    "sample_primary_model_transcript_end",
    "sample_primary_model_3p_distance",
    "sample_primary_candidate_id",
    "sample_primary_transcript_strand",
    "sample_primary_strand_relation",
    "sample_primary_transcript_end",
    "sample_distance_to_primary_transcript_end",
    "sample_signed_distance_to_primary_transcript_end",
    "sample_downstream_distance_to_primary_transcript_end",
    "sample_primary_assignment_reason",
    "sample_primary_assignment_status",
    "atlas_confidence_tier",
    "atlas_region_class",
    *[f"sample_{column}" for column in NEIGHBOR_PROVENANCE_COLUMNS],
]

PAS_ATLAS_SAMPLE_COUNT_COLUMNS = [
    "atlas_pas_id",
    "sample_id",
    "atlas_total_umi_from_sample_pas",
    "atlas_n_cells_from_sample_pas",
    "mapped_sample_pas_count",
    "exact_sample_pas_count",
    "nearby_sample_pas_count",
]

ATLAS_REQUANT_WINDOW_COLUMNS = [
    "atlas_pas_id",
    "chrom",
    "window_start",
    "window_end",
    "strand",
    "gene_id",
    "representative_cut_site",
    "window_source",
    "window_status",
]

ATLAS_REQUANT_STATUS_COLUMNS = [
    "sample_id",
    "atlas_pas_id",
    "chrom",
    "window_start",
    "window_end",
    "strand",
    "window_status",
    "observed_umi_pre_cell_filter",
    "observed_cells_pre_cell_filter",
    "observed_umi_post_cell_filter",
    "observed_cells_post_cell_filter",
    "has_umi",
    "in_matrix",
    "drop_reason",
    "quantification_status",
]


@dataclass
class PasAtlasResult:
    output_dir: str
    sample_evidence: pd.DataFrame
    pas_atlas: pd.DataFrame
    sample_to_atlas_map: pd.DataFrame
    atlas_sample_counts: pd.DataFrame
    atlas_sample_umi_matrix: pd.DataFrame
    qc_summary: pd.DataFrame
    review_examples: pd.DataFrame


@dataclass
class PasAtlasRequantificationResult:
    output_dir: str
    atlas_windows: pd.DataFrame
    quantification_status: pd.DataFrame
    requantified_umi_matrix: pd.DataFrame
    requantified_cell_matrix: pd.DataFrame
    qc_summary: pd.DataFrame


def validate_pas_atlas_build_config(atlas_config: Optional[dict] = None) -> dict:
    """Validate and normalize the atlas-build evidence policy."""
    config = dict(atlas_config or {})
    enabled = config.get("promote_tez_sensitivity_to_rescue", False)
    if not isinstance(enabled, bool):
        raise ValueError("atlas.build.promote_tez_sensitivity_to_rescue must be true or false")
    raw_min = config.get(
        "tez_rescue_min_downstream_distance_bp",
        DEFAULT_TEZ_RESCUE_MIN_DISTANCE_BP,
    )
    raw_max = config.get(
        "tez_rescue_max_downstream_distance_bp",
        DEFAULT_TEZ_RESCUE_MAX_DISTANCE_BP,
    )
    min_distance = _coerce_integer_bp_bound(raw_min)
    max_distance = _coerce_integer_bp_bound(raw_max)
    if min_distance < 0:
        raise ValueError("atlas.build.tez_rescue_min_downstream_distance_bp must be >= 0")
    if max_distance <= min_distance:
        raise ValueError("atlas.build TES rescue maximum distance must be greater than the minimum")
    config["promote_tez_sensitivity_to_rescue"] = enabled
    config["tez_rescue_min_downstream_distance_bp"] = min_distance
    config["tez_rescue_max_downstream_distance_bp"] = max_distance
    return config


def _coerce_integer_bp_bound(raw_value) -> int:
    """Accept integer-valued bounds without silently truncating fractions."""
    error = ValueError("atlas.build TES rescue distance bounds must be integer base-pair values")
    if isinstance(raw_value, (bool, np.bool_)):
        raise error
    if isinstance(raw_value, Integral):
        return int(raw_value)
    if isinstance(raw_value, Real):
        numeric = float(raw_value)
        if not np.isfinite(numeric) or not numeric.is_integer():
            raise error
        return int(numeric)
    if isinstance(raw_value, str):
        text = raw_value.strip()
        if text and text.lstrip("+-").isdigit():
            return int(text)
    raise error


def build_pas_atlas(
    annotation_dirs: Sequence[str],
    output_dir: str,
    sample_ids: Optional[Sequence[str]] = None,
    atlas_config: Optional[dict] = None,
) -> PasAtlasResult:
    """Build a first-pass species-level PAS atlas from finalized annotation directories.

    The MVP treats per-sample PAS features as evidence, clusters usable evidence
    into atlas PAS features, and writes traceable sample-to-atlas mappings. It
    does not re-quantify BAMs; missing projected sample peaks should not be
    interpreted as true zero.
    """
    config = validate_pas_atlas_build_config(atlas_config)
    if sample_ids is not None and len(sample_ids) != len(annotation_dirs):
        raise ValueError("sample_ids must have the same length as annotation_dirs")

    os.makedirs(output_dir, exist_ok=True)
    evidence_frames = []
    for idx, annotation_dir in enumerate(annotation_dirs):
        sample_id = sample_ids[idx] if sample_ids is not None else _sample_id_from_annotation_dir(annotation_dir)
        evidence_frames.append(
            read_pas_sample_evidence(
                annotation_dir,
                sample_id=sample_id,
                atlas_config=config,
            )
        )

    sample_evidence = (
        pd.concat(evidence_frames, ignore_index=True)
        if evidence_frames
        else pd.DataFrame(columns=SAMPLE_EVIDENCE_COLUMNS)
    )
    sample_evidence = _ensure_columns(sample_evidence, SAMPLE_EVIDENCE_COLUMNS)
    identity = sample_evidence[["species", "genome_build"]].drop_duplicates()
    if len(identity) > 1 and not bool(config.get("allow_mixed_species_or_builds", False)):
        observed = identity.astype(str).agg(":".join, axis=1).tolist()
        raise ValueError(f"A PAS atlas must use one species/genome build, observed: {observed}")

    pas_atlas, sample_to_atlas = cluster_pas_atlas(sample_evidence, atlas_config=config)
    atlas_sample_counts, atlas_sample_umi_matrix = build_pas_atlas_sample_counts(pas_atlas, sample_to_atlas)
    qc_summary = summarize_pas_atlas_qc(
        sample_evidence,
        pas_atlas,
        sample_to_atlas,
        include_sensitivity=bool(config.get("include_sensitivity", False)),
    )
    review_examples = build_pas_atlas_review_examples(pas_atlas)

    tables = {
        "pas_sample_evidence.tsv": sample_evidence,
        "pas_atlas.tsv": pas_atlas,
        "pas_sample_to_atlas_map.tsv": sample_to_atlas,
        "pas_atlas_sample_counts.tsv": atlas_sample_counts,
        "pas_atlas_sample_umi_matrix.tsv": atlas_sample_umi_matrix,
        "pas_atlas_qc_summary.tsv": qc_summary,
        "pas_atlas_review_examples.tsv": review_examples,
    }
    for filename, table in tables.items():
        table.to_csv(os.path.join(output_dir, filename), sep="\t", index=False)

    return PasAtlasResult(
        output_dir=output_dir,
        sample_evidence=sample_evidence,
        pas_atlas=pas_atlas,
        sample_to_atlas_map=sample_to_atlas,
        atlas_sample_counts=atlas_sample_counts,
        atlas_sample_umi_matrix=atlas_sample_umi_matrix,
        qc_summary=qc_summary,
        review_examples=review_examples,
    )


def requantify_pas_atlas(
    pas_atlas: pd.DataFrame | str | Path,
    bam_paths: Sequence[str] | dict[str, str],
    output_dir: str,
    sample_ids: Optional[Sequence[str]] = None,
    valid_barcode_paths: Optional[Sequence[Optional[str]] | dict[str, Optional[str]]] = None,
    umi_tag: str = "UB",
    barcode_tag: str = "CB",
    n_jobs: int = 1,
    show_progress: bool = False,
    min_reads_per_cell: int = 50,
    transcript_strand_from_read: str = "same",
    quantification_backend: str = "auto",
    requant_config: Optional[dict] = None,
) -> PasAtlasRequantificationResult:
    """Re-quantify atlas PAS windows against each sample BAM.

    Overlapping windows for distinct PAS of the same gene are partitioned at the
    midpoint between representative cut sites. Remaining overlaps, including
    cross-gene windows or duplicate cut sites, are marked
    `window_overlap_ambiguous` and are not quantified to avoid double-counting.
    """
    config = requant_config or {}
    atlas = _load_pas_atlas_table(pas_atlas)
    os.makedirs(output_dir, exist_ok=True)

    sample_specs = _normalize_requant_sample_inputs(
        bam_paths=bam_paths,
        sample_ids=sample_ids,
        valid_barcode_paths=valid_barcode_paths,
    )
    atlas_windows = build_pas_atlas_quantification_windows(atlas, requant_config=config)

    status_frames = []
    sample_qc_rows = []
    per_sample_root = Path(output_dir) / "per_sample"
    per_sample_root.mkdir(parents=True, exist_ok=True)

    quant_windows = atlas_windows[atlas_windows["window_status"] == "usable"].copy()
    quant_peak_bed = _atlas_windows_to_peak_bed(quant_windows)

    for sample_id, bam_path, barcode_path in sample_specs:
        sample_dir = per_sample_root / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)

        frames = [_nonusable_requant_rows(atlas_windows, sample_id)]
        qc_row = {
            "sample_id": sample_id,
            "n_atlas_pas": len(atlas_windows),
            "n_usable_windows": int((atlas_windows["window_status"] == "usable").sum()),
            "n_window_overlap_ambiguous": int((atlas_windows["window_status"] == "window_overlap_ambiguous").sum()),
            "n_quantified_detected": 0,
            "n_quantified_zero": 0,
            "n_quantified_only_filtered_cells": 0,
            "n_not_quantified_missing_bam": 0,
            "n_not_quantified_missing_barcode": 0,
            "n_not_quantified_error": 0,
        }

        require_valid_barcodes = bool(config.get("require_valid_barcodes", False))
        if not bam_path or not os.path.exists(bam_path):
            missing = _placeholder_requant_rows(atlas_windows, sample_id, default_status="not_quantified_missing_bam")
            frames.append(missing[missing["window_status"] == "usable"].copy())
            qc_row["n_not_quantified_missing_bam"] = len(missing[missing["window_status"] == "usable"])
        elif require_valid_barcodes and (not barcode_path or not os.path.exists(barcode_path)):
            missing = _placeholder_requant_rows(atlas_windows, sample_id, default_status="not_quantified_missing_barcode")
            frames.append(missing[missing["window_status"] == "usable"].copy())
            qc_row["n_not_quantified_missing_barcode"] = len(missing[missing["window_status"] == "usable"])
        elif quant_peak_bed.empty:
            pass
        else:
            valid_barcodes = _read_valid_barcodes(barcode_path)
            try:
                quantify_pas(
                    bam_path=bam_path,
                    peak_bed_df=quant_peak_bed,
                    output_dir=str(sample_dir),
                    umi_tag=umi_tag,
                    barcode_tag=barcode_tag,
                    n_jobs=n_jobs,
                    show_progress=show_progress,
                    valid_barcodes=valid_barcodes,
                    min_reads_per_cell=min_reads_per_cell,
                    transcript_strand_from_read=transcript_strand_from_read,
                    quantification_backend=quantification_backend,
                )
                quantified = _read_requant_status(sample_dir, quant_windows, sample_id)
                frames.append(quantified)
                qc_row.update(_summarize_requant_status(quantified))
            except Exception:
                if bool(config.get("fail_on_error", False)):
                    raise
                errored = _placeholder_requant_rows(atlas_windows, sample_id, default_status="not_quantified_error")
                frames.append(errored[errored["window_status"] == "usable"].copy())
                qc_row["n_not_quantified_error"] = len(errored[errored["window_status"] == "usable"])

        combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=ATLAS_REQUANT_STATUS_COLUMNS)
        status_frames.append(combined[ATLAS_REQUANT_STATUS_COLUMNS].copy())
        sample_qc_rows.append(qc_row)

    quantification_status = (
        pd.concat(status_frames, ignore_index=True)
        if status_frames
        else pd.DataFrame(columns=ATLAS_REQUANT_STATUS_COLUMNS)
    )
    quantification_status = quantification_status[ATLAS_REQUANT_STATUS_COLUMNS].copy()
    requantified_umi_matrix = _build_requant_matrix(quantification_status, "observed_umi_post_cell_filter", sample_specs, atlas_windows)
    requantified_cell_matrix = _build_requant_matrix(quantification_status, "observed_cells_post_cell_filter", sample_specs, atlas_windows)
    qc_summary = pd.DataFrame(sample_qc_rows)

    atlas_windows.to_csv(os.path.join(output_dir, "pas_atlas_quantification_windows.tsv"), sep="\t", index=False)
    quantification_status.to_csv(os.path.join(output_dir, "pas_atlas_quantification_status.tsv"), sep="\t", index=False)
    requantified_umi_matrix.to_csv(os.path.join(output_dir, "pas_atlas_requantified_umi_matrix.tsv"), sep="\t", index=False)
    requantified_cell_matrix.to_csv(os.path.join(output_dir, "pas_atlas_requantified_cell_matrix.tsv"), sep="\t", index=False)
    qc_summary.to_csv(os.path.join(output_dir, "pas_atlas_requant_qc_summary.tsv"), sep="\t", index=False)

    return PasAtlasRequantificationResult(
        output_dir=output_dir,
        atlas_windows=atlas_windows,
        quantification_status=quantification_status,
        requantified_umi_matrix=requantified_umi_matrix,
        requantified_cell_matrix=requantified_cell_matrix,
        qc_summary=qc_summary,
    )


def read_pas_sample_evidence(
    annotation_dir: str,
    sample_id: Optional[str] = None,
    atlas_config: Optional[dict] = None,
) -> pd.DataFrame:
    """Read a finalized single-sample PAS annotation directory as evidence rows."""
    config = validate_pas_atlas_build_config(atlas_config)
    annotation_path = Path(annotation_dir)
    sample_id = sample_id or _sample_id_from_annotation_dir(annotation_dir)
    feature = _read_tsv(annotation_path / "pas_feature_annotation.tsv")
    if feature.empty:
        return pd.DataFrame(columns=SAMPLE_EVIDENCE_COLUMNS)

    confidence = _read_tsv(annotation_path / "pas_confidence_tiers.tsv")
    if not confidence.empty and "confidence_tier" not in feature.columns:
        keep = [col for col in ["pas_id", "confidence_tier", "tier_reason", "failure_reason"] if col in confidence.columns]
        feature = feature.merge(confidence[keep].drop_duplicates("pas_id"), on="pas_id", how="left")

    sequence = _read_tsv(annotation_path / "pas_sequence_qc.tsv")
    if not sequence.empty:
        keep = [
            col
            for col in [
                "pas_id",
                "best_pas_motif",
                "best_pas_motif_search_class",
                "sequence_qc_status",
                "internal_priming_flag",
            ]
            if col in sequence.columns
        ]
        sequence = sequence[keep].drop_duplicates("pas_id")
        suffix_cols = [col for col in keep if col != "pas_id" and col in feature.columns]
        merged = feature.merge(sequence, on="pas_id", how="left", suffixes=("", "_sequence"))
        for col in suffix_cols:
            seq_col = f"{col}_sequence"
            if seq_col in merged.columns:
                merged[col] = merged[col].where(merged[col].astype(str).str.len() > 0, merged[seq_col])
                merged = merged.drop(columns=[seq_col])
        feature = merged

    external = _read_tsv(annotation_path / "pas_external_overlap.tsv")
    if not external.empty and "external_overlap_status" not in feature.columns:
        external_status = (
            external.groupby("pas_id", dropna=False)["overlap_status"]
            .agg(_external_status_summary)
            .reset_index()
            .rename(columns={"overlap_status": "external_overlap_status"})
        )
        feature = feature.merge(external_status, on="pas_id", how="left")

    if bool(config.get("promote_tez_sensitivity_to_rescue", False)):
        required_neighbor_columns = set(NEIGHBOR_PROVENANCE_COLUMNS) | set(
            NEIGHBOR_POLICY_REQUIRED_COLUMNS
        )
        missing_neighbor_columns = sorted(
            required_neighbor_columns - set(feature.columns)
        )
        if missing_neighbor_columns:
            raise ValueError(
                "New atlas input is missing required neighbor provenance columns: "
                f"{missing_neighbor_columns}"
            )

    feature = _ensure_columns(
        feature,
        [
            "pas_id",
            "raw_peak_id",
            "species",
            "genome_build",
            "chrom",
            "peak_start",
            "peak_end",
            "cut_site",
            "strand",
            "primary_gene_id",
            "primary_gene_symbol",
            "primary_gene_biotype",
            "primary_transcript_id",
            "primary_transcript_biotype",
            "primary_model_transcript_end",
            "primary_model_3p_distance",
            "primary_candidate_id",
            "primary_transcript_strand",
            "primary_strand_relation",
            "primary_transcript_end",
            "distance_to_primary_transcript_end",
            "signed_distance_to_primary_transcript_end",
            "downstream_distance_to_primary_transcript_end",
            "primary_assignment_reason",
            "multi_gene_ambiguous",
            "primary_region_class",
            "primary_assignment_status",
            "confidence_tier",
            "internal_priming_flag",
            "sequence_qc_status",
            "best_pas_motif",
            "best_pas_motif_search_class",
            "external_overlap_status",
            "total_pas_umi",
            "n_cells_detected",
            "in_matrix",
            "quantification_status",
            *NEIGHBOR_PROVENANCE_COLUMNS,
        ],
    )

    rows = []
    for row in feature.itertuples(index=False):
        evidence_tier, candidate_status = _derive_evidence_tier(row)
        policy_status = _tez_sensitivity_rescue_policy_status(row, config)
        if policy_status == "promoted_tez_sensitivity_to_rescue":
            evidence_tier = "rescue"
            candidate_status = "usable_tez_sensitivity_rescue"
        main_atlas_eligible, main_atlas_exclusion_reason = _main_atlas_eligibility(
            row
        )
        if not main_atlas_eligible:
            evidence_tier = "sensitivity"
            candidate_status = (
                "usable_sensitivity_" + main_atlas_exclusion_reason
            )
        rows.append(
            {
                "sample_id": sample_id,
                "species": row.species,
                "genome_build": row.genome_build,
                "sample_pas_id": row.pas_id,
                "raw_peak_id": row.raw_peak_id,
                "chrom": row.chrom,
                "start": row.peak_start,
                "end": row.peak_end,
                "cut_site": row.cut_site,
                "strand": row.strand,
                "primary_gene_id": row.primary_gene_id,
                "primary_gene_symbol": row.primary_gene_symbol,
                "primary_gene_biotype": row.primary_gene_biotype,
                "primary_transcript_id": row.primary_transcript_id,
                "primary_transcript_biotype": row.primary_transcript_biotype,
                "primary_model_transcript_end": row.primary_model_transcript_end,
                "primary_model_3p_distance": row.primary_model_3p_distance,
                "primary_candidate_id": row.primary_candidate_id,
                "primary_transcript_strand": row.primary_transcript_strand,
                "primary_strand_relation": row.primary_strand_relation,
                "primary_transcript_end": row.primary_transcript_end,
                "distance_to_primary_transcript_end": row.distance_to_primary_transcript_end,
                "signed_distance_to_primary_transcript_end": row.signed_distance_to_primary_transcript_end,
                "downstream_distance_to_primary_transcript_end": row.downstream_distance_to_primary_transcript_end,
                "primary_assignment_reason": row.primary_assignment_reason,
                "primary_region_class": row.primary_region_class,
                "primary_assignment_status": row.primary_assignment_status,
                "confidence_tier": row.confidence_tier,
                "evidence_tier_for_atlas": evidence_tier,
                "evidence_group": _evidence_group(evidence_tier, row.primary_region_class),
                "evidence_policy_status": policy_status,
                "internal_priming_flag": _bool_string(row.internal_priming_flag),
                "sequence_qc_status": row.sequence_qc_status,
                "best_pas_motif": row.best_pas_motif,
                "best_pas_motif_search_class": row.best_pas_motif_search_class,
                "external_overlap_status": row.external_overlap_status,
                "total_pas_umi": row.total_pas_umi,
                "n_cells_detected": row.n_cells_detected,
                "in_matrix": row.in_matrix,
                "quantification_status": row.quantification_status,
                "candidate_status": candidate_status,
                "main_atlas_eligible": str(main_atlas_eligible).lower(),
                "main_atlas_exclusion_reason": main_atlas_exclusion_reason,
                **{
                    column: getattr(row, column)
                    for column in NEIGHBOR_PROVENANCE_COLUMNS
                },
            }
        )
    return pd.DataFrame(rows, columns=SAMPLE_EVIDENCE_COLUMNS)


def cluster_pas_atlas(sample_evidence: pd.DataFrame, atlas_config: Optional[dict] = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Cluster usable sample evidence rows into atlas PAS features and mappings."""
    config = validate_pas_atlas_build_config(atlas_config)
    include_sensitivity = bool(config.get("include_sensitivity", False))
    usable_tiers = {"strict", "rescue"}
    if include_sensitivity:
        usable_tiers.add("sensitivity")

    evidence = _ensure_columns(sample_evidence.copy(), SAMPLE_EVIDENCE_COLUMNS)
    if bool(config.get("promote_tez_sensitivity_to_rescue", False)):
        for row in evidence.itertuples(index=False):
            derived_eligible, derived_reason = _main_atlas_eligibility(row)
            observed_eligible = str(row.main_atlas_eligible).strip().lower()
            observed_reason = str(row.main_atlas_exclusion_reason).strip()
            if observed_eligible != str(derived_eligible).lower() or observed_reason != derived_reason:
                raise ValueError(
                    "Sample evidence main-atlas eligibility contradicts neighbor policy: "
                    f"sample={row.sample_id}, pas={row.sample_pas_id}"
                )
    evidence["cut_site_num"] = pd.to_numeric(evidence["cut_site"], errors="coerce")
    evidence = evidence.dropna(subset=["cut_site_num"]).copy()
    evidence["cut_site_num"] = evidence["cut_site_num"].astype(int)

    usable = evidence[
        evidence["evidence_tier_for_atlas"].isin(usable_tiers)
        & evidence["main_atlas_eligible"].astype(str).eq("true")
        & evidence["primary_gene_id"].astype(str).str.len().gt(0)
        & evidence["candidate_status"].astype(str).str.startswith("usable")
    ].copy()

    atlas_rows = []
    mapped_rows = []
    group_cols = ["species", "genome_build", "chrom", "strand", "primary_gene_id", "evidence_group"]
    if not usable.empty:
        usable = usable.sort_values(group_cols + ["cut_site_num", "sample_id", "sample_pas_id"]).reset_index(drop=True)
        for _, group in usable.groupby(group_cols, dropna=False, sort=False):
            for cluster in _cluster_group(group, config):
                atlas_row = _atlas_row(cluster, config)
                atlas_rows.append(atlas_row)
                for row in cluster.itertuples(index=False):
                    distance = abs(int(row.cut_site_num) - int(atlas_row["representative_cut_site"]))
                    mapped_rows.append(
                        {
                            "sample_id": row.sample_id,
                            "sample_pas_id": row.sample_pas_id,
                            "raw_peak_id": row.raw_peak_id,
                            "sample_strand": row.strand,
                            "sample_primary_gene_id": row.primary_gene_id,
                            "atlas_pas_id": atlas_row["atlas_pas_id"],
                            "distance_to_atlas_cut_site": distance,
                            "mapping_status": _mapping_status(distance, atlas_row["broad_or_unstable_cut_site"]),
                            "sample_main_atlas_eligible": row.main_atlas_eligible,
                            "sample_main_atlas_exclusion_reason": row.main_atlas_exclusion_reason,
                            "sample_total_umi": row.total_pas_umi,
                            "sample_n_cells": row.n_cells_detected,
                            "sample_confidence_tier": row.confidence_tier,
                            "sample_region_class": row.primary_region_class,
                            "sample_evidence_tier": row.evidence_tier_for_atlas,
                            "sample_evidence_policy_status": row.evidence_policy_status,
                            "sample_primary_transcript_id": row.primary_transcript_id,
                            "sample_primary_gene_biotype": row.primary_gene_biotype,
                            "sample_primary_transcript_biotype": row.primary_transcript_biotype,
                            "sample_primary_model_transcript_end": row.primary_model_transcript_end,
                            "sample_primary_model_3p_distance": row.primary_model_3p_distance,
                            "sample_primary_candidate_id": row.primary_candidate_id,
                            "sample_primary_transcript_strand": row.primary_transcript_strand,
                            "sample_primary_strand_relation": row.primary_strand_relation,
                            "sample_primary_transcript_end": row.primary_transcript_end,
                            "sample_distance_to_primary_transcript_end": row.distance_to_primary_transcript_end,
                            "sample_signed_distance_to_primary_transcript_end": row.signed_distance_to_primary_transcript_end,
                            "sample_downstream_distance_to_primary_transcript_end": row.downstream_distance_to_primary_transcript_end,
                            "sample_primary_assignment_reason": row.primary_assignment_reason,
                            "sample_primary_assignment_status": row.primary_assignment_status,
                            "atlas_confidence_tier": atlas_row["atlas_confidence_tier"],
                            "atlas_region_class": atlas_row["atlas_region_class"],
                            **{
                                f"sample_{column}": getattr(row, column)
                                for column in NEIGHBOR_PROVENANCE_COLUMNS
                            },
                        }
                    )

    mapped_ids = {row["sample_pas_id"] for row in mapped_rows}
    for row in evidence[~evidence["sample_pas_id"].isin(mapped_ids)].itertuples(index=False):
        mapped_rows.append(
            {
                "sample_id": row.sample_id,
                "sample_pas_id": row.sample_pas_id,
                "raw_peak_id": row.raw_peak_id,
                "sample_strand": row.strand,
                "sample_primary_gene_id": row.primary_gene_id,
                "atlas_pas_id": "",
                "distance_to_atlas_cut_site": "",
                "mapping_status": _unmapped_mapping_status(row, include_sensitivity),
                "sample_main_atlas_eligible": row.main_atlas_eligible,
                "sample_main_atlas_exclusion_reason": row.main_atlas_exclusion_reason,
                "sample_total_umi": row.total_pas_umi,
                "sample_n_cells": row.n_cells_detected,
                "sample_confidence_tier": row.confidence_tier,
                "sample_region_class": row.primary_region_class,
                "sample_evidence_tier": row.evidence_tier_for_atlas,
                "sample_evidence_policy_status": row.evidence_policy_status,
                "sample_primary_transcript_id": row.primary_transcript_id,
                "sample_primary_gene_biotype": row.primary_gene_biotype,
                "sample_primary_transcript_biotype": row.primary_transcript_biotype,
                "sample_primary_model_transcript_end": row.primary_model_transcript_end,
                "sample_primary_model_3p_distance": row.primary_model_3p_distance,
                "sample_primary_candidate_id": row.primary_candidate_id,
                "sample_primary_transcript_strand": row.primary_transcript_strand,
                "sample_primary_strand_relation": row.primary_strand_relation,
                "sample_primary_transcript_end": row.primary_transcript_end,
                "sample_distance_to_primary_transcript_end": row.distance_to_primary_transcript_end,
                "sample_signed_distance_to_primary_transcript_end": row.signed_distance_to_primary_transcript_end,
                "sample_downstream_distance_to_primary_transcript_end": row.downstream_distance_to_primary_transcript_end,
                "sample_primary_assignment_reason": row.primary_assignment_reason,
                "sample_primary_assignment_status": row.primary_assignment_status,
                "atlas_confidence_tier": "",
                "atlas_region_class": "",
                **{
                    f"sample_{column}": getattr(row, column)
                    for column in NEIGHBOR_PROVENANCE_COLUMNS
                },
            }
        )

    atlas = pd.DataFrame(atlas_rows, columns=PAS_ATLAS_COLUMNS)
    sample_to_atlas = pd.DataFrame(mapped_rows, columns=SAMPLE_TO_ATLAS_COLUMNS)
    return atlas, sample_to_atlas


def summarize_pas_atlas_qc(
    sample_evidence: pd.DataFrame,
    pas_atlas: pd.DataFrame,
    sample_to_atlas: pd.DataFrame,
    include_sensitivity: bool = False,
) -> pd.DataFrame:
    """Summarize basic atlas construction QC metrics."""
    if sample_evidence.empty:
        usable_mask = pd.Series(dtype=bool)
        clusterable_mask = pd.Series(dtype=bool)
    else:
        tier = sample_evidence["evidence_tier_for_atlas"].astype(str)
        usable_mask = sample_evidence["candidate_status"].astype(str).str.startswith("usable")
        clusterable_tiers = {"strict", "rescue"}
        if include_sensitivity:
            clusterable_tiers.add("sensitivity")
        clusterable_mask = usable_mask & tier.isin(clusterable_tiers)

    rows = [
        ("n_sample_evidence", len(sample_evidence)),
        ("n_samples", sample_evidence["sample_id"].nunique() if "sample_id" in sample_evidence else 0),
        ("n_usable_evidence", int(usable_mask.sum()) if not sample_evidence.empty else 0),
        ("n_clusterable_evidence", int(clusterable_mask.sum()) if not sample_evidence.empty else 0),
        ("n_atlas_pas", len(pas_atlas)),
        ("n_atlas_core", int((pas_atlas.get("atlas_status", pd.Series(dtype=str)) == "core").sum()) if not pas_atlas.empty else 0),
        ("n_atlas_private", int((pas_atlas.get("atlas_status", pd.Series(dtype=str)) == "private").sum()) if not pas_atlas.empty else 0),
        ("n_atlas_broad", _bool_count(pas_atlas.get("broad_or_unstable_cut_site", pd.Series(dtype=str))) if not pas_atlas.empty else 0),
        ("n_mapped_sample_pas", int(sample_to_atlas["atlas_pas_id"].astype(str).str.len().gt(0).sum()) if not sample_to_atlas.empty else 0),
        ("n_unmapped_sample_pas", int(sample_to_atlas["atlas_pas_id"].astype(str).str.len().eq(0).sum()) if not sample_to_atlas.empty else 0),
        ("n_not_clustered_sensitivity", int((sample_to_atlas.get("mapping_status", pd.Series(dtype=str)) == "not_clustered_sensitivity").sum()) if not sample_to_atlas.empty else 0),
    ]
    for tier in ["strict", "rescue", "sensitivity", "exclude"]:
        rows.append(
            (
                f"n_evidence_{tier}",
                int((sample_evidence.get("evidence_tier_for_atlas", pd.Series(dtype=str)) == tier).sum()) if not sample_evidence.empty else 0,
            )
        )
    return pd.DataFrame(rows, columns=["metric", "value"])


def build_pas_atlas_sample_counts(pas_atlas: pd.DataFrame, sample_to_atlas: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate mapped sample-PAS counts to atlas-PAS by sample.

    This is an evidence-level summary from existing sample PAS quantification,
    not BAM-level re-quantification over atlas windows.
    """
    if pas_atlas.empty or sample_to_atlas.empty:
        return (
            pd.DataFrame(columns=PAS_ATLAS_SAMPLE_COUNT_COLUMNS),
            pd.DataFrame(columns=["atlas_pas_id"]),
        )

    mapped = sample_to_atlas[sample_to_atlas["atlas_pas_id"].astype(str).str.len().gt(0)].copy()
    if mapped.empty:
        return (
            pd.DataFrame(columns=PAS_ATLAS_SAMPLE_COUNT_COLUMNS),
            pd.DataFrame({"atlas_pas_id": pas_atlas["atlas_pas_id"].astype(str)}),
        )

    mapped["sample_total_umi_num"] = _numeric(mapped["sample_total_umi"])
    mapped["sample_n_cells_num"] = _numeric(mapped["sample_n_cells"])
    mapped["exact_flag"] = (mapped["mapping_status"].astype(str) == "exact").astype(int)
    mapped["nearby_flag"] = (mapped["mapping_status"].astype(str) == "nearby").astype(int)
    counts = (
        mapped.groupby(["atlas_pas_id", "sample_id"], dropna=False)
        .agg(
            atlas_total_umi_from_sample_pas=("sample_total_umi_num", "sum"),
            atlas_n_cells_from_sample_pas=("sample_n_cells_num", "sum"),
            mapped_sample_pas_count=("sample_pas_id", "count"),
            exact_sample_pas_count=("exact_flag", "sum"),
            nearby_sample_pas_count=("nearby_flag", "sum"),
        )
        .reset_index()
    )
    for col in [
        "atlas_total_umi_from_sample_pas",
        "atlas_n_cells_from_sample_pas",
        "mapped_sample_pas_count",
        "exact_sample_pas_count",
        "nearby_sample_pas_count",
    ]:
        counts[col] = counts[col].round().astype(int)

    atlas_order = pas_atlas["atlas_pas_id"].astype(str).tolist()
    sample_order = sorted(counts["sample_id"].astype(str).unique())
    counts["atlas_pas_id"] = pd.Categorical(counts["atlas_pas_id"].astype(str), categories=atlas_order, ordered=True)
    counts["sample_id"] = pd.Categorical(counts["sample_id"].astype(str), categories=sample_order, ordered=True)
    counts = counts.sort_values(["atlas_pas_id", "sample_id"]).astype({"atlas_pas_id": str, "sample_id": str})
    counts = counts[PAS_ATLAS_SAMPLE_COUNT_COLUMNS]

    matrix = counts.pivot_table(
        index="atlas_pas_id",
        columns="sample_id",
        values="atlas_total_umi_from_sample_pas",
        aggfunc="sum",
        fill_value=0,
        observed=False,
    )
    matrix = matrix.reindex(atlas_order, fill_value=0)
    matrix = matrix.reindex(columns=sample_order, fill_value=0)
    matrix = matrix.reset_index()
    return counts, matrix


def build_pas_atlas_quantification_windows(
    pas_atlas: pd.DataFrame,
    requant_config: Optional[dict] = None,
) -> pd.DataFrame:
    """Build conservative quantification windows for atlas PAS re-quantification."""
    config = requant_config or {}
    atlas = pas_atlas.copy()
    atlas = _ensure_columns(
        atlas,
        [
            "atlas_pas_id",
            "chrom",
            "strand",
            "gene_id",
            "representative_cut_site",
            "atlas_start",
            "atlas_end",
        ],
    )
    window_padding_bp = int(config.get("window_padding_bp", 100))
    window_source = str(config.get("window_source", "atlas_interval"))
    rows = []
    for row in atlas.itertuples(index=False):
        representative_cut_site = pd.to_numeric(getattr(row, "representative_cut_site", ""), errors="coerce")
        atlas_start = pd.to_numeric(getattr(row, "atlas_start", ""), errors="coerce")
        atlas_end = pd.to_numeric(getattr(row, "atlas_end", ""), errors="coerce")
        if window_source == "atlas_interval" and pd.notna(atlas_start) and pd.notna(atlas_end):
            start = int(atlas_start)
            end = int(atlas_end)
            source = "atlas_interval"
        elif pd.notna(representative_cut_site):
            center = int(representative_cut_site)
            start = max(0, center - window_padding_bp)
            end = center + window_padding_bp + 1
            source = "cut_site_padding"
        else:
            start = -1
            end = -1
            source = "invalid"

        rows.append(
            {
                "atlas_pas_id": getattr(row, "atlas_pas_id", ""),
                "chrom": getattr(row, "chrom", ""),
                "window_start": start,
                "window_end": end,
                "strand": getattr(row, "strand", ""),
                "gene_id": getattr(row, "gene_id", ""),
                "representative_cut_site": getattr(row, "representative_cut_site", ""),
                "window_source": source,
                "window_status": "usable" if start >= 0 and end > start else "invalid_window",
            }
        )
    windows = pd.DataFrame(rows, columns=ATLAS_REQUANT_WINDOW_COLUMNS)
    if bool(config.get("split_same_gene_overlaps", True)):
        windows = _split_same_gene_overlapping_quantification_windows(windows)
    return _mark_overlapping_quantification_windows(windows)


def build_pas_atlas_review_examples(pas_atlas: pd.DataFrame, max_per_category: int = 5) -> pd.DataFrame:
    """Return representative atlas rows for manual review."""
    if pas_atlas.empty:
        return pd.DataFrame(columns=["review_category"] + PAS_ATLAS_COLUMNS)

    examples = []
    seen_atlas_ids: set[str] = set()
    categories = [
        ("core", pas_atlas["atlas_status"] == "core"),
        ("private", pas_atlas["atlas_status"] == "private"),
        ("broad", pas_atlas["broad_or_unstable_cut_site"].astype(str).str.lower() == "true"),
        ("rescue", pas_atlas["atlas_confidence_tier"] == "rescue"),
    ]
    for category, mask in categories:
        unseen = ~pas_atlas["atlas_pas_id"].astype(str).isin(seen_atlas_ids)
        subset = pas_atlas[mask & unseen].head(max_per_category).copy()
        if subset.empty:
            continue
        seen_atlas_ids.update(subset["atlas_pas_id"].astype(str))
        subset.insert(0, "review_category", category)
        examples.append(subset)
    if not examples:
        return pd.DataFrame(columns=["review_category"] + PAS_ATLAS_COLUMNS)
    return pd.concat(examples, ignore_index=True)


def _load_pas_atlas_table(pas_atlas: pd.DataFrame | str | Path) -> pd.DataFrame:
    if isinstance(pas_atlas, pd.DataFrame):
        return pas_atlas.copy()
    return pd.read_csv(pas_atlas, sep="\t", dtype=str).fillna("")


def _normalize_requant_sample_inputs(
    bam_paths: Sequence[str] | dict[str, str],
    sample_ids: Optional[Sequence[str]],
    valid_barcode_paths: Optional[Sequence[Optional[str]] | dict[str, Optional[str]]],
) -> list[tuple[str, str, Optional[str]]]:
    if isinstance(bam_paths, dict):
        ordered_sample_ids = list(bam_paths.keys()) if sample_ids is None else list(sample_ids)
        barcode_lookup = valid_barcode_paths if isinstance(valid_barcode_paths, dict) else {}
        return [
            (
                sample_id,
                str(bam_paths.get(sample_id, "")),
                None if not isinstance(barcode_lookup, dict) else barcode_lookup.get(sample_id),
            )
            for sample_id in ordered_sample_ids
        ]

    bam_list = list(bam_paths)
    if sample_ids is None:
        raise ValueError("sample_ids are required when bam_paths is a sequence")
    if len(sample_ids) != len(bam_list):
        raise ValueError("sample_ids must have the same length as bam_paths")
    if valid_barcode_paths is None:
        barcode_list = [None] * len(bam_list)
    else:
        barcode_list = list(valid_barcode_paths)
        if len(barcode_list) != len(bam_list):
            raise ValueError("valid_barcode_paths must have the same length as bam_paths")
    return [(str(sample_id), str(bam_path), barcode_path) for sample_id, bam_path, barcode_path in zip(sample_ids, bam_list, barcode_list)]


def _atlas_windows_to_peak_bed(windows: pd.DataFrame) -> pd.DataFrame:
    if windows.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "name", "score", "strand"])
    bed = windows.copy()
    bed["start"] = pd.to_numeric(bed["window_start"], errors="coerce").fillna(-1).astype(int)
    bed["end"] = pd.to_numeric(bed["window_end"], errors="coerce").fillna(-1).astype(int)
    bed["score"] = pd.to_numeric(bed["representative_cut_site"], errors="coerce").fillna(0).astype(int)
    bed = bed.rename(columns={"atlas_pas_id": "name", "window_start": "_window_start", "window_end": "_window_end"})
    return bed[["chrom", "start", "end", "name", "score", "strand"]].copy()


def _split_same_gene_overlapping_quantification_windows(windows: pd.DataFrame) -> pd.DataFrame:
    """Partition overlapping windows for distinct PAS of the same gene by cut-site midpoint."""
    if windows.empty:
        return windows

    windows = windows.copy()
    usable = windows["window_status"].eq("usable")
    same_gene = usable & windows["gene_id"].astype(str).str.len().gt(0)
    changed_indices = set()
    for _, group in windows[same_gene].groupby(["chrom", "strand", "gene_id"], dropna=False, sort=False):
        group = group.copy()
        group["cut_num"] = pd.to_numeric(group["representative_cut_site"], errors="coerce")
        group["start_num"] = pd.to_numeric(group["window_start"], errors="coerce")
        group["end_num"] = pd.to_numeric(group["window_end"], errors="coerce")
        group = group.dropna(subset=["cut_num", "start_num", "end_num"])
        group = group.sort_values(["cut_num", "start_num", "end_num", "atlas_pas_id"])
        indices = list(group.index)

        for left_index, right_index in zip(indices, indices[1:]):
            left_cut = int(group.at[left_index, "cut_num"])
            right_cut = int(group.at[right_index, "cut_num"])
            if right_cut <= left_cut:
                continue

            left_start = int(windows.at[left_index, "window_start"])
            left_end = int(windows.at[left_index, "window_end"])
            right_start = int(windows.at[right_index, "window_start"])
            right_end = int(windows.at[right_index, "window_end"])
            if right_start >= left_end:
                continue

            boundary = (left_cut + right_cut) // 2
            left_end_new = min(left_end, boundary)
            right_start_new = max(right_start, boundary)
            if left_end_new <= left_start or right_end <= right_start_new:
                continue

            if left_end_new != left_end:
                windows.at[left_index, "window_end"] = left_end_new
                changed_indices.add(left_index)
            if right_start_new != right_start:
                windows.at[right_index, "window_start"] = right_start_new
                changed_indices.add(right_index)

    for index in changed_indices:
        source = str(windows.at[index, "window_source"])
        if not source.endswith(":same_gene_midpoint"):
            windows.at[index, "window_source"] = f"{source}:same_gene_midpoint"
    return windows


def _mark_overlapping_quantification_windows(windows: pd.DataFrame) -> pd.DataFrame:
    if windows.empty:
        return windows
    windows = windows.copy()
    usable = windows["window_status"] == "usable"
    overlap_indices: set[int] = set()
    for _, group in windows[usable].groupby(["chrom", "strand"], dropna=False, sort=False):
        group = group.copy()
        group["window_start_num"] = pd.to_numeric(group["window_start"], errors="coerce").fillna(-1).astype(int)
        group["window_end_num"] = pd.to_numeric(group["window_end"], errors="coerce").fillna(-1).astype(int)
        group = group.sort_values(["window_start_num", "window_end_num", "atlas_pas_id"])

        max_end = -1
        max_end_index: Optional[int] = None
        for row_index, row in group.iterrows():
            start = int(row["window_start_num"])
            end = int(row["window_end_num"])
            if max_end_index is not None and start < max_end:
                overlap_indices.add(int(max_end_index))
                overlap_indices.add(int(row_index))
            if end > max_end:
                max_end = end
                max_end_index = int(row_index)

    mask = windows.index.isin(overlap_indices) & usable
    windows.loc[mask, "window_status"] = "window_overlap_ambiguous"
    return windows


def _nonusable_requant_rows(atlas_windows: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    rows = _placeholder_requant_rows(atlas_windows, sample_id)
    return rows[rows["window_status"] != "usable"].copy()


def _placeholder_requant_rows(
    atlas_windows: pd.DataFrame,
    sample_id: str,
    default_status: str = "window_overlap_ambiguous",
) -> pd.DataFrame:
    rows = []
    for row in atlas_windows.itertuples(index=False):
        status = default_status if row.window_status == "usable" else row.window_status
        rows.append(
            {
                "sample_id": sample_id,
                "atlas_pas_id": row.atlas_pas_id,
                "chrom": row.chrom,
                "window_start": row.window_start,
                "window_end": row.window_end,
                "strand": row.strand,
                "window_status": row.window_status,
                "observed_umi_pre_cell_filter": 0,
                "observed_cells_pre_cell_filter": 0,
                "observed_umi_post_cell_filter": 0,
                "observed_cells_post_cell_filter": 0,
                "has_umi": "false",
                "in_matrix": "false",
                "drop_reason": status,
                "quantification_status": status,
            }
        )
    return pd.DataFrame(rows, columns=ATLAS_REQUANT_STATUS_COLUMNS)


def _read_valid_barcodes(path: Optional[str]) -> Optional[set[str]]:
    if not path or not os.path.exists(path):
        return None
    with open(path, "r", encoding="utf-8") as handle:
        return {line.strip() for line in handle if line.strip()}


def _read_requant_status(sample_dir: Path, quant_windows: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    status_path = sample_dir / "pas_quantification_status.tsv"
    status = pd.read_csv(status_path, sep="\t", dtype=str).fillna("")
    status = status.rename(
        columns={
            "name": "atlas_pas_id",
            "start": "window_start",
            "end": "window_end",
            "umi_count_pre_cell_filter": "observed_umi_pre_cell_filter",
            "cell_count_pre_cell_filter": "observed_cells_pre_cell_filter",
            "umi_count_post_cell_filter": "observed_umi_post_cell_filter",
            "cell_count_post_cell_filter": "observed_cells_post_cell_filter",
        }
    )
    status["window_start"] = _numeric(status["window_start"]).astype(int)
    status["window_end"] = _numeric(status["window_end"]).astype(int)
    quant_windows = quant_windows.copy()
    quant_windows["window_start"] = _numeric(quant_windows["window_start"]).astype(int)
    quant_windows["window_end"] = _numeric(quant_windows["window_end"]).astype(int)
    merged = quant_windows.merge(status, on=["atlas_pas_id", "chrom", "window_start", "window_end", "strand"], how="left")
    merged["sample_id"] = sample_id
    merged["window_status"] = merged["window_status"].fillna("usable")
    merged["observed_umi_pre_cell_filter"] = _numeric(merged["observed_umi_pre_cell_filter"]).astype(int)
    merged["observed_cells_pre_cell_filter"] = _numeric(merged["observed_cells_pre_cell_filter"]).astype(int)
    merged["observed_umi_post_cell_filter"] = _numeric(merged["observed_umi_post_cell_filter"]).astype(int)
    merged["observed_cells_post_cell_filter"] = _numeric(merged["observed_cells_post_cell_filter"]).astype(int)
    merged["has_umi"] = merged["has_umi"].replace("", "false")
    merged["in_matrix"] = merged["in_matrix"].replace("", "false")
    merged["drop_reason"] = merged["drop_reason"].replace("", "not_quantified_error")
    merged["quantification_status"] = merged.apply(_derive_requant_status, axis=1)
    return merged[ATLAS_REQUANT_STATUS_COLUMNS].copy()


def _derive_requant_status(row: pd.Series) -> str:
    if row["window_status"] != "usable":
        return str(row["window_status"])
    drop_reason = str(row.get("drop_reason", ""))
    umi_post = _to_int(row.get("observed_umi_post_cell_filter", 0))
    umi_pre = _to_int(row.get("observed_umi_pre_cell_filter", 0))
    if drop_reason == "none" and umi_post > 0:
        return "quantified_detected"
    if drop_reason == "only_filtered_cells" and umi_pre > 0 and umi_post == 0:
        return "quantified_only_filtered_cells"
    if drop_reason == "no_overlapping_umi" or (umi_pre == 0 and umi_post == 0):
        return "quantified_zero"
    return "not_quantified_error"


def _summarize_requant_status(status: pd.DataFrame) -> dict[str, int]:
    qc = {}
    value_counts = status["quantification_status"].astype(str).value_counts()
    for key in [
        "quantified_detected",
        "quantified_zero",
        "quantified_only_filtered_cells",
        "not_quantified_missing_bam",
        "not_quantified_missing_barcode",
        "not_quantified_error",
    ]:
        qc[f"n_{key}"] = int(value_counts.get(key, 0))
    return qc


def _build_requant_matrix(
    quantification_status: pd.DataFrame,
    value_column: str,
    sample_specs: Sequence[tuple[str, str, Optional[str]]],
    atlas_windows: pd.DataFrame,
) -> pd.DataFrame:
    atlas_order = atlas_windows["atlas_pas_id"].astype(str).tolist()
    sample_order = [sample_id for sample_id, _, _ in sample_specs]
    if quantification_status.empty:
        matrix = pd.DataFrame({"atlas_pas_id": atlas_order})
        for sample_id in sample_order:
            matrix[sample_id] = 0
        return matrix
    status = quantification_status.copy()
    status[value_column] = _numeric(status[value_column]).astype(int)
    matrix = status.pivot_table(
        index="atlas_pas_id",
        columns="sample_id",
        values=value_column,
        aggfunc="sum",
        fill_value=0,
        observed=False,
    )
    matrix = matrix.reindex(atlas_order, fill_value=0)
    matrix = matrix.reindex(columns=sample_order, fill_value=0)
    return matrix.reset_index()


def _to_int(value) -> int:
    numeric = pd.to_numeric(value, errors="coerce")
    if pd.isna(numeric):
        return 0
    return int(numeric)


def _cluster_group(group: pd.DataFrame, config: dict) -> list[pd.DataFrame]:
    clusters = []
    current_rows = []
    for _, row in group.iterrows():
        if not current_rows:
            current_rows = [row]
            continue
        current = pd.DataFrame(current_rows)
        candidate = pd.concat([current, row.to_frame().T], ignore_index=True)
        representative = _representative_cut_site(current)
        window = _cluster_window(row["evidence_tier_for_atlas"], config)
        max_width = _max_cluster_width(row["evidence_tier_for_atlas"], config)
        distance = abs(int(row["cut_site_num"]) - representative)
        width = int(candidate["cut_site_num"].max()) - int(candidate["cut_site_num"].min())
        if distance <= window and width <= max_width:
            current_rows.append(row)
        else:
            clusters.append(current)
            current_rows = [row]
    if current_rows:
        clusters.append(pd.DataFrame(current_rows))
    return clusters


def _atlas_row(cluster: pd.DataFrame, config: dict) -> dict:
    representative = _representative_cut_site(cluster)
    cut_sites = cluster["cut_site_num"].astype(int)
    cut_min = int(cut_sites.min())
    cut_max = int(cut_sites.max())
    cut_mad = float(np.median(np.abs(cut_sites - np.median(cut_sites)))) if len(cut_sites) else 0.0
    total_umi = int(_numeric(cluster["total_pas_umi"]).sum())
    total_cells = int(_numeric(cluster["n_cells_detected"]).sum())
    support_library_count = int(cluster["sample_id"].nunique())
    confidence_tier = _atlas_confidence_tier(cluster, config)
    max_width = _max_cluster_width(confidence_tier, config)
    broad = str((cut_max - cut_min) > max_width).lower()
    peak_starts = pd.to_numeric(cluster["start"], errors="coerce")
    peak_ends = pd.to_numeric(cluster["end"], errors="coerce")
    valid_peaks = peak_starts.notna() & peak_ends.notna() & peak_ends.gt(peak_starts)
    if valid_peaks.any():
        atlas_start = max(0, int(peak_starts[valid_peaks].min()))
        atlas_end = int(peak_ends[valid_peaks].max())
    else:
        padding = int(config.get("quantification_padding_bp", 100))
        atlas_start = max(0, cut_min - padding)
        atlas_end = cut_max + padding + 1
    atlas_status = _atlas_status(support_library_count, config)
    gene_id = _first_nonempty(cluster["primary_gene_id"])
    raw_id = "|".join(
        [
            _first_nonempty(cluster["species"]),
            _first_nonempty(cluster["genome_build"]),
            _first_nonempty(cluster["chrom"]),
            _first_nonempty(cluster["strand"]),
            gene_id,
            str(representative),
            str(atlas_start),
            str(atlas_end),
            confidence_tier,
        ]
    )
    short_hash = hashlib.sha1(raw_id.encode("utf-8")).hexdigest()[:10]
    atlas_pas_id = (
        f"{_first_nonempty(cluster['species'])}:"
        f"{_first_nonempty(cluster['genome_build'])}:"
        f"{_first_nonempty(cluster['chrom'])}:"
        f"{_first_nonempty(cluster['strand'])}:"
        f"{representative}:atlas:{short_hash}"
    )
    return {
        "atlas_pas_id": atlas_pas_id,
        "species": _first_nonempty(cluster["species"]),
        "genome_build": _first_nonempty(cluster["genome_build"]),
        "chrom": _first_nonempty(cluster["chrom"]),
        "strand": _first_nonempty(cluster["strand"]),
        "representative_cut_site": representative,
        "cut_site_min": cut_min,
        "cut_site_max": cut_max,
        "cut_site_mad": _format_float(cut_mad),
        "atlas_start": atlas_start,
        "atlas_end": atlas_end,
        "gene_id": gene_id,
        "gene_symbol": _first_nonempty(cluster["primary_gene_symbol"]),
        "atlas_region_class": _dominant_value(cluster["primary_region_class"]),
        "atlas_confidence_tier": confidence_tier,
        "support_library_count": support_library_count,
        "support_sample_count": support_library_count,
        "member_peak_count": len(cluster),
        "total_member_umi": total_umi,
        "total_member_cells": total_cells,
        "has_motif_support": str(cluster["best_pas_motif"].astype(str).str.len().gt(0).any()).lower(),
        "has_expected_motif_support": str(
            (
                cluster["best_pas_motif"].astype(str).str.len().gt(0)
                & cluster["best_pas_motif_search_class"].astype(str).isin(["", "expected_upstream"])
            ).any()
        ).lower(),
        "has_expanded_motif_support": str(
            (
                cluster["best_pas_motif"].astype(str).str.len().gt(0)
                & cluster["best_pas_motif_search_class"].astype(str).eq("expanded_proxy_window")
            ).any()
        ).lower(),
        "has_external_support": str(_has_external_support(cluster["external_overlap_status"])).lower(),
        "broad_or_unstable_cut_site": broad,
        "atlas_status": atlas_status,
        "has_same_strand_other_gene_intron_overlap": str(
            cluster["same_strand_other_gene_intron_overlap"].astype(str).eq("true").any()
        ).lower(),
        "has_same_strand_exonic_gene_conflict": str(
            cluster["same_strand_exonic_gene_conflict"].astype(str).eq("true").any()
        ).lower(),
        "has_opposite_strand_gene_overlap": str(
            cluster["opposite_strand_gene_overlap"].astype(str).eq("true").any()
        ).lower(),
        "neighbor_gene_decision_rules": _joined_values(
            cluster["neighbor_gene_decision_rule"]
        ),
        "neighbor_gene_terminal_states": _joined_values(
            cluster["neighbor_gene_terminal_state"]
        ),
        "has_transcript_model_conflict": str(
            cluster["transcript_model_conflict"].astype(str).eq("true").any()
        ).lower(),
        "transcript_context_classes": _joined_values(
            cluster["transcript_context_class"]
        ),
        "has_biotype_assignment_conflict": str(
            cluster["biotype_assignment_conflict"].astype(str).eq("true").any()
        ).lower(),
        "biotype_tie_break_statuses": _joined_values(
            cluster["biotype_tie_break_status"]
        ),
    }


def _derive_evidence_tier(row) -> tuple[str, str]:
    region = str(row.primary_region_class)
    status = str(row.primary_assignment_status)
    confidence = str(row.confidence_tier)
    gene_id = str(row.primary_gene_id)
    sequence_status = str(row.sequence_qc_status)
    internal_priming = _is_true(row.internal_priming_flag)

    if not gene_id:
        return "exclude", "excluded_no_gene"
    if status == "ambiguous_multi_gene" or region in {"multi_gene_ambiguous_pas"}:
        return "exclude", "ambiguous_gene_context"
    if region in EXCLUDED_REGION_CLASSES:
        return "exclude", "excluded"
    if internal_priming:
        return "sensitivity", "retained_internal_priming_risk"
    if confidence == "exclude":
        return "exclude", "excluded"
    if confidence == "relaxed" or status == "assigned_sensitivity":
        return "sensitivity", "usable_sensitivity"
    if region in STRICT_REGION_CLASSES and status == "assigned_primary":
        if sequence_status != "computed":
            return "rescue", "usable_sequence_qc_missing"
        has_motif = bool(str(row.best_pas_motif).strip())
        motif_class = str(row.best_pas_motif_search_class).strip()
        has_expected_motif = has_motif and motif_class in {"", "expected_upstream"}
        has_external = _has_external_support(pd.Series([row.external_overlap_status]))
        if has_expected_motif or has_external:
            return "strict", "usable_positive_cleavage_support"
        if has_motif:
            return "rescue", "usable_expanded_motif_support"
        return "rescue", "usable_no_positive_cleavage_support"
    if region in RESCUE_REGION_CLASSES and status == "assigned_primary":
        return "rescue", "usable"
    if status == "assigned_primary":
        return "rescue", "usable_primary_fallback"
    return "exclude", "excluded"


def _main_atlas_eligibility(row) -> tuple[bool, str]:
    """Validate and apply the shared, fail-closed neighboring-gene gate."""
    return validate_neighbor_policy_row(row)


def _tez_sensitivity_rescue_policy_status(row, config: dict) -> str:
    """Return a structured, fail-closed TES2kb promotion decision."""
    if not bool(config.get("promote_tez_sensitivity_to_rescue", False)):
        return "policy_disabled"

    region = str(row.primary_region_class).strip()
    status = str(row.primary_assignment_status).strip()
    if region != "transcript_end_zone_pas":
        return "not_applicable_region_class"
    if status != "assigned_sensitivity":
        return "not_applicable_assignment_status"

    confidence = str(row.confidence_tier).strip()
    if not confidence:
        return "not_promoted_missing_confidence"
    if confidence == "exclude":
        return "not_promoted_excluded_confidence"
    if not str(row.primary_gene_id).strip():
        return "not_promoted_missing_primary_gene"
    ambiguous = str(getattr(row, "multi_gene_ambiguous", "")).strip().lower()
    if not ambiguous:
        return "not_promoted_missing_ambiguity_evidence"
    if ambiguous == "true":
        return "not_promoted_ambiguous_gene_context"
    if ambiguous != "false":
        return "not_promoted_invalid_ambiguity_evidence"
    if not str(row.primary_transcript_id).strip():
        return "not_promoted_missing_primary_transcript"
    if not str(row.primary_candidate_id).strip():
        return "not_promoted_missing_primary_candidate"

    required_text = [
        "candidate_gene_ids",
        "candidate_region_classes",
        "original_tes_candidate_ids",
        "selected_tes_candidate_id",
        "selected_tes_policy_class",
        "selected_tes_policy_reason",
        "neighbor_gene_decision_rule",
        "neighbor_gene_terminal_state",
    ]
    if str(getattr(row, "selected_tes_policy_class", "")).strip() == "qualified_tes2kb_rescue":
        required_text.extend(
            [
                "qualified_same_strand_tes_candidate_ids",
                "qualified_same_strand_tes_gene_ids",
                "qualified_same_strand_tes_transcript_ids",
            ]
        )
    for field in required_text:
        if not str(getattr(row, field, "")).strip():
            return f"not_promoted_missing_{field}"
    if str(row.selected_tes_candidate_id).strip() != str(row.primary_candidate_id).strip():
        return "not_promoted_inconsistent_selected_tes_candidate"

    for field in [
        "same_strand_other_gene_intron_overlap",
        "same_strand_exonic_gene_conflict",
        "opposite_strand_gene_overlap",
    ]:
        value = str(getattr(row, field, "")).strip().lower()
        if not value:
            return f"not_promoted_missing_{field}"
        if value not in {"true", "false"}:
            return f"not_promoted_invalid_{field}"
    tes_candidate_count = pd.to_numeric(
        getattr(row, "qualified_same_strand_tes_candidate_count", ""),
        errors="coerce",
    )
    tes_gene_count = pd.to_numeric(
        getattr(row, "qualified_same_strand_tes_gene_count", ""),
        errors="coerce",
    )
    if pd.isna(tes_candidate_count) or pd.isna(tes_gene_count):
        return "not_promoted_missing_tes_competition_count"
    selected_policy_class = str(
        getattr(row, "selected_tes_policy_class", "")
    ).strip()
    if selected_policy_class == "qualified_tes2kb_rescue" and (
        float(tes_candidate_count) < 1 or float(tes_gene_count) < 1
    ):
        return "not_promoted_invalid_tes_competition_count"

    transcript_end = pd.to_numeric(
        getattr(row, "primary_transcript_end", ""),
        errors="coerce",
    )
    absolute_distance = pd.to_numeric(
        getattr(row, "distance_to_primary_transcript_end", ""),
        errors="coerce",
    )
    signed_distance = pd.to_numeric(
        getattr(row, "signed_distance_to_primary_transcript_end", ""),
        errors="coerce",
    )
    downstream_distance = pd.to_numeric(
        getattr(row, "downstream_distance_to_primary_transcript_end", ""),
        errors="coerce",
    )
    if pd.isna(transcript_end):
        return "not_promoted_missing_primary_tes"
    if pd.isna(absolute_distance) or pd.isna(signed_distance) or pd.isna(downstream_distance):
        return "not_promoted_missing_tes_distance"
    if (
        float(absolute_distance) < 0
        or float(absolute_distance) != abs(float(signed_distance))
        or float(absolute_distance) != abs(float(downstream_distance))
    ):
        return "not_promoted_inconsistent_tes_distance"

    strand_relation = str(row.primary_strand_relation).strip()
    transcript_strand = str(row.primary_transcript_strand).strip()
    pas_strand = str(row.strand).strip()
    if not strand_relation or transcript_strand not in {"+", "-"} or pas_strand not in {"+", "-"}:
        return "not_promoted_missing_strand_evidence"
    if strand_relation != "same_strand" or transcript_strand != pas_strand:
        return "not_promoted_not_same_strand"
    expected_directed = float(signed_distance) if transcript_strand == "+" else -float(signed_distance)
    if expected_directed != float(downstream_distance):
        return "not_promoted_inconsistent_strand_directed_distance"

    internal_priming = str(getattr(row, "internal_priming_flag", "")).strip().lower()
    if not internal_priming:
        return "not_promoted_missing_internal_priming_evidence"
    if internal_priming == "true":
        return "not_promoted_internal_priming"
    if internal_priming != "false":
        return "not_promoted_invalid_internal_priming_evidence"

    min_distance = int(config["tez_rescue_min_downstream_distance_bp"])
    max_distance = int(config["tez_rescue_max_downstream_distance_bp"])
    downstream_distance = float(downstream_distance)
    if downstream_distance <= min_distance:
        return "not_promoted_at_or_below_primary_boundary"
    if downstream_distance > max_distance:
        return "not_promoted_beyond_max_distance"
    if str(row.same_strand_exonic_gene_conflict).strip().lower() == "true":
        return "not_promoted_same_strand_exonic_gene_conflict"
    return "promoted_tez_sensitivity_to_rescue"


def _evidence_group(evidence_tier: str, region: str) -> str:
    if evidence_tier in {"strict", "rescue"} and region in STRICT_REGION_CLASSES:
        return "utr_or_terminal_exon"
    if evidence_tier == "rescue" and region == "transcript_end_zone_pas":
        return "near_tez_rescue"
    if evidence_tier == "rescue":
        return "structural_rescue"
    if evidence_tier == "strict":
        return "utr_or_terminal_exon"
    if evidence_tier == "sensitivity":
        return "sensitivity_or_other"
    return "excluded"


def _cluster_window(evidence_tier: str, config: dict) -> int:
    if evidence_tier == "strict":
        return int(config.get("strict_cluster_window_bp", 25))
    if evidence_tier == "sensitivity":
        return int(config.get("sensitivity_cluster_window_bp", 50))
    return int(config.get("rescue_cluster_window_bp", 50))


def _max_cluster_width(evidence_tier: str, config: dict) -> int:
    if evidence_tier == "strict":
        return int(config.get("max_cluster_width_strict_bp", 50))
    if evidence_tier == "sensitivity":
        return int(config.get("max_cluster_width_sensitivity_bp", 100))
    return int(config.get("max_cluster_width_rescue_bp", 100))


def _representative_cut_site(group: pd.DataFrame) -> int:
    values = group["cut_site_num"].astype(int).to_numpy()
    weights = _numeric(group["total_pas_umi"]).to_numpy()
    return int(_weighted_median(values, weights))


def _weighted_median(values: np.ndarray, weights: np.ndarray) -> int:
    if len(values) == 0:
        return 0
    if len(values) == 1:
        return int(values[0])
    weights = np.nan_to_num(weights.astype(float), nan=0.0)
    if weights.sum() <= 0:
        return int(np.median(values))
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    cutoff = weights.sum() / 2.0
    return int(values[np.searchsorted(np.cumsum(weights), cutoff, side="left")])


def _atlas_confidence_tier(cluster: pd.DataFrame, config: dict) -> str:
    tiers = set(cluster["evidence_tier_for_atlas"].astype(str))
    if "strict" in tiers:
        return "strict"
    if "rescue" in tiers:
        min_reproducible = int(config.get("min_support_libraries_strict", 2))
        reproducible = cluster["sample_id"].nunique() >= min_reproducible
        strict_region_only = cluster["primary_region_class"].astype(str).isin(STRICT_REGION_CLASSES).all()
        sequence_qc_complete = cluster["sequence_qc_status"].astype(str).eq("computed").all()
        if reproducible and strict_region_only and sequence_qc_complete:
            return "strict"
        return "rescue"
    if "sensitivity" in tiers:
        return "sensitivity"
    return "exclude"


def _atlas_status(support_library_count: int, config: dict) -> str:
    min_core = int(config.get("min_support_libraries_core", 2))
    if support_library_count >= min_core:
        return "core"
    if bool(config.get("allow_private_pas", True)):
        return "private"
    return "excluded"


def _mapping_status(distance: int, broad: str) -> str:
    if str(broad).lower() == "true":
        return "broad_cluster"
    if distance == 0:
        return "exact"
    return "nearby"


def _unmapped_mapping_status(row, include_sensitivity: bool) -> str:
    if str(getattr(row, "main_atlas_eligible", "")).lower() == "false":
        reason = str(getattr(row, "main_atlas_exclusion_reason", ""))
        if reason == "same_strand_exonic_gene_conflict":
            return "excluded_sensitivity_same_strand_exonic_gene_conflict"
        if reason:
            return f"excluded_{reason}"
    if str(row.evidence_tier_for_atlas) == "sensitivity" and not include_sensitivity:
        return "not_clustered_sensitivity"
    return row.candidate_status or "no_match"


def _sample_id_from_annotation_dir(annotation_dir: str) -> str:
    path = Path(annotation_dir)
    if path.name == "pas_annotation" and path.parent.name:
        return path.parent.name
    return path.name


def _read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def _ensure_columns(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    df = df.copy()
    for col in columns:
        if col not in df.columns:
            df[col] = ""
    return df


def _numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").fillna(0)


def _is_true(value) -> bool:
    return str(value).strip().lower() == "true"


def _bool_string(value) -> str:
    return str(_is_true(value)).lower()


def _first_nonempty(values) -> str:
    for value in values:
        if isinstance(value, str) and value:
            return value
        if pd.notna(value) and str(value):
            return str(value)
    return ""


def _dominant_value(values: pd.Series) -> str:
    cleaned = values.astype(str)
    cleaned = cleaned[cleaned.str.len() > 0]
    if cleaned.empty:
        return ""
    counts = cleaned.value_counts()
    top = counts[counts == counts.max()].index.tolist()
    if len(top) == 1:
        return top[0]
    return "mixed:" + ",".join(sorted(top))


def _joined_values(values: pd.Series) -> str:
    return "|".join(sorted(set(value for value in values.astype(str) if value)))


def _external_status_summary(values: pd.Series) -> str:
    statuses = [str(value) for value in values if str(value)]
    if not statuses:
        return ""
    informative = [value for value in statuses if value not in {"unavailable", "not_run"}]
    return ";".join(sorted(set(informative or statuses)))


def _has_external_support(values: pd.Series) -> bool:
    positive = {"supported", "overlap", "present", "exact", "nearby"}
    return any(
        token in positive
        for value in values
        for token in str(value).lower().split(";")
    )


def _bool_count(series: pd.Series) -> int:
    if series.empty:
        return 0
    return int(series.astype(str).str.lower().eq("true").sum())


def _format_float(value: float) -> str:
    return f"{value:.6g}"
