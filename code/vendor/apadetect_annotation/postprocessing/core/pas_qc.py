"""Productized PAS annotation and atlas QC tables, reports, and figures."""

from __future__ import annotations

import platform
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


@dataclass
class PasQcResult:
    output_dir: str
    summary: pd.DataFrame
    report_path: str
    figure_status: pd.DataFrame


def write_pas_annotation_qc(
    annotation_dir: str,
    output_dir: Optional[str] = None,
    write_plots: bool = True,
) -> PasQcResult:
    """Write inspectable annotation attrition, evidence, and review artifacts."""
    annotation_path = Path(annotation_dir)
    qc_dir = Path(output_dir) if output_dir else annotation_path / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)

    feature = _read_tsv(annotation_path / "pas_feature_annotation.tsv")
    candidates = _read_tsv(annotation_path / "pas_annotation_candidates.tsv")
    confidence = _read_tsv(annotation_path / "pas_confidence_tiers.tsv")
    sequence = _read_tsv(annotation_path / "pas_sequence_qc.tsv")
    external = _read_tsv(annotation_path / "pas_external_overlap.tsv")
    testability = _read_tsv(annotation_path / "gene_testability.tsv")

    summary = _annotation_summary(feature, confidence, sequence, external, testability)
    region_counts = _count_table(feature, ["species", "genome_build", "primary_region_class"], "n_pas")
    assignment_counts = _count_table(feature, ["primary_assignment_status"], "n_pas")
    confidence_counts = _count_table(confidence, ["confidence_tier"], "n_pas")
    external_counts = _count_table(external, ["source_name", "overlap_status", "distance_bin"], "n_pas_source_rows")
    motif_counts = _motif_internal_counts(sequence, confidence)
    flow = _candidate_primary_flow(candidates, feature)
    funnel = _gene_testability_funnel(feature, confidence, testability)
    examples = _representative_pas_examples(feature, confidence, sequence, external)

    tables = {
        "pas_annotation_qc_summary.tsv": summary,
        "pas_region_class_counts.tsv": region_counts,
        "pas_assignment_status_counts.tsv": assignment_counts,
        "pas_confidence_tier_counts.tsv": confidence_counts,
        "pas_external_support_counts.tsv": external_counts,
        "pas_motif_internal_priming_counts.tsv": motif_counts,
        "pas_candidate_to_primary_flow.tsv": flow,
        "gene_testability_funnel.tsv": funnel,
        "representative_pas_examples.tsv": examples,
    }
    for filename, table in tables.items():
        table.to_csv(qc_dir / filename, sep="\t", index=False)

    figure_status = _write_annotation_figures(
        qc_dir,
        region_counts,
        flow,
        external_counts,
        motif_counts,
        funnel,
        examples,
        enabled=write_plots,
    )
    figure_status.to_csv(qc_dir / "pas_qc_figure_status.tsv", sep="\t", index=False)
    report_path = qc_dir / "pas_annotation_qc_report.md"
    report_path.write_text(
        _annotation_report(annotation_path, summary, tables, figure_status),
        encoding="utf-8",
    )
    return PasQcResult(str(qc_dir), summary, str(report_path), figure_status)


def write_pas_atlas_qc(
    atlas_dir: str,
    requantification_dir: Optional[str] = None,
    output_dir: Optional[str] = None,
    write_plots: bool = True,
) -> PasQcResult:
    """Write species-atlas attrition, mapping, stability, and re-quant QC."""
    atlas_path = Path(atlas_dir)
    qc_dir = Path(output_dir) if output_dir else atlas_path / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)

    evidence = _read_tsv(atlas_path / "pas_sample_evidence.tsv")
    atlas = _read_tsv(atlas_path / "pas_atlas.tsv")
    mapping = _read_tsv(atlas_path / "pas_sample_to_atlas_map.tsv")
    requant_path = Path(requantification_dir) if requantification_dir else None
    requant_status = _read_tsv(requant_path / "pas_atlas_quantification_status.tsv") if requant_path else pd.DataFrame()

    summary = _atlas_summary(evidence, atlas, mapping, requant_status)
    evidence_counts = _count_table(evidence, ["species", "evidence_tier_for_atlas"], "n_sample_pas")
    atlas_counts = _count_table(atlas, ["species", "atlas_confidence_tier", "atlas_status"], "n_atlas_pas")
    mapping_counts = _count_table(mapping, ["mapping_status"], "n_sample_pas")
    support_counts = _count_table(atlas, ["support_library_count"], "n_atlas_pas")
    adjacent = _adjacent_atlas_distances(atlas)
    requant_counts = _count_table(requant_status, ["sample_id", "quantification_status"], "n_atlas_sample_rows")

    tables = {
        "pas_atlas_product_qc_summary.tsv": summary,
        "pas_atlas_evidence_tier_counts.tsv": evidence_counts,
        "pas_atlas_status_counts.tsv": atlas_counts,
        "pas_atlas_mapping_status_counts.tsv": mapping_counts,
        "pas_atlas_support_library_distribution.tsv": support_counts,
        "pas_atlas_adjacent_site_distances.tsv": adjacent,
        "pas_atlas_requant_status_counts.tsv": requant_counts,
    }
    for filename, table in tables.items():
        table.to_csv(qc_dir / filename, sep="\t", index=False)

    figure_status = _write_atlas_figures(
        qc_dir,
        evidence_counts,
        atlas_counts,
        mapping_counts,
        support_counts,
        adjacent,
        requant_counts,
        enabled=write_plots,
    )
    figure_status.to_csv(qc_dir / "pas_atlas_qc_figure_status.tsv", sep="\t", index=False)
    report_path = qc_dir / "pas_atlas_product_qc_report.md"
    report_path.write_text(
        _atlas_report(atlas_path, requant_path, summary, tables, figure_status),
        encoding="utf-8",
    )
    return PasQcResult(str(qc_dir), summary, str(report_path), figure_status)


def _annotation_summary(feature, confidence, sequence, external, testability) -> pd.DataFrame:
    metrics = [
        ("n_pas", len(feature)),
        ("n_assigned_primary", _eq_count(feature, "primary_assignment_status", "assigned_primary")),
        ("n_assigned_sensitivity", _eq_count(feature, "primary_assignment_status", "assigned_sensitivity")),
        ("n_gene_level_exclude", _eq_count(confidence, "confidence_tier", "exclude")),
        ("n_sequence_qc_computed", _eq_count(sequence, "sequence_qc_status", "computed")),
        ("n_expected_upstream_motif", _eq_count(sequence, "best_pas_motif_search_class", "expected_upstream")),
        ("n_expanded_proxy_window_motif", _eq_count(sequence, "best_pas_motif_search_class", "expanded_proxy_window")),
        ("n_internal_priming_risk", _true_count(sequence, "internal_priming_flag")),
        ("n_external_supported_source_rows", _eq_count(external, "overlap_status", "supported")),
        ("n_external_absent_source_rows", _eq_count(external, "overlap_status", "absent")),
        ("n_external_unavailable_source_rows", _eq_count(external, "overlap_status", "unavailable")),
        ("n_genes", len(testability)),
        ("n_apa_testable_primary_genes", _true_count(testability, "is_apa_testable_primary")),
        ("n_apa_testable_relaxed_genes", _true_count(testability, "is_apa_testable_relaxed")),
        ("n_pas_in_matrix", _true_count(feature, "in_matrix")),
    ]
    return pd.DataFrame(metrics, columns=["metric", "value"])


def _atlas_summary(evidence, atlas, mapping, requant_status) -> pd.DataFrame:
    metrics = [
        ("n_sample_pas_evidence", len(evidence)),
        ("n_samples", evidence["sample_id"].nunique() if "sample_id" in evidence else 0),
        ("n_atlas_pas", len(atlas)),
        ("n_atlas_strict", _eq_count(atlas, "atlas_confidence_tier", "strict")),
        ("n_atlas_rescue", _eq_count(atlas, "atlas_confidence_tier", "rescue")),
        ("n_atlas_core", _eq_count(atlas, "atlas_status", "core")),
        ("n_atlas_private", _eq_count(atlas, "atlas_status", "private")),
        ("n_broad_atlas", _true_count(atlas, "broad_or_unstable_cut_site")),
        ("n_mapped_sample_pas", _nonempty_count(mapping, "atlas_pas_id")),
        ("n_unmapped_sample_pas", len(mapping) - _nonempty_count(mapping, "atlas_pas_id")),
        ("n_requantified_detected", _eq_count(requant_status, "quantification_status", "quantified_detected")),
        ("n_requantified_zero", _eq_count(requant_status, "quantification_status", "quantified_zero")),
        ("n_requant_errors", _eq_count(requant_status, "quantification_status", "not_quantified_error")),
        ("n_requant_overlap_ambiguous", _eq_count(requant_status, "quantification_status", "window_overlap_ambiguous")),
    ]
    return pd.DataFrame(metrics, columns=["metric", "value"])


def _count_table(frame: pd.DataFrame, columns: list[str], value_name: str) -> pd.DataFrame:
    if frame.empty or any(column not in frame.columns for column in columns):
        return pd.DataFrame(columns=columns + [value_name])
    return frame.groupby(columns, dropna=False).size().reset_index(name=value_name)


def _motif_internal_counts(sequence: pd.DataFrame, confidence: pd.DataFrame) -> pd.DataFrame:
    if sequence.empty:
        return pd.DataFrame(columns=["confidence_tier", "motif_search_class", "internal_priming_flag", "n_pas"])
    table = sequence.copy()
    if not confidence.empty and "confidence_tier" in confidence.columns:
        table = table.merge(confidence[["pas_id", "confidence_tier"]].drop_duplicates("pas_id"), on="pas_id", how="left")
    if "confidence_tier" not in table:
        table["confidence_tier"] = ""
    table["motif_search_class"] = table.get("best_pas_motif_search_class", "").replace("", "no_motif")
    table["internal_priming_flag"] = table.get("internal_priming_flag", "false").astype(str).str.lower()
    return _count_table(table, ["confidence_tier", "motif_search_class", "internal_priming_flag"], "n_pas")


def _candidate_primary_flow(candidates: pd.DataFrame, feature: pd.DataFrame) -> pd.DataFrame:
    if candidates.empty or feature.empty:
        return pd.DataFrame(columns=["candidate_region_class", "primary_region_class", "n_candidate_links", "n_pas"])
    merged = candidates[["pas_id", "candidate_region_class"]].merge(
        feature[["pas_id", "primary_region_class"]].drop_duplicates("pas_id"), on="pas_id", how="left"
    )
    return (
        merged.groupby(["candidate_region_class", "primary_region_class"], dropna=False)
        .agg(n_candidate_links=("pas_id", "size"), n_pas=("pas_id", "nunique"))
        .reset_index()
    )


def _gene_testability_funnel(feature, confidence, testability) -> pd.DataFrame:
    values = [
        ("pas_total", len(feature)),
        ("pas_assigned_gene", _nonempty_count(feature, "primary_gene_id")),
        ("pas_primary", _eq_count(confidence, "confidence_tier", "primary")),
        ("pas_relaxed", _eq_count(confidence, "confidence_tier", "relaxed")),
        ("genes_with_pas", len(testability)),
        ("genes_multi_pas_primary", _true_count(testability, "has_multi_pas_primary")),
        ("genes_apa_testable_primary", _true_count(testability, "is_apa_testable_primary")),
        ("genes_apa_testable_relaxed", _true_count(testability, "is_apa_testable_relaxed")),
    ]
    return pd.DataFrame(values, columns=["stage", "value"])


def _representative_pas_examples(feature, confidence, sequence, external, max_per_region: int = 2) -> pd.DataFrame:
    if feature.empty:
        return pd.DataFrame()
    table = feature.copy()
    if not confidence.empty:
        table = table.merge(confidence[["pas_id", "confidence_tier"]].drop_duplicates("pas_id"), on="pas_id", how="left")
    if not sequence.empty:
        keep = [column for column in ["pas_id", "best_pas_motif", "best_pas_motif_search_class", "internal_priming_flag"] if column in sequence]
        table = table.merge(sequence[keep].drop_duplicates("pas_id"), on="pas_id", how="left", suffixes=("", "_sequence"))
    if not external.empty:
        supported = external[external["overlap_status"].eq("supported")].groupby("pas_id").size()
        table["n_external_supported_sources"] = table["pas_id"].map(supported).fillna(0).astype(int)
    table["total_pas_umi_num"] = pd.to_numeric(table.get("total_pas_umi", 0), errors="coerce").fillna(0)
    table = table.sort_values(["primary_region_class", "total_pas_umi_num", "pas_id"], ascending=[True, False, True])
    return table.groupby("primary_region_class", dropna=False).head(max_per_region).reset_index(drop=True)


def _adjacent_atlas_distances(atlas: pd.DataFrame) -> pd.DataFrame:
    columns = ["species", "genome_build", "gene_id", "chrom", "strand", "left_atlas_pas_id", "right_atlas_pas_id", "distance_bp"]
    if atlas.empty:
        return pd.DataFrame(columns=columns)
    rows = []
    grouped = atlas.copy()
    grouped["cut_num"] = pd.to_numeric(grouped["representative_cut_site"], errors="coerce")
    grouped = grouped.dropna(subset=["cut_num"])
    for keys, group in grouped.groupby(["species", "genome_build", "gene_id", "chrom", "strand"], dropna=False):
        group = group.sort_values(["cut_num", "atlas_pas_id"])
        records = list(group.itertuples(index=False))
        for left, right in zip(records, records[1:]):
            rows.append((*keys, left.atlas_pas_id, right.atlas_pas_id, int(right.cut_num - left.cut_num)))
    return pd.DataFrame(rows, columns=columns)


def _write_annotation_figures(qc_dir, region, flow, external, motif, funnel, examples, enabled):
    specs = [
        ("region_class_stacked_bar.png", lambda plt: _bar_plot(plt, region, "primary_region_class", "n_pas", "PAS region classes")),
        ("overlap_priority_flow.png", lambda plt: _flow_heatmap(plt, flow)),
        ("external_support_heatmap.png", lambda plt: _external_heatmap(plt, external)),
        ("motif_internal_priming_panel.png", lambda plt: _motif_plot(plt, motif)),
        ("gene_testability_funnel.png", lambda plt: _bar_plot(plt, funnel, "stage", "value", "Gene-PAS testability funnel")),
        ("representative_pas_browser.png", lambda plt: _browser_plot(plt, examples)),
    ]
    return _write_figures(qc_dir, specs, enabled)


def _write_atlas_figures(qc_dir, evidence, atlas, mapping, support, adjacent, requant, enabled):
    specs = [
        ("atlas_evidence_tiers.png", lambda plt: _bar_plot(plt, evidence, "evidence_tier_for_atlas", "n_sample_pas", "Atlas evidence tiers")),
        ("atlas_status.png", lambda plt: _bar_plot(plt, atlas, "atlas_confidence_tier", "n_atlas_pas", "Atlas confidence tiers")),
        ("atlas_mapping_status.png", lambda plt: _bar_plot(plt, mapping, "mapping_status", "n_sample_pas", "Sample-to-atlas mapping")),
        ("atlas_support_libraries.png", lambda plt: _bar_plot(plt, support, "support_library_count", "n_atlas_pas", "Library support per atlas PAS")),
        ("atlas_adjacent_site_distance.png", lambda plt: _hist_plot(plt, adjacent, "distance_bp", "Adjacent same-gene PAS distance")),
        ("atlas_requantification_status.png", lambda plt: _bar_plot(plt, requant, "quantification_status", "n_atlas_sample_rows", "Atlas re-quantification status")),
    ]
    return _write_figures(qc_dir, specs, enabled)


def _write_figures(qc_dir: Path, specs, enabled: bool) -> pd.DataFrame:
    rows = []
    if not enabled:
        return pd.DataFrame([{"figure": name, "status": "disabled", "reason": "write_plots=false"} for name, _ in specs])
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return pd.DataFrame([{"figure": name, "status": "not_generated", "reason": "matplotlib_unavailable"} for name, _ in specs])

    for name, writer in specs:
        try:
            writer(plt)
            plt.tight_layout()
            plt.savefig(qc_dir / name, dpi=160, bbox_inches="tight")
            plt.close()
            rows.append({"figure": name, "status": "generated", "reason": ""})
        except Exception as exc:
            plt.close("all")
            rows.append({"figure": name, "status": "not_generated", "reason": f"{type(exc).__name__}: {exc}"})
    return pd.DataFrame(rows)


def _bar_plot(plt, frame, category, value, title):
    if frame.empty or category not in frame or value not in frame:
        raise ValueError("no data")
    data = frame.groupby(category, dropna=False)[value].sum().sort_values(ascending=False)
    fig, ax = plt.subplots(figsize=(max(6, min(14, len(data) * 0.75)), 4.5))
    ax.bar([str(item) for item in data.index], data.values, color="#4C78A8")
    ax.set_title(title)
    ax.set_ylabel(value)
    ax.tick_params(axis="x", rotation=45)


def _hist_plot(plt, frame, column, title):
    if frame.empty or column not in frame:
        raise ValueError("no data")
    values = pd.to_numeric(frame[column], errors="coerce").dropna()
    if values.empty:
        raise ValueError("no numeric data")
    upper = float(values.quantile(0.99))
    clipped_values = [min(float(value), upper) for value in values]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(clipped_values, bins=40, color="#59A14F")
    ax.set_title(title)
    ax.set_xlabel("bp")
    ax.set_ylabel("count")


def _flow_heatmap(plt, frame):
    if frame.empty:
        raise ValueError("no data")
    matrix = frame.pivot_table(index="candidate_region_class", columns="primary_region_class", values="n_pas", aggfunc="sum", fill_value=0)
    fig, ax = plt.subplots(figsize=(max(6, matrix.shape[1]), max(4, matrix.shape[0] * 0.45)))
    image = ax.imshow(matrix.values, aspect="auto", cmap="Blues")
    ax.set_xticks(range(matrix.shape[1]), matrix.columns, rotation=45, ha="right")
    ax.set_yticks(range(matrix.shape[0]), matrix.index)
    ax.set_title("Candidate region to primary assignment")
    fig.colorbar(image, ax=ax, label="n PAS")


def _external_heatmap(plt, frame):
    if frame.empty:
        raise ValueError("no external resource data")
    matrix = frame.pivot_table(index="source_name", columns="overlap_status", values="n_pas_source_rows", aggfunc="sum", fill_value=0)
    fig, ax = plt.subplots(figsize=(max(6, matrix.shape[1] * 1.2), max(3, matrix.shape[0] * 0.6)))
    image = ax.imshow(matrix.values, aspect="auto", cmap="Purples")
    ax.set_xticks(range(matrix.shape[1]), matrix.columns, rotation=45, ha="right")
    ax.set_yticks(range(matrix.shape[0]), matrix.index)
    ax.set_title("External PAS support by source")
    fig.colorbar(image, ax=ax, label="PAS-source rows")


def _motif_plot(plt, frame):
    if frame.empty:
        raise ValueError("no sequence QC data")
    data = frame.groupby(["motif_search_class", "internal_priming_flag"])["n_pas"].sum().unstack(fill_value=0)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    data.plot(kind="bar", stacked=True, ax=ax, color=["#4C78A8", "#E45756"][: len(data.columns)])
    ax.set_title("Motif position and internal-priming risk")
    ax.set_ylabel("PAS")
    ax.tick_params(axis="x", rotation=30)


def _browser_plot(plt, examples):
    if examples.empty:
        raise ValueError("no representative PAS data")
    data = examples.head(20).copy()
    starts = pd.to_numeric(data["peak_start"], errors="coerce")
    ends = pd.to_numeric(data["peak_end"], errors="coerce")
    cuts = pd.to_numeric(data["cut_site"], errors="coerce")
    fig, ax = plt.subplots(figsize=(10, max(4, len(data) * 0.35)))
    for index, (_, row) in enumerate(data.iterrows()):
        start, end, cut = starts.iloc[index], ends.iloc[index], cuts.iloc[index]
        if pd.isna(start) or pd.isna(end) or pd.isna(cut):
            continue
        ax.hlines(index, start, end, color="#9C755F", linewidth=5)
        ax.vlines(cut, index - 0.25, index + 0.25, color="#E45756", linewidth=1.5)
    labels = [f"{row.get('primary_gene_symbol', '')}:{row.get('primary_region_class', '')}" for _, row in data.iterrows()]
    ax.set_yticks(range(len(data)), labels)
    ax.set_title("Representative PAS peak intervals and unchanged cut proxies")
    ax.set_xlabel("genomic coordinate")


def _annotation_report(annotation_path, summary, tables, figure_status):
    return _report(
        "PAS annotation QC report",
        [f"Annotation directory: `{annotation_path}`"],
        summary,
        tables,
        figure_status,
        "Counts describe annotation/sequence/external evidence; they do not by themselves prove biological PAS validity.",
    )


def _atlas_report(atlas_path, requant_path, summary, tables, figure_status):
    inputs = [f"Atlas directory: `{atlas_path}`"]
    if requant_path:
        inputs.append(f"Re-quantification directory: `{requant_path}`")
    return _report(
        "PAS atlas QC report",
        inputs,
        summary,
        tables,
        figure_status,
        "Projection missingness and BAM-level quantified zero remain separate states; sensitivity evidence is not silently mixed into the primary atlas.",
    )


def _report(title, inputs, summary, tables, figure_status, caveat):
    lines = [f"# {title}", "", "## Inputs", "", *[f"- {item}" for item in inputs], "", "## Summary", "", _markdown_table(summary), "", "## Methods", "", "Deterministic count summaries are derived directly from the module TSV contracts. Figures are generated from the same tables when matplotlib is available.", "", "## Artifacts", ""]
    lines.extend(f"- `{filename}`" for filename in tables)
    lines.extend(["", "## Figure status", "", _markdown_table(figure_status), "", "## Caveat", "", caveat, "", "## Environment", "", f"- Python: {platform.python_version()}", f"- pandas: {pd.__version__}", f"- numpy: {np.__version__}", "- Random seed: not applicable; summaries are deterministic.", ""])
    return "\n".join(lines)


def _markdown_table(frame: pd.DataFrame) -> str:
    if frame.empty:
        return "(none)"
    columns = list(frame.columns)
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in frame.astype(str).to_dict(orient="records"):
        lines.append("| " + " | ".join(row[column].replace("|", "\\|") for column in columns) + " |")
    return "\n".join(lines)


def _read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def _eq_count(frame: pd.DataFrame, column: str, value: str) -> int:
    return int(frame[column].astype(str).eq(value).sum()) if column in frame else 0


def _true_count(frame: pd.DataFrame, column: str) -> int:
    return int(frame[column].astype(str).str.lower().eq("true").sum()) if column in frame else 0


def _nonempty_count(frame: pd.DataFrame, column: str) -> int:
    return int(frame[column].astype(str).str.len().gt(0).sum()) if column in frame else 0
