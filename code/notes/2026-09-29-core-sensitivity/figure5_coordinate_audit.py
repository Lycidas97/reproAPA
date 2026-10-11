"""Audit annotated-PAS representative-coordinate sensitivity for Figure 5."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
NATIVE = ROOT / "analysis_code/release_2026-09-24/sources/F5/full-cohort-v2/data"
ANNOTATION = ROOT / "analysis_code/release_2026-09-24/sources/F5/annotation-pas-global-v1"
OUT = Path(__file__).resolve().parent / "figure5-coordinate-audit"
KEY = ["species", "route", "gene_id", "category", "chrom", "strand", "start0", "end0", "contrast"]
DISPLAYED_SITES = [
    {
        "site_label": "Eif4ebp2_Pcbp1_retained",
        "gene": "Eif4ebp2",
        "category": "RBP",
        "chrom": "chr10",
        "start0": 61269310,
        "end0": 61269318,
    },
    {
        "site_label": "Eif4ebp2_Pcbp1_changed_1",
        "gene": "Eif4ebp2",
        "category": "RBP",
        "chrom": "chr10",
        "start0": 61266098,
        "end0": 61266106,
    },
    {
        "site_label": "Eif4ebp2_Pcbp1_changed_2",
        "gene": "Eif4ebp2",
        "category": "RBP",
        "chrom": "chr10",
        "start0": 61264544,
        "end0": 61264552,
    },
    {
        "site_label": "Hmgn2_miR23ab3p",
        "gene": "Hmgn2",
        "category": "miRNA",
        "chrom": "chr4",
        "start0": 133692265,
        "end0": 133692273,
    },
]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def structural_status(row: pd.Series) -> str:
    if pd.isna(row.mask_json):
        return "absent"
    mask = np.asarray(json.loads(row.mask_json), dtype=int)
    exact = bool(row.exact) and np.array_equal(mask[:, 0], mask[:, 1])
    if not exact:
        return "boundary_uncertain"
    values = mask[:, 0]
    if np.all(values == 1):
        return "constitutive_included"
    if np.all(values == 0):
        return "constitutive_excluded"
    return "apa_variable"


def changed(frame: pd.DataFrame, prefix: str) -> pd.Series:
    a = frame[f"eq20_lo_{prefix}"]
    b = frame[f"eq50_lo_{prefix}"]
    return a.notna() & b.notna() & a.abs().ge(0.1 - 1e-12) & b.abs().ge(0.1 - 1e-12) & (a * b > 0)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mapping_path = ANNOTATION / "pas_annotation_mapping.tsv.gz"
    family_path = ANNOTATION / "family_migration.tsv"
    mapping = pd.read_csv(mapping_path, sep="\t")
    families = pd.read_csv(family_path, sep="\t")
    distance_summary = mapping.groupby(["species", "route"]).agg(
        pas=("atlas_pas_id", "size"),
        genes=("gene_id", "nunique"),
        representative_distance_median=("representative_distance", "median"),
        representative_distance_p90=("representative_distance", lambda x: x.quantile(0.9)),
        representative_distance_max=("representative_distance", "max"),
        representative_distance_gt50=("representative_distance", lambda x: int((x > 50).sum())),
        nearest_source_distance_max=("nearest_source_distance", "max"),
    ).reset_index()
    distance_summary.to_csv(OUT / "pas_representative_distance_summary.tsv", sep="\t", index=False)
    mapping[mapping.representative_distance.gt(50)].to_csv(
        OUT / "pas_representative_distance_gt50.tsv", sep="\t", index=False
    )
    families.groupby(["species", "route"]).agg(
        families=("gene_id", "size"),
        native_pas=("native_pas", "sum"),
        annotated_pas=("annotated_pas", "sum"),
        single_pas_after_merge=("single_pas_after_merge", "sum"),
    ).reset_index().to_csv(OUT / "family_migration_summary.tsv", sep="\t", index=False)

    joined_parts = []
    input_hashes = {str(mapping_path.relative_to(ROOT)): digest(mapping_path), str(family_path.relative_to(ROOT)): digest(family_path)}
    keep = KEY + ["gene_symbol", "mask_json", "exact", "variable", "eq20_lo", "eq50_lo"]
    for species in ["mouse", "human"]:
        native_path = NATIVE / species / "site_effects.tsv.gz"
        annotation_path = ANNOTATION / f"full-cohort-v2/data/{species}/site_effects.tsv.gz"
        input_hashes[str(native_path.relative_to(ROOT))] = digest(native_path)
        input_hashes[str(annotation_path.relative_to(ROOT))] = digest(annotation_path)
        native = pd.read_csv(native_path, sep="\t", usecols=keep, low_memory=False)
        annotation = pd.read_csv(annotation_path, sep="\t", usecols=keep, low_memory=False)
        assert not native.duplicated(KEY).any()
        assert not annotation.duplicated(KEY).any()
        native["structure"] = native.apply(structural_status, axis=1)
        annotation["structure"] = annotation.apply(structural_status, axis=1)
        joined = native.merge(
            annotation,
            on=KEY,
            how="outer",
            suffixes=("_native", "_annotation"),
            indicator=True,
            validate="one_to_one",
        )
        joined["strict_changed_native"] = changed(joined, "native")
        joined["strict_changed_annotation"] = changed(joined, "annotation")
        joined["structure_changed"] = (
            joined._merge.eq("both")
            & joined.structure_native.ne(joined.structure_annotation)
        )
        joined["strict_change_changed"] = (
            joined.strict_changed_native.ne(joined.strict_changed_annotation)
        )
        joined_parts.append(joined)
    joined = pd.concat(joined_parts, ignore_index=True)
    summary = joined.groupby(["species", "route", "contrast", "category"]).agg(
        native_rows=("_merge", lambda x: int(x.isin(["both", "left_only"]).sum())),
        annotation_rows=("_merge", lambda x: int(x.isin(["both", "right_only"]).sum())),
        shared_rows=("_merge", lambda x: int((x == "both").sum())),
        removed_rows=("_merge", lambda x: int((x == "left_only").sum())),
        added_rows=("_merge", lambda x: int((x == "right_only").sum())),
        structure_changed=("structure_changed", "sum"),
        strict_changed_native=("strict_changed_native", "sum"),
        strict_changed_annotation=("strict_changed_annotation", "sum"),
        strict_change_changed=("strict_change_changed", "sum"),
    ).reset_index()
    summary.to_csv(OUT / "site_classification_sensitivity.tsv", sep="\t", index=False)
    joined[joined.structure_changed | joined.strict_change_changed].to_csv(
        OUT / "site_classification_changes.tsv.gz", sep="\t", index=False
    )

    featured_genes = {"Hmgn2", "Eif4ebp2"}
    featured_mapping = mapping[mapping.gene.isin(featured_genes)].copy()
    featured_mapping.to_csv(OUT / "featured_gene_pas_mapping.tsv", sep="\t", index=False)
    featured_sites = joined[
        joined.gene_symbol_native.isin(featured_genes)
        | joined.gene_symbol_annotation.isin(featured_genes)
    ].copy()
    featured_sites.to_csv(OUT / "featured_gene_site_sensitivity.tsv.gz", sep="\t", index=False)
    displayed_parts = []
    for site in DISPLAYED_SITES:
        selected = featured_sites[
            (featured_sites.chrom == site["chrom"])
            & (featured_sites.start0 == site["start0"])
            & (featured_sites.end0 == site["end0"])
            & (featured_sites.category == site["category"])
            & (
                featured_sites.gene_symbol_native.eq(site["gene"])
                | featured_sites.gene_symbol_annotation.eq(site["gene"])
            )
        ].copy()
        if selected.empty:
            raise AssertionError(f"Displayed site missing from audit inputs: {site}")
        selected.insert(0, "site_label", site["site_label"])
        displayed_parts.append(selected)
    displayed = pd.concat(displayed_parts, ignore_index=True)
    displayed.to_csv(OUT / "displayed_site_sensitivity.tsv", sep="\t", index=False)
    displayed_summary = displayed.groupby("site_label").agg(
        routes=("route", "nunique"),
        contrasts=("contrast", "nunique"),
        structure_changed=("structure_changed", "sum"),
        strict_change_changed=("strict_change_changed", "sum"),
    ).reset_index()
    displayed_summary.to_csv(
        OUT / "displayed_site_summary.tsv", sep="\t", index=False
    )
    featured_summary = []
    for gene in sorted(featured_genes):
        m = featured_mapping[featured_mapping.gene.eq(gene)]
        s = featured_sites[
            featured_sites.gene_symbol_native.eq(gene)
            | featured_sites.gene_symbol_annotation.eq(gene)
        ]
        featured_summary.append(
            {
                "gene": gene,
                "mapped_pas": len(m),
                "max_representative_distance": int(m.representative_distance.max()),
                "representative_distance_gt50": int(m.representative_distance.gt(50).sum()),
                "site_contrast_rows": len(s),
                "structure_changed_rows": int(s.structure_changed.sum()),
                "strict_change_changed_rows": int(s.strict_change_changed.sum()),
                "annotation_strict_changed_rows": int(s.strict_changed_annotation.sum()),
            }
        )
    pd.DataFrame(featured_summary).to_csv(
        OUT / "featured_gene_summary.tsv", sep="\t", index=False
    )
    audit = {
        "status": "PASS",
        "mapping_rows": len(mapping),
        "representative_distance_gt50": int(mapping.representative_distance.gt(50).sum()),
        "families_becoming_single_pas": int(families.single_pas_after_merge.sum()),
        "site_contrast_rows_union": len(joined),
        "shared_site_contrast_rows": int(joined._merge.eq("both").sum()),
        "structure_changed_rows": int(joined.structure_changed.sum()),
        "strict_change_changed_rows": int(joined.strict_change_changed.sum()),
        "displayed_sites": int(displayed.site_label.nunique()),
        "displayed_structure_changed_rows": int(displayed.structure_changed.sum()),
        "displayed_strict_change_changed_rows": int(
            displayed.strict_change_changed.sum()
        ),
        "input_sha256": input_hashes,
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(distance_summary.to_string(index=False))
    print(pd.DataFrame(featured_summary).to_string(index=False))
    print(json.dumps({k: v for k, v in audit.items() if k != "input_sha256"}, indent=2))


if __name__ == "__main__":
    main()
