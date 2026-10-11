#!/usr/bin/env python3
"""Nested-background KEGG pathway analysis for recurrent reproAPA genes."""
from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
SRC = ROOT / "sources_kegg"
OUT = ROOT / "outputs_kegg_nested_v2"
OUT.mkdir(exist_ok=True)


def load_nested_module():
    path = ROOT / "scripts" / "run_go_nested.py"
    spec = importlib.util.spec_from_file_location("run_go_nested_base", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


NESTED = load_nested_module()
BASE = NESTED.BASE
THRESHOLD = NESTED.THRESHOLD


SPECIES_CONFIG = {
    "human": {
        "organism": "hsa",
        "gene_info": "Homo_sapiens.gene_info.gz",
        "ensembl_prefix": "ENSG",
        "species_suffix": " - Homo sapiens (human)",
    },
    "mouse": {
        "organism": "mmu",
        "gene_info": "Mus_musculus.gene_info.gz",
        "ensembl_prefix": "ENSMUSG",
        "species_suffix": " - Mus musculus (house mouse)",
    },
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_gene_info(species: str) -> dict[str, set[str]]:
    config = SPECIES_CONFIG[species]
    mapping: dict[str, set[str]] = {}
    with gzip.open(SRC / config["gene_info"], "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene_id = row["GeneID"]
            for value in row["dbXrefs"].split("|"):
                if not value.startswith("Ensembl:"):
                    continue
                ensembl = value.split(":", 1)[1].split(".", 1)[0]
                if ensembl.startswith(config["ensembl_prefix"]):
                    mapping.setdefault(ensembl, set()).add(gene_id)
    return mapping


def parse_two_column(path: Path) -> list[tuple[str, str]]:
    rows = []
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            left, right = line.rstrip("\n").split("\t")[:2]
            rows.append((left, right))
    return rows


def load_kegg(species: str):
    config = SPECIES_CONFIG[species]
    organism = config["organism"]
    ensembl_to_entrez = parse_gene_info(species)

    ncbi_to_kegg = {
        left.removeprefix("ncbi-geneid:"): right
        for left, right in parse_two_column(SRC / f"{organism}_ncbi_conv.tsv")
    }
    kegg_to_pathways: dict[str, set[str]] = {}
    for gene, pathway in parse_two_column(SRC / f"{organism}_gene_pathway.tsv"):
        pathway_id = pathway.removeprefix("path:")
        kegg_to_pathways.setdefault(gene, set()).add(pathway_id)

    annotations: dict[str, set[str]] = {}
    mapping_rows = []
    for ensembl, entrez_ids in ensembl_to_entrez.items():
        pathways = set()
        kegg_ids = set()
        for entrez in entrez_ids:
            kegg = ncbi_to_kegg.get(entrez)
            if kegg:
                kegg_ids.add(kegg)
                pathways.update(kegg_to_pathways.get(kegg, set()))
        if pathways:
            annotations[ensembl] = pathways
        mapping_rows.append(
            {
                "species": species,
                "ensembl_gene_id": ensembl,
                "entrez_gene_ids": ";".join(sorted(entrez_ids)),
                "kegg_gene_ids": ";".join(sorted(kegg_ids)),
                "n_pathways": len(pathways),
            }
        )

    terms = {}
    suffix = config["species_suffix"]
    for pathway_id, name in parse_two_column(SRC / f"{organism}_pathways.tsv"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
        terms[pathway_id] = {"name": name}
    return annotations, terms, mapping_rows


def main() -> None:
    annotations = {}
    terms = {}
    mapping_rows = []
    for species in ["mouse", "human"]:
        annotations[species], terms[species], rows = load_kegg(species)
        mapping_rows.extend(rows)
        print(
            f"loaded KEGG {species}: {len(annotations[species]):,} Ensembl genes, "
            f"{len(terms[species]):,} pathways",
            flush=True,
        )

    g0 = NESTED.add_expression_strata(NESTED.build_g0())
    events = pd.read_csv(DATA / "event_baseline_matrix.tsv.gz", sep="\t")
    recurrent = pd.read_csv(DATA / "recurrent_gene_sets.tsv.gz", sep="\t")
    recurrent = recurrent[
        recurrent.model.eq("original")
        & recurrent.omitted_study.eq("none")
        & recurrent.threshold.eq(THRESHOLD)
    ].copy()

    g1_parts = []
    for (species, contrast), _ in events.groupby(["species", "contrast"]):
        cov = BASE.add_bins(BASE.covariate_table(events, species, contrast))
        cov = cov[cov.eligible_units >= 2].copy()
        cov["species"] = species
        cov["contrast"] = contrast
        g1_parts.append(cov)
    g1 = pd.concat(g1_parts, ignore_index=True)

    summary_rows = []
    attrition_rows = []
    result_frames = []
    diagnostics = []

    for species, contrast in NESTED.DESIGN:
        expr = g0[(g0.species == species) & (g0.contrast == contrast)].copy()
        apa = g1[(g1.species == species) & (g1.contrast == contrast)].copy()
        rec = recurrent[
            recurrent.species.eq(species) & recurrent.contrast.eq(contrast)
        ].copy()

        g0_ids = set(expr.gene_id)
        g1_raw_ids = set(apa.gene_id)
        g2_ids = set(rec.gene_id)
        g1_outside_g0 = g1_raw_ids - g0_ids
        g1_ids = g1_raw_ids & g0_ids
        if not g2_ids <= g1_ids <= g0_ids:
            raise AssertionError(
                {
                    "species": species,
                    "contrast": contrast,
                    "g1_outside_g0": sorted(g1_outside_g0)[:20],
                    "g2_outside_g1_nested": sorted(g2_ids - g1_ids)[:20],
                }
            )
        apa_symbols = apa.set_index("gene_id").gene_symbol.to_dict()
        for gene in sorted(g1_outside_g0):
            attrition_rows.append(
                {
                    "species": species,
                    "contrast": contrast,
                    "gene_id": gene,
                    "gene_symbol": apa_symbols.get(gene, ""),
                    "exclusion": (
                        "APA-testable but absent from the cross-study endpoint "
                        "nonzero-expression universe"
                    ),
                }
            )
        apa = apa[apa.gene_id.isin(g1_ids)].copy()

        ann = annotations[species]
        for layer, ids in [
            ("G0_expression", g0_ids),
            ("G1_APA_testable_raw", g1_raw_ids),
            ("G1_APA_testable_nested", g1_ids),
            ("G2_recurrent", g2_ids),
        ]:
            summary_rows.append(
                {
                    "species": species,
                    "contrast": contrast,
                    "layer": layer,
                    "n_genes": len(ids),
                    "n_kegg_annotated": sum(bool(ann.get(g)) for g in ids),
                }
            )

        targets = {
            "all": g2_ids,
            "proximal": set(rec.loc[rec.direction.eq("negative"), "gene_id"]),
            "distal": set(rec.loc[rec.direction.eq("positive"), "gene_id"]),
        }
        comparisons = []
        for target_name, target_ids in targets.items():
            comparisons.append(
                ("G2_vs_G0", target_name, target_ids, expr, "expression_stratum")
            )
        comparisons.append(
            ("G1_vs_G0", "all", g1_ids, expr, "expression_stratum")
        )
        for target_name, target_ids in targets.items():
            comparisons.append(
                ("G2_vs_G1", target_name, target_ids, apa, "stratum_full")
            )

        for comparison, target_name, target_ids, background, stratum_col in comparisons:
            frame, annotated_bg, fg, family, term_ids, hits = NESTED.scan(
                species,
                contrast,
                comparison,
                target_name,
                target_ids,
                background,
                ann,
                terms[species],
            )
            matched, diagnostic = NESTED.stratified_matched_null(
                frame.analysis_id.iloc[0],
                background,
                annotated_bg,
                fg,
                family,
                term_ids,
                hits,
                stratum_col,
            )
            frame = frame.merge(
                matched, on=["analysis_id", "term_id"], validate="one_to_one"
            )
            result_frames.append(frame)
            diagnostics.append(diagnostic)
            print(
                f"finished {frame.analysis_id.iloc[0]}: N={frame.N.iloc[0]:,}, "
                f"n={frame.n.iloc[0]:,}, pathways={len(frame):,}",
                flush=True,
            )

    results = pd.concat(result_frames, ignore_index=True)
    results = NESTED.add_grid_q(results, "p_hyper", "q_hyper_grid")
    results = NESTED.add_grid_q(results, "p_matched", "q_matched_grid")
    primary_mask = (
        results.contrast.eq("P_S")
        & results.comparison.eq("G2_vs_G0")
        & results.target_set.eq("all")
    )
    results["q_hyper_primary_pair"] = np.nan
    results["q_matched_primary_pair"] = np.nan
    results.loc[primary_mask, "q_hyper_primary_pair"] = BASE.bh(
        results.loc[primary_mask, "p_hyper"]
    )
    results.loc[primary_mask, "q_matched_primary_pair"] = BASE.bh(
        results.loc[primary_mask, "p_matched"]
    )

    top = (
        results.sort_values(["analysis_id", "p_hyper"])
        .groupby("analysis_id", as_index=False, group_keys=False)
        .head(20)
    )
    significant = results[
        results.q_hyper_family.lt(0.05) | results.q_matched_family.lt(0.05)
    ].copy()

    pd.DataFrame(summary_rows).to_csv(
        OUT / "kegg_universe_summary.tsv", sep="\t", index=False
    )
    pd.DataFrame(
        attrition_rows,
        columns=["species", "contrast", "gene_id", "gene_symbol", "exclusion"],
    ).to_csv(OUT / "kegg_universe_attrition.tsv", sep="\t", index=False)
    pd.DataFrame(mapping_rows).to_csv(
        OUT / "kegg_ensembl_mapping.tsv.gz", sep="\t", index=False
    )
    results.to_csv(OUT / "kegg_enrichment.tsv.gz", sep="\t", index=False)
    top.to_csv(OUT / "kegg_top_terms.tsv", sep="\t", index=False)
    significant.to_csv(OUT / "kegg_significant.tsv", sep="\t", index=False)
    pd.DataFrame(diagnostics).to_csv(
        OUT / "kegg_matched_diagnostics.tsv", sep="\t", index=False
    )

    source_files = sorted(path for path in SRC.iterdir() if path.is_file())
    summary = {
        "design": "G2/G0 + G1/G0 + G2/G1",
        "database": "KEGG pathway",
        "threshold": THRESHOLD,
        "n_perm": NESTED.N_PERM,
        "seed": NESTED.SEED,
        "min_pathway": NESTED.MIN_TERM,
        "max_pathway": NESTED.MAX_TERM,
        "n_scans": int(results.analysis_id.nunique()),
        "n_test_rows": int(len(results)),
        "n_family_significant_hyper": int(results.q_hyper_family.lt(0.05).sum()),
        "n_family_significant_matched": int(results.q_matched_family.lt(0.05).sum()),
        "n_grid_significant_hyper": int(results.q_hyper_grid.lt(0.05).sum()),
        "n_grid_significant_matched": int(results.q_matched_grid.lt(0.05).sum()),
        "n_G1_outside_G0": int(len(attrition_rows)),
        "all_resamples_fixed_n": bool(pd.DataFrame(diagnostics).fixed_n.all()),
        "inputs": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in [
                DATA / "event_baseline_matrix.tsv.gz",
                DATA / "recurrent_gene_sets.tsv.gz",
            ]
        },
        "sources": {
            str(path.relative_to(ROOT)): sha256(path) for path in source_files
        },
    }
    (OUT / "run_summary_kegg.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
