#!/usr/bin/env python3
"""Nested-background GO-BP analysis for recurrent reproAPA genes.

The analysis separates three gene universes:

G0: expression-measurable in both endpoints of at least two independent studies,
    without requiring multiple PAS;
G1: APA-testable in at least two independent studies;
G2: frozen same-direction recurrent APA genes at |delta WPP| >= 0.02.

Three comparisons locate where a functional signal enters:

1. G2 / G0: total functional composition relative to the cross-study
   expression-measurable transcriptome;
2. G1 / G0: functional composition associated with APA testability;
3. G2 / G1: recurrence-specific enrichment conditional on APA testability.

Hypergeometric tests are accompanied by covariate-stratified resampling. Each
resample draws without replacement within exact strata and therefore preserves
the target size and stratum counts exactly.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import hypergeom


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
SRC = ROOT / "sources"
OUT = ROOT / "outputs_nested_v2"
OUT.mkdir(exist_ok=True)

THRESHOLD = 0.02
N_PERM = 10_000
CHUNK = 500
SEED = 20260927
MIN_TERM, MAX_TERM = 10, 500


def load_base_module():
    path = ROOT / "scripts" / "run_go.py"
    spec = importlib.util.spec_from_file_location("run_go_base", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


BASE = load_base_module()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stable_seed(label: str) -> int:
    digest = hashlib.sha256(label.encode()).digest()
    return (SEED + int.from_bytes(digest[:4], "big")) % (2**32)


def gex_root() -> Path:
    env = os.environ.get("REPROAPA_GEX_ROOT")
    if env:
        root = Path(env)
    else:
        server = Path(
            "/root/nfsdata/apadetect/analysis_projects/"
            "figure3_apa_gex_state_20260919"
        )
        local = Path(
            "/Users/lycidas/Documents/ChatGPT/reproAPA/analysis_code/"
            "release_2026-09-24/sources/MOUSE/analysis/"
            "figure3-apa-gex-state-analysis-2026-09-19/raw"
        )
        root = server if server.exists() else local
    if not root.exists():
        raise FileNotFoundError(f"full-GEX root not found: {root}")
    return root


STUDY_FOLDERS = {
    "GSE103221": "gse103221_pluripotency_reprogramming_core",
    "GSE115943": "gse115943_pluripotency_reprogramming_core",
    "GSE118258": "gse118258_human_pluripotency_reprogramming",
    "GSE150311": "gse150311_human_pluripotency_reprogramming",
    "GSE178325": "gse178325_human_chemical_pluripotency_reprogramming",
    "GSE221739": "gse221739_human_pluripotency_reprogramming",
}


# Study composition follows the frozen Figure 3 recurrence analysis. Multiple
# Guan series remain one independent study unit.
DESIGN = {
    ("mouse", "P_S"): {
        "states": ("S", "P"),
        "studies": {"GSE103221": ["CIP"], "GSE115943": ["WOT"]},
    },
    ("mouse", "P_M"): {
        "states": ("M", "P"),
        "studies": {"GSE103221": ["CIP"], "GSE115943": ["WOT"]},
    },
    ("human", "P_S"): {
        "states": ("S", "P"),
        "studies": {
            "GSE118258": ["GSE118258"],
            "GSE178325": ["178325_0618", "178325_0809"],
        },
    },
    ("human", "P_M"): {
        "states": ("M", "P"),
        "studies": {
            "GSE118258": ["GSE118258"],
            "GSE150311": ["GSE150311"],
            "GSE221739": ["GSE221739"],
        },
    },
}


def read_full_gex(study: str) -> pd.DataFrame:
    path = gex_root() / STUDY_FOLDERS[study] / "full_gex_states.tsv.gz"
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(
        path,
        sep="\t",
        usecols=["gene_id", "gene_symbol", "route", "state", "GEX"],
    )


def build_g0() -> pd.DataFrame:
    """Build contrast-specific expression-measurable universes."""
    needed = sorted({study for d in DESIGN.values() for study in d["studies"]})
    cache = {}
    for study in needed:
        cache[study] = read_full_gex(study)
        print(f"loaded full GEX {study}: {len(cache[study]):,} rows", flush=True)

    rows = []
    for (species, contrast), design in DESIGN.items():
        state_a, state_b = design["states"]
        study_values = []
        symbols = {}
        for study, routes in design["studies"].items():
            x = cache[study]
            route_values = []
            for route in routes:
                y = x[x.route.eq(route)]
                if y.empty:
                    raise ValueError(f"missing route: {study} {route}")
                # Equal-library state means, then equal-state route means.
                z = (
                    y[y.state.isin([state_a, state_b])]
                    .groupby(["gene_id", "gene_symbol", "state"], as_index=False)
                    .GEX.mean()
                )
                wide = z.pivot_table(
                    index=["gene_id", "gene_symbol"],
                    columns="state",
                    values="GEX",
                    aggfunc="mean",
                ).reset_index()
                if state_a not in wide or state_b not in wide:
                    raise ValueError(f"missing endpoint state: {study} {route}")
                wide = wide.dropna(subset=[state_a, state_b]).copy()
                wide["route_GEX"] = wide[[state_a, state_b]].mean(axis=1)
                wide["route_min_GEX"] = wide[[state_a, state_b]].min(axis=1)
                route_values.append(wide[["gene_id", "gene_symbol", "route_GEX", "route_min_GEX"]])
            rv = pd.concat(route_values, ignore_index=True)
            # Multiple routes are alternative series within one study, not
            # independent recurrence units.
            sv = rv.groupby("gene_id", as_index=False).agg(
                gene_symbol=("gene_symbol", "first"),
                study_GEX=("route_GEX", "mean"),
                study_min_GEX=("route_min_GEX", "min"),
            )
            sv["study"] = study
            study_values.append(sv)
            symbols.update(dict(zip(sv.gene_id, sv.gene_symbol)))

        all_studies = pd.concat(study_values, ignore_index=True)
        summary = all_studies.groupby("gene_id", as_index=False).agg(
            expr_units=("study", "nunique"),
            mean_GEX=("study_GEX", "mean"),
            min_endpoint_GEX=("study_min_GEX", "min"),
        )
        summary = summary[summary.expr_units >= 2].copy()
        summary["gene_symbol"] = summary.gene_id.map(symbols)
        summary["species"] = species
        summary["contrast"] = contrast
        rows.append(summary)
        print(
            f"G0 {species} {contrast}: {len(summary):,} genes from "
            f"{len(design['studies'])} study units",
            flush=True,
        )
    return pd.concat(rows, ignore_index=True)


def add_expression_strata(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["expr_units_bin"] = df.expr_units.astype(int).astype(str)
    # Quintiles reduce expression-detection imbalance without conditioning on
    # PAS number or reference position, which are part of the G1/G0 contrast.
    df["mean_GEX_bin"] = (
        pd.qcut(df.mean_GEX.rank(method="first"), 5, labels=False)
        .astype("Int64")
        .astype(str)
    )
    df["expression_stratum"] = df.expr_units_bin + "|" + df.mean_GEX_bin
    return df


def term_family(bg_ids, annotations):
    members = defaultdict(set)
    for gene in bg_ids:
        for term in annotations.get(gene, set()):
            members[term].add(gene)
    return {
        term: genes
        for term, genes in members.items()
        if MIN_TERM <= len(genes) <= MAX_TERM
    }


def scan(species, contrast, comparison, target_set, target_ids, bg, annotations, terms):
    bg_ids = set(bg.gene_id)
    annotated_bg = {g for g in bg_ids if annotations.get(g)}
    target = set(target_ids)
    if not target <= bg_ids:
        missing = sorted(target - bg_ids)[:10]
        raise AssertionError((species, contrast, comparison, target_set, "target outside background", missing))
    fg = target & annotated_bg
    family = term_family(annotated_bg, annotations)
    term_ids = sorted(family)
    N, n = len(annotated_bg), len(fg)
    K = np.array([len(family[t]) for t in term_ids], dtype=int)
    hits = np.array([len(fg & family[t]) for t in term_ids], dtype=int)
    p = hypergeom.sf(hits - 1, N, K, n)
    q = BASE.bh(p)
    analysis_id = f"{species}|{contrast}|{comparison}|{target_set}"
    rows = []
    for i, term in enumerate(term_ids):
        rows.append(
            dict(
                analysis_id=analysis_id,
                species=species,
                contrast=contrast,
                comparison=comparison,
                target_set=target_set,
                term_id=term,
                term_name=terms[term]["name"],
                N=N,
                n=n,
                K=int(K[i]),
                hits=int(hits[i]),
                fold=(hits[i] / (n * K[i] / N)) if n and K[i] else np.nan,
                p_hyper=float(p[i]),
                q_hyper_family=float(q[i]),
            )
        )
    return pd.DataFrame(rows), annotated_bg, fg, family, term_ids, hits


def stratified_matched_null(
    analysis_id,
    bg,
    annotated_bg,
    fg,
    family,
    term_ids,
    observed_hits,
    stratum_col,
):
    """Sample exact target stratum counts without replacement."""
    bg = bg[bg.gene_id.isin(annotated_bg)].copy().set_index("gene_id")
    bg_ids = list(bg.index)
    gene_index = {g: i for i, g in enumerate(bg_ids)}
    term_index = {term: i for i, term in enumerate(term_ids)}
    # Counts can reach MAX_TERM=500; int8 wraps at 128 and invalidates the null.
    A = sparse.lil_matrix((len(bg_ids), len(term_ids)), dtype=np.int32)
    for term, genes in family.items():
        j = term_index[term]
        for gene in genes:
            A[gene_index[gene], j] = 1
    A = A.tocsr()

    pools = {
        key: np.array([gene_index[g] for g in frame.index], dtype=np.int64)
        for key, frame in bg.groupby(stratum_col, sort=False)
    }
    target_counts = Counter(bg.loc[sorted(fg), stratum_col])
    for key, count in target_counts.items():
        if key not in pools or count > len(pools[key]):
            raise AssertionError((analysis_id, key, count, len(pools.get(key, []))))

    rng = np.random.default_rng(stable_seed(analysis_id))
    exceed = np.zeros(len(term_ids), dtype=np.int64)
    expected_n = len(fg)
    min_n, max_n = expected_n, expected_n
    for start in range(0, N_PERM, CHUNK):
        size = min(CHUNK, N_PERM - start)
        row_parts, col_parts = [], []
        for row in range(size):
            chosen = []
            for key, count in target_counts.items():
                chosen.extend(rng.choice(pools[key], size=count, replace=False).tolist())
            unique_n = len(set(chosen))
            min_n, max_n = min(min_n, unique_n), max(max_n, unique_n)
            if unique_n != expected_n:
                raise AssertionError((analysis_id, expected_n, unique_n))
            row_parts.extend([row] * expected_n)
            col_parts.extend(chosen)
        D = sparse.csr_matrix(
            (
                np.ones(len(row_parts), dtype=np.int32),
                (np.asarray(row_parts), np.asarray(col_parts)),
            ),
            shape=(size, len(bg_ids)),
        )
        H = (D @ A).toarray()
        if np.any(H < 0) or np.any(H > min(expected_n, MAX_TERM)):
            raise AssertionError((analysis_id, "invalid null hit count"))
        exceed += (H >= observed_hits[None, :]).sum(axis=0)
    p = (1.0 + exceed) / (N_PERM + 1.0)
    q = BASE.bh(p)
    result = pd.DataFrame(
        {
            "analysis_id": analysis_id,
            "term_id": term_ids,
            "p_matched": p,
            "q_matched_family": q,
        }
    )
    diag = dict(
        analysis_id=analysis_id,
        stratum_col=stratum_col,
        n_strata=len(target_counts),
        target_n=expected_n,
        min_resample_n=min_n,
        max_resample_n=max_n,
        fixed_n=bool(min_n == expected_n == max_n),
        replacement=False,
        count_dtype=str(H.dtype),
    )
    return result, diag


def add_grid_q(df: pd.DataFrame, p_col: str, q_col: str) -> pd.DataFrame:
    df = df.copy()
    df[q_col] = BASE.bh(df[p_col].to_numpy())
    return df


def main():
    terms, aliases, ancestors = BASE.ontology()
    annotations = {
        species: BASE.load_go_genes(species, terms, aliases, ancestors)
        for species in ["mouse", "human"]
    }
    g0 = add_expression_strata(build_g0())

    events = pd.read_csv(DATA / "event_baseline_matrix.tsv.gz", sep="\t")
    recurrent = pd.read_csv(DATA / "recurrent_gene_sets.tsv.gz", sep="\t")
    recurrent = recurrent[
        recurrent.model.eq("original")
        & recurrent.omitted_study.eq("none")
        & recurrent.threshold.eq(THRESHOLD)
    ].copy()

    g1_parts = []
    for (species, contrast), frame in events.groupby(["species", "contrast"]):
        cov = BASE.add_bins(BASE.covariate_table(events, species, contrast))
        cov = cov[cov.eligible_units >= 2].copy()
        cov["species"] = species
        cov["contrast"] = contrast
        g1_parts.append(cov)
    g1 = pd.concat(g1_parts, ignore_index=True)

    membership_rows = []
    summary_rows = []
    result_frames = []
    matched_frames = []
    diagnostics = []
    attrition_rows = []

    for species, contrast in DESIGN:
        expr = g0[(g0.species == species) & (g0.contrast == contrast)].copy()
        apa = g1[(g1.species == species) & (g1.contrast == contrast)].copy()
        rec = recurrent[
            (recurrent.species == species) & (recurrent.contrast == contrast)
        ].copy()
        g0_ids = set(expr.gene_id)
        g1_raw_ids = set(apa.gene_id)
        g2_ids = set(rec.gene_id)
        g1_outside_g0 = g1_raw_ids - g0_ids
        g1_ids = g1_raw_ids & g0_ids
        if not g2_ids <= g1_ids <= g0_ids:
            raise AssertionError(
                dict(
                    species=species,
                    contrast=contrast,
                    g1_outside_g0=sorted(g1_outside_g0)[:20],
                    g2_outside_g1_raw=sorted(g2_ids - g1_raw_ids)[:20],
                    g2_outside_g0=sorted(g2_ids - g0_ids)[:20],
                )
            )
        apa_symbols = apa.set_index("gene_id").gene_symbol.to_dict()
        for gene in sorted(g1_outside_g0):
            attrition_rows.append(
                dict(
                    species=species,
                    contrast=contrast,
                    gene_id=gene,
                    gene_symbol=apa_symbols.get(gene, ""),
                    exclusion=(
                        "APA-testable but absent from the cross-study endpoint "
                        "nonzero-expression universe"
                    ),
                )
            )
        # Use one explicitly nested universe in all comparisons. PAS/GEX
        # discordances remain visible in the attrition table rather than
        # silently expanding the expression-defined G0 background.
        apa = apa[apa.gene_id.isin(g1_ids)].copy()

        ann = annotations[species]
        for layer, ids in [
            ("G0_expression", g0_ids),
            ("G1_APA_testable_raw", g1_raw_ids),
            ("G1_APA_testable_nested", g1_ids),
            ("G2_recurrent", g2_ids),
        ]:
            summary_rows.append(
                dict(
                    species=species,
                    contrast=contrast,
                    layer=layer,
                    n_genes=len(ids),
                    n_go_annotated=sum(bool(ann.get(g)) for g in ids),
                )
            )
        genes = sorted(g0_ids)
        rec_dir = rec.set_index("gene_id").direction.to_dict()
        symbols = expr.set_index("gene_id").gene_symbol.to_dict()
        for gene in genes:
            membership_rows.append(
                dict(
                    species=species,
                    contrast=contrast,
                    gene_id=gene,
                    gene_symbol=symbols.get(gene, ""),
                    in_G0=True,
                    in_G1_raw=gene in g1_raw_ids,
                    in_G1_nested=gene in g1_ids,
                    in_G2=gene in g2_ids,
                    recurrent_direction=rec_dir.get(gene, ""),
                )
            )

        targets = {
            "all": g2_ids,
            "proximal": set(rec.loc[rec.direction.eq("negative"), "gene_id"]),
            "distal": set(rec.loc[rec.direction.eq("positive"), "gene_id"]),
        }
        comparisons = []
        for target_name, target_ids in targets.items():
            comparisons.append(("G2_vs_G0", target_name, target_ids, expr, "expression_stratum"))
        comparisons.append(("G1_vs_G0", "all", g1_ids, expr, "expression_stratum"))
        for target_name, target_ids in targets.items():
            comparisons.append(("G2_vs_G1", target_name, target_ids, apa, "stratum_full"))

        for comparison, target_name, target_ids, background, stratum_col in comparisons:
            frame, annotated_bg, fg, family, term_ids, hits = scan(
                species,
                contrast,
                comparison,
                target_name,
                target_ids,
                background,
                ann,
                terms,
            )
            matched, diag = stratified_matched_null(
                frame.analysis_id.iloc[0],
                background,
                annotated_bg,
                fg,
                family,
                term_ids,
                hits,
                stratum_col,
            )
            frame = frame.merge(matched, on=["analysis_id", "term_id"], validate="one_to_one")
            result_frames.append(frame)
            diagnostics.append(diag)
            print(
                f"finished {frame.analysis_id.iloc[0]}: N={frame.N.iloc[0]:,}, "
                f"n={frame.n.iloc[0]:,}, terms={len(frame):,}",
                flush=True,
            )

    results = pd.concat(result_frames, ignore_index=True)
    results = add_grid_q(results, "p_hyper", "q_hyper_grid")
    results = add_grid_q(results, "p_matched", "q_matched_grid")
    primary_mask = results.contrast.eq("P_S") & results.comparison.eq("G2_vs_G0") & results.target_set.eq("all")
    results["q_hyper_primary_pair"] = np.nan
    results["q_matched_primary_pair"] = np.nan
    results.loc[primary_mask, "q_hyper_primary_pair"] = BASE.bh(results.loc[primary_mask, "p_hyper"])
    results.loc[primary_mask, "q_matched_primary_pair"] = BASE.bh(results.loc[primary_mask, "p_matched"])

    top = (
        results.sort_values(["analysis_id", "p_hyper"])
        .groupby("analysis_id", as_index=False, group_keys=False)
        .head(20)
    )
    significant = results[
        (results.q_hyper_family < 0.05) | (results.q_matched_family < 0.05)
    ].copy()

    pd.DataFrame(summary_rows).to_csv(OUT / "nested_universe_summary.tsv", sep="\t", index=False)
    pd.DataFrame(
        attrition_rows,
        columns=["species", "contrast", "gene_id", "gene_symbol", "exclusion"],
    ).to_csv(OUT / "nested_universe_attrition.tsv", sep="\t", index=False)
    pd.DataFrame(membership_rows).to_csv(
        OUT / "nested_universe_membership.tsv.gz", sep="\t", index=False
    )
    results.to_csv(OUT / "nested_go_enrichment.tsv.gz", sep="\t", index=False)
    top.to_csv(OUT / "nested_go_top_terms.tsv", sep="\t", index=False)
    significant.to_csv(OUT / "nested_go_significant.tsv", sep="\t", index=False)
    pd.DataFrame(diagnostics).to_csv(
        OUT / "nested_matched_diagnostics.tsv", sep="\t", index=False
    )

    inputs = {
        str(path.relative_to(ROOT)): sha256(path)
        for path in [
            DATA / "event_baseline_matrix.tsv.gz",
            DATA / "recurrent_gene_sets.tsv.gz",
            SRC / "go-basic.obo",
            SRC / "mouse_go_ids.tsv",
            SRC / "human_go_ids.tsv",
        ]
    }
    gex_inputs = {}
    for study, folder in STUDY_FOLDERS.items():
        path = gex_root() / folder / "full_gex_states.tsv.gz"
        gex_inputs[f"{study}:{path}"] = sha256(path)
    summary = dict(
        design="G2/G0 + G1/G0 + G2/G1",
        threshold=THRESHOLD,
        n_perm=N_PERM,
        seed=SEED,
        min_term=MIN_TERM,
        max_term=MAX_TERM,
        n_scans=int(results.analysis_id.nunique()),
        n_test_rows=int(len(results)),
        n_family_significant_hyper=int((results.q_hyper_family < 0.05).sum()),
        n_family_significant_matched=int((results.q_matched_family < 0.05).sum()),
        n_grid_significant_hyper=int((results.q_hyper_grid < 0.05).sum()),
        n_grid_significant_matched=int((results.q_matched_grid < 0.05).sum()),
        n_G1_outside_G0=int(len(attrition_rows)),
        all_resamples_fixed_n=bool(pd.DataFrame(diagnostics).fixed_n.all()),
        inputs=inputs,
        gex_inputs=gex_inputs,
    )
    (OUT / "run_summary_nested.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
