#!/usr/bin/env python3
"""GO-BP functional analysis of reproAPA recurrent APA gene sets.

Question: are recurrent (>=2 independent studies, same-direction) P-S APA genes
concentrated in particular biological processes, and does this differ by response
direction (proximal vs distal)?

Design (pre-declared, no threshold fishing):
- Targets: Figure 3 recurrent sets, per species, per contrast, split by direction.
  Primary contrast = P_S, threshold 0.02 (frozen). P_M reported separately.
- Background universe: genes eligible for cross-study recurrence evaluation,
  taken from the frozen event observability matrix. Primary = eligible in >=2
  studies (the minimum needed for recurrence); sensitivity = all eligible genes.
- Annotation: Ensembl release 116 gene->GO ids + go-basic.obo; GO Biological
  Process only; ancestor propagation; root excluded.
- Test: one-sided hypergeometric over-representation; term family = 10..500
  genes within the tested annotated background; BH per family.
- Testability-matched sensitivity: frequency-matched resampling of the annotated
  background on APA-measurability covariates (eligible_units, mean_GEX,
  APA_coverage, PAS_number, reference_WPP), 10,000 draws, term-wise empirical
  upper-tail P, BH over the same full term family.

No visualisation in this run.
"""
from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import hypergeom, mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
SRC = ROOT / "sources"
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)

CONTRASTS = ["P_S", "P_M"]
PRIMARY_CONTRAST = "P_S"
THRESHOLD = 0.02
N_PERM = 10_000
CHUNK = 1_000
SEED = 20260926
MIN_TERM, MAX_TERM = 10, 500
ROOT_TERM = "GO:0008150"  # biological_process


# ---------------------------------------------------------------- utilities
def bh(p):
    a = np.asarray(p, dtype=float)
    order = np.argsort(a, kind="stable")
    q = np.empty(len(a), dtype=float)
    ranked = a[order] * len(a) / np.arange(1, len(a) + 1)
    q[order] = np.minimum(1.0, np.minimum.accumulate(ranked[::-1])[::-1])
    return q


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ontology():
    """Parse go-basic.obo: term namespace + ancestor closure."""
    source = (SRC / "go-basic.obo").read_text()
    terms, aliases = {}, {}
    for block in source.split("\n[Term]\n")[1:]:
        rec = defaultdict(list)
        for line in block.split("\n[")[0].splitlines():
            if ": " in line:
                k, v = line.split(": ", 1)
                rec[k].append(v)
        if not rec["id"] or rec["is_obsolete"] == ["true"]:
            continue
        term = rec["id"][0]
        parents = [x.split()[0] for x in rec["is_a"]]
        parents += [x.split()[1] for x in rec["relationship"] if x.startswith("part_of ")]
        terms[term] = dict(name=rec["name"][0], namespace=rec["namespace"][0], parents=parents)
        for a in rec["alt_id"]:
            aliases[a] = term

    @lru_cache(None)
    def ancestors(t):
        if t not in terms:
            return frozenset()
        return frozenset([t]).union(*(ancestors(p) for p in terms[t]["parents"]))

    return terms, aliases, ancestors


def load_go_genes(species, terms, aliases, ancestors):
    """gene -> set of BP terms (with ancestors); root excluded."""
    path = SRC / f"{species}_go_ids.tsv"
    raw = pd.read_csv(path, sep="\t", header=None, names=["gene_id", "go_id"], dtype=str,
                      keep_default_na=False)
    raw = raw[raw.gene_id.ne("[success]")]
    out = defaultdict(set)
    for g, t in zip(raw.gene_id, raw.go_id):
        if not t:
            continue
        t = aliases.get(t, t)
        if t not in terms:
            continue
        for a in ancestors(t):
            if terms[a]["namespace"] == "biological_process" and a != ROOT_TERM:
                out[g].add(a)
    return out


def covariate_table(events, species, contrast):
    """Gene-level APA-measurability covariates from the frozen event matrix."""
    e = events[(events.species == species) & (events.contrast == contrast)]
    g = e.groupby("gene_id", as_index=False).agg(
        gene_symbol=("gene_symbol", "first"),
        eligible_units=("eligible_units", "max"),
        mean_GEX=("mean_GEX", "median"),
        APA_coverage=("APA_coverage", "median"),
        PAS_number=("PAS_number", "median"),
        reference_WPP=("reference_WPP", "median"),
    )
    return g


def add_bins(g):
    g = g.copy()
    g["eligible_units_bin"] = np.where(g.eligible_units >= 3, "3plus", "2")
    for c in ["mean_GEX", "APA_coverage", "PAS_number", "reference_WPP"]:
        vals = g[c]
        try:
            g[c + "_bin"] = pd.qcut(vals.rank(method="first"), 3, labels=False).astype("Int64").astype(str)
        except ValueError:
            g[c + "_bin"] = "0"
    g["stratum_full"] = (g.eligible_units_bin + "|" + g.PAS_number_bin + "|" + g.APA_coverage_bin
                         + "|" + g.mean_GEX_bin + "|" + g.reference_WPP_bin)
    g["stratum_relax1"] = (g.eligible_units_bin + "|" + g.PAS_number_bin + "|" + g.APA_coverage_bin
                           + "|" + g.mean_GEX_bin)
    g["stratum_relax2"] = g.eligible_units_bin + "|" + g.PAS_number_bin + "|" + g.APA_coverage_bin
    g["stratum_relax3"] = g.eligible_units_bin + "|" + g.PAS_number_bin
    g["stratum_relax4"] = g.eligible_units_bin
    return g


# ---------------------------------------------------------------- main
def main():
    print("loading ontology ...", flush=True)
    terms, aliases, ancestors = ontology()
    go_genes = {}
    for sp in ["mouse", "human"]:
        go_genes[sp] = load_go_genes(sp, terms, aliases, ancestors)
        print(f"  {sp}: {len(go_genes[sp])} GO-annotated genes", flush=True)

    events = pd.read_csv(DATA / "event_baseline_matrix.tsv.gz", sep="\t")
    rec = pd.read_csv(DATA / "recurrent_gene_sets.tsv.gz", sep="\t")
    rec = rec[(rec.model == "original") & (rec.omitted_study == "none") & (rec.threshold == THRESHOLD)]
    gene_eff = events.groupby(["species", "contrast", "gene_id"]).APA_delta_mean.median()

    # sanity: direction sign matches effect sign
    chk = events.merge(rec[["species", "contrast", "gene_id", "direction"]],
                       on=["species", "contrast", "gene_id"], how="inner")
    for d, sgn in [("negative", -1), ("positive", 1)]:
        m = chk[chk.direction == d].APA_delta_mean
        assert (np.sign(m) == sgn).mean() > 0.98, (d, (np.sign(m) == sgn).mean())

    cov = {}
    for sp in ["mouse", "human"]:
        for ct in CONTRASTS:
            cov[(sp, ct)] = add_bins(covariate_table(events, sp, ct))

    target_sets, balance_rows, universe_rows = [], [], []
    results, member_rows, matched_rows, match_diag = [], [], [], []

    for sp in ["mouse", "human"]:
        for ct in CONTRASTS:
            r = rec[(rec.species == sp) & (rec.contrast == ct)]
            if r.empty:
                continue
            r_meta = r.set_index("gene_id")
            cg = cov[(sp, ct)]
            genes = set(cg.gene_id)
            ann_genes = {g: go_genes[sp].get(g, set()) for g in genes}
            annotated_all = {g for g, s in ann_genes.items() if s}
            universes = {
                "eligible_ge2": cg[cg.eligible_units >= 2],
                "eligible_all": cg,
            }
            # covariate balance figures (target vs primary background)
            for name, sub in universes.items():
                universe_rows.append(dict(species=sp, contrast=ct, universe=name,
                                          n_genes=len(sub), n_annotated=len(set(sub.gene_id) & annotated_all)))

            for direction, label in [("all", "all"), ("negative", "proximal"), ("positive", "distal")]:
                t = r if direction == "all" else r[r.direction == direction]
                tids = list(t.gene_id)
                if not tids:
                    continue
                target_sets.append(dict(species=sp, contrast=ct, target_set=label,
                                        direction=direction, n_target=len(tids),
                                        n_target_annotated=len(set(tids) & annotated_all)))
                for name, sub in universes.items():
                    bg_ann = {g for g in sub.gene_id if ann_genes.get(g)}
                    N = len(bg_ann)
                    fg = set(tids) & bg_ann
                    n = len(fg)
                    assert set(tids) <= set(sub.gene_id), (sp, ct, name, "target not in universe")
                    # term family within this annotated background
                    term_members = defaultdict(set)
                    for g in bg_ann:
                        for term in ann_genes[g]:
                            term_members[term].add(g)
                    fam = {t_: m for t_, m in term_members.items() if MIN_TERM <= len(m) <= MAX_TERM}
                    term_ids = sorted(fam)
                    if not term_ids:
                        continue
                    K = np.array([len(fam[t_]) for t_ in term_ids])
                    hits = np.array([len(fg & fam[t_]) for t_ in term_ids])
                    p = hypergeom.sf(hits - 1, N, K, n)
                    q = bh(p)
                    for i, t_ in enumerate(term_ids):
                        results.append(dict(species=sp, contrast=ct, target_set=label, universe=name,
                                            method="hypergeometric", term_id=t_,
                                            term_name=terms[t_]["name"], N=N, n=n, K=int(K[i]),
                                            hits=int(hits[i]),
                                            fold=(hits[i] / (n * K[i] / N)) if n and K[i] else np.nan,
                                            p=float(p[i]), q=float(q[i])))
                    # member genes for nominally enriched terms (direction resolved later)
                    for i, t_ in enumerate(term_ids):
                        if hits[i] > 0 and (p[i] < 0.005 or q[i] < 0.05):
                            for g in sorted(fg & fam[t_]):
                                member_rows.append(dict(species=sp, contrast=ct, target_set=label,
                                                        universe=name, method="hypergeometric", term_id=t_,
                                                        term_name=terms[t_]["name"], gene_id=g,
                                                        gene_symbol=r_meta.gene_symbol.get(g, ""),
                                                        direction=r_meta.direction.get(g, ""),
                                                        positive_support=r_meta.positive_support.get(g, 0),
                                                        negative_support=r_meta.negative_support.get(g, 0),
                                                        APA_delta_median=float(gene_eff.get((sp, ct, g), np.nan))))
                    # ---- testability-matched resampling null (primary universe only)
                    if name == "eligible_ge2":
                        sub2 = cg[cg.gene_id.isin(bg_ann)].copy()
                        fg_list = sorted(fg)
                        if fg_list and term_ids:
                            matched, fallback_frac = matched_null(sub2, fg_list, term_ids, fam, fg, hits, N, K, n)
                            match_diag.append(dict(species=sp, contrast=ct, target_set=label,
                                                   n_target_annotated=len(fg_list),
                                                   unmatched_fraction=fallback_frac))
                            for i, t_ in enumerate(term_ids):
                                matched_rows.append(dict(species=sp, contrast=ct, target_set=label,
                                                         term_id=t_, term_name=terms[t_]["name"],
                                                         N=N, n=n, K=int(K[i]), hits=int(hits[i]),
                                                         p_hyper=float(p[i]), q_hyper=float(q[i]),
                                                         p_matched=float(matched[i])))
                                if matched[i] < 0.005 and hits[i] > 0:
                                    for g in sorted(fg & fam[t_]):
                                        member_rows.append(dict(species=sp, contrast=ct, target_set=label,
                                                                universe="eligible_ge2_matched", method="matched",
                                                                term_id=t_, term_name=terms[t_]["name"], gene_id=g,
                                                                gene_symbol=r_meta.gene_symbol.get(g, ""),
                                                                direction=r_meta.direction.get(g, ""),
                                                                positive_support=r_meta.positive_support.get(g, 0),
                                                                negative_support=r_meta.negative_support.get(g, 0),
                                                                APA_delta_median=float(gene_eff.get((sp, ct, g), np.nan))))
            # covariate balance for each target set vs ge2 background
            for direction, label in [("all", "all"), ("negative", "proximal"), ("positive", "distal")]:
                t = r if direction == "all" else r[r.direction == direction]
                sub = universes["eligible_ge2"]
                for c in ["eligible_units", "mean_GEX", "APA_coverage", "PAS_number", "reference_WPP"]:
                    a = cg[cg.gene_id.isin(t.gene_id)][c].dropna()
                    b = sub[c].dropna()
                    if len(a) > 2 and len(b) > 2:
                        u, pv = mannwhitneyu(a, b, alternative="two-sided")
                    else:
                        pv = np.nan
                    balance_rows.append(dict(species=sp, contrast=ct, target_set=label, covariate=c,
                                             target_median=float(a.median()) if len(a) else np.nan,
                                             background_median=float(b.median()) if len(b) else np.nan,
                                             p_mannwhitney=float(pv) if pv == pv else np.nan))
            print(f"  {sp} {ct}: done", flush=True)

    # assemble matched table with BH over the full term family per set
    matched = pd.DataFrame(matched_rows)
    if not matched.empty:
        matched["p_matched_full_family"] = 1.0
        out = []
        for (sp, ct, ts), f in matched.groupby(["species", "contrast", "target_set"]):
            # terms with k=0 never enter; family size = hypergeometric family size
            fam_size = sum(1 for x in results if x["species"] == sp and x["contrast"] == ct
                           and x["target_set"] == ts and x["universe"] == "eligible_ge2")
            f = f.copy()
            q = bh(np.concatenate([f.p_matched.to_numpy(), np.ones(max(0, fam_size - len(f)))]))[:len(f)]
            f["q_matched"] = q
            out.append(f)
        matched = pd.concat(out, ignore_index=True)

    # ---- Figure 4 shared gene per-gene annotation
    shared = pd.read_csv(DATA / "figure4_shared_membership.tsv", sep="\t")
    f18 = pd.read_csv(DATA / "frozen18_outcomes.tsv", sep="\t")
    f18 = f18[f18.view == "reapplied"].drop_duplicates("human_gene_id").set_index("human_gene_id")
    f4_rows = []
    res_df = pd.DataFrame(results)
    sig_terms = set(res_df[(res_df.universe == "eligible_ge2") & (res_df.q < 0.05)].term_id) if not res_df.empty else set()
    nom_terms = set(res_df[(res_df.universe == "eligible_ge2") & (res_df.p < 0.005)].term_id) if not res_df.empty else set()
    for row in shared.itertuples():
        gs = go_genes["human"].get(row.human_gene_id, set())
        names = ";".join(terms[t]["name"] for t in sorted(gs))
        f4_rows.append(dict(human_gene_id=row.human_gene_id, human_symbol=row.human_symbol,
                            set_relation=row.set_relation,
                            in_preexisting_eight=bool(f18.in_preexisting_eight.get(row.human_gene_id, False)),
                            n_go_bp_terms=len(gs),
                            go_bp_terms=";".join(sorted(gs)),
                            go_bp_term_names=names,
                            in_significant_recurrent_term=bool(gs & sig_terms),
                            significant_terms=";".join(sorted(gs & sig_terms)),
                            in_nominal_recurrent_term=bool(gs & nom_terms),
                            nominal_terms=";".join(sorted(gs & nom_terms))))
    f4 = pd.DataFrame(f4_rows)

    # ---- top terms per group (hypergeometric primary + matched null)
    top_rows = []
    if not res_df.empty:
        for key, f in res_df.groupby(["species", "contrast", "target_set", "universe"]):
            for rank, r_ in enumerate(f.sort_values("p").head(20).itertuples(), start=1):
                top_rows.append(dict(species=key[0], contrast=key[1], target_set=key[2], universe=key[3],
                                     rank=rank,
                                     term_id=r_.term_id, term_name=r_.term_name, N=r_.N, n=r_.n, K=r_.K,
                                     hits=r_.hits, fold=r_.fold, p=r_.p, q=r_.q))
    pd.DataFrame(top_rows).to_csv(OUT / "go_top_terms.tsv", sep="\t", index=False)

    # ---- write outputs
    pd.DataFrame(target_sets).to_csv(OUT / "target_sets.tsv", sep="\t", index=False)
    pd.DataFrame(universe_rows).to_csv(OUT / "background_universes.tsv", sep="\t", index=False)
    res_df.to_csv(OUT / "go_enrichment_all.tsv", sep="\t", index=False)
    if not res_df.empty:
        res_df[(res_df.universe == "eligible_ge2") & (res_df.q < 0.05)].to_csv(
            OUT / "go_enrichment_significant.tsv", sep="\t", index=False)
    pd.DataFrame(member_rows).to_csv(OUT / "go_term_members_nominal.tsv", sep="\t", index=False)
    matched.to_csv(OUT / "go_matched_null.tsv", sep="\t", index=False)
    pd.DataFrame(balance_rows).to_csv(OUT / "covariate_balance.tsv", sep="\t", index=False)
    pd.DataFrame(match_diag).to_csv(OUT / "matched_null_diagnostics.tsv", sep="\t", index=False)
    f4.to_csv(OUT / "figure4_shared_gene_go.tsv", sep="\t", index=False)

    summary = dict(
        go_data_version=next(l for l in (SRC / "go-basic.obo").read_text().splitlines()
                             if l.startswith("data-version:")),
        ensembl_release=json.loads((SRC / "manifest.json").read_text()).get("ensembl_release"),
        threshold=THRESHOLD, n_perm=N_PERM, seed=SEED,
        min_term=MIN_TERM, max_term=MAX_TERM,
        n_tests=int(len(res_df)),
        n_significant_ge2=int(((res_df.universe == "eligible_ge2") & (res_df.q < 0.05)).sum()) if not res_df.empty else 0,
        n_matched_significant=int((matched.q_matched < 0.05).sum()) if not matched.empty else 0,
        inputs={str(p.relative_to(ROOT)): sha256(p) for p in
                [DATA / "event_baseline_matrix.tsv.gz", DATA / "recurrent_gene_sets.tsv.gz",
                 DATA / "figure4_shared_membership.tsv", DATA / "frozen18_outcomes.tsv"]},
    )
    (OUT / "run_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


def matched_null(sub2, fg_list, term_ids, fam, fg, hits, N, K, n):
    """Frequency-matched resampling null over annotated background genes."""
    rng = np.random.default_rng(SEED)
    bg = sub2.set_index("gene_id")
    term_index = {t_: i for i, t_ in enumerate(term_ids)}
    bg_ids = list(bg.index)
    gi = {g: i for i, g in enumerate(bg_ids)}
    A = sparse.lil_matrix((len(bg_ids), len(term_ids)), dtype=np.int8)
    for t_ in term_ids:
        j = term_index[t_]
        for g in fam[t_]:
            if g in gi:
                A[gi[g], j] = 1
    A = A.tocsr()
    # pools per covariate stratum, hierarchical fallback
    def make_pools(stratum_col):
        pools = defaultdict(list)
        for g, s in zip(bg.index, bg[stratum_col]):
            pools[s].append(gi[g])
        return {k: np.asarray(v) for k, v in pools.items()}
    pool_levels = [make_pools(c) for c in
                   ["stratum_full", "stratum_relax1", "stratum_relax2", "stratum_relax3", "stratum_relax4"]]
    all_pool = np.arange(len(bg_ids))
    gene_pool = {}
    n_fallback = 0
    for g in fg_list:
        pool = None
        for lvl, pools in zip(["stratum_full", "stratum_relax1", "stratum_relax2",
                               "stratum_relax3", "stratum_relax4"], pool_levels):
            s = bg.at[g, lvl]
            if s in pools and len(pools[s]) > 0:
                pool = pools[s]
                break
        if pool is None:
            n_fallback += 1
        gene_pool[g] = all_pool if pool is None else pool
    observed = np.array(hits, dtype=int)
    exceed = np.zeros(len(term_ids), dtype=np.int64)
    for start in range(0, N_PERM, CHUNK):
        c = min(CHUNK, N_PERM - start)
        idx = np.empty((c, len(fg_list)), dtype=np.int64)
        for j, g in enumerate(fg_list):
            pool = gene_pool[g]
            idx[:, j] = pool[rng.integers(0, len(pool), size=c)]
        rows = np.repeat(np.arange(c), len(fg_list))
        D = sparse.csr_matrix((np.ones(idx.size, dtype=np.int8), (rows, idx.ravel())),
                              shape=(c, len(bg_ids)))
        D.sum_duplicates()
        D.data[:] = 1
        H = (D @ A).toarray()
        exceed += (H >= observed[None, :]).sum(axis=0)
    p = (1.0 + exceed) / (N_PERM + 1.0)
    return p, n_fallback / max(1, len(fg_list))


if __name__ == "__main__":
    main()
