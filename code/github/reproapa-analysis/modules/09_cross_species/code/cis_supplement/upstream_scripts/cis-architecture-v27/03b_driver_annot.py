"""v27 (annotation-based) step 3: driver sites and model input on annotated PAS.

Detected study PAS usage is aggregated onto the annotated PAS they matched, then the
usage-gain / usage-loss driver sites are read off in the discovery species. Conservation of
the driver PAS uses the annotated reciprocal one-to-one map.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/root/nfsdata/apadetect/analysis_projects/reproapa_cis_architecture_v27")
IN = ROOT / "inputs"; OUT = ROOT / "outputs"; WORK = ROOT / "work"
HUMAN_CODES = ["H1", "H2", "H3", "H4"]; MOUSE_CODES = ["W", "C"]
CONTRASTS = {"P-S": ("P", "S"), "P-M": ("P", "M")}

usage = pd.read_csv(OUT / "site_usage_by_node.tsv.gz", sep="\t")
d2a = pd.read_csv(OUT / "detected_to_annot.tsv.gz", sep="\t")
d2a = d2a.dropna(subset=["annot_cut"])
d2a = d2a.rename(columns={"study": "code", "cut_site": "detected_cut"})
map_key = d2a[["species", "code", "gene_id", "chrom", "detected_cut", "strand", "annot_cut"]].drop_duplicates()
u = usage.merge(map_key, left_on=["species", "code", "gene_id", "chrom", "cut_site", "strand"],
                right_on=["species", "code", "gene_id", "chrom", "detected_cut", "strand"], how="inner")
au = u.groupby(["species", "code", "gene_id", "chrom", "annot_cut", "strand", "node"], as_index=False).agg(
    usage=("usage", "sum"), n_sites=("detected_cut", "size"), n_libraries=("n_libraries", "max"))
au = au.rename(columns={"annot_cut": "cut_site"})
au.to_csv(OUT / "site_usage_annot.tsv.gz", sep="\t", index=False, compression="gzip")

one = pd.read_csv(OUT / "pas_pairs_one2one_annot.tsv.gz", sep="\t")
h_cons = set(zip(one[one.reciprocal25].human_gene_id, one[one.reciprocal25].human_chrom,
                 one[one.reciprocal25].human_cut, one[one.reciprocal25].human_strand))
m_cons = set(zip(one[one.reciprocal25].mouse_gene_id, one[one.reciprocal25].mouse_chrom,
                 one[one.reciprocal25].mouse_cut, one[one.reciprocal25].mouse_strand))
orth = pd.read_csv(IN / "ortholog_universe.tsv", sep="\t",
                   usecols=["human_gene_id", "mouse_gene_id", "human_gene_symbol", "mouse_gene_symbol"])
m2h = orth.set_index("mouse_gene_id").human_gene_id

driver_rows = []
for contrast, (target, ref) in CONTRASTS.items():
    sub = au[au.node.isin([target, ref])]
    for species, codes in [("human", HUMAN_CODES), ("mouse", MOUSE_CODES)]:
        s = sub[sub.species.eq(species) & sub.code.isin(codes)]
        piv = s.pivot_table(index=["gene_id", "chrom", "cut_site", "strand", "code"],
                            columns="node", values="usage", aggfunc="first").reset_index()
        if target not in piv or ref not in piv:
            continue
        piv = piv.dropna(subset=[target, ref])
        piv["delta"] = piv[target] - piv[ref]
        agg = piv.groupby(["gene_id", "chrom", "cut_site", "strand"], as_index=False).agg(
            delta=("delta", "mean"), n_studies=("code", "nunique"),
            target_usage=(target, "mean"), ref_usage=(ref, "mean"))
        cons = h_cons if species == "human" else m_cons
        agg["is_conserved"] = [int((g, c, s_, st) in cons) for g, c, s_, st in
                               zip(agg.gene_id, agg.chrom, agg.cut_site, agg.strand)]
        for gene, g in agg.groupby("gene_id"):
            gain = g.loc[g.delta.idxmax()]; loss = g.loc[g.delta.idxmin()]
            driver_rows.append(dict(
                human_gene_id=gene if species == "human" else m2h.get(gene),
                species=species, contrast=contrast, gene_id=gene,
                n_discovery_pas=len(g), n_conserved_discovery_pas=int(g.is_conserved.sum()),
                frac_conserved_discovery_pas=float(g.is_conserved.mean()),
                gain_cut=int(gain.cut_site), gain_delta=float(gain.delta), gain_conserved=int(gain.is_conserved),
                gain_n_studies=int(gain.n_studies),
                loss_cut=int(loss.cut_site), loss_delta=float(loss.delta), loss_conserved=int(loss.is_conserved),
                loss_n_studies=int(loss.n_studies),
                driver_conserved=int(gain.is_conserved) + int(loss.is_conserved),
                driver_min_studies=int(min(gain.n_studies, loss.n_studies)),
                driver_same_site=int(gain.cut_site == loss.cut_site)))
driver = pd.DataFrame(driver_rows)
driver.to_csv(OUT / "driver_sites_annot.tsv.gz", sep="\t", index=False, compression="gzip")

uni = pd.read_csv(IN / "reciprocal_gene_universe.tsv.gz", sep="\t")
uni = uni[uni.contrast.isin(CONTRASTS) & uni.discovery.isin(["human", "mouse"])].copy()
arch = pd.read_csv(OUT / "gene_architecture_annot.tsv.gz", sep="\t")
arch = arch[["human_gene_id", "mouse_gene_id", "n_h_pas", "n_m_pas", "n_matched_pas", "A_g",
             "h_span", "m_span", "frac_h_conserved", "proximal_conserved", "distal_conserved",
             "both_anchors_conserved", "h_multi_source", "m_multi_source"]]
mi = uni.merge(driver.rename(columns={"species": "discovery"}),
               on=["human_gene_id", "contrast", "discovery"], how="left")
mi = mi.merge(arch, on="human_gene_id", how="left", suffixes=("", "_arch"))
mi["Y"] = mi.validation_same_direction.astype(int)
mi["abs_effect"] = mi.disc_meta.abs()
mi.to_csv(OUT / "model_input_annot.tsv.gz", sep="\t", index=False, compression="gzip")
qc = {"driver_rows": int(len(driver)), "model_input_rows": int(len(mi)),
      "driver_available": int(mi.driver_conserved.notna().sum()),
      "annot_site_rows": int(len(au)),
      "driver_contrast_species": {"|".join(map(str, k)): int(v) for k, v in driver.groupby(["contrast", "species"]).size().items()}}
(OUT / "driver_annot_qc.json").write_text(json.dumps(qc, indent=2) + "\n")
print(json.dumps(qc, indent=2))
