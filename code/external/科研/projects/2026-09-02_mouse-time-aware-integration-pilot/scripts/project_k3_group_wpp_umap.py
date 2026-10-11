"""Project Figure 2 C1-C3 mean WPP onto the frozen Figure 1 joint UMAP."""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.font_manager import FontProperties, fontManager
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.ndimage import gaussian_filter


WORK = Path(__file__).resolve().parents[3]
PROJECT = Path(__file__).resolve().parents[1]
FIG1 = WORK / "projects/2026-09-10_reproapa-figure1"
HUMAN = WORK / "projects/2026-09-09_human-time-aware-integration/data/results_metacell_pas_v2"
MOUSE = PROJECT / "data/metacell-pas-v7-v1"
HUMAN_SUPPORT = PROJECT / "data/apa-majority-shared-v1/human_annotation"
MOUSE_SUPPORT = PROJECT / "data/pas-annotation-support-v1/studies"
UMAP = WORK / "projects/2026-09-10_human-mouse-time-aware-umap/data/results_raw_cell_species_integrated_v3_faiss2/raw_cell_nodes.tsv.gz"
STATES = UMAP.parent / "state_and_trajectory_display_v13/raw_cell_stateflow_groups.tsv.gz"
FIG1_LABELS = UMAP.parent / "state_and_trajectory_display_v26/state_label_positions.tsv"
MEMBERSHIP = PROJECT / "outputs/figure2-k3-srelative-v2/gene_membership_k3.tsv.gz"
OUT = PROJECT / "outputs/figure2-k3-group-umap-v1"
FONT = Path("/Users/lycidas/Library/Fonts/nimbus-sans-l.regular.otf")
CLASSES = ("C1", "C2", "C3")
FEATURE_COLUMNS = ["feature_index_1based", "atlas_pas_id", "chrom", "cut_site", "strand",
                   "primary_gene_id", "keep_for_apa_primary"]


def study_values(source, support, species, membership):
    meta = pd.read_csv(source / "metacells.tsv.gz", sep="\t", usecols=["metacell_id"])
    features = pd.read_csv(source / "features.tsv.gz", sep="\t", usecols=FEATURE_COLUMNS,
                           low_memory=False)
    assert np.array_equal(features.feature_index_1based.to_numpy(),
                          np.arange(1, len(features) + 1))
    features = features.loc[features.keep_for_apa_primary &
                            features.primary_gene_id.isin(membership.gene_id)].copy()
    features = features.sort_values(["primary_gene_id", "cut_site", "atlas_pas_id"]).reset_index(drop=True)
    reference = pd.read_csv(support / source.name / "primary_pas_reference_matches.tsv.gz",
                            sep="\t", usecols=["atlas_pas_id", "matched_source_count_50nt"])
    features = features.merge(reference, on="atlas_pas_id", validate="one_to_one")
    assert features.matched_source_count_50nt.notna().all()
    genes = pd.Index(features.primary_gene_id.unique())
    gene_index = genes.get_indexer(features.primary_gene_id)
    starts = np.r_[0, np.flatnonzero(np.diff(gene_index)) + 1]
    strand = features.strand.iloc[starts].to_numpy()
    signed_position = features.cut_site.to_numpy(np.float64) * np.where(features.strand.eq("+"), 1., -1.)
    low = np.minimum.reduceat(signed_position, starts)
    high = np.maximum.reduceat(signed_position, starts)
    matched = features.matched_source_count_50nt.ge(1).to_numpy()
    valid = (features.groupby("primary_gene_id", sort=False).chrom.nunique().to_numpy() == 1)
    valid &= (features.groupby("primary_gene_id", sort=False).strand.nunique().to_numpy() == 1)
    valid &= np.isin(strand, ["+", "-"]) & (high > low)
    valid &= np.add.reduceat(matched.astype(np.int32), starts) >= 2
    relative = np.divide(signed_position - low[gene_index], (high - low)[gene_index],
                         out=np.zeros(len(features)), where=(high - low)[gene_index] > 0)

    counts = sparse.load_npz(source / "metacell_by_pas_counts.npz").tocsr()
    assert len(meta) == counts.shape[0] and features.feature_index_1based.max() <= counts.shape[1]
    counts = counts[:, features.feature_index_1based.to_numpy(np.int64) - 1].multiply(matched).tocsr()
    mapping = sparse.csr_matrix((np.ones(len(features), np.float32),
                                 (np.arange(len(features)), gene_index)),
                                shape=(len(features), len(genes)))
    weighting = sparse.csr_matrix((relative.astype(np.float32),
                                   (np.arange(len(features)), gene_index)),
                                  shape=mapping.shape)
    totals = (counts @ mapping).toarray()
    numerators = (counts @ weighting).toarray()
    observable = (totals >= 20) & valid
    wpp = np.divide(numerators, totals, out=np.full(totals.shape, np.nan, np.float32),
                    where=observable)
    direct_checks = 0
    for row in np.linspace(0, len(meta) - 1, min(6, len(meta)), dtype=int):
        for gene in np.flatnonzero(observable[row])[:2]:
            stop = starts[gene + 1] if gene + 1 < len(starts) else len(features)
            pas_counts = counts[row, starts[gene]:stop].toarray().ravel()
            direct = np.dot(pas_counts, relative[starts[gene]:stop]) / pas_counts.sum()
            assert np.isclose(direct, wpp[row, gene], atol=2e-5)
            direct_checks += 1
    gene_mean = np.divide(np.nansum(wpp, axis=0), observable.sum(axis=0),
                          out=np.full(len(genes), np.nan), where=observable.sum(axis=0) > 0)
    cluster = membership.set_index("gene_id").cluster.reindex(genes).to_numpy()
    result = pd.DataFrame({"species": species, "study": source.name,
                           "metacell_id": meta.metacell_id})
    audit = {"species": species, "study": source.name, "metacells": len(meta),
             "direct_pas_checks": direct_checks}
    for klass in CLASSES:
        take = cluster == klass
        count = observable[:, take].sum(axis=1)
        raw = np.divide(np.nansum(wpp[:, take], axis=1), count,
                        out=np.full(len(meta), np.nan), where=count > 0)
        centered = np.divide(np.nansum(wpp[:, take] - gene_mean[take], axis=1), count,
                             out=np.full(len(meta), np.nan), where=count > 0)
        result[f"{klass}_wpp"] = raw
        result[f"{klass}_centered"] = centered
        result[f"{klass}_n"] = count
        audit[klass] = {"candidate_genes": int((membership.cluster == klass).sum()),
                        "valid_genes": int(valid[take].sum()),
                        "metacells_with_any": int((count > 0).sum()),
                        "median_evaluable_genes": float(np.median(count)),
                        "p10_evaluable_genes": float(np.quantile(count, .1)),
                        "min_evaluable_genes": int(count.min())}
    return result, audit


def attach_cells(metacells):
    maps = []
    for species, root in (("human", HUMAN), ("mouse", MOUSE)):
        part = pd.read_csv(root / "all_cell_to_pas_and_metacell.tsv.gz", sep="\t",
                           usecols=["cell_id", "metacell_id"])
        part["species"] = species
        maps.append(part)
    mapping = pd.concat(maps, ignore_index=True)
    cells = pd.read_csv(UMAP, sep="\t",
                        usecols=["cell_id", "species_key", "study_id", "UMAP1", "UMAP2"])
    cells = cells.merge(mapping, on="cell_id", validate="one_to_one")
    cells = cells.merge(metacells, on=["species", "metacell_id"], validate="many_to_one")
    assert len(cells) == 551201 and cells.cell_id.is_unique
    assert cells.species.eq(cells.species_key).all()
    assert cells.study_id.replace({"WOT_GSE115943": "GSE115943"}).str.lower().eq(
        cells.study.str.split("_").str[0]).all()
    return cells


def state_positions(cells):
    states = pd.read_csv(STATES, sep="\t", usecols=["cell_id", "state_group"])
    nodes = cells[["cell_id", "species", "UMAP1", "UMAP2"]].merge(
        states, on="cell_id", validate="one_to_one")
    codes = {"shared_s": "S", "shared_r": "R", "shared_e": "E",
             "shared_tp": "P", "shared_m": "M"}
    reference_labels = pd.read_csv(FIG1_LABELS, sep="\t").set_index("state_group")
    result = {}
    for species in ("mouse", "human"):
        result[species] = []
        for group, label in codes.items():
            part = nodes.loc[nodes.species.eq(species) & nodes.state_group.eq(group)]
            assert len(part) > 0
            hist, xedges, yedges = np.histogram2d(part.UMAP1, part.UMAP2, bins=100)
            peak = np.unravel_index(np.argmax(gaussian_filter(hist, 1.5)), hist.shape)
            x = (xedges[peak[0]] + xedges[peak[0] + 1]) / 2
            y = (yedges[peak[1]] + yedges[peak[1] + 1]) / 2
            # Human R has several UMAP islands. Its absolute density maximum is a
            # small lower island; use the frozen Figure 1 R landmark instead.
            if species == "human" and group == "shared_r":
                x = float(reference_labels.loc[group, "anchor_x"])
                y = float(reference_labels.loc[group, "anchor_y"])
            # The human M peak over-emphasizes one Guan2022 island. The median
            # lies within the main M region and is locally state-pure.
            if species == "human" and group == "shared_m":
                x = float(part.UMAP1.median())
                y = float(part.UMAP2.median())
            result[species].append((label, x, y))
    return result


def render(cells, landmarks, species, metric, filename):
    fontManager.addfont(str(FONT))
    family = FontProperties(fname=str(FONT)).get_name()
    plt.rcParams.update({"font.family": family, "font.size": 8, "pdf.fonttype": 42,
                         "svg.fonttype": "none", "axes.linewidth": .6})
    subset = cells.loc[cells.species.eq(species)]
    values = np.concatenate([subset[f"{klass}_{metric}"].dropna().to_numpy() for klass in CLASSES])
    if metric == "centered":
        limit = float(np.quantile(np.abs(values), .99))
        norm = colors.TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
        cmap = "coolwarm"
        key = "Mean WPP relative to each study/gene mean"
        ticks = [-limit, 0, limit]
    else:
        lower, upper = np.quantile(values, [.01, .99])
        norm = colors.Normalize(vmin=float(lower), vmax=float(upper))
        cmap = "viridis"
        key = "Mean weighted PAS position (WPP)"
        ticks = [float(lower), float((lower + upper) / 2), float(upper)]
    xy = cells[["UMAP1", "UMAP2"]].to_numpy(np.float32)
    low = xy.min(axis=0) - .035 * np.ptp(xy, axis=0)
    high = xy.max(axis=0) + .035 * np.ptp(xy, axis=0)
    back = np.random.default_rng(260924).permutation(len(cells))
    fig, axes = plt.subplots(1, 3, figsize=(160 / 25.4, 73 / 25.4))
    fig.subplots_adjust(left=.025, right=.975, top=.88, bottom=.24, wspace=.02)
    artist = None
    qa = []
    for klass, ax in zip(CLASSES, axes):
        column = f"{klass}_{metric}"
        valid = subset.loc[subset[column].notna()]
        order = np.random.default_rng(260924 + int(klass[1])).permutation(len(valid))
        ax.scatter(xy[back, 0], xy[back, 1], s=.12, color="#A9B3BB", alpha=.075,
                   linewidths=0, rasterized=True)
        artist = ax.scatter(valid.UMAP1.to_numpy()[order], valid.UMAP2.to_numpy()[order],
                            c=valid[column].to_numpy()[order], s=.30, cmap=cmap, norm=norm,
                            alpha=.87, linewidths=0, rasterized=True)
        ax.set_xlim(low[0], high[0]); ax.set_ylim(low[1], high[1])
        ax.set_aspect("equal"); ax.axis("off")
        for label, x, y in landmarks[species]:
            ax.text(x, y, label, ha="center", va="center", fontsize=8,
                    fontweight="bold", color="#253947", zorder=8,
                    bbox=dict(boxstyle="round,pad=.13", facecolor="white",
                              edgecolor="#546774", linewidth=.55, alpha=.96))
        ax.text(.03, .97, klass, transform=ax.transAxes, ha="left", va="top",
                fontsize=10, fontweight="bold", color="#253947")
        qa.append({"class": klass, "cells_colored": len(valid),
                   "metacells_colored": int(valid.metacell_id.nunique()),
                   "cell_median_evaluable_genes": float(valid[f"{klass}_n"].median()),
                   "cell_value_median": float(valid[column].median())})
    bar = fig.add_axes([.27, .14, .46, .035])
    cb = fig.colorbar(artist, cax=bar, orientation="horizontal", ticks=ticks)
    cb.ax.tick_params(labelsize=8, direction="in", length=3, pad=2)
    cb.ax.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    cb.outline.set_linewidth(.6)
    fig.text(.5, .195, key, ha="center", fontsize=8, color="#253947")
    fig.text(.025, .945, species.capitalize(), fontsize=10, color="#253947")
    fig.savefig(OUT / f"{filename}.pdf", dpi=300, facecolor="white")
    fig.savefig(OUT / f"{filename}.png", dpi=300, facecolor="white")
    plt.close(fig)
    return {"species": species, "metric": metric, "color_range": [norm.vmin, norm.vmax],
            "classes": qa}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    membership = pd.read_csv(MEMBERSHIP, sep="\t", usecols=["species", "gene_id", "cluster"])
    assert membership.gene_id.notna().all()
    assert not membership.duplicated(["species", "gene_id"]).any()
    assert set(membership.cluster) == set(CLASSES)
    rows, study_audit = [], []
    for species, root, support in (("human", HUMAN, HUMAN_SUPPORT),
                                   ("mouse", MOUSE, MOUSE_SUPPORT)):
        members = membership.loc[membership.species.eq(species)]
        for source in sorted((root / "studies").iterdir()):
            frame, audit = study_values(source, support, species, members)
            rows.append(frame); study_audit.append(audit)
            print(species, source.name, {k: audit[k]["median_evaluable_genes"] for k in CLASSES}, flush=True)
    metacells = pd.concat(rows, ignore_index=True)
    assert len(metacells) == 11104 and metacells.metacell_id.is_unique
    metacells.to_csv(OUT / "k3_group_wpp_metacells.tsv.gz", sep="\t", index=False,
                     compression="gzip")
    cells = attach_cells(metacells)
    landmarks = state_positions(cells)
    figures = []
    for species in ("mouse", "human"):
        for metric in ("wpp", "centered"):
            figures.append(render(cells, landmarks, species, metric,
                                  f"{species}_C1-C3_mean_{metric}_joint_umap"))
    qa = {"status": "AWAITING_VISUAL_QA", "umap_cells": len(cells),
          "metacells": len(metacells), "study_audit": study_audit, "figures": figures,
          "landmarks": landmarks,
          "source_membership": str(MEMBERSHIP), "source_umap": str(UMAP),
          "interpretation": "Metacell-level group mean WPP assigned to member cells; centered values remove each study/gene metacell mean. No per-cell PAS estimate."}
    (OUT / "qa.json").write_text(json.dumps(qa, indent=2) + "\n")
    print(json.dumps({k: v for k, v in qa.items() if k != "study_audit"}, indent=2))


if __name__ == "__main__":
    main()
