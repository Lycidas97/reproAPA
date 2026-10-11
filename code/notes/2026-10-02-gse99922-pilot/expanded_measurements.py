"""Expand verified gene identities without changing frozen members or PAS rules."""
from collections import defaultdict
from pathlib import Path
import json
import subprocess

import numpy as np
import pandas as pd

from matched_analysis import digest, read

HERE = Path(__file__).resolve().parent
OUT = HERE / "vector_results"
BASE = Path("/root/nfsdata/apadetect/analysis_projects/reproapa_pas_reliability_2026_09_13/data")
STUDIES = ["gse115943_pluripotency_reprogramming_core",
           "gse103221_pluripotency_reprogramming_core",
           "gse171127_pluripotency_vs_trophoblast_bifurcation"]
ENDPOINT = ["MEF_R1", "MEF_R2", "iPSC_R1", "iPSC_R2"]


def save(frame, name):
    frame.to_csv(OUT / name, sep="\t", index=False, float_format="%.12g")


def resolve_gene(entrez, symbol, xrefs, reciprocal, geometry, symbols):
    """Authoritative ambiguity cannot be rescued by a matching symbol."""
    ensembl = xrefs.get(entrez, set())
    if ensembl:
        if len(ensembl) != 1:
            return "", "rejected", "ambiguous_ncbi_xref"
        gene = next(iter(ensembl))
        if len(reciprocal[gene]) != 1:
            return "", "rejected", "nonreciprocal_ncbi_xref"
        if gene not in geometry:
            return "", "rejected", "no_unique_reference_geometry"
        return gene, "reciprocal_ncbi_geneid", "mapped"
    if symbol in symbols:
        return symbols[symbol], "exact_symbol_fallback", "mapped"
    return "", "rejected", "no_verified_identity"


def lift(sites, offset):
    stem = OUT / f"expanded_offset{offset}"
    pd.DataFrame({"chrom": sites.chrom, "start": sites.pos - offset,
        "end": sites.pos - offset + 1, "name": sites.COORDS,
        "score": 0, "strand": sites.source_strand}).to_csv(
            str(stem) + ".mm9.bed", sep="\t", index=False, header=False)
    for source, target, chain in [("mm9", "mm10", "mm9ToMm10.over.chain.gz"),
                                  ("mm10", "mm39", "mm10ToMm39.over.chain.gz")]:
        subprocess.run([str(HERE / "tools/liftOver"), str(stem) + f".{source}.bed",
            str(HERE / "inputs" / chain), str(stem) + f".{target}.bed",
            str(stem) + f".{source}.unmapped.bed"], check=True,
            stdout=subprocess.DEVNULL)
    lifted = read(Path(str(stem) + ".mm39.bed"), header=None,
        names=["mapped_chrom", "start", "end", "COORDS", "score", "orientation"])
    assert lifted.COORDS.is_unique
    lifted["mapped_cut"] = lifted.start + 1
    return lifted[["COORDS", "mapped_chrom", "mapped_cut", "orientation"]]


def profile_metrics(proportions, local_cut, source=(0, 1)):
    target = [i for i in range(4) if i not in source]
    change = proportions[target].mean(axis=0) - proportions[list(source)].mean(axis=0)
    delta_bp = change @ local_cut
    span = np.ptp(local_cut)
    return np.array([np.abs(change).sum() / 2, np.ptp(change),
                     delta_bp / span if span else np.nan, delta_bp])


def run():
    OUT.mkdir(exist_ok=True)
    protocol = HERE / "vector-protocol.json"
    hashes = {"vector-protocol.json": digest(protocol),
        "inputs/frozen_candidates.tsv": digest(HERE / "inputs/frozen_candidates.tsv"),
        "inputs/GSE99922_readcounts-table.txt.gz": digest(HERE / "inputs/GSE99922_readcounts-table.txt.gz"),
        "inputs/Mus_musculus.gene_info.gz": digest(HERE / "inputs/Mus_musculus.gene_info.gz")}
    frames = []
    for study in STUDIES:
        path = BASE / study / "primary_features.tsv.gz"
        hashes[str(path)] = digest(path)
        frames.append(read(path, usecols=["atlas_pas_id", "genome_build", "chrom",
            "cut_site", "strand", "primary_gene_id", "primary_gene_symbol"]))
    ref = pd.concat(frames).drop_duplicates().reset_index(drop=True)
    assert ref.genome_build.eq("GRCm39_2024_A").all() and ref.atlas_pas_id.is_unique
    geometry_frame = ref[["primary_gene_id", "primary_gene_symbol", "chrom", "strand"]].drop_duplicates()
    unique_gene = geometry_frame.groupby("primary_gene_id").size().eq(1)
    geometry_frame = geometry_frame[geometry_frame.primary_gene_id.isin(unique_gene[unique_gene].index)]
    geometry = geometry_frame.set_index("primary_gene_id").to_dict("index")
    unique_symbol = geometry_frame.groupby("primary_gene_symbol").size().eq(1)
    symbols = geometry_frame[geometry_frame.primary_gene_symbol.isin(unique_symbol[unique_symbol].index)].set_index("primary_gene_symbol").primary_gene_id.to_dict()
    info = read(HERE / "inputs/Mus_musculus.gene_info.gz", dtype=str, keep_default_na=False)
    info = info[info["#tax_id"].eq("10090")]
    xrefs, reciprocal = defaultdict(set), defaultdict(set)
    for row in info.itertuples(index=False):
        for token in row.dbXrefs.split("|"):
            if token.startswith("Ensembl:ENSMUSG"):
                gene = token.split(":", 1)[1]
                xrefs[row.GeneID].add(gene)
                reciprocal[gene].add(row.GeneID)
    raw = read(HERE / "inputs/GSE99922_readcounts-table.txt.gz", dtype={"IDs": str})
    count_columns = list(raw.columns[5:])
    duplicates = raw[raw.COORDS.duplicated(False)]
    assert duplicates.groupby("COORDS")[count_columns].nunique().le(1).all().all()
    assert raw[count_columns].ge(0).all().all()
    assert np.equal(raw[count_columns], np.floor(raw[count_columns])).all().all()
    identities = []
    for row in raw[["IDs", "GENES"]].drop_duplicates().itertuples(index=False):
        gene, method, status = resolve_gene(row.IDs, row.GENES, xrefs, reciprocal, geometry, symbols)
        identities.append({"IDs": row.IDs, "GENES": row.GENES, "gene_id": gene,
            "identity_method": method, "identity_status": status})
    identity = pd.DataFrame(identities)
    save(identity, "expanded_id_map.tsv.gz")
    # Keep the prior conservative exclusions, even if aliases might be resolvable.
    bad = raw.groupby("COORDS").GENES.nunique().gt(1)
    sites = raw[~raw.COORDS.isin(bad[bad].index)]
    type_conflict = sites.groupby("COORDS").TYPES.nunique().gt(1)
    sites = sites[~sites.COORDS.isin(type_conflict[type_conflict].index)].drop_duplicates("COORDS")
    sites = sites[sites.TYPES.eq("LastExon")].merge(identity, on=["IDs", "GENES"], validate="many_to_one")
    sites = sites[sites.identity_status.eq("mapped")].copy()
    sites[["chrom", "source_strand", "pos"]] = sites.COORDS.str.extract(r"^(chr.+)([+-])(\d+)$")
    assert sites.pos.notna().all() and sites.COORDS.is_unique
    sites["pos"] = sites.pos.astype(int)
    sites["ref_chrom"] = sites.gene_id.map(lambda g: geometry[g]["chrom"])
    sites["ref_strand"] = sites.gene_id.map(lambda g: geometry[g]["strand"])
    sites["gene_symbol"] = sites.gene_id.map(lambda g: geometry[g]["primary_gene_symbol"])
    mapped = sites.merge(lift(sites, 1), on="COORDS", how="left", validate="one_to_one")
    alternate = lift(sites, 0).rename(columns={"mapped_cut": "alternate_cut", "mapped_chrom": "alternate_chrom", "orientation": "alternate_orientation"})
    mapped = mapped.merge(alternate, on="COORDS", how="left", validate="one_to_one")
    mapped["coordinate_base_sensitive"] = (mapped.mapped_chrom.ne(mapped.alternate_chrom)
        | mapped.orientation.ne(mapped.alternate_orientation)
        | ((mapped.mapped_cut - mapped.alternate_cut).abs() > 1)
        | mapped.mapped_cut.isna() | mapped.alternate_cut.isna())
    by_gene = {g: v for g, v in ref.groupby("primary_gene_id")}
    matches = []
    for row in mapped.itertuples():
        result = {"COORDS": row.COORDS, "atlas_pas_id": "", "distance_nt": np.nan,
                  "compatible50": False, "nearest_tie": False}
        if (row.mapped_chrom == row.ref_chrom and row.orientation == row.ref_strand
            and row.orientation == row.source_strand and not row.coordinate_base_sensitive):
            r = by_gene[row.gene_id]
            d = (r.cut_site - row.mapped_cut).abs()
            other = (r.cut_site - row.alternate_cut).abs()
            result.update(distance_nt=d.min(), nearest_tie=int(d.eq(d.min()).sum()) != 1)
            if (d.min() <= 50 and not result["nearest_tie"] and other.min() <= 50
                and other.eq(other.min()).sum() == 1 and d.idxmin() == other.idxmin()):
                result.update(atlas_pas_id=r.loc[d.idxmin(), "atlas_pas_id"], compatible50=True)
        matches.append(result)
    mapped = mapped.merge(pd.DataFrame(matches), on="COORDS", validate="one_to_one")
    mapped["many_to_one_collision"] = mapped.compatible50 & mapped.atlas_pas_id.duplicated(False)
    mapped["usable"] = mapped.compatible50 & ~mapped.many_to_one_collision
    save(mapped, "expanded_site_mapping.tsv.gz")
    depths = raw.drop_duplicates("COORDS")[ENDPOINT].sum()
    rows, vectors = [], []
    for gene, all_sites in mapped.groupby("gene_id"):
        shared = all_sites[all_sites.usable].copy()
        signed = shared.mapped_cut.to_numpy() * (1 if all_sites.ref_strand.iloc[0] == "+" else -1)
        local = signed - signed.min() if len(signed) else np.array([])
        counts = shared[ENDPOINT].to_numpy().T
        totals = counts.sum(axis=1)
        public_totals = all_sites[ENDPOINT].sum().to_numpy()
        proportions = np.divide(counts, totals[:, None], out=np.full(counts.shape, np.nan), where=totals[:, None] > 0)
        fractions = np.divide(totals, public_totals, out=np.zeros(4), where=public_totals > 0)
        overlapping = len(shared) > 1 and shared.pos.sort_values().diff().min() <= 80
        span = np.ptp(local) if len(local) else 0
        support = int((counts.sum(axis=0) >= 10).sum())
        metrics = profile_metrics(proportions, local) if len(shared) else np.full(4, np.nan)
        pooled = proportions.mean(axis=0)
        record = {"gene_id": gene, "gene_symbol": all_sites.gene_symbol.iloc[0],
            "identity_methods": ";".join(sorted(all_sites.identity_method.unique())),
            "n_external_last_exon_sites": len(all_sites), "n_shared_sites": len(shared),
            "n_collision_sites": int(all_sites.many_to_one_collision.sum()),
            "overlapping_shared_windows": bool(overlapping), "shared_span_bp": span,
            "supported_sites": support, "min_shared_count": totals.min(),
            "tvd": metrics[0], "mpro": metrics[1], "delta_shared_axis": metrics[2], "delta_bp": metrics[3],
            "mean_log2_last_exon_CPM": float(np.log2(public_totals / depths.to_numpy() * 1e6 + .5).mean()),
            "log2_min_shared_count": float(np.log2(max(totals.min(), 1))),
            "log2_external_last_exon_sites": float(np.log2(len(all_sites))),
            "min_shared_fraction": fractions.min(), "log2_shared_span_bp": float(np.log2(max(span, 1))),
            "pooled_shared_axis_usage": (pooled @ local) / span if span else np.nan,
            "pooled_proportion_impurity": 1 - np.square(pooled).sum(),
            "pas_stratum": "2" if len(shared) == 2 else ("3-4" if len(shared) <= 4 else "5+")}
        for threshold in [20, 50, 100]:
            record[f"measurable{threshold}"] = (len(shared) >= 2 and support >= 2 and span >= 100
                and not overlapping and totals.min() >= threshold and fractions.min() >= .8)
        for i, sample in enumerate(ENDPOINT):
            record[f"shared_count__{sample}"] = totals[i]
            record[f"shared_fraction__{sample}"] = fractions[i]
            shared[f"usage__{sample}"] = proportions[i]
        shared["local_cut_bp"] = local
        vectors.append(shared)
        rows.append(record)
    genes = pd.DataFrame(rows)
    save(genes, "all_gene_measurements.tsv.gz")
    save(pd.concat(vectors, ignore_index=True), "shared_pas_vectors.tsv.gz")
    frozen = read(HERE / "inputs/frozen_candidates.tsv")
    assert frozen.groupby("set_name").size().to_dict() == {"primary_fixed130": 130, "secondary_PS488": 488}
    candidates = frozen.merge(genes, on="gene_id", how="left", suffixes=("", "_external"), validate="many_to_one")
    old = read(HERE / "results/frozen_candidate_results.tsv")
    old = old[["set_name", "gene_id", "measurable20__iPSC"]].rename(columns={"measurable20__iPSC": "prior_measurable20"})
    candidates = candidates.merge(old, on=["set_name", "gene_id"], validate="one_to_one")
    save(candidates, "expanded_candidate_results.tsv")
    coverage = []
    for name, group in candidates.groupby("set_name"):
        for threshold in [20, 50, 100]:
            measured = group[f"measurable{threshold}"].eq(True)
            prior = group.prior_measurable20.eq(True)
            coverage.append({"set_name": name, "threshold": threshold, "frozen_n": len(group),
                "identity_mapped_n": int(group.n_shared_sites.notna().sum()),
                "two_shared_sites_n": int(group.n_shared_sites.ge(2).sum()), "measurable_n": int(measured.sum()),
                "newly_measurable_vs_prior_count20": int((measured & ~prior).sum()),
                "previous_count20_missing_here": int((~measured & prior).sum())})
    save(pd.DataFrame(coverage), "expanded_coverage.tsv")
    record = {"status": "completed", "script_sha256": digest(Path(__file__)), "inputs_sha256": hashes,
        "public_gene_identity_statuses": identity.identity_status.value_counts().to_dict(),
        "mapping_methods": identity.identity_method.value_counts().to_dict(),
        "mapped_gene_count": len(genes), "measurable20_gene_count": int(genes.measurable20.sum()),
        "scope": "GeneID expanded measurement, not replacement of parent signed tests"}
    (OUT / "measurement_run.json").write_text(json.dumps(record, indent=2) + "\n")
    print(pd.DataFrame(coverage).to_string(index=False))


if __name__ == "__main__":
    run()
