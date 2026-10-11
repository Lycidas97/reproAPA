"""Verify full predicted intervals and coordinate reciprocity, not binding."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd


BASE = Path(__file__).resolve().parent
SITE = BASE / "sites"
KEYS = ["source_population", "tolerance", "contrast", "human_gene_id"]


def main():
    seal = json.loads((SITE / "mirna_sealed_manifest.json").read_text())
    assert seal["status"] == "SEALED_MIRNA_FOR_INDEPENDENT_REVIEW"
    for name, expected in seal["files"].items():
        assert hashlib.sha256((SITE / name).read_bytes()).hexdigest() == expected, name
    seeds = pd.read_csv(SITE / "strict_case_predicted_seeds.tsv", sep="\t")
    pairs = pd.read_csv(SITE / "reciprocal_predicted_seed_pairs.tsv", sep="\t")
    summary = pd.read_csv(SITE / "strict_case_seed_summary.tsv", sep="\t")
    intervals = pd.read_csv(SITE / "strict_case_intervals.tsv", sep="\t")
    traces = pd.read_csv(SITE / "strict_case_interval_inclusion.tsv", sep="\t")
    queries = pd.read_csv(SITE / "seed_mapping_queries.tsv", sep="\t").set_index("query_id")
    calibration = pd.read_csv(SITE / "chain_calibration_queries.tsv", sep="\t")
    mappings = {x["query_id"]: x for x in map(json.loads, (SITE / "seed_chain_maps.jsonl").read_text().splitlines())}
    assert set(queries.index) | set(calibration.query_id) == set(mappings)
    for row in calibration.itertuples(index=False):
        mapped = mappings[row.query_id]
        if row.expected_status == "unique":
            assert mapped["status"] == "unique" and len(mapped["hits"]) == 1
            hit = mapped["hits"][0]
            cut = hit["end0"] if hit["strand"] == "+" else hit["start0"]
            assert (hit["chrom"], hit["strand"], cut) == (row.expected_chrom, row.expected_strand, row.expected_cut)
        else:
            assert row.expected_status == "unmapped" and not mapped["hits"]
    assert seeds.end0.sub(seeds.start0).isin([7, 8]).all()
    assert seeds.start0.ge(seeds.core_start0).all() and seeds.end0.le(seeds.core_end0).all()
    assert not pairs.duplicated(KEYS + ["human_query_id", "mouse_query_id"]).any()
    for pair in pairs.itertuples(index=False):
        h, m = queries.loc[pair.human_query_id], queries.loc[pair.mouse_query_id]
        assert h.species == "human" and m.species == "mouse" and h.seed == m.seed == pair.seed
        assert h.end0 - h.start0 == m.end0 - m.start0
        for source, target in [(pair.human_query_id, m), (pair.mouse_query_id, h)]:
            hits = mappings[source]["hits"]
            assert len(hits) == 1 and mappings[source]["status"] == "unique"
            assert (hits[0]["chrom"], hits[0]["strand"], hits[0]["start0"], hits[0]["end0"]) == (target.chrom, target.strand, target.start0, target.end0)
        for query_id in [pair.human_query_id, pair.mouse_query_id]:
            occurrence = seeds[(seeds.source_population == pair.source_population) & (seeds.tolerance == pair.tolerance) & (seeds.contrast == pair.contrast) & (seeds.human_gene_id == pair.human_gene_id) & (seeds.query_id == query_id)]
            assert len(occurrence) == 1
    for row in summary.to_dict("records"):
        mask = pd.Series(True, index=seeds.index)
        pmask = pd.Series(True, index=pairs.index)
        for key in KEYS:
            mask &= seeds[key].eq(row[key])
            pmask &= pairs[key].eq(row[key])
        group, paired = seeds[mask], pairs[pmask]
        h = group[group.species.eq("human")]
        m = group[group.species.eq("mouse")]
        assert len(h) == row["human_predicted_physical_sites"] and len(m) == row["mouse_predicted_physical_sites"]
        assert len(set(h.seed) & set(m.seed)) == row["same_seed_groups_both_species"]
        assert len(paired) == row["reciprocal_physical_site_pairs"]
        assert paired.seed.nunique() == row["reciprocal_seed_groups"]
    group_keys = KEYS + ["species"]
    reconstructed = traces.groupby(group_keys + ["study_id"])[["target_inclusion", "reference_inclusion", "delta_inclusion"]].mean().groupby(group_keys).mean().sort_index()
    recorded = intervals.set_index(group_keys)[reconstructed.columns].sort_index()
    np.testing.assert_allclose(reconstructed, recorded, atol=1e-12)
    np.testing.assert_allclose(traces.target_inclusion - traces.reference_inclusion, traces.delta_inclusion, atol=1e-12)
    outcome = {
        "status": "PASS_MIRNA_SUMMARIES_AND_RECIPROCAL_OUTPUT_COORDINATES",
        "sealed_files_verified": len(seal["files"]),
        "mapped_query_intervals": len(queries),
        "PAS_point_calibrations_independently_checked": len(calibration),
        "gene_context_summaries": len(summary),
        "reciprocal_pair_rows_across_tolerance_contexts": len(pairs),
        "unique_physical_pairs_without_tolerance_duplication": len(pairs.drop_duplicates(["human_gene_id", "human_query_id", "mouse_query_id"])),
        "binding_or_functional_inference": False,
        "scope": "All saved coordinate reciprocities and summaries; chain-block source algorithm and full native inclusion audit are in review/.",
    }
    (BASE / "root_mirna_validation.json").write_text(json.dumps(outcome, indent=2) + "\n")
    print(json.dumps(outcome))


if __name__ == "__main__":
    main()
