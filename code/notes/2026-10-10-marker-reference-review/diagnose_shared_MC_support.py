"""Diagnose source-MC double counting in intersecting cross-reference cohorts.

No classification, group effects or direction is changed. Only the existing
joint PAS qualification is compared with a distinct source-MC union.
"""
import csv
import gzip
import json

import numpy as np
import pandas as pd

from verify_APA_delivery import HERE, SOURCE, independent_cohorts, read
from verify_marker_delivery import STUDIES, records, sha


FIELDS = ["study_id", "route", "library_key", "setting", "version", "contrast"]


def joint_counts(counts, a, b):
    ma, mb = counts[a].sum(axis=0), counts[b].sum(axis=0)
    na, nb = (counts[a] > 0).sum(axis=0), (counts[b] > 0).sum(axis=0)
    shared = np.intersect1d(a, b)
    molecules = ma + mb
    nonzero = na + nb
    distinct_molecules = molecules - counts[shared].sum(axis=0)
    distinct_nonzero = nonzero - (counts[shared] > 0).sum(axis=0)
    old_joint = (molecules >= 10) & (nonzero >= 2)
    unique_joint = (distinct_molecules >= 10) & (distinct_nonzero >= 2)
    assert np.all(~unique_joint | old_joint)
    return shared, old_joint, unique_joint


def verify_intent():
    # The same MC alone expressing a distal PAS must not count as two sources,
    # even when both compared groups independently contain three MC.
    counts = np.array([[10, 10], [20, 0], [20, 0], [20, 0], [20, 0]])
    shared, old, distinct = joint_counts(counts, np.array([0, 1, 2]), np.array([0, 3, 4]))
    assert shared.tolist() == [0] and old.tolist() == [True, True]
    assert distinct.tolist() == [True, False]
    # Disjoint original P/M and mutually exclusive P*/M* keep the original gate.
    shared, old, distinct = joint_counts(counts, np.array([0, 1]), np.array([3, 4]))
    assert not len(shared) and np.array_equal(old, distinct)


def main():
    verify_intent()
    z = np.load(SOURCE / "inputs.npz", allow_pickle=False)
    expected = pd.read_csv(HERE / "expected_marker_classification.tsv.gz", sep="\t", float_precision="round_trip")
    groups, metas = independent_cohorts(z, expected)
    eligibility = read("library_apa_eligibility.tsv")
    contexts = {}
    for row in eligibility[eligibility.status.eq("pre_gene_eligible")].to_dict("records"):
        lib = (row["study_id"], row["route"], row["library_key"])
        a, b = groups[lib + (row["target_cohort"],)], groups[lib + (row["reference_cohort"],)]
        shared = np.intersect1d(a, b)
        if len(shared):
            key = tuple(row[f] for f in FIELDS)
            contexts[key] = dict(row, a=a, b=b, shared=shared)
    print("CANDIDATE_OVERLAP_CONTEXTS", len(contexts), flush=True)
    delivered = {}
    with gzip.open(SOURCE / "pair_eligibility.tsv.gz", "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = tuple(row[f] for f in FIELDS)
            if key in contexts:
                k = key + (row["gene_id"],)
                assert k not in delivered
                delivered[k] = (row["status"], row["reason"], int(row["n_joint_PAS"]))
    print("DELIVERED_OVERLAP_GENE_CONTEXTS", len(delivered), flush=True)
    summaries = {key: dict(zip(FIELDS, key), n_candidate_intersection_MC=len(row["shared"]),
        n_genes_checked=0, n_genes_with_qualified_intersection=0, n_original_eligible_pairs=0,
        n_original_eligible_pairs_with_qualified_intersection=0, n_joint_count_changed_genes=0,
        n_eligible_joint_count_changed_pairs=0, n_eligible_pairs_lost_after_dedup=0)
        for key, row in contexts.items()}
    changed_keys, old_eligible_with_intersection = [], []
    fields = FIELDS + ["gene_id", "qualified_MC_intersection", "n_qualified_MC_intersection", "old_status",
        "n_old_joint_PAS", "n_unique_joint_PAS", "n_lost_joint_PAS", "lost_joint_PAS_ids", "eligible_after_dedup"]
    with gzip.open(HERE / "shared_MC_joint_support.tsv.gz", "wt") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fields); writer.writeheader()
        checked = 0
        for study in STUDIES:
            selected_contexts = {key: row for key, row in contexts.items() if key[0] == study}
            if not selected_contexts:
                continue
            meta = metas[study]
            features = records(z, study + "_features_json"); features["column"] = np.arange(len(features))
            data, indices, ptr = (z[study + suffix] for suffix in ["_pas_data", "_pas_indices", "_pas_indptr"])
            for setting in sorted({key[3] for key in selected_contexts}):
                retained = features[features.coordinate_eligible]
                if setting == "supported50":
                    retained = retained[retained.matched_source_count_50nt.ge(1)]
                active = {key: row for key, row in selected_contexts.items() if key[3] == setting}
                for gene, sites in retained.groupby("primary_gene_id", sort=True):
                    if len(sites) < 2:
                        continue
                    counts = np.zeros((len(meta), len(sites)), dtype=data.dtype)
                    for j, column in enumerate(sites.column):
                        start, end = ptr[column:column+2]
                        counts[indices[start:end], j] = data[start:end]
                    depths = counts.sum(axis=1)
                    for key, row in active.items():
                        qa, qb = row["a"][depths[row["a"]] >= 20], row["b"][depths[row["b"]] >= 20]
                        shared, old_joint, distinct = joint_counts(counts, qa, qb)
                        reason = "insufficient_gene_qualified_metacells" if min(len(qa), len(qb)) < 3 else "insufficient_gene_qualified_cells" if min(meta.n_cells.iloc[qa].sum(), meta.n_cells.iloc[qb].sum()) < 30 else "fewer_than_two_joint_supported_PAS" if old_joint.sum() < 2 else ""
                        actual = delivered[key + (gene,)]
                        assert actual == ("not_measurable" if reason else "eligible", reason, int(old_joint.sum())), (key, gene, actual, reason, int(old_joint.sum()))
                        checked += 1
                        summary = summaries[key]; summary["n_genes_checked"] += 1
                        summary["n_original_eligible_pairs"] += int(not reason)
                        if not len(shared):
                            continue
                        summary["n_genes_with_qualified_intersection"] += 1
                        summary["n_original_eligible_pairs_with_qualified_intersection"] += int(not reason)
                        summary["n_joint_count_changed_genes"] += int(old_joint.sum() != distinct.sum())
                        summary["n_eligible_joint_count_changed_pairs"] += int(not reason and old_joint.sum() != distinct.sum())
                        lost = bool(not reason and distinct.sum() < 2)
                        summary["n_eligible_pairs_lost_after_dedup"] += int(lost)
                        if not reason:
                            old_eligible_with_intersection.append(key + (gene,))
                        if lost:
                            changed_keys.append(key + (gene,))
                        writer.writerow(dict(zip(FIELDS, key), gene_id=gene,
                            qualified_MC_intersection=json.dumps(meta.metacell_id.iloc[shared].tolist(), separators=(",", ":")),
                            n_qualified_MC_intersection=len(shared), old_status=actual[0],
                            n_old_joint_PAS=int(old_joint.sum()), n_unique_joint_PAS=int(distinct.sum()),
                            n_lost_joint_PAS=int((old_joint & ~distinct).sum()),
                            lost_joint_PAS_ids=json.dumps(sites.atlas_pas_id[old_joint & ~distinct].tolist(), separators=(",", ":")),
                            eligible_after_dedup=bool(not reason and distinct.sum() >= 2)))
            del data, indices, ptr
            print("STUDY_CHECKED", study, checked, "LOST_PAIRS", len(changed_keys), flush=True)
    assert checked == len(delivered)
    pd.DataFrame(summaries.values()).to_csv(HERE / "shared_MC_joint_support_contexts.tsv", sep="\t", index=False)
    pd.DataFrame(changed_keys, columns=FIELDS + ["gene_id"]).to_csv(HERE / "shared_MC_pairs_lost.tsv", sep="\t", index=False)
    pd.DataFrame(old_eligible_with_intersection, columns=FIELDS + ["gene_id"]).to_csv(HERE / "shared_MC_eligible_intersection_pairs.tsv.gz", sep="\t", index=False)
    report = dict(status="PASS", diagnostic_only=True, classification_or_effects_changed=False,
        n_candidate_overlap_contexts=len(contexts), n_gene_contexts_checked=checked,
        n_eligible_gene_pairs_with_qualified_intersection=len(old_eligible_with_intersection),
        n_eligible_pairs_lost_after_distinct_MC_dedup=len(changed_keys),
        n_gene_contexts_joint_count_changed=sum(s["n_joint_count_changed_genes"] for s in summaries.values()),
        n_eligible_pairs_joint_count_changed=sum(s["n_eligible_joint_count_changed_pairs"] for s in summaries.values()),
        verifier_sha256=sha(__file__), inputs_sha256=sha(SOURCE / "inputs.npz"))
    (HERE / "shared_MC_joint_support_validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
