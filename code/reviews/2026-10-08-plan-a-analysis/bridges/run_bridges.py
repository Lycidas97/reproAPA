"""Plan A items 4–6, on immutable manuscript-cohort tables."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "notes/2026-10-06-supplement-polyapipe/contract-r2/bundle/tables"
OUT = Path(__file__).resolve().parent
INPUTS = {}
TESTS = []


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read(name, **kwargs):
    p = BASE / name
    INPUTS[str(p)] = digest(p)
    return pd.read_csv(p, sep="\t", **kwargs)


def save(frame, name):
    pd.DataFrame(frame).to_csv(OUT / name, sep="\t", index=False)


def recurrent(contrast):
    f = read("recurrence/recurrence_members.tsv.gz")
    f = f[f.metric.eq("WPP") & f.threshold.eq(.02) & f.series_matching &
          f.variant.eq("primary") & f.contrast.eq(contrast)].copy()
    assert not f.duplicated(["species", "gene_id"]).any()
    return f


def human_feasibility(rec, dual):
    inv = read("intermediate/state_library_inventory.tsv")
    inv = inv[inv.setting.eq("supported50_p08") &
              inv.study_id.isin(["GSE178325", "GSE150311"])].copy()
    libs = read("frozen/state_libraries.tsv")
    dates = libs[["study_id", "library_key", "experimental_time"]].drop_duplicates()
    assert not dates.duplicated(["study_id", "library_key"]).any()
    inv = inv.merge(dates, on=["study_id", "library_key"], validate="many_to_one")
    save(inv, "A4-qualified-library-inventory.tsv")
    rows = []
    for sid, g in inv.groupby("study_id"):
        p = g[g.state.eq("P")]
        m = g[g.state.eq("M")]
        same = p.merge(m, on=["study_id", "route", "experimental_time"], suffixes=("_P", "_M"))
        rows.append(dict(study_id=sid, P_libraries=p.library_key.nunique(),
                         M_libraries=m.library_key.nunique(), same_recorded_time_PM_pairs=len(same),
                         same_library_PM_pairs=int(same.library_key_P.eq(same.library_key_M).sum()),
                         documented_source=sid == "GSE178325",
                         conclusion="no_qualified_same_recorded_time_PM" if not len(same) else
                         "observed_S_sensitivity_only_no_documented_source"))
        save(same, f"A4-{sid}-same-date-pairs.tsv")
    save(rows, "A4-feasibility.tsv")
    states = read("intermediate/library_state_values.tsv.gz")
    states = states[states.study_id.eq("GSE150311") & states.setting.eq("supported50_p08") &
                    states.state.isin(list("PMS"))].copy()
    keys = ["species", "study_id", "route", "gene_id", "gene_symbol", "library_key", "sample_unit", "time"]
    assert not states.duplicated(keys + ["state"]).any()
    wide = states.pivot(index=keys, columns="state", values="value").reindex(columns=list("PMS")).dropna().reset_index()
    pairs = read("intermediate/paired_library_effects.tsv.gz")
    pairs = pairs[pairs.study_id.eq("GSE150311") & pairs.setting.eq("supported50_p08") &
                  pairs.contrast.isin(["P_S", "P_M", "M_S"]) & pairs.n_joint_pas.ge(2)]
    qualified = pairs.groupby(keys).contrast.nunique()
    qualified = qualified[qualified.eq(3)].reset_index()[keys]
    wide = wide.merge(qualified, on=keys, validate="one_to_one")
    wide["DWUI_PS"] = wide.P - wide.S
    wide["DWUI_PM"] = wide.P - wide.M
    wide["DWUI_MS"] = wide.M - wide.S
    assert np.max(np.abs(wide.DWUI_PS - wide.DWUI_PM - wide.DWUI_MS)) < 1e-12
    save(wide, "A4-Liu2020-observed-S-library-triangles.tsv.gz")
    gkeys = ["species", "study_id", "route", "gene_id", "gene_symbol"]
    unit = wide.groupby(gkeys + ["time", "sample_unit"])[["DWUI_PS", "DWUI_PM", "DWUI_MS"]].mean().reset_index()
    daily = unit.groupby(gkeys + ["time"])[["DWUI_PS", "DWUI_PM", "DWUI_MS"]].mean().reset_index()
    gene = daily.groupby(gkeys)[["DWUI_PS", "DWUI_PM", "DWUI_MS"]].mean().reset_index()
    r = rec[rec.species.eq("human")]
    gene = r[["species", "gene_id", "gene_symbol", "direction"]].merge(
        gene.drop(columns="gene_symbol"), on=["species", "gene_id"], how="left", validate="one_to_one")
    sign = gene.direction.map({"positive": 1, "negative": -1})
    gene["measurable"] = gene.DWUI_PM.notna()
    gene["aligned_original_direction"] = gene.measurable & gene.DWUI_PM.mul(sign).ge(.02) & gene.DWUI_PS.mul(sign).ge(.02)
    gene["Xing_dual_member"] = gene.gene_id.isin(dual[dual.species.eq("human")].gene_id)
    save(gene, "A4-Liu2020-fixed-human-candidate-sensitivity.tsv")
    summary = [dict(candidates=len(gene), measurable=int(gene.measurable.sum()),
                    aligned=int(gene.aligned_original_direction.sum()),
                    Xing_dual_candidates=int(gene.Xing_dual_member.sum()),
                    Xing_dual_measurable=int((gene.Xing_dual_member & gene.measurable).sum()),
                    Xing_dual_aligned=int((gene.Xing_dual_member & gene.aligned_original_direction).sum()),
                    true_source_replication=False)]
    save(summary, "A4-Liu2020-sensitivity-summary.tsv")


def interval_bridge(rec, dual):
    regions = read("candidate_screen/recurrent_region_pairs.tsv.gz")
    exact = read("site_statistics/exact_site_cross_study.tsv.gz")
    model = read("intervals/model_effects.tsv")
    rows = []
    audits = []
    metrics = []
    for sp, routes in [("mouse", {"CIP", "WOT"}), ("human", {"GSE118258", "GSE150311"})]:
        candidates = rec[rec.species.eq(sp)].copy()
        d = set(dual[dual.species.eq(sp)].gene_id)
        strict = set(regions[regions.species.eq(sp) & regions.contrast.eq("MP")].gene_id)
        assert len(strict & d) == (71 if sp == "mouse" else 5)
        mod = model[model.species.eq(sp) & model.contrast.eq("MP") & model.route.isin(routes)]
        elig = mod[mod.eligible].groupby("gene_id").route.agg(set)
        eligible = set(elig[elig.map(lambda x: x == routes)].index)
        both = exact[exact.species.eq(sp) & exact.contrast.eq("MP")].copy()
        measured = set(both[both.both_observed].gene_id) & eligible
        assert strict & set(candidates.gene_id) <= measured
        for event in candidates.itertuples():
            audits.append(dict(species=sp, gene_id=event.gene_id, dual=event.gene_id in d,
                               modeled_both=event.gene_id in eligible,
                               observed_variable_site_interval=event.gene_id in measured,
                               strict_interval=event.gene_id in strict))
        yes = candidates[candidates.gene_id.isin(d & measured)]
        no = candidates[~candidates.gene_id.isin(d) & candidates.gene_id.isin(measured)]
        a, c = len(set(yes.gene_id) & strict), len(set(no.gene_id) & strict)
        b, e = len(yes) - a, len(no) - c
        if min(len(yes), len(no)) >= 3:
            odds, p = fisher_exact([[a, b], [c, e]])
            status = "tested_conditional_on_existing_union_pool_and_annotation"
        else:
            odds, p, status = np.nan, 1., "insufficient_measurable_groups"
        TESTS.append(dict(family="A5-interval-bridge", species=sp, test="strict_interval_rate", raw_p=p, status=status))
        rows.append(dict(species=sp, dual_all=len(d), dual_strict_all=len(d & strict),
                         nondual_all=len(candidates) - len(d),
                         dual_measurable=len(yes), dual_strict=a,
                         nondual_measurable=len(no), nondual_strict=c,
                         dual_rate=a / len(yes) if len(yes) else np.nan,
                         nondual_rate=c / len(no) if len(no) else np.nan,
                         odds_ratio=odds, raw_p=p, status=status))
        # Shared physical sites are the unit within a gene; genes are the testing unit.
        f = both[both.evaluable_a & both.evaluable_b].copy()
        f = f[f.gene_id.isin(candidates.gene_id)]
        for col in ["effect20_a", "effect50_a", "effect20_b", "effect50_b"]:
            # Constant full retention has an exact zero effect by construction.
            side = col[-1]
            const = f["retained_" + side]
            assert f.loc[~const, col].notna().all()
            assert np.allclose(f.loc[const, col].fillna(0), 0.)
            f.loc[const, col] = 0.
        physical = ["gene_id", "chrom", "strand", "start0", "end0"]
        effects = ["effect20_a", "effect50_a", "effect20_b", "effect50_b"]
        assert f.groupby(physical)[effects].nunique().le(1).all().all()
        f = f.drop_duplicates(physical)
        f["net_retention"] = f[effects].mean(axis=1)
        gm = f.groupby("gene_id").agg(net_retention=("net_retention", "mean"), measured_physical_sites=("start0", "size")).reset_index()
        gm["dual"] = gm.gene_id.isin(d)
        gm["species"] = sp
        metrics.append(gm)
        x, y = gm.loc[gm.dual, "net_retention"], gm.loc[~gm.dual, "net_retention"]
        if min(len(x), len(y)) >= 3:
            p = mannwhitneyu(x, y, alternative="two-sided", method="asymptotic").pvalue
            status = "tested_gene_equal_shared_predicted_site_inclusion"
        else:
            p, status = 1., "insufficient_measurable_groups"
        TESTS.append(dict(family="A5-interval-bridge", species=sp, test="site_net_retention", raw_p=p, status=status))
        rows[-1].update(dual_site_genes=len(x), nondual_site_genes=len(y),
                        dual_mean_site_retention=x.mean(), nondual_mean_site_retention=y.mean(),
                        dual_median_site_retention=x.median(), nondual_median_site_retention=y.median(),
                        site_raw_p=p, site_status=status)
    save(rows, "A5-interval-bridge-summary.tsv")
    save(audits, "A5-candidate-measurability.tsv")
    save(pd.concat(metrics, ignore_index=True), "A5-site-net-retention-per-gene.tsv")


def endpoint_shapes():
    members = read("profiles/gene_membership_k3.tsv.gz")
    events = read("recurrence/study_event_covariates.tsv.gz")
    events = events[events.metric.eq("WPP") & events.threshold.eq(.02) & events.contrast.eq("P_S")].copy()
    rec = recurrent("P_S")
    rows = []
    audits = []
    for sp, units in [("mouse", {"GSE103221", "GSE115943"}), ("human", {"GSE118258", "GSE178325"})]:
        g = members[members.species.eq(sp)].copy()
        ev = events[events.species.eq(sp) & events.unit.isin(units)]
        assert not ev.duplicated(["unit", "gene_id"]).any()
        by_gene = ev.groupby("gene_id").unit.agg(set)
        measurable = set(by_gene[by_gene.map(lambda x: x == units)].index)
        selected = set(rec[rec.species.eq(sp)].gene_id)
        g["measurable"] = g.gene_id.isin(measurable)
        g["recurrent"] = g.gene_id.isin(selected)
        g["measured_units"] = g.gene_id.map(by_gene.map(len)).fillna(0).astype(int)
        g["positive_votes"] = g.gene_id.map(ev[ev.event.eq(1)].groupby("gene_id").size()).fillna(0).astype(int)
        g["negative_votes"] = g.gene_id.map(ev[ev.event.eq(-1)].groupby("gene_id").size()).fillna(0).astype(int)
        assert not (g.recurrent & ~g.measurable).any()
        audits.append(g)
        for cluster, q in g.groupby("cluster"):
            m = q[q.measurable]
            rows.append(dict(species=sp, cluster=cluster, all_candidates=len(q), all_recurrent=int(q.recurrent.sum()),
                             measurable=len(m), recurrent=int(m.recurrent.sum()),
                             measurable_rate=m.recurrent.mean(), missing=len(q) - len(m)))
        for other in ["C1", "C3"]:
            a = g[g.cluster.eq("C2") & g.measurable]
            b = g[g.cluster.eq(other) & g.measurable]
            odds, p = fisher_exact([[int(a.recurrent.sum()), int((~a.recurrent).sum())],
                                   [int(b.recurrent.sum()), int((~b.recurrent).sum())]])
            TESTS.append(dict(family="A6-shape-endpoint", species=sp, test="C2_vs_" + other,
                              raw_p=p, odds_ratio=odds, rate_difference=a.recurrent.mean() - b.recurrent.mean(),
                              status="tested_measurable_in_both_endpoint_studies"))
    save(rows, "A6-measurable-endpoint-recurrence.tsv")
    save(pd.concat(audits, ignore_index=True), "A6-profile-member-audit.tsv")


def main():
    INPUTS[str(OUT / "PROTOCOL.md")] = digest(OUT / "PROTOCOL.md")
    rec = recurrent("P_M")
    dual = read("candidate_screen/dual_reference_candidates.tsv")
    assert dual.groupby("species").size().to_dict() == {"human": 128, "mouse": 222}
    assert rec.groupby("species").size().to_dict() == {"human": 268, "mouse": 420}
    human_feasibility(rec, dual)
    interval_bridge(rec, dual)
    endpoint_shapes()
    save(TESTS, "raw-tests.tsv")
    assert all(digest(Path(p)) == h for p, h in INPUTS.items())
    (OUT / "input-hashes.json").write_text(json.dumps(INPUTS, indent=2) + "\n")
    (OUT / "run-receipt.json").write_text(json.dumps(dict(
        status="complete", primary_source="polyApipe contract-r2", items=[4, 5, 6],
        source_code_sha256=digest(Path(__file__)), inference_tests=len(TESTS),
        python_version=__import__("sys").version, numpy_version=np.__version__, pandas_version=pd.__version__), indent=2) + "\n")
    print(pd.DataFrame(TESTS).to_string(index=False))


if __name__ == "__main__":
    main()
