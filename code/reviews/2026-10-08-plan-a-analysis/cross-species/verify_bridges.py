"""Independent replay of root's A4–A6 denominators and gene-level statistics."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BRIDGES = HERE.parent/'bridges'
BASE = ROOT/'notes/2026-10-06-supplement-polyapipe/contract-r2/bundle/tables'
READS = {}


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read(p):
    READS[str(p)] = digest(p)
    return pd.read_csv(p,sep='\t',float_precision='round_trip')


def close(a,b):
    assert np.isclose(a,b,rtol=1e-12,atol=1e-14,equal_nan=True),(a,b)


def main():
    h = json.loads((BRIDGES/'input-hashes.json').read_text())
    assert all(digest(p)==sha for p,sha in h.items())
    receipt = json.loads((BRIDGES/'run-receipt.json').read_text())
    assert digest(BRIDGES/'run_bridges.py') == receipt['source_code_sha256']
    READS[str(BRIDGES/'run_bridges.py')] = receipt['source_code_sha256']
    READS[str(BRIDGES/'PROTOCOL.md')] = digest(BRIDGES/'PROTOCOL.md')
    tests = read(BRIDGES/'raw-tests.tsv')
    assert tests.groupby('family').size().to_dict() == {'A5-interval-bridge':4,'A6-shape-endpoint':4}
    rec = read(BASE/'recurrence/recurrence_members.tsv.gz')
    rec = rec[rec.variant.eq('primary')&rec.series_matching&rec.metric.eq('WPP')&rec.threshold.eq(.02)]
    dual = read(BASE/'candidate_screen/dual_reference_candidates.tsv')

    # A4: qualification and physical-library triangles; no replacement of S meaning.
    inventory = read(BASE/'intermediate/state_library_inventory.tsv')
    inventory = inventory[inventory.setting.eq('supported50_p08')&inventory.study_id.isin(['GSE150311','GSE178325'])]
    assert inventory.n_mc.ge(3).all() and inventory.n_cells.ge(30).all()
    metadata = read(BASE/'frozen/state_libraries.tsv')
    times = metadata[['study_id','library_key','experimental_time']].drop_duplicates()
    inv = inventory.merge(times,on=['study_id','library_key'],validate='many_to_one')
    feasibility = read(BRIDGES/'A4-feasibility.tsv')
    assert 'same_recorded_time_PM_pairs' in feasibility.columns
    assert 'same_date_PM_pairs' not in feasibility.columns
    for sid,q in inv.groupby('study_id'):
        p,m = q[q.state.eq('P')],q[q.state.eq('M')]
        pairs = p.merge(m,on=['study_id','route','experimental_time'],suffixes=('_P','_M'))
        r = feasibility[feasibility.study_id.eq(sid)].iloc[0]
        assert (p.library_key.nunique(),m.library_key.nunique(),len(pairs),int(pairs.library_key_P.eq(pairs.library_key_M).sum())) == (r.P_libraries,r.M_libraries,r.same_recorded_time_PM_pairs,r.same_library_PM_pairs)
        if sid == 'GSE178325':
            assert r.conclusion == 'no_qualified_same_recorded_time_PM'
    raw_states = read(BASE/'intermediate/library_state_values.tsv.gz')
    raw_states = raw_states[raw_states.study_id.eq('GSE150311')&raw_states.setting.eq('supported50_p08')&raw_states.state.isin(['P','M','S'])]
    assert raw_states.n_mc.ge(3).all() and raw_states.n_cells.ge(30).all() and raw_states.pas_counts.ge(20).all()
    edges = read(BASE/'intermediate/paired_library_effects.tsv.gz')
    edges = edges[edges.study_id.eq('GSE150311')&edges.setting.eq('supported50_p08')&edges.contrast.isin(['P_S','P_M','M_S'])&edges.n_joint_pas.ge(2)]
    keys = ['species','study_id','route','gene_id','gene_symbol','library_key','sample_unit','time']
    state_values = raw_states.pivot(index=keys,columns='state',values='value').reindex(columns=['P','M','S']).dropna()
    edgecounts = edges.groupby(keys).contrast.nunique()
    valid = state_values.loc[state_values.index.intersection(edgecounts[edgecounts.eq(3)].index)]
    triangle = read(BRIDGES/'A4-Liu2020-observed-S-library-triangles.tsv.gz')
    assert set(valid.index.tolist()) == set(pd.MultiIndex.from_frame(triangle[keys]).tolist())
    x = valid.reset_index()
    for contrast,a,b in [('PS','P','S'),('PM','P','M'),('MS','M','S')]:
        x['DWUI_'+contrast] = x[a]-x[b]
    gene_keys = ['species','study_id','route','gene_id','gene_symbol']
    cols = ['DWUI_PS','DWUI_PM','DWUI_MS']
    mean = x.groupby(gene_keys+['time','sample_unit'])[cols].mean().groupby(gene_keys+['time']).mean().groupby(gene_keys).mean().reset_index()
    cand = rec[rec.species.eq('human')&rec.contrast.eq('P_M')][['gene_id','direction']].merge(mean,on='gene_id',how='left',validate='one_to_one')
    sign = cand.direction.map({'positive':1,'negative':-1})
    cand['aligned'] = cand.DWUI_PM.mul(sign).ge(.02)&cand.DWUI_PS.mul(sign).ge(.02)
    cand['dual'] = cand.gene_id.isin(dual[dual.species.eq('human')].gene_id)
    s = read(BRIDGES/'A4-Liu2020-sensitivity-summary.tsv').iloc[0]
    assert (len(cand),int(cand.DWUI_PM.notna().sum()),int(cand.aligned.sum()),int(cand.dual.sum()),
            int((cand.dual&cand.DWUI_PM.notna()).sum()),int((cand.dual&cand.aligned).sum())) == (s.candidates,s.measurable,s.aligned,s.Xing_dual_candidates,s.Xing_dual_measurable,s.Xing_dual_aligned)
    assert not s.true_source_replication
    print('A4: 268 candidates; 122 observed-S measurable; 71 aligned; Xing dual 59 measurable/34 aligned. One physical Liu triangle library; true_source_replication=False.')

    # A5: exclude unmeasured genes, then remove duplicated factor/category site links.
    regions = read(BASE/'candidate_screen/recurrent_region_pairs.tsv.gz')
    exact = read(BASE/'site_statistics/exact_site_cross_study.tsv.gz')
    models = read(BASE/'intervals/model_effects.tsv')
    summary = read(BRIDGES/'A5-interval-bridge-summary.tsv')
    per_gene = read(BRIDGES/'A5-site-net-retention-per-gene.tsv')
    denom_audit = read(BRIDGES/'A5-candidate-measurability.tsv')
    a5_rows = []
    for species,routes in [('mouse',{'CIP','WOT'}),('human',{'GSE118258','GSE150311'})]:
        candidates = set(rec.loc[rec.species.eq(species)&rec.contrast.eq('P_M'),'gene_id'])
        d = set(dual.loc[dual.species.eq(species),'gene_id'])
        strict = set(regions.loc[regions.species.eq(species)&regions.contrast.eq('MP'),'gene_id'])
        mod = models[models.species.eq(species)&models.contrast.eq('MP')&models.route.isin(routes)&models.eligible]
        modroutes = mod.groupby('gene_id').route.agg(set)
        eligible = set(modroutes[modroutes.map(lambda v:v==routes)].index)
        e = exact[exact.species.eq(species)&exact.contrast.eq('MP')]
        assert set(e.route_a.unique())|set(e.route_b.unique()) == routes
        measured = set(e.loc[e.both_observed,'gene_id'])&eligible
        assert strict&candidates <= measured
        yes, no = d&candidates&measured,(candidates-d)&measured
        counts = [[len(yes&strict),len(yes-strict)],[len(no&strict),len(no-strict)]]
        odds,p = fisher_exact(counts)
        r = summary[summary.species.eq(species)].iloc[0]
        assert (len(yes),len(no),counts[0][0],counts[1][0])==(r.dual_measurable,r.nondual_measurable,r.dual_strict,r.nondual_strict)
        close(odds,r.odds_ratio);close(p,r.raw_p)
        close(p,tests[tests.species.eq(species)&tests.test.eq('strict_interval_rate')].raw_p.iloc[0])
        audit = denom_audit[denom_audit.species.eq(species)]
        assert set(audit.gene_id)==candidates
        assert set(audit.loc[audit.observed_variable_site_interval,'gene_id'])==measured&candidates
        f = e[e.evaluable_a&e.evaluable_b&e.gene_id.isin(candidates)].copy()
        cols = ['effect20_a','effect50_a','effect20_b','effect50_b']
        assert np.isfinite(f[cols]).all().all(), 'Measured effects cannot be filled from missing.'
        for side in ['a','b']:
            retained = f['retained_'+side]
            assert f.loc[retained,'status_'+side].eq('constitutive_retained').all()
            assert f.loc[~retained,'status_'+side].eq('variable_observed').all()
            assert np.allclose(f.loc[retained,['effect20_'+side,'effect50_'+side]],0.)
            # A full-retention mask has exactly zero change; remove summation epsilon.
            f.loc[retained,['effect20_'+side,'effect50_'+side]] = 0.
        physical = ['gene_id','chrom','strand','start0','end0']
        assert f.groupby(physical)[cols].nunique().le(1).all().all()
        physical_f = f.drop_duplicates(physical)
        physical_f['site_effect'] = physical_f[cols].mean(axis=1)
        gm = physical_f.groupby('gene_id').agg(net_retention=('site_effect','mean'),measured_physical_sites=('start0','size'))
        published = per_gene[per_gene.species.eq(species)].set_index('gene_id')
        assert set(gm.index)==set(published.index)
        assert np.allclose(gm.net_retention.sort_index(),published.net_retention.sort_index(),rtol=1e-12,atol=1e-14)
        assert np.array_equal(gm.measured_physical_sites.sort_index(),published.measured_physical_sites.sort_index())
        vx = gm.loc[gm.index.isin(d),'net_retention'];vy=gm.loc[~gm.index.isin(d),'net_retention']
        psite = mannwhitneyu(vx,vy,alternative='two-sided',method='asymptotic').pvalue
        close(psite,r.site_raw_p)
        close(psite,tests[tests.species.eq(species)&tests.test.eq('site_net_retention')].raw_p.iloc[0])
        assert (len(vx),len(vy))==(r.dual_site_genes,r.nondual_site_genes)
        a5_rows.append(dict(species=species,rate_2x2=counts,missing_interval_candidates=len(candidates-measured),
            factor_link_rows=len(f),physical_sites=len(physical_f),dual_site_genes=len(vx),nondual_site_genes=len(vy),
            interval_p=p,site_p=psite,measured_retained_missing=0))
        print('A5:',json.dumps(a5_rows[-1]))

    # A6: membership in measured events, including event=0, determines eligibility.
    profile = read(BASE/'profiles/gene_membership_k3.tsv.gz')
    events = read(BASE/'recurrence/study_event_covariates.tsv.gz')
    events = events[events.metric.eq('WPP')&events.threshold.eq(.02)&events.contrast.eq('P_S')]
    endpoints = read(BRIDGES/'A6-measurable-endpoint-recurrence.tsv')
    member_audit = read(BRIDGES/'A6-profile-member-audit.tsv')
    a6_rows = []
    for species,units in [('mouse',{'GSE103221','GSE115943'}),('human',{'GSE118258','GSE178325'})]:
        g = profile[profile.species.eq(species)].copy()
        e = events[events.species.eq(species)&events.unit.isin(units)]
        assert not e.duplicated(['unit','gene_id']).any()
        assert np.isfinite(e[['effect_mean','effect_min','effect_max']]).all().all()
        measurable = set.intersection(*(set(e.loc[e.unit.eq(u),'gene_id']) for u in units))
        selected = set(rec.loc[rec.species.eq(species)&rec.contrast.eq('P_S'),'gene_id'])
        for cluster,q in g.groupby('cluster'):
            ids = set(q.gene_id)&measurable
            r = endpoints[endpoints.species.eq(species)&endpoints.cluster.eq(cluster)].iloc[0]
            assert (len(q),len(ids),len(ids&selected))==(r.all_candidates,r.measurable,r.recurrent)
        for other in ['C1','C3']:
            a = set(g.loc[g.cluster.eq('C2'),'gene_id'])&measurable
            b = set(g.loc[g.cluster.eq(other),'gene_id'])&measurable
            cells = [[len(a&selected),len(a-selected)],[len(b&selected),len(b-selected)]]
            odds,p = fisher_exact(cells)
            r = tests[tests.species.eq(species)&tests.test.eq('C2_vs_'+other)].iloc[0]
            close(p,r.raw_p);close(odds,r.odds_ratio)
            close(len(a&selected)/len(a)-len(b&selected)/len(b),r.rate_difference)
            a6_rows.append(dict(species=species,test='C2_vs_'+other,cells=cells,p=p))
        ma = member_audit[member_audit.species.eq(species)]
        assert set(ma.gene_id)==set(g.gene_id)
        assert set(ma.loc[ma.measurable,'gene_id'])==set(g.gene_id)&measurable
        assert set(ma.loc[ma.recurrent,'gene_id'])==set(g.gene_id)&selected
        print('A6:',species,'all measurable endpoint event=0 entries',int(e.event.eq(0).sum()),'fully measurable profile',int(ma.measurable.sum()))
    assert all(digest(p)==sha for p,sha in READS.items())
    result = dict(status='PASS',tests_recomputed=8,tests_skipped=0,A4_observed_S_is_sensitivity_only=True,
        A4_symbolic_time_labels_are_not_verified_dates=True,
        A4_final_label_same_recorded_time_not_same_date_verified=True,
        A5_missing_not_negative=True,A5_physical_sites_not_factor_links=True,
        A6_measurable_not_threshold_denominator=True,A5=a5_rows,A6=a6_rows,
        inputs_outputs_SHA256=READS,verifier_sha256=digest(Path(__file__)))
    (HERE/'bridges-counts-verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='inputs_outputs_SHA256'},indent=2))


if __name__ == '__main__':
    main()
