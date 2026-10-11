"""Opposite-sign genes (netshift-local, mouse P-M): decomposition and strict-interval inclusion, polyApipe main bundle.
Functions copied from notes/2026-10-07-riok3-example/derive.py (same code path); Riok3 recomputed as a check against Table S7c values."""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[3]
T = ROOT/'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables'
OUT = Path(__file__).resolve().parent
ROUTES = {'CIP': 'GSE103221', 'WOT': 'GSE115943'}
models = {r: json.loads((T/f'intervals/annotation/mouse/{r}/models.json').read_text()) for r in ROUTES}
mask_for = lambda start, end, strand, cuts: (end <= cuts if strand == '+' else start >= cuts-1).astype(int)
def qualified(model, state, cutoff):
    return {a['library']: a for a in model['profiles'].get(state, {}).get('libraries', []) if a['total'] >= cutoff and a['coverage'] >= .5}
def paired_inclusion(gene_id, route, s, e):
    m = models[route][gene_id]; strand = m['family']['strand']; cuts = np.array([a['cut_site'] for a in m['pas']]); mask = mask_for(s, e, strand, cuts)
    out = {}
    for c in (20, 50):
        M, P = qualified(m, 'M', c), qualified(m, 'P', c); libs = sorted(set(M) & set(P))
        d = [float(np.asarray(P[l]['counts'])@mask/P[l]['total'] - np.asarray(M[l]['counts'])@mask/M[l]['total']) for l in libs]
        out[c] = dict(n=len(libs), dI_pp=float(np.mean(d))*100 if d else None)
    return out
eff = pd.read_csv(T/'state_vectors/route_effects.primary.tsv.gz', sep='\t')
eff = eff[eff.setting.eq('supported50_p08') & eff.route.isin(list(ROUTES)) & eff.contrast.eq('P_M') & eff.design.eq('same_library')]
axes = pd.read_csv(T/'state_vectors/PAS_axes.primary.tsv.gz', sep='\t', usecols=['study_id','setting','primary_gene_id','atlas_pas_id','cut_site','strand','native_span'])
def decompose(gid):
    rows = []
    for a in eff[eff.gene_id.eq(gid)].sort_values('route').itertuples():
        ids = json.loads(a.PAS_ids); P = np.array(json.loads(a.target_usage)); M = np.array(json.loads(a.reference_usage))
        ax = axes[axes.study_id.eq(a.study_id) & axes.setting.eq(a.setting) & axes.primary_gene_id.eq(gid)].set_index('atlas_pas_id').loc[ids]
        cuts = ax.cut_site.to_numpy(); strand = ax.strand.iloc[0]; signed = cuts*(1 if strand == '+' else -1); span = float(ax.native_span.iloc[0])
        np.testing.assert_allclose((P-M)@((signed-signed.min())/span), a.delta_DWUI, atol=1e-12)
        r = dict(route=a.route, n_pairs=int(a.n_pairs), dDWUI=a.delta_DWUI, MPRO=a.MPRO, native_pas=len(ids), span_nt=span)
        m = models[a.route].get(gid)
        if m is None: r['family_status'] = 'no_interval_model'; rows.append(r); continue
        fam = m['family']; lo, hi = fam['common_start0'], fam['common_end0']
        retained = cuts+(strand == '-'); mask = (retained > lo) & (retained <= hi)
        origin = lo if strand == '+' else hi; x = (cuts-origin)*(1 if strand == '+' else -1)
        mp, mm = P[mask].sum(), M[mask].sum()
        r.update(family=[lo, hi], terminal_pas=int(mask.sum()), outside_usage_P=float(P[~mask].sum()), outside_usage_M=float(M[~mask].sum()))
        if mask.sum() >= 2 and min(mp, mm) > 0:
            up = P[mask]@x[mask]/mp; um = M[mask]@x[mask]/mm
            within = (mp+mm)/2*(up-um)/span; mass = (mp-mm)*(up+um)/2/span; outside = (P[~mask]-M[~mask])@x[~mask]/span
            np.testing.assert_allclose(within+mass+outside, a.delta_DWUI, atol=1e-12)
            r.update(family_status='measurable', terminal_shift_nt=up-um, within=within, mass=mass, outside=outside)
        else: r['family_status'] = 'fewer_than_two_family_PAS_or_zero_mass'
        rows.append(r)
    return rows
rr = pd.read_csv(T/'candidate_screen/recurrent_region_pairs.tsv.gz', sep='\t')
rr = rr[rr.contrast.eq('MP') & rr.recurrent_strict.astype(str).eq('True')]
rec = pd.read_csv(T/'recurrence/recurrence_members.tsv.gz', sep='\t')
rec = rec[rec.metric.eq('WPP') & rec.species.eq('mouse') & rec.contrast.eq('P_M') & rec.series_matching.astype(str).eq('True') & rec.variant.eq('primary') & np.isclose(rec.threshold, .02)]
GENES = {'Riok3': 'ENSMUSG00000024404', 'Fgfr1op2': 'ENSMUSG00000040242', 'Tmem183a': 'ENSMUSG00000042305',
         'Rbm39': None, 'Fxr1': None, 'Fgfr1op2_': None}
sym2id = dict(zip(rr.gene_symbol_a, rr.gene_id))
res = {}
for g in ['Riok3', 'Fgfr1op2', 'Tmem183a', 'Rbm39', 'Fxr1']:
    gid = GENES.get(g) or sym2id[g]
    iv = rr[rr.gene_id.eq(gid)][['shared_start0', 'shared_end0']].drop_duplicates().sort_values('shared_start0')
    ivs = []
    for s, e in iv.itertuples(index=False):
        ivs.append(dict(interval=[int(s), int(e)], length=int(e-s), **{rt: paired_inclusion(gid, rt, int(s), int(e)) for rt in ROUTES}))
    m = rec[rec.gene_id.eq(gid)]
    res[g] = dict(gene_id=gid, dwui_pm_recurrent=(m.direction.iloc[0] if len(m) else 'not_recurrent'), decomposition=decompose(gid), strict_intervals=ivs)
# p50: recurrence membership and strict intervals
P5 = ROOT/'notes/2026-10-06-supplement-apadetect-p50/contract-r2/bundle/tables'
r5 = pd.read_csv(P5/'recurrence/recurrence_members.tsv.gz', sep='\t')
r5 = r5[r5.metric.eq('WPP') & r5.species.eq('mouse') & r5.contrast.eq('P_M') & r5.series_matching.astype(str).eq('True') & r5.variant.eq('primary') & np.isclose(r5.threshold, .02)]
q5 = pd.read_csv(P5/'candidate_screen/recurrent_region_pairs.tsv.gz', sep='\t'); q5 = q5[q5.contrast.eq('MP')]
e5 = pd.read_csv(P5/'endpoints/route_effects.tsv.gz', sep='\t') if (P5/'endpoints/route_effects.tsv.gz').exists() else None
for g, v in res.items():
    gid = v['gene_id']; z = q5[q5.gene_id.eq(gid)]
    v['p50'] = dict(dwui_pm_recurrent=(r5[r5.gene_id.eq(gid)].direction.iloc[0] if r5.gene_id.eq(gid).any() else 'not_recurrent'),
                    interval_pairs=int(len(z)), strict_pairs=int(z.recurrent_strict.astype(str).eq('True').sum()))
    if e5 is not None:
        s = e5[e5.gene_id.eq(gid) & e5.contrast.eq('P_M') & e5.route.isin(['CIP','WOT'])]
        if 'setting' in s: s = s[s.setting.astype(str).str.contains('supported50_p08')]
        v['p50']['dDWUI'] = {r.route: float(r.delta_DWUI) for r in s.itertuples()} if len(s) else 'none'
json.dump(res, open(OUT/'opposite-sign.json', 'w'), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
