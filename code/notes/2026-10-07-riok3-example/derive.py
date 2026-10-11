"""Riok3 rows for Figure 5f from the polyApipe main bundle (2026-10-07).

The bundle's example/geometry nodes covered only the low-shift cohort (Trip12), so Riok3 rows are derived here with the
same code paths as analysis_code/pas_rerun (downstream_example_inclusion.feature_libraries, downstream_geometry.decompose,
migrate.mask_for). Each function is first run on Trip12 and must reproduce the bundle rows exactly.
Outputs (tables/): example_endpoint_effects.tsv, inclusion_paired_libraries.tsv, inclusion-paired-means.tsv,
axis-decomposition.tsv, native-sites.tsv; receipt.json records input SHA-256 and the Trip12 checks.
"""
import json, hashlib, sys
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[2]
T = ROOT/'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables'
OUT = Path(__file__).resolve().parent/'tables'; OUT.mkdir(exist_ok=True)
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
inputs = {}
def rd(rel, **kw):
    p = T/rel; inputs[rel] = sha(p); return pd.read_csv(p, sep='\t', **kw)
def js(rel):
    p = T/rel; inputs[rel] = sha(p); return json.loads(p.read_text())

TRIP12, RIOK3 = 'ENSMUSG00000026219', 'ENSMUSG00000024404'
ROUTES = {'CIP': 'GSE103221', 'WOT': 'GSE115943'}
models = {r: js(f'intervals/annotation/mouse/{r}/models.json') for r in ROUTES}
mask_for = lambda start, end, strand, cuts: (end <= cuts if strand == '+' else start >= cuts-1).astype(int)

# library -> experimental time, from the bundle's own library-state inclusion table (all featured genes)
ls = rd('examples/inclusion_library_states.tsv.gz', usecols=['route', 'library_key', 'time'])
TIMES = {(a.route, a.library_key): str(a.time) for a in ls.drop_duplicates(['route', 'library_key']).itertuples()}

def qualified(model, state, cutoff):
    return {a['library']: a for a in model['profiles'].get(state, {}).get('libraries', []) if a['total'] >= cutoff and a['coverage'] >= .5}

def paired_inclusion(gene_id, route, start0, end0, label='interval'):
    m = models[route][gene_id]; strand = m['family']['strand']
    cuts = np.array([a['cut_site'] for a in m['pas']]); mask = mask_for(start0, end0, strand, cuts)
    fid = f'{route}|{gene_id}|interval|{start0}-{end0}|{label}'
    rows = []
    for cutoff in (20, 50):
        M, P = qualified(m, 'M', cutoff), qualified(m, 'P', cutoff)
        for lib in sorted(set(M) & set(P)):
            mi = float(np.asarray(M[lib]['counts'])@mask/M[lib]['total']); pi = float(np.asarray(P[lib]['counts'])@mask/P[lib]['total'])
            rows.append(dict(library_key=lib, time=TIMES.get((route, lib), ''), reference_inclusion=mi, target_inclusion=pi, delta=pi-mi,
                             species='mouse', study_id=ROUTES[route], route=route, gene_id=gene_id, gene=m['gene_symbol'].upper(),
                             start0=start0, end0=end0, kind='interval', label=label, feature_id=fid, strand=strand,
                             chrom=m['family']['chrom'], cutoff=cutoff))
    return pd.DataFrame(rows), mask

def paired_means(p):
    return (p.groupby(['feature_id', 'gene', 'route', 'cutoff'])
             .agg(n=('delta', 'size'), reference_mean=('reference_inclusion', 'mean'), target_mean=('target_inclusion', 'mean'), delta_mean=('delta', 'mean'))
             .reset_index())

eff = rd('state_vectors/route_effects.primary.tsv.gz')
eff = eff[eff.setting.eq('supported50_p08') & eff.route.isin(list(ROUTES)) & eff.contrast.eq('P_M') & eff.design.eq('same_library')]
axes = rd('state_vectors/PAS_axes.primary.tsv.gz', usecols=['study_id', 'setting', 'primary_gene_id', 'atlas_pas_id', 'cut_site', 'strand', 'primary_region_class', 'position', 'native_span'])

def decompose(gene_id):
    rows, sites = [], []
    for a in eff[eff.gene_id.eq(gene_id)].sort_values('route').itertuples():
        ids = json.loads(a.PAS_ids); P = np.array(json.loads(a.target_usage)); M = np.array(json.loads(a.reference_usage))
        ax = axes[axes.study_id.eq(a.study_id) & axes.setting.eq(a.setting) & axes.primary_gene_id.eq(a.gene_id)].set_index('atlas_pas_id').loc[ids]
        cuts = ax.cut_site.to_numpy(); strand = ax.strand.iloc[0]; signed = cuts*(1 if strand == '+' else -1); span = float(ax.native_span.iloc[0])
        np.testing.assert_allclose((P-M)@((signed-signed.min())/span), a.delta_DWUI, atol=1e-12)
        row = dict(gene_id=a.gene_id, gene_symbol=a.gene_symbol, route=a.route, study_id=a.study_id, native_pas=len(ids), axis_span_nt=span,
                   supported_pas_span_nt=float(np.ptp(signed)), n_native_pairs=a.n_pairs, native_delta_DWUI=a.delta_DWUI, native_MPRO=a.MPRO)
        m = models[a.route][a.gene_id]; fam = m['family']; lo, hi = fam['common_start0'], fam['common_end0']
        retained = cuts+(strand == '-'); mask = (retained > lo) & (retained <= hi)
        origin = lo if strand == '+' else hi; x = (cuts-origin)*(1 if strand == '+' else -1)
        mp, mm = P[mask].sum(), M[mask].sum()
        row.update(family_id=fam['family_id'], family_start0=lo, family_end0=hi, terminal_native_pas=int(mask.sum()),
                   terminal_pas_span_nt=float(np.ptp(signed[mask])) if mask.any() else np.nan, terminal_mass_P=mp, terminal_mass_M=mm)
        ok = mask.sum() >= 2 and min(mp, mm) > 0
        row['family_status'] = 'measurable_native_family' if ok else 'fewer_than_two_native_family_sites_or_zero_mass'
        if ok:
            up = P[mask]@x[mask]/mp; um = M[mask]@x[mask]/mm
            within = (mp+mm)/2*(up-um)/span; mass = (mp-mm)*(up+um)/2/span; outside = (P[~mask]-M[~mask])@x[~mask]/span
            np.testing.assert_allclose(within+mass+outside, a.delta_DWUI, atol=1e-12)
            row.update(terminal_position_P_nt=up, terminal_position_M_nt=um, terminal_shift_nt=up-um, within_family_DWUI=within,
                       family_mass_DWUI=mass, nonterminal_DWUI=outside, nonterminal_mass_delta=mm-mp, decomposition_residual=within+mass+outside-a.delta_DWUI)
        for thr in (20, 50):
            row[f'interval_paired_libraries_{thr}'] = len(set(qualified(m, 'P', thr)) & set(qualified(m, 'M', thr)))
        rows.append(row)
        for i, pas in enumerate(ids):
            sites.append(dict(gene_id=a.gene_id, route=a.route, pas_id=pas, cut_interbase=int(cuts[i]), strand=strand, terminal_family=bool(mask[i]),
                              relative_to_family_origin_nt=float(x[i]), P_usage=P[i], M_usage=M[i]))
    return pd.DataFrame(rows), pd.DataFrame(sites)

checks = {}
# ---- Trip12 reproduction checks
old_ip = rd('examples/inclusion_paired_libraries.tsv.gz'); old_ip = old_ip[old_ip.gene_id.eq(TRIP12) & old_ip.kind.eq('interval')]
mism = 0; n = 0
for (route, s, e), g in old_ip.groupby(['route', 'start0', 'end0']):
    new, _ = paired_inclusion(TRIP12, route, int(s), int(e))
    a = g.sort_values(['cutoff', 'library_key'])[['library_key', 'cutoff', 'reference_inclusion', 'target_inclusion', 'delta', 'time']].reset_index(drop=True)
    b = new.sort_values(['cutoff', 'library_key'])[['library_key', 'cutoff', 'reference_inclusion', 'target_inclusion', 'delta', 'time']].reset_index(drop=True)
    n += 1
    if not (len(a) == len(b) and (a.library_key.values == b.library_key.values).all() and np.allclose(a[['reference_inclusion', 'target_inclusion', 'delta']], b[['reference_inclusion', 'target_inclusion', 'delta']], atol=1e-12)
            and (a.time.astype(str).values == b.time.astype(str).values).all()):
        mism += 1
checks['trip12_inclusion_intervals'] = dict(intervals=n, mismatched=mism)
old_pm = rd('panel_statistics/inclusion-paired-means.tsv'); old_pm = old_pm[old_pm.gene.eq('TRIP12') & old_pm.feature_id.str.contains('|interval|', regex=False)]
new_pm = paired_means(pd.concat([paired_inclusion(TRIP12, r, int(s), int(e))[0] for (r, s, e), _ in old_ip.groupby(['route', 'start0', 'end0'])]))
mg = old_pm.merge(new_pm, on=['feature_id', 'cutoff'], suffixes=('_o', '_n'))
checks['trip12_paired_means'] = dict(rows=len(old_pm), merged=len(mg), max_abs_diff=float(np.abs(mg[['reference_mean_o', 'target_mean_o', 'delta_mean_o']].values - mg[['reference_mean_n', 'target_mean_n', 'delta_mean_n']].values).max()), n_equal=bool((mg.n_o == mg.n_n).all()))
old_ax = rd('geometry/axis-decomposition.tsv'); old_ax = old_ax[old_ax.gene_id.eq(TRIP12)].sort_values('route').reset_index(drop=True)
new_ax, _ = decompose(TRIP12)
num = ['native_pas', 'axis_span_nt', 'native_delta_DWUI', 'family_start0', 'family_end0', 'terminal_native_pas', 'terminal_shift_nt', 'within_family_DWUI', 'family_mass_DWUI', 'nonterminal_DWUI', 'interval_paired_libraries_20', 'interval_paired_libraries_50']
checks['trip12_decomposition'] = dict(columns=num, max_abs_diff=float(np.nanmax(np.abs(old_ax[num].astype(float).values - new_ax[num].astype(float).values))))
old_ee = rd('examples/example_endpoint_effects.tsv.gz'); old_ee = old_ee[old_ee.gene_id.eq(TRIP12) & old_ee.setting.eq('supported50_p08') & old_ee.contrast.eq('P_M')].sort_values('route')
new_ee = eff[eff.gene_id.eq(TRIP12)].sort_values('route')
cols = [c for c in ['n_pairs', 'MPRO', 'delta_DWUI', 'reference_DWUI', 'target_DWUI'] if c in new_ee and c in old_ee]
checks['trip12_endpoint'] = dict(columns=cols, max_abs_diff=float(np.abs(old_ee[cols].astype(float).values - new_ee[cols].astype(float).values).max()),
                                 PAS_ids_equal=bool((old_ee.PAS_ids.values == new_ee.PAS_ids.values).all()),
                                 usage_equal=bool((old_ee.target_usage.values == new_ee.target_usage.values).all() and (old_ee.reference_usage.values == new_ee.reference_usage.values).all()))
# ---- Riok3
IV889 = (12288586, 12289475); IVPOLY = (12288550, 12289754)
ip_rows, masks = [], {}
for route in ROUTES:
    for (s, e), lab in [(IV889, 'common_interval_889nt'), (IVPOLY, 'polyApipe_strict_interval')]:
        p, mk = paired_inclusion(RIOK3, route, s, e, label=lab); ip_rows.append(p); masks[f'{route}:{s}-{e}'] = mk.tolist()
ip = pd.concat(ip_rows, ignore_index=True); pm = paired_means(ip)
ax, sites = decompose(RIOK3)
ee = eff[eff.gene_id.eq(RIOK3)].sort_values('route').copy(); ee['status'] = 'measured'; ee['reason'] = ''
for name, df in [('example_endpoint_effects.tsv', ee), ('inclusion_paired_libraries.tsv', ip), ('inclusion-paired-means.tsv', pm),
                 ('axis-decomposition.tsv', ax), ('native-sites.tsv', sites)]:
    df.to_csv(OUT/name, sep='\t', index=False)
summary = dict(
    native={r.route: dict(delta_DWUI=r.delta_DWUI, MPRO=r.MPRO, n_pairs=int(r.n_pairs), reference_DWUI=getattr(r, 'reference_DWUI', None), target_DWUI=getattr(r, 'target_DWUI', None)) for r in ee.itertuples()},
    family={r.route: dict(terminal_native_pas=int(r.terminal_native_pas), terminal_shift_nt=r.terminal_shift_nt, within_family_DWUI=r.within_family_DWUI,
                          family_mass_DWUI=r.family_mass_DWUI, nonterminal_DWUI=r.nonterminal_DWUI, terminal_mass_P=r.terminal_mass_P, terminal_mass_M=r.terminal_mass_M,
                          axis_span_nt=r.axis_span_nt, family=[int(r.family_start0), int(r.family_end0)]) for r in ax.itertuples()},
    outside_usage={k: dict(P=float(g.loc[~g.terminal_family, 'P_usage'].sum()), M=float(g.loc[~g.terminal_family, 'M_usage'].sum())) for k, g in sites.groupby('route')},
    inclusion_pp={f'{r.route}|{r.feature_id.split("|")[3]}|{r.cutoff}': dict(n=int(r.n), M=r.reference_mean*100, P=r.target_mean*100, delta=r.delta_mean*100) for r in pm.itertuples()},
    masks=masks)
json.dump(dict(inputs=inputs, checks=checks, riok3=summary), open(Path(__file__).resolve().parent/'receipt.json', 'w'), indent=1, default=float)
print(json.dumps(checks, indent=1)); print(json.dumps(summary, indent=1, default=float))
