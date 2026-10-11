"""Major 2 (two-stage review): where do polyApipe strict recurrent interval pairs stand with the alternative PAS calls?
Unit: matched interval pair (route a/b regions, deduplicated by species, contrast, gene and the two region ids) with
recurrent_strict in polyApipe. Each physical region is looked up in the APADETECT p50 coordinate regions of the same route,
gene and contrast: full retention of [start0, end0) requires a PAS 3' of the region's 3'-most base, so its inclusion under
the other calls equals that of the p50 region containing that base (+ strand: end0-1; - strand: start0).
Classes per pair: reproduced (p50 strict, same sign, both studies); measurable_below_threshold (both measured, same sign or
zero change, not strict in >=1 study); direction_conflict (opposite sign in >=1 study); not_measurable (no p50 region
containing the base, region not eligible, or no effect in >=1 study)."""
import json, hashlib
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[2]; HERE = Path(__file__).resolve().parent
PP = ROOT/'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables'; PA = ROOT/'notes/2026-10-06-dual-main-plot-data/apadetect-p50/main/tables'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
F = [PP/'candidate_screen/recurrent_region_pairs.tsv.gz', PA/'intervals/all_coordinate_regions.tsv.gz']
rp = pd.read_csv(F[0], sep='\t', low_memory=False)
rp = rp[rp.recurrent_strict.astype(str).eq('True')]
pairs = rp.drop_duplicates(['species', 'contrast', 'gene_id', 'region_id_a', 'region_id_b']).copy()
reg = pd.read_csv(F[1], sep='\t', low_memory=False)
me = pd.read_csv(PA/'intervals/model_effects.tsv', sep='\t'); F.append(PA/'intervals/model_effects.tsv')
MODELLED = set(zip(me.route, me.gene_id, me.contrast))
reg = reg[reg.eligible.astype(str).eq('True')]
idx = {k: g.sort_values('start0') for k, g in reg.groupby(['route', 'gene_id', 'contrast'])}
def look(route, gene, contrast, s, e, strand):
    g = idx.get((route, gene, contrast))
    if g is None: return None
    b = e-1 if strand == '+' else s
    h = g[(g.start0 <= b) & (g.end0 > b)]
    return None if len(h) != 1 else h.iloc[0]
rows = []
for p in pairs.itertuples():
    out = dict(species=p.species, contrast=p.contrast, gene_id=p.gene_id, gene=p.gene_symbol_a, strand=p.strand)
    sig, st, meas = [], [], []
    for side in 'ab':
        route = getattr(p, 'route_'+side); s, e = int(getattr(p, 'region_start0_'+side)), int(getattr(p, 'region_end0_'+side))
        e20p = float(getattr(p, 'effect20_'+side))
        h = look(route, p.gene_id, p.contrast, s, e, p.strand)
        out.update({f'route_{side}': route, f'poly_region_{side}': f'{s}-{e}', f'poly_effect20_{side}': e20p})
        if h is None or pd.isna(h.effect20) or pd.isna(h.effect50):
            out.update({f'alt_region_{side}': '', f'alt_effect20_{side}': np.nan, f'alt_effect50_{side}': np.nan, f'alt_strict_{side}': False})
            meas.append(False); sig.append(0); st.append(False); continue
        out.update({f'alt_region_{side}': f'{int(h.start0)}-{int(h.end0)}', f'alt_effect20_{side}': h.effect20, f'alt_effect50_{side}': h.effect50,
                    f'alt_strict_{side}': str(h.strict_change) == 'True'})
        meas.append(True); sig.append(int(np.sign(h.effect20)) * int(np.sign(e20p))); st.append(str(h.strict_change) == 'True')
    if not all(meas):
        c = 'not_measurable'
        modelled = all(((getattr(p, 'route_'+s), p.gene_id, p.contrast) in MODELLED) for s in 'ab')
        out['not_measurable_reason'] = 'outside_alternative_family_or_no_effect' if modelled else 'gene_not_in_alternative_interval_pool'
    elif any(x < 0 for x in sig): c = 'direction_conflict'
    elif all(st) and all(x > 0 for x in sig): c = 'reproduced'
    else: c = 'measurable_below_threshold'
    out['class'] = c
    out['alt_strict_studies'] = sum(st) if all(meas) else np.nan
    rows.append(out)
R = pd.DataFrame(rows)
R.to_csv(HERE/'polyapipe-strict-interval-pairs-under-alternative-calls.tsv', sep='\t', index=False)
order = ['reproduced', 'measurable_below_threshold', 'direction_conflict', 'not_measurable']
S = (R.groupby(['species', 'contrast', 'class']).size().unstack(fill_value=0).reindex(columns=order, fill_value=0))
S['interval_pairs'] = S.sum(axis=1)
S = S.reset_index()
# below threshold: how many were strict in one of the two studies
bt = R[R['class'].eq('measurable_below_threshold')].groupby(['species', 'contrast']).alt_strict_studies.apply(lambda x: int((x == 1).sum())).rename('below_threshold_strict_in_one_study')
S = S.merge(bt.reset_index(), how='left', on=['species', 'contrast']).fillna({'below_threshold_strict_in_one_study': 0})
nr = R[R['class'].eq('not_measurable')].groupby(['species', 'contrast', 'not_measurable_reason']).size().unstack(fill_value=0)
S = S.merge(nr.reset_index(), how='left', on=['species', 'contrast']).fillna(0)
for c in S.columns[2:]: S[c] = S[c].astype(int)
S.to_csv(HERE/'interval-pair-classes.tsv', sep='\t', index=False)
json.dump({str(f.relative_to(ROOT)): sha(f) for f in F}, open(HERE/'inputs.json', 'w'), indent=1)
print(S.to_string())
print(R[R.gene.eq('Trip12')].T.to_string())
