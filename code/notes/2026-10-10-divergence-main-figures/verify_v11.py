"""Independent checks for out-v11/Figure.3.pdf: kept c-e region object-identical to Fig/Figure.4.pdf, a and b values
recomputed from the source tables, Fig/ untouched (sha list written before the run), set copies byte-identical."""
import json, hashlib, sys
from pathlib import Path
import numpy as np, pandas as pd, pdfplumber
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]; OUT = HERE/'out-v11'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
rec = json.loads((OUT/'receipt.json').read_text()); r3 = rec['Figure.3']; res = {}
before = dict((Path(l.split()[1]).name, l.split()[0]) for l in Path(sys.argv[1]).read_text().splitlines())
res['fig_unchanged'] = all(sha(ROOT/'Fig'/k) == v for k, v in before.items())
def objs(pdf, region):
    pg = pdfplumber.open(pdf).pages[0]; x0, y0, x1, y1 = region; out = []
    for kind in ('chars', 'rects', 'lines', 'curves'):
        for o in getattr(pg, kind):
            if o['x1']-o['x0'] > 400 and o['bottom']-o['top'] > 400: continue   # page-size background rects (cleared by crop)
            cx, cy = (o['x0']+o['x1'])/2, (o['top']+o['bottom'])/2
            if x0 <= cx <= x1 and y0 <= cy <= y1:
                out.append((kind, o.get('text', ''), round(o['x0'], 1), round(o['top'], 1), round(o['x1'], 1), round(o['bottom'], 1)))
    return sorted(out)
reg = r3['kept'][0]['region']
a, b = set(objs(ROOT/'Fig/Figure.4.pdf', reg)), set(objs(OUT/'Figure.3.pdf', reg))
res['kept_cde'] = dict(source=len(a), candidate=len(b), missing=len(a-b), extra=len(b-a), examples=sorted(a ^ b)[:6])
top = objs(OUT/'Figure.3.pdf', (0, 0, 453.6, 165.9))
res['top_row_letters'] = sorted(o[1] for o in top if o[0] == 'chars' and o[1] in 'ab' and o[3] < 6)
v = r3['values']
o = pd.read_csv(ROOT/'notes/2026-10-10-recurrence-utr-ortholog/ortholog_3x3.tsv', sep='\t')
o = o[(o.model == 'full_assay') & np.isclose(o.threshold, .02) & o.tested]
res['a_rows_ok'] = f"{sum(v[f'a {r.contrast} human {r.human_direction} mouse {r.mouse_direction}']['obs'] == r.observed and abs(v[f'a {r.contrast} human {r.human_direction} mouse {r.mouse_direction}']['q8']-r.q_eight_per_threshold_model) < 1e-12 for _, r in o.iterrows())}/{len(o)}"
s = pd.read_csv(ROOT/'notes/2026-10-10-ortholog-threshold-scan/ortholog_3x3_threshold_scan.tsv', sep='\t')
s = s[(s.model == 'full_assay') & s.tested & (s.analytic_expected >= 1)]
pairs = [('proximal', 'proximal'), ('distal', 'distal'), ('distal', 'proximal')]
s = s[[(h, m) in pairs for h, m in zip(s.human_direction, s.mouse_direction)]]
ok = 0
for _, r in s.iterrows():
    k = f'b {r.contrast} human {r.human_direction} mouse {r.mouse_direction} t={r.threshold:.3f}'
    ok += k in v and v[k]['obs'] == r.observed and abs(v[k]['fold']-r.observed/r.analytic_expected) < 1e-9 and abs(v[k]['q72']-r.q_scan_72_per_model) < 1e-12
res['b_points_ok'] = f'{ok}/{len(s)}'
sig = s[s.q_scan_72_per_model < .05].groupby(['contrast', 'human_direction', 'mouse_direction']).size()
res['b_significant_points'] = {' '.join(k): int(n) for k, n in sig.items()}

# c: event-class stratification, recomputed from the analysis table
st = pd.read_csv(ROOT/'notes/2026-10-10-utr-event-stratification/overlap/primary_proximal_overlap.tsv', sep='\t')
st = st[st.model == 'full_assay']
okc = []
for cls in ('both_tandem', 'both_ALE'):
    r2 = st[(st.pair_class == cls) & np.isclose(st.threshold, .02)].iloc[0]
    okc.append(v[f'c forest {cls}']['obs02'] == r2.observed and abs(v[f'c forest {cls}']['q4']-r2.q_primary4) < 1e-12)
    for _, r in st[(st.pair_class == cls) & (st.analytic_expected >= 1)].iterrows():
        k = f'c scan {cls} t={r.threshold:.3f}'; okc.append(k in v and v[k]['obs'] == r.observed and abs(v[k]['q36']-r.q_scan36) < 1e-12)
res['c_ok'] = f'{sum(okc)}/{len(okc)}'
# e: usage and tiles recomputed from the bundle tables
sys.path.insert(0, str(ROOT/'notes/2026-10-06-supplement-polyapipe/tools'))
from supdraw import Data
B = Data()
hp = B.df('F4E', 'axis_extensions/homolog_site_pairs.tsv.gz'); hp = hp[(hp.human_gene_id == 'ENSG00000196923') & np.isclose(hp.tolerance, 25)]
e = v['e']; oke = []
for route, pair, side in (('GSE118258', 'H1_C', 'h'), ('CIP', 'H1_C', 'm'), ('WOT', 'H1_W', 'm')):
    col = 'human_id' if side == 'h' else 'mouse_id'
    g = hp[(hp.pair == pair) & (hp.design == 'qualified_P_S')].copy(); g['pos'] = g[col].str.split('|').str[-1].astype(int)
    g = g.sort_values('pos', ascending=False)   # minus strand: proximal first
    pr, di = g.iloc[0], g.iloc[1]
    oke.append(abs(e['usage'][f'{route}|P']['distal_pct']-100*di[f'target_{side}']/(pr[f'target_{side}']+di[f'target_{side}'])) < 1e-9)
mem = B.df('F4E', 'cross_trace/member_endpoint_effects.tsv.gz'); mem = mem[(mem.human_gene_id == 'ENSG00000196923') & (mem.setting == 'supported50_p08')]
for c in e['cells']:
    if c['route'] == '178325_0809' or c['value'] is None: continue
    m = mem[(mem.contrast == c['contrast']) & (mem.route == c['route'])].iloc[0]; oke.append(abs(m.delta-c['value']) < 1e-12)
le = v['e last exons']; oke.append(le['ENST00000355572']['start0'] == e['human_sites'][0] and abs(le['ENST00000355841']['start0']-e['human_sites'][1]) == 0)
res['e_ok'] = f'{sum(oke)}/{len(oke)}'
res['set_copies'] = {k: sha(OUT/k) == d['sha256'] for k, d in rec['set'].items()}
res['checks'] = {c: len(r3['checks'][c]) for c in ('outside', 'text_overlaps', 'letter_zone', 'colour_unnamed', 'line_aspect_bad')}
(OUT/'verify.json').write_text(json.dumps(res, indent=1, default=str)); print(json.dumps(res, indent=1, default=str))
