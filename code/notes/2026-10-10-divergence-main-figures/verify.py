"""Independent checks for the 2026-10-10 divergence candidates: Fig/ untouched, kept regions pixel-identical after the
recorded shift, key plotted values recomputed from the source tables without compose.py helpers."""
import json, hashlib, subprocess, tempfile
from pathlib import Path
import numpy as np, pandas as pd
from PIL import Image
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]; OUT = HERE/'out'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
rec = json.loads((OUT/'receipt.json').read_text()); res = {}
before = dict(l.split()[::-1] for l in (OUT/'fig-sha-before.txt').read_text().splitlines())
res['fig_unchanged'] = {k: before[k] == sha(HERE/k) for k in before}
DPI = 150; S = DPI/72
def render(pdf):
    t = Path(tempfile.mkdtemp())/'r'
    subprocess.run(['pdftoppm', '-r', str(DPI), '-png', '-singlefile', str(pdf), str(t)], check=True, capture_output=True)
    return np.asarray(Image.open(str(t)+'.png').convert('L')).astype(int)
def cmp(src, cand, region, dy):
    a, b = render(src), render(cand); x0, y0, x1, y1 = region
    # shrink by 1 px at region edges to avoid antialias of neighbours
    ra = a[int(y0*S)+1:int(y1*S)-1, int(x0*S)+1:int(x1*S)-1]
    yb = int(round((y0+dy)*S))+1
    rb = b[yb:yb+ra.shape[0], int(x0*S)+1:int(x0*S)+1+ra.shape[1]]
    n = min(ra.shape[0], rb.shape[0]); d = np.abs(ra[:n]-rb[:n])
    return dict(pixels=int(d.size), diff_gt32=int((d > 32).sum()), max=int(d.max()))
r3, r4 = rec['Figure.3'], rec['Figure.4']
import pdfplumber
def objs(pdf, region, dy=0.0):
    pg = pdfplumber.open(pdf).pages[0]; x0, y0, x1, y1 = region; out = []
    for kind in ('chars', 'rects', 'lines', 'curves'):
        for o in getattr(pg, kind):
            if o['x1']-o['x0'] > 400 and o['bottom']-o['top'] > 400: continue
            cx, cy = (o['x0']+o['x1'])/2, (o['top']+o['bottom'])/2-dy
            if x0 <= cx <= x1 and y0 <= cy <= y1:
                out.append((kind, o.get('text', ''), round(o['x0'], 1), round(o['top']-dy, 1), round(o['x1'], 1), round(o['bottom']-dy, 1)))
    return sorted(out)
def cmp(src, cand, region, dy, ignore=()):
    inside = lambda o: any(bx[0] <= (o[2]+o[4])/2 <= bx[2] and bx[1] <= (o[3]+o[5])/2 <= bx[3] for bx in ignore)
    a, b = [o for o in objs(src, region) if not inside(o)], [o for o in objs(cand, region, dy) if not inside(o)]
    sa, sb = set(a), set(b)
    return dict(source_objects=len(a), candidate_objects=len(b), missing=len(sa-sb), extra=len(sb-sa), examples=sorted(sa ^ sb)[:6])
ign = r3.get('e_relabelled', {}).get('cleared', [])
res['Figure.3_kept'] = [cmp(ROOT/r3['source'], ROOT/r3['pdf'], k['region'], k['dy'], ign if k['dy'] else ()) for k in r3['kept']]
res['Figure.3_e_relabelled'] = r3.get('e_relabelled')
res['Figure.4_kept'] = [cmp(ROOT/r4['source'], ROOT/r4['pdf'], k['region'], 0) for k in r4['kept']]
# values
v3, v4 = r3['values'], r4['values']
dr = pd.read_csv(ROOT/'notes/2026-10-10-directional-recurrence/directional_recurrence.tsv', sep='\t')
ok = []
for _, r in dr[np.isclose(dr.threshold, .02)].iterrows():
    v = v3[f"3a {r.species} {r.contrast} {r.direction}"]; ok.append(v['obs02'] == r.observed and abs(v['fold02']-r.observed/r.null_mean) < 1e-9)
res['3a_rows_ok'] = f'{sum(ok)}/{len(ok)}'
# 3c: independent recompute from qualified date effects
q = pd.read_csv(ROOT/'notes/2026-10-06-supplement-polyapipe/contract-r2/bundle/tables/global/qualified_date_effects.tsv.gz', sep='\t')
b = pd.read_csv(ROOT/'notes/2026-10-06-supplement-polyapipe/contract-r2/bundle/tables/global/background_gene_membership.tsv.gz', sep='\t')
q = q[(q.setting == 'supported50_p08') & (q.contrast == 'P_S')]; b = b[(b.setting == 'supported50_p08') & (b.contrast == 'P_S') & (b.background == 'common_one_to_one_ortholog_all_settings')]
arm = {}
for r in ['GSE118258', '178325_0618', '178325_0809', 'CIP', 'WOT']:
    x = q[q.route == r].groupby('gene_id').APA_delta.mean(); x = x[x.index.isin(b[b.route == r].gene_id)]
    arm[r] = (x.clip(lower=0).mean(), x.clip(upper=0).mean(), (x <= -.02).mean(), (x >= .02).mean())
hp = (arm['GSE118258'][0]+(arm['178325_0618'][0]+arm['178325_0809'][0])/2)/2; hn = (arm['GSE118258'][1]+(arm['178325_0618'][1]+arm['178325_0809'][1])/2)/2
mp = (arm['CIP'][0]+arm['WOT'][0])/2; mn = (arm['CIP'][1]+arm['WOT'][1])/2
res['3c_diff'] = dict(recomputed_proximal_share=(hn-mn)/((hp+hn)-(mp+mn)), plotted=v3['3c human-mouse']['proximal_share'])
res['3c_fractions_ok'] = all(abs(arm[r][2]-v3[f'3c {r}']['fp']) < 1e-12 and abs(arm[r][3]-v3[f'3c {r}']['fd']) < 1e-12 for r in arm)
t = pd.read_csv(ROOT/'notes/2026-10-10-recurrence-utr-ortholog/tandem_directional_recurrence.tsv', sep='\t')
t = t[(t.contrast == 'P_S') & np.isclose(t.threshold, .02)]
res['3d_counts'] = {f'{r.species} {r.direction}': [int(r.observed), v3[f'3d {r.species} {r.direction}']['obs']] for _, r in t.iterrows()}
o = pd.read_csv(ROOT/'notes/2026-10-10-recurrence-utr-ortholog/ortholog_3x3.tsv', sep='\t')
o = o[(o.model == 'full_assay') & np.isclose(o.threshold, .02) & o.tested]
res['4b_rows_ok'] = all(v4[f'4b {r.contrast} human {r.human_direction} mouse {r.mouse_direction}']['obs'] == r.observed and
                        abs(v4[f'4b {r.contrast} human {r.human_direction} mouse {r.mouse_direction}']['q8']-r.q_eight_per_threshold_model) < 1e-12 for _, r in o.iterrows())
res['new_content_checks'] = {k: {c: len(rec[k]['checks'][c]) for c in ('outside', 'text_overlaps', 'letter_zone', 'colour_unnamed', 'line_aspect_bad')} for k in ('Figure.3', 'Figure.4')}
(OUT/'verify.json').write_text(json.dumps(res, indent=1, default=float)); print(json.dumps(res, indent=1, default=float))
