"""Supplementary Figure S13 candidate (v2, 2026-10-09). The study-pair sharing panel of old main Figure 4a is REDRAWN with
supdraw (supplement house style: sizes, LW/MS/MEW tokens, named colours, layout checks) into the rectangle that old S13a+b
occupied (333.3 x 112 pt), so the figure keeps the original S13 geometry: row 1 = new a | old c (now b), band = old d, e, f
(now c, d, e), each moved 0 pt and only relettered. Study colours are explained by the shared study key under panel d.
Old S13a (threshold/layer sensitivity) -> Table S5b block S5b.5 only (user, 2026-10-09); old S13b is contained in the new a.
Data: cross_sharing/sharing_summary.tsv (same SHA in the main and supplement bundles). Usage: python3 compose_supp.py"""
import json, math
from pathlib import Path
import pandas as pd
HERE = Path(__file__).resolve().parent
g = {'__file__': str(HERE/'compose_main.py')}
exec((HERE/'compose_main.py').read_text().split("if __name__ == '__main__':")[0], g)
FIG, ROOT, crop, place, new_page, finish, W, letters, letter, supdraw, sha = (g[k] for k in (
    'FIG', 'ROOT', 'crop', 'place', 'new_page', 'finish', 'W', 'letters', 'letter', 'supdraw', 'sha'))
from supdraw import Lin, PROX, DIST, INK, MUTED, AXC, NULLBAR, textwidth
from v5draw import BAND
S13 = FIG/'SuppFig.S13.pdf'; H = 404.3; L13 = letters(S13)
TAB = ROOT/'notes/2026-10-06-supplement-polyapipe/contract-r2/bundle/tables/cross_sharing/sharing_summary.tsv'
dst = new_page(H)
lb = lambda L: (L13[L]['box'][0]-0.6, L13[L]['box'][1]-0.6, L13[L]['box'][2]+0.6, L13[L]['box'][3]+0.6)
KEEP = [((340.0, 0, W, 120.0), 'c', 'b'), ((0, 122.0, 160.0, 296.0), 'd', 'c'), ((160.0, 122.0, W, H), 'e', 'd'), ((0, 296.0, 160.0, H), 'f', 'e')]
for R, old, new in KEEP:
    place(dst, H, crop(S13, H, R, [lb(old)]), R, R[0], R[1])
F = supdraw.Fig(H, 'SuppFig.S13')
for R, old, new in KEEP:
    letter(F.t, new, L13[old]['ox'], L13[old]['base'])
# ---------------------------------------------------------------- a  same-direction sharing by study pair (redrawn)
x, y, w, h = F.panel('a', 4.5, 4.5, 349.1-11.3-4.5, 112.0)
ss = pd.read_csv(TAB, sep='\t')
PAIRS = [('GSE118258', 'GSE103221', 'GSE118258', 'CIP'), ('GSE118258', 'GSE115943', 'GSE118258', 'WOT'),
         ('GSE178325', 'GSE103221', '178325_0618', 'CIP'), ('GSE178325', 'GSE115943', '178325_0618', 'WOT')]
ROW0, PITCH, PGAP = y+21.5, 8.2, 2.0
ry = lambda k: ROW0+k*PITCH+(k//2)*PGAP
for i, (hs, ms, hr, mr) in enumerate(PAIRS):
    if i % 2 == 0:
        F.rect(x, ry(2*i)-4.6, w, ry(2*i+1)-ry(2*i)+9.2, fc=BAND)
    F.study_pair(x+0.4, (ry(2*i)+ry(2*i+1))/2, hr, mr)
RLO, RHI = 0.6, 2.0
Av = {}
for bi, (metric, head) in enumerate((('WPP', 'DWUI'), ('RWUI', 'RWUI'))):
    b0 = x+24.0+bi*159.3
    lx = lambda v, b0=b0: b0+(math.log(v)-math.log(RLO))/(math.log(RHI)-math.log(RLO))*58.0
    sc = Lin(RLO, RHI, b0, b0+58.0)
    xq, xb, bw = b0+58.0+6.0+18.0, b0+58.0+6.0+36.0+4.0, 24.0
    F.text(b0+72.0, y+0.6, f'{head}, same direction', 7.5, ha='center')
    F.text(lx(1.0), y+10.6, 'Obs./null', 7.0, ha='center', color=MUTED)
    F.text(xq, y+10.6, 'q8', 7.0, ha='center', color=MUTED)
    F.text(xb+bw/2+6.0, y+10.6, 'Shared', 7.0, ha='center', color=MUTED)
    F.line([lx(1.0)]*2, [ry(0)-4.6, ry(7)+4.6], lw='ref', color=AXC, ls='ref')
    k = 0
    for hs, ms, hr, mr in PAIRS:
        for thr in (0.02, 0.05):
            q = ss[(ss.metric == metric) & (ss.human_study == hs) & (ss.mouse_study == ms) & ((ss.threshold-thr).abs() < 1e-9)].set_index('test')
            r = q.loc['same_direction']; mean = float(r.null_mean)
            lo, hi, enr = float(r.null_p025)/mean, float(r.null_p975)/mean, float(r.enrichment)
            assert RLO <= lo and hi <= RHI and RLO <= enr <= RHI and abs(enr-float(r.observed)/mean) < 1e-9
            yy = ry(k); k += 1
            F.line([lx(lo), lx(hi)], [yy, yy], lw='bold', color=NULLBAR)
            sig = float(r.q_primary_8) < 0.05; dag = float(r.q_joint_primary_24) < 0.05
            F.marker(lx(enr), yy, 'o' if thr == 0.02 else 'D', 'mark' if thr == 0.02 else 'small', fc=INK if sig else 'white', ec=INK)
            F.text(xq, yy-3.4, f'{float(r.q_primary_8):.4f}'.lstrip('0')+(' †' if dag else ''), 7.0, ha='center')
            neg, pos = int(q.loc['shared_negative'].observed), int(q.loc['shared_positive'].observed)
            assert neg+pos == int(r.observed)
            F.rect(xb, yy-2.1, bw*neg/(neg+pos), 4.2, fc=PROX); F.rect(xb+bw*neg/(neg+pos), yy-2.1, bw*pos/(neg+pos), 4.2, fc=DIST)
            F.text(xb+bw+1.6, yy-3.4, str(int(r.observed)), 7.0, color=MUTED)
            Av[f'{head}|{hr}|{mr}|{thr}'] = dict(observed=int(r.observed), null_mean=mean, ratio=enr, null=[lo, hi], q8=float(r.q_primary_8),
                                                 q24=float(r.q_joint_primary_24), proximal=neg, distal=pos)
    F.xaxis(sc, ry(7)+5.2, [], size=7.0)
    for v, lab in ((0.7, '0.7'), (1.0, '1'), (1.4, '1.4'), (2.0, '2')):
        F.line([lx(v)]*2, [ry(7)+5.2, ry(7)+7.2], lw='axis', color=AXC); F.text(lx(v), ry(7)+8.2, lab, 7.0, ha='center')
items = [('marker', ('o', INK), '0.02'), ('marker', ('D', INK), '0.05'), ('marker', ('o', 'white', INK), 'q8 ≥ .05'),
         ('block', PROX, 'Proximal'), ('block', DIST, 'Distal'), ('line', (NULLBAR, None, 'bold'), 'Null 95%')]
tail = '† q24 < .05'
wk = textwidth('Threshold', 7.0)+5.0+F.key_width(items, 7.0, 6.0)+6.0+textwidth(tail, 7.0)
kx = x+w/2-wk/2; ky = y+h-7.4
kx += F.text(kx, ky, 'Threshold', 7.0)+5.0
kx = F.key_row(kx, ky, items, size=7.0, gap=6.0)+6.0
F.text(kx, ky, tail, 7.0)
F.values['a'] = Av
chk = F.check()
print('checks: outside', len(chk['outside']), 'overlaps', len(chk['text_overlaps']), 'letter-zone', len(chk['letter_zone']), 'colour-unnamed', chk.get('colour_unnamed'))
if chk['text_overlaps']: print(chk['text_overlaps'][:6])
finish(dst, F.t, 'SuppFig.S13', H)
rec = dict(version=2, sources={'SuppFig.S13.pdf': sha(S13), 'sharing_summary.tsv': sha(TAB)}, page_pt=[W, H],
           new_a=dict(rect=[x, y, w, h], values=Av), kept=[dict(region=R, old=o, new=n) for R, o, n in KEEP],
           checks=dict(outside=len(chk['outside']), text_overlaps=len(chk['text_overlaps']), letter_zone=len(chk['letter_zone'])),
           dropped='old S13a -> Table S5b block S5b.5; old S13b contained in the new a; old main Figure 4a values redrawn here')
rc = json.loads((HERE/'out/receipt.json').read_text()); rc['SuppFig.S13'] = {**g['RECEIPT']['SuppFig.S13'], **rec}
(HERE/'out/receipt.json').write_text(json.dumps(rc, indent=1, default=float, ensure_ascii=False))
