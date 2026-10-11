"""Supplementary Figure D17 (printed number from supdraw.SNUM) — the shared proximal candidates (R4), polyApipe v2 primary.

a  the 23 shared P–S candidates (supported layer, |ΔDWUI| >= .02): P–S ΔDWUI in every documented-source route and the
   species means used for admission
b  the 20 shared P–M candidates: paired P–M ΔDWUI per study (same-library estimator) and the species means
PDLIM7 (Figure 4e) in bold. Grey cells (MISSING): not measured (not zero).
"""
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools'))
from supdraw import *   # noqa

FID = 'D17'
D = Data(); OUT = outdir(FID)
X0, RR = MARGIN, PAGE_W-MARGIN
mem = D.df('D17A', 'cross_overlap/new_shared_members.tsv.gz')
mem = mem[(mem.layer == 'supported50') & (abs(mem.threshold-0.02) < 1e-9)]
re_ = D.df('D17A', 'state_vectors/route_effects.primary.tsv.gz', usecols=['gene_id', 'route', 'contrast', 'design', 'setting', 'delta_DWUI'])
re_ = re_[(re_.contrast == 'P_S') & (re_.design == 'true_source') & (re_.setting == 'supported50_p08')].set_index(['route', 'gene_id']).delta_DWUI
pm = D.df('D17B', 'cross_overlap/paired_pm_study_effects.tsv.gz')
pm = pm[pm.layer == 'supported50'].set_index(['code', 'gene_id']).mean_delta
LIM = 0.3
NPS, NPM = int((mem.contrast == 'PS').sum()), int((mem.contrast == 'PM').sum())
rh = 9.0
NAMES = {c: STUDY[c][0] for c in STUDY}
VLW = max(textwidth(NAMES[r], 7.0) for r in ['GSE118258', '178325_0618', '178325_0809', 'GSE150311', 'GSE221739', 'CIP', 'WOT'])
VLW = max(VLW, textwidth('Species mean', 7.0))
HDR = 24.0+VLW+3.0
H1 = HDR+max(NPS, NPM)*rh
H = MARGIN*2+H1
F = Fig(H, figname(FID))
COLS_A = [('route', 'GSE118258', 'h'), ('route', '178325_0618', 'h'), ('route', '178325_0809', 'h'), ('mean', 'h_PS_mean_effect', 'h'),
          ('route', 'CIP', 'm'), ('route', 'WOT', 'm'), ('mean', 'm_PS_mean_effect', 'm')]
COLS_B = [('code', 'H1', 'GSE118258', 'h'), ('code', 'H2', 'GSE150311', 'h'), ('code', 'H4', 'GSE221739', 'h'), ('mean', 'h_PM_mean_effect', None, 'h'),
          ('code', 'C', 'CIP', 'm'), ('code', 'W', 'WOT', 'm'), ('mean', 'm_PM_mean_effect', None, 'm')]


def heat(L, x, w, rows, cols, getter, title, key=False):
    x, y, w, h = F.panel(L, x, MARGIN, w, H1)
    labw = max(textwidth(s, 7.0, 'BolIta') for s in rows.human_gene_symbol)
    hx = x+labw+4.0; cw = (x+w-hx)/len(cols)
    F.text((hx+x+w)/2, y+0.4, title, 7.5, ha='center')
    for (k0, k1, lab) in ((0, 4, 'Human'), (4, 7, 'Mouse')):
        F.text(hx+cw*(k0+k1)/2, y+11.0, lab, 7.5, ha='center')
        F.line([hx+cw*k0+1.2, hx+cw*k1-1.2], [y+21.0]*2, lw='thin', color=INK)
    vals = {}
    for k, c in enumerate(cols):
        cx_ = hx+cw*k+cw/2
        if c[0] == 'mean':
            F.vlabel(cx_, y+HDR-3.0, 'Species mean', 7.0, color=MUTED)
        else:
            r = c[2] if c[0] == 'code' else c[1]
            F.vlabel(cx_, y+HDR-3.0, NAMES[r], 7.0)
    F.line([hx+cw*4]*2, [y+HDR, y+HDR+len(rows)*rh], lw='axis', color='white')
    for i, r in enumerate(rows.itertuples()):
        yy = y+HDR+i*rh
        bold = r.human_gene_symbol.upper() == 'PDLIM7'
        F.text(hx-2.0, yy+rh/2-3.4, r.human_gene_symbol, 7.0, font='BolIta' if bold else 'RegIta', ha='right')
        vals[r.human_gene_symbol] = {}
        for k, c in enumerate(cols):
            v = getter(r, c)
            cx_ = hx+cw*k
            if v is None or v != v:
                F.rect(cx_+0.4, yy+0.4, cw-0.8, rh-0.8, fc=MISSING)
                vals[r.human_gene_symbol][str(c[1])] = None
            else:
                fc = diverge(float(v), LIM)
                F.rect(cx_+0.4, yy+0.4, cw-0.8, rh-0.8, fc=fc, ec=INK if c[0] == 'mean' else 'none', lw='thin' if c[0] == 'mean' else 0)
                if abs(v) >= 0.02: F.marker(cx_+cw/2, yy+rh/2, 'o', 'point', fc=on(fc))
                vals[r.human_gene_symbol][str(c[1])] = float(v)
    if key:
        yb = y+HDR+len(rows)*rh
        F.cbar_key((hx+x+w)/2, yb+8.0, 'ΔDWUI', '−.3', '+.3', lambda f_: diverge(-LIM+2*LIM*f_, LIM), wbar=56)
        F.key_row(None, yb+19.0, [("block", MISSING, "not measured"), ('marker', ('o', INK), '|ΔDWUI| ≥ .02')], center=(hx+x+w)/2)
    F.values[L] = vals


def ga(r, c):
    if c[0] == 'mean': return getattr(r, c[1])
    gid = r.human_gene_id if c[2] == 'h' else r.mouse_gene_id
    return re_.get((c[1], gid), None)


def gb(r, c):
    if c[0] == 'mean': return getattr(r, c[1])
    gid = r.human_gene_id if c[3] == 'h' else r.mouse_gene_id
    return pm.get((c[1], gid), None)


W2 = (RR-X0-GAP)/2
ps = mem[mem.contrast == 'PS'].sort_values('h_PS_mean_effect')
pmr = mem[mem.contrast == 'PM'].sort_values('h_PM_mean_effect')
heat('a', X0, W2, ps, COLS_A, ga, f'Shared P–S candidates ({NPS}), P–S ΔDWUI')
heat('b', X0+W2+GAP, W2, pmr, COLS_B, gb, f'Shared P–M candidates ({NPM}), paired P–M ΔDWUI', key=True)
F.notes.append('colour key under b applies to both panels; outlined cells = species mean used for admission; colour clipped at ±.3')
F.save(OUT, D)
