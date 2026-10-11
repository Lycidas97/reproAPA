"""Write build_figure5_riok3.py: the polyApipe Figure 5 builder with panel F pinned to Riok3 (rows from derive.py)."""
from pathlib import Path
SRC = Path(__file__).resolve().parents[1]/'2026-10-06-main-figures-polyapipe/builders/figure5/build.py'
lines = SRC.read_text().split('\n')
start = next(i for i, l in enumerate(lines) if l.startswith('# ------------------------------------------------------------------ F'))
end = next(i for i, l in enumerate(lines) if l.startswith("ver['F'] = Fd"))
head = '\n'.join(lines[:start]).replace("HERE = Path(__file__).resolve().parent\nsys.path.insert(0, str(HERE.parent/'tools-local'))",
    "HERE = Path(__file__).resolve().parent\nsys.path.insert(0, str(HERE.parent/'2026-10-06-main-figures-polyapipe/builders/tools-local'))", 1)
assert '2026-10-06-main-figures-polyapipe/builders/tools-local' in head
tail = '\n'.join(lines[end+1:])
F = r'''# ------------------------------------------------------------------ F  (2026-10-07: Riok3, cross-workflow example)
# Riok3 is not a low-shift screen hit (|ΔDWUI| > .01 in both routes with polyApipe); it is pinned here because its terminal-family
# redistribution and interval inclusion loss are supported with both PAS workflows. Its native, family and inclusion rows are not
# in this bundle's example/geometry tables (those cover the low-shift cohort only) and come from derive.py, which reproduces the
# bundle's Trip12 rows exactly with the same code paths (receipt.json, checks).
import hashlib as _hl
DER = Path(__file__).resolve().parent/'tables'
_dsha = {p.name: _hl.sha256(p.read_bytes()).hexdigest() for p in sorted(DER.glob('*.tsv'))}
GIDF, GSYMF = 'ENSMUSG00000024404', 'Riok3'
IV889 = (12288586, 12289475)
R.chosen['low_shift@F5F'] = dict(mode='pin', chosen=dict(gene_id=GIDF, gene_symbol=GSYMF), note='cross-workflow example; not a low-shift screen hit')
R.derive('F5F', f'Riok3 native endpoint, terminal-family decomposition and paired inclusion from derived tables {json.dumps(_dsha)} (notes/2026-10-07-riok3-example/derive.py; Trip12 reproduction exact)')
R.derive('F5F', 'shared interval fixed to chr18:[12288586,12289475) (889 nt), the interval supported by both PAS workflows; the polyApipe strict interval [12288550,12289754) has the same complete-retention mask [0,1]')
use('F5F', 'examples/example_endpoint_effects.tsv.gz', lambda r: False); use('F5F', 'examples/inclusion_paired_libraries.tsv.gz', lambda r: False)
use('F5F', 'geometry/axis-decomposition.tsv', lambda r: False); use('F5F', 'candidate_screen/recurrent_region_pairs.tsv.gz', lambda r: False)
use('F5F', 'panel_statistics/inclusion-paired-means.tsv', lambda r: False)
ee = tsv(DER/'example_endpoint_effects.tsv')
ip = [x for x in tsv(DER/'inclusion_paired_libraries.tsv') if (int(x['start0']), int(x['end0'])) == IV889]
pmF = [x for x in tsv(DER/'inclusion-paired-means.tsv') if x['cutoff'] in ('20', '20.0') and f'|{IV889[0]}-{IV889[1]}|' in x['feature_id']]
famF = {r['route']: r for r in tsv(DER/'axis-decomposition.tsv')}
inv = use('F5F', 'site_statistics/site_factor_inventory.tsv.gz', lambda r: r['gene_id'] == GIDF)
PSC_SEEDS = {r['seed'] for r in R.rows('F5F', 'external_mirna/external_seed_expression.tsv',
                                       lambda r: r['species'] == 'mouse' and r['expression_group'] == 'PSC')}
_pex = {}
for r_ in R.rows('F5F', 'rbp_expression/state_expression.tsv', lambda r: r['species'] == 'mouse' and r['state'] == 'P' and r['route'] in ('CIP', 'WOT')):
    _pex.setdefault(r_['factor'], set()).add(r_['route']) if r_['expressed_any_2CPM'] == 'True' else None
P_RBP = {f for f, rs in _pex.items() if rs == {'CIP', 'WOT'}}
RUG_ROWS = [
    ('miRNA', 'miRNA', lambda x: x['factor'] in PSC_SEEDS, 412.0, INK, 0.6),
    ('RBP_direct', 'RBP', lambda x: x['provenance'] == 'RBP', 421.5, '#8C6D31', 0.6),
    ('RBP_transfer', 'RBP', lambda x: x['provenance'] == 'RBP_transfer' and any(n in P_RBP for n in x['names'].split(';')), 431.0, '#C9B68F', 0.5)]
R.derive('F5F', 'transferred-motif RBP row: provenance RBP_transfer and any named RBP expressed >= 2 CPM in P of both CIP and WOT (rbp_expression/state_expression)')
R.derive('F5F', 'site rug limited to PSC-enriched miRNA seeds (external_seed_expression, mouse, expression_group PSC) and same-species RBP motifs (provenance RBP)')
RUG_LABEL = {'miRNA': 'PSC miRNA', 'RBP_direct': 'Direct RBP', 'RBP_transfer': 'Transferred RBP'}
RUG_SHORT = {'RBP_transfer': 'Transf. RBP'}
RUG_FS = 7.5
IVS = [IV889]
inany = lambda a, b: any(a >= s_ and b <= e_ for s_, e_ in IVS)
R.derive('F5F', 'family PAS usage renormalised within the terminal family of geometry/axis-decomposition (GAP-D06)')
t.text_top(24.9, 389.6, GSYMF, 9, font='BolIta')
_x = 92.0
t.line([_x, _x+10.0], [394.0, 394.0], lw=1.0, color='#8A959C', ls=(0, (3, 1.6))); t.text_top(_x+12.5, 390.0, 'M', 8); _x += 12.5+textwidth('M', 8)+9.0
t.line([_x, _x+10.0], [394.0, 394.0], lw=1.2, color=INK); t.text_top(_x+12.5, 390.0, 'P', 8); _x += 12.5+textwidth('P', 8)+12.0
_sl = 'Shared interval'
t.rect(_x, 391.4, 7.0, 6.0, fc='#E4DFEE'); t.text_top(_x+9.5, 390.0, _sl, 8); _x += 9.5+textwidth(_sl, 8)+12.0
for _c, _n in [('#CF4C00', 'Guo2019–CIP'), ('#207450', 'Schiebinger2019–2i')]:
    t.rect(_x, 391.4, 6.0, 6.0, fc=_c); t.text_top(_x+8.5, 390.0, _n, 8); _x += 8.5+textwidth(_n, 8)+12.0
assert _x < 452.0, _x
nsite = {}
for route in ('CIP', 'WOT'):
    rr = [r for r in inv if r['contrast'] == 'MP' and r['route'] == route and r['strict_changed'] == 'True' and inany(int(r['start0']), int(r['end0']))]
    nsite[route] = {c: sum(1 for r in rr if r['category'] == c) for c in ('miRNA', 'RBP')}
    nsite[route]['coords'] = {c: len({(r['start0'], r['end0']) for r in rr if r['category'] == c}) for c in ('miRNA', 'RBP')}
Fd = {}
for k, (route, title, colr) in enumerate([('CIP', 'Guo2019–CIP', '#CF4C00'), ('WOT', 'Schiebinger2019–2i', '#207450')]):
    dx = 230.2*k
    r = [r for r in ee if r['route'] == route]
    if len(r) != 1 or route not in famF:
        R.unavailable('F5F', f'{route}: native P–M usage or terminal family not available; column left empty'); continue
    r = r[0]
    ids = json.loads(r['PAS_ids']); pos = [pas_pos(i_) for i_ in ids]; strand = pas_strand(ids[0])
    tu = json.loads(r['target_usage']); ru = json.loads(r['reference_usage'])
    fs_, fe_ = int(float(famF[route]['family_start0'])), int(float(famF[route]['family_end0']))
    plusF = strand == '+'
    fam = sorted([i for i, p in enumerate(pos) if fs_ <= p <= fe_], key=lambda i: pos[i] if plusF else -pos[i])
    if len(fam) != int(float(famF[route]['terminal_native_pas'])):
        R.flag('F5F', 'consistency', f'{route}: {len(fam)} PAS inside family bounds != terminal_native_pas {famF[route]["terminal_native_pas"]}')
    if not fam:
        R.unavailable('F5F', f'{route}: no PAS inside the terminal family'); continue
    base = pos[fam[0]]
    DF = (lambda p: p-base) if plusF else (lambda p: base-p)
    FMAX = 2200.0
    R.axis('F5F', f'{route} complete family distance', [DF(pos[j]) for j in fam], 0, FMAX)
    XF = lambda nt, dx=dx: 34.6+dx+nt*0.0365*2500/FMAX
    YF = lambda v: 499.2-v*0.528
    t.rect(6.0+dx, 402.0, 212.0, 3.0, fc=colr)
    for iv in IVS:
        xa_, xb_ = sorted((XF(DF(iv[0])), XF(DF(iv[1]))))
        t.rect(xa_, YF(100), xb_-xa_, YF(0)-YF(100), fc='#E4DFEE')
    xl_ = min(XF(DF(v_)) for iv in IVS for v_ in iv); xr_ = max(XF(DF(v_)) for iv in IVS for v_ in iv)
    HDR = [f"MPRO {float(r['MPRO']):.3f}, n = {r['n_pairs']}", 'DWUI M → P',
           f"{float(r['reference_DWUI']):.3f} → {float(r['target_DWUI']):.3f}"]
    for k_, h_ in enumerate(HDR):
        t.text_top(216.0+dx, 408.8+9.6*k_, h_, 8, ha='right')
    OBST = [(216.0+dx-textwidth(h_, 8), 408.8+9.6*k_, 408.8+9.6*k_+8.0) for k_, h_ in enumerate(HDR)]
    OBST += [(153.2+dx-textwidth('100', 8), 440.2, 448.2)]
    rug = {}
    for key, cat, keep_, yr, cc, lw_ in RUG_ROWS:
        rr = [x for x in inv if x['contrast'] == 'MP' and x['route'] == route and x['strict_changed'] == 'True' and x['category'] == cat
              and inany(int(x['start0']), int(x['end0'])) and keep_(x)]
        ds = sorted({DF((int(x['start0'])+int(x['end0']))/2) for x in rr})
        if not ds: R.unavailable('F5F', f'{route} {key}: no site of this row inside the interval')
        for d in ds:
            t.line([XF(d)]*2, [yr, yr+4.2], lw=lw_, color=cc)
        lab_ = f'{RUG_LABEL[key]} ({len(ds)})'
        if xl_-2.0-textwidth(lab_, RUG_FS) > 30.6+dx+6.0:
            t.text_top(xl_-2.0, yr-2.1, lab_, RUG_FS, ha='right', color=MUTED)
        else:
            y0_, y1_ = yr-2.1, yr-2.1+RUG_FS
            lim_ = min([x_ for x_, a_, b_ in OBST if a_ < y1_+1.0 and b_ > y0_-1.0] + [218.0+dx]) - 4.0
            if xr_+2.5+textwidth(lab_, RUG_FS) >= lim_ and key in RUG_SHORT:
                lab_ = f'{RUG_SHORT[key]} ({len(ds)})'
            t.text_top(xr_+2.5, yr-2.1, lab_, RUG_FS, ha='left', color=MUTED)
            assert xr_+2.5+textwidth(lab_, RUG_FS) < lim_, lab_
        rug[key] = dict(rule=RUG_LABEL[key], annotations=len(rr), distinct_sites=len(ds), factors=sorted({x['factor'] for x in rr}),
                        named=sorted({n for x in rr for n in x['names'].split(';')}) if key == 'RBP_transfer' else None)
    for usage, style, col, lw in [(ru, (0, (3, 1.6)), '#8A959C', 1.0), (tu, '-', INK, 1.2)]:
        u = np.array([usage[i] for i in fam]); u = u/u.sum(); cur = 1.0; pts = [(XF(0), YF(100))]
        for i, j in enumerate(fam):
            d = DF(pos[j]); pts.append((XF(d), YF(cur*100))); cur -= u[i]; pts.append((XF(d), YF(max(cur, 0)*100)))
        pts.append((XF(FMAX), YF(max(cur, 0)*100)))
        t.line([p[0] for p in pts], [p[1] for p in pts], lw=lw, color=col, ls=style)
    t.line([30.6+dx, 30.6+dx], [YF(100), YF(0)], lw=0.5, color=AXC); t.line([30.6+dx, XF(FMAX)], [501.6, 501.6], lw=0.5, color=AXC)
    for v in (100, 50, 0):
        t.line([28.6+dx, 30.6+dx], [YF(v)]*2, lw=0.5, color=AXC); t.text_top(27.3+dx, YF(v)-4.0, str(v), 8, ha='right')
    for v in (0, 1000, 2000):
        t.line([XF(v)]*2, [501.6, 503.6], lw=0.5, color=AXC); t.text_top(XF(v), 504.0, str(v), 8, ha='center')
    t.text_top(XF(1000), 513.0, 'Distance (nt)', 8, ha='center')
    t.vtext(9.0+dx, 470.0, 'Family PAS usage (%)', 7.5)
    lo_, hi_ = 0.0, 1.0
    YI = lambda v: 496.2-(min(max(v, lo_), hi_)-lo_)/(hi_-lo_)*52.0
    xm, xp = 168.0+dx, 205.8+dx
    lib = [x for x in ip if x['route'] == route and x['cutoff'] in ('20', '20.0')]
    refs = np.array([float(x['reference_inclusion']) for x in lib]); tars = np.array([float(x['target_inclusion']) for x in lib])
    pm_ = [m for m in pmF if m['route'] == route]
    incl = []
    if not lib or len(pm_) != 1:
        R.unavailable('F5F', f'{route}: inclusion pairs / paired mean not available')
    else:
        pm_ = pm_[0]
        R.axis('F5F', f'{route} interval inclusion', list(refs)+list(tars), lo_, hi_)
        for a_, b_ in zip(refs, tars):
            t.line([xm, xp], [YI(a_), YI(b_)], lw=0.4, color=colr, alpha=0.55)
        t.line([xm, xp], [YI(float(pm_['reference_mean'])), YI(float(pm_['target_mean']))], lw=1.6, color=colr)
        t.text_top(xm, 504.0, 'M', 8, ha='center'); t.text_top(xp, 504.0, 'P', 8, ha='center')
        if int(float(pm_['n'])) != len(lib): R.flag('F5F', 'consistency', f'{route}: paired-mean n {pm_["n"]} != {len(lib)} libraries')
        incl.append(dict(interval=1, start0=IV889[0], end0=IV889[1], length_nt=IV889[1]-IV889[0], n=int(float(pm_['n'])),
                         source_feature_id=pm_['feature_id'], reference_mean=float(pm_['reference_mean']), target_mean=float(pm_['target_mean']),
                         delta_I_pp=float(pm_['delta_mean'])*100))
    t.line([156.5+dx]*2, [YI(1.0), YI(lo_)], lw=0.5, color=AXC)
    for v, s in [(1.0, '100'), (0.5, '50'), (0.0, '0')]:
        t.line([154.5+dx, 156.5+dx], [YI(v)]*2, lw=0.5, color=AXC); t.text_top(153.2+dx, YI(v)-4.0, s, 8, ha='right')
    if incl:
        lab = 'ΔI ' + ' / '.join(num(x['delta_I_pp'], 1) for x in incl) + ' pp; n=' + '/'.join(str(x['n']) for x in incl)
        fs_ = 8 if textwidth(lab, 8) < 216.0-112.0 else 7.5
        t.text_top(216.0+dx, 515.2, lab, fs_, ha='right')
        assert 216.0-textwidth(lab, fs_) > 108.0, lab
    t.vtext(135.0+dx, 470.0, 'Interval\ninclusion (%)', 7.5)
    Fd[route] = dict(MPRO=float(r['MPRO']), n_native_pairs=int(r['n_pairs']), DWUI_M=float(r['reference_DWUI']), DWUI_P=float(r['target_DWUI']),
                     delta_DWUI=float(r['delta_DWUI']), family_shift_nt=float(famF[route]['terminal_shift_nt']),
                     within_family_DWUI=float(famF[route]['within_family_DWUI']), nonterminal_DWUI=float(famF[route]['nonterminal_DWUI']),
                     intervals=incl, sites=nsite[route], site_rug=rug)
ver['F'] = Fd'''
out = Path(__file__).resolve().parent/'build_figure5_riok3.py'
out.write_text(head+'\n'+F+'\n'+tail)
print(out, len(out.read_text().split('\n')))
