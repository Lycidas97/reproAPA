"""Audit all supplement candidates: layout checks from the receipts, page height, embedded fonts, colour-scale record."""
import json, subprocess, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from supdraw import ORDER, SNUM, outdir
MAXH = 190/25.4*72
# S31 (working ID A01) is drawn outside ORDER, in notes/2026-10-06-alt-catalogue-suppfig, from its own pinned tables.
EXTRA = [('A01', 31, Path(__file__).resolve().parents[2]/'2026-10-06-alt-catalogue-suppfig/final/SuppFig.S31.receipt.json')]
HEIGHT_PENDING = set()   # A01 was 202 mm (pending); the 2026-10-06 v4 redesign is 189.97 mm, so no exception is needed
CONS = '--consolidated' in sys.argv     # the 21 consolidated figures (tools/compose.py, consolidated/<plan letter>/)
if CONS:
    from compose import PRINT, PNUM, cdir
    ITEMS = [(c, PNUM[c], cdir(c)/f'SuppFig.S{PNUM[c]}.receipt.json') for c in PRINT]
else:
    ITEMS = [(f, SNUM[f], outdir(f)/f'SuppFig.S{SNUM[f]}.receipt.json') for f in ORDER]+EXTRA
rows, bad = [], []
for fid, sn, rc in ITEMS:
    if not rc.exists():
        rows.append(dict(id=fid, S=sn, status='not redrawn in this series')); continue
    r = json.loads(rc.read_text()); c = r['checks']
    pdf = rc.with_suffix('').with_suffix('.pdf')
    fonts = subprocess.run(['pdffonts', str(pdf)], capture_output=True, text=True).stdout.splitlines()[2:]
    not_emb = [f for f in fonts if f.split()[-5] != 'yes'] if fonts else []
    slack = {L: v for L, v in c['slack'].items() if max(v.values()) > 1.5}
    strips = {L: v[0] for d in ('internal_vertical_strip', 'internal_horizontal_strip') for L, v in c[d].items() if v[0] >= 8.0}
    gaps = {k: v for k, v in c['ink_gaps'].items() if v < 11.3-1.0}
    row = dict(id=fid, S=sn, height_mm=round(r['page_pt'][1]/72*25.4, 1), panels=len(r['panels']),
               outside=len(c['outside']), overlaps=len(c['text_overlaps']), letter_zone=len(c['letter_zone']),
               line_aspect_bad=len(c['line_aspect_bad']), colour_unnamed=c.get('colour_unnamed', 'not checked'), sizes=c['sizes'], fonts_not_embedded=len(not_emb),
               slack_gt_1p5=slack, strips_ge_8=strips, ink_gaps_lt_10p3=gaps,
               diverging=r.get('colour_scales', {}).get('diverging'), tables=len(r.get('tables', {})),
               external=sorted(r.get('external_sources', {})), interface_exceptions=sorted(k for k, v in r.get('tables', {}).items() if v['mode'].startswith('interface')))
    tall = r['page_pt'][1] > MAXH+0.01
    if tall and fid in HEIGHT_PENDING: row['height_over_190mm_pending_decision'] = True; tall = False
    hard = row['outside']+row['overlaps']+row['letter_zone']+row['line_aspect_bad']+row['fonts_not_embedded']+tall+len(c.get('colour_unnamed', []))
    row['hard_failures'] = int(hard)
    if CONS and 'consolidated' in r:
        cr = r['consolidated']; row['letter_map'] = cr.get('letter_map'); row['dropped_panels'] = cr.get('dropped_panels')
        row['values_identical_to_standalone'] = all(v['kept_values_identical'] for v in cr.get('values_check', {}).values()) if 'values_check' in cr else 'copied'
        row['layout_warnings'] = cr.get('layout_warnings')
    if hard: bad.append(fid)
    rows.append(row)
out = Path(__file__).resolve().parent.parent/'verification'; out.mkdir(exist_ok=True)
(out/('audit-consolidated.json' if CONS else 'audit.json')).write_text(json.dumps(rows, indent=1, ensure_ascii=False))
for r in rows:
    if 'status' in r: print(f"S{r['S']:<3}{r['id']:<5}{r['status']}"); continue
    print(f"S{r['S']:<3}{r['id']:<5}{r['height_mm']:>6} mm  panels {r['panels']}  hard {r['hard_failures']}  slack>{len(r['slack_gt_1p5'])} strips>{len(r['strips_ge_8'])} gaps<{len(r['ink_gaps_lt_10p3'])}  div={r['diverging']}")
print('HARD FAILURES:', bad)
