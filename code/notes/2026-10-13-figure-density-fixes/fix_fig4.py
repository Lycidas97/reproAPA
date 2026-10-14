"""Figure 4 density fixes (2026-10-13).
4b  122 dots (and the two key dots) in a beeswarm (15 rows at 1.8 pt) of 3.6 pt dots with a same-colour 0.5 pt edge: neighbours merged into
    one mass. Same centres, same footprint (3.6 pt), but a 0.45 pt white edge, so touching dots stay separable.
4c  75% of the 51 genes have MPRO 0.10–0.21 and sat in the bottom quarter of a linear axis. The y axis becomes
    logarithmic over 0.1–0.5 (ticks .1 .2 .3 .5); each point's MPRO is recovered exactly from the linear mapping
    (ticks .1 at y 161.16 and .5 at y 62.16) and re-placed; x positions, marker styles, the shading and the axis frame
    are unchanged; the Trip12 ring and leader line follow its point."""
import json, math, pdfplumber
from fixlib import Job, chars_of, ROOT
J = Job(4)
pg = pdfplumber.open(J.src).pages[0]
def is_marker(o, lo, hi):
    w, h = o['x1']-o['x0'], o['bottom']-o['top']
    return lo < w < hi and abs(w-h) < 0.6
hexc = lambda c: '#%02X%02X%02X' % tuple(round(v*255) for v in c[:3])

# ---- 4b
B = [o for o in pg.curves if 15 < o['x0'] < 335 and 124 < o['top'] < 168 and is_marker(o, 2.5, 6)]
assert len(B) == 122, len(B)
for o in B: J.exact.append(((o['x0'], o['top'], o['x1'], o['bottom']), '4b dot'))
dots_b = [((o['x0']+o['x1'])/2, (o['top']+o['bottom'])/2, hexc(o['non_stroking_color'])) for o in B]
KEY = [o for o in pg.curves if 250 < o['x0'] < 300 and 105 < o['top'] < 118 and is_marker(o, 2.5, 6)]   # Lost / Gained key
assert len(KEY) == 2, len(KEY)
for o in KEY: J.exact.append(((o['x0'], o['top'], o['x1'], o['bottom']), '4b key dot'))
dots_b += [((o['x0']+o['x1'])/2, (o['top']+o['bottom'])/2, hexc(o['non_stroking_color'])) for o in KEY]

# ---- 4c
C = [o for o in pg.curves if 345 < o['x0'] < 453.5 and 55 < o['top'] < 165 and is_marker(o, 2.0, 4.2)]
assert len(C) == 51, len(C)
Y01, Y05 = 161.16, 62.16
mpro = lambda y: 0.1+(Y01-y)/(Y01-Y05)*0.4
ylog = lambda m: Y01-math.log(m/0.1)/math.log(5)*(Y01-Y05)
pts_c = []
for o in C:
    J.exact.append(((o['x0'], o['top'], o['x1'], o['bottom']), '4c marker'))
    cy = (o['top']+o['bottom'])/2; m = mpro(cy)
    pts_c.append(dict(x=(o['x0']+o['x1'])/2, y_old=cy, mpro=m, y_new=ylog(m), d=o['x1']-o['x0'],
                      filled=o['non_stroking_color'] not in (1.0, (1.0,), [1.0], (1, 1, 1)), lw=o.get('linewidth'),
                      fc=o['non_stroking_color'], ec=o['stroking_color']))
ring = [o for o in pg.curves if 390 < o['x0'] < 391 and 120 < o['top'] < 121.5]; assert len(ring) == 1
lead = [o for o in pg.lines if abs(o['x0']-366.19) < 0.05 and abs(o['top']-78.18) < 0.05]; assert len(lead) == 1
for o in ring + lead: J.exact.append(((o['x0'], o['top'], o['x1'], o['bottom']), '4c Trip12 call-out'))
mid = [o for o in pg.lines if 353 < o['x0'] < 354 and abs(o['top']-111.66) < 0.05]; assert len(mid) == 1
J.exact.append(((mid[0]['x0'], mid[0]['top'], mid[0]['x1'], mid[0]['bottom']), '4c linear .3 tick'))
lab1 = chars_of(J.src, '.1', (345, 156, 353, 165))
J.drop_word('.3', (345, 106.5, 353, 115.6), note='4c linear .3 tick label')
trip = min(pts_c, key=lambda p: math.dist((p['x'], p['y_old']), (393.91, 124.13)))
ty = trip['y_new']; r_ring = (ring[0]['x1']-ring[0]['x0'])/2


def draw(t):
    for col in ('#0072B2', '#D55E00'):                       # blue first, orange (net gain) on top, as before
        for x, y, c in dots_b:
            if c == col: t.marker(x, y, 'o', 3.6, fc=c, ec='white', mew=0.45)
    for p in sorted(pts_c, key=lambda p: p['filled']):       # open markers first, filled on top
        if p['filled']: t.marker(p['x'], p['y_new'], 'o', p['d']-p['lw'], fc=hexc(p['fc']), ec=hexc(p['ec']), mew=p['lw'])
        else: t.marker(p['x'], p['y_new'], 'o', p['d']-p['lw'], fc='white', ec=hexc(p['ec']), mew=p['lw'])
    ink = hexc(ring[0]['stroking_color'])
    t.line([366.19, trip['x']], [78.18, ty-r_ring-0.2], lw=0.5, color=ink)
    t.ax.add_patch(__import__('matplotlib.patches', fromlist=['Circle']).Circle((trip['x'], t.H-ty), r_ring, fill=False, ec=ink, lw=1.0))
    for v, s in ((0.2, '.2'), (0.3, '.3')):
        y = ylog(v); t.line([353.38, 355.38], [y, y], lw=0.5, color='#87969E')
        t.text(lab1['x1'], lab1['base']+(y-Y01), s, lab1['size'], color=lab1['color'], ha='right')


rec = J.run(draw=draw)
rec['4c_points'] = [dict(x=round(p['x'], 3), mpro=round(p['mpro'], 5), y_old=round(p['y_old'], 3), y_new=round(p['y_new'], 3), filled=p['filled']) for p in pts_c]
rec['4c_ticks'] = {s: round(ylog(v), 3) for v, s in ((0.1, '.1'), (0.2, '.2'), (0.3, '.3'), (0.5, '.5'))}
json.dump(rec, open('out/receipt.fig4.json', 'w'), indent=1, ensure_ascii=False)
print('exact', rec['exact_removed'], 'removed', rec['words_removed'], 'added', rec['words_added'])
print('MPRO range', round(min(p['mpro'] for p in pts_c), 4), round(max(p['mpro'] for p in pts_c), 4), 'Trip12', round(trip['mpro'], 4), round(ty, 2))
