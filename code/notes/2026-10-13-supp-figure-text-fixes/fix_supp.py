"""Supplementary figures: in-place text fixes after reviews/2026-10-13-figure-recheck/REPORT.zh.md (items 1–13).
PDF-level edits on the production Fig/SuppFig.S*.pdf; plotted values are untouched. Each figure writes
out/SuppFig.S*.pdf and out/receipt.S*.json. The builder scripts that draw these labels are listed in README.zh.md."""
import json, sys
from fixlib import Job, chars_of, textwidth, ROOT, INK
SUB_SIZE, SUB_DROP = 6.8, 1.45          # subscript as in Figure 3e (q8): 6.8 pt, baseline 1.45 pt lower
recs = {}


def job(n):
    return Job(n, src=ROOT/f'Fig/SuppFig.S{n}.pdf', name=f'SuppFig.S{n}')


def finish(n, J, **kw):
    rec = J.run(**kw); recs[f'S{n}'] = rec
    json.dump(rec, open(f'out/receipt.S{n}.json', 'w'), indent=1, ensure_ascii=False)
    print(f"S{n}: removed {rec['words_removed']} added {rec['words_added']} sha {rec['sha256'][:12]}")


def subscript_q(J, box, sub, tail, tail_box):
    """The legend run 'q<sub> <tail>' is one text object: remove it and redraw q, a lowered smaller subscript, and the
    tail at its original position (font, size and colour of the original)."""
    c = J.drop_word('q'+sub, box, note=f'q{sub} -> q with subscript {sub}; run tail {tail!r} redrawn in place')
    t0 = chars_of(J.src, tail.split()[0], tail_box)
    q_w = textwidth('q', c['size'])
    return lambda t: (t.text(c['origin'], c['base'], 'q', c['size'], color=c['color']),
                      t.text(c['origin']+q_w, c['base']+SUB_DROP, sub, SUB_SIZE, color=c['color']),
                      t.text(t0['origin'], t0['base'], tail, t0['size'], color=t0['color']))


only = set(sys.argv[1:])
def want(n): return not only or f'S{n}' in only

# ---- S1e, S2e (item 9): the three program tracks share one scale; keep the top tick label only
for n, lab, tops in [(1, '3.5', (455.3, 477.6)), (2, '2.5', (392.5, 414.8))]:
    if not want(n): continue
    J = job(n)
    for top in tops:
        J.drop_word(lab, (38, top+1, 52, top+6), note=f'e: repeated top tick label {lab} beside a track name')
    finish(n, J)

# ---- S3f (item 5): 'original' -> 'primary fit' (the caption is changed to match); legend lines re-centred
if want(3):
    J = job(3)
    o1 = chars_of(J.src, 'original', (329, 443, 357, 449)); rv = chars_of(J.src, 'refit', (367, 443, 383, 449))
    o2 = chars_of(J.src, 'original', (327, 453, 354, 459))
    J.remove((330.7, 442.2, 410.0, 450.2), 'f legend line 1 text (redrawn)')
    J.remove((328.3, 452.2, 417.2, 460.2), 'f legend line 2 text (redrawn)')
    size, col = o1['size'], o1['color']
    t1, t2, t3 = 'primary fit', 'refit variants', 'primary-fit null 95%, centered'
    # line 1: [diamond 5.1] 2.2 t1 7.0 [dot 3.6] 3.0 t2 ; centre kept at the old centre
    w1 = 5.1+2.2+textwidth(t1, size)+7.0+3.6+3.0+textwidth(t2, size); c1 = (323.6+409.5)/2; x1 = c1-w1/2
    dx_d = x1-323.6; x_t1 = x1+5.1+2.2; x_dot = x_t1+textwidth(t1, size)+7.0; dx_dot = x_dot-362.1; x_t2 = x_dot+3.6+3.0
    # line 2: [bar 10.0] 2.5 t3 ; centre kept
    w2 = 10.0+2.5+textwidth(t3, size); c2 = (316.0+416.9)/2; x2 = c2-w2/2; dx_b = x2-316.0
    def draw3(t):
        t.text(x_t1, o1['base'], t1, size, color=col); t.text(x_t2, rv['base'], t2, size, color=col)
        t.text(x2+12.5, o2['base'], t3, size, color=col)
    finish(3, J, draw=draw3, shifts=[((322.0, 442.5, 330.2, 450.0), dx_d, 0.0, 'f diamond key'),
                                     ((360.5, 443.0, 367.2, 449.5), dx_dot, 0.0, 'f refit dot key'),
                                     ((315.0, 454.5, 327.5, 458.0), dx_b, 0.0, 'f null-range bar key')])
    recs['S3']['layout'] = dict(dx_diamond=dx_d, dx_dot=dx_dot, dx_bar=dx_b)

# ---- S8b (item 13): legend 'Tested class (outlined)' -> 'Tested class', re-centred
if want(8):
    J = job(8)
    tc = chars_of(J.src, 'Tested', (352, 128.9, 378, 137.0))
    J.remove((353.5, 128.8, 430.0, 137.2), 'b legend run Tested class (outlined): removed, redrawn without the parenthesis')
    dx8 = (429.3-396.0)/2                  # keep the key centred where it was
    finish(8, J, draw=lambda t: t.text(tc['origin']+dx8, tc['base'], 'Tested class', tc['size'], color=tc['color']),
           shifts=[((345.8, 129.6, 352.4, 136.2), dx8, 0.0, 'b legend key square re-centred')])

# ---- S9a (item 4) subscripts; S9d (item 1) corner tick label
if want(9):
    J = job(9)
    d1 = subscript_q(J, (52, 113.5, 66.5, 121.5), '16', '< .05', (67, 113.5, 72.5, 121.5))
    d2 = subscript_q(J, (97.5, 113.5, 111.5, 121.5), '16', '≥ .05', (112, 113.5, 117.5, 121.5))
    J.drop_word('−.6', (20, 258.0, 28.5, 265.6), note='d y tick label at the corner (touched the x tick label −.6)')
    finish(9, J, draw=lambda t: (d1(t), d2(t)))

# ---- S11d (items 4, 6): subscript; sub-labels clear of the labels above
if want(11):
    J = job(11)
    d = subscript_q(J, (211, 248, 220.5, 255.8), '8', '< .05', (221, 248, 226, 255.8))
    finish(11, J, draw=d, shifts=[((4.0, 175.45, 64.5, 183.2), 1.2, 'd sub-label Mean of Xing pairs (1)'),
                                  ((4.0, 197.45, 64.5, 205.2), 1.2, 'd sub-label Mean of Xing pairs (2)'),
                                  ((4.0, 219.45, 24.5, 226.6), 1.2, 'd sub-label Liu2020 × Guo2019–CIP keys')])

# ---- S13f, S14b (item 1) corner tick labels; S14c (item 10) explanatory legend entry
if want(13):
    J = job(13)
    J.drop_word('0', (328, 401.8, 332.4, 410.3), note='f y tick label at the corner (touched the x tick label −4)')
    finish(13, J)
if want(14):
    J = job(14)
    J.drop_word('−600', (245, 133.4, 265.0, 141.0), note='b y tick label at the corner (touched the x tick label −600)')
    J.remove((255.0, 322.2, 261.0, 328.2), 'c legend NA key (empty square)')
    J.remove((262.6, 321.2, 382.0, 329.1), 'c legend NA explanation (in the caption; in-plot NA label stays)')
    finish(14, J, shifts=[((60.0, 320.2, 250.0, 330.4), (381.8-248.2)/2, 0.0, 'c legend line 2 re-centred under line 1')])

# ---- S15f (item 12): complete axis titles
if want(15):
    J = job(15)
    a = chars_of(J.src, 'from', (28, 330, 46, 337)); b = chars_of(J.src, 'from', (164, 330, 182, 337))
    J.remove((28.9, 329.5, 112.0, 337.9), 'f left axis title (redrawn complete)')
    J.remove((164.5, 329.5, 265.4, 337.9), 'f right axis title (redrawn complete)')
    L, R = 'Distance from proximal PAS (kb)', 'Distance from the terminal-exon start (kb)'
    cl, cr = (29.1+111.8)/2, (164.7+265.2)/2
    wl, wr = textwidth(L, a['size']), textwidth(R, b['size'])
    gap = (cr-wr/2)-(cl+wl/2)
    shift_r = max(0.0, 8.0-gap)           # keep >= 8 pt between the two titles
    print('S15 titles: left', round(cl-wl/2, 1), round(cl+wl/2, 1), 'right', round(cr+shift_r-wr/2, 1), round(cr+shift_r+wr/2, 1), 'shift', round(shift_r, 1))
    finish(15, J, draw=lambda t: (t.text(cl, a['base'], L, a['size'], color=a['color'], ha='center'),
                                  t.text(cr+shift_r, b['base'], R, b['size'], color=b['color'], ha='center')))
    recs['S15']['title_extents'] = dict(left=[cl-wl/2, cl+wl/2], right=[cr+shift_r-wr/2, cr+shift_r+wr/2])

# ---- S17a (item 8): legend row clear of the x-axis title
if want(17):
    J = job(17)
    finish(17, J, shifts=[((0.0, 134.95, 150.0, 143.8), 1.5, 'a legend row')])

# ---- S19e (item 7): not edited. The 'q' head does not touch the first row (an earlier probe crop had clipped it);
# the six row labels sit at a 7.0 pt pitch with no vertical room above (letter row) or below (one panel gap to f/g),
# so wider rows need a redraw from the builder (see README).

# ---- S20d (item 11): axis title without the definition (in the caption)
if want(20):
    J = job(20)
    c = chars_of(J.src, 'Endpoint', (300, 462, 335, 470))
    J.remove((301.7, 461.9, 409.2, 470.3), 'd axis title (redrawn without = 1 − T)')
    finish(20, J, draw=lambda t: t.text((301.9+409.0)/2, c['base'], 'Endpoint retention D/L', c['size'], color=c['color'], ha='center'))
