"""Production main-figure revision 2026-10-12 (in-place PDF edits on the 10-11 production Fig/*.pdf).

Fig 2  b: subheader shortened to 'Same-state halves' (median recurrent/observed moved to caption);
          'Wells, mouse P–M' -> 'Well-level DWUI, mouse P–M (q < .05)'.
Fig 3  d: species titles above the PDLIM7 panel removed (species named in the caption).
Fig 4  c: funnel box 'strict-interval' -> 'strict recurrent' (same term as panel b and the caption).
Caption-only changes (Figure 3 b and d, Figure 4 c and e) are in main/results/*.tex, not here.
"""
import json
from fixlib import Job, chars_of, textwidth, INK, MUTED
recs = {}

# ---- Figure 2
J = Job(2)
sub = chars_of(J.src, 'Same-state', (249, 110, 290, 120))
J.remove((249.5, 110.6, 440.0, 119.4), '2b subheader shortened (caption carries the description)')
J.remove((249.5, 83.8, 400.0, 93.6), '2b wells label (replaced by Well-level wording)')
wl = chars_of(J.src, 'Wells,', (248, 84, 275, 93))
def draw2(t):
    t.text(sub['origin'], sub['base'], 'Same-state halves', 7.5, color=MUTED)
    t.text(wl['origin'], wl['base'], 'Well-level DWUI, mouse P–M (q < .05)', 7.5, color=MUTED)
recs['fig2'] = J.run(draw=draw2)

# ---- Figure 3
J = Job(3)
J.remove((19.5, 295.6, 81.5, 305.7), '3d species title (human) removed; species in caption')
J.remove((134.5, 295.6, 193.5, 305.7), '3d species title (mouse) removed; species in caption')
recs['fig3'] = J.run()

# ---- Figure 4
J = Job(4)
fw = chars_of(J.src, 'strict-interval', (358.0, 38.5, 411.0, 48.0))
def draw4(t):
    # the old word is a glyph run that clear() does not remove here; cover it with the box fill, then redraw
    t.rect(357.2, 38.6, 54.6, 9.6, fc='#F3F5F6')
    t.text(fw['origin'], fw['base'], 'strict recurrent', 8.0, font='Bol', color=INK)
recs['fig4'] = J.run(draw=draw4)

json.dump({k: {kk: r.get(kk) for kk in ('edits', 'removed_kinds', 'words_removed', 'words_added', 'sha256', 'out')} for k, r in recs.items()},
          open('out/receipt.rev-main.json', 'w'), indent=1, ensure_ascii=False)
print(json.dumps({k: (r['out'], r['sha256'][:12], r.get('words_added')) for k, r in recs.items()}, ensure_ascii=False, indent=1))
