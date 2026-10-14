"""Figure 1b (2026-10-13): three arrowheads entered P's left side within 6.8 pt (E→P brown at y 357.3, E→P dark at
361.7, and the dashed M→P 'unresolved interval' at 364.1, 1.7 pt above P's bottom-left corner). The dashed M→P
connector is re-routed: same start on M (125.13, 379.17), right to x 184.5, then up into P's bottom edge, where it
ends with the same arrowhead as the dash-dot M→P connector 8.2 pt to its right. Style (white 1.26 pt halo under a
0.55 pt ink line, dash 2.27/1.26 pt) is unchanged; the new route crosses no other connector."""
import json
from fixlib import Job
J = Job(1)
H = J.H
# old dashed path (halo + ink share one bbox) and its right-pointing arrowhead (halo + ink)
J.exact += [((125.13, 364.061, 177.563, 379.168), '1b dashed M→P path'),
            ((176.187, 363.613, 177.503, 364.509), '1b dashed M→P arrowhead')]
INK = '#243746'   # (0.1411, 0.2156, 0.2744) as the original connector
X_IN = 184.5
ROUTE = [(125.13, 379.168), (X_IN, 379.168), (X_IN, 366.037)]
HEAD = [(X_IN-0.658, 366.546), (X_IN, 366.059), (X_IN+0.658, 366.546)]   # as the dash-dot arrowhead at x 192.657


def draw(t):
    xs, ys = [p[0] for p in ROUTE], [p[1] for p in ROUTE]
    for lw, col in ((1.258543, 'white'), (0.550613, INK)):
        t.line(xs, ys, lw=lw, color=col, dashes=(2.271571/lw, 1.261983/lw), solid_joinstyle='miter', dash_joinstyle='miter')
    for lw, col in ((1.258543, 'white'), (0.550614, INK)):
        t.ax.fill([p[0] for p in HEAD], [H-p[1] for p in HEAD], fc=INK, ec=col, lw=lw, joinstyle='miter', zorder=4 if col == INK else 3)


rec = J.run(draw=draw)
json.dump(rec, open('out/receipt.fig1.json', 'w'), indent=1, ensure_ascii=False)
print('exact', rec['exact_removed'], 'removed', rec['words_removed'], 'added', rec['words_added'])
