TARGET = 'notes/2026-10-10-supplement-v11/figures/strat.py'
EDITS = [
 ("    F.text(cx0+2*cw-0.5, yk+10.5, '(n) shared proximal, |ΔDWUI| ≥ .02', 7.0, ha='center', color=MUTED)\n    return yk+10.5+0.93*7.0",
  "    # 2026-10-12: the '(n) shared proximal' note moved to the caption\n    return yk+0.93*7.5"),
 ("    F.text(x+20.0, y+0.6, 'Shared proximal P–M pairs by event class', 7.5)\n", ""),
 ("            F.text(px+bw+gin/2, y+0.6, f'{\"Both tandem-3′ UTR\" if cls == \"both_tandem\" else \"Both ALE-compatible\"}, {n} pairs', 7.5, ha='center')",
  "            F.text(px+bw+gin/2, y+0.6, 'Both tandem-3′ UTR' if cls == 'both_tandem' else 'Both ALE-compatible', 7.5, ha='center')   # pair count n: caption"),
 ("    F.text(x+16.0+supdraw.textwidth('DYNLL2/Dynll2', 8.0, 'BolIta')+4.0, y+0.9, 'source to P', 7.5, color=MUTED)\n", ""),
 ("    F.text((si.p0+si.p1)/2, y+11.0, 'Segment inclusion', 7.5, ha='center')\n", ""),
 ("    F.text(x+16.0, y+4.0, 'Corresponding switches', 7.5)\n", ""),
]
