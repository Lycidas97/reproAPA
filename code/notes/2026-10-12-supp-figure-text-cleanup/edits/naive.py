TARGET = 'notes/2026-10-10-supplement-v11/figures/naive.py'
EDITS = [
 ("TOP = 4.5+24.0", "TOP = 4.5+19.0   # 2026-10-12: panel titles moved to axis titles or the caption"),
 ("H = YB+24.6", "H = YB+34.2"),
 ("F.text((sa.p0+sa.p1)/2, 4.5+0.6, 'Overlap with mouse set', 7.5, ha='center')\nF.text((sa.p0+sa.p1)/2, 4.5+10.6, 'Obs./expected, |ΔDWUI| ≥ .02', 7.0, ha='center', color=MUTED)\nF.text(xa+WA, 4.5+10.6, 'Pairs', 7.0, ha='right', color=MUTED)",
  "F.text(xa+WA, 4.5+3.0, 'Pairs', 7.0, ha='right', color=MUTED)"),
 ("ya = F.xaxis(sa, YB, [1, 1.5], ['1', '1.5'])", "ya = F.xaxis(sa, YB, [1, 1.5], ['1', '1.5'], title='Obs./expected overlap')"),
 ("    F.text(x0+16.0, 4.5+0.6, f'Naive − control', 7.5)\n    F.text(x0+16.0, 4.5+10.6, f'Fold difference, {thr:.2f}'.replace('0.', '.'), 7.0, color=MUTED)\n    F.text(x0+w, 4.5+10.6, 'q', 7.0, ha='right', color=MUTED)",
  "    F.text(x0+w, 4.5+3.0, 'q', 7.0, ha='right', color=MUTED)"),
 ("    return F.xaxis(s, YB, ticks, labs)", "    return F.xaxis(s, YB, ticks, labs, title=f'Fold difference ({thr:.2f})'.replace('0.', '.'))"),
]
