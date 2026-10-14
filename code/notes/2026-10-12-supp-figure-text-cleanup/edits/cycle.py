TARGET = 'notes/2026-10-10-supplement-v11/figures/cycle.py'
EDITS = [
 ("F.text(x+16.0, y+0.6, 'Cell-cycle score difference', 7.5)\nF.text((s1.p0+s1.p1)/2, y+11.0, 'P–M, same libraries', 7.0, ha='center', color=MUTED)\nF.text((s2.p0+s2.p1)/2, y+11.0, 'P–S, documented source', 7.0, ha='center', color=MUTED)\ntop, pitch = y+24.0, 9.6",
  "# 2026-10-12: panel titles moved to axis titles or the caption; facet labels shortened (same libraries / documented source: caption)\nF.text((s1.p0+s1.p1)/2, y+3.0, 'P–M', 7.5, ha='center')\nF.text((s2.p0+s2.p1)/2, y+3.0, 'P–S', 7.5, ha='center')\ntop, pitch = y+20.4, 9.6"),
 ("    ya = F.xaxis(s, yb, [0, 0.2, 0.4, 0.6], ['0', '.2', '.4', '.6'])\nF.key_row(0, ya+2.5,",
  "    ya = F.xaxis(s, yb, [0, 0.2, 0.4, 0.6], ['0', '.2', '.4', '.6'])\nF.text((s1.p0+s2.p1)/2, ya+1.6, 'Cell-cycle score difference', 8.0, ha='center'); ya += 1.6+0.93*8.0\nF.key_row(0, ya+2.5,"),
 ("F.text(x+16.0, y+0.6, 'P–S, documented source', 7.5)\n", ""),
 ("F.text(x+16.0, y+0.6, 'CPA-factor expression, Δ mean log2(CPM + 1)', 7.5)\n", ""),
 ("F.text(CX(0)+2.5*cw, y+11.4, 'P–M, same libraries', 7.0, ha='center', color=MUTED)\nF.text(CX(5)+2.5*cw, y+11.4, 'P–S, documented source', 7.0, ha='center', color=MUTED)\nht = y+21.0",
  "F.text(CX(0)+2.5*cw, y+3.0, 'P–M', 7.5, ha='center')\nF.text(CX(5)+2.5*cw, y+3.0, 'P–S', 7.5, ha='center')\nht = y+12.6"),
 ("F.cbar_key(cx0+(CX(9)+cw-cx0)/2, yk, 'Δ', '−1.5', '+1.5',", "F.cbar_key(cx0+(CX(9)+cw-cx0)/2, yk, 'Δ mean log2(CPM + 1)', '−1.5', '+1.5',"),
 ("F.text(x+16.0, y+0.6, 'Upstream UGUA, distal − proximal', 7.5)\n", ""),
 ("F.text(XR-0.2, y+11.0, 'q', 7.0, ha='right', color=MUTED)\ntop, pitch = y+28.0, 11.2", "F.text(XR-0.2, y+3.0, 'q', 7.0, ha='right', color=MUTED)\ntop, pitch = y+28.0, 11.2"),
 ("    F.text(x+0.2, yh, bl+f', {int(p[p.background == bg].n_genes.iloc[0]) if bg == \"paired_strict\" else \"\"}'.rstrip(', '), 7.0, color=MUTED)",
  "    F.text(x+0.2, yh, bl, 7.0, color=MUTED)   # common-ortholog count: caption"),
 ("F.text(x+16.0, y+0.6, 'Within-state partial correlation of cell-cycle score with mean DWUI', 7.5)\n", ""),
 ("top, pitch = y+24.0, 8.4", "top, pitch = y+20.0, 8.4"),
 ("    F.text((s.p0+s.p1)/2, y+12.0, st, 7.5, ha='center', color=SC[st])", "    F.text((s.p0+s.p1)/2, y+3.0, st, 7.5, ha='center', color=SC[st])"),
]
# follow-up (sed): H1, H3 107.3, 108.1 -> 113.3, 105.1
# follow-up (sed): a top 20.4 -> 20.8 (letter zone), H1 113.7
