TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/D22.py'
EDITS = [
 ("H1, H2 = 101.0, 106.0", "H1, H2 = 108.0, 106.0   # 2026-10-12: panel titles moved to axis titles or the caption"),
 ("    F.text((px.p0+px.p1)/2, y+3.0, f'{cat}, partial retention (%)', 7.5, ha='center')", "    F.text((px.p0+px.p1)/2, y+3.0, cat, 7.5, ha='center')"),
 ("    F.xaxis(px, ry_[ORD_A[-1]]+5.0, [0, 0.3, 0.6], ['0', '30', '60'])", "    ab_ = F.xaxis(px, ry_[ORD_A[-1]]+5.0, [0, 0.3, 0.6], ['0', '30', '60'])\nF.text(x+LBA+sw_+6.0, ab_+1.6, 'Partial retention (%)', 8.0, ha='center')"),
 ("F.text((bx.p0+bx.p1)/2, y+3.0, 'Partial in both studies', 7.5, ha='center')\n", ""),
 ("    F.text(x+16.0+textwidth(sym, 8.0, 'BolIta')+3.0, y+1.0, f'{fac.upper()} sites, P–M, {chrom} ({strand.replace(\"-\", \"−\")})', 7.0, color=MUTED)\n", "    # factor, contrast, chromosome and strand: caption\n"),
]
