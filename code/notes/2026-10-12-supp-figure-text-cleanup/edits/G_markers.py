TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/G_markers.py'
EDITS = [
 ("TOT = 2*HDR+GAP+LEGH+nmax*PITCH_LONG+0.8\npitch = {'left': (TOT-2*HDR-GAP-LEGH-0.8)/nL, 'right': (TOT-2*HDR-GAP-LEGH-0.8)/nR}",
  "# 2026-10-12: J/F and bold/flag notes moved to the caption, so only the left column (size/colour key) reserves LEGH\n"
  "TOT = 2*HDR+GAP+0.8+max(nL*PITCH_LONG+LEGH, nR*PITCH_LONG)\npitch = {'left': (TOT-2*HDR-GAP-LEGH-0.8)/nL, 'right': (TOT-2*HDR-GAP-0.8)/nR}"),
 ("""def leg_flags(x, w, gx0, xj, yl):
    xc = x+w/2
    F.text(xc, yl, 'J, Jaccard; F, display-family compatibility (9-partition medians)', 7.5, ha='center')
    F.text(xc, yl+10.0, f'Bold: < {LOWJ:.2f}'.replace('0.', '.')+' · * candidate or mixed'
           + (' · † < 20 cells in every library' if SP == 'mouse' else ''), 7.5, ha='center')
""", ""),
 ("panel('d', PANELS[3][1], PANELS[3][2], PANELS[3][3], XR, yb+GAP, WL, pitch['right'], leg_flags)",
  "panel('d', PANELS[3][1], PANELS[3][2], PANELS[3][3], XR, yb+GAP, WL, pitch['right'], None)   # J/F/bold/flag notes: caption"),
]
