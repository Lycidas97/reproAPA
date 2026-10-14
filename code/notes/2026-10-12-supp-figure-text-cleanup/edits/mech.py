TARGET = 'notes/2026-10-10-supplement-v11/figures/mech.py'
EDITS = [
 ("""S15 = ROOT/'notes/2026-10-06-supplement-polyapipe/consolidated/O/SuppFig.S15.pdf'
H15 = 534.9
OLD = {'a': (4.5, 4.5), 'b': (315.8, 4.5), 'c': (4.5, 113.8), 'd': (277.8, 113.8)}
H1, H2 = 88.6, 160.5
Y_OLD = 4.5+H1+GAP+H2+GAP                  # top of the kept former-S15 rows
OLD_REG = (0.0, 0.0, PAGE_W, 113.8+96.0+0.6)  # former S15 rows 1-2 incl. keys""",
  """import os, pdfplumber
# 2026-10-12: f–i are cropped from a composed page holding only former S15 a–d (consolidated O, bands D19 ab / D20 ab),
# rebuilt after the in-figure text cleanup (V11_S19_FI); letters are read from that page, the whole page is kept.
S15 = Path(os.environ['V11_S19_FI']) if os.environ.get('V11_S19_FI') else ROOT/'notes/2026-10-06-supplement-polyapipe/consolidated/O/SuppFig.S15.pdf'
with pdfplumber.open(S15) as _pdf:
    _pg = _pdf.pages[0]; H15 = float(_pg.height)
    _L = {w_['text']: (w_['x0'], w_['top']) for w_ in _pg.extract_words(extra_attrs=['size']) if w_['size'] > 10.5 and len(w_['text']) == 1}
OLD = {k_: (round(_L[k_][0], 1), round(_L[k_][1]-0.4, 1)) for k_ in 'abcd'}
H1, H2 = 88.6, 160.5
Y_OLD = 4.5+H1+GAP+H2+GAP                  # top of the kept former-S15 rows
OLD_REG = (0.0, 0.0, PAGE_W, H15-4.5)        # whole composed page above its bottom margin"""),
 ("F.text(x+16.0, y+0.6, 'S–M on the shared core', 7.5)\n", "# 2026-10-12: panel titles of a–e moved to axis titles or the caption\n"),
 ("F.text(x+16.0, y+0.6, 'Schiebinger2019–2i, within dates', 7.5)\n", ""),
 ("F.text((sb.p0+sb.p1)/2, y+11.0, 'r, cycle vs proximal strength', 7.0, ha='center', color=MUTED)\nF.text(XR-0.2, y+11.0, 'q', 7.0, ha='right', color=MUTED)",
  "F.text(XR-0.2, y+3.0, 'q', 7.0, ha='right', color=MUTED)"),
 ("yb_ = F.xaxis(sb, yb, [-1, 0, 1], ['−1', '0', '1'])", "yb_ = F.xaxis(sb, yb, [-1, 0, 1], ['−1', '0', '1'], title='r, cycle vs proximal strength')"),
 ("    ('Long-isoform clearance, 67 pairs',", "    ('Long-isoform clearance',"),
 ("    ('External targets, 605 genes',", "    ('External targets',"),
 ("    ('Membrane functions, 751 pairs',", "    ('Membrane functions',"),
 ("allg = GROUPS+[('RNA regulation, 48 shared orthologs', [(lab, v) for lab, v in rna])]", "allg = GROUPS+[('RNA regulation', [(lab, v) for lab, v in rna])]   # background sizes: caption"),
 ("F.text(x+16.0, y+0.6, 'Candidate mechanisms and functions', 7.5)\n", ""),
 ("F.text((scc.p0+scc.p1)/2, y+11.0, 'Observed / expected', 7.0, ha='center', color=MUTED)\nF.text(x+WC-18.0, y+11.0, 'n/exp.', 7.0, ha='right', color=MUTED); F.text(x+WC, y+11.0, 'q', 7.0, ha='right', color=MUTED)",
  "F.text(x+WC-18.0, y+3.0, 'n/exp.', 7.0, ha='right', color=MUTED); F.text(x+WC, y+3.0, 'q', 7.0, ha='right', color=MUTED)"),
 ("yc = F.xaxis(scc, yb, [0, 1, 2, 4, 6], ['0', '1', '2', '4', '6'])", "yc = F.xaxis(scc, yb, [0, 1, 2, 4, 6], ['0', '1', '2', '4', '6'], title='Observed / expected')"),
 ("HD = 74.0", "HD = 86.0"),
 ("F.text(x+16.0, y+0.6, 'U1 perturbation', 7.5)\n", ""),
 ("F.text((sd.p0+sd.p1)/2, y+11.0, 'Reciprocal response, fraction', 7.0, ha='center', color=MUTED)\n", ""),
 ("yd = F.xaxis(sd, yb, [0, 0.5, 1], ['0', '.5', '1'])", "yd = F.xaxis(sd, yb, [0, 0.5, 1], ['0', '.5', '1'], title='Reciprocal response (fraction)')"),
 ("F.text(x+16.0, y+0.6, 'PAS anchor architecture', 7.5)\n", ""),
 ("F.text((se.p0+se.p1)/2, y+11.0, 'Core − expected', 7.0, ha='center', color=MUTED)\nF.text(XR-0.2, y+11.0, 'q', 7.0, ha='right', color=MUTED)",
  "F.text(XR-0.2, y+3.0, 'q', 7.0, ha='right', color=MUTED)"),
 ("yee = F.xaxis(se, yb, [-0.5, 0, 0.5], ['−.5', '0', '.5'])", "yee = F.xaxis(se, yb, [-0.5, 0, 0.5], ['−.5', '0', '.5'], title='Core − expected')"),
]
