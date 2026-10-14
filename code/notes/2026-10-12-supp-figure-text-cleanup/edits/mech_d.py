TARGET = 'notes/2026-10-10-supplement-v11/figures/mech.py'
EDITS = [
 ("    ('Long-isoform clearance',", "    ('Long-isoform clearance, 67 pairs',"),
 ("    ('External targets',", "    ('External targets, 605 genes',"),
 ("    ('Membrane functions',", "    ('Membrane functions, 751 pairs',"),
 ("allg = GROUPS+[('RNA regulation', [(lab, v) for lab, v in rna])]   # background sizes: caption", "allg = GROUPS+[('RNA regulation, 48 shared orthologs', [(lab, v) for lab, v in rna])]"),
 ("    F.text(x+0.2, yy-3.0, gl, 7.0, color=MUTED); yy += pitch", "    F.text(x+0.2, yy-3.0, gl.split(',')[0], 7.0, color=MUTED); yy += pitch   # 2026-10-12: background sizes in the caption; values keys keep the full label"),
]
