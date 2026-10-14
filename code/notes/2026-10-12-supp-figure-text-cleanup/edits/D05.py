TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/D05.py'
EDITS = [
 ("H1, H2, H3, H4 = 110.0, 100.0, 86.0, 150.0", "H1, H2, H3, H4 = 120.0, 105.0, 92.0, 150.0   # 2026-10-12: panel titles -> x-axis titles"),
 ("F.text((xa.p0+xa.p1)/2, y+3.0, 'P–S route mean ΔDWUI, 16 settings', 7.5, ha='center')\nF.text((xd.p0+xd.p1)/2, y+3.0, 'Supported → expanded change', 7.5, ha='center')\n", ""),
 ("['−.04', '−.03', '−.02', '−.01', '0', '.01', '.02'])", "['−.04', '−.03', '−.02', '−.01', '0', '.01', '.02'], title='P–S route mean ΔDWUI')"),
 ("F.xaxis(xd, ay, [-0.008, -0.004, 0, 0.004], ['−.008', '−.004', '0', '.004'])", "F.xaxis(xd, ay, [-0.008, -0.004, 0, 0.004], ['−.008', '−.004', '0', '.004'], title='Expanded − supported')"),
 ("F.text((xb.p0+xb.p1)/2, y+3.0, 'Gene P–S ΔDWUI', 7.5, ha='center')\nF.text((xs.p0+xs.p1)/2, y+3.0, 'Route mean by source', 7.5, ha='center')\n", ""),
 ("['−.15', '−.1', '−.05', '0', '.05', '.1'])", "['−.15', '−.1', '−.05', '0', '.05', '.1'], title='Gene P–S ΔDWUI')"),
 ("F.xaxis(xs, ay, [-0.04, -0.02, 0, 0.02], ['−.04', '−.02', '0', '.02'])", "F.xaxis(xs, ay, [-0.04, -0.02, 0, 0.02], ['−.04', '−.02', '0', '.02'], title='Route mean P–S ΔDWUI')"),
 ("F.text((xc.p0+xc.p1)/2, y+3.0, 'Mean P–S ΔDWUI', 7.5, ha='center')\n", ""),
 ("F.xaxis(xc, ay, [-0.06, -0.04, -0.02, 0, 0.02], ['−.06', '−.04', '−.02', '0', '.02'])", "F.xaxis(xc, ay, [-0.06, -0.04, -0.02, 0, 0.02], ['−.06', '−.04', '−.02', '0', '.02'], title='Mean P–S ΔDWUI')"),
 ("        F.text((sx.p0+sx.p1)/2, y+3.0, 'hN − P ΔDWUI, one RSeT library', 7.5, ha='center')\n", ""),
 ("        F.text((sx.p0+sx.p1)/2, y+3.0, 'Same library, primary setting', 7.5, ha='center')\n", ""),
 ("['−.06', '−.04', '−.02', '0', '.02', '.04'])", "['−.06', '−.04', '−.02', '0', '.02', '.04'], title='hN − P ΔDWUI' if L == 'd' else 'ΔDWUI')"),
]
