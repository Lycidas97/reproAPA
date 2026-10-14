TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/D09.py'
# 2026-10-12 (author): S6f '700' and '1,400' were 4.1 pt apart on a 39 pt axis; ticks at 0 and 1,000
EDITS = [("F.xaxis(fx, yy+30.0, [0, 700, 1400], ['0', '700', '1,400'], title='Recurrent genes')",
          "F.xaxis(fx, yy+30.0, [0, 1000], ['0', '1,000'], title='Recurrent genes')")]
