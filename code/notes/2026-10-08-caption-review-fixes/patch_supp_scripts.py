"""Label fixes in the supplementary figure scripts (2026-10-08 caption review). Each old string must match exactly once."""
from pathlib import Path
S = Path(__file__).resolve().parent.parent/'2026-10-06-supplement-polyapipe'
EDITS = [
 ('tools/compose.py', "lab = name if textwidth(name, 7.0)+7.6 < span+4.0 else", "lab = name if textwidth(name, 7.0)+7.6 < span+12.0 else"),   # S2e: Guo2019–OSK in full (A9)
 ('figures/D12.py', "'Responders in both studies', 7.5", "'Comparable genes', 7.5"),                                                  # S12g header (A7)
 ('figures/D16.py', "[('marker', ('o', SC[s]), f'{s}–S') for s in 'REPM']", "[('marker', ('o', SC[s]), f'{s}–S') for s in 'REPM' if (gp.state == s).any()]"),  # S13d: no M–S key (C3)
 ('figures/D22.py', "f'{cat}, partial (%)'", "f'{cat}, partial retention (%)'"),                                                      # S17c (A6)
 ('figures/D22.py', "f'{fac} sites, P–M,", "f'{fac.upper()} sites, P–M,"),                                                              # S17e CARHSP1 (B8)
 ('figures/G3.py', "[('Source / ECM', 'Source / ECM'), ('Epithelial', 'Epithelial'), ('Pluripotency', 'Pluripotency')]",
                   "[('Source / ECM', 'ECM'), ('Epithelial', 'Epithelial'), ('Pluripotency', 'Pluripotency')]"),                      # S2h: fixed ECM module as in S2e (B14)
]
for f, old, new in EDITS:
    p = S/f; t = p.read_text()
    assert t.count(old) == 1, (f, old, t.count(old))
    p.write_text(t.replace(old, new)); print('ok', f)
