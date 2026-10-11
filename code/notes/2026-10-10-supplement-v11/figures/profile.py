"""Supplementary Figure S20 (v11): state-resolved APA profile classes and partial endpoint retention.
a–c former main Figure 2 a–c (Fig/Figure.2.pdf), vector crop, letters unchanged
d   former Supplementary Figure S7c (Fig/SuppFig.S7.pdf = consolidated G), vector crop, relettered c -> d
Former main Figure 2 d–f and former S6, S7 a,b,d,e, S8 and S11 are not carried into the v11 supplement (author decision
2026-10-10: one condensed profile figure supporting the Discussion statement on endpoint retention).
No data are redrawn; the receipt records source SHA-256 and regions."""
import json, subprocess, tempfile
from pathlib import Path
from common import *                                     # noqa
from common import supdraw, name, OUT, ROOT, crop, place, new_page, sha
import pikepdf

F2 = ROOT/'revision-backups/20261011-before-v11-apply/Fig/Figure.2.pdf'; H2 = 556.74
S7 = ROOT/'revision-backups/20261011-before-v11-apply/Fig/SuppFig.S7.pdf'; H7 = 384.4
R2 = (0.0, 0.0, PAGE_W, 322.0)                  # former Figure 2 rows a–c incl. their keys
R7 = (0.0, 125.8, PAGE_W, 273.6)                # former S7 c
DY7 = 322.0+GAP-126.3
H = R7[3]+DY7+4.5
dst = new_page(H)
place(dst, H, crop(F2, H2, R2), R2, 0.0, 0.0)
place(dst, H, crop(S7, H7, R7, extra_clear=[(4.0, 126.0, 13.5, 138.5)]), R7, 0.0, R7[1]+DY7)
F = supdraw.Fig(H, name('profile'))
F.t.text_top(4.5+0.2, 126.3+0.4+DY7, 'd', 11.0, font='Bol')
tmp = Path(tempfile.mkdtemp()); F.t.save(tmp/'ov.pdf')
ov = pikepdf.open(tmp/'ov.pdf')
dst.pages[0].add_overlay(ov.pages[0], pikepdf.Rectangle(0, 0, PAGE_W, H))
OUT.mkdir(exist_ok=True)
out = OUT/f'{F.name}.pdf'; dst.save(out)
try:
    from pdfclean import finalize
    finalize(out, label=F.name, log=OUT/'sizes.tsv')
except Exception as e:
    print('pdfclean', e)
subprocess.run(['pdftoppm', '-r', '150', '-png', '-singlefile', str(out), str(OUT/F.name)], check=True, capture_output=True)
rec = dict(figure=F.name, page_pt=[PAGE_W, H], pdf_sha256=sha(out), redrawn=False,
           kept=[dict(source='Fig/Figure.2.pdf', source_sha=sha(F2), region=R2, dy=0.0, letters='a, b, c unchanged'),
                 dict(source='Fig/SuppFig.S7.pdf', source_sha=sha(S7), region=R7, dy=DY7, letters={'c': 'd'})],
           not_carried=['former Figure 2 d–f', 'former S6', 'former S7 a, b, d, e', 'former S8', 'former S11'])
(OUT/f'{F.name}.receipt.json').write_text(json.dumps(rec, indent=1))
print(F.name, 'H', round(H, 1), 'mm', round(H/72*25.4, 1))
