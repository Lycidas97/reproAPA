"""Assemble the v11 supplementary figure set S1–S21 into ../set/ (candidate only; Fig/ is never written).
Copies: production Fig/SuppFig.S<old>.pdf under the new number (no change).
Splits: new S7 = old S9 d–h (relettered a–e), new S8 = old S12 f–i (relettered a–d), vector crops moved to the top.
Text fix: new S15 (old S20) in-figure reference 'Figure 5f' -> 'Figure 4d' (v11 main numbering).
New: S10, S12, S18, S19, S20 from ../out/. Writes set/manifest.tsv with source and SHA-256."""
import json, shutil, subprocess, tempfile
from pathlib import Path
from common import *                                     # noqa
from common import supdraw, ROOT, OUT, crop, place, new_page, sha, clear
import pikepdf, pdfplumber

FIG = ROOT/'revision-backups/20261011-before-v11-apply/Fig'
SET = OUT.parent/'set'; SET.mkdir(exist_ok=True)
COPY = {1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 10, 9: 13, 11: 14, 13: 16, 14: 19, 16: 17, 17: 18, 21: 21}
NEW = {10: 'strat', 12: 'naive', 18: 'cycle', 19: 'mech', 20: 'profile'}
rows = []


def page_h(p):
    with pikepdf.open(p) as pdf:
        return float(pdf.pages[0].MediaBox[3])


def letters(p):
    with pdfplumber.open(p) as pdf:
        return {w['text']: (w['x0'], w['top']) for w in pdf.pages[0].extract_words(extra_attrs=['size']) if w['size'] > 10.5 and len(w['text']) == 1}


def overlay_save(dst, H, draw, name_):
    F = supdraw.Fig(H, name_)
    draw(F)
    tmp = Path(tempfile.mkdtemp()); F.t.save(tmp/'ov.pdf')
    ov = pikepdf.open(tmp/'ov.pdf')
    dst.pages[0].add_overlay(ov.pages[0], pikepdf.Rectangle(0, 0, PAGE_W, H))
    out = SET/f'{name_}.pdf'; dst.save(out)
    return out


def split(new, old, y_top, y_bot, relabel):
    src = FIG/f'SuppFig.S{old}.pdf'; Hs = page_h(src); L = letters(src)
    reg = (0.0, y_top, PAGE_W, y_bot)
    dy = 4.5-(y_top+0.5)
    H = y_bot+dy+4.5
    lb = [(L[o][0]-0.6, L[o][1]-0.6, L[o][0]+8.0, L[o][1]+11.6) for o in relabel]
    dst = new_page(H)
    place(dst, H, crop(src, Hs, reg, extra_clear=lb), reg, 0.0, y_top+dy)
    def draw(F):
        for o, n in relabel.items():
            F.t.text_top(L[o][0]-0.2+0.0, L[o][1]-0.4+dy, n, 11.0, font='Bol')
    out = overlay_save(dst, H, draw, f'SuppFig.S{new}')
    rows.append((new, f'Fig/SuppFig.S{old}.pdf', sha(src), f'split: region y {y_top}-{y_bot}, letters {relabel}', sha(out)))


for new, old in COPY.items():
    src = FIG/f'SuppFig.S{old}.pdf'; dst = SET/f'SuppFig.S{new}.pdf'
    shutil.copyfile(src, dst)
    rows.append((new, f'Fig/SuppFig.S{old}.pdf', sha(src), 'copy' if new != old else 'copy (same number)', sha(dst)))

split(7, 9, 102.7, 515.3, {'d': 'a', 'e': 'b', 'f': 'c', 'g': 'd', 'h': 'e'})
split(8, 12, 251.6, 507.9, {'f': 'a', 'g': 'b', 'h': 'c', 'i': 'd'})

# S15 = old S20 with the main-figure reference updated: the whole note line is removed from the content stream
# (one BT..ET block, font size 6.8, starting with 'Upper') and redrawn with the v11 panel reference.
src = FIG/'SuppFig.S20.pdf'; Hs = page_h(src)
with pdfplumber.open(src) as pdf:
    ws = pdf.pages[0].extract_words(extra_attrs=['size'])
    w0 = [x for x in ws if x['text'] == 'Upper' and abs(x['size']-6.8) < 0.05][0]
    line = ' '.join(x['text'] for x in ws if abs(x['top']-w0['top']) < 0.5 and abs(x['size']-6.8) < 0.05)
assert line.startswith('Upper three tiers = Figure 5f rug'), line
blk = pikepdf.open(src); pg = blk.pages[0]
ops = pikepdf.parse_content_stream(pg)
bts = [i for i, (a_, o) in enumerate(ops) if str(o) == 'BT']
hit = []
for i in bts:
    j = i
    while str(ops[j][1]) != 'ET': j += 1
    tf = [o for o in ops[i:j] if str(o[1]) == 'Tf']
    tj = [o for o in ops[i:j] if str(o[1]) == 'TJ']
    txt = ''.join(bytes(x).decode('latin-1') for x in tj[0][0][0] if isinstance(x, pikepdf.String)).replace('\x00', '') if tj else ''
    if tf and abs(float(tf[0][0][1])-6.8) < 0.01 and txt.startswith('Upper three'):
        hit.append((i, j))
assert len(hit) == 1, hit
i, j = hit[0]
ops = ops[:i]+ops[j+1:]
pg.Contents = blk.make_stream(pikepdf.unparse_content_stream(ops))
tmp = Path(tempfile.mkdtemp()); blk.save(tmp/'s20.pdf')
dst = pikepdf.open(tmp/'s20.pdf')
new_line = line.replace('Figure 5f', 'Figure 4d')
out = overlay_save(dst, Hs, lambda F: F.t.text_top(w0['x0'], w0['top'], new_line, 6.8, color=MUTED), 'SuppFig.S15')
rows.append((15, 'Fig/SuppFig.S20.pdf', sha(src), "in-figure note 'Figure 5f' -> 'Figure 4d' (line removed and redrawn)", sha(out)))

for new, code in NEW.items():
    s = OUT/f'SuppFig.S{new}.pdf'; d = SET/f'SuppFig.S{new}.pdf'
    shutil.copyfile(s, d); rows.append((new, f'notes/2026-10-10-supplement-v11/out/SuppFig.S{new}.pdf ({code})', sha(s), 'new', sha(d)))

rows.sort()
with open(SET/'manifest.tsv', 'w') as f:
    f.write('figure\tsource\tsource_sha256\taction\tsha256\tpage_mm\n')
    for new, s, ss, a, ds in rows:
        f.write(f'S{new}\t{s}\t{ss}\t{a}\t{ds}\t{page_h(SET/f"SuppFig.S{new}.pdf")/72*25.4:.1f}\n')
for n in range(1, 22):
    subprocess.run(['pdftoppm', '-r', '110', '-png', '-singlefile', str(SET/f'SuppFig.S{n}.pdf'), str(SET/f'SuppFig.S{n}')], check=True)
subprocess.run(['qpdf', '--empty', '--pages']+[str(SET/f'SuppFig.S{n}.pdf') for n in range(1, 22)]+['--', str(SET/'Supplement-v11-S1-S21.pdf')], check=True)
print(open(SET/'manifest.tsv').read())
