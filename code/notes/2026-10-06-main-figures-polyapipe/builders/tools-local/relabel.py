"""Final stage for the main figures: panel letters in lowercase (house style from 2026-10-06, same as the supplement).

The builders and retile.py keep the uppercase master letters (retile.letters() locates panels by them). This stage
takes the retiled final PDF, removes each panel letter glyph and draws the lowercase letter, 11 pt Nimbus Sans L Bold,
on the same baseline and origin x, in the same colour. Nothing else is touched (the removal is checked: exactly one
glyph per letter, no other object inside the letter box). Usage: relabel.py IN.pdf OUT.pdf FIG
"""
import sys, json, tempfile
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from v5draw import Top, clear, INK
from pdfedit import overlay
import pikepdf, pdfplumber

EXPECT = {1: 'ABCD', 2: 'ABCDEF', 3: 'ABCDE', 4: 'ABCDE', 5: 'ABCDEF'}
SIZE = 11.0


def panel_letters(pdf):
    pg = pdfplumber.open(pdf).pages[0]
    return pg, [c for c in pg.chars if c['size'] >= 9.5 and 'Bol' in c['fontname'] and 'Ita' not in c['fontname']
                and 'Obl' not in c['fontname'] and len(c['text']) == 1 and c['text'].isalpha()]


def relabel(src, dst, fig):
    pg, L = panel_letters(src); W, H = float(pg.width), float(pg.height)
    got = ''.join(sorted(c['text'] for c in L))
    assert got == EXPECT[fig], (fig, got)
    nchars0 = len(pg.chars); rects = []; rec = {}
    t = Top(W, H)
    for c in L:
        r = (c['x0']-0.5, c['top']-0.5, c['x1']+0.5, c['bottom']+0.5); rects.append(r)
        ox, base_pdf = float(c['matrix'][4]), float(c['matrix'][5])
        col = c.get('non_stroking_color')
        colr = INK if not col else (tuple(col[:3]) if len(col) >= 3 else (col[0],)*3)
        t.text(ox, H-base_pdf, c['text'].lower(), SIZE, font='Bol', color=colr)
        rec[c['text'].lower()] = dict(old=c['text'], old_size=round(c['size'], 2), origin_x=round(ox, 2), baseline_top=round(H-base_pdf, 2))
    pdf = pikepdf.open(src)
    removed = clear(pdf, H, rects)
    tmp = Path(tempfile.mkdtemp())
    pdf.save(tmp/'cleared.pdf')
    pg1 = pdfplumber.open(tmp/'cleared.pdf').pages[0]
    assert len(pg1.chars) == nchars0-len(L), ('removal touched other glyphs', nchars0, len(pg1.chars), len(L))
    t.save(tmp/'letters.pdf')
    pdf = pikepdf.open(tmp/'cleared.pdf'); overlay(pdf, tmp/'letters.pdf', 0, 0)
    pdf.save(dst)
    pg2, L2 = panel_letters(dst)
    assert ''.join(sorted(c['text'] for c in L2)) == EXPECT[fig].lower(), [c['text'] for c in L2]
    return dict(figure=fig, source=str(src), removed_objects=removed if isinstance(removed, int) else len(removed), letters=rec)


if __name__ == '__main__':
    print(json.dumps(relabel(sys.argv[1], sys.argv[2], int(sys.argv[3])), indent=1, default=str))
