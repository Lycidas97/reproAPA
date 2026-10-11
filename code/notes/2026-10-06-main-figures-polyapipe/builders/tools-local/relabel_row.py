"""Replace one row label of a final main figure in place (2026-10-07: Figure 3b 'Within-species axis' -> 'Common PAS axis').

Same method as relabel.py: remove exactly the glyphs of the old label, draw the new label in the same font, size,
colour and baseline, right-aligned at the old right edge (the label column is right-aligned). Checks: only the old
label's glyphs are removed, the new label reads back, and the 150-dpi rendering outside the label box is unchanged.
Usage: relabel_row.py IN.pdf OUT.pdf 'old label' 'new label'
"""
import sys, json, tempfile, subprocess
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from v5draw import Top, clear, INK
from pdfedit import overlay
import pikepdf, pdfplumber
from PIL import Image, ImageChops


def find(pdf, label):
    pg = pdfplumber.open(pdf).pages[0]
    words = pg.extract_words(extra_attrs=['fontname', 'size'], keep_blank_chars=False)
    toks = label.split()
    hits = [i for i in range(len(words) - len(toks) + 1)
            if [w['text'] for w in words[i:i+len(toks)]] == toks and len({round(w['top'], 1) for w in words[i:i+len(toks)]}) == 1]
    assert len(hits) == 1, (label, hits)
    ws = words[hits[0]:hits[0]+len(toks)]
    box = (min(w['x0'] for w in ws), min(w['top'] for w in ws), max(w['x1'] for w in ws), max(w['bottom'] for w in ws))
    chars = [c for c in pg.chars if c['x0'] >= box[0]-0.01 and c['x1'] <= box[2]+0.01 and c['top'] >= box[1]-0.01 and c['bottom'] <= box[3]+0.01]
    return pg, box, chars


def render(pdf, png):
    subprocess.run(['pdftoppm', '-r', '150', '-png', '-singlefile', str(pdf), str(png)[:-4]], check=True)
    return Image.open(png).convert('RGB')


def main(src, dst, old, new):
    pg, box, chars = find(src, old); W, H = float(pg.width), float(pg.height)
    assert ''.join(c['text'] for c in chars).replace(' ', '') == old.replace(' ', ''), [c['text'] for c in chars]
    chars_v = [c for c in chars if c['text'].strip()]; c0 = chars_v[0]; size = round(float(c0['size']), 2); base = H - float(c0['matrix'][5])
    assert all(abs(float(c['matrix'][5]) - float(c0['matrix'][5])) < 1e-3 for c in chars_v)
    font = 'Bol' if 'Bol' in c0['fontname'] else 'Reg'
    col = c0.get('non_stroking_color'); colr = INK if not col else (tuple(col[:3]) if len(col) >= 3 else (col[0],)*3)
    n0 = len(pg.chars); right = box[2]
    rect = (box[0]-0.5, box[1]-0.5, box[2]+0.5, box[3]+0.5)
    pdf = pikepdf.open(src); removed = clear(pdf, H, [rect])
    tmp = Path(tempfile.mkdtemp()); pdf.save(tmp/'cleared.pdf')
    assert len(pdfplumber.open(tmp/'cleared.pdf').pages[0].chars) == n0 - len(chars), 'removal touched other glyphs'
    t = Top(W, H); t.text(right, base, new, size, font=font, color=colr, ha='right'); t.save(tmp/'label.pdf')
    pdf = pikepdf.open(tmp/'cleared.pdf'); overlay(pdf, tmp/'label.pdf', 0, 0); pdf.save(dst)
    pg2, box2, chars2 = find(dst, new)
    assert box2[0] >= box[0] - 0.5 and abs(box2[2] - right) < 0.6, (box, box2)
    a, b = render(src, tmp/'a.png'), render(dst, tmp/'b.png'); k = 150/72
    diff = ImageChops.difference(a, b)
    x0, y0, x1, y1 = [int(v*k) for v in (rect[0]-1, rect[1]-1, rect[2]+1, rect[3]+1)]
    mask = Image.new('L', a.size, 255); mask.paste(0, (x0, y0, x1, y1))
    outside = ImageChops.multiply(diff.convert('L'), mask).getbbox()
    assert outside is None, ('pixels changed outside the label box', outside)
    return dict(source=str(src), output=str(dst), old=old, new=new, font=c0['fontname'], size=size,
                baseline_top=round(base, 2), right_edge=round(right, 2), removed_glyphs=len(chars),
                new_box=[round(v, 2) for v in box2], changed_pixels_outside_box=0)


if __name__ == '__main__':
    print(json.dumps(main(*sys.argv[1:5]), indent=1))
