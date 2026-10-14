import subprocess, tempfile
from pathlib import Path
from PIL import Image, ImageDraw
HERE = Path(__file__).resolve().parents[1]; ROOT = HERE.parents[1]
DPI = 220; k = DPI/72
REG = [('Figure.1', 0, 325, 232, 415), ('Figure.4', 0, 100, 453.5, 190), ('SuppFig.S19', 290, 200, 453.54, 312)]
rows = []
for f, x0, y0, x1, y1 in REG:
    pair = []
    for p in (ROOT/f'Fig/{f}.pdf', HERE/f'out/{f}.pdf'):
        o = Path(tempfile.mkdtemp())/'c'
        subprocess.run(['pdftoppm', '-r', str(DPI), '-png', '-singlefile', '-x', str(int(x0*k)), '-y', str(int(y0*k)), '-W', str(int((x1-x0)*k)), '-H', str(int((y1-y0)*k)), str(p), str(o)], check=True, capture_output=True)
        pair.append(Image.open(str(o)+'.png').convert('RGB'))
    rows.append((f, pair))
blocks = []
for f, (a, b) in rows:
    if a.width > 900:   # wide crop: stack before over after
        im = Image.new('RGB', (a.width, a.height+b.height+30), 'white'); im.paste(a, (0, 0)); im.paste(b, (0, a.height+30))
        ImageDraw.Draw(im).line([(0, a.height+15), (a.width, a.height+15)], fill=(200, 0, 0), width=3)
    else:
        im = Image.new('RGB', (a.width+b.width+30, max(a.height, b.height)), 'white'); im.paste(a, (0, 0)); im.paste(b, (a.width+30, 0))
        ImageDraw.Draw(im).line([(a.width+15, 0), (a.width+15, im.height)], fill=(200, 0, 0), width=3)
    blocks.append(im)
W = max(b.width for b in blocks); H = sum(b.height+20 for b in blocks)
sheet = Image.new('RGB', (W, H), 'white'); y = 0
for b in blocks: sheet.paste(b, (0, y)); y += b.height+20
sheet.save(HERE/'work/contact.jpg', quality=82); print(sheet.size)
