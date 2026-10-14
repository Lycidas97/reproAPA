"""Before/after crops of each edited region (250 dpi) stacked into contact sheets."""
import subprocess, tempfile, sys
from pathlib import Path
from PIL import Image, ImageDraw
ROOT = Path(__file__).resolve().parents[3]; HERE = Path(__file__).resolve().parents[1]
DPI = 250; k = DPI/72
REG = [('S1', 0, 425, 200, 500), ('S2', 0, 360, 200, 445), ('S3', 300, 435, 453, 465), ('S8', 320, 115, 453, 145),
       ('S9', 0, 92, 200, 125), ('S9', 0, 235, 80, 285), ('S11', 0, 160, 260, 262), ('S13', 300, 395, 453, 425),
       ('S14', 230, 120, 300, 155), ('S14', 0, 305, 453, 333), ('S15', 0, 315, 330, 342), ('S17', 0, 110, 160, 150),
       ('S20', 280, 450, 453, 475)]
tmp = Path(tempfile.mkdtemp()); rows = []
for i, (f, x0, y0, x1, y1) in enumerate(REG):
    pair = []
    for tag, p in [('before', ROOT/f'Fig/SuppFig.{f}.pdf'), ('after', HERE/f'out/SuppFig.{f}.pdf')]:
        o = tmp/f'{i}{tag}'
        subprocess.run(['pdftoppm', '-r', str(DPI), '-png', '-singlefile', '-x', str(int(x0*k)), '-y', str(int(y0*k)),
                        '-W', str(int((x1-x0)*k)), '-H', str(int((y1-y0)*k)), str(p), str(o)], check=True, capture_output=True)
        pair.append(Image.open(str(o)+'.png').convert('RGB'))
    rows.append((f, pair))
sheets = [rows[:7], rows[7:]]
for j, rs in enumerate(sheets):
    W = max(a.width+b.width for _, (a, b) in rs)+30; H = sum(max(a.height, b.height)+28 for _, (a, b) in rs)
    im = Image.new('RGB', (W, H), 'white'); d = ImageDraw.Draw(im); y = 0
    for f, (a, b) in rs:
        d.text((4, y+4), f'{f}  before | after', fill=(200, 0, 0)); y += 22
        im.paste(a, (0, y)); im.paste(b, (a.width+30, y)); d.line([(a.width+15, y), (a.width+15, y+max(a.height, b.height))], fill=(200, 0, 0), width=2)
        y += max(a.height, b.height)+6
    im.save(HERE/f'work/contact{j+1}.png'); print(im.size)
