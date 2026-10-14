import sys, pdfplumber, json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
for n in (1,2,3,4):
    pg = pdfplumber.open(ROOT/f'Fig/Figure.{n}.pdf').pages[0]
    ws = pg.extract_words(x_tolerance=1.5, keep_blank_chars=True, extra_attrs=['size','fontname'])
    # merge words on same line (top within 0.6) and gap < 4pt
    ws.sort(key=lambda w:(round(w['top']),w['x0']))
    lines=[]
    for w in ws:
        if lines and abs(lines[-1]['top']-w['top'])<0.6 and w['x0']-lines[-1]['x1']<4.5 and abs(lines[-1]['size']-w['size'])<0.3:
            L=lines[-1]; L['text']+=' '+w['text']; L['x1']=w['x1']
        else: lines.append(dict(text=w['text'],x0=w['x0'],x1=w['x1'],top=w['top'],size=w['size'],font=w['fontname'].split('+')[-1]))
    print(f'##### Figure {n}  {pg.width:.1f}x{pg.height:.1f}')
    for L in lines:
        print(f"{L['top']:6.1f} {L['x0']:6.1f}-{L['x1']:6.1f} {L['size']:4.1f} {L['font'][-10:]:>10} | {L['text']}")
