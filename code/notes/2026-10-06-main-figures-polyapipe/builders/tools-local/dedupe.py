"""Lossless size reduction for the re-tiled finals: identical streams (fonts, images, forms copied once per
panel segment) are merged into one object, then the file is rewritten with object streams.
Usage: python3 dedupe.py [figure numbers]"""
import sys, hashlib
from pathlib import Path
import pikepdf
HERE = Path(__file__).resolve().parent


def key(o):
    # stream dictionary serialised with indirect references as 'n g R', so two forms are equal only if they
    # point at the same (already merged) fonts, images and child forms
    parts = []
    for k in sorted(o.keys()):
        if k == '/Length': continue
        v = o[k]
        parts.append(k+'='+(f'{v.objgen[0]} {v.objgen[1]} R' if isinstance(v, pikepdf.Object) and v.is_indirect
                             else (v.unparse().decode('latin-1') if isinstance(v, pikepdf.Object) else repr(v))))
    return hashlib.sha256(o.read_raw_bytes()+'|'.join(parts).encode('latin-1')).hexdigest()


def dedupe(path):
    pdf = pikepdf.open(path, allow_overwriting_input=True)
    for _ in range(8):                      # forms that differ only by duplicated children merge on the next pass
        first, dup = {}, {}
        for o in pdf.objects:
            if isinstance(o, pikepdf.Stream):
                k = key(o)
                if k in first: dup[o.objgen] = first[k]
                else: first[k] = o
        if not dup: break
        def fix(c):
            if isinstance(c, pikepdf.Dictionary) or isinstance(c, pikepdf.Stream):
                for k_ in list(c.keys()):
                    v = c[k_]
                    if not isinstance(v, pikepdf.Object): continue
                    if v.is_indirect and v.objgen in dup: c[k_] = dup[v.objgen]
                    elif isinstance(v, (pikepdf.Dictionary, pikepdf.Array)) and not v.is_indirect: fix(v)
            elif isinstance(c, pikepdf.Array):
                for i in range(len(c)):
                    v = c[i]
                    if not isinstance(v, pikepdf.Object): continue
                    if v.is_indirect and v.objgen in dup: c[i] = dup[v.objgen]
                    elif isinstance(v, (pikepdf.Dictionary, pikepdf.Array)) and not v.is_indirect: fix(v)
        for o in list(pdf.objects):
            if isinstance(o, (pikepdf.Dictionary, pikepdf.Stream, pikepdf.Array)): fix(o)
        fix(pdf.trailer)
        pdf.remove_unreferenced_resources() if False else None
    pdf.save(path, object_stream_mode=pikepdf.ObjectStreamMode.generate, compress_streams=True, recompress_flate=True)


if __name__ == '__main__':
    for f in [int(v) for v in sys.argv[1:]] or [1, 2, 3, 4, 5]:
        p = HERE.parent/f'figure{f}'/f'Figure.{f}.pdf'
        b = p.stat().st_size; dedupe(p); print(f, b, '->', p.stat().st_size)
