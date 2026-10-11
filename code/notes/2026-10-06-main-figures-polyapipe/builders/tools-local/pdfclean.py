"""Final-save clean-up for every figure PDF (standing rule; called by retile.py, never per figure by hand).

finalize(path):
  1. drop resource entries (Font, XObject, ExtGState, Shading) that no content stream uses, page and nested forms;
  2. merge identical streams (font programs, images, forms) into one shared object (dedupe.py logic, references
     compared as object numbers so only truly identical objects merge), repeated until stable;
  3. write with object streams and recompressed flate; only reachable objects are written, so anything removed
     in 1-2 (and anything left over from earlier passes) leaves the file;
  4. check that the page renders pixel-identical before and after (150 dpi) - the clean-up may never change the figure;
  5. append size and object counts to sizes.tsv and warn when the file grew > 15 % since the last entry.
"""
import sys, json, hashlib, subprocess, tempfile, time, collections
from pathlib import Path
import numpy as np
from PIL import Image
import pikepdf
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from dedupe import key

OPS = {'Tf': '/Font', 'Do': '/XObject', 'gs': '/ExtGState', 'sh': '/Shading'}


PAINT = {'f', 'F', 'f*', 'B', 'B*', 'b', 'b*', 'S', 's', 'sh', 'Tj', 'TJ', "'", '"', 'BI', 'INLINE IMAGE'}


def strip_invisible(pdf):
    """remove leftovers that draw nothing: text blocks without any shown text, and Do calls of forms that
    (recursively) paint nothing; returns counts"""
    memo, n = {}, collections.Counter()
    def paints(form, pres=None):
        if form.objgen in memo: return memo[form.objgen]
        memo[form.objgen] = True                      # guard against cycles
        res = form.get('/Resources', pres); xo = res.get('/XObject') if res is not None else None
        try: ops = pikepdf.parse_content_stream(form)
        except Exception: return True
        out = False
        for operands, op in ops:
            o = str(op)
            if o in PAINT: out = True; break
            if o == 'Do' and xo is not None and str(operands[0]) in xo:
                x = xo[str(operands[0])]
                if x.get('/Subtype') != '/Form' or paints(x, res): out = True; break
        memo[form.objgen] = out
        return out
    TEXT_STATE = {'Tf', 'Tc', 'Tw', 'Tz', 'TL', 'Tr', 'Ts'}
    POSITION = {'BT', 'ET', 'Td', 'TD', 'Tm', 'T*'}
    def rewrite(obj, pres=None):
        res = obj.get('/Resources', pres); xo = res.get('/XObject') if res is not None else None
        try: ops = pikepdf.parse_content_stream(obj)
        except Exception: return
        # state set inside a removed (empty) text block still applies to what follows, so it is carried forward
        # and re-emitted just before the next operator that could depend on it
        new, block, changed = [], None, False
        pend = [{}]
        def flush(text):
            for k_ in list(pend[-1].keys()):
                if text or k_ not in TEXT_STATE:
                    new.append(pend[-1].pop(k_))
        for operands, op in ops:
            o = str(op)
            if o == 'Do' and xo is not None and str(operands[0]) in xo:
                x = xo[str(operands[0])]
                if x.get('/Subtype') == '/Form' and not paints(x, res):
                    changed = True; n['empty_form_draws'] += 1; continue
            if o == 'BT':
                block = [(operands, op)]; continue
            if block is not None:
                block.append((operands, op))
                if o == 'ET':
                    if any(str(b_[1]) in ('Tj', 'TJ', "'", '"') for b_ in block):
                        flush(True); new.extend(block)
                    else:
                        changed = True; n['empty_text_blocks'] += 1
                        for opr_, op_ in block:
                            if str(op_) not in POSITION:
                                pend[-1].pop(str(op_), None); pend[-1][str(op_)] = (opr_, op_)
                    block = None
                continue
            if o == 'q':
                flush(False); new.append((operands, op)); pend.append(dict(pend[-1])); continue
            if o == 'Q':
                pend.pop() if len(pend) > 1 else None; new.append((operands, op)); continue
            flush(o == 'Do')
            new.append((operands, op))
        if block: new.extend(block)
        if changed:
            data = pikepdf.unparse_content_stream(new)
            if isinstance(obj, pikepdf.Stream): obj.write(data)
            else: obj.Contents = pdf.make_stream(data)
        if xo is not None:
            for k in xo.keys():
                x = xo[k]
                if x.get('/Subtype') == '/Form' and x.objgen not in done:
                    done.add(x.objgen); rewrite(x, res)
    done = set()
    for pg in pdf.pages: rewrite(pg.obj)
    return dict(n)


def prune_resources(pdf):
    """drop resource entries no content uses; a form without its own /Resources uses its caller's, so its names
    count as used by the caller"""
    seen = set(); removed = collections.Counter()
    def names(obj, res):
        used = collections.defaultdict(set)
        try:
            ops = pikepdf.parse_content_stream(obj)
        except Exception:
            return None
        xo = res.get('/XObject') if res is not None else None
        for operands, op in ops:
            cat = OPS.get(str(op))
            if cat and operands:
                used[cat].add(str(operands[0]))
            if str(op) == 'Do' and xo is not None and str(operands[0]) in xo:
                x = xo[str(operands[0])]
                if x.get('/Subtype') == '/Form' and '/Resources' not in x:
                    sub = names(x, res)
                    if sub is None: return None
                    for k_, v_ in sub.items(): used[k_] |= v_
        return used
    # resource dictionaries can be shared between forms (after merging), so usage is pooled per dictionary
    pool, dicts, bad = collections.defaultdict(set), {}, set()
    def dkey(owner, res, cat):
        d = res[cat]
        if d.is_indirect: return ('i', d.objgen)
        if res.is_indirect: return ('r', res.objgen, cat)
        return ('o', owner.objgen, cat)
    def walk(obj):
        if obj.objgen != (0, 0):
            if obj.objgen in seen: return
            seen.add(obj.objgen)
        res = obj.get('/Resources')
        if res is None: return
        used = names(obj, res)
        for cat in OPS.values():
            if cat in res:
                k_ = dkey(obj, res, cat); dicts[k_] = (res[cat], cat)
                if used is None or obj.objgen == (0, 0) and False: bad.add(k_)
                if used is None: bad.add(k_)
                else: pool[k_] |= used[cat]
        if '/XObject' in res:
            for k in res['/XObject'].keys():
                x = res['/XObject'][k]
                if x.get('/Subtype') == '/Form': walk(x)
    for pg in pdf.pages: walk(pg.obj)
    for k_, (d, cat) in dicts.items():
        if k_ in bad: continue
        for name in list(d.keys()):
            if name not in pool[k_]:
                del d[name]; removed[cat] += 1
    return dict(removed)


def ser(v, top=False):
    """serialise with indirect references as 'n g R' (the top object itself is expanded)"""
    if isinstance(v, pikepdf.Object) and v.is_indirect and not top:
        return f'{v.objgen[0]} {v.objgen[1]} R'
    if isinstance(v, (pikepdf.Dictionary, pikepdf.Stream)):
        return '<<'+' '.join(k+' '+ser(v[k]) for k in sorted(v.keys()) if k != '/Length')+'>>'
    if isinstance(v, pikepdf.Array):
        return '['+' '.join(ser(x) for x in v)+']'
    if isinstance(v, pikepdf.Object):
        return v.unparse().decode('latin-1')
    return repr(v)


SKIP_TYPES = {'/Page', '/Pages', '/Catalog'}


def merge_identical(pdf):
    """merge identical streams and identical indirect dictionaries/arrays (fonts, widths, resources, forms)"""
    merged = 0
    for _ in range(12):
        first, dup = {}, {}
        for o in pdf.objects:
            if isinstance(o, pikepdf.Stream):
                k = 'S'+key(o)
            elif isinstance(o, (pikepdf.Dictionary, pikepdf.Array)) and o.is_indirect:
                if isinstance(o, pikepdf.Dictionary) and str(o.get('/Type')) in SKIP_TYPES: continue
                k = 'D'+hashlib.sha256(ser(o, top=True).encode('latin-1')).hexdigest()
            else:
                continue
            if k in first: dup[o.objgen] = first[k]
            else: first[k] = o
        if not dup: break
        merged += len(dup)
        def fix(c):
            items = list(c.items()) if isinstance(c, (pikepdf.Dictionary, pikepdf.Stream)) else list(enumerate(c))
            for k_, v in items:
                if not isinstance(v, pikepdf.Object): continue
                if v.is_indirect and v.objgen in dup: c[k_] = dup[v.objgen]
                elif not v.is_indirect and isinstance(v, (pikepdf.Dictionary, pikepdf.Array)): fix(v)
        for o in list(pdf.objects):
            if isinstance(o, (pikepdf.Dictionary, pikepdf.Stream, pikepdf.Array)): fix(o)
    return merged


def render(path, dpi=150):
    tmp = Path(tempfile.mkdtemp())/'r'
    subprocess.run(['pdftoppm', '-r', str(dpi), '-png', '-singlefile', str(path), str(tmp)], check=True, capture_output=True)
    return np.asarray(Image.open(str(tmp)+'.png').convert('RGB'), dtype=int)


def clean_to(path, out, steps=('strip', 'prune', 'merge')):
    pdf = pikepdf.open(path)
    # merge first (cheap and it collapses the per-panel copies), then strip and prune, then merge what became equal
    merged = merge_identical(pdf) if 'merge' in steps else 0
    if merged:
        buf = __import__('io').BytesIO(); pdf.save(buf); pdf = pikepdf.open(__import__('io').BytesIO(buf.getvalue()))
    stripped = strip_invisible(pdf) if 'strip' in steps else {}
    pruned = prune_resources(pdf) if 'prune' in steps else {}
    merged += merge_identical(pdf) if 'merge' in steps else 0
    pdf.save(out, object_stream_mode=pikepdf.ObjectStreamMode.generate, compress_streams=True, recompress_flate=True)
    return stripped, pruned, merged


def finalize(path, label=None, log=None):
    path = Path(path); before = render(path); b0 = path.stat().st_size
    tmp = Path(tempfile.mkdtemp())/'clean.pdf'            # outside the project folder
    stripped, pruned, merged = clean_to(path, tmp)
    after = render(tmp)
    if not (before.shape == after.shape and int(np.abs(before-after).max()) == 0):
        raise AssertionError(f'clean-up would change the render of {path}; file left unchanged')
    import shutil; shutil.copyfile(tmp, path)
    p2 = pikepdf.open(path)
    streams = [o for o in p2.objects if isinstance(o, pikepdf.Stream)]
    n_img = sum(1 for o in streams if o.get('/Subtype') == '/Image')
    n_font = sum(1 for o in streams if any(k in o for k in ('/Length1', '/Length2', '/Length3')) or o.get('/Subtype') in ('/Type1C', '/CIDFontType0C', '/OpenType'))
    row = dict(time=time.strftime('%Y-%m-%d %H:%M'), file=label or path.name, bytes=path.stat().st_size, bytes_before_clean=b0,
               objects=len(p2.objects), streams=len(streams), images=n_img, font_programs=n_font,
               stripped_invisible=sum(stripped.values()), pruned_resource_entries=sum(pruned.values()), merged_objects=merged)
    log = Path(log or HERE.parent/'sizes.tsv'); prev = None
    if log.exists():
        rows = [l.rstrip('\n').split('\t') for l in open(log)]
        hdr = rows[0]; same = [dict(zip(hdr, r)) for r in rows[1:] if dict(zip(hdr, r)).get('file') == row['file']]
        prev = int(same[-1]['bytes']) if same else None
    else:
        open(log, 'w').write('\t'.join(row)+'\n')
    open(log, 'a').write('\t'.join(str(v) for v in row.values())+'\n')
    row['grew_vs_last'] = None if prev is None else round(row['bytes']/prev-1, 3)
    if prev is not None and row['bytes'] > 1.15*prev:
        print(f'WARNING {row["file"]}: {prev} -> {row["bytes"]} bytes (+{100*(row["bytes"]/prev-1):.0f} %) since the last export', file=sys.stderr)
    return row


if __name__ == '__main__':
    for p in sys.argv[1:]:
        print(json.dumps(finalize(p)))
