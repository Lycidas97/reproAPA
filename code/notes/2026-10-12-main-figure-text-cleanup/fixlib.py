"""Small in-place edits on the production main figures.

Copied from notes/2026-10-11-caption-figure-fixes/fixlib.py on 2026-10-12 (main-figure text cleanup) so that
outputs land in this folder; additions: Job.move_word (re-place one word by dy), Job.run(shifts=...) (move a
vector block via v5draw.shift_block after the removals), and enlarge_bottom (add white space under a page).

Same drawing path as the earlier in-place relabels: pdfedit.remove_regions (via v5draw.clear) deletes the glyphs
being replaced, v5draw.Top draws the new text/marks in Nimbus Sans L, pdfedit.overlay places them 1:1, and
pdfclean.finalize writes the file. recolor() changes the colour of marks inside a region without moving anything:
it wraps the affected paint operators (and geometry-only marker forms) with a colour change and restores the
original colour right after, so nothing outside the region and no geometry changes.
"""
import sys, json, hashlib, tempfile, subprocess, collections
from pathlib import Path
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'notes/2026-10-06-main-figures-polyapipe/builders/tools-local'))
from v5draw import Top, clear, shift_block, INK, MUTED, AXC, NULLBAR, textwidth   # noqa
from pdfedit import overlay, mul, apply, bbox, PATH_BUILD, PAINT, SHOW   # noqa
import pikepdf, pdfplumber
from pikepdf import Operator, Name

sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
OUT = HERE/'out'; OUT.mkdir(exist_ok=True)
PROX = '#0072B2'; DIST = '#D55E00'


def hex2rgb(h):
    h = h.lstrip('#'); return tuple(int(h[i:i+2], 16)/255 for i in (0, 2, 4))


def words(path, box=None):
    pg = pdfplumber.open(path).pages[0]
    if box: pg = pg.crop(box, strict=False)
    return [(w['text'], round(w['x0'], 1), round(w['top'], 1)) for w in pg.extract_words(x_tolerance=1.5)]


def chars_of(path, text, box):
    """Chars of the first word equal to `text` inside box (top-origin x0, top, x1, bottom)."""
    pg = pdfplumber.open(path).pages[0]
    c = pg.crop(box, strict=False)
    ws = [w for w in c.extract_words(x_tolerance=1.5, extra_attrs=['size', 'fontname']) if w['text'] == text]
    assert len(ws) == 1, (text, box, [w['text'] for w in c.extract_words(x_tolerance=1.5)])
    w = ws[0]
    ch = [k for k in c.chars if k['x0'] >= w['x0']-0.01 and k['x1'] <= w['x1']+0.01 and abs(k['top']-w['top']) < 0.5]
    first = min(ch, key=lambda k: k['x0'])
    H = float(pg.height)
    col = first.get('non_stroking_color')
    return dict(x0=w['x0'], x1=w['x1'], top=w['top'], bottom=w['bottom'], size=round(float(first['size']), 2),
                font=first['fontname'], origin=float(first['matrix'][4]), base=H-float(first['matrix'][5]),
                color=tuple(col[:3]) if col and len(col) >= 3 else INK)


def font_key(fontname):
    f = fontname.split('+')[-1]
    if 'BolIta' in f or 'BoldItal' in f: return 'BolIta'
    if 'RegIta' in f or 'Ital' in f: return 'RegIta'
    if 'Bol' in f or 'Bold' in f: return 'Bol'
    return 'Reg'


class Recolor:
    COLOR_OPS = {'rg', 'RG', 'g', 'G', 'k', 'K', 'sc', 'scn', 'SC', 'SCN', 'cs', 'CS'}

    def __init__(self, pdf, rects, mapping, tol=0.004):
        self.pdf, self.rects, self.tol = pdf, rects, tol
        self.map = [(hex2rgb(a), hex2rgb(b)) for a, b in mapping.items()]
        self.log = collections.Counter()

    def target(self, c):
        if c is None: return None
        for a, b in self.map:
            if all(abs(x-y) <= self.tol for x, y in zip(c, a)): return b
        return None

    def hit(self, bb):
        cx, cy = (bb[0]+bb[2])/2, (bb[1]+bb[3])/2
        return any(r[0] <= cx <= r[2] and r[1] <= cy <= r[3] for r in self.rects)

    def overlaps(self, bb):
        return any(not (bb[2] < r[0] or bb[0] > r[2] or bb[3] < r[1] or bb[1] > r[3]) for r in self.rects)

    def pure(self, xo):
        try: ops = pikepdf.parse_content_stream(xo)
        except Exception: return False
        return not any(str(o) in self.COLOR_OPS or str(o) in ('Do', 'BT') for _, o in ops)

    def edit(self, ops, resources, ctm, fc, sc):
        out = []; stack = []; path = []; pathops = []
        tm = lm = (1, 0, 0, 1, 0, 0); tl = 0
        for operands, op in ops:
            o = str(op)
            if o == 'q': stack.append((ctm, fc, sc)); out.append((operands, op)); continue
            if o == 'Q':
                if stack: ctm, fc, sc = stack.pop()
                out.append((operands, op)); continue
            if o == 'cm': ctm = mul(tuple(map(float, operands)), ctm); out.append((operands, op)); continue
            if o in self.COLOR_OPS:
                v = [float(x) for x in operands if isinstance(x, (int, float)) or hasattr(x, '__float__')] if o not in ('cs', 'CS') else []
                col = tuple(v) if len(v) == 3 else ((v[0],)*3 if len(v) == 1 and o in ('g', 'G') else None)
                if o in ('cs',): col = (0.0, 0.0, 0.0)
                if o in ('CS',): col = (0.0, 0.0, 0.0)
                if o in ('rg', 'g', 'k', 'sc', 'scn', 'cs'): fc = col
                else: sc = col
                out.append((operands, op)); continue
            if o in PATH_BUILD:
                pathops.append((operands, op)); v = list(map(float, operands))
                if o == 're':
                    x, y, w, h = v; path += [apply(ctm, x, y), apply(ctm, x+w, y+h)]
                else:
                    path += [apply(ctm, v[k], v[k+1]) for k in range(0, len(v), 2)]
                continue
            if o in ('W', 'W*'): pathops.append((operands, op)); continue
            if o in PAINT:
                pathops.append((operands, op))
                pre, post = [], []
                if path and o != 'n' and self.hit(bbox(path)):
                    tf = self.target(fc) if o in ('f', 'F', 'f*', 'B', 'B*', 'b', 'b*') else None
                    ts = self.target(sc) if o in ('S', 's', 'B', 'B*', 'b', 'b*') else None
                    if tf: pre.append((list(tf), Operator('rg'))); post.append((list(fc), Operator('rg'))); self.log['path fill'] += 1
                    if ts: pre.append((list(ts), Operator('RG'))); post.append((list(sc), Operator('RG'))); self.log['path stroke'] += 1
                out += pre + pathops + post; path = []; pathops = []; continue
            if o == 'BT': tm = lm = (1, 0, 0, 1, 0, 0); out.append((operands, op)); continue
            if o == 'TL': tl = float(operands[0]); out.append((operands, op)); continue
            if o == 'Tm': tm = lm = tuple(map(float, operands)); out.append((operands, op)); continue
            if o in ('Td', 'TD'):
                tx, ty = map(float, operands)
                if o == 'TD': tl = -ty
                lm = mul((1, 0, 0, 1, tx, ty), lm); tm = lm; out.append((operands, op)); continue
            if o == 'T*': lm = mul((1, 0, 0, 1, 0, -tl), lm); tm = lm; out.append((operands, op)); continue
            if o in SHOW:
                x0, y0 = apply(mul(tm, ctm), 0, 0)
                tf = self.target(fc) if self.hit((x0, y0, x0, y0)) else None
                if tf:
                    out += [(list(tf), Operator('rg')), (operands, op), (list(fc), Operator('rg'))]; self.log['text'] += 1
                else:
                    out.append((operands, op))
                continue
            if o == 'Do':
                name = operands[0]; xo = resources.XObject[name]
                if str(xo.get('/Subtype')) != '/Form': out.append((operands, op)); continue
                m = tuple(map(float, xo.get('/Matrix', [1, 0, 0, 1, 0, 0]))); fctm = mul(m, ctm)
                b = list(map(float, xo.BBox))
                bb = bbox([apply(fctm, b[0], b[1]), apply(fctm, b[2], b[3]), apply(fctm, b[0], b[3]), apply(fctm, b[2], b[1])])
                if self.pure(xo):
                    tf, ts = (self.target(fc), self.target(sc)) if self.hit(bb) else (None, None)
                    pre = ([(list(tf), Operator('rg'))] if tf else []) + ([(list(ts), Operator('RG'))] if ts else [])
                    post = ([(list(fc), Operator('rg'))] if tf else []) + ([(list(sc), Operator('RG'))] if ts else [])
                    if pre: self.log['form mark'] += 1
                    out += pre + [(operands, op)] + post; continue
                if not self.overlaps(bb): out.append((operands, op)); continue
                res = xo.get('/Resources', resources)
                new_ops = self.edit(pikepdf.parse_content_stream(xo), res, fctm, fc, sc)
                nx = pikepdf.Stream(self.pdf, pikepdf.unparse_content_stream(new_ops))
                for k, v in xo.items():
                    if k not in ('/Length', '/Filter', '/DecodeParms'): nx[k] = v
                nn = Name(str(name)+'r')
                while nn in resources.XObject: nn = Name(str(nn)+'r')
                resources.XObject[nn] = self.pdf.make_indirect(nx)
                out.append(([nn], op)); continue
            out.append((operands, op))
        return out


def recolor(pdf, H, rects_top, mapping):
    rects = [(x0, H-y1, x1, H-y0) for x0, y0, x1, y1 in rects_top]
    rc = Recolor(pdf, rects, mapping)
    page = pdf.pages[0]
    new = rc.edit(pikepdf.parse_content_stream(page), page.Resources, (1, 0, 0, 1, 0, 0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
    page.Contents = pdf.make_stream(pikepdf.unparse_content_stream(new))
    return dict(rc.log)


class Job:
    """One figure: remove listed glyph runs, recolour regions, draw new content, write out/Figure.N.pdf."""
    def __init__(self, fig, src=None, name=None):
        self.fig = fig; self.src = Path(src) if src else ROOT/f'Fig/Figure.{fig}.pdf'
        self.name = name or f'Figure.{fig}'
        self.pdf = pikepdf.open(self.src)
        pg = pdfplumber.open(self.src).pages[0]; self.W, self.H = float(pg.width), float(pg.height)
        self.t = Top(self.W, self.H); self.clear_rects = []; self.rec = dict(source=str(self.src.relative_to(ROOT)), source_sha=sha(self.src), edits=[])
        self.expect_gone = []; self.expect_new = []; self._moves = []

    def replace(self, old, box, new, align='left', size=None, font=None, color=None, note=''):
        c = chars_of(self.src, old, box)
        self.clear_rects.append((c['x0']-0.2, c['top']-0.2, c['x1']+0.2, c['bottom']+0.2))
        size = size or c['size']; font = font or font_key(c['font']); color = color or c['color']
        if align == 'left': x = c['origin']
        elif align == 'center': x = (c['x0']+c['x1'])/2
        else: x = c['x1']
        self.t.text(x, c['base'], new, size, font=font, color=color, ha=align)
        self.rec['edits'].append(dict(kind='text', old=old, new=new, at=[round(c['x0'], 2), round(c['top'], 2)], size=size, font=font, note=note))
        self.expect_gone.append(old)
        return c

    def drop_word(self, text, box, note=''):
        c = chars_of(self.src, text, box)
        self.clear_rects.append((c['x0']-0.2, c['top']-0.2, c['x1']+0.2, c['bottom']+0.2))
        self.rec['edits'].append(dict(kind='drop text', old=text, at=[round(c['x0'], 2), round(c['top'], 2)], note=note))
        return c

    def move_word(self, text, box, dy, new=None, align='left', note=''):
        """Re-place one word (or a new text) dy pt lower, same size, font and colour."""
        c = chars_of(self.src, text, box)
        self.clear_rects.append((c['x0']-0.2, c['top']-0.2, c['x1']+0.2, c['bottom']+0.2))
        x = c['origin'] if align == 'left' else c['x1']
        self._moves.append((x, c['base']+dy, new or text, c['size'], font_key(c['font']), c['color'], align))
        self.rec['edits'].append(dict(kind='move text', old=text, new=new or text, dy=dy, at=[round(c['x0'], 2), round(c['top'], 2)], note=note))
        if new and new != text: self.expect_gone.append(text)
        return c

    def remove(self, rect, note):
        self.clear_rects.append(rect); self.rec['edits'].append(dict(kind='remove', rect=rect, note=note))

    def run(self, recolors=(), draw=None, center=True, shifts=()):
        before = collections.Counter(w for w, _, _ in words(self.src))
        if self.clear_rects:
            removed = clear(self.pdf, self.H, self.clear_rects, center=center)
            kinds = collections.Counter(r[0] for r in removed)
            self.rec['removed_kinds'] = dict(kinds)
        for rects, mapping, note in recolors:
            log = recolor(self.pdf, self.H, rects, mapping)
            self.rec['edits'].append(dict(kind='recolor', rects=rects, mapping=mapping, changed=log, note=note))
        for region, dy, note in shifts:
            n = shift_block(self.pdf, self.src, self.H, self.W, None, region, dy=dy)
            self.rec['edits'].append(dict(kind='shift block', region=list(region), dy=dy, removed_on_page=n, note=note))
        for x, base, s_, size, font, color, ha in self._moves:
            self.t.text(x, base, s_, size, font=font, color=color, ha=ha)
        if draw: draw(self.t)
        tmp = Path(tempfile.mkdtemp()); ov = tmp/'overlay.pdf'; self.t.save(ov)
        overlay(self.pdf, ov, 0, 0)
        out = OUT/f'{self.name}.pdf'; self.pdf.save(out)
        from pdfclean import finalize
        finalize(out, label=self.name, log=OUT/'sizes.tsv')
        after = collections.Counter(w for w, _, _ in words(out))
        self.rec['words_removed'] = dict(before-after); self.rec['words_added'] = dict(after-before)
        self.rec['out'] = str(out.relative_to(ROOT)); self.rec['sha256'] = sha(out)
        subprocess.run(['pdftoppm', '-r', '150', '-png', '-singlefile', str(out), str(OUT/self.name)], check=True, capture_output=True)
        return self.rec


def reorder_rows(src, out, moves):
    """Rebuild a page from vector crops: moves = [(region_top_origin, dy), ...]; regions must not overlap after the move.
    Each crop keeps only objects touching its region (objects entirely outside are deleted) and is clipped to it."""
    pg = pdfplumber.open(src).pages[0]; W, H = float(pg.width), float(pg.height)
    dst = pikepdf.new(); dst.add_blank_page(page_size=(W, H))
    placed = []
    for (x0, y0, x1, y1), dy in moves:
        blk = pikepdf.open(src)
        clear(blk, H, [(-20, -20, x0-0.05, H+20), (x1+0.05, -20, W+20, H+20), (-20, -20, W+20, y0-0.05), (-20, y1+0.05, W+20, H+20)])
        p = blk.pages[0]; box = [x0, H-y1, x1, H-y0]; p.MediaBox = box; p.CropBox = box
        ny0, ny1 = y0+dy, y1+dy
        assert ny0 >= -0.01 and ny1 <= H+0.01, (ny0, ny1)
        for a, b in placed: assert ny1 <= a+1e-6 or ny0 >= b-1e-6, ('overlap', (ny0, ny1), (a, b))
        placed.append((ny0, ny1))
        dst.pages[0].add_overlay(blk.pages[0], pikepdf.Rectangle(x0, H-ny1, x1, H-ny0))
    dst.save(out)
    return placed


def enlarge_bottom(path, extra, label):
    """Add `extra` pt of white space below the page content; top-origin coordinates stay unchanged."""
    pdf = pikepdf.open(path, allow_overwriting_input=True); pg = pdf.pages[0]
    x0, y0, x1, y1 = map(float, pg.MediaBox)
    pg.contents_coalesce()
    pre = pdf.make_stream(f'q 1 0 0 1 0 {extra} cm\n'.encode()); post = pdf.make_stream(b'\nQ\n')
    pg.Contents = pikepdf.Array([pre, pg.Contents, post])
    pg.MediaBox = [x0, y0, x1, y1+extra]
    if '/CropBox' in pg: pg.CropBox = [x0, y0, x1, y1+extra]
    pdf.save(path)
    from pdfclean import finalize
    finalize(path, label=label, log=OUT/'sizes.tsv')
    return [x0, y0, x1, y1+extra]
