"""Region-scoped vector editing for single-page figure PDFs (pikepdf only).

remove_regions(pdf, page, rects): drops painted paths, images and text-show operators whose
page-space extent lies inside any rect (x0, y0, x1, y1 in PDF points, origin bottom-left).
Form XObjects that straddle a rect are copied and edited recursively, so nothing outside
the rects changes. Clip-only paths are kept. No rasterization.
"""
import pikepdf
from pikepdf import Operator, Name

PATH_BUILD = {'m', 'l', 'c', 'v', 'y', 'h', 're'}
PAINT = {'S', 's', 'f', 'F', 'f*', 'B', 'B*', 'b', 'b*', 'n'}
SHOW = {'Tj', 'TJ', "'", '"'}


def mul(a, b):  # a then b (row-vector convention, PDF style)
    return (a[0]*b[0]+a[1]*b[2], a[0]*b[1]+a[1]*b[3],
            a[2]*b[0]+a[3]*b[2], a[2]*b[1]+a[3]*b[3],
            a[4]*b[0]+a[5]*b[2]+b[4], a[4]*b[1]+a[5]*b[3]+b[5])


def apply(m, x, y):
    return (m[0]*x+m[2]*y+m[4], m[1]*x+m[3]*y+m[5])


def bbox(pts):
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    return (min(xs), min(ys), max(xs), max(ys))


def inside(bb, rects, tol=0.6):
    return any(bb[0] >= r[0]-tol and bb[1] >= r[1]-tol and bb[2] <= r[2]+tol and bb[3] <= r[3]+tol for r in rects)


def overlaps(bb, rects):
    return any(not (bb[2] < r[0] or bb[0] > r[2] or bb[3] < r[1] or bb[1] > r[3]) for r in rects)


def f(x):
    return float(x)


def intersect(a, b):
    return (max(a[0], b[0]), max(a[1], b[1]), min(a[2], b[2]), min(a[3], b[3]))


BIG = (-1e6, -1e6, 1e6, 1e6)


class Editor:
    def __init__(self, pdf, rects, keep=None, log=None, center=False):
        self.pdf = pdf; self.rects = rects; self.keep = keep; self.removed = []
        self.log = log; self.center = center

    def drop(self, bb):
        if self.center:
            cx, cy = (bb[0]+bb[2])/2, (bb[1]+bb[3])/2
            hit = any(r[0] <= cx < r[2] and r[1] <= cy < r[3] for r in (self.keep if self.keep is not None else self.rects))
            return (not hit) if self.keep is not None else hit
        if self.keep is not None:
            return not inside(bb, self.keep)
        return inside(bb, self.rects)

    def touches(self, bb):
        if self.keep is not None:
            return True
        return overlaps(bb, self.rects)

    def edit_stream(self, ops, resources, base_ctm):
        out = []; ctm = base_ctm; stack = []; clip = getattr(self, '_clip', BIG)
        path = []; pathops = []
        tm = lm = (1, 0, 0, 1, 0, 0); tl = 0; fs = 1
        i = 0
        while i < len(ops):
            operands, op = ops[i]; o = str(op)
            if o == 'q':
                stack.append((ctm, clip)); out.append((operands, op))
            elif o == 'Q':
                ctm, clip = stack.pop() if stack else (base_ctm, BIG); out.append((operands, op))
            elif o == 'cm':
                ctm = mul(tuple(map(f, operands)), ctm); out.append((operands, op))
            elif o in PATH_BUILD:
                pathops.append((operands, op))
                v = list(map(f, operands))
                if o == 're':
                    x, y, w, h = v; path += [apply(ctm, x, y), apply(ctm, x+w, y+h), apply(ctm, x+w, y), apply(ctm, x, y+h)]
                else:
                    for k in range(0, len(v), 2):
                        path.append(apply(ctm, v[k], v[k+1]))
            elif o in ('W', 'W*'):
                pathops.append((operands, op))
            elif o in PAINT:
                pathops.append((operands, op))
                isclip = any(str(p[1]) in ('W', 'W*') for p in pathops)
                eff = intersect(bbox(path), clip) if path else None
                if path and o != 'n' and (eff[2] < eff[0] or eff[3] < eff[1] or self.drop(eff)):
                    self.removed.append(('path', eff))
                    if isclip:  # keep clipping effect, drop painting
                        out += [p for p in pathops[:-1]] + [([], Operator('n'))]
                else:
                    out += pathops
                if isclip and path:
                    clip = intersect(clip, bbox(path))
                path = []; pathops = []
            elif o == 'BT':
                tm = lm = (1, 0, 0, 1, 0, 0); out.append((operands, op))
            elif o == 'Tf':
                fs = f(operands[1]); out.append((operands, op))
            elif o == 'TL':
                tl = f(operands[0]); out.append((operands, op))
            elif o == 'Tm':
                tm = lm = tuple(map(f, operands)); out.append((operands, op))
            elif o in ('Td', 'TD'):
                tx, ty = map(f, operands)
                if o == 'TD': tl = -ty
                lm = mul((1, 0, 0, 1, tx, ty), lm); tm = lm; out.append((operands, op))
            elif o == 'T*':
                lm = mul((1, 0, 0, 1, 0, -tl), lm); tm = lm; out.append((operands, op))
            elif o in SHOW:
                if o in ("'", '"'):
                    lm = mul((1, 0, 0, 1, 0, -tl), lm); tm = lm
                trm = mul(tm, ctm)
                x0, y0 = apply(trm, 0, 0)
                n = 0
                for a in operands:
                    if isinstance(a, pikepdf.Array):
                        n += sum(len(bytes(e)) for e in a if isinstance(e, pikepdf.String))
                    elif isinstance(a, pikepdf.String):
                        n += len(bytes(a))
                xa, ya = apply(trm, 0.45*max(n, 1)*(0.5 if n > 1 and self._twobyte(resources) else 1), 0.7)
                bb = bbox([(x0, y0), (xa, ya)])
                if self.drop((x0, y0, x0, y0)):
                    self.removed.append(('text', bb, str(operands)[:60]))
                    if o in ("'", '"'):
                        out.append(([], Operator('T*')))
                else:
                    out.append((operands, op))
                # advance approximation (only matters for consecutive shows)
                tm = mul((1, 0, 0, 1, 0.5*max(n, 1)*fs, 0), tm)
            elif o == 'Do':
                name = operands[0]; xo = resources.XObject[name]
                st = str(xo.get('/Subtype'))
                if st == '/Image':
                    bb = bbox([apply(ctm, x, y) for x, y in ((0, 0), (1, 0), (0, 1), (1, 1))])
                    effb = intersect(bb, clip)
                    if self.drop(effb):
                        self.removed.append(('image', bb))
                    else:
                        if self.keep is None and overlaps(bb, self.rects) and not inside(bb, self.rects):
                            self.removed.append(('image-partial-kept', bb))
                        out.append((operands, op))
                else:
                    m = tuple(map(f, xo.get('/Matrix', [1, 0, 0, 1, 0, 0])))
                    fctm = mul(m, ctm)
                    b = list(map(f, xo.BBox))
                    bb = bbox([apply(fctm, b[0], b[1]), apply(fctm, b[2], b[3]), apply(fctm, b[0], b[3]), apply(fctm, b[2], b[1])])
                    effb = intersect(bb, clip)
                    if self.center:
                        reg = self.keep if self.keep is not None else self.rects
                        if not overlaps(effb, reg):
                            if self.keep is not None: self.removed.append(('form', bb))
                            else: out.append((operands, op))
                            i += 1; continue
                        self._clip = clip
                        newname = self._edit_form(resources, name, xo, fctm)
                        self._clip = BIG
                        out.append(([newname], op))
                    elif self.keep is not None:
                        if inside(effb, self.keep):
                            out.append((operands, op)); i += 1; continue
                        if not overlaps(effb, self.keep):
                            self.removed.append(('form', bb)); i += 1; continue
                        self._clip = clip
                        newname = self._edit_form(resources, name, xo, fctm)
                        self._clip = BIG
                        out.append(([newname], op))
                    elif self.drop(effb):
                        self.removed.append(('form', bb))
                    elif self.touches(effb):
                        self._clip = clip
                        newname = self._edit_form(resources, name, xo, fctm)
                        self._clip = BIG
                        out.append(([newname], op))
                    else:
                        out.append((operands, op))
            elif o == 'INLINE IMAGE':
                bb = bbox([apply(ctm, x, y) for x, y in ((0, 0), (1, 0), (0, 1), (1, 1))])
                if self.drop(intersect(bb, clip)):
                    self.removed.append(('inline-image', bb))
                else:
                    out.append((operands, op))
            else:
                out.append((operands, op))
            i += 1
        return out

    def _twobyte(self, resources):
        return False

    def _edit_form(self, resources, name, xo, fctm):
        res = xo.get('/Resources', resources)
        ops = pikepdf.parse_content_stream(xo)
        new_ops = self.edit_stream(ops, res, fctm)
        data = pikepdf.unparse_content_stream(new_ops)
        newxo = pikepdf.Stream(self.pdf, data)
        for k, v in xo.items():
            if k not in ('/Length', '/Filter', '/DecodeParms'):
                newxo[k] = v
        newname = Name(str(name) + 'e')
        while newname in resources.XObject:
            newname = Name(str(newname) + 'e')
        resources.XObject[newname] = self.pdf.make_indirect(newxo)
        return newname


def keep_only(pdf, rect, page_index=0, center=False):
    page = pdf.pages[page_index]
    ed = Editor(pdf, [], keep=(rect if isinstance(rect[0], (list, tuple)) else [rect]), center=center)
    ops = pikepdf.parse_content_stream(page)
    new = ed.edit_stream(ops, page.Resources, (1, 0, 0, 1, 0, 0))
    page.Contents = pdf.make_stream(pikepdf.unparse_content_stream(new))
    return ed.removed


def overlay_shifted(pdf, part_pdf_path, dx, dy, page_index=0):
    src = pikepdf.open(part_pdf_path); sp = src.pages[0]; mb = [float(v) for v in sp.MediaBox]
    pdf.pages[page_index].add_overlay(sp, pikepdf.Rectangle(mb[0]+dx, mb[1]+dy, mb[2]+dx, mb[3]+dy))


def remove_regions(pdf, rects, page_index=0, center=False):
    page = pdf.pages[page_index]
    ed = Editor(pdf, rects, center=center)
    ops = pikepdf.parse_content_stream(page)
    new = ed.edit_stream(ops, page.Resources, (1, 0, 0, 1, 0, 0))
    page.Contents = pdf.make_stream(pikepdf.unparse_content_stream(new))
    return ed.removed


def overlay(pdf, panel_pdf_path, x0, y0, page_index=0):
    """Place page 1 of panel_pdf_path with its lower-left at (x0, y0) at 1:1 scale."""
    src = pikepdf.open(panel_pdf_path)
    sp = src.pages[0]; mb = [float(v) for v in sp.MediaBox]
    w = mb[2]-mb[0]; h = mb[3]-mb[1]
    page = pdf.pages[page_index]
    page.add_overlay(sp, pikepdf.Rectangle(x0, y0, x0+w, y0+h))
    return (x0, y0, x0+w, y0+h)


def recanvas(src_path, out_path, grow_top=0.0, grow_bottom=0.0, shift_x=0.0):
    """Place the page on a new canvas: +grow_top adds space above, +grow_bottom adds space below
    (negative values crop). Vector content is preserved via a form XObject."""
    src = pikepdf.open(src_path); sp = src.pages[0]; mb = [float(v) for v in sp.MediaBox]
    w, h = mb[2]-mb[0], mb[3]-mb[1]
    new = pikepdf.new(); page = new.add_blank_page(page_size=(w, h+grow_top+grow_bottom))
    pad = 200.0   # widen the source box so content moved beyond the old canvas is not clipped
    sp.MediaBox = [mb[0]-pad, mb[1]-pad, mb[2]+pad, mb[3]+pad]
    if '/CropBox' in sp: del sp['/CropBox']
    page.add_overlay(sp, pikepdf.Rectangle(shift_x-pad, grow_bottom-pad, shift_x+w+pad, grow_bottom+h+pad), shrink=False, expand=False)
    new.save(out_path)
