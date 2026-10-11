"""Re-tile a figure so every visible gutter is the same: G between adjacent panels, M at the page edge.

Input  figureN/Figure.N.draft.pdf (the build output, one slot per panel as in audit.PANELS).
Output figureN/Figure.N.pdf, figureN/Figure.N.panels.json (final panel boxes = ink box + G/2), retile.json.

Each panel is copied (only its own objects; everything else is removed from the copy), clipped to its
ink box + 0.8 pt and translated. Horizontal and vertical gaps are set between ink boxes of panels that
share a slot boundary; the page width is unchanged, the page height follows the content.
Values (same for all five figures): G = 11.3 pt (4.0 mm), M = 4.5 pt (1.6 mm).
"""
import sys, json, subprocess, tempfile
from pathlib import Path
import numpy as np
from PIL import Image
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from audit import PANELS
from v5draw import clear, prune_xobjects
import pikepdf, pdfplumber

G, M = 11.3, 4.5
# Panel letters (round 5): the letter is not part of the box used for alignment. Each panel's content box
# (ink without the letter) is tiled with gaps G and equal heights; the letter is then placed with its ink
# top-left on the content box's top-left corner, so it never sticks out of the panel. The build must keep that
# corner free (LETTER_SLOT); retile.py --check lists every panel whose content enters it.
SLOT_PAD = (2.5, 2.0)
# panels that must share top and bottom edges: (stack one side, stack other side); the shorter side is filled by
# widening its existing all-white horizontal bands (content is translated, never scaled)
EQUAL = {1: [(['B', 'C'], ['D'])], 2: [(['B'], ['C']), (['E'], ['F'])], 3: [(['A'], ['B']), (['D'], ['E'])],
         4: [(['A'], ['B']), (['C', 'E'], ['D'])], 5: [(['A', 'B'], ['C']), (['D'], ['E'])]}
NBANDS = 4          # at most this many bands share one panel's extra height


def slot_of_center(r, P):
    cx, cy = (r[0]+r[2])/2, (r[1]+r[3])/2
    hit = [q for q, (qa, qb, qc, qd) in P.items() if qa <= cx <= qc and qb <= cy <= qd]
    return hit[0] if hit else None


def letters(pdf, P):
    """panel letter of each slot: the bold glyph >= 9.5 pt equal to the panel id nearest the slot's top-left"""
    pg = pdfplumber.open(pdf).pages[0]; out = {}
    for p, (a, b, c, d) in P.items():
        cand = [ch for ch in pg.chars if ch['text'] == p and ch['size'] >= 9.5 and 'Bol' in ch['fontname']
                and a-1 <= ch['x0'] <= c and b-14 <= ch['top'] <= d]
        assert cand, (pdf, p)
        ch = min(cand, key=lambda ch: ch['top']+ch['x0'])
        out[p] = (float(ch['x0']), float(ch['top']), float(ch['x1']), float(ch['bottom']))
    return out


def inkboxes(pdf, P, L, dpi=300):
    tmp = Path(tempfile.mkdtemp())/'r'
    subprocess.run(['pdftoppm', '-r', str(dpi), '-png', '-singlefile', str(pdf), str(tmp)], check=True, capture_output=True)
    im = np.asarray(Image.open(str(tmp)+'.png').convert('L'), dtype=float) < 245; k = dpi/72
    out = {}; bands = {}; lbox = {}; conflict = {}
    # grow each letter's rectangle to its whole glyph (master letters render beyond their reported char box)
    for p in list(L):
        x0_, y0_, x1_, y1_ = L[p]
        for _ in range(30):
            grew = False
            for side in range(4):
                r = [x0_, y0_, x1_, y1_]
                ring = {0: (x0_-1.0, y0_, x0_, y1_), 1: (x0_, y0_-1.0, x1_, y0_), 2: (x1_, y0_, x1_+1.0, y1_), 3: (x0_, y1_, x1_, y1_+1.0)}[side]
                reg = im[int(ring[1]*k):int(np.ceil(ring[3]*k)), int(ring[0]*k):int(np.ceil(ring[2]*k))]
                if reg.any():
                    grew = True
                    if side == 0: x0_ -= 0.5
                    elif side == 1: y0_ -= 0.5
                    elif side == 2: x1_ += 0.5
                    else: y1_ += 0.5
            if not grew: break
        L[p] = (x0_, y0_, x1_, y1_)
    allm = np.zeros_like(im)
    for lx0, ly0, lx1, ly1 in L.values():
        allm[int((ly0-0.9)*k):int(np.ceil((ly1+0.9)*k)), int((lx0-0.9)*k):int(np.ceil((lx1+0.9)*k))] = True
    for p in P:
        lx0, ly0, lx1, ly1 = L[p]
        reg = im[int((ly0-0.9)*k):int(np.ceil((ly1+0.9)*k)), int((lx0-0.9)*k):int(np.ceil((lx1+0.9)*k))]
        ly, lx = np.where(reg)
        ox, oy = int((lx0-0.9)*k)/k, int((ly0-0.9)*k)/k
        lbox[p] = (ox+lx.min()/k, oy+ly.min()/k, ox+(lx.max()+1)/k, oy+(ly.max()+1)/k)
    im = im & ~allm
    for p, (a, b, c, d) in P.items():
        y0, y1, x0, x1 = max(0, int(round(b*k))), min(im.shape[0], int(round(d*k))), max(0, int(round(a*k))), min(im.shape[1], int(round(c*k)))
        sub = im[y0:y1, x0:x1].copy()
        ys, xs = np.where(sub)
        out[p] = tuple(float(v) for v in (x0/k+xs.min()/k, y0/k+ys.min()/k, x0/k+(xs.max()+1)/k, y0/k+(ys.max()+1)/k))
        # the letter's slot at the content box's top-left corner must be empty
        lw, lh = lbox[p][2]-lbox[p][0], lbox[p][3]-lbox[p][1]
        sx0, sy0 = max(0, int(round((lbox[p][0]-x0/k)*k))), int(round((out[p][1]-y0/k)*k))   # slot: letter's x, content's top
        sl = sub[sy0:sy0+int((lh+SLOT_PAD[1])*k), sx0:sx0+int((lw+SLOT_PAD[0])*k)]
        if sl.any():
            yy, xx = np.where(sl)
            conflict[p] = [round(x0/k+(sx0+xx.min())/k, 1), round(y0/k+(sy0+yy.min())/k, 1), round(x0/k+(sx0+xx.max()+1)/k, 1),
                           round(y0/k+(sy0+yy.max()+1)/k, 1)]
        crop = sub[int(round((out[p][1]-y0/k)*k)):int(round((out[p][3]-y0/k)*k)), int(round((out[p][0]-x0/k)*k)):int(round((out[p][2]-x0/k)*k))]
        blank = ~crop.any(axis=1); runs = []; i = 0
        while i < len(blank):
            if blank[i]:
                j = i
                while j < len(blank) and blank[j]: j += 1
                runs.append((out[p][1]+i/k, out[p][1]+j/k)); i = j
            else:
                i += 1
        bands[p] = runs
    return out, bands, lbox, conflict


def neighbours(P):
    left, up = {p: [] for p in P}, {p: [] for p in P}
    for p, (a, b, c, d) in P.items():
        for q, (a2, b2, c2, d2) in P.items():
            if p == q: continue
            if abs(c2-a) < 0.6 and min(d, d2)-max(b, b2) > 5: left[p].append(q)
            if abs(d2-b) < 0.6 and min(c, c2)-max(a, a2) > 5: up[p].append(q)
    return left, up


def check(fig, P=None, d=None):
    d = Path(d) if d else HERE.parent/f'figure{fig}'; draft = d/f'Figure.{fig}.draft.pdf'; P = P or PANELS[fig]
    L = letters(draft, P); box, bands, lbox, conflict = inkboxes(draft, P, L)
    return conflict


def retile(fig, P=None, clean=True, parts=None, assemble=True, d=None):
    d = Path(d) if d else HERE.parent/f'figure{fig}'
    draft = d/f'Figure.{fig}.draft.pdf'
    P = P or PANELS[fig]
    pg = pdfplumber.open(draft).pages[0]; W, H = float(pg.width), float(pg.height)
    L = letters(draft, P)
    box, bands, lbox, conflict = inkboxes(draft, P, L)
    hgt = {p: box[p][3]-box[p][1] for p in P}; ins = {p: [] for p in P}
    for sa, sb in EQUAL.get(fig, []):
        ha = sum(hgt[q] for q in sa)+G*(len(sa)-1); hb = sum(hgt[q] for q in sb)+G*(len(sb)-1)
        if abs(ha-hb) < 0.3: continue
        short, delta = (sa, hb-ha) if ha < hb else (sb, ha-hb)
        cand = sorted(((b1-b0, q, (b0+b1)/2) for q in short for b0, b1 in bands[q] if b1-b0 >= 3.0), reverse=True)[:NBANDS]
        assert cand, (fig, short)
        for w_, q, yc in cand:
            ins[q].append((yc, delta/len(cand)))
        for w_, q, yc in cand: hgt[q] += delta/len(cand)
    left, up = neighbours(P)
    # horizontally the letter stays where the build put it (inside the panel, at its left); vertically it is
    # aligned with the content top, and only the content counts for vertical alignment
    hx0 = {p: min(box[p][0], lbox[p][0]) for p in P}; hx1 = {p: max(box[p][2], lbox[p][2]) for p in P}
    nx, ny = {}, {}
    for p in sorted(P, key=lambda q: P[q][0]):
        nx[p] = M if not left[p] else max(nx[q]+(hx1[q]-hx0[q]) for q in left[p])+G
    for p in sorted(P, key=lambda q: P[q][1]):
        ny[p] = M if not up[p] else max(ny[q]+hgt[q] for q in up[p])+G
    newH = float(max(ny[p]+hgt[p] for p in P)+M)
    right = {p: round(W-(nx[p]+hx1[p]-hx0[p]), 2) for p in P}
    import io, hashlib
    cache = d/'_retile_cache'; cache.mkdir(parents=True, exist_ok=True)
    tag = hashlib.sha256(open(draft, 'rb').read()).hexdigest()[:16]
    todo = list(P) if parts is None else parts
    final = {}
    # one heavy clear per slot (everything outside the slot removed); segments and letters are cut from that
    reduced = {}
    for q, (qa, qb, qc, qd) in P.items():
        if q not in todo and not any(slot_of_center(lbox[t_], P) == q for t_ in todo): continue
        base = pikepdf.open(draft)
        clear(base, H, [(-20, -20, qa-0.05, H+20), (qc+0.05, -20, W+20, H+20), (-20, -20, W+20, qb-0.05), (-20, qd+0.05, W+20, H+20)])
        prune_xobjects(base); buf = io.BytesIO(); base.save(buf); reduced[q] = buf.getvalue()
    def slot_of(pt):
        return [q for q, (qa, qb, qc, qd) in P.items() if qa <= pt[0] <= qc and qb <= pt[1] <= qd][0]
    for p in P:
        final[p] = [round(v, 2) for v in (nx[p]-G/2, ny[p]-G/2, nx[p]+hx1[p]-hx0[p]+G/2, ny[p]+hgt[p]+G/2)]
    for p in todo:
        part = pikepdf.new(); part.add_blank_page(page_size=(W, newH)); page = part.pages[0]
        a, b, c, dd = P[p]; x0, y0, x1, y1 = box[p]
        cx0, cy0, cx1, cy1 = x0-0.8, y0-0.8, x1+0.8, y1+0.8
        cuts = sorted(ins[p]); edges = [cy0]+[yc for yc, _ in cuts]+[cy1]; offs = [0.0]
        for _, e in cuts: offs.append(offs[-1]+e)
        for (s0, s1), off in zip(zip(edges[:-1], edges[1:]), offs):
            src = pikepdf.open(io.BytesIO(reduced[p]))
            # only this segment's objects (nothing hidden by the clip is carried along), and no panel letters
            clear(src, H, [(-20, -20, W+20, s0-0.05), (-20, s1+0.05, W+20, H+20)]+[(q0-0.8, q1-0.8, q2+0.8, q3+0.8) for q0, q1, q2, q3 in L.values()])
            prune_xobjects(src)
            sp = src.pages[0]
            sp.MediaBox = [cx0, H-s1, cx1, H-s0]; sp.CropBox = [cx0, H-s1, cx1, H-s0]
            dx, dy = nx[p]-hx0[p], ny[p]-y0+off
            page.add_overlay(sp, pikepdf.Rectangle(cx0+dx, newH-(s1+dy), cx1+dx, newH-(s0+dy)))
        # the letter: own copy, ink top on the content top, x as drawn
        lx0, ly0, lx1, ly1 = L[p]; r0, r1, r2, r3 = lx0-0.8, ly0-0.8, lx1+0.8, ly1+0.8
        src = pikepdf.open(io.BytesIO(reduced[slot_of(((lbox[p][0]+lbox[p][2])/2, (lbox[p][1]+lbox[p][3])/2))]))
        clear(src, H, [(-20, -20, r0-0.05, H+20), (r2+0.05, -20, W+20, H+20), (-20, -20, W+20, r1-0.05), (-20, r3+0.05, W+20, H+20)])
        prune_xobjects(src)
        sp = src.pages[0]; sp.MediaBox = [r0, H-r3, r2, H-r1]; sp.CropBox = [r0, H-r3, r2, H-r1]
        dx, dy = nx[p]-hx0[p], ny[p]-lbox[p][1]
        page.add_overlay(sp, pikepdf.Rectangle(r0+dx, newH-(r3+dy), r2+dx, newH-(r1+dy)))
        prune_xobjects(part); part.save(cache/f'F{fig}_{p}_{tag}.pdf')
    if not assemble:
        return None
    out = pikepdf.new(); out.add_blank_page(page_size=(W, newH)); page = out.pages[0]
    keep_open = []
    for p in P:
        cp = cache/f'F{fig}_{p}_{tag}.pdf'
        assert cp.exists(), f'missing part {cp.name}: run retile.py {fig} --parts={p} first'
        keep_open.append(pikepdf.open(cp))
        page.add_overlay(keep_open[-1].pages[0], pikepdf.Rectangle(0, 0, W, newH))
    prune_xobjects(out)
    out.save(d/f'Figure.{fig}.pdf', object_stream_mode=pikepdf.ObjectStreamMode.generate, compress_streams=True, recompress_flate=True)
    if clean:
        from pdfclean import finalize
        clean = finalize(d/f'Figure.{fig}.pdf', label=f'Figure.{fig}.pdf', log=d/'sizes.tsv')
    else:
        clean = 'pending: run pdfclean.py on the final'
    assert not conflict, f'Figure {fig}: content inside the letter slot {conflict}' 
    json.dump(final, open(d/f'Figure.{fig}.panels.json', 'w'), indent=1)
    rep = dict(G=G, M=M, page=[W, H], new_height=round(newH, 2), shift={p: [round(nx[p]-hx0[p], 2), round(ny[p]-box[p][1], 2)] for p in P},
               draft_inkboxes={p: [round(v, 2) for v in box[p]] for p in P}, right_margin=right,
               inserted_space={p: [[round(y_, 1), round(e_, 2)] for y_, e_ in v] for p, v in ins.items() if v},
               final_heights={p: round(hgt[p], 2) for p in P}, clean=clean, letter_ink=lbox, letter_slot_conflicts=conflict)
    json.dump(rep, open(d/'retile.json', 'w'), indent=1)
    return rep


if __name__ == '__main__':
    args = [v for v in sys.argv[1:] if not v.startswith('--')]
    if '--check' in sys.argv:
        for f in [int(v) for v in args] or [1, 2, 3, 4, 5]:
            print(f, check(f))
        sys.exit(0)
    parts = None
    for v in sys.argv[1:]:
        if v.startswith('--parts='): parts = v.split('=', 1)[1].split(',')
        if v == '--assemble-only': parts = []
    for f in [int(v) for v in args] or [1, 2, 3, 4, 5]:
        r = retile(f, clean='--no-clean' not in sys.argv, parts=parts, assemble=parts is None or '--assemble' in sys.argv or '--assemble-only' in sys.argv)
        if r is None:
            print(f'Figure {f}: parts {parts} written'); continue
        bad = {p: v for p, v in r['right_margin'].items() if v < 8 and not (3.0 <= v <= 7.0)}
        print(f'Figure {f}: height {r["page"][1]:.1f} -> {r["new_height"]}, right margins {r["right_margin"]}', 'CHECK' if bad else '')
