"""Draw exact-coordinate vector annotations (PDF points, origin bottom-left) as a transparent page."""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from pathlib import Path
FD = Path(__file__).resolve().parent/'fonts'
FONTS = {k: FontProperties(fname=str(FD/f'NimbusSanL-{k}.otf')) for k in ['Reg', 'Bol', 'RegIta', 'BolIta']}
plt.rcParams.update({'pdf.fonttype': 42, 'svg.fonttype': 'none'})
INK = '#263238'


class Sheet:
    def __init__(self, w, h):
        self.w, self.h = w, h
        self.fig = plt.figure(figsize=(w/72, h/72))
        self.fig.patch.set_alpha(0)
        self.ax = self.fig.add_axes([0, 0, 1, 1]); self.ax.set_xlim(0, w); self.ax.set_ylim(0, h); self.ax.axis('off')
        self.ax.patch.set_alpha(0)

    def text(self, x, y, s, size=8, font='Reg', color=INK, ha='left', va='baseline', **kw):
        return self.ax.text(x, y, s, fontproperties=FONTS[font], fontsize=size, color=color, ha=ha, va=va, **kw)

    def line(self, xs, ys, lw=0.5, color=INK, **kw):
        self.ax.plot(xs, ys, lw=lw, color=color, solid_capstyle='butt', **kw)

    def rect(self, x, y, w, h, fc='white', ec='none', lw=0, **kw):
        from matplotlib.patches import Rectangle
        self.ax.add_patch(Rectangle((x, y), w, h, fc=fc, ec=ec, lw=lw, **kw))

    def save(self, path):
        self.fig.savefig(path, transparent=True); plt.close(self.fig)


_TT = {}
def textwidth(s, size, font='Reg'):
    from fontTools.ttLib import TTFont
    if font not in _TT:
        t = TTFont(str(FD/f'NimbusSanL-{font}.otf')); _TT[font] = (t.getBestCmap(), t['hmtx'], t['head'].unitsPerEm)
    cmap, hmtx, upm = _TT[font]
    return sum(hmtx[cmap.get(ord(ch), '.notdef')][0] for ch in s)*size/upm
