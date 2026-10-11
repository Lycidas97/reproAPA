"""Dataset swatches and full-name keys for the existing supplement renderer."""
import re
import numpy as np
from matplotlib.text import Text
from matplotlib.patches import Rectangle
from matplotlib.transforms import Affine2D

# Names follow the manuscript inventory and its route/series identities.
DATASETS={
 'Xing2020':('#753399',['Xing','H1','GSE118258']),
 'Liu2020':('#648290',['Liu','H2','GSE150311']),
 'Guan2022':('#8B1B56',['Guan','H3','GSE178325']),
 'Guan2022–0330':('#CB6489',['178325_0330']),
 'Guan2022–0605':('#662249',['Guan0605']),
 'Guan2022–0618':('#A82C60',['Guan0618','Guan–0618','178325_0618']),
 'Guan2022–0809':('#E0829F',['Guan0809','Guan–0809','178325_0809']),
 'Panariello2023':('#30363D',['Panariello','H4','GSE221739']),
 'Guo2019':('#B26E00',['Guo','GSE103221']),
 'Guo2019–CIP':('#CF4C00',['CIP','C']),
 'Guo2019–OSK':('#255FBC',['OSK']),
 'Schiebinger2019':('#4A8F72',['Schiebinger','GSE115943']),
 'Schiebinger2019–2i':('#207450',['WOT','W','WOT_GSE115943']),
 'Jaber2022':('#999933',['Jaber','GSE171127']),
 'Jaber2022–D3 GETM':('#17A8B8',['D3 GETM']),
 'Jaber2022–D3 OSKM':('#A5A51D',['D3 OSKM']),
 'Jaber2022–D6 pooled':('#8C564B',['D6']),
 'Jaber2022–D12 pooled':('#8C6D31',['D12']),
}
ALIASES={alias:key for key,(_,aliases) in DATASETS.items() for alias in aliases+[key]}
PATTERN=re.compile(r'(?<![A-Za-z0-9_])('+ '|'.join(re.escape(k) for k in sorted(ALIASES,key=len,reverse=True))+r')(?![A-Za-z0-9_])')
SLOT='\u2002\u2002'  # Two en spaces keep a measurable, fixed-size square slot.
def matches(s):return list(PATTERN.finditer(str(s)))
def dataset_color(s,default=None):
 m=matches(s)
 return DATASETS[ALIASES[m[0].group()]][0] if len(m)==1 else default

class DatasetText(Text):
    """Regular Text geometry, with vector squares occupying reserved text slots."""
    def set_text(self,s):
        # Axis formatters refresh tick text on every draw; retain the swatch replacement.
        if hasattr(self,'dataset_slots') and matches(s):
            rendered,slots=swatch_string(str(s));self.dataset_original=str(s);self.dataset_slots=slots
            return super().set_text(rendered)
        return super().set_text(s)
    def draw(self,renderer):
        super().draw(renderer)
        if not self.get_visible():return
        _,lines,_=self._get_layout(renderer)
        pos=self.get_transform().transform((float(self.convert_xunits(self._x)),float(self.convert_yunits(self._y))))
        if not np.isfinite(pos).all():return
        rotation=self.get_rotation();unit=renderer.points_to_pixels(1);side=self.get_fontsize()*.68*unit
        offset=0
        for line,_,(lx,ly) in lines:
            for start,key in self.dataset_slots:
                if offset<=start<offset+len(line):
                    width=renderer.get_text_width_height_descent(line[:start-offset],self._fontproperties,False)[0]
                    slot=renderer.get_text_width_height_descent(SLOT,self._fontproperties,False)[0]
                    tr=Affine2D().rotate_deg(rotation).translate(pos[0]+lx,pos[1]+ly)
                    mark=Rectangle((width+(slot-side)/2,.35*unit),side,side,transform=tr,facecolor=DATASETS[key][0],edgecolor='none',clip_on=False)
                    mark.set_figure(self.figure);mark.draw(renderer)
            offset+=len(line)+1

def swatch_string(original):
    ms=matches(original)
    pieces=[];slots=[];end=0
    for m in ms:
        pieces.append(original[end:m.start()]);start=sum(map(len,pieces));key=ALIASES[m.group()];slots.append((start,key));pieces.append(SLOT);end=m.end()
    pieces.append(original[end:]);return ''.join(pieces),slots

def convert_text(t):
    original=t.get_text()
    if not matches(original):return []
    if type(t) is not Text:raise TypeError(('Dataset text needs explicit handling',type(t),original))
    rendered,slots=swatch_string(original);t.set_text(rendered);t.__class__=DatasetText;t.dataset_slots=slots;t.dataset_original=original
    return [k for _,k in slots]

def prepare_dataset_labels(r):
    """Collect dataset identities by owning axis before replacing their text."""
    r.fig.canvas.draw();used=set();inventory=[];r.dataset_axis_names={}
    owned_text=set()
    for ax in r.fig.axes:
        local=set();lg=ax.get_legend()
        if lg is not None:
            labels=[t.get_text() for t in lg.get_texts()]
            for handle,label in zip(lg.legend_handles,labels):
                col=dataset_color(label)
                if col and hasattr(handle,'set_color'):handle.set_color(col)
            for artist in list(ax.lines)+list(ax.collections):
                local.update(ALIASES[m.group()] for m in matches(artist.get_label()))
                col=dataset_color(artist.get_label())
                if col and hasattr(artist,'set_color'):artist.set_color(col)
            # Dataset-only labels are reproduced inside this axis's final panel.
            if labels and all(len(matches(s))==1 and PATTERN.sub('',s).strip()=='' for s in labels):
                local.update(ALIASES[matches(s)[0].group()] for s in labels);lg.remove()
        for t in ax.findobj(match=lambda x:isinstance(x,Text)):
            if t in owned_text or not t.get_visible():continue
            owned_text.add(t);old=t.get_text();keys=convert_text(t)
            if keys:local.update(keys);inventory.append(dict(original=old,datasets=keys,axis_index=r.fig.axes.index(ax)))
        r.dataset_axis_names[ax]=[k for k in DATASETS if k in local];used.update(local)
    # A figure-level dataset identity would recreate an ambiguous shared key.
    for t in r.fig.texts:
        if t.get_visible() and matches(t.get_text()):raise ValueError(('Dataset text has no panel owner',t.get_text()))
    r.dataset_names=[k for k in DATASETS if k in used];r.dataset_label_inventory=inventory
