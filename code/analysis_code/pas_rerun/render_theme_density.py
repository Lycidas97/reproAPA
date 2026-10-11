"""Count rendered data encodings, separately from area and measured pixel coverage."""
import re,json
import numpy as np
from matplotlib.collections import PathCollection,LineCollection,PolyCollection
from matplotlib.container import BarContainer
from matplotlib.patches import Rectangle,PathPatch
from matplotlib import colors

def census(r):
    rows=[];details=[]
    for panel,axes,rect,keys in r.strict_groups:
        records=[];families=[];schematic=False
        for ax in sorted(axes,key=lambda a:r.fig.axes.index(a)):
            meta=next(m for m in r.axes_meta if m['axis'] is ax);family=meta['family'];families.append(family)
            def add(kind,n,channels=1,basis=''):
                if n:records.append(dict(axis=r.fig.axes.index(ax),mark_type=kind,marks=int(n),channels=int(channels),information_count=int(n)*int(channels),basis=basis))
            if family=='table':
                n=sum(bool(re.search(r'\d',t.get_text())) for t in ax.texts);add('table numeric cell',n,basis='Numeric body text; headers and categorical identities excluded');continue
            if family=='structure' and not ax.patches and not ax.collections:
                vals=sum(len(re.findall(r'(?<![A-Za-z])\d+(?:\.\d+)?',t.get_text())) for t in ax.texts)
                if vals:add('count schematic',vals,basis='Numeric count labels, excluding conceptual labels')
                else:schematic=True
                continue
            bars=set()
            for container in ax.containers:
                if isinstance(container,BarContainer):
                    for p in container.patches:
                        bars.add(p)
                        fc=np.array(p.get_facecolor())[:3]
                        if np.min(fc)>.85 and np.isclose(abs(p.get_width()),1):continue # Fixed 100% track.
                        add('bar',1,basis='One measured segment length, including an explicit zero')
            for im in ax.images:add('matrix cell',np.ma.count(np.ma.masked_invalid(im.get_array())),basis='Finite, unmasked color-encoded values; missing cells excluded')
            for c in ax.collections:
                if isinstance(c,PathCollection):
                    xy=np.ma.asarray(c.get_offsets());good=~np.ma.getmaskarray(xy).any(axis=1)&np.isfinite(np.ma.filled(xy,np.nan)).all(axis=1);n=int(good.sum())
                    channels=2 if family in ['scatter','embedding'] else 1
                    arr=c.get_array();sizes=c.get_sizes();fc=c.get_facecolors()
                    if arr is not None and len(arr)>1:channels+=1
                    elif len(fc)>1 and len(np.unique(fc,axis=0))>1:channels+=1
                    if len(sizes)>1 and len(np.unique(sizes))>1:channels+=1
                    add('point',n,channels,'Two coordinate channels for scatter/embedding; one for categorical-row points; variable color/area add channels')
                elif isinstance(c,LineCollection):
                    seg=c.get_segments();color=c.get_colors()
                    if len(color) and np.min(color[:,:3])>.80:continue # Row guides.
                    if c.get_transform()==ax.transData:add('interval/track segment',len(seg),2,'Two data endpoints; axes-fraction reference rules excluded')
                elif isinstance(c,PolyCollection):
                    n=sum(len(np.unique(p.vertices[:,0])) for p in c.get_paths())
                    add('density/band sample',n,2,'Finite x samples, upper and lower profile/band values; not extra biological replicates')
            for line in ax.lines:
                if line.get_transform()!=ax.transData:continue
                xy=np.asarray(line.get_xydata());n=np.isfinite(xy).all(axis=1).sum()
                if n and colors.to_rgba(line.get_color())[:3]==(0.8666666666666667,0.8901960784313725,0.9058823529411765):continue
                add('line sample',n,basis='Samples stored by the plotted line; axes-fraction references excluded')
            for t in ax.texts:
                if re.fullmatch(r'(?:n\s*=\s*)?[+−-]?\d+(?:\.\d+)?(?:M|%|/\d+)?(?:\s*\(\d+\))?',t.get_text().strip()):
                    add('numeric annotation',len(re.findall(r'\d+(?:\.\d+)?',t.get_text())),basis='Additional printed numeric encodings; not interpreted as independent observations')
            for p in ax.patches:
                if p in bars or p.get_zorder()<0:continue
                if isinstance(p,PathPatch):add('flow/track feature',1,basis='One encoded connection or structural feature')
                elif isinstance(p,Rectangle) and p.get_transform().contains_branch(ax.transData):
                    add('category/qualification glyph',1,basis='State strip or qualification outline; background row bands excluded')
            if family=='structure' and not records:schematic=True
        count=sum(a['information_count'] for a in records)
        if any(getattr(a,'_concept_schematic',False) for a in axes):count=None;schematic=True
        area=(rect[2]-rect[0])*(rect[3]-rect[1])/100
        rows.append(dict(panel=panel['panel'],mark_type=' + '.join(dict.fromkeys(a['mark_type'] for a in records)) or 'concept schematic',families=' + '.join(dict.fromkeys(families)),information_count=count if count is not None else None,panel_area_cm2=float(area),density_values_per_cm2=count/area if count is not None else None,rectangle_mm=panel['rectangle_mm'],schematic=schematic))
        details.append(dict(panel=panel['panel'],definition='Count visible data encodings, not independent observations or unique biological facts; repeated forms encode repeated values.',marks=records))
    (r.out/'verification'/f'{r.fid}.information.json').write_text(json.dumps(dict(panels=rows,artist_census=details),ensure_ascii=False,indent=2))
