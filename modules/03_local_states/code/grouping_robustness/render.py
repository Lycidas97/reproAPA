"""One compact evidence page; all detailed and adverse results remain in tables."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.transforms import Bbox

ROOT=Path(__file__).resolve().parents[1];D=ROOT/'data';OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
for name in ['Arial.ttf','Arial-Bold.ttf']:font_manager.fontManager.addfont(str(ROOT/'inputs'/name))
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.labelsize':8,'xtick.labelsize':7,
    'ytick.labelsize':7,'axes.linewidth':.45,'pdf.fonttype':42,'svg.fonttype':'path',
    'figure.facecolor':'white','savefig.facecolor':'white','axes.spines.top':False,'axes.spines.right':False})
W,H=180,252
fig=plt.figure(figsize=(W/25.4,H/25.4),dpi=150)
def axmm(x,y,w,h):return fig.add_axes([x/W,(H-y-h)/H,w/W,h/H])
def textmm(x,y,s,**kw):return fig.text(x/W,1-y/H,s,va='top',**kw)
def header(letter,y,title,note):
    return [textmm(8,y,letter,fontsize=10,fontweight='bold'),textmm(14,y+.5,title,fontsize=8,color='#263D4A'),
        textmm(14,y+6,note,fontsize=6.5,color='#52626D')]
def legend(y,items):
    arts=[]
    for x,c,label in items:
        arts.extend([textmm(x,y,'■',fontsize=7,color=c),textmm(x+3,y,label,fontsize=7)])
    return arts
textmm(8,6,'SMC4  State grouping and APA membership sensitivity',fontsize=9,color='#263D4A')
groups=[];plotted=[]
run=pd.read_csv(D/'partition_runs.tsv',sep='\t');run=run[run.role.eq('sensitivity')]
units=['gse118258','gse150311','gse178325','gse221739','wot_common_to_2i','gse103221_cip','gse103221_osk',
    'gse171127_d3_getm','gse171127_d3_oskm','gse171127_d6_pooled','gse171127_d12_pooled']
names=['GSE118258','GSE150311','GSE178325','GSE221739','WOT','GSE103221 CIP','GSE103221 OSK',
    'GSE171127 D3 GETM','GSE171127 D3 OSKM','GSE171127 D6','GSE171127 D12']
ha=header('A',18,'Partition stability and compatibility with original families',
    'Three resolutions × three seeds per original graph. Points: medians; lines: full ranges.')
a=axmm(45,31,119,49)
for j,unit in enumerate(units):
    q=run[run.unit.eq(unit)];assert len(q)==9
    for offset,metric,color,marker in [(-.15,'ARI','#0072B2','o'),(.15,'family_Bcubed_precision','#96753C','D')]:
        v=q[metric];a.plot([v.min(),v.max()],[j+offset]*2,c=color,lw=.85)
        a.scatter(v.median(),j+offset,c=color,marker=marker,s=13,zorder=3)
        plotted.append(dict(panel='A',unit=unit,metric=metric,minimum=v.min(),median=v.median(),maximum=v.max(),denominator=9))
a.set(xlim=(.3,1.02),ylim=(10.6,-.6));a.set_xticks([.4,.6,.8,1]);a.set_yticks(range(11),names)
a.tick_params(axis='y',length=0,pad=4);a.tick_params(axis='x',length=2,pad=2);a.spines['left'].set_visible(False)
la=legend(88,[(45,'#0072B2','Original-state ARI'),(92,'#96753C','Original-family compatibility')])
notea=textmm(14,95,'Compatibility measures original-family purity within alternative clusters; it is not annotation accuracy.',fontsize=6.5,color='#52626D')
groups.append(ha+[a,*la,notea])

effect=pd.read_csv(D/'local_effect_summary.tsv',sep='\t')
cases=[('WOT','exclude_WOT_L00','SR','WOT −L00: S–R'),('WOT','exclude_WOT_L00','RE','WOT −L00: R–E'),
    ('CIP','exclude_CIP_C012','RE','CIP −C012: R–E'),('CIP','exclude_CIP_C012','EP','CIP −C012: E–P'),
    ('GSE118258','exclude_human_C005','RE','Human −C005: R–E'),('GSE118258','exclude_human_C005','EP','Human −C005: E–P')]
hb=header('B',106,'Local APA effects after removing the indicated state',
    'Denominator: original effects with |ΔPA| ≥0.01. Purity and coverage thresholds reapplied.')
b=axmm(47,120,117,35)
cols=['same_direction','weakened','opposite_direction','new_unavailable'];colors=['#0072B2','#D4A64B','#A16A89','#D6DBDF']
for i,(route,variant,contrast,label) in enumerate(cases):
    q=effect[(effect.route==route)&(effect.variant==variant)&(effect.contrast==contrast)&(effect.policy=='reapply_purity')].iloc[0]
    n=int(q.baseline_effects_ge001);assert sum(int(q[c]) for c in cols)==n
    left=0
    for col,color in zip(cols,colors):
        value=int(q[col]);width=100*value/n
        b.barh(i,width,left=left,height=.66,color=color,edgecolor='white',lw=.3);left+=width
        plotted.append(dict(panel='B',unit=label,metric=col,count=value,denominator=n))
    b.text(102,i,str(n),va='center',fontsize=7)
b.set(xlim=(0,100),ylim=(5.65,-.65));b.set_xticks([0,50,100],['0','50','100%']);b.set_yticks(range(6),[x[-1] for x in cases])
b.text(1.015,1.045,'n',transform=b.transAxes,fontsize=7)
b.tick_params(axis='y',length=0,pad=4);b.tick_params(axis='x',length=2,pad=2);b.spines['left'].set_visible(False)
lb=legend(164,[(28,colors[0],'Same direction'),(72,colors[1],'Weakened'),(109,colors[2],'Opposite'),(141,colors[3],'Unavailable')])
groups.append(hb+[b,*lb])

shape=pd.read_csv(D/'shape_retention.tsv',sep='\t')
cases=[('mouse','exclude_WOT_L00','Mouse / −L00'),('mouse','exclude_CIP_C012','Mouse / −C012'),
    ('human','exclude_human_C005','Human / −C005'),('mouse','exclude_three','Mouse / −three states'),
    ('human','exclude_three','Human / −three states'),('human','exclude_human_candidates','Human / −candidate R/E')]
hc=header('C',181,'Retention of the original APA shape assignment',
    'Original four centers and full gene denominators fixed; refitted assignments are supplied in the table.')
c=axmm(51,195,113,35)
for i,(sp,variant,label) in enumerate(cases):
    q=shape[(shape.species==sp)&(shape.variant==variant)&(shape.policy=='reapply_purity')].iloc[0]
    n=int(q.fixed_genes);values=[int(q.same_fixed_center),int(q.evaluable-q.same_fixed_center),int(q.unavailable)]
    assert sum(values)==n
    left=0
    for value,color,metric in zip(values,[colors[0],colors[1],colors[3]],['same_original_shape','changed_shape','unavailable']):
        width=100*value/n;c.barh(i,width,left=left,height=.66,color=color,edgecolor='white',lw=.3);left+=width
        plotted.append(dict(panel='C',unit=label,metric=metric,count=value,denominator=n))
    c.text(102,i,str(n),va='center',fontsize=7)
c.set(xlim=(0,100),ylim=(5.65,-.65));c.set_xticks([0,50,100],['0','50','100%']);c.set_yticks(range(6),[x[-1] for x in cases])
c.text(1.015,1.045,'n',transform=c.transAxes,fontsize=7)
c.tick_params(axis='y',length=0,pad=4);c.tick_params(axis='x',length=2,pad=2);c.spines['left'].set_visible(False)
lc=legend(240,[(45,colors[0],'Same shape'),(89,colors[1],'Changed'),(125,colors[3],'Unavailable')])
groups.append(hc+[c,*lc])
fig.canvas.draw();renderer=fig.canvas.get_renderer();pw,ph=fig.canvas.get_width_height();scale=25.4/fig.dpi
boxes=[a.get_tightbbox(renderer) for a in fig.axes]+[t.get_window_extent(renderer) for t in fig.texts]
box=Bbox.union(boxes);rowboxes=[Bbox.union([a.get_tightbbox(renderer) for a in group]) for group in groups]
gaps=[(a.y0-b.y1)*scale for a,b in zip(rowboxes,rowboxes[1:])]
assert min(gaps)>1 and box.x0>=0 and box.y0>=0 and box.x1<=pw and box.y1<=ph,(gaps,box)
for ext in ['pdf','svg','png']:fig.savefig(OUT/f'SMC4_grouping_robustness.{ext}',dpi=300)
pd.DataFrame(plotted).to_csv(OUT/'plotted_values.tsv',sep='\t',index=False)
(OUT/'render_qa.json').write_text(json.dumps(dict(width_mm=W,height_mm=H,font='Arial; Nimbus Sans L unavailable',
    minimum_text_pt=6.5,all_rendered_bounds_inside=True,panel_gutters_mm=gaps,
    visible_margins_mm=dict(left=box.x0*scale,right=(pw-box.x1)*scale,top=(ph-box.y1)*scale,bottom=box.y0*scale)),indent=2)+'\n')
print('FIGURE_COMPLETE',gaps,flush=True)
