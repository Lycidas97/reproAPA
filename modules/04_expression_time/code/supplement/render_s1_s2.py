"""Render native marker evidence and actual-edge integration support; no model fitting."""
from pathlib import Path
import json,subprocess
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.transforms import Bbox
from PIL import Image
R=Path(__file__).resolve().parent;O=R/'figures';O.mkdir(exist_ok=True)
P=json.loads((R/'palette.json').read_text());CFG=json.loads((R/'marker_selection.json').read_text());COL=P['state_groups'];INK='#263D4A';BLUE='#377EB8';MUTED='#70777C'
CM=LinearSegmentedColormap.from_list('expression',['#EAF1F6','#89B7D5',BLUE,'#173F65'])
plt.rcParams.update({'font.family':'Nimbus Sans L','font.size':8.05,'axes.labelsize':8.05,'axes.titlesize':8,'xtick.labelsize':8.05,'ytick.labelsize':8.05,'axes.linewidth':.45,'axes.edgecolor':INK,'text.color':INK,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'path'})
NAMES={'GSE118258':'GSE118258','GSE150311':'GSE150311','GSE178325':'GSE178325','GSE221739':'GSE221739','gse103221_osk':'GSE103221 / OSK','gse103221_cip':'GSE103221 / CIP','gse171127_d3_getm':'GSE171127 / D3 GETM','gse171127_d3_oskm':'GSE171127 / D3 OSKM','gse171127_d6_pooled':'GSE171127 / D6','gse171127_d12_pooled':'GSE171127 / D12','wot_common_to_2i':'GSE115943 / 2i'}
FAMILY={'shared_s':'S: source / fibroblast','shared_r':'R: remodeling candidate','shared_e':'E: epithelial / mixed','shared_tp':'P: pluripotency-associated','shared_m':'M: ECM-retaining','shared_k':'K: keratin-associated','human_x':'hX: endoderm-associated','human_n':'hN: naive-associated','human_u':'hU: unresolved','mouse_b':'mB: trophoblast-associated','mouse_n':'mN: neural / neurokeratin','mouse_i':'mI: interferon','mouse_y':'mY: myeloid-associated','mouse_u':'mU: unresolved'}
CODE={k:v.split(':')[0] for k,v in FAMILY.items()}

STUDY={'GSE118258':'Xing2020','GSE150311':'Liu2020','GSE178325':'Guan2022','GSE221739':'Panariello2023','GSE103221':'Guo2019','GSE171127':'Jaber2022','GSE115943':'Schiebinger2019'}
def display_name(label):
 for accession,study in STUDY.items():label=label.replace(accession,study)
 return label.removeprefix('WOT_')

STUDY_COL={'Xing2020':'#6F2C91','Liu2020':'#607D8B','Guan2022':'#A66B00','Panariello2023':'#30363D','Guo2019':'#255FBC','Jaber2022':'#78688F','Schiebinger2019':'#146B4A'}
def study_of(unit):
 if unit in STUDY:return STUDY[unit]
 for prefix,study in [('GSE178325','Guan2022'),('gse103221','Guo2019'),('gse171127','Jaber2022'),('wot_','Schiebinger2019')]:
  if unit.startswith(prefix):return study
 raise ValueError(unit)
def keyed_text(page,x,y,text,study,**kw):
 t=page.text(x,y,text,**kw);page.fig.canvas.draw()
 bb=t.get_window_extent(page.fig.canvas.get_renderer());left=bb.x0*72/page.fig.dpi
 page.text(left-8,y,'■',color=STUDY_COL[study],g=kw.get('g'))
 return t
def key_axis_labels(page,ax,studies):
 page.fig.canvas.draw();renderer=page.fig.canvas.get_renderer()
 for label,study in zip(ax.get_yticklabels(),studies):
  bb=label.get_window_extent(renderer);page.text(bb.x0*72/page.fig.dpi-8,page.h-(bb.y0+bb.y1)/2*72/page.fig.dpi,'■',color=STUDY_COL[study],va='center')

def read(n):return pd.read_csv(R/'source_data'/n,sep='\t',low_memory=False)
class Page:
 def __init__(self,h):
  self.w=160/25.4*72;self.h=h;self.fig=plt.figure(figsize=(self.w/72,h/72));self.gutters=[];self.groups={}
 def text(self,x,y,s,g=None,**kw):
  s=display_name(s);kw.setdefault('va','top');t=self.fig.text(x/self.w,1-y/self.h,s,**kw)
  if g:self.groups.setdefault(g,[]).append(t)
  return t
 def ax(self,x,y,w,h,g=None,**kw):
  a=self.fig.add_axes([x/self.w,1-(y+h)/self.h,w/self.w,h/self.h],**kw)
  a.spines[['top','right']].set_visible(False);a.tick_params(length=2,width=.45,pad=2)
  if g:self.groups.setdefault(g,[]).append(a)
  return a
 def save(self,name):
  self.fig.canvas.draw();r=self.fig.canvas.get_renderer();factor=72/self.fig.dpi
  # Visible union includes glyphs and ticks; center horizontally without resizing.
  bb=Bbox.union([a.get_tightbbox(r) for a in self.fig.axes]+[t.get_window_extent(r) for t in self.fig.texts])
  dx=(self.fig.bbox.width-bb.x1-bb.x0)/2/self.fig.bbox.width
  for a in self.fig.axes:
   b=a.get_position();a.set_position([b.x0+dx,b.y0,b.width,b.height])
  for t in self.fig.texts:x,y=t.get_position();t.set_position((x+dx,y))
  self.fig.canvas.draw();r=self.fig.canvas.get_renderer()
  boxes={k:Bbox.union([a.get_tightbbox(r) if hasattr(a,'get_tightbbox') else a.get_window_extent(r) for a in v]) for k,v in self.groups.items()}
  gaps={a+' / '+b:float((boxes[b].x0-boxes[a].x1 if orient=='h' else boxes[a].y0-boxes[b].y1)*factor) for a,b,orient in self.gutters}
  assert all(v>=2 for v in gaps.values()),(name,gaps)
  allbb=Bbox.union([a.get_tightbbox(r) for a in self.fig.axes]+[t.get_window_extent(r) for t in self.fig.texts]);assert allbb.x0>0 and allbb.x1<self.fig.bbox.width and allbb.y0>0 and allbb.y1<self.fig.bbox.height,(name,allbb)
  raw=O/(name+'_raw.pdf');pdf=O/(name+'.pdf')
  self.fig.savefig(raw,dpi=300);self.fig.savefig(O/(name+'.svg'),dpi=300)
  subprocess.run(['pdftocairo','-pdf',str(raw),str(pdf)],check=True,capture_output=True);raw.unlink()
  # Cairo fixes CFF embedding declarations; restore the exact artwork page size.
  from pypdf import PdfReader,PdfWriter
  reader=PdfReader(pdf);writer=PdfWriter();pg=reader.pages[0];pg.mediabox.upper_right=(self.w,self.h);pg.cropbox.upper_right=(self.w,self.h);writer.add_page(pg)
  with pdf.open('wb') as handle:writer.write(handle)
  subprocess.run(['pdftoppm','-r','300','-singlefile','-png',str(O/(name+'.pdf')),str(O/name)],check=True,capture_output=True)
  im=Image.open(O/(name+'.png'));im.convert('L').save(O/(name+'_gray.png'))
  pixels=np.array(im.convert('RGB'));ys,xs=np.where(pixels.min(2)<245)
  qa={'width_mm':160,'height_mm':self.h/72*25.4,'left_margin_mm':float(xs.min()/300*25.4),'right_margin_mm':float((pixels.shape[1]-1-xs.max())/300*25.4),'top_margin_mm':float(ys.min()/300*25.4),'bottom_margin_mm':float((pixels.shape[0]-1-ys.max())/300*25.4),'measured_gutters_pt':gaps,'font':'Nimbus Sans L','minimum_font_pt':8.05}
  assert abs(qa['left_margin_mm']-qa['right_margin_mm'])<.8
  plt.close(self.fig);return qa

UNIT_COL={**P['routes'],'GSE150311':'#607D8B','GSE178325':'#A66B00',
 'gse171127_d3_getm':'#53406D','gse171127_d3_oskm':'#78688F',
 'gse171127_d6_pooled':'#A396B7','gse171127_d12_pooled':'#C7BED3'}
UNIT_LABEL={'GSE118258':'Xing2020','GSE150311':'Liu2020','GSE178325':'Guan2022','GSE221739':'Panariello2023',
 'gse103221_osk':'Guo2019 · OSK','gse103221_cip':'Guo2019 · CIP',
 'gse171127_d3_getm':'Jaber2022 · D3 GETM','gse171127_d3_oskm':'Jaber2022 · D3 OSKM',
 'gse171127_d6_pooled':'Jaber2022 · D6','gse171127_d12_pooled':'Jaber2022 · D12',
 'wot_common_to_2i':'Schiebinger2019 · 2i',
 'GSE178325::donor_0605':'Guan2022 · Adult skin fibroblasts',
 'GSE178325::donor_0618':'Guan2022 · Adipose stromal cells\nTime course'}
def colour_key(page,x,y,label,colour,g='unit_key'):
 page.text(x,y,'■',color=colour,g=g)
 page.text(x+10,y,label,g=g,linespacing=1.12)
def marker_figure(sp,name):
 s=read('states.tsv');s=s[s.species.eq(sp)].copy().fillna('');g=read('marker_by_state.tsv');g=g[g.species.eq(sp)];mods=read('program_by_state_library.tsv.gz');mods=mods[mods.species.eq(sp)&mods.n_cells.ge(20)]
 genes=[v for vals in CFG['markers'][sp].values() for v in vals];n=len(s);page=Page(566 if sp=='mouse' else 546)
 x=62;w=380;top=28;h=252;bottom=top+h
 page.text(12,8,'A',fontsize=10,weight='bold')
 a=page.ax(x,top,w,h,'A');a.set_xlim(-.5,n-.5);a.set_ylim(len(genes)-.5,-.5);a.set_yticks(range(len(genes)),genes,fontstyle='italic');a.set_xticks([]);a.tick_params(length=0,pad=3);a.spines[:].set_visible(False)
 ix=dict(zip(s.state_key,range(n)));iy=dict(zip(genes,range(len(genes))))
 a.scatter(g.state_key.map(ix),g.gene.map(iy),s=17*g.detected_fraction,c=g.relative_expression,cmap=CM,vmin=0,vmax=1,edgecolors='none',zorder=2)
 f=g.fallback_small;a.scatter(g.loc[f,'state_key'].map(ix),g.loc[f,'gene'].map(iy),s=22,facecolors='none',edgecolors='#6E7377',linewidths=.35)
 # Two aligned annotation tracks replace detached source/route headings.
 for yy,values,palette,label in [(8,s.unit,UNIT_COL,'Unit'),(18,s.family,COL,'Family')]:
  bar=page.ax(x,yy,w,6,'Atracks');bar.set_xlim(-.5,n-.5);bar.set_ylim(0,1);bar.axis('off')
  for i,value in enumerate(values):bar.add_patch(Rectangle((i-.5,0),1,1,color=palette[value],lw=0))
  page.text(x-5,yy-1,label,ha='right',g='Atracks')
 for i,row in enumerate(s.itertuples()):
  label=str(int(row.short_id[1:])) if row.short_id.startswith('C') else row.short_id
  page.text(x+(i+.5)/n*w,bottom+4+(i%2)*11,label,ha='center',g='Aflags')
  if row.review_flag:page.text(x+(i+.5)/n*w,307+(i%2)*10,row.review_flag,ha='center',g='Aflags')
 starts=0
 for unit,q in s.groupby('unit',sort=False):
  if starts:a.axvline(starts-.5,color='#C7CDD1',lw=.5)
  starts+=len(q)
 cb=page.ax(67,330,80,5,'key');plt.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,1),cmap=CM),cax=cb,orientation='horizontal',ticks=[0,.5,1]);cb.set_xlabel('Gene-relative expression',labelpad=2)
 da=page.ax(260,326,145,20,'sizekey');da.set_xlim(0,145);da.set_ylim(0,20);da.axis('off')
 for xx,v in zip([9,49,89],[.25,.5,1]):da.scatter(xx,12,s=17*v,color=BLUE);da.text(xx+6,12,f'{int(v*100)}%',va='center')
 da.text(9,0,'Detected cells (library mean)',va='top')
 families=[k for k in FAMILY if k in set(s.family)]
 for i,k in enumerate(families):
  xx=12+(i%3)*145;yy=362+(i//3)*11
  page.text(xx,yy,'●',color=COL[k],g='A_families');page.text(xx+9,yy,FAMILY[k],g='A_families')
 btop=408
 page.text(12,btop,'B',fontsize=10,weight='bold');ba=page.ax(x,btop,w,27,'B');ba.set_xlim(-.5,n-.5);ba.set_ylim(2.5,-.5)
 libmax=20 if sp=='human' else 80
 vals=np.stack([np.log10(s.n_cells.astype(float))/5,s.n_libraries_ge20.astype(float)/libmax,s.max_library_fraction.astype(float)])
 assert (vals>=0).all() and (vals<=1).all()
 ba.imshow(vals,cmap='Greys',vmin=0,vmax=1,aspect='auto',interpolation='nearest');ba.set_yticks([0,1,2],['log10 n','Lib. ≥20','Max lib.']);ba.set_xticks([]);ba.tick_params(length=0,pad=3);ba.spines[:].set_visible(False)
 page.text(x,btop+31,f'Ranges: 0–5; 0–{libmax}; 0–100% (dark = high).',g='Brange')
 ctop=464
 page.text(12,452,'C',fontsize=10,weight='bold');units=s.unit.drop_duplicates().tolist();maxy=np.ceil(mods.mean_log1p_cp10k.max()*2)/2
 for j,program in enumerate(['ECM','Epithelial','Pluripotency']):
  ax=page.ax(43+j*139,ctop,117,35,'C'+str(j));q=mods[mods.program.eq(program)].sort_values(['state_key','library_key']);rng=np.random.default_rng(1700+j)
  xx=q.unit.map(dict(zip(units,range(len(units))))).to_numpy()+rng.uniform(-.23,.23,len(q))
  ax.scatter(xx,q.mean_log1p_cp10k,s=4.5,c=q.family.map(COL),alpha=.65,linewidths=0,rasterized=True)
  ax.set_xlim(-.5,len(units)-.5);ax.set_ylim(0,maxy*1.04);ax.set_xticks([]);ax.set_ylabel(program,labelpad=3)
  if j:ax.tick_params(labelleft=False)
  bar=page.ax(43+j*139,502,117,5,'Cstrip');bar.set_xlim(-.5,len(units)-.5);bar.set_ylim(0,1);bar.axis('off')
  for i,unit in enumerate(units):bar.add_patch(Rectangle((i-.44,0),.88,1,color=UNIT_COL[unit],lw=0))
 # Full author/year labels live in a dedicated key, never on cramped plot ticks.
 for i,unit in enumerate(units):
  col=i//4 if sp=='mouse' else i//2;row=i%4 if sp=='mouse' else i%2
  colour_key(page,62+col*199,515+row*12,UNIT_LABEL[unit],UNIT_COL[unit])
 page.gutters=[('Aflags','key','v'),('key','A_families','v'),('sizekey','A_families','v'),('A_families','B','v'),('Brange','C0','v'),('C0','C1','h'),('C1','C2','h'),('Cstrip','unit_key','v')]
 return page.save(name)

ROUTES=['GSE118258','GSE178325::donor_0605','GSE178325::donor_0618','GSE221739','gse103221_osk','gse103221_cip','wot_common_to_2i']
RL=['GSE118258','GSE178325\nAdult skin fibroblasts','GSE178325\nAdipose stromal cells\nTime course','GSE221739','GSE103221 / OSK','GSE103221 / CIP','GSE115943 / 2i']
def scatter_map(ax,d,c,**kw):
 ax.scatter(d.UMAP1,d.UMAP2,s=.7,c=c,linewidths=0,rasterized=True,**kw);ax.set_xlim(-5.3,15.1);ax.set_ylim(-4.4,15.2);ax.set_aspect('equal');ax.axis('off')
def integration_figure():
 page=Page(566);s=read('states.tsv');cells=read('display_cells.tsv.gz');expr=read('display_expression.tsv.gz')
 assert set(cells.cell_id)==set(expr.cell_id);cells=cells.merge(expr.drop(columns='species_key'),on='cell_id',validate='one_to_one')
 page.text(13,8,'A',fontsize=10,weight='bold');page.text(258,8,'B',fontsize=10,weight='bold')
 for i,sp in enumerate(['human','mouse']):
  ax=page.ax(14+i*113,20,105,90,'A'+str(i));q=cells[cells.species_key.eq(sp)].sample(frac=1,random_state=27);scatter_map(ax,q,q.reviewed_group.map(COL));page.text(66+i*113,112,sp.capitalize(),ha='center')
 # The categorical legend belongs to A; B retains its separate quantitative scale.
 short_family={'shared_s':'S: source / fibroblast','shared_r':'R: remodeling','shared_e':'E: epithelial / mixed','shared_tp':'P: pluripotency-associated','shared_m':'M: ECM-retaining','shared_k':'K: keratin-associated','human_x':'hX: endoderm-associated','human_n':'hN: naive-associated','human_u':'hU: unresolved','mouse_b':'mB: trophoblast-associated','mouse_n':'mN: neural / neurokeratin','mouse_i':'mI: interferon','mouse_y':'mY: myeloid-associated','mouse_u':'mU: unresolved'}
 for i,(family,label) in enumerate(short_family.items()):
  xx=12+(i//7)*127;yy=134+(i%7)*11
  page.text(xx,yy,'●',color=COL[family],g='Akey');page.text(xx+8,yy,label,fontsize=8.05,g='Akey')
 edges=read('actual_state_pair_edges.tsv.gz');hs=s[s.species.eq('human')];ms=s[s.species.eq('mouse')]
 mat=edges.pivot(index='human_state',columns='mouse_state',values='human_row_share').reindex(index=hs.state_key,columns=ms.state_key)
 bx,by,bw,bh=327,23,110,114
 ax=page.ax(bx,by,bw,bh,'B');masked=np.ma.masked_invalid(mat.to_numpy());im=ax.imshow(masked,cmap=CM,vmin=0,vmax=1,aspect='auto',interpolation='nearest')
 # Study labels and adjacent family strips identify the actual original-state pairs.
 def blocks(frame):
  out=[];start=0
  for study,part in frame.groupby('study_id',sort=False):
   out.append((start+(len(part)-1)/2,display_name(study)));start+=len(part)
  return out
 hb=blocks(hs);mb=blocks(ms)
 ax.set_yticks([a for a,b in hb],[b for a,b in hb]);ax.tick_params(axis='y',length=0,pad=10)
 ax.set_xticks([a for a,b in mb],[b[:-4]+'\n'+b[-4:] for a,b in mb],rotation=40,ha='right');ax.tick_params(axis='x',length=0,pad=3)
 for positions,direction in [(hs.groupby('unit',sort=False).size().cumsum().to_numpy()[:-1], 'h'),(ms.groupby('unit',sort=False).size().cumsum().to_numpy()[:-1], 'v')]:
  for n in positions:(ax.axhline(n-.5,color='#8A939A',lw=.4) if direction=='h' else ax.axvline(n-.5,color='#8A939A',lw=.4))
 for frame,orient in [(hs,'y'),(ms,'x')]:
  bar=page.ax(bx-6,by,4,bh,'Bstrip') if orient=='y' else page.ax(bx,by-6,bw,4,'Bstrip')
  bar.set_xlim(-.5,len(frame)-.5) if orient=='x' else bar.set_xlim(0,1)
  bar.set_ylim(0,1) if orient=='x' else bar.set_ylim(len(frame)-.5,-.5)
  bar.axis('off')
  for i,family in enumerate(frame.family):bar.add_patch(Rectangle((i-.5,0) if orient=='x' else (0,i-.5),1,1,color=COL[family],lw=0))
 for frame,orient in [(hs,'y'),(ms,'x')]:
  bar=page.ax(bx-12,by,4,bh) if orient=='y' else page.ax(bx,by-12,bw,4)
  bar.set_xlim(0,1) if orient=='y' else bar.set_xlim(-.5,len(frame)-.5)
  bar.set_ylim(len(frame)-.5,-.5) if orient=='y' else bar.set_ylim(0,1);bar.axis('off')
  for i,unit in enumerate(frame.unit):bar.add_patch(Rectangle((0,i-.5) if orient=='y' else (i-.5,0),1,1,color=STUDY_COL[study_of(unit)],lw=0))
 ax.tick_params(axis='y',pad=17)
 state_index=s[['species','figure_row','study_id','unit','short_id','state_key','family']].copy()
 state_index['study_label']=state_index.study_id.map(display_name)
 state_index.to_csv(R/'state_index.tsv',sep='\t',index=False)
 # Mark original states without any outgoing edge as NE rather than zero.
 for idx in np.flatnonzero(mat.isna().all(axis=1)):ax.axhspan(idx-.5,idx+.5,color='#D9D9D9',lw=0)
 cb=page.ax(305,186,122,5,'Bkey');plt.colorbar(im,cax=cb,orientation='horizontal',ticks=[0,.5,1]);cb.set_xlabel('Row-normalized edge weight',labelpad=2)
 page.text(12,215,'C',fontsize=10,weight='bold')
 scales=[]
 for j,program in enumerate(['ECM','Epithelial','Pluripotency']):
  ax=page.ax(12+j*148,218,133,70,'C'+str(j));c=cells.copy();v=np.zeros(len(c))
  for sp in ['human','mouse']:
   mask=c.species_key.eq(sp);cap=float(np.quantile(c.loc[mask,program],.99));v[mask]=np.minimum(c.loc[mask,program]/cap,1);scales.append(dict(species=sp,program=program,percentile_99=cap))
  c['relative_program']=v;c=c.sort_values('relative_program');scatter_map(ax,c,c.relative_program,cmap=CM,vmin=0,vmax=1)
  page.text(78+j*148,290,program,ha='center',g='Cvars')
 cb=page.ax(172,310,108,5,'Ckey');plt.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,1),cmap=CM),cax=cb,orientation='horizontal',ticks=[0,.5,1]);cb.set_xlabel('Relative program expression',labelpad=2)
 pd.DataFrame(scales).to_csv(R/'source_data/display_program_scales.tsv',sep='\t',index=False)
 page.text(12,332,'D',fontsize=10,weight='bold');page.text(308,332,'E',fontsize=10,weight='bold')
 ax=page.ax(86,372,144,102,'D');edge=read('S2C_route_expression.tsv');fullpoints=[]
 for i,(h,m) in enumerate((h,m) for h in ROUTES[:4] for m in ROUTES[4:]):
  q=edge[edge.human.eq(h)&edge.mouse.eq(m)];orig=q[q.setting.eq('Original')].iloc[0];full=q[q.setting.eq('Full variants')];valid=full[full.status.eq('OK')]
  ax.plot([orig.null_q025-orig.null_mean,orig.null_q975-orig.null_mean],[i,i],color='#A8ADB2',lw=2,zorder=1)
  jitter=np.linspace(-.17,.17,len(valid));ax.scatter(valid.excess,i+jitter,s=9,c=BLUE,alpha=.7,lw=0)
  ax.scatter(orig.excess,i,s=16,marker='D',color=INK,zorder=4)
  passed=int(((valid.excess>0)&(valid.observed>valid.null_q975)).sum());ne=int(full.status.ne('OK').sum())
  ax.text(1.04,i,str(passed)+'/12'+('*' if ne else ''),transform=ax.get_yaxis_transform(),va='center')
  fullpoints.append(dict(human=h,mouse=m,full_pass=passed,full_total=12,full_NE=ne,original_excess=orig.excess))
  if i in [2,5,8]:ax.axhline(i+.5,color='#E1E5E8',lw=.4)
 ax.axvline(0,color='#A8ADB2',lw=.5,ls='--');ax.set_yticks([]);ax.set_ylim(11.5,-.5);ax.set_xlim(-.006,.058);ax.set_xticks([0,.02,.04]);ax.set_xlabel('Excess held-out correlation',labelpad=3)
 page.text(14,358,'Human');page.text(56,358,'Mouse');page.text(234,358,'Support')
 for xx,human in [(26,True),(68,False)]:
  bar=page.ax(xx,372,10,102,'Dtracks');bar.set_xlim(0,1);bar.set_ylim(11.5,-.5);bar.axis('off')
  for i in range(12):
   route=ROUTES[i//3] if human else ROUTES[4+i%3]
   bar.add_patch(Rectangle((0,i-.5),1,1,color=UNIT_COL[route],lw=0))
   if (not human) or i in [2,5,8]:bar.axhline(i+.5,color='white',lw=.7)
 pd.DataFrame(fullpoints).to_csv(R/'source_data/route_support_counts.tsv',sep='\t',index=False)
 geo=read('S2A_state_geometry.tsv');near=read('S2B_neighbor_all_scales.tsv');geo=geo[geo.setting.eq('Original')];near=near[near.setting.eq('Original')&near.k.eq(30)]
 for j in range(2):
  aa=page.ax(339+j*62,372,42,102,'E'+str(j))
  for i,r in enumerate(ROUTES):
   v=geo.loc[geo.facet_id.eq(r),'pair_distance_spearman'].iloc[0] if j==0 else near.loc[near.route.eq(r),'recall_equal_library'].iloc[0]*100
   aa.scatter(v,i,s=14,color=INK,zorder=3)
   if j==1:aa.scatter(near.loc[near.route.eq(r),'chance'].iloc[0]*100,i,s=12,marker='x',color=MUTED,lw=.6)
  aa.set_xlim((0,1.04) if j==0 else (0,25));aa.set_ylim(6.5,-.5);aa.set_xticks([0,1] if j==0 else [0,20]);aa.set_yticks([]);aa.set_xlabel('Distance ρ' if j==0 else 'Recall %',labelpad=3)
 bar=page.ax(322,372,8,102,'Etrack');bar.set_xlim(0,1);bar.set_ylim(6.5,-.5);bar.axis('off')
 for i,r in enumerate(ROUTES):bar.add_patch(Rectangle((0,i-.45),1,.9,color=UNIT_COL[r],lw=0))
 legend_ax=page.ax(13,344,237,9,'Dkey');legend_ax.set_xlim(0,237);legend_ax.set_ylim(0,1);legend_ax.axis('off')
 legend_ax.scatter(3,.5,marker='D',s=16,color=INK);legend_ax.text(10,.5,'Original',va='center')
 legend_ax.scatter(58,.5,s=9,color=BLUE);legend_ax.text(65,.5,'Refit',va='center')
 legend_ax.plot([105,116],[.5,.5],color='#A8ADB2',lw=2);legend_ax.text(120,.5,'Null central 95%',va='center')
 page.text(369,344,'× Chance',g='Ekey')
 # One full author/year + route key serves both levels of D and all rows of E.
 for route,yy in zip(ROUTES[:4],[500,513,526,551]):colour_key(page,13,yy,UNIT_LABEL[route],UNIT_COL[route])
 for route,yy in zip(ROUTES[4:],[500,513,526]):colour_key(page,285,yy,UNIT_LABEL[route],UNIT_COL[route])
 page.text(285,551,'*2 refits not evaluable.',g='Ddefinition')
 page.gutters=[('A0','Akey','v'),('Akey','B','h'),('Akey','C0','v'),('B','Bkey','v'),('Bkey','C1','v'),('Cvars','Ckey','v'),('C0','C1','h'),('C1','C2','h'),('Ckey','Dkey','v'),('D','Etrack','h'),('E0','E1','h'),('D','unit_key','v'),('E0','unit_key','v')]
 return page.save('Figure_S3')
if __name__=='__main__':
 qa={'Figure_S1':marker_figure('human','Figure_S1'),'Figure_S2':marker_figure('mouse','Figure_S2')}
 if (R/'source_data/display_expression.tsv.gz').exists():qa['Figure_S3']=integration_figure()
 (R/'layout_qa.json').write_text(json.dumps(qa,indent=2)+'\n');print(json.dumps(qa,indent=2))
