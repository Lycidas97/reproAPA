"""Compound supplement via render_candidate: fixed type, evidence and height gates."""
from pathlib import Path
from collections import Counter
import copy,hashlib,json,sys,textwrap,subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager,colors,ticker
from matplotlib.patches import Rectangle,PathPatch
from matplotlib.path import Path as MPath
from matplotlib.lines import Line2D
import pymupdf as fitz
from render_common import Renderer,MM
from render_theme_studies import prepare_dataset_labels,dataset_color

ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/'notes/2026-10-03-polyapipe-integrated-reference'
DESIGN=BASE/'supplement-redraw-design-v3-20261004'
INK='#263D4A';MUTED='#586870';BLUE='#0072B2';ORANGE='#D55E00';GREY='#A9B3BB'
PAL=[BLUE,ORANGE,'#8766A4','#009E73','#CC79A7','#5D6D7E']
STATE={'S':'#E69F00','R':BLUE,'E':'#009E73','P':'#CC79A7','M':'#8C6D31'}
STUDY={'CIP':'#CF4C00','WOT':'#207450','OSK':'#255FBC','GSE118258':'#753399','GSE150311':'#648290','GSE221739':'#30363D','178325_0618':'#A82C60','178325_0809':'#E0829F'}
SHORT={'GSE118258':'Xing','GSE150311':'Liu','GSE103221':'Guo2019','GSE115943':'Schiebinger2019','WOT_GSE115943':'WOT','GSE171127':'Jaber','GSE178325':'Guan','178325_0618':'Guan–0618','178325_0809':'Guan–0809','GSE221739':'Panariello','H1':'Xing','H2':'Liu','H3':'Guan','H4':'Panariello','C':'CIP','W':'WOT'}
CX=[('human','SP'),('human','MP'),('mouse','SP'),('mouse','MP')]
CTX=[('human','P_S'),('human','P_M'),('mouse','P_S'),('mouse','P_M')]
CL=['Human P–S','Human P–M','Mouse P–S','Mouse P–M']
FONT=Path('/Users/lycidas/Library/Fonts')
for p in FONT.glob('nimbus-sans-l.*.otf'):font_manager.fontManager.addfont(p)
FM=fitz.Font(fontfile=str(FONT/'nimbus-sans-l.regular.otf'))
plt.rcParams.update({'font.family':'Nimbus Sans L','font.size':8,'axes.labelsize':8,'axes.titlesize':8,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'pdf.fonttype':42,'axes.linewidth':.6,'xtick.major.size':2,'ytick.major.size':2,'xtick.major.pad':2,'ytick.major.pad':2,'axes.spines.top':False,'axes.spines.right':False,'text.color':INK,'axes.labelcolor':INK,'axes.edgecolor':INK,'xtick.color':INK,'ytick.color':INK,'figure.max_open_warning':0})
DIV=colors.LinearSegmentedColormap.from_list('theme_delta',[BLUE,'#F7F7F7',ORANGE]);DIV.set_bad('#DEDEDE')
# Heights follow rows/labels, not a full-page target. Fixed typography is unchanged.
HEIGHTS={1:[136.5,39],2:[133,42],3:[36,40,48],4:[34,36,41],5:[34,35,32],6:[34,42,52],7:[37,49,36],8:[77,32,47],9:[88,37],10:[34,34,26],11:[32,36,39],12:[51,37,42],13:[26,40,37],14:[24,34,36],15:[30,91,34],16:[31,66,33],17:[96,79],18:[21,37,34,31],19:[30,41,77],20:[29,78,42],21:[30,40,42,31],22:[42,32,40],23:[39,40,47],24:[34,36,39],25:[52,34,34],26:[45,47,34]}

# Packed rectangles are complete content boxes, including axis labels and legends.
# Entries: logical block, component, x, top, width, height (mm).
PACKED={
'D06':[('A',0,1.5,5.5,45,35),('A',1,50,5.5,59,35),('B',1,112.5,5.5,46,35),('B',0,1.5,48,65,53),('C',0,70,48,56,53),('C',1,129.5,48,29,53)],
'D10':[('A',0,1.5,5.5,54,40),('A',1,59,5.5,46,40),('C',1,108.5,5.5,50,40),('B',0,1.5,53,99,38),('B',1,104,53,54.5,38)],
'D14':[('A',1,1.5,5.5,50,34),('B',1,55,5.5,50,34),('C',1,108.5,5.5,50,34),('B',0,1.5,47,99,44),('C',0,104,47,54.5,44)],
'D16':[('A',0,1.5,5.5,66,39),('A',1,71,5.5,38,39),('C',0,112.5,5.5,46,39),('B',0,1.5,52,52,66),('B',1,57,52,55,66),('C',1,115.5,52,43,66)],
'D18':[('A',1,1.5,5.5,42,36),('B',1,47,5.5,51,36),('C',1,101.5,5.5,57,36),('B',0,1.5,49,104,39),('D',0,109,49,49.5,39),('C',0,1.5,95.5,82,38),('D',1,87,95.5,71.5,38)],
'D20':[('B',0,1.5,5.5,105,83),('A',1,110,5.5,48.5,32),('B',1,110,45,48.5,43.5),('C',0,1.5,96,157,43)],
'D21':[('A',0,1.5,5.5,74,33),('A',1,79,5.5,39,33),('D',1,121.5,5.5,37,33),('B',0,1.5,46,104,50.5),('B',1,109,46,49.5,50.5),('C',0,1.5,104,69,45),('C',1,74,104,37.5,45),('D',0,115,104,43.5,45)],
'D22':[('A',0,1.5,5.5,70,44),('A',1,75,5.5,43,44),('C',1,121.5,5.5,37,44),('B',0,1.5,57,64.5,44),('C',0,69.5,57,89,44)],
'D23':[('A',0,1.5,5.5,67,43),('A',1,72,5.5,29.5,43),('C',1,105,5.5,53.5,43),('B',0,1.5,56,43,52),('B',1,48,56,49,52),('C',0,100.5,56,58,52)],
'D24':[('A',0,1.5,5.5,50,43),('A',1,55,5.5,47,43),('B',0,105.5,5.5,53,43),('B',1,1.5,56,50,45),('C',0,55,56,60,45),('C',1,118.5,56,40,45)],
'D26':[('A',0,1.5,5.5,80,47),('A',1,85,5.5,35,47),('B',1,123.5,5.5,35,47),('B',0,1.5,60,99,49),('C',0,104,60,54.5,49)]}
CONTENT_FIT=set(PACKED)|{'D07','D12','D15','D17'}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def short(x):return SHORT.get(str(x),str(x))
def context(sp,co):return ('H ' if sp=='human' else 'M ')+co.replace('_','–').replace('-','–').replace('SP','P–S').replace('MP','P–M')
def endpoint(x):return x[(x.contrast.eq('P_S')&x.design.eq('true_source'))|(x.contrast.eq('P_M')&x.design.eq('same_library'))]
def mmwidth(s,pt=8):return max((FM.text_length(t,fontsize=pt) for t in str(s).split('\n')),default=0)*25.4/72
def wrap(s,w,pt=8):
    lines=[]
    for part in str(s).split('\n'):
        line=''
        for word in part.split():
            test=(line+' '+word).strip()
            if line and mmwidth(test,pt)>w:lines.append(line);line=word
            else:line=test
        lines.append(line)
    return '\n'.join(lines)
def resolve(ref):
    for prefix,root in [('root:',ROOT),('base:',BASE),('hyp:',BASE/'validation/hypothesis-extensions-ABCDE'),('axis:',BASE/'validation/axis-exploration'),('prepared:',DESIGN/'prepared-inputs')]:
        if ref.startswith(prefix):return root/ref[len(prefix):]
    return BASE/'statistics_results/v2_primary'/ref
def legend(ax,**kw):return ax.legend(frameon=False,handlelength=1.1,handletextpad=.35,columnspacing=.65,labelspacing=.2,borderpad=0,**kw)
def zeros(ax,both=False):
    ax.axvline(0,color=GREY,ls=':',lw=.6,zorder=0)
    if both:ax.axhline(0,color=GREY,ls=':',lw=.6,zorder=0)

class ThemeRenderer(Renderer):
    mmwidth=staticmethod(mmwidth)
    wrap=staticmethod(wrap)
    def __init__(self,config,out):
        super().__init__(config,out)
        self.contract=json.loads((DESIGN/'redraw-contract.json').read_text())
        if (self.out/'review-contract.json').exists():self.contract=json.loads((self.out/'review-contract.json').read_text())
        self.defs={f['id']:copy.deepcopy(f) for f in self.contract['figures']}
        self.pins={str(Path(r.path).resolve()):r.sha256 for r in pd.read_csv(DESIGN/'source-interface.tsv',sep='\t').itertuples()}
        self.pins.update({v.source_path:v.recorded_sha256 for v in pd.read_csv(ROOT/'notes/2026-10-05-supplement-source-map/source-inventory.tsv',sep='\t').itertuples() if str(v.reference).startswith('root:notes/2026-09-24-integration-top-labels/')})
        self.cache={};self.generated=[]
        self.bundle=None
        if self.config.get('result_records'):
            from render_theme_bundle import SupplementBundle
            self.bundle=SupplementBundle(self.config['result_records'],self)
        for d in ['Fig','panel-data','captions','verification']: (self.out/d).mkdir(exist_ok=True)
        snap=self.out/'before'/'protected-input-hashes.json'
        if not snap.exists():
            protected=[str(p) for folder in ['Fig','main','supplement'] for p in (ROOT/folder).rglob('*') if p.is_file() and p.suffix in ['.pdf','.tex']]
            snap.parent.mkdir(exist_ok=True);snap.write_text(json.dumps({p:sha(p) for p in protected},indent=2))
        self.protected=json.loads(snap.read_text())
    def begin(self,fid):
        self.fid=fid;self.f=copy.deepcopy(self.defs[fid]);self.frame_counts=Counter();self.axes_meta=[];self.evidence={};self.notes={};self.active=None;self.sources=set();self.panels=[];self.panel_axes={};self.panel_titles={};self.y_links=[];self.colorbars=[]
        hs=HEIGHTS.get(int(fid[1:]),[45]*len(self.f['panels']));hs=hs+[45]*(len(self.f['panels'])-len(hs));assert len(hs)==len(self.f['panels']);y=1.5
        for p,h in zip(self.f['panels'],hs):
            p['rectangle_mm']=[1.5,y,157,h+4];cs=p['components'];gap=3.5
            if p['component_layout']=='vertical':
                # Marker grids retain >=2.8mm row pitch; support uses its own short strip.
                hh=[h-gap-21,17];yy=y+4
                for c,ch in zip(cs,hh):c['rectangle_mm']=[1.5,yy,157,ch];yy+=ch+gap+4
            elif p['component_layout']=='left_tall':
                cs[0]['rectangle_mm']=[1.5,y+4,92,h];cs[1]['rectangle_mm']=[97,y+4,61.5,39];cs[2]['rectangle_mm']=[97,y+4+46.5,61.5,h-46.5]
            else:
                ww=np.array([c['min_width_mm'] for c in cs],float);ww*= (157-gap*(len(cs)-1))/ww.sum();x=1.5
                for c,w in zip(cs,ww):c['rectangle_mm']=[x,y+4,float(w),h];x+=w+gap
            y+=h+7.5
        if fid in PACKED:
            assigned=set()
            for letter,k,x0,y0,ww,hh in PACKED[fid]:
                pp=next(p for p in self.f['panels'] if p['id']==fid+letter);pp['components'][k]['rectangle_mm']=[x0,y0,ww,hh];assigned.add((letter,k))
            for pp in self.f['panels']:
                cs=[c['rectangle_mm'] for k,c in enumerate(pp['components']) if (pp['id'][-1],k) in assigned]
                if cs:
                    x0=min(a[0] for a in cs);y0=min(a[1] for a in cs)-4;x1=max(a[0]+a[2] for a in cs);y1=max(a[1]+a[3] for a in cs);pp['rectangle_mm']=[x0,y0,x1-x0,y1-y0]
            y=max(row[3]+row[5] for row in PACKED[fid])+3.5
        self.compact_plan=None
        if (self.out/'compact-row-plan.json').exists():
            plans=json.loads((self.out/'compact-row-plan.json').read_text())
            self.compact_plan=plans.get(fid)
        self.height=max(v[1]+v[3] for v in self.compact_plan.values())+1.5 if self.compact_plan else y-2.0
        self.width=160;assert self.height<=190,(fid,self.height)
        self.f['height_mm']=self.height;self.number='Figure '+fid;self.fig=plt.figure(figsize=(160*MM,self.height*MM))
        self.pdefs={p['id']:p for p in self.f['panels']}
    def panel(self,letter):
        self.active=self.fid+letter;assert self.active in self.pdefs
        self.panels.append(self.active);self.panel_axes[self.active]=[];self.evidence[self.active]={'sources':{},'exports':[]};self.notes[self.active]=[]
        p=self.pdefs[self.active];x,y,w,h=p['rectangle_mm'];# Final letters are assigned to actual shared-y groups, not logical topic rows.
        return self
    def read(self,ref,**kw):
        if self.bundle:return self.bundle.read(ref,**kw)
        path=resolve(ref).resolve();assert path.is_file(),path
        if str(path) not in self.inputs:
            digest=sha(path)
            if str(path) in self.pins:assert digest==self.pins[str(path)],ref
            self.inputs[str(path)]=digest
        self.sources.add(ref);self.evidence[self.active]['sources'][ref]=dict(path=str(path),sha256=self.inputs[str(path)])
        key=(str(path),str(kw))
        if key not in self.cache:self.cache[key]=pd.read_csv(path,sep=kw.pop('sep','\t'),float_precision='round_trip',**kw)
        return self.cache[key].copy()
    def npz(self,ref):
        p=self.bundle.path(ref) if self.bundle else resolve(ref).resolve();digest=sha(p)
        if str(p) in self.pins:assert digest==self.pins[str(p)]
        self.inputs[str(p)]=digest
        if not self.bundle:self.evidence[self.active]['sources'][ref]=dict(path=str(p),sha256=digest)
        return np.load(p)
    def json(self,ref):
        p=self.bundle.path(ref) if self.bundle else resolve(ref).resolve();digest=sha(p)
        if str(p) in self.pins:assert digest==self.pins[str(p)]
        self.inputs[str(p)]=digest
        if not self.bundle:self.evidence[self.active]['sources'][ref]=dict(path=str(p),sha256=digest)
        if str(p) not in self.cache:self.cache[str(p)]=json.loads(p.read_text())
        return self.cache[str(p)]
    def record(self,name,df,note='Source rows or explicit descriptive transformation; no new inference'):
        path=self.out/'panel-data'/self.active;path.mkdir(exist_ok=True);dest=path/(name+'.tsv.gz')
        df.to_csv(dest,sep='\t',index=True)
        self.evidence[self.active]['exports'].append(dict(path=str(dest.relative_to(self.out)),rows=len(df),sha256=sha(dest),operation=note))
        return df
    def note(self,s):self.notes[self.active].append(s)
    def ax(self,k,sub=0,n=None,labels=None,kind=None,pad=None,vertical=False):
        c=self.pdefs[self.active]['components'][k];x,y,w,h=c['rectangle_mm'];n=c['axes'] if n is None else n
        family=kind or c['canonical_family'];assert family in ['matrix','bar','distribution','scatter','line','point_interval','genomic','embedding','structure','table']
        self.frame_counts[family]+=1
        if n>1:
            if vertical:h=(h-4*(n-1))/n;y+=sub*(h+4)
            elif self.active=='D20B' and k==0:
                total=w-5;w=total*(.65 if sub==0 else .35);x+=0 if sub==0 else total*.65+5
            else:w=(w-5*(n-1))/n;x+=sub*(w+5)
        if pad is None:
            left=13 if family in ['scatter','line','distribution'] else 7
            if labels:left=min(max(mmwidth(t) for t in labels)+2,w*.62)
            if family in ['structure','table','embedding']:pad=(1,1,1,1)
            else:pad=(left,3,2,10)
        left,top,right,bottom=pad;assert w-left-right>5 and h-top-bottom>3,(self.active,k,w,h,pad)
        box=[x+left,y+top,w-left-right,h-top-bottom];ax=self.fig.add_axes([box[0]/160,1-(box[1]+box[3])/self.height,box[2]/160,box[3]/self.height])
        if family in ['structure','table','embedding']:ax.set_axis_off()
        else:
            ax.xaxis.set_major_locator(ticker.MaxNLocator(4));ax.yaxis.set_major_locator(ticker.MaxNLocator(4))
        self.panel_axes[self.active].append(ax);self.axes_meta.append(dict(panel=self.active,component=k,family=family,axis=ax,component_rect=[x,y,w,h]))
        return ax
    def table(self,k,headers,rows,widths=None):
        ax=self.ax(k,kind='table');ax.set(xlim=(0,1),ylim=(len(rows)+1,0));widths=np.ones(len(headers)) if widths is None else np.array(widths);xx=np.cumsum(np.r_[0,widths/widths.sum()])
        for j,s in enumerate(headers):ax.text(xx[j]+.015,.5,wrap(str(s),ax.get_position().width*160*(xx[j+1]-xx[j])-.8),va='center',fontsize=8,fontweight='bold')
        for i,row in enumerate(rows):
            if i%2==0:ax.axhspan(i+1,i+2,color='#F0F3F4')
            for j,s in enumerate(row):ax.text(xx[j]+.015,i+1.5,str(s),va='center',fontsize=8)
        self.record(f'table-{k}',pd.DataFrame(rows,columns=headers));return ax
    def diagram(self,k,labels,counts=None,vertical=False,arrows=True):
        ax=self.ax(k,kind='structure');n=len(labels);ax.set(xlim=(0,1),ylim=(0,1))
        for i,lab in enumerate(labels):
            xx=.5 if vertical else (i+.5)/n;yy=1-(i+.5)/n if vertical else .55
            ax.text(xx,yy,str(lab)+('' if counts is None else '\n'+str(counts[i])),ha='center',va='center',fontsize=8)
            if i<n-1 and arrows:
                if vertical:ax.annotate('',xy=(.5,yy-.75/n),xytext=(.5,yy-.3/n),arrowprops=dict(arrowstyle='->',color=GREY,lw=.7))
                else:ax.annotate('',xy=(xx+.64/n,yy),xytext=(xx+.35/n,yy),arrowprops=dict(arrowstyle='->',color=GREY,lw=.7))
        self.record(f'definition-{k}',pd.DataFrame({'label':labels,'value':counts if counts is not None else ['']*n}),'Concept/count schematic; arrows do not imply biological causation');return ax
    def dots(self,k,df,value,labels,lo=None,hi=None,series=None,xlabel='',reference=0,q=None,pad=None,sub=0,n=None):
        ax=self.ax(k,sub,n,labels=labels,kind='point_interval',pad=pad);groups=list(df[series].drop_duplicates()) if series else [None]
        for i in range(len(labels)):
            if i%2==0:ax.axhspan(i-.5,i+.5,color='#F2F4F5',zorder=-3)
        for j,g in enumerate(groups):
            x=df[df[series].eq(g)] if series else df;y=x['_row'].to_numpy()+(j-(len(groups)-1)/2)*min(.22,.7/max(len(groups),1));v=x[value].to_numpy(float);good=np.isfinite(v);col=PAL[j%len(PAL)] if series else INK
            if lo and hi:
                low=x[lo].to_numpy(float);high=x[hi].to_numpy(float);ok=good&np.isfinite(low)&np.isfinite(high);ax.hlines(y[ok],low[ok],high[ok],color=col if series else GREY,lw=1)
            fill=np.where(x[q].to_numpy()<.05,col,'white') if q else col
            ax.scatter(v[good],y[good],s=12,marker=['o','s','^'][j%3],facecolors=fill[good] if isinstance(fill,np.ndarray) else fill,edgecolors=col,lw=.7,label=str(g) if g is not None else None,zorder=3)
            for yy in y[~good]:ax.text(.98,yy,'NA',transform=ax.get_yaxis_transform(),ha='right',va='center',color=MUTED,fontsize=7)
        ax.set_yticks(range(len(labels)),labels,fontsize=8);ax.set_ylim(len(labels)-.45,-.55);ax.set_xlabel(xlabel)
        if reference is not None:ax.axvline(reference,color=GREY,ls=':',lw=.6)
        if series:legend(ax,loc='best',ncol=min(3,len(groups)))
        self.record(f'dots-{k}-{sub}',df);return ax
    def bars(self,k,frame,xlabel='',colorset=None,fraction=False,pad=None):
        labels=list(map(str,frame.index));ax=self.ax(k,labels=labels,kind='bar',pad=pad);left=np.zeros(len(frame))
        for j,col in enumerate(frame):
            v=frame[col].to_numpy(float);ax.barh(range(len(frame)),v,left=left,color=(colorset or PAL)[j%len(colorset or PAL)],height=.65,label=str(col));left+=np.nan_to_num(v)
        ax.set_yticks(range(len(frame)),labels,fontsize=8);ax.set_ylim(len(frame)-.4,-.6);ax.set_xlabel(xlabel);ax.set_xlim(0,1 if fraction else max(1,left.max())*1.1)
        if len(frame.columns)>1:legend(ax,loc='lower right',bbox_to_anchor=(1,1.01),ncol=min(3,len(frame.columns)))
        for i,v in enumerate(left):
            if v==0:ax.text(.015,i,'0',transform=ax.get_yaxis_transform(),va='center',fontsize=7)
        self.record(f'bars-{k}',frame);return ax
    def ecdf(self,k,groups,xlabel='',log=False,pad=None):
        ax=self.ax(k,kind='distribution',pad=pad)
        for j,(lab,v) in enumerate(groups):
            v=np.sort(np.asarray(v,float));v=v[np.isfinite(v)];assert len(v)>0
            ax.step(v,np.arange(1,len(v)+1)/len(v),where='post',color=PAL[j%len(PAL)],ls=['-','--',':'][j//len(PAL)%3],label=lab,lw=.9);self.record(f'ecdf-{k}-{j}',pd.DataFrame({'value':v}))
        ax.set(xlabel=xlabel,ylabel='Cumulative fraction',ylim=(0,1),yticks=[0,.5,1]);
        if log:ax.set_xscale('log')
        legend(ax,loc='best');return ax
    def matrix(self,k,frame,label='',limit=None,cmap=None,sub=0,n=None,pad=None,qframe=None):
        ax=self.ax(k,sub,n,labels=list(frame.index),kind='matrix',pad=pad);v=frame.to_numpy(float);limit=limit or np.nanmax(abs(v));vmin=0 if cmap else -limit
        im=ax.imshow(np.ma.masked_invalid(v),aspect='auto',interpolation='nearest',cmap=cmap or DIV,vmin=vmin,vmax=limit)
        ax.set_xticks(range(len(frame.columns)),frame.columns,fontsize=7);ax.xaxis.tick_top();ax.tick_params(length=0,pad=2);ax.set_yticks(range(len(frame)),frame.index,fontsize=8)
        for s in ax.spines.values():s.set_visible(False)
        for i,j in zip(*np.where(~np.isfinite(v))):ax.text(j,i,'×',ha='center',va='center',color=MUTED,fontsize=7)
        if qframe is not None:
            for i,j in zip(*np.where(qframe.to_numpy(float)<.05)):ax.add_patch(Rectangle((j-.48,i-.48),.96,.96,fill=False,ec=INK,lw=.6))
        pos=ax.get_position();cb=self.fig.add_axes([pos.x0,pos.y0-4/self.height,min(pos.width*.4,.12),1.4/self.height]);self.fig.colorbar(im,cax=cb,orientation='horizontal',ticks=[vmin,limit]);cb.tick_params(labelsize=7,length=1,pad=1);cb.text(1.12,.5,label,transform=cb.transAxes,fontsize=7,va='center');self.panel_axes[self.active].append(cb);self.colorbars.append((ax,cb))
        self.record(f'matrix-{k}-{sub}',frame);return ax
    def share_y(self,axes,basis,hide=True):
        axes=list(axes);assert len(axes)>1
        lo=min(min(a.get_ylim()) for a in axes);hi=max(max(a.get_ylim()) for a in axes);rev=axes[0].yaxis_inverted()
        for a in axes:a.set_ylim((hi,lo) if rev else (lo,hi))
        if hide:
            for a in axes[1:]:
                a.tick_params(axis='y',labelleft=False);a.set_ylabel('')
                m=next(m for m in self.axes_meta if m['axis'] is a);x,y,w,h=m['component_rect'];p=a.get_position();left=(x+2)/160;a.set_position([left,p.y0,p.x1-left,p.height])
        self.y_links.append(dict(axes=axes,basis=basis))
    def row_distribution(self,k,groups,xlabel,pad=None):
        ax=self.ax(k,kind='distribution',labels=[lab for lab,_ in groups],pad=pad)
        for i,(label,values) in enumerate(groups):
            v=np.asarray(values,float);v=v[np.isfinite(v)];assert len(v)
            # Histogram silhouettes preserve multimodality and exact zeros; no KDE inference.
            edges=np.histogram_bin_edges(v,bins=min(32,max(8,int(np.sqrt(len(v))))))
            freq,_=np.histogram(v,edges);cent=(edges[:-1]+edges[1:])/2;half=.28*freq/max(freq.max(),1)
            ax.fill_between(cent,i-half,i+half,color=BLUE,alpha=.25,lw=.4)
            ax.scatter(v,i+np.random.default_rng(20261004+i).uniform(-.14,.14,len(v)),s=2,c=BLUE,alpha=.3,lw=0,rasterized=True)
            qq=np.quantile(v,[.25,.5,.75]);ax.plot(qq[[0,2]],[i,i],color=INK,lw=1);ax.plot(qq[1],i,'|',color=INK,ms=4)
            self.record(f'row-distribution-{k}-{i}',pd.DataFrame({'value':v}),'All finite observations; histogram silhouette scaled within row, raw marks, median/IQR. Width is not between-row sample size.')
        ax.set(yticks=range(len(groups)),yticklabels=[lab for lab,_ in groups],ylim=(len(groups)-.5,-.5),xlabel=xlabel);ax.axvline(0,color=GREY,ls=':',lw=.6);return ax
    def _sync_y_geometry(self):
        for link in self.y_links:
            aa=link['axes'];bottom=max(a.get_position().y0 for a in aa);top=min(a.get_position().y1 for a in aa)
            assert (top-bottom)*self.height>3
            for a in aa:
                p=a.get_position();a.set_position([p.x0,bottom,p.width,top-bottom])
        for parent,cb in self.colorbars:
            p=parent.get_position();old=cb.get_position();cb.set_position([p.x0,p.y0-4/self.height,min(p.width*.4,.12),old.height])
    def _letter_groups(self,draw=True,with_axes=False):
        allaxes=[m['axis'] for m in self.axes_meta];groups=[{a} for a in allaxes]
        for link in self.y_links:
            union=set(link['axes']);hits=[g for g in groups if g&union]
            for g in hits:union|=g;groups.remove(g)
            groups.append(union)
        def box(g):
            mm=[m['component_rect'] for m in self.axes_meta if m['axis'] in g];return [min(v[0] for v in mm),min(v[1] for v in mm),max(v[0]+v[2] for v in mm),max(v[1]+v[3] for v in mm)]
        groups.sort(key=lambda g:(round(box(g)[1],1),box(g)[0]));mapping=[]
        for i,g in enumerate(groups):
            letter=chr(65+i);x,y,right,bottom=box(g)
            if draw:self.fig.text(x/160,1-(y-3.6)/self.height,letter,fontsize=11,fontweight='bold',va='top')
            metas=[m for m in self.axes_meta if m['axis'] in g];links=[a['basis'] for a in self.y_links if set(a['axes'])<=g]
            mapping.append(dict(panel=self.fid+letter,source_logical_panels=list(dict.fromkeys(m['panel'] for m in metas)),components=[dict(logical_panel=m['panel'],component=m['component'],family=m['family']) for m in metas],composite=len(g)>1,shared_y_basis=links,rectangle_mm=[x,y-4,right-x,bottom-y+4]))
        return list(zip(mapping,groups)) if with_axes else mapping
    def _align_row_geometry(self):
        # Equal physical frames do not imply a shared variable or a compound panel.
        groups=[]
        for m in self.axes_meta:
            a=m['axis']
            if not (a.axison or getattr(a,'_row_alignable',False)):continue
            _,y,_,h=m['component_rect'];key=(round(y,4),round(h,4))
            hit=next((g for k,g in groups if k==key),None)
            if hit is None:groups.append((key,{a}))
            else:hit.add(a)
        united=[g for _,g in groups]
        for link in self.y_links:
            union=set(link['axes']);hits=[g for g in united if g&union]
            for g in hits:union|=g;united.remove(g)
            united.append(union)
        self.row_align_checks=[]
        for axes in united:
            if len(axes)<2:continue
            bottom=max(a.get_position().y0 for a in axes);top=min(a.get_position().y1 for a in axes)
            assert (top-bottom)*self.height>3,(self.fid,'row too short')
            for a in axes:
                p=a.get_position();a.set_position([p.x0,bottom,p.width,top-bottom])
            actual=[a.get_position() for a in axes]
            error=max(max(p.y0 for p in actual)-min(p.y0 for p in actual),max(p.y1 for p in actual)-min(p.y1 for p in actual))*self.height
            self.row_align_checks.append(dict(components=[f"{m['panel']}:{m['component']}" for m in self.axes_meta if m['axis'] in axes],top_mm=(1-top)*self.height,bottom_mm=(1-bottom)*self.height,max_y_edge_difference_mm=error))
    def finish(self):
        assert set(self.panels)==set(self.pdefs)
        assert max(self.frame_counts.values())<=2,(self.fid,self.frame_counts)
        # Multi-axis components become composites only where their y dimension agrees.
        for pid in self.panels:
            for k in range(len(self.pdefs[pid]['components'])):
                aa=[m['axis'] for m in self.axes_meta if m['panel']==pid and m['component']==k]
                mm=[m for m in self.axes_meta if m['panel']==pid and m['component']==k]
                if len(aa)>1 and len({m['family'] for m in mm})==1 and mm[0]['family']!='embedding' and not any(set(aa)&set(l['axes']) for l in self.y_links):
                    ylabel=aa[0].get_ylabel();ticks=aa[0].get_yticklabels()
                    if len(set(a.get_ylabel() for a in aa))==1 or mm[0]['family']=='point_interval':self.share_y(aa,ylabel or 'Matched object rows')
        # Determine shared-variable identity before replacing names by equal-width swatches.
        prepare_dataset_labels(self)
        from render_theme_layout import layout
        mapping=layout(self,self._letter_groups(draw=False,with_axes=True))
        assert len(self.fig.legends)==len(self.panel_keys), 'A legend has no panel owner'
        from render_theme_density import census
        census(self)
        flags=self.strict_violations;textboxes=self.strict_textboxes
        (self.out/'verification'/f'{self.fid}.textboxes.json').write_text(json.dumps(textboxes,ensure_ascii=False,indent=2))
        pdf=self.out/'Fig'/f'SuppFig.{self.fid}.pdf';self.fig.savefig(pdf,dpi=600)
        # Cairo embeds the Nimbus CFF outlines using the correct PDF font subtype.
        normalized=pdf.with_suffix('.normalized.pdf')
        proc=subprocess.run(['/opt/homebrew/bin/pdftocairo','-pdf','-origpagesizes',str(pdf),str(normalized)],capture_output=True,text=True)
        assert proc.returncode==0 and normalized.is_file(),proc.stderr
        normalized.replace(pdf)
        # Export default: deduplicate identical PDF objects/streams after font normalization.
        compact=pdf.with_suffix('.compact.pdf')
        with fitz.open(pdf) as cleaned:cleaned.save(compact,garbage=4,deflate=True,clean=True)
        compact.replace(pdf)
        doc=fitz.open(pdf);assert len(doc)==1
        assert abs(doc[0].rect.width*25.4/72-160)<.01
        assert abs(doc[0].rect.height*25.4/72-self.height)<.01
        doc[0].get_pixmap(dpi=150).save(pdf.with_suffix('.png'));doc.close()
        # A separate tiling preview makes the actual compound-panel ownership visible.
        for p in self.f['panels']:
            x,y,w,h=p['rectangle_mm'];self.fig.add_artist(Rectangle((x/160,1-(y+h)/self.height),w/160,h/self.height,transform=self.fig.transFigure,fill=False,ec='#D55E00',lw=.5))
        self.fig.savefig(pdf.with_suffix('.tiling.png'),dpi=150);plt.close(self.fig)
        caption=[]
        for p in self.f['panels']:
            pid=p['id']
            if self.bundle:
                p['sources']=[v['path'] for v in self.evidence[pid]['sources'].values()]
                p['selection']='Current source: '+self.bundle.reader.identity['run_id']+'; selections and denominators are in the exported source rows and notes.'
            caption.append(pid+' '+p['question']+'\n'+p['selection']+'\n'+'\n'.join(self.notes[pid]));self.evidence[pid]['contract_sources']=p['sources']
        (self.out/'captions'/f'{self.fid}.zh.txt').write_text('Actual panel map (shared-y definition)\n'+ '\n'.join(v['panel']+' ← '+','.join(v['source_logical_panels'])+' components '+','.join(str(a['component']) for a in v['components'])+(' · shared y '+ '; '.join(v['shared_y_basis']) if v['composite'] else ' · independent subplot') for v in mapping)+'\n\nSource/evidence notes by logical block\n'+'\n\n'.join(caption))
        receipt=dict(figure=self.fid,width_mm=160,height_mm=self.height,maximum_native_height_mm=190,pdf_sha256=sha(pdf),frame_counts=dict(self.frame_counts),panels=self.evidence,panel_mapping=mapping,dataset_names=self.dataset_names,dataset_label_inventory=self.dataset_label_inventory,legend_ownership=self.legend_ownership,panel_gap_mm=4.0,internal_gap_mm=2.8,pdf_bytes=pdf.stat().st_size,shared_y_checks=[dict(basis=v['basis'],max_y_edge_difference_mm=max(abs(a.get_position().y0-v['axes'][0].get_position().y0)+abs(a.get_position().y1-v['axes'][0].get_position().y1) for a in v['axes'])*self.height,max_y_limit_difference=max(np.max(abs(np.array(a.get_ylim())-v['axes'][0].get_ylim())) for a in v['axes'])) for v in self.y_links],geometry_flags=flags,visual_review=False)
        (self.out/'verification'/f'{self.fid}.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2,default=str))
        (self.out/'verification'/f'{self.fid}.row-alignment.json').write_text(json.dumps(self.row_align_checks,ensure_ascii=False,indent=2))
        (self.out/'verification'/f'{self.fid}.layout.json').write_text(json.dumps(self.f,ensure_ascii=False,indent=2))
        self.generated.append(self.fid);print(f'{self.fid}: {self.height:g} mm; {len(flags)} bounds flags; {dict(self.frame_counts)}',flush=True)

def run(config,out,figures):
    if figures=='package':
        from package_theme import package
        package(out);return
    from render_theme_panels import DRAW
    r=ThemeRenderer(config,out);requested=list(r.defs) if figures=='all' else figures.split(',')
    if (r.out/'review-contract.json').exists():
        from render_theme_review import extend
        DRAW=extend(DRAW)
    for fid in requested:r.begin(fid);DRAW[fid](r);r.finish()
    assert all(Path(p).is_file() and sha(p)==h for p,h in r.protected.items()),'Protected prior artifact changed'
    if r.bundle:r.bundle.write_bindings()
    (r.out/'source-hashes.json').write_text(json.dumps(r.inputs,indent=2))
    print('Completed requested figures. Per-page visual review remains required.',flush=True)
