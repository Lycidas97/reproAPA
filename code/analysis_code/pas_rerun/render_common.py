"""Figure rendering from registered outputs of one native PAS rerun."""
import json
from pathlib import Path
import shutil
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter

from polyapipe_prepare import digest

MM=1/25.4
BLUE='#377EB8';ORANGE='#CC8026';INK='#233D4D';GREY='#A3ADB5'
CLASSES=['C1','C2','C3'];COLORS=[BLUE,'#7664A0',ORANGE]
ROUTES=['GSE118258','GSE150311','178325_0330','178325_0605','178325_0618','178325_0809','GSE221739','OSK','CIP','WOT']
LABELS=dict(zip(ROUTES,['Xing2020','Liu2020','Guan2022–0330','Guan2022–0605','Guan2022–0618','Guan2022–0809','Panariello2023','Guo2019–OSK','Guo2019–CIP','Schiebinger2019–2i']))
DIVERGE=LinearSegmentedColormap.from_list('APA',[BLUE,'#F7F7F7',ORANGE]);DIVERGE.set_bad('#D8DCE0')
HOUSE_INK='#263D4A'
HOUSE_MAP=LinearSegmentedColormap.from_list('APA_house',['#0072B2','#F7F7F7','#D55E00'])
HOUSE_MAP.set_bad('#DEDEDE')
STATE_COLORS=dict(S='#E69F00',R='#0072B2',E='#009E73',P='#CC79A7',M='#8C6D31')
STUDY_STYLE={'WOT':('#207450','o'),'CIP':('#CF4C00','s'),'OSK':('#255FBC','^'),
    'GSE118258':('#753399','o'),'GSE150311':('#648290','s'),'GSE221739':('#30363D','^'),
    '178325_0330':('#795923','D'),'178325_0605':('#B17800','P'),
    '178325_0618':('#A82C60','X'),'178325_0809':('#AB775B','>')}
for family in ['Nimbus Sans L','Arial','DejaVu Sans']:
    try:font_manager.findfont(family,fallback_to_default=False);break
    except ValueError:continue
plt.rcParams.update({'font.family':family,'font.size':7,'axes.labelsize':8,'axes.titlesize':8,
    'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':6.5,'pdf.fonttype':42,'ps.fonttype':42,
    'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.5,'lines.linewidth':.8,
    'text.color':INK,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK,'savefig.facecolor':'white'})


class Renderer:
    def __init__(self,config,out):
        self.config=json.loads(Path(config).read_text());self.root=Path(self.config['analysis_root'])
        self.out=Path(out);self.out.mkdir(parents=True,exist_ok=True)
        self.inputs={str(Path(config).resolve()):digest(config)};self.receipts={};self.frames={};self.records=[]

    def artifact(self,relative):
        node=Path(relative).parts[0]
        if node not in self.receipts:
            p=self.root/node/'result.json';r=json.loads(p.read_text())
            assert r['node']==node and r['status']=='complete' and all(r['checks'].values())
            assert r['run_id']==self.config['run_id'] and r['input_fingerprint']==self.config['input_fingerprint']
            self.receipts[node]={a['path']:a['sha256'] for a in r['artifacts']};self.inputs[str(p)]=digest(p)
        path=self.root/relative;short=str(Path(relative).relative_to(node));assert short in self.receipts[node],relative
        if str(path) not in self.inputs:
            assert digest(path)==self.receipts[node][short],relative
            self.inputs[str(path)]=self.receipts[node][short]
        self.sources.add(relative)
        return path

    def read(self,relative,**kwargs):
        path=self.artifact(relative)
        if relative not in self.frames:self.frames[relative]=pd.read_csv(path,sep='\t',float_precision='round_trip',**kwargs)
        return self.frames[relative].copy()

    def resource(self,name):
        item=self.config['resources'][name];p=Path(item['path']);assert digest(p)==item['sha256'],name
        copy=self.out/'resources'/name/p.name;copy.parent.mkdir(parents=True,exist_ok=True)
        if copy.exists():assert digest(copy)==item['sha256'],name
        else:shutil.copy2(p,copy)
        self.inputs[str(p)]=item['sha256'];self.sources.add('resource:'+name);return p

    def start(self,number,width,height):
        self.number=number;self.sources=set();self.panels=[];self.width=width;self.height=height;self.panel_axes={};self.panel_titles={}
        self.main_layout=None;self.panel_texts={}
        self.fig=plt.figure(figsize=(width*MM,height*MM));return self.fig

    def main(self,n):
        """Use the approved apadetect page and panel anchors for current main figures."""
        path=Path(self.config['layout_lock']);self.inputs[str(path)]=digest(path)
        layout=json.loads(path.read_text())[n-1]
        for font in Path(self.config['font_directory']).glob('*.otf'):font_manager.fontManager.addfont(str(font))
        font_manager.findfont('Nimbus Sans L',fallback_to_default=False)
        plt.rcParams.update({'font.family':'Nimbus Sans L','font.size':8,'axes.labelsize':8,
            'xtick.labelsize':7.2,'ytick.labelsize':7.2,'legend.fontsize':7.2,
            'axes.linewidth':.6,'xtick.major.size':2,'ytick.major.size':2,
            'text.color':HOUSE_INK,'axes.labelcolor':HOUSE_INK,'xtick.color':HOUSE_INK,'ytick.color':HOUSE_INK,
            'savefig.dpi':600})
        self.start(f'Figure {n}',*layout['native_size_mm']);self.main_layout=layout
        grid_path=Path(self.config['layout_grids'])/f'figure{n}/grid.json'
        self.inputs[str(grid_path)]=digest(grid_path)
        self.grid=json.loads(grid_path.read_text())['panels']
        self.number_id=n

    def text(self,x,y,value,**kwargs):
        kwargs.setdefault('va','top');kwargs.setdefault('fontsize',8)
        artist=self.fig.text(x/self.width,1-y/self.height,value,**kwargs)
        self.panel_texts.setdefault(self.panels[-1],[]).append(artist)
        return artist

    def display(self,name,frame):
        folder=self.out/'source-data';folder.mkdir(exist_ok=True)
        frame.to_csv(folder/f'Figure{self.number_id}-{name}.tsv',sep='\t',index=False)

    def legend(self,x,y,handles,ncol=1):
        artist=self.fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(x/self.width,1-y/self.height),
            ncol=ncol,frameon=False,borderaxespad=0,handlelength=1,handletextpad=.4,columnspacing=.8,labelspacing=.3)
        self.panel_texts.setdefault(self.panels[-1],[]).append(artist)
        return artist

    def geometry(self):
        self.fig.canvas.draw();renderer=self.fig.canvas.get_renderer();result={}
        for panel in self.panels:
            bounds=[]
            for ax in self.panel_axes.get(panel,[]):
                box=ax.get_tightbbox(renderer)
                if box is not None:bounds.append(box)
            bounds.extend(t.get_window_extent(renderer) for t in self.panel_texts.get(panel,[]) if t.get_visible())
            if not bounds:continue
            from matplotlib.transforms import Bbox
            b=Bbox.union(bounds);dpi=self.fig.dpi
            actual=[b.x0/dpi*25.4,self.height-b.y1/dpi*25.4,b.width/dpi*25.4,b.height/dpi*25.4]
            target=np.asarray(self.grid[panel]['box_pt'])*25.4/72
            result[panel]=dict(actual_ink_box_mm=actual,master_panel_box_mm=target.tolist(),
                overflow_mm=[max(0,target[0]-actual[0]),max(0,target[1]-actual[1]),
                    max(0,actual[0]+actual[2]-target[0]-target[2]),max(0,actual[1]+actual[3]-target[1]-target[3])])
        return result

    def axes(self,box):
        x,y,w,h=box;ax=self.fig.add_axes([x/self.width,1-(y+h)/self.height,w/self.width,h/self.height])
        panel=self.panels[-1] if self.panels else 'shared';self.panel_axes.setdefault(panel,[]).append(ax)
        return ax

    def label(self,letter,x,y,title=''):
        assert letter not in self.panels;self.panels.append(letter)
        self.panel_titles[letter]=title
        if self.main_layout:
            anchor=next(a for a in self.main_layout['panel_letter_anchors'] if a['text']==letter)
            x,y=np.asarray(anchor['origin'])*25.4/72
            self.fig.text(x/self.width,1-y/self.height,letter,fontsize=anchor['size'],fontweight='bold',va='baseline')
        else:self.fig.text(x/self.width,1-y/self.height,letter,fontsize=10,fontweight='bold',va='top')
        if title:self.fig.text((x+5)/self.width,1-y/self.height,title,fontsize=8,va='top')

    def save(self,target,description,reuse_figure1=None):
        path=self.out/target;path.parent.mkdir(parents=True,exist_ok=True)
        frozen_pixel_checks=None
        # The manuscript includes F1 at 150 mm and F3 at 146 mm. Protect the
        # minimum readable type size after those existing publication scalings.
        printed_width={'Figure 1':150,'Figure 3':146}.get(self.number,self.width)
        from matplotlib.text import Text
        minimum=6.5*self.width/printed_width
        for text in self.fig.findobj(Text):
            if text.get_text() and text.get_fontsize()<minimum:text.set_fontsize(minimum)
        numerical=self.export_panel_data(description)
        layout_qa=self.geometry() if self.main_layout else None
        self.fig.savefig(path,transparent=reuse_figure1 is not None)
        if reuse_figure1 is not None:
            import pymupdf
            overlay=pymupdf.open(path);base=pymupdf.open(reuse_figure1)
            assert len(base)==1 and abs(base[0].rect.width-overlay[0].rect.width)<.02
            # Remove all old PAS panels before adding their newly generated replacements.
            for rect in [(0,458,221,599.761),(221,313,453.544,599.761)]:
                base[0].add_redact_annot(pymupdf.Rect(rect),fill=(1,1,1))
            base[0].apply_redactions(images=2,graphics=2,text=0)
            base[0].show_pdf_page(base[0].rect,overlay,0)
            merged=path.with_suffix('.merged.pdf');base.save(merged);base.close();overlay.close();merged.replace(path)
            document=pymupdf.open(path);pix=document[0].get_pixmap(dpi=300);pix.save(path.with_suffix('.png'))
            original=pymupdf.open(reuse_figure1);old=original[0].get_pixmap(dpi=300)
            before=np.frombuffer(old.samples,np.uint8).reshape(old.height,old.width,old.n)
            after=np.frombuffer(pix.samples,np.uint8).reshape(pix.height,pix.width,pix.n)
            assert before.shape==after.shape
            scale=300/72
            # Exclude the two-pixel antialiasing boundary of each redaction.
            protected=[(slice(0,int(313*scale)-2),slice(None)),
                (slice(int(313*scale)+2,int(458*scale)-2),slice(0,int(221*scale)-2))]
            frozen_pixel_checks=[]
            for y,x in protected:
                difference=np.max(abs(before[y,x].astype(np.int16)-after[y,x].astype(np.int16)),axis=2)
                # PDF path serialization introduces sparse antialiasing changes
                # (30 pixels in the validated source). Keep a narrow pixel guard,
                # while recording rather than asserting bit-identical rasterization.
                changed=int(np.count_nonzero(difference));maximum=int(difference.max())
                fraction=changed/difference.size
                frozen_pixel_checks.append(dict(changed_pixels=changed,pixels=difference.size,
                    maximum_channel_difference=maximum,changed_fraction=fraction,
                    exactly_equal=not changed,passed=maximum<=8 and fraction<=.0001))
            assert all(c['passed'] for c in frozen_pixel_checks),'Frozen GEX panels changed beyond PDF antialiasing tolerance'
            original.close();document.close()
        else:self.fig.savefig(path.with_suffix('.png'),dpi=150)
        if self.main_layout:
            import pymupdf
            from PIL import ImageDraw
            doc=pymupdf.open(path);pix=doc[0].get_pixmap(dpi=150)
            pix.save(path.with_suffix('.png'))
            from PIL import Image
            tiled=Image.open(path.with_suffix('.png')).convert('RGB');draw=ImageDraw.Draw(tiled)
            for panel,item in self.grid.items():
                x,y,w,h=np.asarray(item['box_pt'])*150/72
                draw.rectangle((x,y,x+w,y+h),outline='#BB3377',width=2)
            tiled.save(path.with_suffix('.tiling.png'));doc.close()
        plt.close(self.fig)
        code={str(Path(m.__file__).resolve()):digest(m.__file__) for m in tuple(sys.modules.values())
            if getattr(m,'__file__',None) and Path(m.__file__).resolve().parent==Path(__file__).resolve().parent and Path(m.__file__).suffix=='.py'}
        code_snapshots=[]
        for source,sha in code.items():
            copy=self.out/'render-code'/(sha[:16]+'__'+Path(source).name);copy.parent.mkdir(exist_ok=True)
            if copy.exists():assert digest(copy)==sha
            else:shutil.copy2(source,copy)
            code_snapshots.append(dict(path=str(copy.relative_to(self.out)),sha256=sha))
        record=dict(figure=self.number,target=target,panels=self.panels,source_files=sorted(self.sources),
            description=description,width_mm=self.width,height_mm=self.height,printed_width_mm=printed_width,
            minimum_new_text_pt_at_print=6.5,sha256=digest(path),
            frozen_GEX_pixel_checks=frozen_pixel_checks,
            layout_geometry=layout_qa,
            run_id=self.config['run_id'],input_fingerprint=self.config['input_fingerprint'],
            numerical_panel_exports=numerical,render_code_sha256=code,render_code_snapshots=code_snapshots,visual_review_complete=False)
        (path.with_suffix('.sources.json')).write_text(json.dumps(record,indent=2)+'\n')
        self.records=[]
        for source in sorted(self.out.glob('Fig/*.sources.json')):
            previous=json.loads(source.read_text())
            assert previous['run_id']==self.config['run_id'] and previous['input_fingerprint']==self.config['input_fingerprint']
            assert digest(self.out/previous['target'])==previous['sha256']
            self.records.append(previous)
        (self.out/'render-progress.json').write_text(json.dumps(dict(run_id=self.config['run_id'],
            input_fingerprint=self.config['input_fingerprint'],figures=self.records,inputs=self.inputs,
            status='partial_until_all_figures_and_visual_review'),indent=2)+'\n')

    def export_panel_data(self,description):
        folder=self.out/'panel-data'/self.number.replace(' ','_');folder.mkdir(parents=True,exist_ok=True)
        exports=[]
        for panel in self.panels:
            arrays={};rows=[]
            for i,ax in enumerate(self.panel_axes.get(panel,[])):
                base=dict(figure=self.number,panel=panel,panel_title=self.panel_titles.get(panel,''),axes=i,x_label=ax.get_xlabel(),y_label=ax.get_ylabel(),
                    x_tick_labels=json.dumps([t.get_text() for t in ax.get_xticklabels()]),
                    y_tick_labels=json.dumps([t.get_text() for t in ax.get_yticklabels()]),
                    source_tables=json.dumps(sorted(self.sources)),protocol=description)
                for kind,artists in [('line',ax.lines),('collection',ax.collections),('image',ax.images),('patch',ax.patches)]:
                    for j,artist in enumerate(artists):
                        key=f'axes{i}_{kind}{j}';keys=[]
                        if kind=='line':data=dict(x=artist.get_xdata(orig=False),y=artist.get_ydata(orig=False))
                        elif kind=='image':data=dict(values=artist.get_array(),extent=artist.get_extent())
                        elif kind=='collection':
                            data=dict(offsets=artist.get_offsets())
                            if artist.get_array() is not None:data['values']=artist.get_array()
                            if hasattr(artist,'get_segments'):
                                # Segment lengths can differ, so preserve each independently.
                                data.update({f'segment{k}':v for k,v in enumerate(artist.get_segments())})
                        elif all(hasattr(artist,k) for k in ['get_x','get_y','get_width','get_height']):
                            data=dict(rectangle=[artist.get_x(),artist.get_y(),artist.get_width(),artist.get_height()])
                        else:continue
                        for name,value in data.items():
                            array=np.asarray(np.ma.filled(value,np.nan))
                            assert array.dtype.kind in 'biufc', (key,name,array.dtype)
                            field=key+'_'+name;arrays[field]=array;keys.append(field)
                        rows.append(dict(base,primitive=key,artist_label=artist.get_label(),array_keys=json.dumps(keys)))
            if not rows:
                rows=[dict(figure=self.number,panel=panel,primitive='retained_frozen_GEX' if self.number=='Figure 1' and panel in ['A','B'] else 'explicit_empty_or_text_panel',
                    source_tables=json.dumps(sorted(self.sources)),protocol=description)]
            tsv=folder/f'{panel}.tsv';pd.DataFrame(rows).to_csv(tsv,sep='\t',index=False)
            values=folder/f'{panel}.npz';np.savez_compressed(values,**arrays)
            exports.extend(dict(path=str(p.relative_to(self.out)),sha256=digest(p)) for p in [tsv,values])
        return exports


def empty(ax,message='No qualifying measurements'):
    ax.text(.5,.5,message,transform=ax.transAxes,ha='center',va='center',fontsize=7,wrap=True)
    ax.set_xticks([]);ax.set_yticks([])


def join_frozen_scores(cells,scores):
    aliases={'WOT_GSE115943':'GSE115943'}
    left=cells.assign(_score_study=cells.study_id.replace(aliases))
    right=scores.rename(columns={'study_id':'_score_study'}).copy()
    right['_score_study']=right._score_study.replace(aliases)
    joined=left.merge(right,on=['metacell_id','_score_study'],how='left',sort=False,
        validate='many_to_one',indicator=True)
    assert joined._merge.eq('both').all(),'A frozen cell has no matching metacell-score record'
    assert joined.cell_id.tolist()==cells.cell_id.tolist(),'Frozen cells were reordered, duplicated or lost'
    return joined.drop(columns=['_score_study','_merge'])


def zero(ax,horizontal=True):
    (ax.axhline if horizontal else ax.axvline)(0,color=GREY,lw=.5,zorder=0)


def matrix(ax,frame,limit,label,cbar=True):
    if frame.empty:empty(ax);return
    im=ax.imshow(np.ma.masked_invalid(frame.to_numpy(float)),cmap=DIVERGE,vmin=-limit,vmax=limit,aspect='auto',interpolation='nearest')
    ax.set_xticks(range(len(frame.columns)),frame.columns);ax.set_yticks(range(len(frame)),frame.index)
    if cbar:ax.figure.colorbar(im,ax=ax,fraction=.035,pad=.03,label=label)
    return im


def ecdf(ax,values,label,color,complement=False,linestyle='-'):
    x=np.sort(np.asarray(values,float));x=x[np.isfinite(x)]
    if not len(x):return
    y=np.arange(1,len(x)+1)/len(x)
    ax.step(x,1-y if complement else y,where='post',color=color,label=f'{label} (n={len(x)})',ls=linestyle)


def medians(ax,frame,group,value,order=None,colors=None):
    order=list(frame[group].drop_duplicates()) if order is None else order
    for i,key in enumerate(order):
        x=frame.loc[frame[group].eq(key),value].dropna()
        if len(x):
            color=colors[i%len(colors)] if colors else INK
            lo,mid,hi=x.quantile([.25,.5,.75]);ax.errorbar(i,mid,yerr=[[mid-lo],[hi-mid]],fmt='o',ms=3,capsize=2,color=color)
        else:ax.plot(i,0,'x',color=GREY,ms=4)
    ax.set_xticks(range(len(order)),order)
    return order


def embedding(ax,cells,value,limit):
    x=cells[value].to_numpy(float);order=np.argsort(np.nan_to_num(np.abs(x),nan=-1))
    ax.scatter(cells.UMAP1,cells.UMAP2,s=.025,c='#DADFE3',rasterized=True,linewidths=0)
    good=order[np.isfinite(x[order])]
    im=ax.scatter(cells.UMAP1.to_numpy()[good],cells.UMAP2.to_numpy()[good],c=x[good],s=.025,
        cmap=DIVERGE,vmin=-limit,vmax=limit,rasterized=True,linewidths=0)
    ax.set_aspect('equal');ax.set_axis_off();return im
