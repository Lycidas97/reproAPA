"""Render S21--S24 from frozen results; no statistical analysis is rerun."""
from pathlib import Path
from collections import Counter
import json
import subprocess
import tempfile
import hashlib
import shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle
from matplotlib.transforms import Bbox
from matplotlib.text import Text
from PIL import Image
import pymupdf

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / 'figures'
STYLE = json.loads((ROOT / 'notes/2026-09-13-pas-reliability-supplement/visual-refinement/figure1_style.json').read_text())
for key in ['font_regular', 'font_bold']:
    font_manager.fontManager.addfont(STYLE[key])
FONT = font_manager.FontProperties(fname=STYLE['font_regular']).get_name()
INK, GRAY, LIGHT = [STYLE[k] for k in ['ink', 'muted', 'light']]
BLUE, ORANGE = STYLE['quantitative'].values()
PURPLE = STYLE['schematic']['purple']
COLORS = dict(DWUI=BLUE, RWUI=ORANGE, PDUI=GRAY, PPUI='#328C98', MPRO=PURPLE)
LABELS = {'GSE118258': 'Xing2020', 'GSE150311': 'Liu2020', '178325_0618': 'Guan2022–0618',
          '178325_0809': 'Guan2022–0809', 'GSE221739': 'Panariello2023',
          'CIP': 'Guo2019–CIP', 'WOT': 'Schiebinger2019–2i'}
METRICS = ['DWUI', 'RWUI', 'PDUI', 'PPUI', 'MPRO']
NATIVE_WIDTH_MM, PRINT_WIDTH_MM = 160, 150
PRINT_SCALE = PRINT_WIDTH_MM / NATIVE_WIDTH_MM
plt.rcParams.update({'font.family': FONT, 'font.size': 8, 'text.color': INK,
    'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK,
    'axes.labelsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
    'axes.edgecolor': '#AAB4BB', 'axes.linewidth': .5, 'axes.spines.top': False,
    'axes.spines.right': False, 'pdf.fonttype': 42, 'ps.fonttype': 42,
    'svg.fonttype': 'none', 'legend.fontsize': 8, 'legend.frameon': False,
    'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'lines.linewidth': .8})
QA = {}


def read(name):
    return pd.read_csv(HERE / name, sep='\t')


class Page:
    """Physical millimetres, top-down placement and complete-panel QA."""
    def __init__(self, name, height):
        self.name, self.height = name, height
        self.fig = plt.figure(figsize=(NATIVE_WIDTH_MM / 25.4, height / 25.4), dpi=150)
        self.groups = {}

    def text(self, panel, x, y, s, size=8, **kwargs):
        t = self.fig.text(x / NATIVE_WIDTH_MM, 1 - y / self.height, s, fontsize=size,
                          va='top', linespacing=1.2, **kwargs)
        self.groups.setdefault(panel, []).append(t)
        return t

    def label(self, panel, x, y):
        self.text(panel, x, y, panel, 11, fontweight='bold')

    def axes(self, panel, x, y, w, h):
        ax = self.fig.add_axes([x / NATIVE_WIDTH_MM, 1 - (y + h) / self.height,
                                w / NATIVE_WIDTH_MM, h / self.height])
        self.groups.setdefault(panel, []).append(ax)
        return ax

    def export(self, pairs):
        self.fig.canvas.draw()
        renderer = self.fig.canvas.get_renderer()
        texts = [t for t in self.fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
        text_boxes = [(t.get_text(), t.get_window_extent(renderer)) for t in texts]
        collisions = []
        clearances = []
        for i,(a,aa) in enumerate(text_boxes):
            for b,bb in text_boxes[i+1:]:
                dx = max(aa.x0-bb.x1, bb.x0-aa.x1)
                dy = max(aa.y0-bb.y1, bb.y0-aa.y1)
                clearances.append(np.hypot(max(dx,0),max(dy,0))/self.fig.dpi*72)
                if dx < 0 and dy < 0:
                    collisions.append([a,b])
        assert not collisions, (self.name, collisions)
        assert min(clearances) * PRINT_SCALE >= 2, (self.name, min(clearances) * PRINT_SCALE)
        boxes = {}
        for key, objects in self.groups.items():
            boxes[key] = Bbox.union([o.get_tightbbox(renderer) if hasattr(o, 'get_tightbbox')
                                    else o.get_window_extent(renderer) for o in objects])
        gutters = {}
        for a, b, direction in pairs:
            gap = boxes[b].x0 - boxes[a].x1 if direction == 'horizontal' else boxes[a].y0 - boxes[b].y1
            gutters[a + '/' + b] = gap / self.fig.dpi * 72
            assert gutters[a + '/' + b] * PRINT_SCALE >= 3, (self.name, a, b, gutters)
        for ext in ['pdf', 'svg']:
            self.fig.savefig(OUT / f'{self.name}.{ext}')
        pdf = OUT / f'{self.name}.pdf'
        fonts = subprocess.run(['pdffonts', str(pdf)], capture_output=True, text=True, check=True)
        if 'Mismatch between font type and embedded font file' in fonts.stderr:
            original_text = subprocess.check_output(['pdftotext', str(pdf), '-'], stderr=subprocess.DEVNULL)
            with tempfile.TemporaryDirectory() as tmp:
                normalized = Path(tmp) / pdf.name
                subprocess.run(['gs', '-q', '-dBATCH', '-dNOPAUSE', '-sDEVICE=pdfwrite',
                    '-dCompatibilityLevel=1.7', '-dEmbedAllFonts=true', '-dSubsetFonts=true',
                    '-dAutoRotatePages=/None', '-dDownsampleColorImages=false',
                    '-dDownsampleGrayImages=false', '-sOutputFile='+str(normalized), str(pdf)],
                    check=True, capture_output=True)
                assert Counter(original_text.split()) == Counter(subprocess.check_output(['pdftotext', str(normalized), '-']).split())
                normalized.replace(pdf)
            fonts = subprocess.run(['pdffonts', str(pdf)], capture_output=True, text=True, check=True)
        assert not fonts.stderr.strip(), fonts.stderr
        subprocess.run(['pdftoppm', '-r', '300', '-singlefile', '-png', str(pdf),
                        str(OUT / self.name)], check=True, capture_output=True)
        pixels = np.asarray(Image.open(OUT / f'{self.name}.png').convert('RGB'))
        yy, xx = np.where(np.any(pixels < 250, axis=2))
        margins = dict(left=int(xx.min()), right=int(pixels.shape[1] - 1 - xx.max()),
                       top=int(yy.min()), bottom=int(pixels.shape[0] - 1 - yy.max()))
        assert min(margins.values()) >= 12, (self.name, margins)
        # Raster-visible unions in the renderer-established complete panel regions.
        visible = {}
        for key, box in boxes.items():
            scale = 300 / self.fig.dpi
            x0, x1 = max(0, int(box.x0 * scale) - 3), min(pixels.shape[1], int(box.x1 * scale) + 4)
            y0, y1 = max(0, int(pixels.shape[0] - box.y1 * scale) - 3), min(pixels.shape[0], int(pixels.shape[0] - box.y0 * scale) + 4)
            y, x = np.where(np.any(pixels[y0:y1, x0:x1] < 250, axis=2))
            visible[key] = [int(x.min() + x0), int(y.min() + y0), int(x.max() + x0), int(y.max() + y0)]
        raster_gaps = {a+'/'+b: (visible[b][0]-visible[a][2] if d == 'horizontal'
                        else visible[b][1]-visible[a][3]) for a,b,d in pairs}
        assert min(raster_gaps.values()) > 0
        (OUT / f'{self.name}.fonts.txt').write_text(fonts.stdout)
        assert all('yes yes' in line for line in fonts.stdout.splitlines()[2:]), fonts.stdout
        with pymupdf.open(pdf) as doc:
            spans = [s for block in doc[0].get_text('dict')['blocks'] if 'lines' in block
                     for line in block['lines'] for s in line['spans'] if s['text'].strip()]
        native_sizes = sorted({round(s['size'], 3) for s in spans})
        assert native_sizes == [8, 11], (self.name, native_sizes)
        assert all('NimbusSanL' in s['font'] for s in spans)
        QA[self.name] = dict(font=FONT, width_mm=NATIVE_WIDTH_MM, height_mm=self.height,
            print_width_mm=PRINT_WIDTH_MM, print_scale=PRINT_SCALE,
            native_font_sizes_pt=native_sizes, minimum_text_pt=min(native_sizes) * PRINT_SCALE,
            panel_text_pt=11 * PRINT_SCALE, axis_text_pt=8 * PRINT_SCALE, raster_margins_px=margins,
            complete_panel_boxes_px=visible, complete_panel_gutters_pt=gutters,
            raster_gutters_px=raster_gaps, all_fonts_embedded=True, text_collisions=collisions,
            minimum_text_clearance_pt=min(clearances) * PRINT_SCALE,
            complete_panel_gutters_final_pt={k: v * PRINT_SCALE for k, v in gutters.items()},
            panel_titles=False)
        plt.close(self.fig)


def overview():
    p = Page('multimetric_overview', 133)
    p.label('A', 4, 4)
    ax = p.axes('A', 30, 16, 46, 42)
    data = read('global_metric_summary.tsv')
    data = data[data.contrast.eq('P_S') & data.subset.eq('all')]
    order = ['GSE118258', '178325_0618', '178325_0809', 'CIP', 'WOT']
    for j, (metric, marker) in enumerate(zip(METRICS[:4], ['o', 's', '^', 'D'])):
        q = data[data.metric.eq(metric)].set_index('route').loc[order]
        ax.scatter(q['mean'] * (-1 if metric == 'PPUI' else 1), np.arange(5) + (j-1.5)*.14,
                   s=13, marker=marker, color=COLORS[metric], label='−PPUI' if metric == 'PPUI' else metric, zorder=3)
    ax.axvline(0, color='#AAB4BB', lw=.6)
    ax.set(yticks=np.arange(5), yticklabels=[LABELS[x] for x in order], ylim=(4.5, -.5),
           xlim=(-.017, .043), xticks=[-.01, .01, .03], xlabel='Mean index change')
    ax.tick_params(axis='y', length=0)
    ax.legend(ncol=2, loc='lower left', bbox_to_anchor=(-.06, 1.01), columnspacing=.8,
              handletextpad=.25, borderaxespad=0, labelspacing=.3)
    p.label('B', 86, 4)
    ax = p.axes('B', 116, 16, 39, 42)
    data = read('hidden_redistribution_summary.tsv')
    data = data[data.subset.eq('PAS_ge3') & data.MPRO_threshold.eq(.1)]
    order = list(LABELS)
    for contrast, marker, offset in [('P_S','o',-.12), ('P_M','s',.12)]:
        q = data[data.contrast.eq(contrast)].set_index('route')
        for i, route in enumerate(order):
            if route in q.index:
                r = q.loc[route]
                ax.scatter(r.hidden_fraction*100, i+offset, marker=marker, s=16,
                           facecolor=PURPLE if contrast=='P_S' else 'white', edgecolor=PURPLE, lw=.8)
    ax.set(yticks=np.arange(7), yticklabels=[LABELS[x] for x in order], ylim=(6.5,-.5),
           xlim=(0,40), xticks=[0,10,20,30,40], xlabel='Candidate fraction (%)')
    ax.tick_params(axis='y', length=0)
    for marker, label, filled in [('o','P–S',True),('s','P–M',False)]:
        ax.scatter([],[],marker=marker,s=16,facecolor=PURPLE if filled else 'white',edgecolor=PURPLE,label=label)
    ax.legend(ncol=2,loc='lower left',bbox_to_anchor=(0,1.02),handletextpad=.2,columnspacing=.8,borderaxespad=0)
    axes_table = read('frozen_PAS_axes.tsv.gz')
    for letter, x, row in zip(['C','D'], [4,86], read('composition_examples.tsv').itertuples()):
        p.label(letter, x, 73)
        ax = p.axes(letter, x+14, 85, 55, 34)
        ids = json.loads(row.PAS_ids)
        weights = axes_table[axes_table.study_id.eq(row.study_id)].set_index('atlas_pas_id').loc[ids].position.to_numpy()
        a, b = np.array(json.loads(row.target_usage)), np.array(json.loads(row.reference_usage))
        for values, marker, color, label in [(b,'o',GRAY,'S reference'),(a,'s',BLUE,'P target')]:
            ax.vlines(weights,0,values,color=color,lw=.7,alpha=.7)
            ax.scatter(weights,values,marker=marker,s=18,facecolor='white' if marker=='o' else color,
                       edgecolor=color,lw=.8,label=label,zorder=3)
        ax.set(xlim=(-.045,1.045),ylim=(-.025,1.02),xticks=[0,.5,1],yticks=[0,.5,1],
               xlabel='Relative PAS position',ylabel='PAS proportion')
        ax.legend(loc='lower left',bbox_to_anchor=(-.06,1.03),ncol=2,handletextpad=.2,columnspacing=.8,borderaxespad=0)
    p.export([('A','B','horizontal'),('C','D','horizontal'),('A','C','vertical'),('B','D','vertical')])


def set_figure():
    p = Page('metric_set_robustness', 122)
    cmap = LinearSegmentedColormap.from_list('enrichment', ['#F3F6F7','#BED4E3',BLUE])
    norm = Normalize(1,4.1)
    settings = [('human','P_S'),('mouse','P_S'),('human','P_M'),('mouse','P_M')]
    pairs = [('GSE118258','GSE103221'),('GSE118258','GSE115943'),('GSE178325','GSE103221'),('GSE178325','GSE115943')]
    rows = [['Human P–S','Mouse P–S','Human P–M','Mouse P–M'],
            ['Xing / Guo CIP','Xing / Schiebinger','Guan / Guo CIP','Guan / Schiebinger']]
    all_cells = []
    for letter, y, source, keys, labels in zip(['A','B'],[4,54],['nulls/recurrence_summary.tsv','nulls/sharing_summary.tsv'],[settings,pairs],rows):
        p.label(letter,5.3,y)
        data=read(source)
        if letter=='B': data=data[data.test.isin(['same_direction','shared_unsigned'])]
        ax=p.axes(letter,47.3,y+12,108,32)
        ax.set(xlim=(-.5,5.5),ylim=(3.5,-.5),yticks=range(4),yticklabels=labels,
               xticks=range(6),xticklabels=['0.02','0.05','0.02','0.05','0.10','0.20'])
        ax.xaxis.tick_top(); ax.tick_params(length=0,pad=3)
        for sp in ax.spines.values():sp.set_visible(False)
        for j,metric in enumerate(['WPP','RWUI','MPRO']):
            ax.text(j*2+.5,-1.12,'DWUI' if metric=='WPP' else metric,ha='center',va='center',fontsize=8)
            for k,threshold in enumerate([.10,.20] if metric=='MPRO' else [.02,.05]):
                q=data[data.metric.eq(metric) & np.isclose(data.threshold,threshold)]
                index=['species','contrast'] if letter=='A' else ['human_study','mouse_study']
                q=q.set_index(index).loc[keys]
                for i,r in enumerate(q.itertuples()):
                    x=j*2+k
                    supported=(r.q_all_24 if letter=='A' else r.q_joint_primary_24)<.05
                    ax.add_patch(Rectangle((x-.46,i-.44),.92,.88,facecolor=cmap(norm(r.enrichment)),
                                          edgecolor=INK if supported else 'white',lw=.8 if supported else .3))
                    ax.text(x,i,f'{r.enrichment:.2f}',ha='center',va='center',fontsize=8,
                            color='white' if r.enrichment>3.4 else INK)
                    all_cells.append(dict(panel=letter,metric=metric,threshold=threshold,row=labels[i],
                                          enrichment=r.enrichment,q=r.q_all_24 if letter=='A' else r.q_joint_primary_24))
    p.text('B',47.3,103,'Outline: q < 0.05 (24 tests per panel)')
    cax=p.axes('B',47.3,110,50,2)
    plt.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation='horizontal',ticks=[1,2,3,4])
    cax.tick_params(labelsize=8,pad=1)
    p.text('B',103.3,109,'Observed / null mean')
    pd.DataFrame(all_cells).to_csv(OUT/'metric_set_robustness.plotted.tsv',sep='\t',index=False)
    p.export([('A','B','vertical')])


def inference_figure():
    p=Page('replicate_supported_inference',100)
    data=read('metric_inference/gene_metric_tests.tsv.gz')
    p.label('A',5.9,4)
    ax=p.axes('A',22.9,13,50,42)
    counts=data.groupby('metric').joint_q05.sum().reindex(METRICS)
    ax.barh(range(5),counts,color=[COLORS[m] for m in METRICS],height=.58)
    for i,n in enumerate(counts):ax.text(n+5,i,str(n),va='center',fontsize=8)
    ax.set(yticks=range(5),yticklabels=METRICS,ylim=(4.6,-.6),xlim=(0,230),xticks=[0,100,200],xlabel='Genes with joint q < 0.05')
    ax.tick_params(axis='y',length=0)
    p.label('B',87.9,4)
    ax=p.axes('B',101.9,13,52,42)
    q=data.pivot(index='gene_id',columns='metric',values='q_BH_joint')
    assert q.shape == (268,5) and not q.isna().any().any()
    x,y=-np.log10(q.DWUI),-np.log10(q.MPRO)
    cut=-np.log10(.05)
    categories=[((q.DWUI>=.05)&(q.MPRO>=.05),'#BCC4C8'),((q.DWUI<.05)&(q.MPRO<.05),BLUE),
                ((q.DWUI<.05)&(q.MPRO>=.05),ORANGE),((q.DWUI>=.05)&(q.MPRO<.05),PURPLE)]
    assert [int(mask.sum()) for mask,_ in categories] == [68,171,15,14]
    for mask,color in categories:ax.scatter(x[mask],y[mask],s=10,color=color,alpha=.75,edgecolors='white',linewidths=.15)
    for value,axis in [(cut,'x'),(cut,'y')]:
        (ax.axvline if axis=='x' else ax.axhline)(value,color=GRAY,lw=.65,ls='--')
    ax.set(xlim=(-.12,4.15),ylim=(-.12,4.15),xticks=[0,1,2,3,4],yticks=[0,1,2,3,4],
           xlabel='DWUI −log10(q)',ylabel='MPRO −log10(q)')
    for xx,yy,txt,color in [(.12,3.8,'14',PURPLE),(1.55,3.8,'171',BLUE),
                            (1.55,.8,'15',ORANGE),(.1,.8,'68',GRAY)]:
        ax.text(xx,yy,txt,fontsize=8,color=color,va='top',linespacing=1.15,
                bbox=dict(facecolor='white',edgecolor='none',alpha=.88,pad=.5))
    p.label('C',5.9,72)
    ax=p.axes('C',67.9,74,89,22)
    candidate = read('metric_inference/candidate_detection_summary.tsv').set_index('metric').loc['DWUI']
    values=[int(candidate[k]) for k in ['original_candidates','inference_covered_genes',
        'significant_original_direction','joint_q05_original_direction_effect_ge_0p02']]
    assert values == [216,62,62,60]
    labels=['Original mouse candidates','Eligible for well-level tests','Same-direction DWUI q < 0.05','Also |ΔDWUI| ≥ 0.02']
    ax.barh(range(4),values,color=[LIGHT,GRAY,BLUE,BLUE],height=.57)
    for i,n in enumerate(values):ax.text(n+4,i,str(n),va='center',fontsize=8)
    ax.set(yticks=range(4),yticklabels=labels,ylim=(3.5,-.5),xlim=(0,240),xticks=[])
    ax.tick_params(axis='y',length=0)
    for sp in ax.spines.values():sp.set_visible(False)
    ax.text(120,1.05,'154 not covered',fontsize=8,color=GRAY,va='center')
    p.export([('A','B','horizontal'),('A','C','vertical'),('B','C','vertical')])


def diagnostics():
    p=Page('inference_design_limits',83)
    p.label('A',5.2,4)
    ax=p.axes('A',19.2,12,56,46)
    data=read('inference/source_endpoint_WPP_tests.tsv.gz')
    sig=data.equal_MC_WPP_q_BH_combined.lt(.05)
    zero=data.equal_MC_WPP_source_zero_observed_variance
    assert [int((sig & ~zero).sum()),int((sig & zero).sum())] == [115,37]
    x=data.equal_MC_WPP_estimate; y=-np.log10(data.equal_MC_WPP_q_BH_combined)
    ax.scatter(x[~sig],y[~sig],s=7,color='#C5CCCF',alpha=.6,edgecolors='none')
    ax.scatter(x[sig & ~zero],y[sig & ~zero],s=10,color=BLUE,alpha=.8,edgecolors='none',label='Supported; source variance > 0 (115)')
    ax.scatter(x[sig & zero],y[sig & zero],s=15,marker='o',facecolors='none',edgecolors=ORANGE,lw=.6,label='Supported; source variance = 0 (37)')
    ax.axhline(-np.log10(.05),ls='--',color=GRAY,lw=.6)
    ax.axvline(0,color='#AAB4BB',lw=.5)
    ax.set(xlabel='Date-balanced ΔDWUI',ylabel='−log10(joint q)',xticks=[-.6,0,.6],yticks=[0,2,4,6])
    ax.legend(loc='upper left',bbox_to_anchor=(-.2,-.25),borderaxespad=0,handletextpad=.3,labelspacing=.4)
    p.label('B',87.2,4)
    ax=p.axes('B',103.2,12,51,46)
    d=read('pool_inference/gene_metric_tests.tsv.gz').query('metric == "DWUI"')
    minimum=np.sort(d.minimum_two_sided_p.to_numpy());rank=np.arange(1,len(d)+1)
    ax.step(rank,minimum,where='post',color=PURPLE,lw=1.3,label='Minimum attainable P')
    ax.plot(rank,.05*rank/len(d),color=INK,lw=1,ls='--',label='BH rejection boundary')
    best=float(np.min(minimum*len(d)/rank))
    assert np.isclose(best,.06312751004016064)
    ax.text(35,.137,f'Min. attainable q = {best:.4f}',fontsize=8,va='top')
    ax.set(xlim=(0,1050),ylim=(0,.15),xticks=[0,500,1000],yticks=[0,.05,.1,.15],
           xlabel='Gene rank (attainable P)',ylabel='P value')
    ax.legend(loc='upper left',bbox_to_anchor=(-.23,-.25),borderaxespad=0,handletextpad=.5,labelspacing=.4)
    p.export([('A','B','horizontal')])


def main():
    OUT.mkdir(exist_ok=True)
    overview();set_figure();inference_figure();diagnostics()
    (HERE/'figure_QA.json').write_text(json.dumps(QA,indent=2)+'\n')
    # Promote the verified PDFs through the same asset manifest as S1--S20.
    manifest=ROOT/'Fig/figure-manifest.tsv'
    mapping=[('multimetric_overview','Position summaries and PAS redistribution'),
        ('metric_set_robustness','Metric dependence of recurrence and cross-species sharing'),
        ('replicate_supported_inference','Matched five-metric independent-well inference'),
        ('inference_design_limits','Source replication and biological-unit resolution')]
    lines=[line for line in manifest.read_text().splitlines()
           if line.split('\t')[0] not in {f'Figure S{n}' for n in range(21,25)}]
    for n,(name,title) in enumerate(mapping,21):
        src=OUT/f'{name}.pdf'; dest=ROOT/f'Fig/SuppFig.S{n}.pdf'
        shutil.copy2(src,dest)
        lines.append('\t'.join([f'Figure S{n}',dest.name,'../'+str(src.relative_to(ROOT)),
                                 hashlib.sha256(src.read_bytes()).hexdigest(),title]))
    manifest.write_text('\n'.join(lines)+'\n')
    print(json.dumps(QA,indent=2))


if __name__=='__main__':
    main()
