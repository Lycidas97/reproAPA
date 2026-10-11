"""Render the coarse-state overview with reviewed shortcut qualifications."""
from pathlib import Path
import itertools
import json
import math
import subprocess
import numpy as np
import pandas as pd
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

P=Path(__file__).resolve().parents[1]
D=P/'data/state-flow-v1'
OUT=P/'outputs/figures/human-coarse-state-flow-v2'
TABLES=P/'data/state-flow-v2'
W,H=2600,1930
SCALE=300/25.4*72/W
INK,MUTED,FLOW,CHEM='#263D4A','#536575','#41647C','#31776F'


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    TABLES.mkdir(parents=True,exist_ok=True)
    mapping=pd.read_csv(D/'local_to_coarse_state.tsv',sep='\t')
    flows=pd.read_csv(D/'all_coarse_time_flows.tsv',sep='\t')
    selected=flows[flows.display_major & ~flows.coarse_source.eq('U') &
        ~(flows.coarse_source.eq('M') & flows.coarse_target.eq('E'))].copy()
    rare=flows[flows.edge_id.eq('H0045')].copy()
    assert len(rare)==1 and rare.iloc[0].coarse_source=='E' and rare.iloc[0].coarse_target=='K'
    selected['figure_selection']='major'
    rare['figure_selection']='explicit_minor_keratin_branch'
    selected=pd.concat([selected,rare],ignore_index=True)
    original=pd.read_csv(D/'figure_temporal_rows.tsv',sep='\t')
    pd.testing.assert_frame_equal(selected[original.columns].sort_values('edge_id').reset_index(drop=True),
        original.sort_values('edge_id').reset_index(drop=True))
    selected['display_style']=selected.edge_id.map({
        'H0120':'unresolved_interval','H0128':'unresolved_interval',
        'H0255':'endpoint_conditional','H0259':'endpoint_conditional',
        'H0045':'minor'}).fillna('main')
    selected['display_qualification']=selected.display_style.map({
        'main':'研究内相邻观测的主要 OT 联系；不是已验证谱系',
        'unresolved_interval':'中间过程未解析；不是一步转换证据',
        'endpoint_conditional':'目标 P 为 7214/7218；条件运输比例不是成功率',
        'minor':'低于主要显示阈值的次要联系；不是实测转换率'})
    omitted=flows[flows.display_major & ~flows.edge_id.isin(selected.edge_id)].copy()
    omitted['omission_reason']='overview simplification; retained in complete table; not evidence of absence'
    units={
      'S':(80,560,370,220,'成纤维来源 / 来源样','#F7EFE0',
           ['四研究的来源相关群','221739：从 D3 来源样开始']),
      'R':(700,560,400,220,'来源退出与重塑','#EAF0F5',
           ['早中期重塑与可塑性窗口','118258 / 化学路线']),
      'E':(1380,560,420,220,'上皮化 / 混合过渡','#E8F1F7',
           ['上皮、来源 / 目标混合候选','118258 / 150311']),
      'P':(2110,560,410,220,'多能性程序获得','#E2F1EB',
           ['内源支持程序及不同亚型','化学终点 / 成熟参考分列']),
      'N':(2110,170,410,230,'naive 相关候选','#E2F1EB',
           ['150311；与其他目标亚型分列','仅表达归类，未建立时间入边','含 naive / RSeT 培养背景']),
      'M':(660,1170,790,260,'间质 / ECM 保留','#F5EFDF',
           ['来源保留、基质重塑及晚期混合群','保留后续目标联系；未确证转化命运','118258 / 150311 / 化学 / 221739']),
      'K':(1430,950,370,180,'角质细胞样候选','#F0EAF6',
           ['118258 C008 / 150311 C006']),
      'X':(1750,1310,550,280,'内胚层相关分泌程序','#EEEAF7',
           ['化学 C009：TTR / APOA1 / AFP','伴 GATA6；目标核心程序弱','不沿用旧 naive 标签，不确证 XEN']),
      'U':(80,1620,740,170,'混合程序 / 未明确集合','#F6F7F9',
           ['150311：神经 / HAND1 混合候选；221739：低 RNA 群']),
    }
    # Coordinates are schematic; no embedding, time or probability scale is implied.
    paths={
      ('S','R','main'):([(450,675),(700,675)],FLOW),
      ('R','E','main'):([(1100,675),(1380,675)],FLOW),
      ('E','P','main'):([(1800,675),(2110,675)],FLOW),
      ('S','P','unresolved_interval'):([(265,560),(265,445),(2340,445),(2340,560)],FLOW),
      ('S','M','main'):([(265,780),(265,1300),(660,1300)],FLOW),
      ('R','M','main'):([(890,780),(890,1170)],CHEM),
      ('E','M','main'):([(1380,740),(1240,740),(1240,1170)],FLOW),
      ('M','P','unresolved_interval'):([(1450,1220),(1930,1220),(1930,720),(2110,720)],FLOW),
      ('M','P','endpoint_conditional'):([(1450,1270),(2370,1270),(2370,780)],CHEM),
      ('X','P','endpoint_conditional'):([(2300,1360),(2440,1360),(2440,780)],CHEM),
      ('M','X','main'):([(1450,1370),(1750,1370)],CHEM),
      ('X','M','main'):([(1800,1590),(1800,1640),(1360,1640),(1360,1430)],CHEM),
      ('E','K','minor'):([(1590,780),(1590,950)],FLOW),
    }
    assert set(paths)==set(zip(selected.coarse_source,selected.coarse_target,selected.display_style))
    dashes={'main':[],'unresolved_interval':[16,10],
        'endpoint_conditional':[16,7,3,7],'minor':[3,7]}
    pdfmetrics.registerFont(TTFont('Latin','/System/Library/Fonts/Supplemental/Arial.ttf'))
    pdfmetrics.registerFont(TTFont('Chinese','/System/Library/Fonts/STHeiti Light.ttc',subfontIndex=0))
    c=canvas.Canvas(str(OUT/'human_coarse_state_flow.pdf'),pagesize=(W*SCALE,H*SCALE),initialFontName='Latin')
    c.setTitle('人类重编程：整合粗状态流（直连解释修订）');c.scale(SCALE,SCALE)
    bounds=[]

    def text(x,y,s,size=26,color=MUTED,center=False,maxw=None):
        runs=[]
        for ch in s:
            f='Latin' if ord(ch)<128 else 'Chinese'
            assert ord(ch) in pdfmetrics.getFont(f).face.charToGlyph,(ch,f)
            if runs and runs[-1][0]==f:runs[-1]=(f,runs[-1][1]+ch)
            else:runs.append((f,ch))
        width=sum(pdfmetrics.stringWidth(t,f,size) for f,t in runs)
        assert maxw is None or width<=maxw,(s,width,maxw)
        xx=x-width/2 if center else x
        bounds.append(dict(text=s,x=xx,y=y,w=width,h=size*1.1))
        c.setFillColor(color)
        for f,t in runs:
            c.setFont(f,size);c.drawString(xx,H-y-size*.84,t);xx+=pdfmetrics.stringWidth(t,f,size)

    for (_,_,style),(pts,color) in paths.items():
        c.setLineWidth(4);c.setStrokeColor(color);c.setDash(dashes[style])
        path=c.beginPath();path.moveTo(pts[0][0],H-pts[0][1])
        for x,y in pts[1:]:path.lineTo(x,H-y)
        c.drawPath(path);c.setDash([])
        x,y=pts[-1];dx=x-pts[-2][0];dy=y-pts[-2][1];length=math.hypot(dx,dy);dx/=length;dy/=length
        path=c.beginPath();path.moveTo(x,H-y)
        for sign in [1,-1]:path.lineTo(x-15*dx+sign*7*dy,H-(y-15*dy-sign*7*dx))
        path.close();c.setFillColor(color);c.drawPath(path,fill=1,stroke=0)
    for key,(x,y,w,h,title,fill,lines) in units.items():
        c.setFillColor(fill);c.setStrokeColor('#9DAFBC');c.setLineWidth(2)
        c.setDash([10,7] if key in ['R','E','N','K','X','U'] else [])
        c.roundRect(x,H-y-h,w,h,20,fill=1,stroke=1);c.setDash([])
        text(x+w/2,y+23,title,size=33 if key!='S' else 31,color=INK,center=True,maxw=w-24)
        for j,line in enumerate(lines):text(x+w/2,y+80+37*j,line,size=24,center=True,maxw=w-24)
        m=mapping[mapping.coarse_state.eq(key)]
        text(x+18,y+h-32,f'{len(m)} 个原群 · {m.n_cells.sum():,} cells',size=23,maxw=w-36)
    text(80,45,'人类重编程：整合粗状态流',size=52,color=INK)
    text(82,123,'62 个原群 → 8 个粗状态 + 1 个未明确集合；215,114 cells 完整保留。',size=28)
    text(80,220,'按人类表达证据归并；同名粗状态保留研究内亚型。',size=27)
    text(80,265,'时间箭头来自研究内 OT；粗状态不等于必须依次通过的阶段。',size=27)
    text(700,350,'221739 · D9 来源样 → D11 多能性程序',size=25)
    text(700,390,'中间过程未解析；不代表从 D0 一步转化',size=25)
    text(488,610,'118258 / 化学',size=25)
    text(1140,610,'118258',size=25)
    text(1835,605,'118258',size=25)
    text(1835,642,'D12 → D16',size=22)
    text(310,1130,'118258 · D0 → D2',size=23)
    text(310,1170,'118258 · D2 → D8',size=23)
    text(310,1210,'221739 · D9 → D11',size=23)
    text(930,970,'化学 0618',size=24,color=CHEM)
    text(930,1010,'II d20 → III d4',size=23,color=CHEM)
    text(1270,865,'118258',size=23)
    text(1270,902,'D12 → D16',size=23)
    text(1620,820,'118258 · D12 → D16',size=22)
    text(1620,856,'3.1% 去向（次要联系）',size=22)
    text(1970,845,'221739 · D11 → D13',size=24)
    text(1970,885,'ECM → 目标程序',size=24)
    text(1970,925,'中间过程未解析',size=24)
    text(1970,1010,'化学 0618 · IV d10',size=23,color=CHEM)
    text(1970,1050,'→ endpoint：两类条件联系',size=23,color=CHEM)
    text(1970,1110,'终点 P：7,214 / 7,218',size=23,color=CHEM)
    text(1970,1150,'高比例受终点组成限制',size=23,color=CHEM)
    text(1480,1285,'0618 · IV d2 → IV d4',size=23,color=CHEM)
    text(1420,1680,'0618 · III d12 → IV d1；IV d4 → IV d10',size=23,color=CHEM)
    text(1420,1720,'不同采样段的联系；不解释为可逆命运或细胞周期',size=23)
    text(80,1820,'箭头：实线为主要联系；长虚线为中间过程未解析；点划线为化学终点条件联系；点线为次要联系。',size=24)
    text(80,1860,'所有箭头均为 OT 推断，不代表已验证谱系或转换率；仅沿同一研究 / 供者读取。完整时间流见附表。',size=24)
    collisions=[]
    for (source,target,style),(pts,_) in paths.items():
        pair=(source,target)
        samples=np.concatenate([np.linspace(a,b,1001) for a,b in zip(pts[:-1],pts[1:])])
        for key,(x,y,w,h,*_) in units.items():
            if key not in pair and ((samples[:,0]>x-2)&(samples[:,0]<x+w+2)&(samples[:,1]>y-2)&(samples[:,1]<y+h+2)).any():
                collisions.append([pair,'node',key])
        for b in bounds:
            if ((samples[:,0]>b['x'])&(samples[:,0]<b['x']+b['w'])&(samples[:,1]>b['y'])&(samples[:,1]<b['y']+b['h'])).any():
                collisions.append([pair,'text',b['text']])
    for a,b in itertools.combinations(bounds,2):
        if a['x']<b['x']+b['w'] and b['x']<a['x']+a['w'] and a['y']<b['y']+b['h'] and b['y']<a['y']+a['h']:
            collisions.append(['text_overlap',a['text'],b['text']])
    assert not collisions,collisions
    c.showPage();c.save()
    selected.to_csv(TABLES/'figure_temporal_rows.tsv',sep='\t',index=False)
    omitted.to_csv(TABLES/'major_rows_not_drawn.tsv',sep='\t',index=False)
    bundles=[]
    for (a,b,style),g in selected.groupby(['coarse_source','coarse_target','display_style']):
        bundles.append(dict(source=a,target=b,edge_ids=';'.join(g.edge_id),routes=';'.join(sorted(g.route_id.unique())),
                            display_style=style,across_pair_probability_aggregation='none'))
    pd.DataFrame(bundles).to_csv(TABLES/'figure_edge_bundles.tsv',sep='\t',index=False)
    poppler=Path('/Users/lycidas/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/poppler/bin')
    subprocess.run([str(poppler/'pdftoppm'),'-r','300','-singlefile','-png',str(OUT/'human_coarse_state_flow.pdf'),str(OUT/'human_coarse_state_flow')],check=True)
    subprocess.run([str(poppler/'pdftocairo'),'-svg',str(OUT/'human_coarse_state_flow.pdf'),str(OUT/'human_coarse_state_flow.svg')],check=True)
    review=dict(status='GEOMETRY_PASS_PENDING_VISUAL',original_states=62,cells=215114,coarse_states=8,
        unresolved_collections=1,displayed_time_rows=len(selected),arrow_bundles=len(bundles),
        major_rows_not_drawn=omitted.edge_id.tolist(),geometry_collisions=collisions,
        min_font_pt=22*SCALE,width_mm=300,font='Arial + STHeiti Light',original_coordinates_modified=False,
        original_selected_rows_unchanged=True,styles=dashes)
    (OUT/'review.json').write_text(json.dumps(review,ensure_ascii=False,indent=2)+'\n')
    (OUT/'text_bounds.json').write_text(json.dumps(bounds,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(review,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
