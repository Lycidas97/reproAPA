"""Native recurrence, contributor matching and independent-unit evidence displays."""
import numpy as np
import pandas as pd

from render_common import BLUE,ORANGE,INK,GREY,LABELS,empty,zero


def primary(frame,metric='WPP',threshold=None):
    f=frame[frame.metric.eq(metric)]
    if 'series_matching' in f:f=f[f.series_matching.astype(str).str.lower().eq('true')]
    if threshold is not None:f=f[np.isclose(f.threshold,threshold)]
    return f.copy()


def null_forest(ax,frame,labels,observed='observed',q='q_metric_8',ratio=False):
    assert {observed,q,'null_mean','null_p025','null_p975'}<=set(frame),frame.columns.tolist()
    ax.set_yticks(range(len(frame)),labels);ax.invert_yaxis()
    for i,row in enumerate(frame.to_dict('records')):
        if not np.isfinite(row.get('null_mean',np.nan)):
            ax.text(.01,i,'Unavailable',transform=ax.get_yaxis_transform(),fontsize=6.5,color=GREY);continue
        if ratio and row['null_mean']==0:
            ax.text(.01,i,f"O={row[observed]:g}; null mean=0; ratio undefined",transform=ax.get_yaxis_transform(),fontsize=6.5,color=INK);continue
        assert row['null_mean']>=0
        scale=row['null_mean'] if ratio else 1.
        ax.plot([row['null_p025']/scale,row['null_p975']/scale],[i,i],color=GREY,lw=2.5)
        ax.plot(row['null_mean']/scale,i,'|',color=INK,ms=7)
        color=BLUE if row.get('species','human')=='human' else ORANGE
        ax.plot(row[observed]/scale,i,'o',ms=3.5,mec=color,mfc=color if row[q]<.05 else 'white')
    ax.set_ylim(len(frame)-.5,-.5)
    ax.set_xlabel('Observed / null mean' if ratio else 'Recurrent genes')


def ladder_frame(r):
    frame=r.read('recurrence_controls/recurrence_ladder.tsv').rename(columns={'expected':'null_mean'})
    with np.load(r.artifact('recurrence_controls/recurrence_draws.npz')) as arrays:
        for i,row in frame.iterrows():
            key='|'.join(map(str,[row.species,row.contrast,row.threshold,row.model]));values=arrays[key]
            assert np.isclose(values.mean(),row.null_mean)
            frame.loc[i,['null_p025','null_p975']]=np.quantile(values,[.025,.975])
    return frame


def grouped_null(ax,frame,group,q,ratio=False):
    assert {group,q,'null_mean','null_p025','null_p975'}<=set(frame)
    order=list(frame[group].drop_duplicates());positions=[]
    contexts=[('human','P_S'),('human','P_M'),('mouse','P_S'),('mouse','P_M')]
    for row in frame.to_dict('records'):
        j=contexts.index((row['species'],row['contrast']));y=order.index(row[group])+(j-1.5)*.19;positions.append(y)
        scale=row['null_mean'] if ratio else 1.
        if not np.isfinite(scale) or not np.isfinite(row['null_mean']):
            ax.plot(.015,y,'x',transform=ax.get_yaxis_transform(),color=GREY,ms=3);continue
        if ratio and scale==0:
            ax.text(.01,y,f"O={row['observed']:g}; null=0",transform=ax.get_yaxis_transform(),fontsize=6.5,color=INK);continue
        assert scale>0
        color=BLUE if row['species']=='human' else ORANGE
        ax.plot([row['null_p025']/scale,row['null_p975']/scale],[y,y],color=GREY,lw=1.4)
        ax.plot(row['observed']/scale,y,'o' if row['contrast']=='P_S' else 's',ms=2.7,mec=color,mfc=color if row[q]<.05 else 'white')
    ax.set_yticks(range(len(order)),order);ax.set_ylim(len(order)-.5,-.5);ax.set_xlabel('Observed / null mean' if ratio else 'Recurrent genes')
    return positions


def matched_flow(members,matched,threshold=.02):
    roster=primary(members,threshold=threshold);roster=roster[roster.contrast.eq('P_M')]
    matched=matched[matched.setting.eq('supported50_p08')]
    detail=[];flow=[]
    for sp in ['human','mouse']:
        routes=['GSE118258'] if sp=='human' else ['CIP','WOT'];required=len(routes)
        f=roster[roster.species.eq(sp)].drop_duplicates('gene_id')
        for row in f.itertuples():
            q=matched[matched.species.eq(sp)&matched.gene_id.eq(row.gene_id)&matched.route.isin(routes)]
            votes=[]
            for sid,g in q.groupby('study_id'):
                votes.append(dict(study_id=sid,**{co:1 if g['DWUI_'+co].ge(threshold).all() else -1 if g['DWUI_'+co].le(-threshold).all() else 0 for co in ['PM','PS']}))
            sign=1 if row.direction=='positive' else -1
            complete=len(votes)==required
            retained=complete and all(v['PM']==sign for v in votes)
            ps=complete and all(v['PS']==sign for v in votes)
            detail.append(dict(species=sp,gene_id=row.gene_id,matched_studies=len(votes),required_studies=required,complete=complete,
                retained_PM=retained,additional_PS=retained and ps,original_PM_direction=sign))
        d=pd.DataFrame([x for x in detail if x['species']==sp])
        counts=[len(f),int(d.complete.sum()) if len(d) else 0,int(d.retained_PM.sum()) if len(d) else 0,int(d.additional_PS.sum()) if len(d) else 0]
        flow.extend(dict(species=sp,stage=label,n=n) for label,n in zip(['P–M recurrent','Matched P/M/S','P–M retained','Also aligned P–S'],counts))
    return pd.DataFrame(flow),pd.DataFrame(detail)


def example_panel(r,ax,symbol):
    means=r.read('matched_pms/matched_state_vectors.tsv.gz')
    f=means[means.setting.eq('supported50_p08')&means.gene_symbol.str.upper().eq(symbol.upper())]
    if not len(f):empty(ax,f'{symbol}: no qualified matched P/M/S');return
    mc=r.read('examples/example_metacell_scores.tsv.gz')
    pairs=r.read('matched_pms/weighted_triplets.tsv.gz')
    for j,route in enumerate(['CIP','WOT']):
        mean=f[f.route.eq(route)];p=pairs[pairs.setting.eq('supported50_p08')&pairs.route.eq(route)&pairs.gene_symbol.str.upper().eq(symbol.upper())]
        if not len(mean):
            ax.plot([],[],ls='',marker='x',color=GREY,label=route+': unavailable');continue
        assert len(mean)==1
        allowed=set(p.library_key_a)|set(p.library_key_b)
        q=mc[mc.setting.eq('supported50_p08')&mc.gene_symbol.str.upper().eq(symbol.upper())&mc.route.eq(route)&mc.library_key.isin(allowed)&mc.qualified]
        for i,state in enumerate('SPM'):
            g=q[q.node.eq(state)].sort_values('metacell_id');x=i+(j-.5)*.25
            jitter=np.linspace(-.06,.06,len(g))
            ax.scatter(x+jitter,g.DWUI,s=2,alpha=.3,color=[BLUE,ORANGE][j],rasterized=True)
            lib=g.groupby('library_key').DWUI.mean();ax.scatter(np.full(len(lib),x),lib,s=10,facecolors='white',edgecolors=[BLUE,ORANGE][j],linewidths=.6)
            ax.scatter(x,mean[state+'_DWUI'].iloc[0],s=17,marker='D',color=[BLUE,ORANGE][j],zorder=4,label=route if i==0 else None)
        ax.plot(np.arange(3)+(j-.5)*.25,[mean[s+'_DWUI'].iloc[0] for s in 'SPM'],color=[BLUE,ORANGE][j],lw=.7)
    ax.set_xticks(range(3),list('SPM'));ax.set_ylabel('DWUI');ax.set_title(symbol,fontsize=8,fontstyle='italic');ax.legend(frameon=False,fontsize=6.5)


def figure3(r):
    import json
    from matplotlib.patches import Rectangle
    from matplotlib.lines import Line2D
    from render_common import HOUSE_INK,STATE_COLORS,STUDY_STYLE
    r.main(3);rec=primary(r.read('recurrence/recurrence_summary.tsv'));members=primary(r.read('recurrence/recurrence_members.tsv.gz'),threshold=.02)
    contexts=[('human','P_S'),('mouse','P_S'),('human','P_M'),('mouse','P_M')]
    r.label('A',0,0);r.text(31,5,'Obs. / null median',ha='center',fontsize=7.2);r.text(60,5,'Obs./median',ha='center',fontsize=7.2);r.text(75,5,'Studies',ha='center',fontsize=7.2)
    ax=r.axes((23,12,27,29));shown=[]
    for i,(sp,co) in enumerate(contexts):
        r.text(2,12+i*7.5,f'{sp.title()}\n{co.replace("_","–")}',fontsize=7.2)
        for j,tau in enumerate([.02,.05]):
            q=rec[rec.species.eq(sp)&rec.contrast.eq(co)&rec.threshold.eq(tau)];assert len(q)==1;a=q.iloc[0];y=i+(j-.5)*.25
            ax.plot([a.null_p025/a.null_median,a.null_p975/a.null_median],[y,y],color='#C4CDD3',lw=1.7)
            ax.plot(a.observed/a.null_median,y,['o','D'][j],color=HOUSE_INK,ms=3.4)
            r.text(60,12+(y+.35)*7.25,f'{int(a.observed)}/{int(a.null_median)}',ha='center',fontsize=7.2)
            shown.append(a)
        r.text(75,13+i*7.25,str(int(a.n_studies)),ha='center',fontsize=7.2)
    ax.axvline(1,color=GREY,lw=.6,ls=':');ax.set_ylim(3.6,-.6);ax.set_yticks([]);ax.set_xticks([1,2,3,4]);ax.set_xlim(.6,max(4,rec.observed.div(rec.null_median).max()+.1))
    r.legend(2,44,[Line2D([],[],color=HOUSE_INK,marker=m,ls='',ms=3,label=f'≥{t}') for m,t in [('o','0.02'),('D','0.05')]]+[Line2D([],[],color='#C4CDD3',lw=1.7,label='95% null range')],ncol=3);r.display('A-tests',pd.DataFrame(shown))
    r.label('B',0,0);r.text(92,4,'P–S at 0.02',fontsize=7.2)
    controls=primary(r.read('measurement/nulls/recurrence_summary.tsv'),threshold=.02);controls=controls[controls.contrast.eq('P_S')]
    ladder=ladder_frame(r);ladder=ladder[ladder.contrast.eq('P_S')&ladder.threshold.eq(.02)]
    kinds=[('interval_supported','Bootstrap interval'),('depth20','20 PAS UMIs'),('cycle_matched','Cell-cycle matched'),('common','Common PAS axis'),('nested','Nested background')]
    ax=r.axes((120,9,34,17));control_rows=[]
    for i,(kind,label) in enumerate(kinds):
        for j,sp in enumerate(['human','mouse']):
            f=ladder[ladder.species.eq(sp)].sort_values('model').tail(1) if kind=='nested' else controls[controls.species.eq(sp)&controls.variant.eq(kind)]
            assert len(f)==1;a=f.iloc[0];y=i+(j-.5)*.22;scale=a.get('null_median',a.null_mean)
            ax.plot([a.null_p025/scale,a.null_p975/scale],[y,y],color='#C4CDD3',lw=1.5);ax.plot(a.observed/scale,y,'o',mfc=HOUSE_INK if j==0 else 'white',mec=HOUSE_INK,ms=3)
            control_rows.append(dict(a,display_control=kind,display_scale=scale))
    ax.set_yticks(range(5),[v[1] for v in kinds],fontsize=7.2);ax.set_ylim(4.6,-.6);ax.axvline(1,color=GREY,ls=':',lw=.5);ax.set_xticks([1,1.5,2]);ax.set_xlabel('Observed / null',fontsize=7.2,labelpad=1)
    r.text(91,32,'● Human   ○ Mouse',fontsize=7.2)
    mm=primary(r.read('measurement/nulls/recurrence_members.tsv.gz'),threshold=.02);mm=mm[mm.contrast.eq('P_S')]
    bars=[]
    for j,sp in enumerate(['human','mouse']):
        native=mm[mm.species.eq(sp)&mm.variant.eq('native')][['gene_id','direction']].drop_duplicates()
        for i,(variant,label) in enumerate([('interval_supported','Interval excludes 0'),('common','Common-axis identity')]):
            f=mm[mm.species.eq(sp)&mm.variant.eq(variant)][['gene_id','direction']].drop_duplicates();n=len(native.merge(f,on=['gene_id','direction']));den=len(native)
            x=117+j*21;y=38+i*4;ax=r.axes((x,y,9,1.5));ax.barh(0,1,color='#E1E7EB');ax.barh(0,n/den,color=HOUSE_INK);ax.set_xlim(0,1);ax.axis('off');r.text(x+10,y-.2,f'{n}/{den}',fontsize=7.2)
            bars.append(dict(species=sp,check=variant,n=n,denominator=den))
            if j==0:r.text(84,y-.2,label,fontsize=7.2)
    r.text(84,46,'Well P–M DWUI   145/199 eligible genes',fontsize=7.2)
    wells=r.read('well_inference/gene_metric_tests.tsv.gz');w=wells[wells.metric.eq('DWUI')];assert len(w)==199 and w.q_BH_joint.lt(.05).sum()==145
    r.display('B-control-tests',pd.DataFrame(control_rows));r.display('B-retention',pd.DataFrame(bars));r.display('B-well-tests',w)
    r.label('C',0,0);matched=r.read('matched_pms/matched_state_vectors.tsv.gz');matched=matched[matched.setting.eq('supported50_p08')]
    flow,detail=matched_flow(members,matched);r.display('C-flow',flow);r.display('C-membership',detail)
    native=r.read('state_vectors/route_effects.primary.tsv.gz');native=native[native.contrast.eq('P_M')&native.design.eq('same_library')]
    selections=[]
    for si,(sp,routes,count,top) in enumerate([('mouse',['CIP','WOT'],10,64),('human',['GSE118258'],8,106)]):
        info=detail[detail.species.eq(sp)&detail.complete].copy();coverage=matched[matched.species.eq(sp)&matched.route.isin(routes)].groupby('gene_id').agg(n_dates=('n_target_dates','sum'),n_libraries=('n_target_libraries','sum'),gene_symbol=('gene_symbol','first'))
        info=info.merge(coverage,on='gene_id').sort_values(['additional_PS','n_dates','n_libraries','gene_id'],ascending=[False,False,False,True])
        if sp=='mouse':info=pd.concat([info[info.gene_symbol.eq('Btg2')],info[~info.gene_symbol.eq('Btg2')]])
        chosen=info.head(count);selections.append(chosen.assign(selection='Btg2 first if present; then dual-reference support, target dates/libraries descending; gene ID tie-break'))
        r.text(8,top-8,sp.title());r.text(39,top-8,'Native P–M',ha='center',fontsize=7.2)
        for co,x in [('PS',61),('PM',90)]:r.text(x+10,top-8,'P–'+co[-1]+' ΔDWUI',ha='center',fontsize=7.2)
        for i,row in enumerate(chosen.itertuples()):
            yy=top+i*3.0
            if i%2==0:r.fig.add_artist(Rectangle((2/r.width,1-(yy+2.9)/r.height),106/r.width,3.2/r.height,transform=r.fig.transFigure,fc='#F0F3F5',ec='none',zorder=-1))
            r.text(2,yy,row.gene_symbol,fontstyle='italic',fontsize=7.2)
            voteroutes=['CIP','WOT'] if sp=='mouse' else ['GSE118258','GSE150311','GSE221739']
            for j,route in enumerate(voteroutes):
                f=native[native.gene_id.eq(row.gene_id)&native.route.eq(route)];v=f.delta_DWUI.iloc[0] if len(f)==1 else np.nan
                glyph='·' if not np.isfinite(v) else '+' if v>=.02 else '−' if v<=-.02 else '0'
                r.text(35+j*5,yy,glyph,color=STUDY_STYLE[route][0],fontsize=7.2)
        for co,x in [('PS',58),('PM',86)]:
            ax=r.axes((x,top-1,21,count*3.0));allvals=[]
            for j,route in enumerate(routes):
                f=matched[matched.route.eq(route)].set_index('gene_id').reindex(chosen.gene_id);vals=f['DWUI_'+co];allvals.extend(vals.dropna());color,marker=STUDY_STYLE[route]
                ax.scatter(vals,np.arange(count)+(j-.5)*.16,s=10,color=color,marker=marker,zorder=3)
            ax.set_ylim(count-.25,-.75);ax.set_yticks([]);ax.set_xlim(min(-.5,min(allvals)-.03),max(.25,max(allvals)+.03));ax.set_xticks([-.5,0,.25]);ax.xaxis.tick_top();ax.axvline(0,color=GREY,lw=.5);ax.spines['bottom'].set_visible(False);ax.spines['left'].set_visible(False)
    for i,route in enumerate(['CIP','WOT']):
        ax=r.axes((120,63+i*42,36,24));f=matched[matched.route.eq(route)].merge(detail[detail.species.eq('mouse')][['gene_id','additional_PS']],on='gene_id',how='left')
        cohort=f.additional_PS.notna();dual=f.additional_PS.eq(True)
        ax.scatter(f.DWUI_PS,f.DWUI_PM,s=3,color='#D1D9DE',lw=0,rasterized=True)
        ax.scatter(f.loc[cohort,'DWUI_PS'],f.loc[cohort,'DWUI_PM'],s=11,facecolors='white',edgecolors=HOUSE_INK,lw=.4)
        ax.scatter(f.loc[dual,'DWUI_PS'],f.loc[dual,'DWUI_PM'],s=11,color=HOUSE_INK,lw=0)
        ax.axhline(0,color=GREY,lw=.4,ls=':');ax.axvline(0,color=GREY,lw=.4,ls=':');ax.set_xlabel('P–S ΔDWUI',fontsize=7.2,labelpad=1);ax.set_ylabel('P–M ΔDWUI',fontsize=7.2,labelpad=1);ax.set_xticks([-1,0,1]);ax.set_yticks([-.5,0,.5]);ax.set_xlim(-1.03,1.03);ax.set_ylim(-1.03,1.03)
        r.text(138,58+i*42,LABELS[route],ha='center',fontsize=7.2)
    r.text(114,139,'Grey  All matched genes\n○ P–M cohort   ● Dual reference',fontsize=7.2)
    r.legend(2,132,[Line2D([],[],color=STUDY_STYLE[v][0],marker=STUDY_STYLE[v][1],ls='',ms=3,label=LABELS[v]) for v in ['CIP','WOT']],ncol=2)
    r.legend(2,136,[Line2D([],[],color=STUDY_STYLE[v][0],marker=STUDY_STYLE[v][1],ls='',ms=3,label=LABELS[v]) for v in ['GSE118258','GSE150311','GSE221739']],ncol=3)
    for i,sp in enumerate(['mouse','human']):
        f=flow[flow.species.eq(sp)].set_index('stage').n
        r.text(3,142+i*4,f'{sp.title()}   {f["Matched P/M/S"]}/{f["P–M recurrent"]} measurable;  {f["Also aligned P–S"]}/{f["Matched P/M/S"]} dual-reference',fontsize=7.2)
    r.display('C-example-rows',pd.concat(selections));r.display('C-all-matched',matched)
    r.label('D',0,0);r.text(12,154,'Btg2',fontstyle='italic',fontweight='bold')
    mc=pd.read_csv(r.resource('Btg2_metacells'),sep='\t');pairs=r.read('matched_pms/weighted_triplets.tsv.gz');axes=r.read('state_vectors/PAS_axes.primary.tsv.gz');checks=[]
    for i,route in enumerate(['CIP','WOT']):
        a=matched[matched.route.eq(route)&matched.gene_symbol.eq('Btg2')].iloc[0];pp=pairs[pairs.route.eq(route)&pairs.gene_symbol.eq('Btg2')&pairs.setting.eq('supported50_p08')]
        allowed=set(pp.library_key_a)|set(pp.library_key_b);f=mc[mc.route.eq(route)&mc.library_key.isin(allowed)&mc.qualified]
        x=13+i*36;ax=r.axes((x,165,24,18));means=[]
        for j,st in enumerate('SPM'):
            v=f.loc[f.node.eq(st),'DWUI'];lo,hi=v.quantile([.1,.9]);mean=a[st+'_DWUI'];means.append(mean)
            ax.plot([j,j],[lo,hi],color=STATE_COLORS[st],lw=.7);ax.plot(j,mean,marker=['^','o','D'][j],color=STATE_COLORS[st],ms=3.5)
            checks.append(dict(route=route,state=st,mean=mean,q10=lo,q90=hi,n_metacells=len(v)))
        ax.plot([0,1],means[:2],color=GREY,lw=.7);ax.axvline(1.5,color=GREY,ls=':',lw=.4);ax.set_xticks(range(3),list('SPM'));ax.set_xlim(-.3,2.3);ax.set_ylim(max(0,min(f.DWUI.min(),min(means))-.03),min(1,max(f.DWUI.max(),max(means))+.03))
        if i==0:ax.set_ylabel('DWUI',labelpad=1)
        r.text(x+12,161,LABELS[route],ha='center',fontsize=7.2)
        r.text(x+12,188,f'P–S {a.DWUI_PS:+.3f}\nP–M {a.DWUI_PM:+.3f}',ha='center',fontsize=7.2)
        ax=r.axes((x,204,24,15));ids=json.loads(a.PAS_ids);pa=axes[axes.study_id.eq(a.study_id)].drop_duplicates('atlas_pas_id').set_index('atlas_pas_id').loc[ids];signed=pa.cut_site.to_numpy()*np.where(pa.strand.eq('+'),1,-1);order=np.argsort(signed);pos=(signed[order]-signed.min())/1000
        for st,ls in [('S','--'),('P','-'),('M','-.')]:
            u=np.asarray(json.loads(a[st+'_usage']))[order];xx=np.r_[0,pos,pos[-1]+.02];yy=np.r_[100,100*(1-np.cumsum(u)),0];ax.step(xx,yy,where='post',color=STATE_COLORS[st],ls=ls,lw=.9)
            assert np.isclose(np.dot(json.loads(a[st+'_usage']),pa.position),a[st+'_DWUI'])
        ax.set_ylim(0,105);ax.set_yticks([0,50,100]);ax.set_xlim(0,pos[-1]+.02)
        if i==0:ax.set_ylabel('Distal fraction (%)',fontsize=7.2,labelpad=1)
    r.text(18,197,'– – S       — P       – · M',fontsize=7.2);r.text(40,223,'Distance from proximal PAS (kb)',ha='center',fontsize=7.2);r.display('D-state-values',pd.DataFrame(checks))
    r.label('E',0,0);r.text(91,154,'APA–GEX conditional tests',fontsize=8)
    nulls=r.read('apa_gex/four_class_nulls.tsv.gz');nulls=nulls[nulls.threshold.eq(.02)&nulls.baseline_bins.eq(5)&nulls.cell.eq('B')]
    tests=r.read('apa_gex/prespecified_eight_test_family.tsv');tests=tests[tests.threshold.eq(.02)&tests.baseline_bins.eq(5)]
    for j,model in enumerate(['response_identity','direction_specific']):
        for i,sp in enumerate(['human','mouse']):
            f=nulls[nulls.species.eq(sp)&nulls.null_model.eq(model)].sort_values('shuffled_study');test=tests[tests.species.eq(sp)&tests.null_model.eq(model)].iloc[0]
            ax=r.axes((99+i*34,166+j*32,20,12));
            for k,row in enumerate(f.itertuples()):
                ax.plot([row.null_p025,row.null_p975],[k,k],color='#B9C6CD',lw=2);ax.plot(row.null_mean,k,'|',ms=5,color=HOUSE_INK);ax.plot(row.observed,k,'o',ms=3.4,color=HOUSE_INK)
            ax.set_yticks([0,1],['Study 1','Study 2'],fontsize=7.2);ax.set_ylim(1.5,-.5);ax.set_xlim(0,max(f.observed.max(),f.null_p975.max())*1.12);ax.set_xticks([0,40,80]);ax.set_xlabel('')
            if j==0:r.text(107.5+i*34,161,sp.title(),ha='center')
            r.text(109+i*34,186+j*32,f'n={test.observed}; q={test.q_family:.4f}',ha='center',fontsize=7.2)
        r.text(84,166+j*32,'Identity' if j==0 else 'Direction',fontsize=7.2,rotation=90,ha='center')
    r.text(88,222,'Genes   ● Observed   | Mean   — 95% null',fontsize=7.2)
    r.display('E-null-values',nulls);r.display('E-family-tests',tests)
    r.save('Fig/Figure.3.pdf','Original A–E layout. Current recurrence and matched contributors replace all PAS values. B removes historical six-rule counts; bars are independent identity checks, not a sequential funnel. C scatter denominators are complete; table rows are explicitly selected examples. D Btg2 uses common-weight S/P/M means and observed metacell 10th–90th ranges. E displays same-APA/opposite-GEX cell B under both shuffled-study nulls, with the existing two-direction family test q; null intervals are not effect confidence intervals.')

def supplement11(r):
    r.start('Figure S11',150,196)
    r.label('A',2,2,'Metric-specific recurrence on fresh backgrounds')
    rec=r.read('recurrence/recurrence_summary.tsv');rec=rec[rec.series_matching]
    ax=r.axes((46,14,94,65));f=rec.sort_values(['metric','contrast','threshold','species'])
    null_forest(ax,f,[f'{"DWUI" if x.metric=="WPP" else x.metric} {x.species[0].upper()} {x.contrast.replace("_","–")} {x.threshold:g}' for x in f.itertuples()],q='q_all_24',ratio=True)
    r.label('B',2,90,'DWUI / MPRO membership')
    inventory=r.read('matched_pms/metric_membership.tsv.gz')
    planned=pd.MultiIndex.from_product([['human','mouse'],['P_S','P_M'],[.02,.05]],names=['species','contrast','DWUI_threshold'])
    counts=inventory.groupby(['species','contrast','DWUI_threshold','cohort']).size().unstack(fill_value=0).reindex(index=planned,columns=['DWUI_only','shared','MPRO_only'],fill_value=0)
    ax=r.axes((34,102,34,30));left=np.zeros(len(counts))
    for j,key in enumerate(counts):ax.barh(np.arange(len(counts)),counts[key],left=left,color=[BLUE,GREY,ORANGE][j],label=key.replace('_',' '));left+=counts[key]
    ax.set_yticks(range(len(counts)),[f'{a[0].upper()} {b.replace("_","–")} {c:g}/{.1 if c==.02 else .2:g}' for a,b,c in counts.index],fontsize=6.5)
    ax.invert_yaxis();ax.set_xlabel('Genes; thresholds DWUI / MPRO');ax.legend(frameon=False,loc='upper left',bbox_to_anchor=(-.95,-.31),ncol=3,fontsize=6.5,handlelength=.7,columnspacing=.5)
    r.label('C',78,90,'Matched P/M/S membership')
    flow,detail=matched_flow(r.read('recurrence/recurrence_members.tsv.gz'),r.read('matched_pms/matched_state_vectors.tsv.gz'))
    ax=r.axes((89,102,51,25))
    for sp,color in [('human',BLUE),('mouse',ORANGE)]:ax.plot(range(4),flow[flow.species.eq(sp)].n,'o-',color=color,ms=3,label=sp)
    ax.set_xticks(range(4),['P–M','Complete','Retain\nP–M','Also\nP–S']);ax.set_ylabel('Genes');ax.legend(frameon=False)
    # A full-width metacell panel gives the original example its actual new denominator.
    r.label('D',2,148,'Slc38a2: metacells and common-weight state means')
    ax=r.axes((22,159,117,28));example_panel(r,ax,'Slc38a2');ax.set_title('')
    r.save('Fig/SuppFig.S10.pdf','Fresh metric cohorts, 24-test correction and matched-reference admission. The Slc38a2 metacell count is recomputed and is not fixed at the historical value.')


def supplement10(r):
    r.start('Figure S10',150,194)
    tests=r.read('well_inference/gene_metric_tests.tsv.gz');rec=r.read('recurrence/recurrence_members.tsv.gz')
    tests['strict']=tests.q_BH_joint.lt(.05)&tests.effect.abs().ge(np.where(tests.metric.eq('MPRO'),.1,.02))
    r.label('A',2,2,'Joint five-metric independent-well tests')
    ax=r.axes((22,14,47,30));summary=tests.groupby('metric').agg(tested=('gene_id','size'),q=('joint_q05','sum'),strict=('strict','sum'))
    ax.barh(range(len(summary)),summary.tested,color='#E3E8EB');ax.barh(range(len(summary)),summary.q,color=GREY);ax.barh(range(len(summary)),summary.strict,color=BLUE)
    ax.set_yticks(range(len(summary)),summary.index);ax.set_xlabel('Genes')
    r.label('B',78,2,'DWUI / MPRO joint-q correspondence')
    ax=r.axes((93,14,47,30));wide=tests.pivot(index='gene_id',columns='metric',values='q_BH_joint')
    if len(wide):ax.scatter(-np.log10(wide.DWUI.clip(lower=1e-300)),-np.log10(wide.MPRO.clip(lower=1e-300)),s=4,color=INK,rasterized=True)
    else:empty(ax)
    for fn in [ax.axhline,ax.axvline]:fn(-np.log10(.05),color=GREY,ls='--',lw=.6)
    ax.set(xlabel='DWUI −log10(q)',ylabel='MPRO −log10(q)')
    r.label('C',2,61,'Current mouse P–M recurrence cohort')
    members=primary(rec,threshold=.02);members=members[members.species.eq('mouse')&members.contrast.eq('P_M')]
    rows=[]
    for metric in ['DWUI','RWUI','PDUI','PPUI','MPRO']:
        f=members[['gene_id','direction']].merge(tests[tests.metric.eq(metric)],on='gene_id',how='left',validate='one_to_one')
        concordant=np.sign(f.effect).eq(np.where(f.direction.eq('positive'),1,-1)*(-1 if metric=='PPUI' else 1)) if metric!='MPRO' else f.effect.notna()
        q=f.q_BH_joint.lt(.05);threshold=.1 if metric=='MPRO' else .02
        strict=q&concordant&f.effect.abs().ge(threshold)
        rows.append(dict(metric=metric,strict=int(strict.sum()),q_only=int((q&~strict).sum()),tested_other=int((f.effect.notna()&~q).sum()),unmeasured=int(f.effect.isna().sum())))
    cohort=pd.DataFrame(rows).set_index('metric');ax=r.axes((22,73,47,30));left=np.zeros(len(cohort))
    for key,color in [('strict',BLUE),('q_only',ORANGE),('tested_other',GREY),('unmeasured','#E7EBEE')]:
        ax.barh(range(len(cohort)),cohort[key],left=left,color=color,label=key.replace('_',' '));left+=cohort[key]
    ax.set_yticks(range(len(cohort)),cohort.index);ax.set_xlabel('All current recurrent genes');ax.legend(frameon=False,fontsize=6.5,loc='upper left',bbox_to_anchor=(-.35,-.3),ncol=2)
    r.label('D',78,61,'Documented-source endpoint inference')
    source=r.read('well_inference/source_endpoint_WPP_tests.tsv.gz');f=source[source.metric.eq('equal_MC_WPP')]
    ax=r.axes((93,73,47,30))
    for positive,marker in [(True,'o'),(False,'s')]:
        g=f[f.source_variance.gt(0).eq(positive)];ax.scatter(g.effect,-np.log10(g.q_BH_joint.clip(lower=1e-300)),s=5,marker=marker,facecolors=BLUE if positive else 'white',edgecolors=BLUE,linewidths=.4)
    zero(ax,False);ax.axhline(-np.log10(.05),color=GREY,ls='--');ax.set(xlabel='Equal-MC ΔDWUI',ylabel='−log10(joint q)')
    r.label('E',2,142,'Pool-level probability resolution')
    pool=r.read('pool_inference/gene_metric_tests.tsv.gz');ax=r.axes((17,154,51,28))
    p=pool.sort_values('minimum_two_sided_p')
    if len(p):
        ranks=np.arange(1,len(p)+1);ax.plot(ranks,p.minimum_two_sided_p,color=INK,label='Minimum attainable p');ax.plot(ranks,.05*ranks/(len(pool)),color=ORANGE,label='Joint BH boundary');ax.set_yscale('log')
    else:empty(ax)
    ax.set(xlabel='Ranked metric–gene tests',ylabel='Probability');ax.legend(frameon=False,fontsize=6.5)
    r.label('F',78,142,'Current mouse P–S candidates')
    roster=primary(rec,threshold=.02);roster=roster[roster.species.eq('mouse')&roster.contrast.eq('P_S')]
    f=roster[['gene_id','direction']].merge(source[source.metric.eq('equal_MC_WPP')],on='gene_id',how='left',validate='one_to_one')
    ax=r.axes((93,154,47,28));ax.scatter(f.effect,-np.log10(f.q_BH_joint.clip(lower=1e-300)),s=5,color=INK,rasterized=True)
    ax.set_title(f'{f.effect.notna().sum()} / {len(f)} tested',fontsize=7);zero(ax,False);ax.axhline(-np.log10(.05),color=GREY,ls='--');ax.set(xlabel='Equal-MC ΔDWUI',ylabel='−log10(joint q)')
    r.save('Fig/SuppFig.S9.pdf','Fresh gene-family correction and current recurrence rosters; wells/pools are inference units. Untested members remain in cohort denominators. Variance-zero sources are displayed separately.')


def supplement9(r):
    r.start('Figure S9',150,224)
    r.label('A',2,2,'Native measurement controls')
    f=primary(r.read('measurement/nulls/recurrence_summary.tsv'),threshold=.02)
    ax=r.axes((33,14,36,46));positions=grouped_null(ax,f,'variant','q_all_followup_tests')
    ax.set_yticklabels([{'native':'Native','interval_supported':'Bootstrap','depth20':'20 UMIs','depth40':'40 UMIs','cycle_matched':'Cycle match','common':'Common PAS'}[v] for v in f.variant.drop_duplicates()])
    thinning=primary(r.read('measurement/nulls/thinning_recurrence_all_draws.tsv'),threshold=.02)
    for i,row in enumerate(f.itertuples()):
        d=thinning[thinning.variant.eq(row.variant)&thinning.species.eq(row.species)&thinning.contrast.eq(row.contrast)]
        if len(d):ax.plot([d.recurrent_genes.min(),d.recurrent_genes.max()],[positions[i],positions[i]],color=INK,lw=1.3)
    r.label('B',78,2,'Historical six-rule consensus')
    historical=r.read('consensus/historical/figure3-joint-consensus-null/joint_consensus_summary.tsv')
    ax=r.axes((111,14,30,46));null_forest(ax,historical,[x.view.replace('common_contributors','Common').replace('fixed_genes','Fixed')+' '+x.species[0].upper()+f' {x.contrast[-1]} {x.threshold:g}' for x in historical.itertuples()],observed='observed_consensus_same_direction',q='q_within_view_four_tests')
    r.label('C',2,75,'Same-state half-split controls')
    half=primary(r.read('measurement/same-state/recurrence_summary.tsv'),threshold=.02)
    ax=r.axes((35,88,33,29))
    for i,row in enumerate(half.itertuples()):
        ax.plot([row.pseudo_median,row.pseudo_p95],[i,i],color=GREY,lw=2)
        ax.plot(row.observed_same_background_median,i,'o',color=BLUE if row.species=='human' else ORANGE,ms=3)
    ax.set_yticks(range(len(half)),[x.pseudo_kind+' '+x.species[0].upper()+' '+x.contrast[-1] for x in half.itertuples()]);ax.invert_yaxis();ax.set_xlabel('Recurrent genes')
    r.label('D',78,75,'Nested conditional backgrounds')
    f=ladder_frame(r);f=f[f.threshold.eq(.02)]
    f=f.copy();f['model']=f.model.replace({'L0_feasible_series':'L0 feasible','L1_measurement':'L1 measure','L2_baseline':'L2 baseline','L2_PAS_first':'L2 PAS','L3_full':'L3 full'})
    ax=r.axes((106,88,35,29));grouped_null(ax,f,'model','q_eight',ratio=True)
    r.label('E',2,132,'Exact human P–M study support')
    rec=primary(r.read('recurrence/recurrence_summary.tsv'));rec=rec[rec.species.eq('human')&rec.contrast.eq('P_M')]
    ax=r.axes((17,145,51,29))
    for j,row in enumerate(rec.to_dict('records')):
        for n in range(1,int(row.get('n_studies',0))+1):
            mean=row[f'exact_{n}_null_median'];lo=row[f'exact_{n}_null_p025'];hi=row[f'exact_{n}_null_p975'];obs=row[f'exact_{n}_observed'];x=n+(j-.5)*.15
            ax.plot([x,x],[max(1,lo),max(1,hi)],color=GREY,lw=2);ax.plot(x,max(1,mean),'_',color=INK)
            ax.plot(x,max(1,obs),['o','s'][j],color=[BLUE,ORANGE][j],ms=3,label=f'{row["threshold"]:g}' if n==1 else None)
    ax.set(xlabel='Exact number of supporting studies',ylabel='Genes (values ≤ 1 at 1)',yscale='log');ax.legend(frameon=False)
    ax.set_xticks(range(1,int(rec.n_studies.max())+1))
    r.label('F',78,132,'Human P–M study omissions')
    f=primary(r.read('recurrence/leave_one_study_out_nulls.tsv'));f=f[f.species.eq('human')&f.contrast.eq('P_M')&f.status.ne('omitted_study_not_evaluable')]
    ax=r.axes((93,145,47,29));order=list(f.omitted_study.drop_duplicates())
    for j,t in enumerate(sorted(f.threshold.unique())):
        g=f[f.threshold.eq(t)].set_index('omitted_study').reindex(order)
        for i,row in enumerate(g.itertuples()):
            x=i+(j-.5)*.15
            if row.status!='tested':ax.text(x,0,'×',ha='center',color=GREY);continue
            ax.plot([x,x],[row.null_p025,row.null_p975],color=GREY,lw=2)
            ax.plot(x,row.null_median,'_',color=INK)
            ax.plot(x,row.observed,'o',mec=[BLUE,ORANGE][j],mfc=[BLUE,ORANGE][j] if row.q_metric_4<.05 else 'white',
                ms=3,label=f'{t:g}' if i==0 else None)
    ax.set_xticks(range(len(order)),[v.replace('GSE','') for v in order],rotation=20);ax.set_ylabel('Recurrent genes');ax.set_xlabel('Omitted study');ax.legend(frameon=False)
    r.label('G',2,190,'Guan series as separate units: sensitivity only')
    native=primary(r.read('recurrence/recurrence_summary.tsv'),threshold=.02);native=native[native.species.eq('human')&native.contrast.eq('P_S')]
    series=primary(r.read('recurrence_controls/guan_series_unit_sensitivity.tsv'),threshold=.02);series=series[series.species.eq('human')&series.contrast.eq('P_S')]
    f=pd.concat([native,series]);ax=r.axes((47,202,93,12));null_forest(ax,f,['Guan one study','Guan separate series'])
    r.save('Fig/SuppFig.S8.pdf','All native controls recomputed. B alone retains historical original-pipeline consensus and is excluded from native inference. A shows all thinning draws; E preserves exact support counts; G is a pseudounit sensitivity, not additional independent replication.')
