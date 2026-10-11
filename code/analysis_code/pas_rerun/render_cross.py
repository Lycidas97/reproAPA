"""New cross-species cohorts, complete test families and state-context figures."""
import json
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

from render_common import BLUE,ORANGE,INK,GREY,DIVERGE,LABELS,empty,zero
from render_recurrence import null_forest

CODES=['H1','H2','H3','H4','C','W']
CODE_LABELS={'H1':'Xing','H2':'Liu','H3':'Guan–0618','H4':'Panariello','C':'Guo–CIP','W':'Schiebinger'}


def sharing_labels(frame):
    return [LABELS.get(x.human_study,x.human_study).replace('2020','').replace('GSE178325','Guan')+' / '+
            ('Guo–CIP' if x.mouse_study=='GSE103221' else 'Schiebinger')+f' · {x.threshold:g}' for x in frame.itertuples()]


def interval_rows(ax,frame,labels,value,lo,hi,q=None):
    ax.set_yticks(range(len(frame)),labels);ax.set_ylim(len(frame)-.5,-.5);zero(ax,False)
    for i,row in enumerate(frame.to_dict('records')):
        if not np.isfinite(row.get(value,np.nan)):
            ax.text(.01,i,'NA',transform=ax.get_yaxis_transform(),fontsize=6.5,color=GREY);continue
        color=BLUE if row.get('discovery','human')=='human' else ORANGE
        if np.isfinite(row.get(lo,np.nan)) and np.isfinite(row.get(hi,np.nan)):
            ax.plot([row[lo],row[hi]],[i,i],color=color,lw=1)
        ax.plot(row[value],i,'o',ms=3,mec=color,mfc=color if q is None or row[q]<.05 else 'white')


def endpoint_example(r,gene,box,contrast):
    x,y,w,h=box
    vectors=r.read('state_vectors/route_effects.primary.tsv.gz')
    all_axes=r.read('state_vectors/PAS_axes.primary.tsv.gz')
    for i,(sid,route) in enumerate([('GSE118258','GSE118258'),('GSE103221','CIP')]):
        ax=r.axes((x+i*w/2,y,w/2-9,h))
        q=vectors[vectors.gene_symbol.str.upper().eq(gene)&vectors.study_id.eq(sid)&vectors.route.eq(route)&
            vectors.contrast.eq(contrast)&vectors.design.eq('true_source' if contrast=='P_S' else 'same_library')]
        assert len(q)<=1
        if not len(q):empty(ax,'No qualified '+contrast.replace('_','–')+' vector');continue
        a=q.iloc[0];ids=json.loads(a.PAS_ids)
        axes=all_axes[all_axes.study_id.eq(sid)].drop_duplicates('atlas_pas_id').set_index('atlas_pas_id').loc[ids]
        pos=axes.position.to_numpy()
        for field,color,marker,label in [('reference_usage',GREY,'o',contrast[-1]),('target_usage',BLUE,'s','P')]:
            usage=np.asarray(json.loads(a[field]));ax.vlines(pos,0,usage,color=color,lw=.6)
            ax.scatter(pos,usage,s=9,facecolors='white' if label!='P' else color,edgecolors=color,marker=marker,label=label)
        other=vectors[vectors.study_id.eq(sid)&vectors.route.eq(route)&vectors.gene_id.eq(a.gene_id)&
            ((vectors.contrast.eq('P_S')&vectors.design.eq('true_source'))|(vectors.contrast.eq('P_M')&vectors.design.eq('same_library')))]
        effect='; '.join(z.contrast.replace('_','–')+f' {z.delta_DWUI:+.3f}' for z in other.itertuples())
        ax.set(title=CODE_LABELS['H1' if i==0 else 'C'],xlim=(-.03,1.03),ylim=(0,1.03),xlabel='PAS position',ylabel='PAS usage' if i==0 else '')
        ax.text(.5,1.2,effect,transform=ax.transAxes,ha='center',fontsize=6.5)
        if i==0:ax.legend(frameon=False,ncol=2,fontsize=6.5)


def figure4(r):
    from render_common import HOUSE_MAP,HOUSE_INK,STATE_COLORS,STUDY_STYLE
    r.main(4);r.label('A',0,0)
    sharing=r.read('cross_sharing/sharing_summary.tsv');w=sharing[sharing.metric.isin(['WPP','RWUI'])&sharing.test.eq('same_direction')]
    pairs=[('GSE118258','GSE103221'),('GSE118258','GSE115943'),('GSE178325','GSE103221'),('GSE178325','GSE115943')]
    pairnames=['Xing / Guo','Xing / Sch.','Guan / Guo','Guan / Sch.']
    for j,metric in enumerate(['WPP','RWUI']):
        ax=r.axes((29+j*46,11,22,32));r.text(40+j*46,2,'DWUI' if metric=='WPP' else 'RWUI',ha='center');r.text(60+j*40,7,'q24',ha='center',fontsize=7.2)
        for i,(hs,ms) in enumerate(pairs):
            for k,tau in enumerate([.02,.05]):
                q=w[w.metric.eq(metric)&w.human_study.eq(hs)&w.mouse_study.eq(ms)&w.threshold.eq(tau)];assert len(q)==1;a=q.iloc[0];y=i+(k-.5)*.3
                ax.plot([a.null_p025/a.null_mean,a.null_p975/a.null_mean],[y,y],color='#BAC6CE',lw=1.6)
                ax.plot(a.enrichment,y,['o','D'][k],ms=3.3,mec=HOUSE_INK,mfc=HOUSE_INK if a.q_joint_primary_24<.05 else 'white')
                r.text(60+j*40,12+(y+.4)*7.6,f'{a.q_joint_primary_24:.4f}',ha='center',fontsize=7.2)
        ax.set_ylim(3.6,-.6);ax.set_xlim(.7,1.42);ax.set_xticks([.8,1,1.3]);ax.axvline(1,color=GREY,ls=':',lw=.5);ax.set_yticks(range(4),pairnames if j==0 else ['']*4,fontsize=7.2);ax.set_xlabel('Obs. / null',fontsize=7.2,labelpad=1)
    r.legend(3,51,[Line2D([],[],color=HOUSE_INK,mfc='white',marker=m,ls='',ms=3,label=t) for m,t in [('o','0.02'),('D','0.05')]]+[Line2D([],[],color=HOUSE_INK,marker='o',ls='',ms=3,label='q24 < .05')],ncol=3);r.display('A-sharing',w)
    r.label('B',0,0);q=r.read('cross_continuous/continuous_gene_block_summary.tsv');q=q[q.layer.eq('supported50')&q.state.isin(list('REP'))]
    ax=r.axes((123,8,34,10))
    for i,a in enumerate(q.itertuples()):ax.plot([a.lo,a.hi],[i,i],color='#758791',lw=1);ax.plot(a.rho,i,'o',color=HOUSE_INK,ms=3)
    ax.set_ylim(2.6,-.6);ax.set_yticks(range(3),q.state+'–S');ax.axvline(0,color=GREY,ls=':',lw=.5);ax.set_xticks([-.04,0,.04]);ax.set_xlim(-.05,.065);ax.set_xlabel('Mean Spearman ρ',labelpad=1,fontsize=7.2)
    r.display('B-continuous',q)
    r.label('C',0,0);sets=r.read('cross_overlap/overlap_null_results.tsv');sets=sets[sets.model.str.startswith('independent_')]
    for j,co in enumerate(['PS','PM']):
        f=sets[sets.contrast.eq(co)].sort_values('model');assert len(f)==2;ax=r.axes((119+j*22,37,17,9))
        for k,a in enumerate(f.itertuples()):ax.plot([k,k],[a.null_q025,a.null_q975],color=GREY,lw=1);ax.plot(k,a.null_mean,['o','s'][k],mfc='white',mec='#758791',ms=3)
        obs=f.observed.iloc[0];ax.axhline(obs,color=HOUSE_INK,lw=.8);ax.text(.5,obs+1,str(obs),ha='center',fontsize=7.2)
        ax.set_xlim(-.5,1.5);ax.set_ylim(0,max(obs+5,f.null_q975.max()+2));ax.set_xticks([0,1],['Q5','R10']);ax.set_yticks([0,15,30] if j==0 else []);r.text(127.5+j*22,29,'P–'+co[-1],ha='center');r.text(127.5+j*22,33,'observed S' if co=='PS' else 'paired',ha='center',fontsize=7.2)
        for k,a in enumerate(f.itertuples()):r.text(123.2+j*22+k*8.5,51,f'{a.q_within_model_two_contrasts:.3f}',ha='center',fontsize=7.2)
        if j==0:ax.set_ylabel('Genes',labelpad=1,fontsize=7.2)
    r.display('C-set-tests',sets)
    r.label('D',0,0)
    roster=r.read('cross_trace/member_migration.tsv');roster=roster[roster.new_member&roster.contrast.eq('P_S')]
    trace=r.read('cross_trace/member_background_trace.tsv.gz');trace=trace[trace.new_member&trace.contrast.eq('P_S')&trace.background.eq('available')&trace.setting.eq('supported50_p08')]
    stable=r.read('cross_trace/member_setting_source_stability.tsv.gz');stable=stable[stable.new_member&stable.contrast.eq('P_S')]
    effects=r.read('cross_trace/member_endpoint_effects.tsv.gz');effects=effects[effects.new_member&effects.contrast.eq('P_S')&effects.setting.eq('supported50_p08')]
    native_axes=r.read('state_vectors/PAS_axes.primary.tsv.gz');names=native_axes[['primary_gene_id','primary_gene_symbol']].drop_duplicates().groupby('primary_gene_id').primary_gene_symbol.first()
    roster=roster.assign(symbol=roster.human_gene_id.map(names).fillna(roster.human_gene_id))
    keys=['human_gene_id','species','study_id','route'];opp=trace.loc[trace.member_proximal_background_distal,keys].drop_duplicates();matched=stable.merge(opp,on=keys,validate='many_to_one')
    opposite=set(opp.human_gene_id);checked=set(matched.loc[matched.new_endpoint_support&matched.mean_source_checked_sign.gt(0),'human_gene_id'])
    four=matched.groupby(keys).new_endpoint_support.agg(['size','all']);robust=set(four[(four['size']==4)&four['all']].index.get_level_values('human_gene_id'))
    roster=roster.assign(opposite=roster.human_gene_id.isin(opposite),source_checked=roster.human_gene_id.isin(checked),four_settings=roster.human_gene_id.isin(robust)).sort_values(['four_settings','source_checked','opposite','symbol'],ascending=[False,False,False,True])
    routes=['GSE118258','178325_0618','178325_0809','CIP','WOT'];assert len(roster)==26
    values=effects.pivot(index='human_gene_id',columns='route',values='delta').reindex(index=roster.human_gene_id,columns=routes)
    supports=effects.pivot(index='human_gene_id',columns='route',values='new_endpoint_support').reindex(index=roster.human_gene_id,columns=routes).fillna(False)
    globalmean=r.read('global/global_direction_summary.tsv');globalmean=globalmean[globalmean.setting.eq('supported50_p08')&globalmean.background.eq('available')&globalmean.contrast.eq('P_S')].set_index('route')['mean']
    for j,label in enumerate(['Xing','Guan\n0618','Guan\n0809','Guo\nCIP','Sch.\n2i']):r.text(30.6+j*8.4,59,label,ha='center',fontsize=7.2)
    r.text(2,72,'Mean sign',fontsize=7.2)
    for j,route in enumerate(routes):r.text(30.6+j*8.4,72,'+' if globalmean[route]>0 else '−',ha='center',color=HOUSE_MAP(.85 if globalmean[route]>0 else .15))
    r.text(76.4,77,'1  2  3',ha='center',fontsize=7.2)
    ax=r.axes((26.4,83,42,91));limit=.10
    im=ax.imshow(np.ma.masked_invalid(values.to_numpy()),cmap=HOUSE_MAP,vmin=-limit,vmax=limit,aspect='auto');ax.set_yticks(range(26),roster.symbol,fontsize=8,fontstyle='italic');ax.set_xticks([]);ax.tick_params(length=0,pad=2)
    for i in range(26):
        for j in range(5):
            val=values.iloc[i,j]
            if np.isfinite(val) and abs(val)>=.02:ax.text(j,i,'+' if val>0 else '−',ha='center',va='center',fontsize=7.2)
            if supports.iloc[i,j]:ax.add_patch(Rectangle((j-.47,i-.47),.94,.94,fill=False,lw=.6,ec=HOUSE_INK))
        for k,field in enumerate(['opposite','source_checked','four_settings']):ax.plot(5+k*.45,i,'o',ms=2.5,mfc=HOUSE_INK if roster.iloc[i][field] else 'white',mec=HOUSE_INK,clip_on=False,mew=.5)
    ax.set_xlim(-.5,4.5)
    for sp in ax.spines.values():sp.set_visible(False)
    cb=r.axes((27,178,40,1.4));r.fig.colorbar(im,cax=cb,orientation='horizontal',ticks=[-limit,0,limit],extend='both');r.text(43,185,'ΔDWUI; grey NA; ± at |Δ| ≥ 0.02',ha='center',fontsize=7.2)
    r.text(3,190,'□ Proximal support in defining routes',fontsize=7.2)
    r.text(3,195,f'1  Proximal in distal-mean route    {len(opposite)}/26',fontsize=7.2)
    r.text(3,199,f'2  + Source-checked background    {len(checked)}/26',fontsize=7.2)
    r.text(3,203,f'3  Four settings in such a route    {len(robust)}/26',fontsize=7.2)
    r.display('D-current-roster',roster.drop(columns=['old_member','new_member','migration']));r.display('D-effect-matrix',values.reset_index());r.display('D-background-means',globalmean.reset_index())
    r.label('E',0,0);r.text(97,58,'SZRD1',fontstyle='italic',fontweight='bold')
    matched=r.read('matched_pms/matched_state_vectors.tsv.gz');triplets=r.read('matched_pms/weighted_triplets.tsv.gz');mc=r.read('examples/example_metacell_scores.tsv.gz');case=[]
    for i,(sid,route) in enumerate([('GSE118258','GSE118258'),('GSE103221','CIP')]):
        f=matched[matched.setting.eq('supported50_p08')&matched.gene_symbol.str.upper().eq('SZRD1')&matched.study_id.eq(sid)&matched.route.eq(route)];assert len(f)==1;a=f.iloc[0]
        tt=triplets[triplets.setting.eq('supported50_p08')&triplets.study_id.eq(sid)&triplets.route.eq(route)&triplets.gene_id.eq(a.gene_id)];allowed=set(tt.library_key_a)|set(tt.library_key_b)
        x=99+i*35;ax=r.axes((x,72,25,16));vals=[a[s+'_DWUI'] for s in 'SPM']
        fmc=mc[mc.setting.eq('supported50_p08')&mc.study_id.eq(sid)&mc.route.eq(route)&mc.gene_symbol.str.upper().eq('SZRD1')&mc.qualified&mc.library_key.isin(allowed)]
        for j,st in enumerate('SPM'):
            v=fmc.loc[fmc.node.eq(st),'DWUI'];lo,hi=v.quantile([.1,.9]);ax.plot([j,j],[lo,hi],color=STATE_COLORS[st],lw=.6);ax.plot(j,vals[j],marker=['^','o','D'][j],color=STATE_COLORS[st],ms=3.5)
        ax.plot([0,1],vals[:2],color=GREY,lw=.7);ax.axvline(1.5,color=GREY,ls=':',lw=.4);ax.set_xticks(range(3),list('SPM'));ax.set_xlim(-.3,2.3);ax.set_ylim(0,1)
        if i==0:ax.set_ylabel('DWUI',labelpad=1)
        r.text(x+12.5,66,LABELS[route],ha='center',fontsize=7.2);r.text(x+12.5,93,f'P–S ΔDWUI {a.DWUI_PS:+.3f}',ha='center',fontsize=7.2)
        ids=json.loads(a.PAS_ids);pa=native_axes[native_axes.study_id.eq(sid)].drop_duplicates('atlas_pas_id').set_index('atlas_pas_id').loc[ids];signed=pa.cut_site.to_numpy()*np.where(pa.strand.eq('+'),1,-1);order=np.argsort(signed);pos=(signed[order]-signed.min())/1000
        ax=r.axes((x,105,25,18))
        for field,st,ls in [('S_usage','S','--'),('P_usage','P','-'),('M_usage','M','-.')]:
            u=np.asarray(json.loads(a[field]))[order];assert np.isclose(np.dot(json.loads(a[field]),pa.position),a[st+'_DWUI']);ax.step(np.r_[0,pos,pos[-1]+.02],np.r_[100,100*(1-np.cumsum(u)),0],where='post',color=STATE_COLORS[st],ls=ls,lw=.9)
        ax.set_ylim(0,105);ax.set_yticks([0,50,100]);ax.set_xlim(0,pos[-1]+.02)
        if i==0:ax.set_ylabel('Distal fraction (%)',fontsize=7.2,labelpad=1)
        case.append(a)
    r.text(107,98,'– – S     — P     – · M',fontsize=7.2);r.text(124,129,'Distance from proximal PAS (kb)',ha='center',fontsize=7.2);r.display('E-case-estimators',pd.DataFrame(case))
    r.label('F',0,0);gains=r.read('cis/paired_cv_uncertainty.tsv');gains=gains[gains.layer.eq('supported50')&gains.tau.eq(.02)&gains.model.eq('All cis')&gains.reference.eq('M1')]
    tests=r.read('cis/recomputed_models.tsv');tests=tests[tests.layer.eq('supported50')&tests.tau.eq(.02)]
    for i,discovery in enumerate(['human','mouse']):
        x=99+i*35;r.text(x+12,139,'Human → mouse' if i==0 else 'Mouse → human',ha='center',fontsize=7.2)
        for j,metric in enumerate(['auc','logloss']):
            ax=r.axes((x,149+j*30,24,16));f=gains[gains.discovery.eq(discovery)&gains.metric.eq(metric)].set_index('contrast').reindex(['P-S','P-M'])
            for k,a in enumerate(f.itertuples()):ax.plot([a.lo,a.hi],[k,k],color='#758791',lw=1);ax.plot(a.delta,k,'o' if k==0 else 'D',color=HOUSE_INK,ms=3.5)
            ax.set_ylim(1.5,-.5);ax.set_yticks([0,1],['P–S','P–M']);ax.axvline(0,color=GREY,lw=.5,ls=':')
            ax.set_xlim((-.13,.23) if metric=='auc' else (-.07,.06));ax.set_xticks([-.1,0,.2] if metric=='auc' else [-.05,0,.05]);ax.set_xlabel('Δ AUC' if j==0 else 'Δ log-loss',fontsize=7.2,labelpad=1)
    r.text(89,203,'All cis vs M1; conditional 95% intervals',fontsize=7.2)
    r.display('F-cis-gains',gains);r.display('F-joint-cis-tests',tests)
    r.save('Fig/Figure.4.pdf','Original A–F panel footprint. A retains all four study pairs and both thresholds, now DWUI/RWUI in the original two columns with primary q24 and Q5 nulls. B continuous effects and C separately selected proximal sets remain distinct. D displays all 26 current P–S members, including conflicting and missing contexts; no old/new subset interpretation. E uses one common-contributor P/M/S summary and all three PAS curves. F shows both predictive directions and both performance metrics for All cis vs M1; conditional intervals are not joint cis test p values.')

def supplement15(r):
    r.start('Figure S15',150,214)
    pairs=r.read('cross_continuous/continuous_pair_summary.tsv');pairs=pairs[pairs.layer.eq('supported50')]
    r.label('A',2,2,'Individual study-pair correlations · ortholog-block 95% intervals')
    pair_order=pd.MultiIndex.from_tuples([(h,m) for m in ['C','W'] for h in ['H1','H2','H3','H4']],names=['human','mouse'])
    short={'H1':'Xing','H2':'Liu','H3':'Guan0618','H4':'Panariello','C':'CIP','W':'WOT'}
    for j,state in enumerate('REPM'):
        q=pairs[pairs.state.eq(state)].set_index(['human','mouse']).reindex(pair_order).reset_index();ax=r.axes((31+j*29,15,20,45))
        labels=[short[a.human]+'/'+short[a.mouse] for a in q.itertuples()]
        interval_rows(ax,q,labels if j==0 else ['']*len(q),'rho','lo','hi');ax.set_title(state+'–S')
        # Preserve every pair's actual denominator, independent of symbol encoding.
        for i,a in enumerate(q.itertuples()):
            if np.isfinite(a.rho):ax.scatter(a.rho,i,s=5+a.n/100,facecolors='none',edgecolors=INK,lw=.4)
    r.fig.legend(handles=[ax.scatter([],[],s=5+n/100,facecolors='none',edgecolors=INK,lw=.4,label=f'n={n:,}') for n in [500,1500,2500]],
        loc='upper center',bbox_to_anchor=(.54,1-66/r.height),ncol=3,frameon=False)
    r.label('B',2,75,'Matched-state correlation differences · permutation ranges')
    tests=r.read('cross_continuous/state_and_matched_permutation.tsv')
    q=tests[tests.layer.eq('supported50')&tests.matched&tests.metric.eq('spearman')]
    ax=r.axes((36,88,104,18));interval_rows(ax,q,q.comparison.str.replace('_','–'),'observed','null_lo','null_hi','q_within_layer');ax.set_xlabel('P-state minus R/E/M-state Spearman ρ')
    r.label('C',2,121,'Recurrence-associated aligned effect · unique-gene 95% intervals')
    rec=r.read('cross_continuous/reciprocal_unique_gene_tests.tsv');q=rec[rec.primary]
    labels=[('Hs→Mm' if a.discovery=='human' else 'Mm→Hs')+' '+a.contrast for a in q.itertuples()]
    ax=r.axes((36,134,104,20));interval_rows(ax,q,labels,'observed_continuous_difference','stratified_gene_bootstrap_continuous_lo','stratified_gene_bootstrap_continuous_hi','exact_q_primary_four');ax.set_xlabel('Aligned ΔDWUI: recurrent minus matched controls')
    r.label('D',2,169,'Complete PAS-layer / threshold sensitivity')
    for j,layer in enumerate(['supported50','expanded_primary']):
        ax=r.axes((36+j*63,182,43,16));g=rec[rec.layer.eq(layer)]
        contexts=[('human','P-S'),('human','P-M'),('mouse','P-S'),('mouse','P-M')]
        for i,(sp,co) in enumerate(contexts):
            for k,tau in enumerate([.01,.02,.03,.05]):
                a=g[g.discovery.eq(sp)&g.contrast.eq(co)&g.tau.eq(tau)];assert len(a)==1
                v=a.iloc[0].get('observed_continuous_difference',np.nan)
                if np.isfinite(v):ax.plot(v,i+(k-1.5)*.18,['o','s','D','^'][k],ms=3,color=BLUE if v<0 else ORANGE,label=f'{tau:g}' if i==0 else None)
                else:ax.text(.02,i+(k-1.5)*.18,'×',transform=ax.get_yaxis_transform(),color=GREY,fontsize=6.5)
        zero(ax,False);ax.set_yticks(range(4),['Hs→Mm P–S','Hs→Mm P–M','Mm→Hs P–S','Mm→Hs P–M'] if j==0 else ['']*4,fontsize=6.5);ax.invert_yaxis();ax.set_title(layer.replace('_',' '),fontsize=7)
        if j==0:
            handles=[Line2D([],[],color=INK,marker=marker,ls='',ms=3,label=f'{tau:g}') for tau,marker in zip([.01,.02,.03,.05],['o','s','D','^'])]
            r.fig.legend(handles=handles,frameon=False,ncol=4,loc='upper center',bbox_to_anchor=(.53,1-207/r.height),fontsize=6.5)
    r.save('Fig/SuppFig.S14.pdf','A preserves each study-pair denominator; B compares jointly measurable state effects with its twelve-test family. C/D display recurrence-minus-control aligned effects, not unadjusted means; C intervals resample each unique ortholog once, retaining reused-control weights. All 32 positions, including unavailable settings, remain in the source tables.')


def supplement16(r):
    r.start('Figure S16',150,219)
    sharing=r.read('cross_sharing/sharing_summary.tsv');primary=sharing[sharing.test.isin(['same_direction','shared_unsigned'])]
    r.label('A',2,2,'All metrics · observed / null sharing and metric-specific q8')
    for i,metric in enumerate(['WPP','RWUI','MPRO']):
        q=primary[primary.metric.eq(metric)];ax=r.axes((47+i*34,15,26,33))
        null_forest(ax,q,sharing_labels(q) if i==0 else ['']*len(q),q='q_primary_8',ratio=True);ax.set_title('DWUI' if metric=='WPP' else metric);ax.set_xlabel('O / null mean')
    r.label('B',2,62,'Correction scope on the same current 24 primary tests')
    ax=r.axes((16,75,123,28))
    for a in primary.itertuples():
        ax.plot([0,1,2],[a.q_primary_8,a.q_joint_primary_24,a.q_joint_all_56],color={'WPP':BLUE,'RWUI':GREY,'MPRO':ORANGE}[a.metric],lw=.5,alpha=.65)
    ax.set_yscale('log');ax.axhline(.05,color=INK,ls=':',lw=.6);ax.set_xticks([0,1,2],['Within metric · 8','Across metrics · 24','Signed components · 56']);ax.set_ylabel('BH q')
    r.label('C',2,117,'Less-conditioned proximal-set null · full 16-position family')
    q=r.read('cross_overlap/overlap_threshold_layer_sensitivity.tsv');ax=r.axes((31,130,108,26))
    for j,(layer,co) in enumerate([(l,c) for l in ['supported50','expanded_primary'] for c in ['PS','PM']]):
        g=q[q.layer.eq(layer)&q.contrast.eq(co)]
        for k,a in enumerate(g.itertuples()):
            y=j+(k-1.5)*.17
            if a.null_mean<=0:ax.text(.02,y,f'O={a.observed:g}, E=0',transform=ax.get_yaxis_transform(),fontsize=6.5);continue
            ax.plot([a.null_q025/a.null_mean,a.null_q975/a.null_mean],[y,y],color=GREY,lw=1)
            ax.plot(a.observed/a.null_mean,y,['o','s','D','^'][k],ms=3,mec=BLUE if co=='PS' else ORANGE,mfc=(BLUE if co=='PS' else ORANGE) if a.q_16<.05 else 'white')
    ax.set_yticks(range(4),['Ref P–S','Ref P–M','Expanded P–S','Expanded P–M']);ax.invert_yaxis();ax.set_xlabel('Observed / null mean')
    r.label('D',2,171,'Nested discrete null distributions')
    ladder=r.read('cross_controls/cross_species_ladder.tsv')
    with np.load(r.artifact('cross_controls/cross_species_draws.npz')) as draws:
        for j,co in enumerate(['PS','PM']):
            ax=r.axes((14,184+j*12,57,9))
            for model,color in [('L1_capacity',GREY),('L2_baseline',BLUE),('L3_full',ORANGE)]:
                vals=draws[co+'|'+model];count=np.bincount(vals.astype(int));ax.step(range(len(count)),count/count.sum(),where='mid',color=color,label=model.split('_',1)[-1])
            q=ladder[ladder.contrast.eq(co)];assert q.observed.nunique()==1;ax.axvline(q.observed.iloc[0],color=INK,lw=1);ax.set_ylabel('P–'+co[-1]);ax.set_xlabel('Shared genes' if j else '')
            if j==0:ax.tick_params(axis='x',labelbottom=False)
            if j==0:ax.legend(frameon=False,ncol=3,fontsize=6.5,loc='upper left',bbox_to_anchor=(-.1,1.9))
    r.label('E',79,171,'Conditional detection power')
    power=r.read('cross_controls/cross_species_power.tsv');power=power[power.experiment.eq('conditional_injection')&power.alpha.eq(.05)]
    ax=r.axes((94,184,46,23))
    for co,color in [('PS',BLUE),('PM',ORANGE)]:
        for model,ls in [('L1_capacity','--'),('L3_full','-')]:
            q=power[power.contrast.eq(co)&power.tested_model.eq(model)].sort_values('odds_ratio');ax.plot(q.odds_ratio,q.rejection_rate,ls=ls,marker='o',ms=2,color=color,label='P–'+co[-1]+' '+model.split('_',1)[-1])
    ax.set(ylim=(0,1),xlabel='Injected odds ratio',ylabel='Detection probability');ax.legend(frameon=False,fontsize=6.5)
    r.save('Fig/SuppFig.S15.pdf','Every panel is recomputed with current native PASs. The former less-conditioned analysis is rerun in C, not copied from its historical PDF. Narrow and broad sharing corrections retain distinct families. D uses discrete permutation distributions; E uses conditional injection and does not reinterpret null ranges as confidence intervals.')


def supplement17(r):
    r.start('Figure S17',150,160)
    roster=r.read('cross_trace/member_migration.tsv');selected=roster[roster.new_member]
    genes=sorted(selected.human_gene_id.unique());height=max(95,52+len(genes)*3.3)
    if height>250:raise ValueError('New selected ortholog set requires multi-page S17 reflow; no members dropped')
    r.height=height;r.fig.set_size_inches(150/25.4,height/25.4)
    r.label('A',2,2,f'State context of all {len(genes)} selected orthologs')
    if not genes:
        empty(r.axes((15,16,125,height-28)),'No ortholog passes current shared-set criteria');r.save('Fig/SuppFig.S16.pdf','The current shared ortholog set is empty; the historical roster is retained only in migration tables.');return
    nodes=r.read('cross_overlap/all_mapped_nodes.tsv.gz');nodes=nodes[nodes.layer.eq('supported50')&nodes.measurement_pass&nodes.human_gene_id.isin(genes)]
    wide=nodes.pivot(index=['human_gene_id','code'],columns='node',values='position_index').reindex(columns=list('SREPM'))
    effects=wide.sub(wide.S,axis=0).drop(columns='S');frame=effects.stack().unstack(['code','node'])
    columns=pd.MultiIndex.from_tuples([(c,s) for c in CODES for s in 'REPM']);frame=frame.reindex(index=genes,columns=columns)
    ax=r.axes((31,14,92,height-52));im=ax.imshow(np.ma.masked_invalid(frame.to_numpy()),cmap=DIVERGE,vmin=-.1,vmax=.1,aspect='auto')
    shared=r.read('cross_overlap/new_shared_members.tsv.gz')
    symbols=shared.drop_duplicates('human_gene_id').set_index('human_gene_id').human_gene_symbol.to_dict()
    ax.set_xticks(range(len(columns)),[CODE_LABELS[c]+' '+s for c,s in columns],rotation=90,fontsize=6.5);ax.set_yticks([])
    for i,j in zip(*np.where(abs(frame.to_numpy())>=.02)):ax.plot(j,i,'.',color=INK,ms=1.5)
    for boundary in np.arange(3.5,len(columns),4):ax.axvline(boundary,color='white',lw=1.2)
    membership=selected.pivot(index='human_gene_id',columns='contrast',values='new_member').reindex(index=genes,columns=['P_S','P_M']).fillna(False).astype(int)
    for i,gene in enumerate(genes):
        for code in CODES:
            j=columns.get_loc((code,'P'))
            if membership.loc[gene,'P_S'] and frame.iloc[i,j]<=-.02:
                ax.add_patch(Rectangle((j-.47,i-.47),.94,.94,fill=False,edgecolor=INK,lw=.6))
    strip=r.axes((25,14,5,height-52));strip.imshow(membership.to_numpy(),cmap='Greys',vmin=0,vmax=1,aspect='auto');strip.set_xticks([0,1],['P–S','P–M'],rotation=90,fontsize=6.5)
    strip.set_yticks(range(len(genes)),[symbols.get(g,g) for g in genes],fontsize=6.5)
    pm=r.read('cross_overlap/paired_pm_study_effects.tsv.gz');pm=pm[pm.layer.eq('supported50')&pm.quality_pass&pm.human_gene_id.isin(genes)]
    pm=pm[(pm.code.isin(['W','C'])&pm.n_sample_units.ge(2))|(pm.code.isin(['H1','H2','H3','H4'])&pm.n_sample_units.ge(1))]
    paired=pm.pivot(index='human_gene_id',columns='code',values='mean_delta').reindex(index=genes,columns=CODES)
    axp=r.axes((126,14,21,height-52));axp.imshow(np.ma.masked_invalid(paired.to_numpy()),cmap=DIVERGE,vmin=-.1,vmax=.1,aspect='auto');axp.set_xticks(range(6),[CODE_LABELS[c] for c in CODES],rotation=90,fontsize=6.5);axp.set_yticks([]);axp.set_title('Paired P–M',fontsize=7)
    for i,gene in enumerate(genes):
        for j,code in enumerate(CODES):
            a=paired.loc[gene,code];row=wide.reindex(pd.MultiIndex.from_tuples([(gene,code)]));b=(row.P-row.M).iloc[0]
            if membership.loc[gene,'P_M'] and a<=-.02:axp.add_patch(Rectangle((j-.47,i-.47),.94,.94,fill=False,edgecolor=INK,lw=.6))
            if np.isfinite(a) and np.isfinite(b) and np.sign(a)!=np.sign(b):axp.plot([j-.3,j+.3],[i+.3,i-.3],color=INK,lw=.6)
    cb=r.fig.add_axes([.33,.067,.32,.014]);r.fig.colorbar(im,cax=cb,orientation='horizontal',ticks=[-.1,0,.1]);cb.set_xlabel('Observed-S ΔDWUI / paired P–M ΔDWUI',fontsize=7)
    r.save('Fig/SuppFig.S16.pdf','The complete newly selected set is shown; old-only genes remain in migration tables. Grey denotes unavailable measurements and dots mark absolute effect ≥0.02. Outlines mark negative effects ≤−0.02 for selected P–S genes in P columns and selected P–M genes in the paired block. Slashes identify paired-versus-separate P–M direction disagreement. Selection and plotted observed-state estimators remain distinct.')
