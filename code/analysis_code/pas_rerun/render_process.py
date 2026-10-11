"""Common-reference process, APA/GEX and nested functional-test supplements."""
import json
import numpy as np
import pandas as pd

from render_common import BLUE,ORANGE,INK,GREY,LABELS,empty,zero

CLASS_NAMES={'endpoint_aligned':'Both endpoint-aligned','reverse':'Both reversed','within_threshold':'Both below threshold',
    'direction_difference':'Opposite directions','mixed_or_threshold_undetermined':'Mixed / threshold-undetermined'}
CLASS_COLORS=[BLUE,ORANGE,'#CBD4DA','#8064A2','#758B87']


def supplement13(r):
    r.start('Figure S13',150,218)
    flow=r.read('intermediate/process/coverage_flow.tsv')
    main=flow[flow.variant.eq('shared_p08')&flow.scope.eq('requalified')]
    r.label('A',2,2,'Endpoint admission and parallel state coverage')
    ax=r.axes((16,14,124,26));stages=['original_candidates','endpoint_measurable','endpoint_retained','R','E','complete_required_states']
    for sp,color,marker in [('human',BLUE,'o'),('mouse',ORANGE,'s')]:
        f=main[main.species.eq(sp)];assert len(f)<=1
        if not len(f):continue
        values=f[stages].iloc[0].to_numpy(float);ax.plot(range(3),values[:3],marker=marker,ms=3,color=color,label=sp)
        ax.scatter(np.arange(3,6),values[3:],marker=marker,s=12,color=color)
        for i,v in enumerate(values):
            if np.isfinite(v):ax.annotate(str(int(v)),(i,v),xytext=(2,-7 if sp=='human' else 6),textcoords='offset points',fontsize=6.5,color=color)
    ax.margins(y=.22);ax.text(4,.10,'n/a',transform=ax.get_xaxis_transform(),ha='center',color=BLUE,fontsize=6.5)
    ax.set_xticks(range(6),['Current P–S','Endpoint\nmeasured','Endpoint\nretained','R measured','E measured','Required\nstates complete']);ax.set_ylabel('Genes');ax.legend(frameon=False)
    pairs=r.read('intermediate/process/same_series_process_classes.tsv.gz')
    base=pairs[pairs.variant.eq('shared_p08')&pairs.scope.eq('requalified')&pairs.endpoint_retained&pairs.pair_measurable]
    r.label('B',2,56,'Intermediate shifts aligned to each study’s endpoint')
    for i,(sp,state) in enumerate([('human','R'),('mouse','R'),('mouse','E')]):
        ax=r.axes((14+i*47,68,31,27));f=base[base.species.eq(sp)&base.state.eq(state)]
        if len(f):
            for (name,color) in zip(CLASS_NAMES,CLASS_COLORS):
                g=f[f.pair_class.eq(name)];ax.scatter(g.aligned_delta_1,g.aligned_delta_2,s=3,color=color,alpha=.65,rasterized=True)
        else:empty(ax)
        for v in [-.02,.02]:ax.axhline(v,color=GREY,ls=':',lw=.5);ax.axvline(v,color=GREY,ls=':',lw=.5)
        ax.set(title=f'{sp.title()} {state} · n={len(f)}',xlabel='Study 1 aligned ΔDWUI',ylabel='Study 2 aligned ΔDWUI')
    r.label('C',2,111,'Five exhaustive intermediate classes')
    ax=r.axes((27,123,42,23));order=[('human','R'),('mouse','R'),('mouse','E')];left=np.zeros(3)
    for key,color in zip(CLASS_NAMES,CLASS_COLORS):
        values=[int((base.species.eq(sp)&base.state.eq(state)&base.pair_class.eq(key)).sum()) for sp,state in order]
        ax.barh(range(3),values,left=left,color=color,label=CLASS_NAMES[key]);left+=values
    ax.set_yticks(range(3),['Human R','Mouse R','Mouse E']);ax.set_xlabel('Qualified genes');ax.legend(frameon=False,loc='upper left',bbox_to_anchor=(-.5,-.28),fontsize=6.5)
    r.label('D',78,111,'Opposite-direction fraction across settings')
    q=pairs[pairs.scope.eq('requalified')&pairs.endpoint_retained&pairs.pair_measurable]
    totals=q.groupby(['species','state','variant']).size().rename('n');opposite=q[q.pair_class.eq('direction_difference')].groupby(['species','state','variant']).size()
    summary=pd.DataFrame({'n':totals,'opposite':opposite.reindex(totals.index,fill_value=0)}).reset_index()
    variants=list(pairs.variant.drop_duplicates());ax=r.axes((93,123,47,23))
    for i,(sp,state) in enumerate(order):
        g=summary[summary.species.eq(sp)&summary.state.eq(state)].set_index('variant').reindex(variants)
        ax.plot(range(len(g)),g.opposite/g.n,marker=['o','s','^'][i],ms=2.5,color=[BLUE,ORANGE,'#8064A2'][i],label=f'{sp[0].upper()} {state}')
    ax.set_xticks(range(len(variants)),[v.replace('shared_','').replace('exclude_','−').replace('native_PAS_coordinates','native').replace('support_','GEX ') for v in variants],rotation=60,ha='right',fontsize=6.5);ax.set_ylabel('Opposite / qualified');ax.legend(frameon=False,fontsize=6.5)
    summary.to_csv(r.out/'S13-process-sensitivity-denominators.tsv',sep='\t',index=False)
    r.label('E',2,175,'Frozen expression-support retention')
    support=r.read('intermediate/process/frozen_expression_support_summary.tsv');support=support[support.node.isin(list('SREP'))]
    ax=r.axes((16,187,124,11));x=np.arange(len(support));ax.bar(x,support.support_q95/support.metacells,color=GREY)
    ax.scatter(x,support.support_q90/support.metacells,color=INK,marker='_',s=12)
    ax.set_xticks(x,[s.study_id.replace('GSE','')+' '+s.node for s in support.itertuples()],rotation=45,ha='right',fontsize=6.5);ax.set_ylim(0,1);ax.set_ylabel('Retained / MCs')
    r.save('Fig/SuppFig.S12.pdf','Fresh endpoint admission, fixed-source R/E effects and exhaustive process classes. Human E is unavailable in the specified same-route pair. D retains each setting’s own denominator; fixed-cohort attrition is separately archived. E reuses unchanged GEX support only.')


def supplement12(r):
    r.start('Figure S12',150,192)
    tests=r.read('apa_gex/prespecified_eight_test_family.tsv');four=r.read('apa_gex/four_class_nulls.tsv.gz')
    r.label('A',2,2,'Same-APA / opposite-GEX conditional nulls')
    t=tests[tests.baseline_bins.eq(5)].copy();ax=r.axes((53,14,86,39))
    for i,row in enumerate(t.itertuples()):
        g=four[four.species.eq(row.species)&four.threshold.eq(row.threshold)&four.null_model.eq(row.null_model)&four.baseline_bins.eq(5)&four.cell.eq('B')]
        for j,b in enumerate(g.itertuples()):ax.plot([b.null_p025,b.null_p975],[i+(j-.5)*.2]*2,color=GREY,lw=2)
        color=BLUE if row.species=='human' else ORANGE
        ax.plot(row.observed,i,'o',mec=color,mfc=color if row.q_family<.05 else 'white',ms=3)
    ax.set_yticks(range(len(t)),[x.species.title()+' '+x.null_model.replace('_',' ')+f' {x.threshold:g}' for x in t.itertuples()],fontsize=6.5);ax.invert_yaxis();ax.set_xlabel('Genes; ranges summarize two null orientations')
    r.label('B',2,68,'All four APA/GEX sign combinations')
    for i,sp in enumerate(['human','mouse']):
        for j,threshold in enumerate([.02,.05]):
            ax=r.axes((15+j*64,80+i*22,50,15));f=four[four.species.eq(sp)&four.threshold.eq(threshold)&four.null_model.eq('response_identity')&four.baseline_bins.eq(5)&four.shuffled_study.eq(1)].set_index('cell').reindex(list('ABCD'))
            if len(f):ax.bar(range(4),f.observed,color=[GREY,BLUE,ORANGE,'#8064A2'])
            ax.set_xticks(range(4),['Same/Same','Same/Opp.','Opp./Same','Opp./Opp.'],rotation=20,ha='right',fontsize=6.5);ax.set_title(f'{sp.title()} · {threshold:g}',fontsize=7);ax.set_ylabel('Genes')
    r.label('C',2,139,'New endpoint-weighted PAS fractions')
    vectors=r.read('state_vectors/route_effects.primary.tsv.gz')
    intended=(vectors.species.eq('mouse')&vectors.gene_symbol.str.upper().eq('RTF2'))|(vectors.species.eq('human')&vectors.gene_symbol.str.upper().eq('AP1S2'))
    q=vectors[vectors.contrast.eq('P_S')&vectors.design.eq('true_source')&intended]
    axes=r.read('state_vectors/PAS_axes.primary.tsv.gz')
    ax=r.axes((24,151,47,29));labels=[];y=0
    for row in q.itertuples():
        position=axes[axes.study_id.eq(row.study_id)].drop_duplicates('atlas_pas_id').set_index('atlas_pas_id').loc[json.loads(row.PAS_ids)].position
        order=np.argsort(position.to_numpy(),kind='stable')
        for state,column in [('S','reference_usage'),('P','target_usage')]:
            usage=np.array(json.loads(getattr(row,column)))[order];left=0.
            for k,value in enumerate(usage):ax.barh(y,value,left=left,color=str(.2+.65*k/max(1,len(usage)-1)),height=.8);left+=value
            labels.append(row.gene_symbol+' '+LABELS.get(row.route,row.route)+' '+state);y+=1
    ax.set_yticks(range(y),labels,fontsize=6.5);ax.invert_yaxis();ax.set_xlim(0,1);ax.set_xlabel('PAS usage, proximal → distal')
    if not y:empty(ax,'Neither example has a qualified P–S vector')
    r.label('D',78,139,'GEX on the same\nPAS-qualified contributors')
    ax=r.axes((110,151,31,29))
    for i,row in enumerate(q.itertuples()):ax.plot([0,row.GEX_delta],[i,i],color=BLUE if row.GEX_delta<0 else ORANGE);ax.plot(row.GEX_delta,i,'o',color=INK,ms=3)
    ax.set_yticks(range(len(q)),[row.gene_symbol+' '+LABELS.get(row.route,row.route) for row in q.itertuples()],fontsize=6.5);ax.invert_yaxis();zero(ax,False);ax.set_xlabel('Δ mean log2(CPM+1)')
    r.save('Fig/SuppFig.S11.pdf','Both conditional-null definitions and both orientations are rebuilt. C uses full current PAS vectors and D the same qualified contributors. Missing example routes remain absent with explicit eligibility in source data; no historical effect labels are carried forward.')


def supplement14(r):
    r.start('Figure S14',150,218)
    membership=r.read('functional/nested_universe_membership.tsv.gz');tables=[]
    for db in ['go','kegg']:tables.append(r.read(f'functional/{db}_ladder.tsv.gz'))
    all_terms=pd.concat(tables,ignore_index=True);exact=all_terms[all_terms.model.eq('L3_full')]
    r.label('A',2,2,'Nested expression / APA / recurrence universes')
    ax=r.axes((25,14,115,24));counts=membership.groupby(['species','contrast'])[['in_G0','in_G1_nested','in_G2']].sum()
    for j,key in enumerate(['in_G0','in_G1_nested','in_G2']):ax.bar(np.arange(len(counts))+(j-1)*.24,counts[key],.22,color=[GREY,BLUE,ORANGE][j],label=['G0 expressed','G1 APA-testable','G2 recurrent'][j])
    ax.set_xticks(range(len(counts)),[a.title()+' '+b.replace('_','–') for a,b in counts.index]);ax.set_ylabel('Genes');ax.legend(frameon=False,ncol=3)
    r.label('B',2,53,'Complete comparison-level exact-test support')
    selected=exact[exact.q_comparison.lt(.05)].sort_values(['database','q_comparison','analysis_id','term_id'])
    assert len(selected)<=40,'Functional support exceeds the complete two-column layout; reflow required'
    for j,indices in enumerate(np.array_split(np.arange(len(selected)),2)):
        ax=r.axes((43+j*72,66,26,70));panel=selected.iloc[indices]
        if panel.empty:empty(ax,'No supported term');continue
        for i,row in enumerate(panel.itertuples()):
            ax.plot([row.expected,row.hits],[i,i],color=GREY,lw=.8);ax.plot(row.expected,i,'o',mfc='white',mec=INK,ms=3);ax.plot(row.hits,i,'o',color=BLUE if row.q_database_model<.05 else ORANGE,ms=3)
        ax.set_yticks(range(len(panel)),[x.term_id+' '+x.species[0].upper()+' '+x.contrast.replace('_','–') for x in panel.itertuples()],fontsize=6.5)
        ax.set_ylim(len(panel)-.5,-.5);ax.set_xlabel('Expected / observed hits')
    r.label('C',2,160,'Support across conditioning levels')
    count=all_terms.groupby(['database','model']).q_comparison.apply(lambda x:int(x.lt(.05).sum())).reset_index(name='supported')
    ax=r.axes((17,173,51,32));models=list(count.model.drop_duplicates())
    for db,color in [('GO',BLUE),('KEGG',ORANGE)]:
        f=count[count.database.eq(db)].set_index('model').reindex(models);ax.plot(range(len(models)),f.supported,'o-',color=color,ms=3,label=db)
    ax.set_xticks(range(len(models)),[m.replace('_','\n') for m in models],fontsize=6.5);ax.set_ylabel('Supported term positions');ax.legend(frameon=False)
    r.label('D',78,160,'Exact conditional detection power')
    power=r.read('functional/functional_power.tsv.gz');power=power[power.model.eq('L3_full')&power.analysis_id.str.contains('G2_vs_G1')]
    ax=r.axes((92,173,49,32))
    for label,color in [('nominal',BLUE),('bonferroni_terms',ORANGE)]:
        f=power[power.threshold.eq(label)].groupby('odds_ratio').rejection_probability
        if not len(f):continue
        mid=f.median();lo=f.quantile(.25);hi=f.quantile(.75);ax.plot(mid.index,mid,'o-',color=color,ms=3,label=label.replace('_',' '));ax.fill_between(mid.index,lo,hi,color=color,alpha=.15)
    ax.set(xlabel='Injected term odds ratio',ylabel='Detection probability',ylim=(0,1));ax.legend(frameon=False,fontsize=6.5)
    r.save('Fig/SuppFig.S13.pdf','All current GO/KEGG families and nested universes are retained. Supported terms are selected only after full family correction; no former term list is fixed. Power bands are IQRs across analysis positions, not confidence intervals.')
