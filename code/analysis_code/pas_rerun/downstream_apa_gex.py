"""Fresh shared-axis APA/GEX labels and response-identity/direction nulls."""
import argparse
import itertools
import json

import numpy as np
import pandas as pd

from downstream_runtime import Context,bh,seed
from downstream_cross import qcodes
from downstream_shared_states import BASE
from polyapipe_effects import weights
from downstream_parallel import ordered_map

AXES={'mouse':'GSE103221__GSE115943','human':'GSE118258__GSE178325'}
STUDIES={'mouse':['GSE103221','GSE115943'],'human':['GSE118258','GSE178325']}


def sign(low,high,t):return np.where(low>=t,1,np.where(high<=-t,-1,0))


def four_counts(x,y,g):
    return np.array([((x*y>0)&(g==1)).sum(),((x*y>0)&(g==-1)).sum(),((x*y<0)&(g==1)).sum(),((x*y<0)&(g==-1)).sum()],int)


def inputs(ctx):
    root=ctx.analysis/'intermediate/shared';states=ctx.read(root/'state_values.tsv.gz');pairs=ctx.read(root/'matched_library_effects.tsv.gz')
    pairs=pairs[pairs.setting.eq('shared_p08')&pairs.contrast.eq('P_S')&pairs.design.eq('true_source')&pairs.axis_id.isin(AXES.values())].copy()
    st=states[states.setting.eq('shared_p08')]
    for suffix,node in [('a','P'),('b','S')]:
        values=['APA','GEX_CPM','site_counts'];cols=BASE+['library_key']
        z=st[st.state.eq(node)][cols+values].rename(columns={c:c+'_'+suffix for c in ['library_key',*values]})
        pairs=pairs.merge(z,on=BASE+['library_key_'+suffix],validate='many_to_one')
    route=[]
    for key,q in pairs.groupby(BASE+['contrast','design']):
        w=weights(q);reference=float(w@q.APA_b)
        route.append(dict(zip(BASE+['contrast','design'],key),APA=float(w@q.APA_delta),GEX=float(w@q.GEX_delta),
            GEX_all=float(w@q.GEX_all_delta),meas=bool(q.GEX_measurable.all()),reference_WPP=reference,
            depth=float(w@((q.pas_counts_a/q.n_mc_a+q.pas_counts_b/q.n_mc_b)/2)),
            PAS=len(json.loads(q.site_counts_a.iloc[0]))))
    r=pd.DataFrame(route,columns=BASE+['contrast','design','APA','GEX','GEX_all','meas','reference_WPP','depth','PAS'])
    rows=[]
    for key,q in r.groupby(['axis_id','species','study_id','setting','gene_id','gene_symbol']):
        rows.append(dict(zip(['axis_id','species','study_id','setting','gene_id','gene_symbol'],key),APA=q.APA.mean(),alo=q.APA.min(),ahi=q.APA.max(),
            GEX=q.GEX.mean(),glo=q.GEX.min(),ghi=q.GEX.max(),meas=q.meas.all(),depth=q.depth.mean(),PAS=q.PAS.max(),
            reference_WPP=q.reference_WPP.median(),reference_WPP_min=q.reference_WPP.min(),reference_WPP_max=q.reference_WPP.max(),routes=';'.join(sorted(q.route))))
    study=pd.DataFrame(rows,columns=['axis_id','species','study_id','setting','gene_id','gene_symbol','APA','alo','ahi','GEX','glo','ghi','meas','depth','PAS','reference_WPP','reference_WPP_min','reference_WPP_max','routes'])
    paired=[]
    for sp,ids in STUDIES.items():
        q=study[study.species.eq(sp)&study.axis_id.eq(AXES[sp])];a=q[q.study_id.eq(ids[0])];b=q[q.study_id.eq(ids[1])]
        z=a.merge(b,on=['axis_id','species','setting','gene_id','gene_symbol'],suffixes=('_1','_2'),validate='one_to_one')
        z['guan_series_pattern']='not_applicable_mouse'
        if sp=='human':
            def pattern(routes):
                x=set(routes.split(';'));a='178325_0618' in x;b='178325_0809' in x
                return 'both_measurable' if a and b else '0618_only' if a else '0809_only' if b else 'neither_prespecified_series'
            z['guan_series_pattern']=z.routes_2.map(pattern)
            assert not z.guan_series_pattern.eq('neither_prespecified_series').any()
        paired.append(z)
    ctx.table(r,'shared_route_covariates.tsv.gz');ctx.table(study,'shared_study_covariates.tsv.gz')
    frame=pd.concat(paired,ignore_index=True);ctx.table(frame,'permutation_input.tsv.gz');return frame


def permutations(frame,sp,t,model,bins,draws=10000):
    q=frame[frame.species.eq(sp)].reset_index(drop=True);labels=[sign(q['alo_'+s],q['ahi_'+s],t) for s in ['1','2']]
    gx=[np.where(q['meas_'+s],sign(q['glo_'+s],q['ghi_'+s],.25),0) for s in ['1','2']];g=gx[0]*gx[1]
    observed=four_counts(*labels,g);rows=[];strata=[]
    for shuffled in [1,2]:
        suffix=str(shuffled);lab=labels[shuffled-1];lo=q['reference_WPP_min_'+suffix].ge(t-1e-9);hi=q['reference_WPP_max_'+suffix].le(1-t+1e-9)
        assert not ((lab==-1)&~lo).any() and not ((lab==1)&~hi).any()
        f=pd.DataFrame(dict(gex_response=pd.cut(q['GEX_'+suffix],[-np.inf,-.25,.25,np.inf],labels=['down','flat','up']).astype(str),
            depth_bin=qcodes(q['depth_'+suffix],3),pas_cap=q['PAS_'+suffix].clip(upper=4),baseline=qcodes(q['reference_WPP_'+suffix],bins),
            can_negative=lo,can_positive=hi,series=q.guan_series_pattern))
        if model=='direction_specific':f['responder']=lab!=0
        key=f.astype(str).agg('|'.join,axis=1) if len(f) else pd.Series(dtype=str)
        groups=list(f.groupby(key).indices.values());rng=np.random.default_rng(seed('apa_gex',sp,t,model,bins,shuffled))
        null=np.empty((draws,4),int);movable=0
        for ix in groups:
            count=np.array([(lab[ix]==v).sum() for v in [-1,0,1]]);p=count/len(ix);mixed=np.count_nonzero(count)>1;movable+=len(ix)*mixed
            strata.append(dict(species=sp,threshold=t,null_model=model,baseline_bins=bins,shuffled_study=shuffled,
                stratum=key.iloc[ix[0]],n_genes=len(ix),negative=count[0],non_event=count[1],positive=count[2],mixed_labels=mixed,
                expected_labels_changed=len(ix)*(1-p@p)))
        for b in range(draws):
            shuffled_labels=lab.copy()
            for ix in groups:shuffled_labels[ix]=rng.permutation(lab[ix])
            null[b]=four_counts(shuffled_labels,labels[1],g) if shuffled==1 else four_counts(labels[0],shuffled_labels,g)
        if model=='direction_specific':
            assert (null[:,0]+null[:,2]==observed[0]+observed[2]).all() and (null[:,1]+null[:,3]==observed[1]+observed[3]).all()
        for i,cell in enumerate('ABCD'):
            mean=null[:,i].mean();rows.append(dict(species=sp,axis_id=AXES[sp],threshold=t,null_model=model,baseline_bins=bins,shuffled_study=shuffled,
                cell=cell,observed=observed[i],null_mean=mean,null_p025=np.quantile(null[:,i],.025),null_p975=np.quantile(null[:,i],.975),
                enrichment=observed[i]/mean if mean else np.nan,p_enrich=(1+(null[:,i]>=observed[i]).sum())/(draws+1),
                p_deplete=(1+(null[:,i]<=observed[i]).sum())/(draws+1),comparable_genes=len(q),shuffleable_genes=movable,permutations=draws,
                status='measured' if len(q) else 'empty_background'))
    return rows,strata


def candidate_relations(ctx):
    effects=ctx.read(ctx.analysis/'intermediate/shared/matched_route_effects.tsv.gz');native=ctx.read(ctx.analysis/'recurrence/recurrence_members.tsv.gz')
    members=native[native.metric.eq('WPP')&native.threshold.eq(.02)&native.series_matching.astype(bool)].drop_duplicates(['species','contrast','gene_id'])
    axes_by_species={sp:q.axis_id.unique() for sp,q in effects.groupby('species',sort=False)}
    settings_by_axis={axis:q.setting.unique() for axis,q in effects.groupby('axis_id',sort=False)}
    gene_rows=effects.groupby(['axis_id','setting','gene_id'],sort=False).indices
    rows=[]
    for rec in members.itertuples(index=False):
        direction=1 if rec.direction=='positive' else -1
        for axis in axes_by_species.get(rec.species,[]):
            for setting in settings_by_axis[axis]:
                q=effects.iloc[gene_rows.get((axis,setting,rec.gene_id),[])]
                for contrast in ['P_S','P_M','R_S','E_S','M_S','E_R','P_E','P_R','E_M','R_M']:
                    design='true_source' if contrast.endswith('_S') else 'same_library';z=q[q.contrast.eq(contrast)&q.design.eq(design)]
                    by=z.groupby('study_id').agg(lo=('APA_delta','min'),hi=('APA_delta','max'),glo=('GEX_delta','min'),ghi=('GEX_delta','max'),meas=('GEX_all_pairs_measurable','all'))
                    support=by.lo.ge(.02) if direction==1 else by.hi.le(-.02);opp=by.hi.le(-.02) if direction==1 else by.lo.ge(.02)
                    rows.append(dict(species=rec.species,entrance=rec.contrast,gene_id=rec.gene_id,gene_symbol=rec.gene_symbol,axis_id=axis,setting=setting,
                        contrast=contrast,design=design,direction=direction,evaluable_studies=len(by),supporting_studies=int(support.sum()),opposite_studies=int(opp.sum()),
                        recurrent=bool(support.sum()>=2),GEX_opposite=bool(len(by)>=2 and by.meas.all() and by.glo.ge(.25).any() and by.ghi.le(-.25).any()),
                        source_interpretation_allowed=rec.contrast=='P_S',status='measured' if len(by)>=2 else 'fewer_than_two_studies'))
    ctx.table(rows,'current_candidate_shared_state_relations.tsv.gz',columns=['species','entrance','gene_id','gene_symbol','axis_id','setting','contrast','design','direction','evaluable_studies','supporting_studies','opposite_studies','recurrent','GEX_opposite','source_interpretation_allowed','status'])


def permutation_task(ctx,item):
    frame=ctx.read(ctx.out/'permutation_input.tsv.gz')
    return permutations(frame,*item)


def run(ctx):
    frame=inputs(ctx);rows=[];diagnostics=[]
    tasks=itertools.product(['human','mouse'],[.02,.05],['response_identity','direction_specific'],[5,3])
    for a,b in ordered_map(ctx,permutation_task,tasks):rows.extend(a);diagnostics.extend(b)
    four=pd.DataFrame(rows);ctx.table(four,'four_class_nulls.tsv.gz')
    ctx.table(diagnostics,'exchangeability.tsv.gz',columns=['species','threshold','null_model','baseline_bins','shuffled_study','stratum','n_genes','mixed_labels'])
    tests=four[four.cell.eq('B')].groupby(['species','axis_id','threshold','null_model','baseline_bins']).agg(observed=('observed','first'),
        p_both=('p_enrich','max'),enrichment_min=('enrichment','min'),comparable_genes=('comparable_genes','first'),shuffleable_genes_min=('shuffleable_genes','min')).reset_index()
    assert tests.groupby('baseline_bins').size().eq(8).all();tests['q_family']=tests.groupby('baseline_bins').p_both.transform(bh)
    ctx.table(tests,'prespecified_eight_test_family.tsv');candidate_relations(ctx)
    ctx.finish(dict(new_common_PAS_and_MC_background=True,all_evaluable_Guan_series=True,mouse_null_recomputed=True,
        response_identity_and_direction_distinct=True,larger_of_two_shuffle_P=True,primary_and_q3_eight_test_families=True,
        candidate_states_retraced=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);a=p.parse_args();run(Context(a.config,'apa_gex',a.out))
