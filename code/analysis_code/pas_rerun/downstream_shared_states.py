"""Fresh cross-study PAS axes, common-state expression and conditional models."""
import hashlib
import itertools
import json

import numpy as np
import pandas as pd

from downstream_runtime import Study,local_state_labels

BASE=['axis_id','species','study_id','route','setting','gene_id','gene_symbol']
CONTRASTS=[('P','M'),('P','S'),('R','S'),('E','S'),('M','S'),('E','R'),('P','E'),('P','R'),('E','M'),('R','M')]
SOURCE='github/reproapa-analysis/modules/08_recurrence/code/apa_gex_state/analyze.py'
VARIANTS={'shared_p08':(.8,20,{}),'shared_p09':(.9,20,{}),'shared_depth40':(.8,40,{}),
    'exclude_C012':(.8,20,{'GSE103221':('E','C012')}),'exclude_L00':(.8,20,{'GSE115943':('R','L00')}),
    'exclude_both':(.8,20,{'GSE103221':('E','C012'),'GSE115943':('R','L00')}),
    'support_q95':(.8,20,{}),'support_q90':(.8,20,{})}
MARKERS=[['MCM5','PCNA','TYMS','FEN1','MCM2','MCM4','RRM1','UNG','GINS2','MCM6'],['CDK1','TOP2A','NUSAP1','CDC20','CDC45','UBE2C','BIRC5','TPX2']]


def clusters(axes,studies):
    f=pd.concat(axes,ignore_index=True)
    if f.chrom.nunique()!=1 or f.strand.nunique()!=1:return []
    groups=[]
    for row in f.sort_values(['cut_site','study_id','atlas_pas_id']).to_dict('records'):
        if not groups or row['cut_site']-groups[-1][0]['cut_site']>50:groups.append([])
        groups[-1].append(row)
    return [g for g in groups if {r['study_id'] for r in g}==set(studies)]


def conditional(x,values):
    rank=np.linalg.matrix_rank(x);condition=np.linalg.cond(x)
    good=rank==x.shape[1] and condition<1e6 and len(x)>x.shape[1]
    return dict(adjusted=np.linalg.lstsq(x,values,rcond=None)[0][1] if good else np.nan,rank=rank,condition=condition,
        residual_df=len(x)-rank,status='estimated' if good else 'rank_condition_or_df_failure')


def state_groups(meta,sid,local,support):
    """Freeze metadata-only groups; gene-specific PAS depth is applied afterwards."""
    groups={}
    for setting,(purity,depth,exclusions) in VARIANTS.items():
        mask=meta.purity.ge(purity)&meta.node.isin(list('SREPM'))&~meta.route.str.contains('reference')
        if sid in exclusions:
            state,local_id=exclusions[sid];mask &= ~(meta.node.eq(state)&(local==local_id))
        if setting.startswith('support_q'):
            mask &= meta.metacell_id.map(support[setting]).to_numpy(dtype=bool,na_value=False)
        groups[setting]=[(key,part.index.to_numpy()) for key,part in meta[mask].groupby(['route','library_key','node'])]
    return groups


def extract(ctx):
    studies={sid:Study(ctx,sid) for sid in ['GSE103221','GSE115943','GSE118258','GSE150311','GSE178325','GSE221739']}
    cycles={};local={};designs={};groups_by_setting={};cell_counts={};metacell_ids={}
    support=ctx.read(ctx.resource('state_expression_support')).set_index('metacell_id')
    assert support.index.is_unique
    for sid,s in studies.items():
        s.load_gex();symbols=s.gex_symbols().str.upper().str.split(';')
        cols=[]
        for markers in MARKERS:
            keep=symbols.map(lambda x:bool(set(x)&set(markers))).to_numpy();assert keep.sum()>=5
            cols.append(np.log2(1+s.gex_counts[:,keep].toarray()/s.gex_total[:,None]*1e6).mean(1))
        cycles[sid]=np.column_stack(cols)
        if sid in ['GSE103221','GSE115943']:
            full=ctx.read(s.root/'metacells.tsv.gz');local[sid]=local_state_labels(full,s.meta.metacell_id)
        designs[sid]=s.design[['route','library_key','time','sample_unit','replicate_status']].drop_duplicates().set_index(['route','library_key']).to_dict('index')
        groups_by_setting[sid]=state_groups(s.meta,sid,local.get(sid),support)
        cell_counts[sid]=s.meta.n_cells.to_numpy();metacell_ids[sid]=s.meta.metacell_id.to_numpy()
    rows=[];regions=[];coverage=[];adjust=[]
    for sp,sids in [('mouse',['GSE103221','GSE115943']),('human',['GSE118258','GSE150311','GSE178325','GSE221739'])]:
        for pair in itertools.combinations(sids,2):
            axis_id='__'.join(pair);shared=set(studies[pair[0]].genes)&set(studies[pair[1]].genes)
            for gid in sorted(shared):
                axes=[studies[sid].axis(gid) for sid in pair];groups=clusters(axes,pair)
                coverage.append(dict(axis_id=axis_id,species=sp,gene_id=gid,n_shared_regions=len(groups)))
                if len(groups)<2:continue
                cuts=np.array([np.mean([r['cut_site'] for r in g]) for g in groups]);signed=cuts*(1 if axes[0].strand.iloc[0]=='+' else -1)
                if np.ptp(signed)==0:continue
                positions=(signed-signed.min())/np.ptp(signed);symbol=axes[0].primary_gene_symbol.iloc[0]
                for sid in pair:
                    s=studies[sid]
                    if gid not in s.gex_index:continue
                    counts=np.column_stack([np.asarray(s.counts[:,[r['feature_index_1based']-1 for r in g if r['study_id']==sid]].sum(1)).ravel() for g in groups])
                    totals=counts.sum(1);apa=np.divide(counts@positions,totals,out=np.zeros(len(totals)),where=totals>0)
                    gcounts=s.gex_counts[:,s.gex_index.get_loc(gid)].toarray().ravel();gl=np.log2(1+gcounts/s.gex_total*1e6)
                    for i,g in enumerate(groups):
                        for r in g:
                            if r['study_id']==sid:regions.append(dict(axis_id=axis_id,species=sp,study_id=sid,gene_id=gid,region=i,
                                chrom=r['chrom'],strand=r['strand'],cut_site=r['cut_site'],atlas_pas_id=r['atlas_pas_id'],axis=positions[i]))
                    for setting,(purity,depth,exclusions) in VARIANTS.items():
                        if sp=='human' and exclusions:continue
                        accepted={}
                        for (route,lib,st),allrows in groups_by_setting[sid][setting]:
                            ii=allrows[totals[allrows]>=depth]
                            n_cells=cell_counts[sid][ii].sum()
                            if len(ii)<3 or n_cells<30:continue
                            d=designs[sid][route,lib];accepted[route,lib,st]=ii
                            rows.append(dict(axis_id=axis_id,species=sp,study_id=sid,route=route,setting=setting,gene_id=gid,gene_symbol=symbol,
                                library_key=lib,state=st,**d,APA=apa[ii].mean(),GEX=gl[ii].mean(),GEX_all=gl[allrows].mean(),
                                GEX_count=int(gcounts[ii].sum()),GEX_CPM=gcounts[ii].sum()/s.gex_total[ii].sum()*1e6,
                                n_mc=len(ii),n_cells=int(n_cells),pas_counts=int(totals[ii].sum()),
                                site_counts=json.dumps(counts[ii].sum(0).tolist()),site_mc=json.dumps((counts[ii]>0).sum(0).tolist()),
                                usage=json.dumps((counts[ii]/totals[ii,None]).mean(0).tolist()),
                                matched_MC_SHA256=hashlib.sha256(';'.join(metacell_ids[sid][ii]).encode()).hexdigest()))
                        if setting!='shared_p08':continue
                        for route,lib in sorted({k[:2] for k in accepted}):
                            for a,b in CONTRASTS:
                                if (route,lib,a) not in accepted or (route,lib,b) not in accepted:continue
                                ia,ib=accepted[route,lib,a],accepted[route,lib,b];ii=np.r_[ia,ib]
                                if ((counts[ii].sum(0)>=10)&((counts[ii]>0).sum(0)>=2)).sum()<2:continue
                                def z(x):return (x-x.mean(0))/np.maximum(x.std(0),1e-8)
                                state=np.r_[np.ones(len(ia)),np.zeros(len(ib))];cy=z(cycles[sid][ii]);ge=z(gl[ii]);dp=z(np.log1p(totals[ii]))
                                for model in ['cycle','cycle_and_gex','cycle_gex_depth']:
                                    columns=[np.ones(len(ii)),state,cy[:,0],cy[:,1]]
                                    if model!='cycle':columns.append(ge)
                                    if model=='cycle_gex_depth':columns.append(dp)
                                    adjust.append(dict(axis_id=axis_id,species=sp,study_id=sid,route=route,library_key=lib,
                                        **designs[sid][route,lib],gene_id=gid,gene_symbol=symbol,contrast=a+'_'+b,model=model,
                                        unadjusted=apa[ia].mean()-apa[ib].mean(),n_mc=len(ii),**conditional(np.column_stack(columns),apa[ii])))
            print('shared_axis',axis_id,'rows',len(rows),flush=True)
    state=pd.DataFrame(rows,columns=BASE+['library_key','state','time','sample_unit','replicate_status','APA','GEX','GEX_all','GEX_count','GEX_CPM','n_mc','n_cells','pas_counts','site_counts','site_mc','usage','matched_MC_SHA256'])
    ctx.table(state,'shared/state_values.tsv.gz')
    ctx.table(regions,'shared/PAS_regions.tsv.gz',columns=['axis_id','species','study_id','gene_id','region','chrom','strand','cut_site','atlas_pas_id','axis'])
    ctx.table(coverage,'shared/PAS_coverage.tsv.gz',columns=['axis_id','species','gene_id','n_shared_regions'])
    ctx.table(adjust,'shared/conditional_sensitivity.tsv.gz',columns=['axis_id','species','study_id','route','library_key','time','sample_unit','gene_id','contrast','model','unadjusted','adjusted','rank','condition','n_mc'])
    return state


def compare(ctx,state):
    f=ctx.functions(SOURCE,['pair','collapse_window'],dict(keys=BASE,values=['APA','GEX','GEX_all']))
    parts=[]
    for axis,w0 in state.groupby('axis_id'):
        for a,b in CONTRASTS:
            z=f.pair(w0[w0.state.eq(a)],w0[w0.state.eq(b)],BASE+['library_key','time','sample_unit','replicate_status'],a+'_'+b,'same_library')
            if len(z):parts.append(z)
        for route,w in w0.groupby('route'):
            source=w[w.state.eq('S')&w.time.astype(str).isin(['D0','source','Stage_I_day_0'])];design='true_source'
            if source.empty and route=='GSE221739':source=w[w.state.eq('S')&w.time.astype(str).eq('D3')];design='induced_D3_reference'
            if source.empty:continue
            for target in ['R','E','P','M']:
                z=f.pair(w[w.state.eq(target)],source,BASE,target+'_S',design)
                if len(z):parts.append(z)
            for b in ['E','R']:
                inter=w[w.state.eq(b)]
                if inter.empty or w[w.state.eq('P')].empty or set(inter.library_key)&set(w[w.state.eq('P')].library_key):continue
                z=f.pair(w[w.state.eq('P')],f.collapse_window(inter),BASE,'P_'+b,'target_vs_intermediate_window')
                if len(z):parts.append(z)
    cols=BASE+['contrast','design','n_joint_PAS','GEX_measurable']+[v+'_delta' for v in ['APA','GEX','GEX_all']]+[c+s for c in ['library_key','time','sample_unit','replicate_status','n_mc','n_cells','pas_counts'] for s in ['_a','_b']]
    pairs=pd.concat(parts,ignore_index=True) if parts else pd.DataFrame(columns=cols);ctx.table(pairs,'shared/matched_library_effects.tsv.gz')
    k=BASE+['contrast','design'];values=['APA_delta','GEX_delta','GEX_all_delta']
    units=pairs.groupby(k+['time_a','sample_unit_a','time_b','sample_unit_b'],dropna=False)[values].mean()
    units=units.groupby(k+['time_a','sample_unit_a','time_b']).mean().groupby(k+['time_a','sample_unit_a']).mean()
    dates=units.groupby(k+['time_a']).mean().reset_index();effects=dates.groupby(k)[values].mean().reset_index()
    effects=effects.merge(pairs.groupby(k).GEX_measurable.all().rename('GEX_all_pairs_measurable').reset_index(),on=k,validate='one_to_one')
    ctx.table(dates,'shared/matched_date_effects.tsv.gz');ctx.table(effects,'shared/matched_route_effects.tsv.gz')
    # Triangle comparisons share the physical target library, source library and all PAS-qualified contributors.
    keys=BASE+['library_key_a','library_key_b','time_a','time_b','sample_unit_a','sample_unit_b']
    ps=pairs[pairs.design.eq('true_source')&pairs.contrast.eq('P_S')];ms=pairs[pairs.design.eq('true_source')&pairs.contrast.eq('M_S')]
    tri=ps.merge(ms,on=keys,suffixes=('_PS','_MS'),validate='one_to_one')
    pm=pairs[pairs.design.eq('same_library')&pairs.contrast.eq('P_M')]
    tri=tri.merge(pm[BASE+['library_key_a','GEX_measurable']+values],on=BASE+['library_key_a'],validate='many_to_one')
    for val in values:
        tri[val+'_triangle_error']=tri[val+'_PS']-tri[val+'_MS']-tri[val]
        assert np.allclose(tri[val+'_triangle_error'],0,rtol=0,atol=1e-12)
    ctx.table(tri,'shared/common_window_triangle.tsv.gz')
    return effects
