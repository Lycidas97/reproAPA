"""Frozen contributor matched descriptive cycle/source-state attribution."""
from pathlib import Path
import ast
import hashlib
import json
import numpy as np
import pandas as pd
H=Path(__file__).resolve().parent; R=H.parents[2]
B=R/'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables'
SETTING='supported50_p08'
SOURCES={}
WEIGHT_CACHE={}

def read(p,**kwargs):
    SOURCES[str(p.relative_to(R))]=hashlib.sha256(p.read_bytes()).hexdigest()
    return pd.read_csv(p,sep='\t',float_precision='round_trip',**kwargs)

def save(f,n):
    pd.DataFrame(f).to_csv(H/n,sep='\t',index=False,na_rep='NA',float_format='%.17g',compression={'method':'gzip','mtime':0} if n.endswith('.gz') else None)

def weights(f):
    key=tuple(zip(f.time_a,f.sample_unit_a,f.time_b,f.sample_unit_b))
    if key in WEIGHT_CACHE:return WEIGHT_CACHE[key]
    levels=[['time_a','sample_unit_a','time_b','sample_unit_b'],['time_a','sample_unit_a','time_b'],['time_a','sample_unit_a'],['time_a'],[]]
    q=f.reset_index(drop=True)
    out=1/q.groupby(levels[0],dropna=False).time_a.transform('size').to_numpy(float)
    for child,parent in zip(levels[:-1],levels[1:]):
        distinct=q[child].drop_duplicates()
        if not parent: out/=len(distinct)
        else:
            n=distinct.groupby(parent,dropna=False).size().rename('_n').reset_index()
            out/=q[parent].merge(n,on=parent,how='left',validate='many_to_one')._n.to_numpy()
    assert abs(out.sum()-1)<1e-12
    WEIGHT_CACHE[key]=out
    return out

KEY=['species','study_id','route','contrast','design','gene_id']
VAR=['exact_Cycle','exact_Adjusted_Cycle','cohort_Cycle','cohort_Adjusted_Cycle']

def cycle_pairs(p,s):
    lut={(row.study_id,row.metacell_id):(row.Cycle,row.Adjusted_Cycle) for row in s.itertuples()}
    full=s.groupby(['study_id','route','library_key','node'])[['Cycle','Adjusted_Cycle']].mean()
    cache={}
    for side in ['a','b']:
        values=[]; cohorts=[]
        for row in p.itertuples():
            ids=getattr(row,'metacell_ids_'+side)
            key=(row.study_id,ids)
            if key not in cache:
                cache[key]=np.mean([lut[(row.study_id,i)] for i in json.loads(ids)],axis=0)
            values.append(cache[key])
            node='P' if side=='a' else row.contrast[-1]
            cohorts.append(full.loc[(row.study_id,row.route,getattr(row,'library_key_'+side),node)].to_numpy(float))
        values=np.asarray(values);cohorts=np.asarray(cohorts)
        for idx,col in enumerate(['Cycle','Adjusted_Cycle']):
            p['exact_'+col+'_'+side]=values[:,idx];p['cohort_'+col+'_'+side]=cohorts[:,idx]
    for col in VAR:p[col+'_delta']=p[col+'_a']-p[col+'_b']
    return p,full,lut

def collapse(p):
    rows=[];err=0
    for key,q in p.groupby(KEY,sort=False):
        w=weights(q);err=max(err,float(np.abs(w-q.weight.to_numpy()).max()))
        row=dict(zip(KEY,key))|dict(n_pairs=len(q),n_target_libraries=q.library_key_a.nunique(),n_source_libraries=q.library_key_b.nunique(),delta_DWUI=float(w@q.delta_DWUI.to_numpy()))
        row.update({col:float(w@q[col].to_numpy()) for col in p if col.startswith(('exact_','cohort_'))})
        rows.append(row)
    assert err<1e-12
    return pd.DataFrame(rows),err

def summarize(g,o):
    rows=[]
    universes={}
    for co in ['P_S','P_M']:
        q=g[g.contrast.eq(co)]
        common=o.copy()
        for (sp,rt),part in q.groupby(['species','route']):
            common=common[common[sp+'_gene_id'].isin(part.gene_id)]
        universes[co]=common
        save(common,'common_'+co+'_orthologs.tsv')
    for key,q in g.groupby(['species','study_id','route','contrast','design'],sort=True):
        sp,st,rt,co,d=key
        for label,part in [('native',q),('common_ortholog',q[q.gene_id.isin(universes[co][sp+'_gene_id'])])]:
            assert len(part)
            context=dict(zip(['species','study_id','route','contrast','design'],key))|dict(background=label,genes=len(part),mean_DWUI=float(part.delta_DWUI.mean()),median_DWUI=float(part.delta_DWUI.median()),proximal=int(part.delta_DWUI.le(-.02).sum()),distal=int(part.delta_DWUI.ge(.02).sum()))
            for col in VAR:
                rows.append(context|dict(cycle_variant=col,mean_Cycle_delta=float(part[col+'_delta'].mean()),mean_Cycle_target=float(part[col+'_a'].mean()),mean_Cycle_reference=float(part[col+'_b'].mean())))
    return pd.DataFrame(rows),universes

def fit_summary(frame):
    fits=[];residuals=[];loso=[]
    for (background,variant),f in frame[frame.contrast.eq('P_S')].groupby(['background','cycle_variant']):
        f=f.sort_values('route').copy();f['study_weight']=np.where(f.study_id.eq('GSE178325'),.5,1.)
        x=f.mean_Cycle_delta.to_numpy(float);y=f.mean_DWUI.to_numpy(float);mouse=f.species.eq('mouse').to_numpy(float);w=f.study_weight.to_numpy(float)
        assert len(f)==5 and len(set(f.study_id))==4 and abs(w.sum()-4)<1e-12
        groups={sp:float(np.average(f.loc[f.species.eq(sp),'mean_DWUI'],weights=f.loc[f.species.eq(sp),'study_weight'])) for sp in ['human','mouse']}
        group_cycle={sp:float(np.average(f.loc[f.species.eq(sp),'mean_Cycle_delta'],weights=f.loc[f.species.eq(sp),'study_weight'])) for sp in ['human','mouse']}
        for model in ['cycle_only','cycle_plus_mouse']:
            X=np.column_stack([np.ones(len(x)),x]+([mouse] if model=='cycle_plus_mouse' else []))
            beta,_,rank,sval=np.linalg.lstsq(X*np.sqrt(w)[:,None],y*np.sqrt(w),rcond=None)
            predicted=X@beta;err=y-predicted;ym=np.average(y,weights=w);sse=float(w@(err*err));sst=float(w@((y-ym)**2))
            row=dict(background=background,cycle_variant=variant,model=model,n_routes=5,n_studies=4,total_study_weight=4.,rank=int(rank),parameters=X.shape[1],residual_df=len(f)-int(rank),condition_number=float(sval[0]/sval[-1]),intercept=float(beta[0]),cycle_beta=float(beta[1]),mouse_beta=float(beta[2]) if len(beta)>2 else np.nan,weighted_R2=1-sse/sst,weighted_RMSE=np.sqrt(sse/w.sum()),unadjusted_mouse_minus_human=groups['mouse']-groups['human'],mouse_minus_human_Cycle=group_cycle['mouse']-group_cycle['human'],p=np.nan)
            fits.append(row)
            for i,rec in enumerate(f.to_dict('records')):residuals.append(rec|dict(model=model,predicted_DWUI=float(predicted[i]),residual_DWUI=float(err[i])))
            for study in sorted(set(f.study_id)):
                mask=f.study_id.ne(study).to_numpy();XX=X[mask];yy=y[mask];ww=w[mask]
                bb,_,rr,ss=np.linalg.lstsq(XX*np.sqrt(ww)[:,None],yy*np.sqrt(ww),rcond=None)
                identified=rr==X.shape[1]
                for i in np.flatnonzero(~mask):
                    pred=float(X[i]@bb) if identified else np.nan
                    loso.append(dict(background=background,cycle_variant=variant,model=model,omitted_study=study,heldout_route=f.iloc[i].route,training_routes=int(mask.sum()),training_studies=len(set(f.loc[mask,'study_id'])),rank=int(rr),parameters=X.shape[1],training_residual_df=int(mask.sum())-int(rr),zero_residual_df=int(mask.sum())==int(rr),identified=identified,condition_number=float(ss[0]/ss[-1]),intercept=float(bb[0]) if identified else np.nan,cycle_beta=float(bb[1]) if identified else np.nan,mouse_beta=float(bb[2]) if identified and len(bb)>2 else np.nan,predicted_DWUI=pred,error_DWUI=float(y[i]-pred),heldout_Cycle=float(x[i]),training_Cycle_min=float(x[mask].min()),training_Cycle_max=float(x[mask].max()),cycle_extrapolation=bool(x[i]<x[mask].min() or x[i]>x[mask].max())))
    return pd.DataFrame(fits),pd.DataFrame(residuals),pd.DataFrame(loso)

def triangles(p,s,full,lut):
    t=read(B/'matched_pms/weighted_triplets.tsv.gz',usecols=['species','study_id','route','setting','gene_id','library_key_a','library_key_b','time_a','time_b','sample_unit_a','sample_unit_b','metacell_ids_a','metacell_ids_b','triplet_weight'])
    t=t[t.setting.eq(SETTING)].copy()
    pm=p[p.contrast.eq('P_M')].set_index(['study_id','route','gene_id','library_key_a'])
    assert pm.index.is_unique
    t=t.merge(pm[['metacell_ids_b','exact_Cycle_b','exact_Adjusted_Cycle_b','cohort_Cycle_b','cohort_Adjusted_Cycle_b']].rename(columns={'metacell_ids_b':'M_metacell_ids','exact_Cycle_b':'exact_Cycle_M','exact_Adjusted_Cycle_b':'exact_Adjusted_Cycle_M','cohort_Cycle_b':'cohort_Cycle_M','cohort_Adjusted_Cycle_b':'cohort_Adjusted_Cycle_M'}),left_on=['study_id','route','gene_id','library_key_a'],right_index=True,validate='many_to_one',how='left')
    assert t.M_metacell_ids.notna().all()
    cache={}
    for state,side in [('P','a'),('S','b')]:
        v=[];c=[]
        for row in t.itertuples():
            key=(row.study_id,getattr(row,'metacell_ids_'+side))
            if key not in cache:cache[key]=np.mean([lut[(row.study_id,i)] for i in json.loads(key[1])],axis=0)
            v.append(cache[key]);c.append(full.loc[(row.study_id,row.route,getattr(row,'library_key_'+side),state)].to_numpy(float))
        v=np.array(v);c=np.array(c)
        for idx,col in enumerate(['Cycle','Adjusted_Cycle']):t['exact_'+col+'_'+state]=v[:,idx];t['cohort_'+col+'_'+state]=c[:,idx]
    rows=[];weight_error=0
    for key,q in t.groupby(['species','study_id','route','gene_id'],sort=False):
        w=weights(q);weight_error=max(weight_error,float(np.abs(w-q.triplet_weight.to_numpy()).max()))
        row=dict(zip(['species','study_id','route','gene_id'],key))|dict(n_triplets=len(q))
        for col in VAR:
            for st in ['P','M','S']:row[col+'_'+st]=float(w@q[col+'_'+st].to_numpy())
            for co,a,b in [('PS','P','S'),('MS','M','S'),('PM','P','M')]:row[col+'_'+co]=row[col+'_'+a]-row[col+'_'+b]
        rows.append(row)
    assert weight_error<1e-12
    g=pd.DataFrame(rows)
    m=read(B/'matched_pms/matched_state_vectors.tsv.gz',usecols=['species','study_id','route','setting','gene_id','DWUI_PS','DWUI_MS','DWUI_PM']);m=m[m.setting.eq(SETTING)].drop(columns='setting')
    g=g.merge(m,on=['species','study_id','route','gene_id'],validate='one_to_one');assert len(g)==len(m)
    ids=read(R/'notes/2026-10-10-route-step-decomposition/three_route_ortholog_gene_effects.tsv.gz',usecols=['human_gene_id','mouse_gene_id','setting'])
    ids=ids[ids.setting.eq(SETTING)]
    assert len(ids)==678 and ids.human_gene_id.is_unique and ids.mouse_gene_id.is_unique
    summaries=[]
    for (sp,st,rt),q in g.groupby(['species','study_id','route']):
        for label,z in [('native',q),('matched_678',q[q.gene_id.isin(ids[sp+'_gene_id'])])]:
            if label=='matched_678':assert len(z)==678
            for col in VAR:
                for co in ['PS','MS','PM']:summaries.append(dict(species=sp,study_id=st,route=rt,background=label,genes=len(z),contrast=co,cycle_variant=col,mean_DWUI=float(z['DWUI_'+co].mean()),mean_Cycle_delta=float(z[col+'_'+co].mean())))
    summary=pd.DataFrame(summaries);cor=[]
    for (bg,v),z in summary.groupby(['background','cycle_variant']):
        cor.append(dict(background=bg,cycle_variant=v,n_algebraically_dependent_points=len(z),independent_routes=3,Pearson_r=float(np.corrcoef(z.mean_Cycle_delta,z.mean_DWUI)[0,1]),p=np.nan,scope='descriptive_only_PS_equals_MS_plus_PM'))
    save(g,'matched_triangle_gene_cycles.tsv.gz');save(summary,'matched_triangle_summary.tsv');save(cor,'matched_triangle_descriptive_correlation.tsv')
    cycle_error=max(float(np.abs(g[v+'_PS']-g[v+'_MS']-g[v+'_PM']).max()) for v in VAR)
    apa_error=float(np.abs(g.DWUI_PS-g.DWUI_MS-g.DWUI_PM).max())
    assert max(cycle_error,apa_error)<1e-12
    return dict(triplet_weight_max_error=weight_error,triangle_Cycle_additivity_max_error=cycle_error,triangle_APA_additivity_max_error=apa_error)

def main():
    protocol=R/'notes/2026-10-10-cycle-source-attribution/PROTOCOL.zh.md'
    protocol_hash=hashlib.sha256(protocol.read_bytes()).hexdigest()
    authorization=json.loads((H/'authorization.json').read_text());assert authorization['protocol_sha256']==protocol_hash
    s=read(R/'notes/2026-10-10-proliferation-apa/metacell_cycle_scores.tsv.gz')
    cols=['species','study_id','route','contrast','design','setting','gene_id','library_key_a','library_key_b','time_a','time_b','sample_unit_a','sample_unit_b','metacell_ids_a','metacell_ids_b','n_mc_a','n_mc_b','weight','delta_DWUI']
    p=read(B/'state_vectors/pair_effects.primary.tsv.gz',usecols=cols);p=p[p.setting.eq(SETTING)].copy()
    p,full,lut=cycle_pairs(p,s);g,werr=collapse(p)
    frozen=read(B/'state_vectors/route_effects.primary.tsv.gz',usecols=KEY+['setting','delta_DWUI']);frozen=frozen[frozen.setting.eq(SETTING)]
    audit=g.merge(frozen,on=KEY,validate='one_to_one',suffixes=['_replay','_frozen']);assert len(audit)==len(g)==len(frozen)
    err=float(np.abs(audit.delta_DWUI_replay-audit.delta_DWUI_frozen).max());assert err<1e-12
    o=read(B/'cross_overlap/ortholog_matched_features.tsv.gz',usecols=['human_gene_id','mouse_gene_id','homology_type','is_high_confidence','human_partner_count','mouse_partner_count'])
    assert o.human_gene_id.is_unique and o.mouse_gene_id.is_unique
    assert o.homology_type.eq('ortholog_one2one').all() and o.is_high_confidence.eq(1).all() and o.human_partner_count.eq(1).all() and o.mouse_partner_count.eq(1).all()
    summary,univ=summarize(g,o[['human_gene_id','mouse_gene_id']]);assert len(univ['P_S'])==1499 and len(univ['P_M'])==463
    fits,residuals,loso=fit_summary(summary)
    save(p,'matched_pair_cycles.tsv.gz');save(g,'route_gene_cycles.tsv.gz');save(summary,'route_contrast_cycle_summary.tsv');save(fits,'source_PS_descriptive_models.tsv');save(residuals,'source_PS_model_residuals.tsv');save(loso,'source_PS_leave_one_study_out.tsv')
    tri=triangles(p,s,full,lut)
    validation=dict(status='PASS',protocol_sha256=protocol_hash,pair_rows=len(p),route_genes=len(g),route_effect_replay_max_error=err,weight_replay_max_error=werr,true_source_PS_common=1499,same_library_PM_common=463,all_score_members_covered=True,no_p_or_species_population_CI=True,observed_source_augmentation='pending_sm_core_output',**tri)
    (H/'validation.json').write_text(json.dumps(validation,indent=2)+'\n');(H/'input_sha256.json').write_text(json.dumps(SOURCES,indent=2)+'\n')
    print(summary[(summary.background.eq('common_ortholog'))&summary.cycle_variant.eq('exact_Cycle')].to_string(index=False));print(fits[(fits.background.eq('common_ortholog'))&fits.cycle_variant.eq('exact_Cycle')].to_string(index=False))

if __name__=='__main__':main()
