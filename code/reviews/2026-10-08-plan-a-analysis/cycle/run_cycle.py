"""Fixed-cohort, common-cycle-bin sensitivity on frozen PAS counts."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
BASES={
    'polyapipe':ROOT/'notes/2026-10-06-supplement-polyapipe/contract-r2/bundle/tables',
    'apadetect_p50':ROOT/'notes/2026-10-06-dual-main-plot-data/apadetect-p50/supplement/tables',
}
SCOPE={'mouse':['CIP','WOT'],'human':['GSE118258']}
INPUTS={}

def register(path):
    p=Path(path).resolve();h=hashlib.sha256(p.read_bytes()).hexdigest()
    if str(p) in INPUTS:assert INPUTS[str(p)]==h
    INPUTS[str(p)]=h;return p

def read(path,**kwargs):
    kwargs.setdefault('low_memory',False)
    return pd.read_csv(register(path),sep='\t',float_precision='round_trip',**kwargs)

def save(frame,name):
    if not isinstance(frame,pd.DataFrame):frame=pd.DataFrame(frame)
    path=HERE/name;path.parent.mkdir(exist_ok=True,parents=True)
    frame.to_csv(path,sep='\t',index=False,na_rep='NA',float_format='%.12g',
        compression={'method':'gzip','mtime':0} if path.suffix=='.gz' else None)

def rng_for(*keys):
    return np.random.default_rng(20261008+int(hashlib.sha256('|'.join(map(str,keys)).encode()).hexdigest()[:12],16))

def inference(values,*keys):
    v=np.asarray(values,float);v=v[np.isfinite(v)];n=len(v)
    out=dict(n_libraries=n,mean=float(v.mean()) if n else np.nan,ci_low=np.nan,ci_high=np.nan,p=np.nan,
        status='estimated' if n>=2 else 'fewer_than_two_libraries',null_method='exact_library_sign_flip' if n<=20 else 'monte_carlo_library_sign_flip')
    if n<2:return out
    r=rng_for(*keys);boot=v[r.integers(0,n,(10000,n))].mean(1)
    out['ci_low'],out['ci_high']=np.quantile(boot,[.025,.975])
    if n<=20:
        draws=np.array(list(itertools.product([-1,1],repeat=n)),dtype=float)@v/n
        out['p']=float((np.abs(draws)>=abs(v.mean())-1e-12).mean())
    else:
        draws=r.choice([-1,1],(10000,n))@v/n
        out['p']=float((1+(np.abs(draws)>=abs(v.mean())-1e-12).sum())/10001)
    return out

def cycle_archive():
    path=register(ROOT/'notes/2026-09-26-state-comparability-audit/expression_inputs.npz')
    assert INPUTS[str(path)]=='dd47800f7b3867d1ddbff5cb012a25f790d34fff2b26fddd0d3f0b4679bc8d4b'
    result={}
    with np.load(path,allow_pickle=False) as a:
        for sid in ['GSE103221','GSE115943','GSE118258','GSE178325']:
            m=pd.DataFrame(json.loads(str(a[sid+'_meta_json'])))
            species=m.species.iloc[0];symbols=a[species+'_marker_symbols'].tolist()
            programs=json.loads(str(a[species+'_program_json']));scores=a[sid+'_marker_logcp10k']
            for name in ['Cycle_S','Cycle_G2M']:
                m[name]=scores[:,[symbols.index(g.upper()) for g in programs[name].split()]].mean(1)
            for _,q in m.groupby('route'):
                for name in ['Cycle_S','Cycle_G2M']:
                    edges=np.unique(np.quantile(q[name],[1/3,2/3]));m.loc[q.index,name+'_bin']=np.searchsorted(edges,q[name],side='right')
            m['cycle_bin']=(m.Cycle_S_bin.astype(int)*3+m.Cycle_G2M_bin.astype(int)).astype(int)
            result[sid]=m
    additional=HERE/'additional_cycle_scores.npz'
    if additional.exists() and additional.with_suffix('.complete.json').exists():
        receipt=json.loads(register(additional.with_suffix('.complete.json')).read_text())
        assert receipt['status']=='complete' and receipt['sha256']==hashlib.sha256(additional.read_bytes()).hexdigest()
        with np.load(register(additional),allow_pickle=False) as a:
            for sid in ['GSE150311','GSE221739','GSE171127']:
                m=pd.DataFrame(json.loads(str(a[sid+'_meta_json'])))
                for score in ['Cycle_S','Cycle_G2M']:m[score]=a[sid+'_'+score]
                result[sid]=m
    return result

def cycle_scores(archive):
    means=[];diffs=[]
    for sid,m in archive.items():
        for (route,lib),q in m[m.purity.ge(.8)].groupby(['route','library_key']):
            by={s:q[q.node.eq(s)] for s in 'PMS'}
            for state,f in by.items():
                for score in ['Cycle_S','Cycle_G2M']:
                    means.append(dict(species=m.species.iloc[0],study_id=sid,route=route,library_key=lib,state=state,score=score,
                        n_metacells=len(f),n_cells=int(f.n_cells.sum()),mean=float(f[score].mean()),
                        eligible=len(f)>=3 and f.n_cells.sum()>=30))
            for contrast,a,b in [('P_M','P','M'),('P_S','P','S'),('M_S','M','S')]:
                for score in ['Cycle_S','Cycle_G2M']:
                    aa,bb=by[a],by[b];good=min(len(aa),len(bb))>=3 and min(aa.n_cells.sum(),bb.n_cells.sum())>=30
                    diffs.append(dict(species=m.species.iloc[0],study_id=sid,route=route,library_key=lib,contrast=contrast,score=score,
                        target_mean=float(aa[score].mean()),reference_mean=float(bb[score].mean()),
                        delta=float(aa[score].mean()-bb[score].mean()) if good else np.nan,
                        n_target_metacells=len(aa),n_reference_metacells=len(bb),
                        n_target_cells=int(aa.n_cells.sum()),n_reference_cells=int(bb.n_cells.sum()),
                        status='descriptive_library_unit' if good else 'missing_or_insufficient_state'))
    means=pd.DataFrame(means);diffs=pd.DataFrame(diffs);save(means,'cycle_score_library_states.tsv');save(diffs,'cycle_score_library_differences.tsv')
    tests=[]
    for key,q in diffs.groupby(['species','study_id','route','contrast','score']):
        tests.append(dict(zip(['species','study_id','route','contrast','score'],key))|
            dict(family='A2-cycle-scores',test_id='cycle_score|'+'|'.join(key),**inference(q.delta,*key)))
    save(tests,'cycle_score_tests.tsv')
    coverage={sid:dict(species=m.species.iloc[0],n_metacells=len(m),n_libraries=int(m.library_key.nunique()),
        source='frozen_expression_archive' if sid in ['GSE103221','GSE115943','GSE118258','GSE178325'] else 'additional_frozen_gex_extraction') for sid,m in archive.items()}
    (HERE/'score_coverage.json').write_text(json.dumps(dict(studies=coverage,all_seven_studies_covered=len(coverage)==7,
        excluded_studies=sorted(set(['GSE103221','GSE115943','GSE118258','GSE178325','GSE150311','GSE221739','GSE171127'])-set(coverage))),indent=2)+'\n')

class Study:
    def __init__(self,caller,sid,archive,pack,axes):
        key=caller+'__'+sid+'__';self.sid=sid;self.caller=caller
        meta=pd.DataFrame(json.loads(str(pack[key+'meta_json'])))
        cov=archive[sid].set_index('metacell_id')
        assert meta.metacell_id.is_unique and meta.metacell_id.isin(cov.index).all()
        for col in ['route','node','library_key','purity','n_cells','matrix_row']:
            assert np.array_equal(meta[col],cov.loc[meta.metacell_id,col]),(caller,sid,col)
        self.meta=meta.merge(cov[['Cycle_S','Cycle_G2M','cycle_bin']],left_on='metacell_id',right_index=True,validate='one_to_one')
        self.counts=sparse.csr_matrix((pack[key+'data'],pack[key+'indices'],pack[key+'indptr']),shape=tuple(pack[key+'shape'])).tocsc()
        self.axes=axes[axes.study_id.eq(sid)].reset_index(drop=True)
        assert np.array_equal(self.axes.feature_index_1based-1,pack[key+'columns'])
        self.groups={key:q.index.to_numpy() for key,q in self.meta[self.meta.purity.ge(.8)].groupby(['route','library_key','node'])}
        self.gene_indices={g:q.index.to_numpy() for g,q in self.axes.groupby('primary_gene_id')}
        self.cache={}

    def gene(self,gene):
        if gene not in self.cache:
            idx=self.gene_indices[gene];a=self.axes.iloc[idx]
            self.cache={gene:(self.counts[:,idx].toarray(),a)}
        return self.cache[gene]

    def state(self,data,route,library,state):
        rows=self.groups.get((route,library,state),np.array([],dtype=int));rows=rows[data[rows].sum(1)>=20]
        c=data[rows];total=c.sum(1);p=c/total[:,None] if len(c) else np.zeros(c.shape)
        return dict(rows=rows,c=c,p=p,cells=self.meta.n_cells.iloc[rows].to_numpy(),bins=self.meta.cycle_bin.iloc[rows].to_numpy(),
            mean=p.mean(0) if len(c) else np.zeros(data.shape[1]),
            Cycle_S=self.meta.Cycle_S.iloc[rows].to_numpy(),Cycle_G2M=self.meta.Cycle_G2M.iloc[rows].to_numpy())

def matched(states):
    common=sorted(set.intersection(*[set(s['bins']) for s in states.values()]))
    if not common:return None,'no_common_cycle_bin',{}
    kept={name:np.isin(s['bins'],common) for name,s in states.items()};audit={'common_bins':len(common)}
    for name,s in states.items():
        mask=kept[name];audit[name+'_n_mc']=int(mask.sum());audit[name+'_n_cells']=int(s['cells'][mask].sum())
    if min(audit[n+'_n_mc'] for n in states)<3:return None,'insufficient_common_metacells',audit
    if min(audit[n+'_n_cells'] for n in states)<30:return None,'insufficient_common_cells',audit
    for a,b in itertools.combinations(states,2):
        ca,cb=states[a]['c'][kept[a]],states[b]['c'][kept[b]]
        good=((ca.sum(0)+cb.sum(0)>=10)&((ca>0).sum(0)+(cb>0).sum(0)>=2)).sum()
        if good<2:return None,'fewer_than_two_joint_PAS_'+a+b,audit
    w=np.array([min(int((s['bins']==z).sum()) for s in states.values()) for z in common],float);w/=w.sum()
    usage={name:w@np.array([s['p'][s['bins']==z].mean(0) for z in common]) for name,s in states.items()}
    for name,s in states.items():
        for score in ['Cycle_S','Cycle_G2M']:
            audit['raw_'+name+'_'+score]=float(s[score].mean())
            audit['matched_'+name+'_'+score]=float(w@np.array([s[score][s['bins']==z].mean() for z in common]))
    assert all(np.isclose(v.sum(),1,rtol=0,atol=1e-12) for v in usage.values())
    audit['bin_weights']=json.dumps(dict(zip(map(str,common),w.tolist())))
    return usage,'measurable',audit

def candidates(base):
    p=read(base/'supplement_statistics/dual-reference-all-measurable.tsv.gz')
    assert not p.duplicated(['species','route','gene_id']).any()
    return p

def route_call(ps,pm,sign):
    return 'aligned' if min(ps*sign,pm*sign)>=.02-1e-12 else 'flipped' if min(-ps*sign,-pm*sign)>=.02-1e-12 else 'unaligned'

def do_dual(caller,base,studies,own,primary):
    triplets=read(base/'matched_pms/weighted_triplets.tsv.gz');triplets=triplets[triplets.setting.eq('supported50_p08')]
    wanted=set(own.gene_id)|set(primary.gene_id);triplets=triplets[triplets.gene_id.isin(wanted)&triplets.route.isin(sum(SCOPE.values(),[]))]
    result=[];max_error=0.
    for (sid,gene),q in triplets.groupby(['study_id','gene_id'],sort=False):
        study=studies[sid];data,axes=study.gene(gene);linear=axes.set_index('atlas_pas_id').WPP
        for r in q.itertuples():
            ids=json.loads(r.PAS_ids);cols=pd.Index(axes.atlas_pas_id).get_indexer(ids);assert (cols>=0).all()
            d=data[:,cols];weight=linear.loc[ids].to_numpy()
            states={name:study.state(d,r.route,lib,name) for name,lib in [('P',r.library_key_a),('M',r.library_key_a),('S',r.library_key_b)]}
            for name,col in [('P','target_usage'),('M','M_usage'),('S','reference_usage')]:
                err=float(np.max(np.abs(states[name]['mean']-json.loads(getattr(r,col)))));max_error=max(max_error,err);assert err<1e-10,(caller,sid,gene,name,err)
            raw={name:float(s['mean']@weight) for name,s in states.items()};u,status,audit=matched(states)
            row=dict(caller=caller,species=r.species,study_id=sid,route=r.route,gene_id=gene,gene_symbol=r.gene_symbol,
                library_key=r.library_key_a,source_library_key=r.library_key_b,triplet_weight=r.triplet_weight,status=status,
                raw_PS=raw['P']-raw['S'],raw_PM=raw['P']-raw['M'],raw_MS=raw['M']-raw['S'],**audit)
            for co,a,b in [('PS','P','S'),('PM','P','M'),('MS','M','S')]:row['matched_'+co]=float((u[a]-u[b])@weight) if u else np.nan
            if u:assert abs(row['matched_PS']-row['matched_PM']-row['matched_MS'])<1e-12
            result.append(row)
    t=pd.DataFrame(result);save(t,caller+'/dual_triplets.tsv.gz');routes=[];summaries=[];library_tests=[];global_rows=[]
    balance=[]
    for (sp,route),q in t[t.status.eq('measurable')].groupby(['species','route']):
        for score in ['Cycle_S','Cycle_G2M']:
            balance.append(dict(species=sp,route=route,score=score,n_triplets=len(q),
                raw_PM_mean=float((q['raw_P_'+score]-q['raw_M_'+score]).mean()),
                matched_PM_mean=float((q['matched_P_'+score]-q['matched_M_'+score]).mean()),
                raw_PS_mean=float((q['raw_P_'+score]-q['raw_S_'+score]).mean()),
                matched_PS_mean=float((q['matched_P_'+score]-q['matched_S_'+score]).mean())))
    save(balance,caller+'/dual_cycle_balance_descriptive.tsv')
    cohorts={'caller_own':own}
    if caller!='polyapipe':cohorts['fixed_polyapipe']=primary
    for cohort,points in cohorts.items():
        records=[]
        for r in points.itertuples():
            q=t[t.species.eq(r.species)&t.route.eq(r.route)&t.gene_id.eq(r.gene_id)];ok=q[q.status.eq('measurable')]
            sign=1 if r.recurrent_direction=='positive' else -1
            row=dict(caller=caller,cohort=cohort,species=r.species,study_id=r.study_id,route=r.route,gene_id=r.gene_id,gene_symbol=r.gene_symbol,
                recurrent_direction=r.recurrent_direction,baseline_dual=bool(r.dual_qualified_gene),n_original_triplets=len(q),n_matched_triplets=len(ok),
                original_call=route_call(r.DWUI_PS,r.DWUI_PM,sign),matched_call='unmeasurable',raw_all_PS=np.nan,raw_all_PM=np.nan,
                raw_common_PS=np.nan,raw_common_PM=np.nan,matched_PS=np.nan,matched_PM=np.nan)
            if len(q):
                w=q.triplet_weight.to_numpy(copy=True);w/=w.sum()
                for co in ['PS','PM']:row['raw_all_'+co]=float(w@q['raw_'+co])
                row['raw_caller_call']=route_call(row['raw_all_PS'],row['raw_all_PM'],sign)
                if cohort=='caller_own':
                    assert abs(row['raw_all_PS']-r.DWUI_PS)<1e-10 and abs(row['raw_all_PM']-r.DWUI_PM)<1e-10
            else:row['raw_caller_call']='unmeasurable'
            if len(ok):
                w=ok.triplet_weight.to_numpy(copy=True);w/=w.sum()
                for co in ['PS','PM']:
                    row['raw_common_'+co]=float(w@ok['raw_'+co]);row['matched_'+co]=float(w@ok['matched_'+co])
                row['matched_call']=route_call(row['matched_PS'],row['matched_PM'],sign)
            records.append(row)
        frame=pd.DataFrame(records);routes.append(frame)
        for (sp,route),f in frame.groupby(['species','route']):
            for label,z in [('all_candidates',f),('baseline_dual',f[f.baseline_dual]),('original_route_aligned',f[f.original_call.eq('aligned')])]:
                summaries.append(dict(caller=caller,cohort=cohort,species=sp,route=route,subset=label,candidates=len(z),
                    measurable=int(z.matched_call.ne('unmeasurable').sum()),retained=int(z.matched_call.eq('aligned').sum()),
                    flipped=int(z.matched_call.eq('flipped').sum()),lost_alignment=int(z.matched_call.eq('unaligned').sum()),
                    unmeasurable=int(z.matched_call.eq('unmeasurable').sum())))
            # Candidate-gene proportions first calculated per physical target library.
            q=t[t.species.eq(sp)&t.route.eq(route)&t.gene_id.isin(f.gene_id)&t.status.eq('measurable')]
            signs=f.set_index('gene_id').recurrent_direction.map({'positive':1,'negative':-1})
            librows=[]
            for (lib,gene),z in q.groupby(['library_key','gene_id']):
                w=z.triplet_weight.to_numpy(copy=True);w/=w.sum();sg=signs.loc[gene]
                original=route_call(float(w@z.raw_PS),float(w@z.raw_PM),sg)=='aligned'
                current=route_call(float(w@z.matched_PS),float(w@z.matched_PM),sg)=='aligned'
                librows.append(dict(library_key=lib,gene_id=gene,raw_aligned=int(original),matched_aligned=int(current)))
            if librows:
                li=pd.DataFrame(librows).groupby('library_key').agg(n_genes=('gene_id','size'),raw_fraction=('raw_aligned','mean'),matched_fraction=('matched_aligned','mean')).reset_index()
                li['delta']=li.matched_fraction-li.raw_fraction;li['caller']=caller;li['cohort']=cohort;li['species']=sp;li['route']=route
                save(li,caller+'/'+cohort+'-'+sp+'-'+route+'-library_alignment.tsv')
                library_tests.append(dict(caller=caller,cohort=cohort,species=sp,route=route,family='A2-cycle-alignment',
                    test_id='alignment_change|'+caller+'|'+cohort+'|'+sp+'|'+route,**inference(li.delta,caller,cohort,sp,route)))
            else:library_tests.append(dict(caller=caller,cohort=cohort,species=sp,route=route,family='A2-cycle-alignment',
                test_id='alignment_change|'+caller+'|'+cohort+'|'+sp+'|'+route,**inference([],caller,cohort,sp,route)))
        for (sp,gene),z in frame.groupby(['species','gene_id']):
            assert set(z.route)==set(SCOPE[sp])
            status='unmeasurable' if z.matched_call.eq('unmeasurable').any() else 'aligned' if z.matched_call.eq('aligned').all() else 'flipped' if z.matched_call.eq('flipped').all() else 'unaligned'
            global_rows.append(dict(caller=caller,cohort=cohort,species=sp,gene_id=gene,gene_symbol=z.gene_symbol.iloc[0],baseline_dual=bool(z.baseline_dual.all()),
                matched_call=status,routes=';'.join(z.route),measurable_routes=int(z.matched_call.ne('unmeasurable').sum())))
    save(pd.concat(routes),caller+'/dual_route_calls.tsv');save(summaries,caller+'/dual_route_summary.tsv');save(library_tests,caller+'/dual_alignment_tests.tsv')
    global_frame=pd.DataFrame(global_rows);save(global_frame,caller+'/dual_species_calls.tsv');species_summary=[]
    for key,z in global_frame.groupby(['caller','cohort','species']):
        for subset,f in [('all_candidates',z),('baseline_dual',z[z.baseline_dual])]:
            species_summary.append(dict(zip(['caller','cohort','species'],key))|dict(subset=subset,candidates=len(f),
                measurable=int(f.matched_call.ne('unmeasurable').sum()),retained=int(f.matched_call.eq('aligned').sum()),
                flipped=int(f.matched_call.eq('flipped').sum()),lost_alignment=int(f.matched_call.eq('unaligned').sum()),
                unmeasurable=int(f.matched_call.eq('unmeasurable').sum())))
    save(species_summary,caller+'/dual_species_summary.tsv');return max_error

def do_route_means(caller,base,studies):
    pairs=read(base/'state_vectors/pair_effects.primary.tsv.gz');pairs=pairs[pairs.contrast.eq('P_M')&pairs.design.eq('same_library')&pairs.study_id.isin(studies)]
    rows=[];max_error=0.
    for (sid,gene),q in pairs.groupby(['study_id','gene_id'],sort=False):
        study=studies[sid];data,axes=study.gene(gene)
        for r in q.itertuples():
            ids=json.loads(r.PAS_ids);ix=pd.Index(axes.atlas_pas_id).get_indexer(ids);assert (ix>=0).all()
            d=data[:,ix];linear=axes.WPP.to_numpy()[ix]
            states={state:study.state(d,r.route,r.library_key_a,state) for state in 'PM'}
            raw=float((states['P']['mean']-states['M']['mean'])@linear);err=abs(raw-r.delta_DWUI);max_error=max(max_error,err);assert err<1e-10
            u,status,audit=matched(states)
            rows.append(dict(caller=caller,species=r.species,study_id=sid,route=r.route,library_key=r.library_key_a,gene_id=gene,gene_symbol=r.gene_symbol,
                raw_PM=raw,matched_PM=float((u['P']-u['M'])@linear) if u else np.nan,status=status,**audit))
    frame=pd.DataFrame(rows);assert not frame.duplicated(['study_id','route','library_key','gene_id']).any();save(frame,caller+'/pm_gene_library_effects.tsv.gz')
    libs=[];tests=[]
    for (sp,sid,route,lib),q in frame.groupby(['species','study_id','route','library_key']):
        common=q[q.status.eq('measurable')]
        for version,part,col in [('raw_all',q,'raw_PM'),('raw_common',common,'raw_PM'),('cycle_matched',common,'matched_PM')]:
            libs.append(dict(caller=caller,species=sp,study_id=sid,route=route,library_key=lib,version=version,n_genes=len(part),mean=float(part[col].mean())))
    libs=pd.DataFrame(libs);save(libs,caller+'/pm_library_means.tsv')
    for key,q in libs.groupby(['species','study_id','route','version']):
        tests.append(dict(caller=caller,**dict(zip(['species','study_id','route','version'],key)),family='A8-route-means',
            test_id='pm_mean|'+caller+'|'+'|'.join(key),**inference(q['mean'],caller,*key)))
    save(tests,caller+'/pm_route_mean_tests.tsv');return max_error

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--caller',choices=list(BASES));parser.add_argument('--skip-means',action='store_true');parser.add_argument('--scores-only',action='store_true');args=parser.parse_args()
    archive=cycle_archive();cycle_scores(archive)
    register(HERE/'run_cycle.py');register(HERE/'PROTOCOL.md')
    if args.scores_only:
        path=HERE/'input_sha256.json';previous=json.loads(path.read_text()) if path.exists() else {}
        previous.update(INPUTS);path.write_text(json.dumps(previous,indent=2)+'\n')
        audit_path=HERE/'audit.json'
        if audit_path.exists():
            audit=json.loads(audit_path.read_text());audit['score_studies']=sorted(archive)
            audit['all_seven_cycle_score_studies_covered']=len(archive)==7
            audit['additional_cycle_score_sha256']=INPUTS.get(str((HERE/'additional_cycle_scores.npz').resolve()))
            audit_path.write_text(json.dumps(audit,indent=2)+'\n')
        print('Cycle scores complete for',sorted(archive),flush=True);return
    primary=candidates(BASES['polyapipe'])
    assert primary.groupby('species').gene_id.nunique().to_dict()=={'human':250,'mouse':395}
    assert primary[primary.dual_qualified_gene].groupby('species').gene_id.nunique().to_dict()=={'human':128,'mouse':222}
    with np.load(register(HERE/'pas_subsets.npz'),allow_pickle=False) as pack:
        extraction=json.loads(register(HERE/'extraction_manifest.json').read_text())
        assert json.loads(str(pack['provenance_json']))=={r[k+'_path']:r[k+'_sha256'] for r in extraction['remote'] for k in ['counts','meta','full_meta']}
        audits=[]
        for caller in [args.caller] if args.caller else BASES:
            base=BASES[caller];own=candidates(base);axes=read(base/'state_vectors/PAS_axes.primary.tsv.gz')
            studies={sid:Study(caller,sid,archive,pack,axes) for sid in ['GSE103221','GSE115943','GSE118258']}
            dual_error=do_dual(caller,base,studies,own,primary);print(caller,'dual complete',dual_error,flush=True)
            pm_error=do_route_means(caller,base,studies) if not args.skip_means else None
            audits.append(dict(caller=caller,max_raw_dual_usage_error=dual_error,max_raw_pm_effect_error=pm_error,
                primary_counts_verified=True,shared_bin_additivity_verified=True,library_inference_only=True))
    for path,value in INPUTS.items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==value
    (HERE/'input_sha256.json').write_text(json.dumps(INPUTS,indent=2)+'\n')
    (HERE/'audit.json').write_text(json.dumps(dict(checks=audits,seed=20261008,bootstrap=10000,protocol_sha256=hashlib.sha256((HERE/'PROTOCOL.md').read_bytes()).hexdigest()),indent=2)+'\n')
    print('complete',audits,flush=True)

if __name__=='__main__':main()
