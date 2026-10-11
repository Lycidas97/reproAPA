"""Frozen-context continuous cycle balance and descriptive fixed 2x2 sensitivity."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.linalg import qr
from scipy.optimize import linprog, minimize

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
OLD=ROOT/'reviews/2026-10-08-plan-a-analysis/cycle'
BASES={
    'polyapipe':ROOT/'notes/2026-10-06-supplement-polyapipe/contract-r2/bundle/tables',
    'apadetect_p50':ROOT/'notes/2026-10-06-dual-main-plot-data/apadetect-p50/supplement/tables',
}
STUDIES=['GSE103221','GSE115943','GSE118258','GSE150311','GSE221739']
ROUTES=[('mouse','GSE103221','CIP'),('mouse','GSE115943','WOT'),('human','GSE118258','GSE118258'),
    ('human','GSE150311','GSE150311'),('human','GSE221739','GSE221739')]
SCORES=['Cycle_S','Cycle_G2M']
METHODS=['continuous','fixed_2x2']
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
    if not len(frame.columns):
        identity=['caller','species','study_id','route','library_key','source_library_key','gene_id','method','scope','status']
        state_cols=identity+['state','n_raw_mc','n_mc','n_cells','ESS','ESS_fraction','max_weight']+[s+'_'+z for s in SCORES for z in ['raw_mean','raw_variance','mean','variance','q10','q50','q90']]
        schemas={
            'pm_candidate_library_changes.tsv':['caller','species','study_id','route','library_key','method','cohort','n_genes','mean_signed_change','mean_support_trim','mean_reweight','raw_context_common_signed','adjusted_signed'],
            'dual_candidate_gene_library_margins.tsv':['caller','species','study_id','route','library_key','gene_id','method','cohort','raw_context_common_margin','raw_support_common_margin','adjusted_margin','margin_change','support_trim','reweight'],
            'pm_gene_state_balance.tsv.gz':state_cols,'dual_gene_state_balance.tsv.gz':state_cols,
        }
        for scope,prefixes in [('pm',['PM']),('dual',['PM','PS','MS'])]:
            schemas[scope+'_gene_pair_balance.tsv.gz']=identity+[a+'_'+s+'_'+z for a in prefixes for s in SCORES for z in ['SMD','variance_ratio','KS','raw_pooled_sd','library_scale_SMD']]+['PM_joint_PAS']+(['PS_joint_PAS'] if scope=='dual' else [])
        basename=Path(name).name;assert basename in schemas,('Empty output needs schema',name)
        frame=pd.DataFrame(columns=schemas[basename])
    path=HERE/name;path.parent.mkdir(exist_ok=True,parents=True)
    frame.to_csv(path,sep='\t',index=False,na_rep='NA',float_format='%.12g',
        compression={'method':'gzip','mtime':0} if path.suffix=='.gz' else None)

def rng_for(*keys):
    return np.random.default_rng(20261008+int(hashlib.sha256('|'.join(map(str,keys)).encode()).hexdigest()[:12],16))

def inference(values,*keys):
    v=np.asarray(values,float);v=v[np.isfinite(v)];n=len(v)
    out=dict(n_libraries=n,mean=float(v.mean()) if n else np.nan,ci_low=np.nan,ci_high=np.nan,p=np.nan,
        status='estimated' if n>=2 else 'fewer_than_two_libraries',generic_signflip_lower_bound=2/(2**n) if n>=2 else np.nan,min_attainable_p=np.nan,
        null_method='exact_library_sign_flip' if n<=20 else 'monte_carlo_library_sign_flip')
    if n<2:return out
    r=rng_for(*keys);boot=v[r.integers(0,n,(10000,n))].mean(1)
    out['ci_low'],out['ci_high']=np.quantile(boot,[.025,.975])
    if n<=20:
        draws=np.array(list(itertools.product([-1,1],repeat=n)),dtype=float)@v/n
        out['p']=float((np.abs(draws)>=abs(v.mean())-1e-12).mean())
        out['min_attainable_p']=float((np.abs(draws)>=np.max(np.abs(draws))-1e-12).mean())
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
            sp=m.species.iloc[0];symbols=a[sp+'_marker_symbols'].tolist();programs=json.loads(str(a[sp+'_program_json']))
            for name in SCORES:m[name]=a[sid+'_marker_logcp10k'][:,[symbols.index(g.upper()) for g in programs[name].split()]].mean(1)
            result[sid]=m
    path=register(OLD/'additional_cycle_scores.npz');receipt=json.loads(register(OLD/'additional_cycle_scores.complete.json').read_text())
    assert receipt['status']=='complete' and receipt['sha256']==INPUTS[str(path)]
    with np.load(path,allow_pickle=False) as a:
        for sid in ['GSE150311','GSE221739','GSE171127']:
            m=pd.DataFrame(json.loads(str(a[sid+'_meta_json'])))
            for score in SCORES:m[score]=a[sid+'_'+score]
            result[sid]=m
    for m in result.values():
        m['cycle_bin']=0
        for _,q in m.groupby('route'):
            eligible=q[q.purity.ge(.8)]
            for k,score in enumerate(SCORES):
                edge=float(eligible[score].median());m.loc[q.index,'cycle_bin']+=np.where(q[score]>=edge,2 if k==0 else 1,0)
        m['cycle_bin']=m.cycle_bin.astype(int)
    return result

def weighted_quantiles(x,w):
    ix=np.argsort(x,kind='stable');xx=x[ix];cdf=np.cumsum(w[ix])
    return xx[np.minimum(np.searchsorted(cdf,[.1,.5,.9],side='left'),len(xx)-1)]

def weighted_ks(x,wx,y,wy):
    grid=np.unique(np.r_[x,y]);ix=np.argsort(x,kind='stable');iy=np.argsort(y,kind='stable')
    cx=np.r_[0,np.cumsum(wx[ix])];cy=np.r_[0,np.cumsum(wy[iy])]
    return float(np.max(np.abs(cx[np.searchsorted(x[ix],grid,side='right')]-cy[np.searchsorted(y[iy],grid,side='right')])))

def balance(raw,weights,library_sd=None):
    """All frozen gates; raw groups define SMD scale, weights define estimand."""
    stats={};fail=[]
    for name,s in raw.items():
        w=np.asarray(weights[name],float).copy();w[w<=1e-8]=0
        if w.sum()>0:w/=w.sum()
        mask=w>0;nn=int(mask.sum());cells=int(s['cells'][mask].sum());ess=float(1/(w@w)) if nn else 0.
        rec=dict(n_raw_mc=len(w),n_mc=nn,n_cells=cells,ESS=ess,ESS_fraction=ess/len(w) if len(w) else 0.,max_weight=float(w.max()) if len(w) else np.nan)
        if nn<3:fail.append('insufficient_metacells_'+name)
        if cells<30:fail.append('insufficient_cells_'+name)
        if ess<3-1e-10:fail.append('insufficient_ESS_'+name)
        for score in SCORES:
            x=s[score];mu=float(w@x) if nn else np.nan;va=float(w@((x-mu)**2)) if nn else np.nan
            qq=weighted_quantiles(x[mask],w[mask]) if nn else [np.nan]*3
            rec[score+'_raw_mean']=float(x.mean()) if len(x) else np.nan
            rec[score+'_raw_variance']=float(x.var()) if len(x) else np.nan
            rec[score+'_mean']=mu;rec[score+'_variance']=va
            for label,value in zip(['q10','q50','q90'],qq):rec[score+'_'+label]=float(value)
        stats[name]=rec
    pair_stats={}
    for a,b in itertools.combinations(raw,2):
        wa=np.asarray(weights[a],float).copy();wb=np.asarray(weights[b],float).copy()
        wa[wa<=1e-8]=0;wb[wb<=1e-8]=0
        if wa.sum()>0:wa/=wa.sum()
        if wb.sum()>0:wb/=wb.sum()
        for score in SCORES:
            x=raw[a][score];y=raw[b][score];mx=stats[a][score+'_mean'];my=stats[b][score+'_mean']
            vx=stats[a][score+'_variance'];vy=stats[b][score+'_variance']
            sd=float(np.sqrt((x.var()+y.var())/2)) if len(x) and len(y) else np.nan
            smd=abs(mx-my)/sd if sd>0 else 0. if abs(mx-my)<=1e-12 else np.inf
            ratio=max(vx,vy)/min(vx,vy) if min(vx,vy)>1e-14 else 1. if max(vx,vy)<=1e-14 and abs(mx-my)<=1e-12 else np.inf
            ks=weighted_ks(x,wa,y,wb) if wa.sum()>0 and wb.sum()>0 else np.nan
            key=a+b+'_'+score
            pair_stats[key+'_SMD']=float(smd);pair_stats[key+'_variance_ratio']=float(ratio);pair_stats[key+'_KS']=float(ks);pair_stats[key+'_raw_pooled_sd']=sd
            if library_sd is not None:
                z=library_sd[key];pair_stats[key+'_library_scale_SMD']=abs(mx-my)/z if z>0 else 0. if abs(mx-my)<=1e-12 else np.inf
            if not np.isfinite(smd) or smd>.1+1e-10:fail.append('SMD_'+key)
            if not np.isfinite(ratio) or ratio>2+1e-10:fail.append('variance_ratio_'+key)
            if not np.isfinite(ks) or ks>.2+1e-10:fail.append('KS_'+key)
    return not fail,stats,pair_stats,sorted(set(fail))

def optimize(raw):
    names=list(raw);sizes=[len(raw[n]['rows']) for n in names];offset=np.cumsum([0]+sizes);n=offset[-1]
    if min(sizes,default=0)<3 or min((s['cells'].sum() for s in raw.values()),default=0)<30:
        return None,dict(solver_status='insufficient_metadata_state',constraint_error=np.nan,objective=np.nan)
    A=[];b=[];scale={}
    for j,name in enumerate(names):
        row=np.zeros(n);row[offset[j]:offset[j+1]]=1;A.append(row);b.append(1.)
    for score in SCORES:
        all_x=np.concatenate([raw[name][score] for name in names]);center=float(all_x.mean());sd=float(all_x.std());scale[score]=(center,sd)
        if sd==0:continue
        for j in range(1,len(names)):
            row=np.zeros(n);row[:offset[1]]=(raw[names[0]][score]-center)/sd
            row[offset[j]:offset[j+1]]=-(raw[names[j]][score]-center)/sd;A.append(row);b.append(0.)
    fullA=np.array(A);fullb=np.array(b)
    _,R,piv=qr(fullA.T,mode='economic',pivoting=True);rank=int((np.abs(np.diag(R))>1e-10).sum());keep=np.sort(piv[:rank]);A=fullA[keep];b=fullb[keep]
    equality_point=np.linalg.lstsq(A,b,rcond=1e-10)[0]
    if np.max(np.abs(fullA@equality_point-fullb))>1e-8:
        return None,dict(solver_status='no_common_mean_support',constraint_error=np.nan,objective=np.nan,message='inconsistent_full_equalities_after_rank_reduction')
    feasible=linprog(np.zeros(n),A_eq=A,b_eq=b,bounds=(0,None),method='highs')
    if not feasible.success:return None,dict(solver_status='no_common_mean_support' if feasible.status==2 else 'linprog_failure',constraint_error=np.nan,objective=np.nan,message=feasible.message)
    mult=np.repeat(sizes,sizes);uniform=np.concatenate([np.repeat(1/z,z) for z in sizes])
    fun=lambda w:float(np.sum(mult*(w-uniform)**2))
    jac=lambda w:2*mult*(w-uniform)
    sol=minimize(fun,feasible.x,jac=jac,method='SLSQP',bounds=[(0,None)]*n,
        constraints={'type':'eq','fun':lambda w:A@w-b,'jac':lambda w:A},options={'ftol':1e-10,'maxiter':300})
    error=float(np.max(np.abs(fullA@sol.x-fullb)))
    info=dict(solver_status='optimal' if sol.success and error<=1e-8 else 'SLSQP_failure',constraint_error=error,
        objective=fun(sol.x),iterations=int(sol.nit),message=sol.message,n_constraints=len(fullA),rank=rank)
    if info['solver_status']!='optimal':return None,info
    result={name:sol.x[offset[j]:offset[j+1]].copy() for j,name in enumerate(names)}
    for name,w in result.items():
        w[w<=1e-8]=0;assert w.sum()>0;w/=w.sum()
    mean_error=max(abs(result[a]@raw[a][s]-result[b]@raw[b][s]) for a,b in itertools.combinations(names,2) for s in SCORES)
    info['post_zero_mean_error']=float(mean_error)
    if mean_error>1e-8:info['solver_status']='zero_trim_constraint_failure';return None,info
    return result,info

def bins_weight(raw):
    common=sorted(set.intersection(*[set(s['bins']) for s in raw.values()]))
    if not common:return None,dict(solver_status='no_common_2x2_bin',constraint_error=np.nan,objective=np.nan)
    target=np.array([min(int((s['bins']==z).sum()) for s in raw.values()) for z in common],float);target/=target.sum()
    result={name:np.zeros(len(s['rows'])) for name,s in raw.items()}
    for name,s in raw.items():
        for z,t in zip(common,target):result[name][s['bins']==z]=t/(s['bins']==z).sum()
    return result,dict(solver_status='fixed_bins',constraint_error=np.nan,objective=np.nan,common_bins=json.dumps(list(map(int,common))),bin_weights=json.dumps(target.tolist()))

class Study:
    def __init__(self,caller,sid,archive,pack,axes):
        key=caller+'__'+sid+'__';self.sid=sid;self.caller=caller
        meta=pd.DataFrame(json.loads(str(pack[key+'meta_json'])));cov=archive[sid].set_index('metacell_id')
        assert meta.metacell_id.is_unique and meta.metacell_id.isin(cov.index).all()
        for col in ['route','node','library_key','purity','n_cells','matrix_row']:assert np.array_equal(meta[col],cov.loc[meta.metacell_id,col]),(caller,sid,col)
        self.meta=meta.merge(cov[SCORES+['cycle_bin']],left_on='metacell_id',right_index=True,validate='one_to_one')
        self.counts=sparse.csr_matrix((pack[key+'data'],pack[key+'indices'],pack[key+'indptr']),shape=tuple(pack[key+'shape'])).tocsc()
        self.axes=axes[axes.study_id.eq(sid)].reset_index(drop=True);assert np.array_equal(self.axes.feature_index_1based-1,pack[key+'columns'])
        self.groups={key:q.index.to_numpy() for key,q in self.meta[self.meta.purity.ge(.8)].groupby(['route','library_key','node'])}
        self.gene_indices={g:q.index.to_numpy() for g,q in self.axes.groupby('primary_gene_id')};self.cache={};self.contexts={}
        self.weight_rows=[];self.context_rows=[];self.context_state_rows=[];self.context_pair_rows=[]

    def gene(self,gene):
        if gene not in self.cache:
            ix=self.gene_indices[gene];self.cache={gene:(self.counts[:,ix].toarray(),self.axes.iloc[ix])}
        return self.cache[gene]

    def state(self,route,library,name,data=None):
        rows=self.groups.get((route,library,name),np.array([],dtype=int))
        if data is not None:rows=rows[data[rows].sum(1)>=20]
        q=self.meta.iloc[rows];result=dict(rows=rows,cells=q.n_cells.to_numpy(),bins=q.cycle_bin.to_numpy())
        for score in SCORES:result[score]=q[score].to_numpy()
        if data is not None:
            c=data[rows];total=c.sum(1);result['c']=c;result['p']=c/total[:,None] if len(c) else np.zeros(c.shape)
            result['mean']=result['p'].mean(0) if len(c) else np.zeros(data.shape[1])
        return result

    def context(self,route,library,source,scope,method):
        key=(route,library,source,scope,method)
        if key not in self.contexts:
            raw={name:self.state(route,lib,name) for name,lib in [('P',library),('M',library)]+([('S',source)] if scope=='PMS' else [])}
            weights,info=optimize(raw) if method=='continuous' else bins_weight(raw)
            ident=dict(caller=self.caller,study_id=self.sid,route=route,library_key=library,source_library_key=source,scope=scope,method=method,
                context_id='|'.join(map(str,[self.caller,self.sid,*key])))
            sd={a+b+'_'+s:float(np.sqrt((raw[a][s].var()+raw[b][s].var())/2)) if len(raw[a][s]) and len(raw[b][s]) else np.nan for a,b in itertools.combinations(raw,2) for s in SCORES}
            if weights is None:passed=False;ss={};pp={};fail=[info['solver_status']]
            else:passed,ss,pp,fail=balance(raw,weights)
            full={name:np.zeros(len(self.meta)) for name in raw}
            if weights is not None:
                for name,w in weights.items():
                    full[name][raw[name]['rows']]=w
                    for ix,value in zip(raw[name]['rows'],w):self.weight_rows.append(ident|dict(state=name,metacell_id=self.meta.metacell_id.iloc[ix],row=int(ix),weight=float(value),
                        Cycle_S=float(self.meta.Cycle_S.iloc[ix]),Cycle_G2M=float(self.meta.Cycle_G2M.iloc[ix]),n_cells=int(self.meta.n_cells.iloc[ix])))
                for name,z in ss.items():self.context_state_rows.append(ident|dict(state=name,**z))
                self.context_pair_rows.append(ident|pp)
            row=ident|info|dict(status='balanced' if passed else 'failed',failures=json.dumps(fail),n_metadata=int(sum(len(s['rows']) for s in raw.values())))
            self.context_rows.append(row);self.contexts[key]=dict(weights=full,passed=passed,status=row['status'],failures=fail,sd=sd,ident=ident)
        return self.contexts[key]

def gene_adjust(states,context):
    if not context['passed']:return None,'context_failed',dict(failures=json.dumps(context['failures'])),[],[]
    weights={name:context['weights'][name][s['rows']].copy() for name,s in states.items()}
    for name,w in weights.items():
        w[w<=1e-8]=0
        if w.sum()>0:w/=w.sum()
    passed,ss,pp,fail=balance(states,weights,context['sd'])
    for a,b in [('P','M')]+([('P','S')] if 'S' in states else []):
        ca=states[a]['c'][weights[a]>1e-8];cb=states[b]['c'][weights[b]>1e-8]
        good=int(((ca.sum(0)+cb.sum(0)>=10)&((ca>0).sum(0)+(cb>0).sum(0)>=2)).sum())
        pp[a+b+'_joint_PAS']=good
        if good<2:fail.append('fewer_than_two_joint_PAS_'+a+b)
    fail=sorted(set(fail));status='measurable' if not fail else 'gene_failed'
    if fail:return None,status,dict(failures=json.dumps(fail)),ss,pp
    usage={name:weights[name]@s['p'] for name,s in states.items()}
    support={name:s['p'][weights[name]>1e-8].mean(0) for name,s in states.items()}
    assert all(np.isclose(v.sum(),1,rtol=0,atol=1e-12) for v in usage.values())
    return dict(adjusted=usage,raw_support_common=support),status,dict(failures='[]'),ss,pp

def recurrent(base):
    p=read(base/'recurrence/recurrence_members.tsv.gz')
    p=p[p.metric.eq('WPP')&p.contrast.eq('P_M')&p.threshold.eq(.02)&p.variant.eq('primary')&~p.series_matching]
    assert not p.duplicated(['species','gene_id']).any();return p

def dual_candidates(base):
    p=read(base/'supplement_statistics/dual-reference-all-measurable.tsv.gz')
    assert not p.duplicated(['species','route','gene_id']).any();return p

def call(ps,pm,sign):
    return 'aligned' if min(ps*sign,pm*sign)>=.02-1e-12 else 'flipped' if min(-ps*sign,-pm*sign)>=.02-1e-12 else 'unaligned'

def pm_call(pm,sign):
    return 'aligned' if pm*sign>=.02-1e-12 else 'flipped' if -pm*sign>=.02-1e-12 else 'unaligned'

def add_diagnostics(row,ss,pp,states,rows_state,rows_pair):
    ident={k:row[k] for k in ['caller','species','study_id','route','library_key','source_library_key','gene_id','method','scope','status']}
    if ss:
        for name,z in ss.items():rows_state.append(ident|dict(state=name,**z))
        rows_pair.append(ident|pp)

def do_pm(caller,base,studies,primary,own):
    pairs=read(base/'state_vectors/pair_effects.primary.tsv.gz');pairs=pairs[pairs.contrast.eq('P_M')&pairs.design.eq('same_library')&pairs.setting.eq('supported50_p08')&pairs.study_id.isin(studies)]
    assert not pairs.duplicated(['study_id','route','library_key_a','gene_id']).any()
    rows=[];state_diag=[];pair_diag=[];max_error=0.
    for (sid,gene),q in pairs.groupby(['study_id','gene_id'],sort=False):
        study=studies[sid];data,axes=study.gene(gene)
        for r in q.itertuples():
            ids=json.loads(r.PAS_ids);ix=pd.Index(axes.atlas_pas_id).get_indexer(ids);assert (ix>=0).all()
            d=data[:,ix];linear=axes.WPP.to_numpy()[ix];states={s:study.state(r.route,r.library_key_a,s,d) for s in 'PM'}
            raw=float((states['P']['mean']-states['M']['mean'])@linear);err=abs(raw-r.delta_DWUI);max_error=max(max_error,err);assert err<1e-10,(caller,sid,gene,err)
            for method in METHODS:
                context=study.context(r.route,r.library_key_a,'','PM',method);u,status,audit,ss,pp=gene_adjust(states,context)
                row=dict(caller=caller,species=r.species,study_id=sid,route=r.route,library_key=r.library_key_a,source_library_key='',gene_id=gene,gene_symbol=r.gene_symbol,
                    method=method,scope='PM',status=status,original_pair_weight=r.weight,raw_all_PM=raw,raw_context_common_PM=raw if u else np.nan,
                    raw_support_common_PM=float((u['raw_support_common']['P']-u['raw_support_common']['M'])@linear) if u else np.nan,
                    adjusted_PM=float((u['adjusted']['P']-u['adjusted']['M'])@linear) if u else np.nan,**audit)
                rows.append(row);add_diagnostics(row,ss,pp,states,state_diag,pair_diag)
    frame=pd.DataFrame(rows);save(frame,caller+'/pm_gene_library_effects.tsv.gz');save(state_diag,caller+'/pm_gene_state_balance.tsv.gz');save(pair_diag,caller+'/pm_gene_pair_balance.tsv.gz')
    libs=[];formal=[];description=[];candidate_libs=[];candidate_calls=[]
    for method in METHODS:
        fm=frame[frame.method.eq(method)]
        for sp,sid,route in ROUTES:
            q=fm[fm.study_id.eq(sid)]
            for lib,z in q.groupby('library_key'):
                good=z[z.status.eq('measurable')]
                for version,part in [('raw_all',z),('raw_context_common',good),('raw_support_common',good),('adjusted',good)]:
                    libs.append(dict(caller=caller,species=sp,study_id=sid,route=route,library_key=lib,method=method,version=version,n_genes=len(part),mean=float(part[version+'_PM'].mean())))
            for cohort,points in [('fixed_primary',primary),('caller_own',own)]:
                p=points[points.species.eq(sp)];signs=p.set_index('gene_id').direction.map({'positive':1,'negative':-1})
                z=q[q.gene_id.isin(signs.index)]
                for lib,g in z[z.status.eq('measurable')].groupby('library_key'):
                    sg=g.gene_id.map(signs).to_numpy()
                    candidate_libs.append(dict(caller=caller,species=sp,study_id=sid,route=route,library_key=lib,method=method,cohort=cohort,n_genes=len(g),
                        mean_signed_change=float(np.mean(sg*(g.adjusted_PM-g.raw_context_common_PM))),
                        mean_support_trim=float(np.mean(sg*(g.raw_support_common_PM-g.raw_context_common_PM))),
                        mean_reweight=float(np.mean(sg*(g.adjusted_PM-g.raw_support_common_PM))),
                        raw_context_common_signed=float(np.mean(sg*g.raw_context_common_PM)),adjusted_signed=float(np.mean(sg*g.adjusted_PM))))
                for r in p.itertuples():
                    zz=z[z.gene_id.eq(r.gene_id)];gg=zz[zz.status.eq('measurable')];sg=1 if r.direction=='positive' else -1
                    candidate_calls.append(dict(caller=caller,species=sp,study_id=sid,route=route,gene_id=r.gene_id,gene_symbol=r.gene_symbol,method=method,cohort=cohort,
                        n_original_libraries=len(zz),n_adjusted_libraries=len(gg),raw_all_PM=float(zz.raw_all_PM.mean()),
                        raw_context_common_PM=float(gg.raw_context_common_PM.mean()),raw_support_common_PM=float(gg.raw_support_common_PM.mean()),adjusted_PM=float(gg.adjusted_PM.mean()),
                        raw_all_call=pm_call(float(zz.raw_all_PM.mean()),sg) if len(zz) else 'unmeasurable',
                        raw_context_common_call=pm_call(float(gg.raw_context_common_PM.mean()),sg) if len(gg) else 'unmeasurable',
                        adjusted_call=pm_call(float(gg.adjusted_PM.mean()),sg) if len(gg) else 'unmeasurable'))
    libs=pd.DataFrame(libs);cl=pd.DataFrame(candidate_libs);cc=pd.DataFrame(candidate_calls)
    save(libs,caller+'/pm_library_means.tsv');save(cl,caller+'/pm_candidate_library_changes.tsv');save(cc,caller+'/pm_candidate_route_calls.tsv')
    for sp,sid,route in ROUTES:
        z=cl[(cl.study_id.eq(sid))&cl.method.eq('continuous')&cl.cohort.eq('fixed_primary')] if len(cl) else pd.DataFrame()
        formal.append(dict(test_id='PM_change|'+caller+'|'+sid+'|'+route,family='A2-cycle-alignment-redesign',pas_source=caller,scope='PM',species=sp,study_id=sid,route=route,
            method='continuous',cohort='fixed_primary',statistic='signed_PM_adjusted_minus_raw_context_common',**inference(z.mean_signed_change if len(z) else [],caller,'PM',sid,route)))
        z=libs[libs.study_id.eq(sid)&libs.method.eq('continuous')&libs.version.eq('adjusted')]
        formal.append(dict(test_id='PM_mean|'+caller+'|'+sid+'|'+route,family='A8-cycle-PM-redesign',pas_source=caller,scope='PM_all_original_genes',species=sp,study_id=sid,route=route,
            method='continuous',cohort='all_original_measurable',statistic='adjusted_PM_library_mean_vs_zero',**inference(z['mean'],caller,'A8',sid,route)))
    for key,z in libs.groupby(['species','study_id','route','method','version']):
        result=inference(z['mean'],caller,'descriptive',*key);result.pop('p');description.append(dict(caller=caller,**dict(zip(['species','study_id','route','method','version'],key)),**result))
    save(description,caller+'/pm_route_means_descriptive.tsv')
    summaries=[]
    for key,z in cc.groupby(['species','study_id','route','method','cohort']):
        for subset,f in [('all_candidates',z),('raw_all_aligned',z[z.raw_all_call.eq('aligned')])]:
            summaries.append(dict(caller=caller,**dict(zip(['species','study_id','route','method','cohort'],key)),subset=subset,candidates=len(f),
                measurable=int(f.adjusted_call.ne('unmeasurable').sum()),retained=int(f.adjusted_call.eq('aligned').sum()),flipped=int(f.adjusted_call.eq('flipped').sum()),
                lost_alignment=int(f.adjusted_call.eq('unaligned').sum()),unmeasurable=int(f.adjusted_call.eq('unmeasurable').sum())))
    save(summaries,caller+'/pm_candidate_summary.tsv')
    identities=pairs[['species','study_id','route','library_key_a','sample_unit_a','replicate_status_a','time_a']].drop_duplicates()
    save(identities,caller+'/pm_library_identities.tsv')
    return formal,max_error

def do_dual(caller,base,studies,primary,own):
    trip=read(base/'matched_pms/weighted_triplets.tsv.gz');trip=trip[trip.setting.eq('supported50_p08')&trip.study_id.isin(STUDIES[:3])&trip.gene_id.isin(set(primary.gene_id)|set(own.gene_id))]
    rows=[];state_diag=[];pair_diag=[];max_error=0.
    for (sid,gene),q in trip.groupby(['study_id','gene_id'],sort=False):
        study=studies[sid];data,axes=study.gene(gene)
        for r in q.itertuples():
            ix=pd.Index(axes.atlas_pas_id).get_indexer(json.loads(r.PAS_ids));assert (ix>=0).all();d=data[:,ix];linear=axes.WPP.to_numpy()[ix]
            states={name:study.state(r.route,lib,name,d) for name,lib in [('P',r.library_key_a),('M',r.library_key_a),('S',r.library_key_b)]}
            for name,col in [('P','target_usage'),('M','M_usage'),('S','reference_usage')]:
                err=float(np.max(np.abs(states[name]['mean']-json.loads(getattr(r,col)))));max_error=max(max_error,err);assert err<1e-10
            raw={n:float(s['mean']@linear) for n,s in states.items()}
            for method in METHODS:
                context=study.context(r.route,r.library_key_a,r.library_key_b,'PMS',method);u,status,audit,ss,pp=gene_adjust(states,context)
                row=dict(caller=caller,species=r.species,study_id=sid,route=r.route,library_key=r.library_key_a,source_library_key=r.library_key_b,gene_id=gene,gene_symbol=r.gene_symbol,
                    triplet_weight=r.triplet_weight,method=method,scope='PMS',status=status,**audit)
                for co,a,b in [('PS','P','S'),('PM','P','M'),('MS','M','S')]:
                    row['raw_all_'+co]=raw[a]-raw[b];row['raw_context_common_'+co]=raw[a]-raw[b] if u else np.nan
                    for version in ['raw_support_common','adjusted']:row[version+'_'+co]=float((u[version][a]-u[version][b])@linear) if u else np.nan
                if u:
                    for version in ['raw_context_common','raw_support_common','adjusted']:assert abs(row[version+'_PS']-row[version+'_PM']-row[version+'_MS'])<1e-12
                rows.append(row);add_diagnostics(row,ss,pp,states,state_diag,pair_diag)
    frame=pd.DataFrame(rows);save(frame,caller+'/dual_triplet_effects.tsv.gz');save(state_diag,caller+'/dual_gene_state_balance.tsv.gz');save(pair_diag,caller+'/dual_gene_pair_balance.tsv.gz')
    candidate_routes=[];candidate_libs=[];formal=[];summary=[];global_calls=[]
    for method in METHODS:
        fm=frame[frame.method.eq(method)]
        for cohort,points in [('fixed_primary',primary),('caller_own',own)]:
            for r in points.itertuples():
                z=fm[fm.species.eq(r.species)&fm.route.eq(r.route)&fm.gene_id.eq(r.gene_id)];good=z[z.status.eq('measurable')];sg=1 if r.recurrent_direction=='positive' else -1
                out=dict(caller=caller,species=r.species,study_id=r.study_id,route=r.route,gene_id=r.gene_id,gene_symbol=r.gene_symbol,method=method,cohort=cohort,
                    baseline_dual=bool(r.dual_qualified_gene),n_original_triplets=len(z),n_adjusted_triplets=len(good),original_call=call(r.DWUI_PS,r.DWUI_PM,sg))
                for version,part in [('raw_all',z),('raw_context_common',good),('raw_support_common',good),('adjusted',good)]:
                    w=part.triplet_weight.to_numpy(copy=True);w/=w.sum() if len(w) else 1
                    for co in ['PS','PM','MS']:out[version+'_'+co]=float(w@part[version+'_'+co]) if len(part) else np.nan
                    out[version+'_call']=call(out[version+'_PS'],out[version+'_PM'],sg) if len(part) else 'unmeasurable'
                if cohort=='caller_own':assert abs(out['raw_all_PS']-r.DWUI_PS)<1e-10 and abs(out['raw_all_PM']-r.DWUI_PM)<1e-10
                candidate_routes.append(out)
                for lib,g in good.groupby('library_key'):
                    w=g.triplet_weight.to_numpy(copy=True);w/=w.sum();rec=dict(caller=caller,species=r.species,study_id=r.study_id,route=r.route,library_key=lib,gene_id=r.gene_id,method=method,cohort=cohort)
                    margins={}
                    for version in ['raw_context_common','raw_support_common','adjusted']:
                        ps=float(w@g[version+'_PS']);pm=float(w@g[version+'_PM']);margins[version]=min(sg*ps,sg*pm)
                        rec[version+'_margin']=margins[version]
                    rec['margin_change']=margins['adjusted']-margins['raw_context_common'];rec['support_trim']=margins['raw_support_common']-margins['raw_context_common'];rec['reweight']=margins['adjusted']-margins['raw_support_common']
                    candidate_libs.append(rec)
    cr=pd.DataFrame(candidate_routes);cl=pd.DataFrame(candidate_libs)
    save(cr,caller+'/dual_candidate_route_calls.tsv');save(cl,caller+'/dual_candidate_gene_library_margins.tsv')
    libmeans=cl.groupby(['caller','species','study_id','route','library_key','method','cohort']).agg(n_genes=('gene_id','size'),mean_margin_change=('margin_change','mean'),
        mean_support_trim=('support_trim','mean'),mean_reweight=('reweight','mean'),raw_context_common_margin=('raw_context_common_margin','mean'),adjusted_margin=('adjusted_margin','mean')).reset_index() if len(cl) else pd.DataFrame(columns=['caller','species','study_id','route','library_key','method','cohort','n_genes','mean_margin_change','mean_support_trim','mean_reweight','raw_context_common_margin','adjusted_margin'])
    save(libmeans,caller+'/dual_library_margin_changes.tsv')
    for sp,sid,route in ROUTES[:3]:
        z=libmeans[libmeans.study_id.eq(sid)&libmeans.method.eq('continuous')&libmeans.cohort.eq('fixed_primary')]
        formal.append(dict(test_id='PMS_margin_change|'+caller+'|'+sid+'|'+route,family='A2-cycle-alignment-redesign',pas_source=caller,scope='PMS',species=sp,study_id=sid,route=route,
            method='continuous',cohort='fixed_primary',statistic='signed_min_PS_PM_adjusted_minus_raw_context_common',**inference(z.mean_margin_change,caller,'PMS',sid,route)))
    for key,z in cr.groupby(['species','study_id','route','method','cohort']):
        for subset,f in [('all_candidates',z),('baseline_dual',z[z.baseline_dual]),('original_route_aligned',z[z.original_call.eq('aligned')])]:
            summary.append(dict(caller=caller,**dict(zip(['species','study_id','route','method','cohort'],key)),subset=subset,candidates=len(f),measurable=int(f.adjusted_call.ne('unmeasurable').sum()),
                retained=int(f.adjusted_call.eq('aligned').sum()),flipped=int(f.adjusted_call.eq('flipped').sum()),lost_alignment=int(f.adjusted_call.eq('unaligned').sum()),unmeasurable=int(f.adjusted_call.eq('unmeasurable').sum())))
    for key,z in cr.groupby(['species','gene_id','method','cohort']):
        sp,gene,method,cohort=key;expected={'CIP','WOT'} if sp=='mouse' else {'GSE118258'};assert set(z.route)==expected
        rec=dict(caller=caller,species=sp,gene_id=gene,gene_symbol=z.gene_symbol.iloc[0],method=method,cohort=cohort,baseline_dual=bool(z.baseline_dual.all()))
        for version in ['raw_all','raw_context_common','raw_support_common','adjusted']:
            c=z[version+'_call'];rec[version+'_call']='unmeasurable' if c.eq('unmeasurable').any() else 'aligned' if c.eq('aligned').all() else 'flipped' if c.eq('flipped').all() else 'unaligned'
        global_calls.append(rec)
    save(summary,caller+'/dual_route_summary.tsv');gc=pd.DataFrame(global_calls);save(gc,caller+'/dual_species_calls.tsv')
    species_summary=[]
    for key,z in gc.groupby(['species','method','cohort']):
        for subset,f in [('all_candidates',z),('baseline_dual',z[z.baseline_dual])]:
            species_summary.append(dict(caller=caller,**dict(zip(['species','method','cohort'],key)),subset=subset,candidates=len(f),measurable=int(f.adjusted_call.ne('unmeasurable').sum()),
                retained=int(f.adjusted_call.eq('aligned').sum()),flipped=int(f.adjusted_call.eq('flipped').sum()),lost_alignment=int(f.adjusted_call.eq('unaligned').sum()),unmeasurable=int(f.adjusted_call.eq('unmeasurable').sum())))
    save(species_summary,caller+'/dual_species_summary.tsv');return formal,max_error

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--caller',choices=list(BASES));args=parser.parse_args()
    for p in [HERE/'PROTOCOL.md',HERE.parent/'PROTOCOL.md',HERE.parent/'CYCLE-FREEZE-ADDENDUM.md',HERE/'run_cycle.py',HERE/'freeze-receipt.json',HERE/'execution-receipt.json']:register(p)
    receipt=json.loads((HERE/'execution-receipt.json').read_text())
    for p,h in receipt['sha256'].items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
    archive=cycle_archive();primary_pm=recurrent(BASES['polyapipe']);primary_dual=dual_candidates(BASES['polyapipe'])
    assert primary_pm.groupby('species').gene_id.nunique().to_dict()=={'human':268,'mouse':420}
    assert primary_dual.groupby('species').gene_id.nunique().to_dict()=={'human':250,'mouse':395}
    tests=[];audit=[]
    with np.load(register(OLD/'pas_subsets.npz'),allow_pickle=False) as original,np.load(register(HERE/'additional_pas_subsets.npz'),allow_pickle=False) as additional:
        for pack,path in [(original,OLD/'extraction_manifest.json'),(additional,HERE/'extraction_manifest.json')]:
            manifest=json.loads(register(path).read_text());assert json.loads(str(pack['provenance_json']))=={r[k+'_path']:r[k+'_sha256'] for r in manifest['remote'] for k in ['counts','meta','full_meta']}
        for caller in [args.caller] if args.caller else BASES:
            base=BASES[caller];axes=read(base/'state_vectors/PAS_axes.primary.tsv.gz');own_pm=recurrent(base);own_dual=dual_candidates(base)
            studies={sid:Study(caller,sid,archive,original if sid in STUDIES[:3] else additional,axes) for sid in STUDIES}
            rows,pm_error=do_pm(caller,base,studies,primary_pm,own_pm);tests+=rows;print(caller,'PM complete',pm_error,flush=True)
            rows,dual_error=do_dual(caller,base,studies,primary_dual,own_dual);tests+=rows;print(caller,'PMS complete',dual_error,flush=True)
            for name,attr in [('context_weights.tsv.gz','weight_rows'),('context_solvers.tsv','context_rows'),('context_state_balance.tsv','context_state_rows'),('context_pair_balance.tsv','context_pair_rows')]:
                save(sum([getattr(s,attr) for s in studies.values()],[]),caller+'/'+name)
            audit.append(dict(caller=caller,max_original_PM_error=pm_error,max_original_PMS_usage_error=dual_error,context_count=sum(len(s.contexts) for s in studies.values())))
            del studies
    save(tests,'raw-tests.tsv');assert pd.DataFrame(tests).groupby('family').size().to_dict()==({'A2-cycle-alignment-redesign':8,'A8-cycle-PM-redesign':5} if args.caller else {'A2-cycle-alignment-redesign':16,'A8-cycle-PM-redesign':10})
    for p,h in INPUTS.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
    (HERE/'input_sha256.json').write_text(json.dumps(INPUTS,indent=2)+'\n')
    (HERE/'audit.json').write_text(json.dumps(dict(checks=audit,seed=20261008,bootstrap=10000,formal_positions=len(tests),post_hoc_method_redesign=True),indent=2)+'\n')
    print('Complete',audit,flush=True)

if __name__=='__main__':main()
