"""Independent audit from frozen raw tables; no production implementation imports."""
from pathlib import Path
from itertools import combinations, product
from collections import Counter, defaultdict
from fractions import Fraction
import math, json, hashlib, argparse
import numpy as np
import pandas as pd

BASE = Path('notes/2026-10-10-conserved-response')
OUT = BASE / 'review'
RES = Path('notes/2026-10-03-polyapipe-rerun/statistics_package_v1/resources')
ROOTS = {'GO:0003723':'molecular_function','GO:0006397':'biological_process','GO:0006417':'biological_process'}
CHECKS=[]
def check(label,condition,detail=None):
    CHECKS.append({'check':label,'pass':bool(condition),'detail':detail})
    if not condition: raise AssertionError(label+': '+str(detail))
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p): return pd.read_csv(p,sep='\t',low_memory=False)
def writejson(p,x): p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n')
def bh(p):
    a=np.asarray(p,dtype=float); order=np.argsort(a,kind='stable'); ranks=np.arange(1,len(a)+1)
    v=np.minimum.accumulate((a[order]*len(a)/ranks)[::-1])[::-1]; q=np.empty(len(a));q[order]=np.minimum(v,1);return q

def pmf_float(N,K,n):
    v=np.zeros(n+1); den=math.comb(N,n)
    for k in range(max(0,n-(N-K)),min(n,K)+1):v[k]=math.comb(K,k)*math.comb(N-K,n-k)/den
    return v

def exact_distribution(layers):
    pmf=np.array([1.0])
    for N,K,n in layers: pmf=np.convolve(pmf,pmf_float(N,K,n))
    return pmf

def toy():
    trials=[[(4,2,2)],[(4,2,2),(3,1,1)],[(3,3,2),(4,1,2),(1,0,1)],[(2,1,0),(1,1,1)],[(5,0,3),(2,2,2)]]
    tested=[]
    for layers in trials:
        choices=[]
        for N,K,n in layers:
            choices.append([sum(i<K for i in c) for c in combinations(range(N),n)])
        outcomes=Counter(sum(x) for x in product(*choices)); den=math.prod(map(len,choices))
        exhaustive={k:Fraction(v,den) for k,v in outcomes.items()}
        pmf=exact_distribution(layers)
        check('toy exhaustive PMF '+str(layers),all(abs(pmf[k]-float(exhaustive.get(k,0)))<2e-15 for k in range(len(pmf))))
        mean=sum(n*K/N for N,K,n in layers)
        var=sum(n*(K/N)*(1-K/N)*(N-n)/(N-1) if N>1 else 0 for N,K,n in layers)
        check('toy mean '+str(layers),abs(np.dot(np.arange(len(pmf)),pmf)-mean)<2e-14)
        check('toy variance '+str(layers),abs(np.dot((np.arange(len(pmf))-mean)**2,pmf)-var)<2e-14)
        for k in range(len(pmf)):
            check('toy tail '+str(layers)+' '+str(k),abs(pmf[k:].sum()-float(sum(v for j,v in exhaustive.items() if j>=k)))<2e-14)
        tested.append({'layers':layers,'exhaustive_labels':den,'pmf':{str(k):str(v) for k,v in exhaustive.items()}})
    check('BH includes NE p=1',np.allclose(bh([.001,.03,1,.06]),[.004,.06,1,.08]))
    check('BH ties',np.allclose(bh([.02,.02,.3,1]),[.04,.04,.4,1]))
    # A namespace-unknown case must leave the pre-fixed bins intact.
    bins=[0,0,1,1,2]; known=[True,False,True,True,False]
    check('unknown filter keeps original strata',[b for b,k in zip(bins,known) if k]==[0,1,1])
    check('union has one vote',len(set(['g1','g2'])|set(['g2','g3']))==3)
    # Opposite and cross-contrast two-side recurrence are not single-side controls.
    for pair,expected in [([('proximal','distal')],False),([('proximal','nonrecurrent'),('nonrecurrent','proximal')],False),([('proximal','nonrecurrent')],True)]:
        check('control eligibility toy '+str(pair),(any(h!='nonrecurrent' for h,m in pair)^any(m!='nonrecurrent' for h,m in pair))==expected)
    writejson(OUT/'toy_validation.json',{'status':'PASS','trials':tested,'checks':CHECKS.copy()})

def reconstruct_main(threshold):
    m=read('notes/2026-10-10-recurrence-utr-ortholog/measured_ortholog_members.tsv.gz')
    m=m[m.threshold==threshold]
    check('unique measured pair contrast '+str(threshold),not m.duplicated(['human_gene_id','mouse_gene_id','contrast']).any())
    f=read('notes/2026-10-10-recurrence-utr-ortholog/ortholog_assay_features.tsv.gz'); f=f[f.threshold==threshold]
    check('unique features '+str(threshold),not f.duplicated(['gene_id','species','contrast']).any())
    fl=f.set_index(['gene_id','species','contrast'])
    result=[]
    traits=['mean_GEX','APA_coverage','PAS_number','baseline','proximal_capacity','distal_capacity']
    for (hg,mg),g in m.groupby(['human_gene_id','mouse_gene_id'],sort=True):
        row={'human_gene_id':hg,'mouse_gene_id':mg,'human_gene_symbol':g.human_gene_symbol.iloc[0],'mouse_gene_symbol':g.mouse_gene_symbol.iloc[0],'threshold':threshold}
        ds=[];hr=[];mr=[]
        for c in ['P_S','P_M']:
            q=g[g.contrast==c]; row[c+'_measurable']=len(q)==1
            if len(q):
                h=str(q.human_direction.iloc[0]);mm=str(q.mouse_direction.iloc[0]);ds.append((c,h,mm));hr.append(h);mr.append(mm)
            else:h=mm=None
            row[c+'_human_direction']=h;row[c+'_mouse_direction']=mm
        row['case_all_four']=any(h==mm and h!='nonrecurrent' for c,h,mm in ds)
        row['case_first_three']=any(h==mm and h!='nonrecurrent' and not(c=='P_M' and h=='distal') for c,h,mm in ds)
        row['human_recurrent_any']=any(x!='nonrecurrent' for x in hr);row['mouse_recurrent_any']=any(x!='nonrecurrent' for x in mr)
        row['control']=row['human_recurrent_any'] ^ row['mouse_recurrent_any']
        row['control_side']='human' if row['control'] and row['human_recurrent_any'] else ('mouse' if row['control'] else None)
        dirs=set(hr+mr)-{'nonrecurrent'};row['direction_class']='both' if len(dirs)==2 else (next(iter(dirs)) if dirs else 'neither')
        row['measurable_mask']=''.join(str(int(row[c+'_measurable'])) for c in ['P_S','P_M'])
        for sp,gid in [('human',hg),('mouse',mg)]:
            crows=[fl.loc[(gid,sp,c)] for c,h,mm in ds]
            for t in traits: row[sp+'_'+t]=np.nanmedian([x[t] for x in crows])
            for c in ['P_S','P_M']:row[sp+'_'+c+'_measured_pattern']=str(fl.loc[(gid,sp,c),'measured_pattern']) if row[c+'_measurable'] else 'NE'
        result.append(row)
    allrows=pd.DataFrame(result)
    return allrows

def add_bins(d,is_naive):
    d=d.copy()
    for t in ['mean_GEX','APA_coverage','PAS_number']:
        h=d['human_'+t].rank(method='average',pct=True);m=d['mouse_'+t].rank(method='average',pct=True)
        for op,rank in [('min',pd.concat([h,m],axis=1).min(axis=1,skipna=False)),('max',pd.concat([h,m],axis=1).max(axis=1,skipna=False))]:
            d[op+'_'+t+'_rank']=rank
            cuts=np.quantile(rank.dropna(),[1/3,2/3],method='linear')
            d[op+'_'+t+'_bin']=[int(np.searchsorted(cuts,x,side='left')) if pd.notna(x) else -1 for x in rank]
    cols=['measurable_mask','direction_class']+['min_'+t+'_bin' for t in ['mean_GEX','APA_coverage','PAS_number']]
    d['coarse_stratum']=d[cols].astype(str).agg('|'.join,axis=1)
    extra=['max_'+t+'_bin' for t in ['mean_GEX','APA_coverage','PAS_number']] if is_naive else ['human_P_S_measured_pattern','human_P_M_measured_pattern','mouse_P_S_measured_pattern','mouse_P_M_measured_pattern']
    d['fine_stratum']=d[['coarse_stratum']+extra].astype(str).agg('|'.join,axis=1)
    return d

def reconstruct_naive(threshold,set_name):
    d=read('notes/2026-10-10-proximal-module-specificity/joint_membership_and_covariates.tsv.gz')
    target='P_M_recurrent' if set_name=='hN79' else 'P_M_PS_same_contributor_dual'
    d=d[(d.block=='naive_direction')&(d.mouse_target==target)&(d.threshold==threshold)].copy()
    check('unique naive raw '+set_name+str(threshold),not d.duplicated(['human_gene_id','mouse_gene_id']).any())
    hp=d.pooled_hN_P7_proximal;hd=d.pooled_hN_P7_distal;mp=d.proximal;md=d.distal
    check('known naive flags '+set_name+str(threshold),not pd.concat([hp,hd,mp,md],axis=1).isna().any().any())
    case=hp&mp;control=(hp & ~mp & ~md)| (mp & ~hp & ~hd)
    d['case']=case;d['control']=control;d=d[case|control].copy()
    d['human_gene_symbol']=d.pooled_hN_P7_h_gene_symbol;d['mouse_gene_symbol']=d.m_gene_symbol
    rename={'mean_GEX':('pooled_hN_P7_GEX','m_GEX'),'APA_coverage':('pooled_hN_P7_coverage','m_coverage'),'PAS_number':('pooled_hN_P7_PAS','m_PAS'),'baseline':('pooled_hN_P7_baseline','m_baseline'),'proximal_capacity':('pooled_hN_P7_capacity_proximal','m_capacity_proximal'),'distal_capacity':('pooled_hN_P7_capacity_distal','m_capacity_distal')}
    for t,(h,m) in rename.items():d['human_'+t]=d[h];d['mouse_'+t]=d[m]
    d['direction_class']='proximal';d['measurable_mask']='hN_PM' if set_name=='hN79' else 'hN_PM_PS';d['control_side']=np.where(d.control,np.where(d.pooled_hN_P7_proximal,'human','mouse'),None)
    d['role']=np.where(d.case,'case','control');d['set_name']=set_name
    return add_bins(d,True)

def reconstruct_pools():
    rows=[];allrows=[]
    for t in [.02,.05]:
        a=reconstruct_main(t);allrows.append(a)
        for name,c in [('main48','case_all_four'),('list45','case_first_three')]:
            d=a[a[c]|a.control].copy();d['case']=d[c];d['role']=np.where(d.case,'case','control');d['set_name']=name;rows.append(add_bins(d,False))
        for name in ['hN79','dual46']:rows.append(reconstruct_naive(t,name))
    return pd.concat(rows,ignore_index=True),pd.concat(allrows,ignore_index=True)

def compare_frame(a,b,keys,label):
    a=a.set_index(keys).sort_index();b=b.set_index(keys).sort_index()
    check(label+' index',a.index.equals(b.index),{'independent':len(a),'production':len(b)})
    errors=[]
    for c in sorted(set(a.columns)&set(b.columns)):
        x=a[c];y=b[c]
        if pd.api.types.is_numeric_dtype(x) and not pd.api.types.is_bool_dtype(x):
            ok=np.isclose(pd.to_numeric(x),pd.to_numeric(y),rtol=1e-12,atol=1e-12,equal_nan=True)
        else:ok=x.fillna('__NA__').astype(str).values==y.fillna('__NA__').astype(str).values
        if not np.all(ok):errors.append({'column':c,'n':int((~ok).sum()),'examples':[str(v) for v in a.index[~ok][:3]]})
    check(label+' all shared columns',not errors,errors)

def feasibility(pool):
    result=[]
    for (name,t),g in pool.groupby(['set_name','threshold']):
        for model in ['coarse','fine']:
            grp=g.groupby(model+'_stratum').case.agg(['size','sum']);mobile=grp[(grp['sum']>0)&(grp['sum']<grp['size'])]
            result.append({'set_name':name,'threshold':t,'model':model,'N':len(g),'cases':int(g.case.sum()),'controls':int((~g.case).sum()),'strata':len(grp),'mobile_strata':len(mobile),'mobile_cases':int(mobile['sum'].sum()),'mobile_case_fraction':mobile['sum'].sum()/g.case.sum(),'expected_changed_fraction':sum(2*x['sum']*(x['size']-x['sum'])/x['size'] for _,x in grp.iterrows())/len(g)})
    return pd.DataFrame(result)

def load_go():
    terms={};cur=None
    for line in (RES/'go-basic.obo').read_text().splitlines()+['[End]']:
        if line.startswith('['):
            if cur and cur.get('id'):terms[cur['id']]=cur
            cur={'parents':[],'alt':[]} if line=='[Term]' else None
        elif cur is not None:
            if line.startswith('id: '):cur['id']=line[4:]
            elif line.startswith('namespace: '):cur['namespace']=line[11:]
            elif line.startswith('name: '):cur['name']=line[6:]
            elif line.startswith('is_a: '):cur['parents'].append(line.split()[1])
            elif line.startswith('relationship: part_of '):cur['parents'].append(line.split()[2])
            elif line.startswith('alt_id: '):cur['alt'].append(line[8:])
            elif line=='is_obsolete: true':cur['obsolete']=True
    valid={k:v for k,v in terms.items() if not v.get('obsolete',False)};aliases={a:k for k,v in valid.items() for a in v['alt']};memo={}
    def ancestors(t):
        t=aliases.get(t,t)
        if t not in valid:return frozenset()
        if t not in memo:memo[t]=frozenset({t}).union(*(ancestors(p) for p in valid[t]['parents']))
        return memo[t]
    annotations={};counts={}
    for sp in ['human','mouse']:
        genes=defaultdict(set);invalid=blank=0;artifact_rows=[]
        for line in (RES/(sp+'_go_ids.tsv')).read_text().splitlines():
            if '\t' not in line:
                artifact_rows.append(line)
                check('GO export footer exact '+sp,line=='[success]')
                continue
            gid,raw=line.split('\t',1);t=aliases.get(raw,raw)
            if not raw:blank+=1;continue
            if t not in valid:invalid+=1;continue
            genes[gid].update(ancestors(t))
        annotations[sp]=genes;counts[sp]={'annotated_genes':len(genes),'blank_rows':blank,'invalid_or_obsolete_rows':invalid,'non_data_export_footer_rows':artifact_rows}
    for r,ns in ROOTS.items():check('GO root namespace '+r,valid[r]['namespace']==ns)
    # Literal descendant checks distinguish regulates from is_a/part_of propagation.
    excluded_examples=[]
    for k,v in valid.items():
        if any(p not in valid for p in v['parents']):excluded_examples.append(k)
    writejson(OUT/'go_resource_validation.json',{'valid_terms':len(valid),'obsolete_terms':len(terms)-len(valid),'alt_ids':len(aliases),'genes':counts,'parents_not_valid':excluded_examples})
    return valid,annotations

def gene_go_rows(pool,terms,ann):
    pairs=pool[['human_gene_id','mouse_gene_id']].drop_duplicates();rows=[]
    for row in pairs.itertuples(index=False):
        for root,ns in ROOTS.items():
            entry={'human_gene_id':row.human_gene_id,'mouse_gene_id':row.mouse_gene_id,'category':root}
            for sp,gid in [('human',row.human_gene_id),('mouse',row.mouse_gene_id)]:
                gs=ann[sp].get(gid,set());entry[sp+'_known']=any(terms[t]['namespace']==ns for t in gs);entry[sp+'_member']=root in gs
            rows.append(entry)
    return pd.DataFrame(rows)

def input_audit():
    p=json.loads((BASE/'protocol.json').read_text());
    check('protocol text frozen sha',sha(BASE/'PROTOCOL.zh.md')=='3884344b47e085974979f6f8efa364dd92c8b403593b8bccbb739f479cd65e3c')
    check('protocol metadata frozen sha',sha(BASE/'protocol.json')=='f70aa9ad1cb80c9db020426c2ab76e1969da9fe3ece494f236d825f4668a4b48')
    for path,s in p['inputs_sha256'].items():check('frozen input '+path,sha(path)==s)
    pool,allrows=reconstruct_pools()
    compare_frame(allrows,read(BASE/'function/all_measured_ortholog_classification.tsv.gz').assign(measurable_mask=lambda d:d.measurable_mask.astype(str).str.zfill(2)),['threshold','human_gene_id','mouse_gene_id'],'all measured classification')
    prod=pd.read_csv(BASE/'function/comparison_membership_features.tsv.gz',sep='\t',low_memory=False,dtype={'measurable_mask':str})
    compare_frame(pool,prod,['set_name','threshold','human_gene_id','mouse_gene_id'],'comparison membership/covariates/bins')
    f=feasibility(pool);compare_frame(f,read(BASE/'function/matching_feasibility.tsv'),['set_name','threshold','model'],'matching feasibility')
    pool.to_csv(OUT/'independent_membership_features.tsv.gz',sep='\t',index=False);f.to_csv(OUT/'independent_matching_feasibility.tsv',sep='\t',index=False)
    terms,ann=load_go();go=gene_go_rows(pool,terms,ann);go.to_csv(OUT/'independent_gene_go.tsv.gz',sep='\t',index=False)
    summary=pool.groupby(['set_name','threshold']).case.agg(['size','sum']).reset_index().to_dict('records')
    writejson(OUT/'input_validation.json',{'status':'PASS','summary':summary,'checks':CHECKS.copy(),'production_membership_sha256':sha(BASE/'function/comparison_membership_features.tsv.gz')})
    return pool,terms,ann,go

NAMES={'GO:0003723':'RNA binding','GO:0006397':'mRNA processing','GO:0006417':'regulation of translation'}
def conditional(g,known,member,stratum_col):
    use=g.loc[known].copy();use['member']=member.loc[known].astype(bool)
    margins=[]
    for s,layer in use.groupby(stratum_col,sort=True):
        margins.append({'stratum':s,'N':len(layer),'K':int(layer.member.sum()),'n':int(layer.case.sum()),'hits':int((layer.member&layer.case).sum())})
    layers=[(r['N'],r['K'],r['n']) for r in margins];pmf=exact_distribution(layers)
    check('PMF normalized',abs(pmf.sum()-1)<3e-13)
    hits=int((use.member&use.case).sum());expected=sum(n*K/N for N,K,n in layers)
    var=sum(n*(K/N)*(1-K/N)*(N-n)/(N-1) if N>1 else 0 for N,K,n in layers)
    check('real PMF expectation',abs(np.dot(np.arange(len(pmf)),pmf)-expected)<2e-11)
    check('real PMF variance',abs(np.dot((np.arange(len(pmf))-expected)**2,pmf)-var)<2e-11)
    mobile=[(N,K,n) for N,K,n in layers if 0<n<N]
    status='TESTED' if use.case.sum()>0 and 0<use.member.sum()<len(use) and var>1e-14 else 'NE_no_conditional_variation_or_empty'
    p=min(1.,float(pmf[hits:].sum())) if status=='TESTED' else 1.
    cdf=np.cumsum(pmf)
    result={'N_before_annotation':len(g),'cases_before_annotation':int(g.case.sum()),'N':len(use),'cases':int(use.case.sum()),'controls':int((~use.case).sum()),'annotation_unknown_cases':int((g.case&~known).sum()),'annotation_unknown_controls':int((~g.case&~known).sum()),'K':int(use.member.sum()),'hits':hits,'expected':expected,'matched_fold':hits/expected if expected else np.nan,'excess':hits-expected,'null_variance':var,'null_p025':int(np.searchsorted(cdf,.025)),'null_p975':int(np.searchsorted(cdf,.975)),'p':p,'status':status,'strata':len(layers),'mobile_cases':sum(n for N,K,n in mobile),'mobile_case_fraction':sum(n for N,K,n in mobile)/use.case.sum() if use.case.sum() else np.nan,'expected_changed_fraction':sum(2*n*(N-n)/N for N,K,n in layers)/len(use) if len(use) else np.nan}
    return result,margins,pmf

def final_audit(pool,terms,ann,go):
    # Every production ancestor membership is compared with independent OBO closure.
    published=read(BASE/'function/annotation_gene_membership.tsv.gz')
    for r in published.itertuples(index=False):
        own=ann[r.species].get(r.gene_id,set())
        check('GO propagated terms '+r.species+r.gene_id,own==set(json.loads(r.propagated_terms)))
        for ns in ['molecular_function','biological_process']:
            check('GO known '+r.species+r.gene_id+ns,any(terms[t]['namespace']==ns for t in own)==getattr(r,ns+'_known'))
    check('annotation gene count',len(published)==2*pool[['human_gene_id','mouse_gene_id']].drop_duplicates().shape[0])
    # Reverse the graph independently to enumerate all root descendants.
    children=defaultdict(set)
    for t,v in terms.items():
        for p in v['parents']:children[p].add(t)
    desc=[]
    for root in ROOTS:
        seen=set();todo=[root]
        while todo:
            t=todo.pop()
            if t not in seen:seen.add(t);todo.extend(children[t])
        for t in sorted(seen):desc.append({'root_id':root,'descendant_id':t,'descendant_name':terms[t]['name'],'namespace':terms[t]['namespace']})
    compare_frame(pd.DataFrame(desc),read(BASE/'function/root_descendants.tsv'),['root_id','descendant_id'],'all root descendants')
    rows=[];members=[];marginrows=[];pmfrows=[];balance=[]
    for (name,t),raw in pool.groupby(['set_name','threshold'],sort=True):
        for annotation,root,model in product(['human','mouse','AND'],ROOTS,['coarse','fine']):
            g=raw.merge(go[go.category==root],on=['human_gene_id','mouse_gene_id'],validate='one_to_one')
            if annotation=='AND':known=g.human_known&g.mouse_known;member=g.human_member&g.mouse_member
            else:known=g[annotation+'_known'];member=g[annotation+'_member']
            stats,margins,pmf=conditional(g,known,member,model+'_stratum')
            ident='|'.join([name,str(t),annotation,root,NAMES[root],model])
            stats.update(set_name=name,threshold=t,annotation=annotation,category=root,category_name=NAMES[root],model=model,test_id=ident,primary=name in ['main48','hN79'] and t==.02 and annotation=='human' and model=='coarse')
            rows.append(stats)
            use=g.loc[known].copy(); weights=use.groupby(model+'_stratum').case.transform('mean')
            for species,trait in product(['human','mouse'],['mean_GEX','APA_coverage','PAS_number','baseline','proximal_capacity','distal_capacity']):
                x=use[species+'_'+trait];finite=x.notna();x=x.loc[finite];w=weights.loc[finite];case=use.case.loc[finite]
                obs=x.loc[case].mean();ctrl=x.loc[~case].mean();exp=(x*w).sum()/w.sum();sd=x.std(ddof=1)
                balance.append({'set_name':name,'threshold':t,'annotation':annotation,'category':root,'category_name':NAMES[root],'model':model,'species':species,'covariate':trait,'observed_case_mean':obs,'original_control_mean':ctrl,'conditional_expected_case_mean':exp,'unmatched_case_control_SMD':(obs-ctrl)/sd if sd>0 else np.nan,'conditional_residual_SMD':(obs-exp)/sd if sd>0 else np.nan,'finite_genes':len(x)})
            out=g[['human_gene_id','mouse_gene_id','case',model+'_stratum']].copy();out['annotation_known']=known;out['category_member']=member
            out['set_name']=name;out['threshold']=t;out['annotation']=annotation;out['category']=root;out['model']=model;members.append(out)
            for x in margins:marginrows.append({'test_id':ident,**x})
            for k in np.flatnonzero(pmf):pmfrows.append({'test_id':ident,'hit_count':int(k),'probability':float(pmf[k])})
    d=pd.DataFrame(rows);d['q_primary6']=np.nan;d['q_sensitivity138']=np.nan
    check('family 6 and 138',d.primary.sum()==6 and (~d.primary).sum()==138)
    d.loc[d.primary,'q_primary6']=bh(d.loc[d.primary,'p']);d.loc[~d.primary,'q_sensitivity138']=bh(d.loc[~d.primary,'p']);d['q_joint144']=bh(d.p)
    compare_frame(d,read(BASE/'function/category_results.tsv'),['test_id'],'144 statistics and 6/138/144 BH')
    compare_frame(pd.concat(members,ignore_index=True),read(BASE/'function/category_membership.tsv.gz'),['set_name','threshold','annotation','category','model','human_gene_id','mouse_gene_id'],'all category memberships/unknown/filter strata')
    compare_frame(pd.DataFrame(marginrows),read(BASE/'function/category_stratum_margins.tsv.gz'),['test_id','stratum'],'all conditional layer margins')
    compare_frame(pd.DataFrame(pmfrows),read(BASE/'function/category_exact_pmf.tsv.gz'),['test_id','hit_count'],'all exact PMF support and probabilities')
    compare_frame(pd.DataFrame(balance),read(BASE/'function/covariate_balance.tsv.gz'),['set_name','threshold','annotation','category','model','species','covariate'],'all original covariate means/conditional balance/residuals')
    pd.DataFrame(balance).to_csv(OUT/'independent_covariate_balance.tsv.gz',sep='\t',index=False)
    d.to_csv(OUT/'independent_category_results.tsv',sep='\t',index=False)
    pd.DataFrame(pmfrows).to_csv(OUT/'independent_category_exact_pmf.tsv.gz',sep='\t',index=False)
    # Exploratory scan is reconstructed from all valid terms; eligibility uses background K.
    scan=[]
    for name in ['main48','hN79']:
        g=pool[(pool.set_name==name)&(pool.threshold==.02)].copy()
        sets=[ann['human'].get(gid,set()) for gid in g.human_gene_id]
        for ns in ['molecular_function','biological_process']:
            known=pd.Series([any(terms[t]['namespace']==ns for t in gs) for gs in sets],index=g.index)
            counts=Counter(t for gs,k in zip(sets,known) if k for t in gs if terms[t]['namespace']==ns)
            for term,K in sorted(counts.items()):
                if not 10<=K<=500:continue
                member=pd.Series([term in gs for gs in sets],index=g.index)
                stats,_,_=conditional(g,known,member,'coarse_stratum');stats.update(set_name=name,term_id=term,term_name=terms[term]['name'],namespace=ns);scan.append(stats)
    scan=pd.DataFrame(scan);scan['q_set']=np.nan
    for name in ['main48','hN79']:sel=scan.set_name==name;scan.loc[sel,'q_set']=bh(scan.loc[sel,'p'])
    scan['q_joint_scan']=bh(scan.p)
    compare_frame(scan,read(BASE/'function/exploratory_GO.tsv.gz'),['set_name','term_id'],'exploratory BP/MF eligibility/background K/statistics/BH')
    scan.to_csv(OUT/'independent_exploratory_GO.tsv.gz',sep='\t',index=False)
    checks_file=OUT/'final_validation.json'
    hashes={str(p.relative_to(BASE)):sha(p) for p in sorted((BASE/'function').glob('*')) if p.is_file()}
    writejson(checks_file,{'status':'PASS','independent_algorithm':'stdlib math.comb and float convolution; exhaustive Fraction toy; independent OBO parser and reverse graph; no production imports','n_functional_tests':len(d),'n_primary':int(d.primary.sum()),'n_sensitivity':int((~d.primary).sum()),'n_exploratory_terms':len(scan),'checks':CHECKS,'sealed_function_file_sha256':hashes})
    return d,scan

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--phase',choices=['inputs','final'],default='inputs');args=parser.parse_args()
    OUT.mkdir(exist_ok=True);toy();pool,terms,ann,go=input_audit()
    if args.phase=='final':final_audit(pool,terms,ann,go)
    print(json.dumps({'status':'PASS','phase':args.phase,'checks':len(CHECKS),'pool_rows':len(pool),'orthologs':pool[['human_gene_id','mouse_gene_id']].drop_duplicates().shape[0]}))
