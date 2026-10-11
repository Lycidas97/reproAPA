"""Fresh matched marker-reference ortholog nulls under the root-frozen protocol."""
import argparse,concurrent.futures,hashlib,importlib.util,json,sys
from pathlib import Path
import numpy as np,pandas as pd
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];OLD=ROOT/'notes/2026-10-10-ortholog-threshold-scan'
def module(name,path):
 spec=importlib.util.spec_from_file_location(name,path);out=importlib.util.module_from_spec(spec);spec.loader.exec_module(out);return out
A=module('marker_frozen_ortholog',ROOT/'notes/2026-10-10-recurrence-utr-ortholog/analyze.py')
S=module('marker_frozen_scan',OLD/'analyze.py')
LABELS=A.LABELS
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p,**k:pd.read_csv(p,sep='\t',float_precision='round_trip',**k)
def save(f,n):f.to_csv(HERE/n,sep='\t',index=False,na_rep='NA',float_format='%.17g',compression={'method':'gzip','mtime':0} if n.endswith('.gz') else None)

def one_null(task):
 definition,threshold,model,part,pairs,feats,labels,n=task
 cells={sp:A.partitions(feats[sp],model) for sp in ['human','mouse']}
 prob={sp:A.class_probabilities(labels[sp],cells[sp]) for sp in ['human','mouse']}
 expected=prob['human'].T@prob['mouse']
 seed=A.seed('ortholog_direction','P_M',threshold,model)
 null=A.cross_draws(labels['human'],labels['mouse'],cells['human'],cells['mouse'],np.random.default_rng(seed),n=n)
 context=dict(definition=definition,contrast='P_M',threshold=threshold,model=model,measured_orthologs=len(pairs))
 diag=[];mob=[]
 for sp in ['human','mouse']:
  d,m=S.feature_diagnostics(feats[sp],labels[sp],cells[sp],context|dict(species=sp));diag.append(d);mob.append(m)
 rows=[]
 for i,hd in enumerate(LABELS):
  for j,md in enumerate(LABELS):
   row=context|dict(human_direction=hd,mouse_direction=md,seed=seed,**A.null_row(int(((labels['human']==i)&(labels['mouse']==j)).sum()),null[:,i,j],expected[i,j]))
   row.update(S.tail_precision(row['tail_count'],n));row['tested']=i<2 and j<2;row['raw_p_enrichment']=row['p_enrichment']
   row['inferential_status']='estimable' if row['inferentially_mobile'] else 'not_estimable_zero_null_variance'
   row['p_for_correction']=row['raw_p_enrichment'] if row['inferentially_mobile'] else 1.
   if not row['inferentially_mobile'] or not row['tested']:row['p_enrichment']=np.nan
   rows.append(row)
 return definition,threshold,model,pd.DataFrame(rows),pd.concat(diag,ignore_index=True),pd.DataFrame(mob),null

def corrections(frame):
 f=frame.copy();formal=f[f.tested]
 for n in ['q8','q72','q144','q288']:f[n]=np.nan
 for _,ii in formal.groupby(['definition','threshold','model']).groups.items():
  assert len(ii)==8;f.loc[ii,'q8']=A.bh(f.loc[ii,'p_for_correction'])
 for _,ii in formal.groupby(['definition','model']).groups.items():
  assert len(ii)==72;f.loc[ii,'q72']=A.bh(f.loc[ii,'p_for_correction'])
 for _,ii in formal.groupby('definition').groups.items():
  assert len(ii)==144;f.loc[ii,'q144']=A.bh(f.loc[ii,'p_for_correction'])
 unique=[]
 for key,p in formal.groupby(['contrast','threshold','model','human_direction','mouse_direction'],sort=True):
  if key[0]=='P_S':
   assert len(p)==3 and p.p_for_correction.nunique()==1;unique.append((tuple(key)+('fixed_P_S',),float(p.p_for_correction.iloc[0])))
  else:
   assert len(p)==3
   unique.extend((tuple(key)+(r.definition,),float(r.p_for_correction)) for r in p.itertuples())
 assert len(unique)==288
 q=A.bh([v for _,v in unique]);qm={k:float(v) for (k,_),v in zip(unique,q)}
 for ii,r in formal.iterrows():
  key=(r.contrast,r.threshold,r.model,r.human_direction,r.mouse_direction,'fixed_P_S' if r.contrast=='P_S' else r.definition);f.loc[ii,'q288']=qm[key]
 return f

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=4);args=parser.parse_args()
 protocol=json.loads((HERE.parent/'protocol.json').read_text());pf=json.loads((HERE.parent/'protocol_freeze_sha256.json').read_text())
 for n,h in pf.items():assert sha(HERE.parent/n)==h,n
 covfreeze=json.loads((HERE/'covariate_freeze.json').read_text());assert covfreeze['status']=='ACTUAL_COVARIATES_AND_EFFECTS_FROZEN_BEFORE_LABELS'
 for n,h in covfreeze['output_sha256'].items():assert sha(HERE/n)==h,n
 paths=[HERE.parent/'PROTOCOL.zh.md',HERE.parent/'protocol.json',HERE/'covariate_freeze.json',HERE/'study_covariates.tsv.gz',HERE/'fixed_background.tsv',OLD/'ortholog_3x3_threshold_scan.tsv',OLD/'null_draws.npz',OLD/'measured_ortholog_members.tsv.gz',ROOT/'notes/2026-10-10-recurrence-utr-ortholog/analyze.py',OLD/'analyze.py']
 hashes={str(p.relative_to(ROOT)):sha(p) for p in paths}
 studies=read(HERE/'study_covariates.tsv.gz');bg=read(HERE/'fixed_background.tsv')
 assert len(bg)==protocol['background_n']==728
 pairs=bg[['human_gene_id','mouse_gene_id','human_gene_symbol','mouse_gene_symbol','is_high_confidence','homology_type']].copy()
 assert pairs.human_gene_id.is_unique and pairs.mouse_gene_id.is_unique
 tasks=[];events=[];features=[];members=[]
 for definition,sources in protocol['definitions'].items():
  base=studies[[r.source_contrast==sources[r.species] for r in studies.itertuples()]].copy().sort_values(['species','contrast','unit','gene_id']).reset_index(drop=True)
  base['eligible_units']=base.groupby(['species','gene_id']).unit.transform('nunique')
  assert base.eligible_units.ge(2).all()
  fixed={}
  for sp in ['human','mouse']:
   source=base[base.species.eq(sp)].copy();source['can_negative']=False;source['can_positive']=False
   fixed[sp]=A.gene_features(source,pairs[sp+'_gene_id'].tolist())
  for t in protocol['thresholds']:
   part=S.threshold_events(base,t);events.append(part.assign(definition=definition))
   labels={};feats={}
   for sp in ['human','mouse']:
    sf=part[part.species.eq(sp)];ids=pairs[sp+'_gene_id'].tolist();feats[sp]=fixed[sp].copy()
    caps=sf.groupby('gene_id')[['can_negative','can_positive']].sum().reindex(ids)
    feats[sp]['proximal_capacity']=caps.can_negative.to_numpy(int);feats[sp]['distal_capacity']=caps.can_positive.to_numpy(int)
    # Cached quantiles and measured-pattern are identical to the original helper.
    check=A.gene_features(sf,ids)
    pd.testing.assert_frame_equal(feats[sp],check,check_exact=True)
    labels[sp]=A.class_labels(sf,ids)
    features.append(feats[sp].assign(definition=definition,contrast='P_M',species=sp,threshold=t))
   mem=S.paired_roster(pairs,part,labels,t).assign(definition=definition,contrast='P_M');members.append(mem)
   for model in protocol['models']:tasks.append((definition,t,model,part,pairs,feats,labels,protocol['permutations']))
 save(pd.concat(events,ignore_index=True),'study_events.tsv.gz');save(pd.concat(features,ignore_index=True),'ortholog_assay_features.tsv.gz');mem=pd.concat(members,ignore_index=True);save(mem,'measured_ortholog_members.tsv.gz')
 # Report the primary threshold first without changing any per-context seed or family.
 order={d:i for i,d in enumerate(protocol['definitions'])}
 tasks.sort(key=lambda task:(task[1]!=.02,order[task[0]],task[1],protocol['models'].index(task[2])))
 print('LABELS_AND_FEATURES_COMPLETE',len(tasks),'fresh_null_contexts',flush=True)
 rows=[];diagnostics=[];mobility=[];archive={}
 with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
  for definition,t,model,result,diag,mob,null in pool.map(one_null,tasks):
   rows.append(result);diagnostics.append(diag);mobility.append(mob);archive['|'.join(map(str,('ortholog','P_M',t,model,definition)))]=null
   print('FRESH_NULL_COMPLETE',definition,t,model,flush=True)
 # Fixed P–S observed/expected/raw p/null values are copied exactly; q is new-family-specific.
 ps=read(OLD/'ortholog_3x3_threshold_scan.tsv');ps=ps[ps.contrast.eq('P_S')].copy()
 assert len(ps)==9*2*9
 ps=ps.drop(columns=[n for n in ps if n.startswith('q_')])
 ps['raw_p_enrichment']=ps.p_enrichment;ps['p_for_correction']=ps.p_enrichment.fillna(1.);ps['inferential_status']='unchanged_fixed_calibration'
 for d in protocol['definitions']:rows.append(ps.assign(definition=d))
 full=corrections(pd.concat(rows,ignore_index=True));save(full,'ortholog_3x3_threshold_scan.tsv');save(full[full.tested],'directional_tests.tsv');save(full[full.contrast.eq('P_M')&full.threshold.eq(.02)],'PM_3x3_primary_threshold.tsv')
 save(pd.concat(diagnostics,ignore_index=True),'ortholog_stratum_diagnostics.tsv.gz');save(pd.concat(mobility,ignore_index=True),'mobility_summary.tsv')
 with np.load(OLD/'null_draws.npz') as old:
  for k in old.files:
   if k.split('|')[1]=='P_S':archive[k]=old[k]
 np.savez_compressed(HERE/'null_draws.npz',**archive)
 oldmem=read(OLD/'measured_ortholog_members.tsv.gz');core=oldmem[oldmem.contrast.eq('P_M')&oldmem.threshold.eq(.02)&oldmem.human_direction.eq('proximal')&oldmem.mouse_direction.eq('proximal')][['human_gene_id','mouse_gene_id','human_gene_symbol','mouse_gene_symbol']]
 assert len(core)==27
 flows=[]
 for (d,t),part in mem.groupby(['definition','threshold']):
  flow=core.merge(part.drop(columns=['human_gene_symbol','mouse_gene_symbol','is_high_confidence','homology_type']),on=['human_gene_id','mouse_gene_id'],how='left',validate='one_to_one')
  flow['definition']=d;flow['threshold']=t;flow['measurable']=flow.human_direction.notna()
  flow['status']=np.where(~flow.measurable,'not_estimable_background_loss',np.where(flow.human_direction.eq('proximal')&flow.mouse_direction.eq('proximal'),'retained_shared_proximal','measurable_changed_class'))
  flows.append(flow)
 save(pd.concat(flows,ignore_index=True),'original_core27_flow.tsv')
 primary=mem[mem.threshold.eq(.02)&mem.human_direction.eq('proximal')&mem.mouse_direction.eq('proximal')].copy();primary['original_core27']=primary.human_gene_id.isin(core.human_gene_id);save(primary,'shared_proximal_primary_members.tsv')
 for p,h in hashes.items():assert sha(ROOT/p)==h,p
 outputs=['study_events.tsv.gz','ortholog_assay_features.tsv.gz','measured_ortholog_members.tsv.gz','ortholog_3x3_threshold_scan.tsv','directional_tests.tsv','PM_3x3_primary_threshold.tsv','ortholog_stratum_diagnostics.tsv.gz','mobility_summary.tsv','null_draws.npz','original_core27_flow.tsv','shared_proximal_primary_members.tsv']
 manifest={'status':'COMPLETE_AWAITING_INDEPENDENT_QA','input_sha256':hashes,'output_sha256':{n:sha(HERE/n) for n in outputs},'script_sha256':sha(Path(__file__)),'permutations':protocol['permutations'],'fresh_PM_contexts':len(tasks),'unique_formal_family':288,'P_S_raw_values_and_draws_unchanged':True,'full_assay_actual_GEX':True,'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__}
 (HERE/'analysis_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print('ANALYSIS_COMPLETE_AWAITING_QA',flush=True)
if __name__=='__main__':main()
