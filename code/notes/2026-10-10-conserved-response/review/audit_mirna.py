"""Independent strict native-core/inclusion and full physical-seed qualification."""
from pathlib import Path
from collections import defaultdict
import hashlib,json
import numpy as np
import pandas as pd
B=Path('notes/2026-10-10-conserved-response');S=B/'sites';O=B/'review';T=Path('notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables')
checks=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check(n,ok):
 checks.append({'check':n,'pass':bool(ok)})
 if not ok:raise AssertionError(n)
def read(p):return pd.read_csv(p,sep='\t',float_precision='round_trip',low_memory=False)
manifest=S/'mirna_sealed_manifest.json';check('announced miRNA manifest',sha(manifest)=='3187281ee84266da1b2357cacc0486f293b321c239a0f38bc0238832846b2668')
for f,h in json.loads(manifest.read_text())['files'].items():check('miRNA seal '+f,sha(S/f)==h)
intervals=read(S/'strict_case_intervals.tsv');effects=read(S/'strict_case_interval_inclusion.tsv');traces=read(S/'site_usage_traces.tsv.gz');contexts=read(S/'site_event_contexts.tsv.gz');mapping=read(T/'cis/detected_to_annotation.tsv.gz').set_index('native_pas_id').id.to_dict();tx=read(T/'tandem/annotation_terminal_utrs.tsv.gz')
routes=read(T/'endpoints/route_effects.tsv.gz');routes=routes[routes.setting.eq('supported50_p08')&((routes.contrast.eq('P_S')&routes.design.eq('true_source'))|(routes.contrast.eq('P_M')&routes.design.eq('same_library')))];routes=routes.set_index(['species','study_id','gene_id','contrast'],drop=False).sort_index()
naive=read(T/'naive_state/library_pair_effects.tsv.gz');naive=naive[naive.version.eq('v2_primary')&naive.setting.eq('supported50_p08')&naive.node_a.eq('N')&naive.node_b.eq('P')&naive.library_a.str.endswith(('6','8'))&naive.library_b.str.endswith('7')]
npz=Path('notes/2026-10-10-marker-reference-sensitivity/inputs.npz');check('native feature metadata immutable',sha(npz)=='6652eb94bc70ad1e361992d2991a4e24df9c7d0fe83e694db155ce156bbaccb6');native={}
with np.load(npz,allow_pickle=False) as a:
 for study in sorted(set(effects.study_id)):
  f=json.loads(str(a[study+'_features_json']));native[study]={x['atlas_pas_id']:x for x in f}
maxerr=0.;observedseed=[];strand_counts=intervals.strand.value_counts().to_dict()
for row in intervals.itertuples(index=False):
 match=contexts[contexts.source_population.eq(row.source_population)&contexts.gene_threshold.eq(.02)&contexts.tolerance.eq(row.tolerance)&contexts.contrast.eq(row.contrast)&contexts.human_gene_id.eq(row.human_gene_id)]
 check('interval actual shared strict case',len(match)==1 and match.parent_case.iloc[0] and match.context_role.iloc[0]=='shared_recurrent' and match.strict_double_driver.iloc[0])
 ctx=match.iloc[0];gain=ctx[row.species+'_gain_id'];loss=ctx[row.species+'_loss_id'];drivercuts=defaultdict(list);profiles=[]
 source=traces[traces.source_population.eq(row.source_population)&traces.gene_threshold.eq(.02)&traces.tolerance.eq(row.tolerance)&traces.contrast.eq(row.contrast)&traces.human_gene_id.eq(row.human_gene_id)&traces.species.eq(row.species)]
 for p in source.itertuples():
  if row.source_population=='hN79' and row.species=='human':
   raw=naive[naive.gene_id.eq(row.gene_id)&naive.library_a.eq(p.series)];check('native hN pair unique',len(raw)==1);raw=raw.iloc[0];u=np.asarray(json.loads(raw.usage_a));v=np.asarray(json.loads(raw.usage_b))
  else:
   c='P_M' if row.source_population=='hN79' else row.contrast;raw=routes.loc[(row.species,p.study_id,row.gene_id,c)]
   if isinstance(raw,pd.Series):raw=raw.to_frame().T
   if p.study_id=='GSE178325':raw=raw[raw.route.str.endswith(p.series.removeprefix('178325_'))]
   elif p.study_id in ['GSE103221','GSE115943']:raw=raw[raw.route.eq(p.series)]
   check('native actual endpoint unique',len(raw)==1);raw=raw.iloc[0];u=np.asarray(json.loads(raw.target_usage));v=np.asarray(json.loads(raw.reference_usage))
  ids=json.loads(raw.PAS_ids);cuts=[native[p.study_id][k]['cut_site'] for k in ids]
  for k,cut in zip(ids,cuts):
   if mapping.get(k) in [gain,loss]:drivercuts[mapping[k]].append(cut)
  profiles.append((p.study_id,p.series,np.asarray(cuts),u,v))
 # Identify which complete driver collection is proximal from transcription direction.
 a,b=drivercuts[gain],drivercuts[loss];check('native cuts for both maximum drivers',bool(a) and bool(b))
 aprox=(np.mean(a)<np.mean(b)) if row.strand=='+' else (np.mean(a)>np.mean(b));prox,dist=(a,b) if aprox else (b,a)
 lo,hi=(max(prox),min(dist)) if row.strand=='+' else (max(dist),min(prox))
 check('actual full native core boundaries',lo==row.start0 and hi==row.end0 and hi-lo==row.length and lo<hi)
 check('native boundary lists literal',sorted(set(prox))==json.loads(row.proximal_native_cuts) and sorted(set(dist))==json.loads(row.distal_native_cuts))
 fits=tx[tx.species.eq(row.species)&tx.gene_id.eq(row.gene_id)&tx.transcript_id.isin(row.common_transcripts.split(';'))]
 check('core entirely inside all actual common terminal UTRs',len(fits)==len(row.common_transcripts.split(';')) and fits.utr_start.sub(1).le(lo).all() and fits.utr_end.ge(hi).all())
 bystudy=defaultdict(list)
 for study,series,cuts,u,v in profiles:
  keep=cuts>=hi if row.strand=='+' else cuts<=lo;x=float(u[keep].sum());y=float(v[keep].sum());e=effects[effects.source_population.eq(row.source_population)&effects.tolerance.eq(row.tolerance)&effects.contrast.eq(row.contrast)&effects.human_gene_id.eq(row.human_gene_id)&effects.species.eq(row.species)&effects.study_id.eq(study)&effects.series.eq(series)]
  check('actual series inclusion row unique',len(e)==1);e=e.iloc[0];err=max(abs(x-e.target_inclusion),abs(y-e.reference_inclusion),abs(x-y-e.delta_inclusion));maxerr=max(maxerr,err);check('native unrenormalized per-series inclusion',err<2e-12);bystudy[study].append((x,y,x-y))
 mean=np.asarray([np.asarray(v).mean(axis=0) for v in bystudy.values()]).mean(axis=0);check('equal-study native inclusion summary',np.allclose(mean,[row.target_inclusion,row.reference_inclusion,row.delta_inclusion],rtol=0,atol=2e-12))
 physical=read(T/f'external_mirna/{row.species}_physical_seed_sites.tsv.gz');physical=physical[physical.gene_id.eq(row.gene_id)&physical.chrom.eq(row.chrom)&physical.strand.eq(row.strand)&physical.start0.ge(lo)&physical.end0.le(hi)]
 for h in physical.itertuples():observedseed.append((row.source_population,row.tolerance,row.contrast,row.human_gene_id,row.species,row.gene_id,h.seed,h.chrom,h.strand,h.start0,h.end0))
seeds=read(S/'strict_case_predicted_seeds.tsv');cols=['source_population','tolerance','contrast','human_gene_id','species','gene_id','seed','chrom','strand','start0','end0'];actual=list(seeds[cols].itertuples(index=False,name=None))
check('all complete frozen physical predictions in actual cores',sorted(observedseed)==sorted(actual));check('full physical widths 7 and 8 retained',(seeds.end0-seeds.start0).isin([7,8]).all());check('seed intervals not truncated',(seeds.start0>=seeds.core_start0).all() and (seeds.end0<=seeds.core_end0).all())
# Mapping is the previously sealed read-only raw chain result, not independently re-fetched here.
records={x['query_id']:x for x in map(json.loads,(S/'seed_chain_maps.jsonl').read_text().splitlines())};links=read(S/'reciprocal_predicted_seed_pairs.tsv');unique_seed=seeds.drop_duplicates('query_id').set_index('query_id')
for r in links.itertuples():
 h=unique_seed.loc[r.human_query_id];m=unique_seed.loc[r.mouse_query_id];forward=records[r.human_query_id];back=records[r.mouse_query_id]
 check('both maps individually unique/full physical width',forward['status']==back['status']=='unique' and len(forward['hits'])==len(back['hits'])==1)
 for hit,dest in [(forward['hits'][0],m),(back['hits'][0],h)]:check('whole exact reciprocal coordinate and strand',(hit['chrom'],hit['strand'],hit['start0'],hit['end0'])==(dest.chrom,dest.strand,int(dest.start0),int(dest.end0)))
 check('same predicted DNA group at corresponding whole physical sites',h.seed==m.seed==r.seed and h.end0-h.start0==m.end0-m.start0)
result={'status':'PASS','native_interval_rows':len(intervals),'native_series_inclusion_rows':len(effects),'native_interval_strands':strand_counts,'complete_predicted_seed_context_rows':len(seeds),'physical_widths':(seeds.end0-seeds.start0).value_counts().to_dict(),'reciprocal_pair_context_rows_checked':len(links),'max_actual_native_inclusion_error':maxerr,'chain_toy_validation_sha256':sha(O/'chain_toy_validation.json'),'sealed_mirna_sha256':sha(manifest),'scope_limitation':'No independent remote full-chain re-fetch; standalone production mapper tested on 11 independently expected fake-chain cases, read-only raw maps verified against exact full physical sites and actual native cores. Root independently verifies mapping closure/calibration/pair-summary counts. Predictions are not binding or functional validation.','checks':checks}
(O/'mirna_validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='checks'}))
