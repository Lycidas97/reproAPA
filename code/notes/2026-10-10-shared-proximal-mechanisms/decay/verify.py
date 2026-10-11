"""Independent reconstruction: regex hits, interval sets, analytic null and exact discrete tails."""
from pathlib import Path
import hashlib,json,math,re
from collections import defaultdict
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
TABLES=ROOT/'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables'

def sha(path):
 with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def hitstarts(s,w):return [m.start() for m in re.finditer('(?='+re.escape(w)+')',s)]
def independent_class(a,b,chrom,strand,utrs,model):
 # Convert interbase cuts to the retained 1-based nucleotide; classify using 1-based inclusive intervals.
 pos=[int(a) if strand=='+' else int(a)+1,int(b) if strand=='+' else int(b)+1];lo,hi=min(pos),max(pos)
 us=utrs.get((chrom,strand),[]);m=model.get((chrom,strand),[])
 if any(start<=lo and hi<=end for _,start,end in us):return 'strict_tandem_UTR'
 ex=[(t,s+1,e) for t,k,s,e in m if k=='exon'];tx=[(t,s+1,e) for t,k,s,e in m if k=='transcript']
 bytx=defaultdict(list)
 for t,s,e in ex:bytx[t].append((s,e))
 term=[]
 for t,v in bytx.items():
  selected=max(v) if strand=='+' else min(v);term.append((t,*selected))
 if not m or not term:return 'unknown'
 if any(s<=lo and hi<=e for _,s,e in term):return 'terminal_exon_nonUTR'
 anyboth=any(s<=lo and hi<=e for _,s,e in ex)
 allterminal=all(any(s<=p<=e for _,s,e in term) for p in pos)
 if allterminal and not anyboth:return 'alternative_terminal_exon_compatible'
 if any(not any(s<=p<=e for _,s,e in ex) and any(s<=p<=e for _,s,e in tx) for p in pos):return 'intronic_compatible'
 if any(any(s<=p<=e for _,s,e in term) for p in pos) or anyboth:return 'mixed'
 return 'unknown'

def exact_binary(g,nforeground,end):
 dist=np.array([1.])
 for _,sub in g.groupby('stratum_'+end[0]):
  n=len(sub);k=int(sub.foreground_002.sum());h=int(sub[end[1]].sum());component=np.zeros(k+1)
  for a in range(max(0,k-(n-h)),min(k,h)+1):component[a]=math.comb(h,a)*math.comb(n-h,k-a)/math.comb(n,k)
  dist=np.convolve(dist,component)
 return dist

def main():
 protocol=json.loads((HERE/'protocol.json').read_text());assert sha(HERE/'PROTOCOL.zh.md')==protocol['protocol_sha256']
 inputs=json.loads((HERE/'input_hashes.json').read_text());inputs.update(json.loads((HERE/'additional_source_hashes.json').read_text()));assert all(sha(ROOT/k)==v for k,v in inputs.items())
 seq=pd.read_csv(HERE/'full_aUTR_sequences.tsv',sep='\t');features=pd.read_csv(HERE/'species_motif_features.tsv.gz',sep='\t').set_index(['species','gene_id']);hits=pd.read_csv(HERE/'motif_hits.tsv.gz',sep='\t')
 motifs={'ESCC_7mer_m8':'AGCACTT','ESCC_7mer_A1':'GCACTTA','ESCC_8mer':'AGCACTTA','U4':'TTTT'};ind_hits=[]
 for row in seq.to_dict('records'):
  s=row['sequence'];f=features.loc[(row['species'],row['gene_id'])];assert len(s)==row['end0']-row['start0'];assert re.fullmatch('[ACGT]+',s)
  assert abs(f.GC-(s.count('G')+s.count('C'))/len(s))<1e-15
  for key,w in motifs.items():
   offsets=hitstarts(s,w);assert int(f[key+'_count'])==len(offsets)
   for off in offsets:
    start=row['start0']+off if row['strand']=='+' else row['end0']-off-len(w)
    ind_hits.append((row['species'],row['gene_id'],key,off,int(start),int(start+len(w))))
  assert f.ESCC_ge7==(len(hitstarts(s,'AGCACTT'))+len(hitstarts(s,'GCACTTA'))>0)
  assert f.ESCC_8mer==(len(hitstarts(s,'AGCACTTA'))>0)
  assert abs(f.U4_density-len(hitstarts(s,'TTTT'))/max(1,len(s)-3))<1e-15
 stored=set(tuple(x) for x in hits[['species','gene_id','motif','aUTR_offset0','genomic_start0','genomic_end0']].itertuples(index=False,name=None));assert set(ind_hits)==stored and len(ind_hits)==len(stored)
 for sp,name in [('human','hsa-miR-302b-3p'),('human','hsa-miR-372-3p'),('mouse','mmu-miR-294-3p')]:
  file=ROOT/'notes/2026-10-03-polyapipe-rerun/statistics_package_v2r1/resources'/f'{sp}_mirna_annotations__{sp}_mirna_annotations.tsv';a=pd.read_csv(file,sep='\t');z=a[a.mirbase_exact_aliases.fillna('').str.split(';').apply(lambda v:name in v)];assert len(z)==1 and z.iloc[0].mature_sequence[1:8]=='AAGUGCU' and z.iloc[0].seed_2_8=='AAGUGCU'
 events=pd.read_csv(HERE/'route_event_inventory.tsv.gz',sep='\t');models=pd.read_csv(HERE/'background_transcript_exon_CDS.tsv.gz',sep='\t');u=pd.read_csv(TABLES/'tandem/annotation_terminal_utrs.tsv.gz',sep='\t')
 bymodel=defaultdict(lambda:defaultdict(list));byutr=defaultdict(lambda:defaultdict(list))
 for row in models.itertuples(index=False):bymodel[(row.species,row.gene_id)][(str(row.chrom).removeprefix('chr'),row.strand)].append((row.transcript_id,row.kind,row.start0,row.end0))
 for row in u.itertuples(index=False):byutr[(row.species,row.gene_id)][(str(row.chrom).removeprefix('chr'),row.strand)].append((row.transcript_id,row.utr_start,row.utr_end))
 for row in events.itertuples(index=False):
  k=(row.species,row.gene_id);cls=independent_class(row.gain_cut,row.loss_cut,str(row.chrom).removeprefix('chr'),row.strand,byutr[k],bymodel[k]) if row.gain_delta>0 and row.loss_delta<0 else 'unknown';assert row.driver_class==cls,(row.gene_id,row.study_id,cls,row.driver_class)
 effects=pd.read_csv(TABLES/'state_vectors/route_effects.primary.tsv.gz',sep='\t');effects=effects[effects.setting.eq('supported50_p08')&effects.contrast.eq('P_M')].set_index(['species','study_id','route','gene_id'])
 for row in events.itertuples(index=False):
  e=effects.loc[(row.species,row.study_id,row.route,row.gene_id)];vals=dict(zip(json.loads(e.PAS_ids),json.loads(e.delta_usage)));mx=max(vals.values());mn=min(vals.values());assert row.gain_pas==sorted(k for k,v in vals.items() if v==mx)[0];assert row.loss_pas==sorted(k for k,v in vals.items() if v==mn)[0]
 b=pd.read_csv(HERE/'paired_motif_background.tsv.gz',sep='\t');coverage=pd.read_csv(HERE/'paired_coverage.tsv.gz',sep='\t');assert len(coverage)==751 and coverage.foreground_002.sum()==27;assert len(b)==67 and b.foreground_002.sum()==8
 eligible=coverage.human_status.eq('eligible')&coverage.mouse_status.eq('eligible');assert set(coverage[eligible].human_gene_id)==set(b.human_gene_id);assert all(b.human_all_drivers_strict&b.mouse_all_drivers_strict)
 r=pd.read_csv(HERE/'enrichment_tests.tsv',sep='\t');draws=np.load(HERE/'null_draws.npz');fg=b.foreground_002;nr=int(fg.sum());endpoints=list(r[r.model.eq('coarse')].endpoint);comparisons=[]
 for model in ('coarse','full'):
  actual=b.loc[fg,endpoints].astype(float).mean().to_numpy();weights=b.groupby('stratum_'+model).foreground_002.transform('mean').to_numpy();Y=b[endpoints].astype(float).to_numpy();analytic=(weights[:,None]*Y).sum(axis=0)/nr
  for j,name in enumerate(endpoints):
   row=r[r.model.eq(model)&r.endpoint.eq(name)].iloc[0];assert abs(row.observed-actual[j])<1e-14;null=draws[model][:,j];assert abs(row.expected-null.mean())<1e-14
   assert abs(null.mean()-analytic[j])<=max(5*null.std()/math.sqrt(len(null)),1e-12)
   p=(1+sum(null>=actual[j]-1e-15))/(len(null)+1);assert abs(row.p_upper-p)<1e-14
   exactp=None
   if 'ESCC' in name:
    exact=exact_binary(b,nr,(model,name));count=int(b.loc[fg,name].sum());exactp=float(exact[count:].sum());assert abs(exact.sum()-1)<1e-12
    assert abs(p-exactp)<=max(5*math.sqrt(exactp*(1-exactp)/len(null)),2/(len(null)+1))
   comparisons.append(dict(model=model,endpoint=name,analytic_expectation=float(analytic[j]),stored_expected=row.expected,exact_binary_upper_p=exactp,stored_mc_p=row.p_upper))
 p=r.p_upper.to_numpy();expectedq=[]
 for value in p:
  possible=[v*len(p)/sum(p<=v) for v in p if v>=value];expectedq.append(min(1.,min(possible)))
 assert np.max(np.abs(r.q18-np.array(expectedq)))<1e-14
 assert r.q18.eq(1).all() and not r.positive_evidence_gate.any()
 pd.DataFrame(comparisons).to_csv(HERE/'independent_numerical_comparison.tsv',sep='\t',index=False,na_rep='NA',float_format='%.17g')
 reference=json.loads((HERE/'independent_reference_validation.json').read_text());assert reference['status']=='passed' and reference['full_sequences_compared']==342
 result={'status':'passed','inputs_sha_checks':len(inputs),'events_reclassified':len(events),'driver_choices_recomputed':len(events),'complete_sequence_motif_reconstructions':len(seq),'physical_motif_hit_checks':len(hits),'exact_binary_nulls':12,'planned_tests_recomputed':18,'coverage_background':751,'coverage_core':27,'eligible_background':67,'eligible_core':8,'independent_reference_receipt':'independent_reference_validation.json','skipped_checks':[]}
 (HERE/'validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
