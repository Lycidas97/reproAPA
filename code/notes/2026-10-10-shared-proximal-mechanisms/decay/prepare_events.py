"""Classify native route drivers and freeze full terminal-UTR intervals before seed scans."""
from pathlib import Path
import gzip,hashlib,json,re
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
TABLES=ROOT/'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables'
CONFIG=ROOT/'notes/2026-10-06-dual-main-plot-data/provenance/polyapipe-pipeline.config.json'
CODES={'H1':'GSE118258','H2':'GSE150311','H3':'GSE178325','H4':'GSE221739','W':'GSE115943','C':'GSE103221'}

def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def save(d,name):
 d.to_csv(HERE/name,sep='\t',index=False,na_rep='NA',float_format='%.17g',compression={'method':'gzip','mtime':0} if name.endswith('.gz') else None)
def retained(cut,strand):return int(cut)-1 if strand=='+' else int(cut)
def contains(lo,hi,base):return int(lo)<=base<int(hi)

def parse_models(path,genes,sp):
 rows=[];gid_re=re.compile(r'gene_id "([^".]+)(?:\.\d+)?"');tid_re=re.compile(r'transcript_id "([^".]+)(?:\.\d+)?"')
 with gzip.open(path,'rt') as stream:
  for line in stream:
   if line.startswith('#'):continue
   c=line.rstrip().split('\t');
   if c[2] not in ('transcript','exon','CDS'):continue
   m=gid_re.search(c[8]);
   if m is None or m.group(1) not in genes:continue
   t=tid_re.search(c[8]);
   if t is None:continue
   rows.append(dict(species=sp,gene_id=m.group(1),transcript_id=t.group(1),chrom=c[0].removeprefix('chr'),strand=c[6],start0=int(c[3])-1,end0=int(c[4]),kind=c[2]))
 return pd.DataFrame(rows)

def classify(cut1,cut2,chrom,strand,utrs,models):
 p1,p2=retained(cut1,strand),retained(cut2,strand)
 u=utrs[utrs.chrom.astype(str).str.removeprefix('chr').eq(str(chrom).removeprefix('chr'))&utrs.strand.eq(strand)]
 common=u[u.utr_start.sub(1).le(min(p1,p2))&u.utr_end.gt(max(p1,p2))]
 if not common.empty:return 'strict_tandem_UTR',';'.join(sorted(common.transcript_id.unique()))
 m=models[models.chrom.eq(str(chrom).removeprefix('chr'))&models.strand.eq(strand)]
 if m.empty:return 'unknown',''
 ex=m[m.kind.eq('exon')];tx=m[m.kind.eq('transcript')]
 terms=[]
 for tid,g in ex.groupby('transcript_id',sort=True):
  e=g.sort_values(['start0','end0']).iloc[-1 if strand=='+' else 0];terms.append(e)
 term=pd.DataFrame(terms)
 if term.empty:return 'unknown',''
 same=term[term.start0.le(min(p1,p2))&term.end0.gt(max(p1,p2))]
 if not same.empty:return 'terminal_exon_nonUTR',';'.join(sorted(same.transcript_id.unique()))
 t1=term[term.start0.le(p1)&term.end0.gt(p1)];t2=term[term.start0.le(p2)&term.end0.gt(p2)]
 together=ex[ex.start0.le(min(p1,p2))&ex.end0.gt(max(p1,p2))]
 if not t1.empty and not t2.empty and together.empty:return 'alternative_terminal_exon_compatible',';'.join(sorted(set(t1.transcript_id)|set(t2.transcript_id)))
 unex=[p for p in (p1,p2) if ex[ex.start0.le(p)&ex.end0.gt(p)].empty and not tx[tx.start0.le(p)&tx.end0.gt(p)].empty]
 if unex:return 'intronic_compatible',''
 if not t1.empty or not t2.empty or not together.empty:return 'mixed',''
 return 'unknown',''

def main():
 rosterpath=ROOT/'notes/2026-10-10-rna-module-enrichment/ortholog_roster.tsv.gz'
 fgpath=ROOT/'notes/2026-10-10-rna-module-enrichment/target_members.tsv'
 roster=pd.read_csv(rosterpath,sep='\t');roster=roster[roster.contrast.eq('P_M')].copy();assert len(roster)==751 and roster.foreground_002.sum()==27
 fg=pd.read_csv(fgpath,sep='\t');assert set(fg[fg.contrast.eq('P_M')].human_gene_id)==set(roster[roster.foreground_002].human_gene_id)
 effpath=TABLES/'state_vectors/route_effects.primary.tsv.gz'
 e=pd.read_csv(effpath,sep='\t');e=e[e.setting.eq('supported50_p08')&e.contrast.eq('P_M')];e=e[(e.species.eq('human')&e.gene_id.isin(roster.human_gene_id))|(e.species.eq('mouse')&e.gene_id.isin(roster.mouse_gene_id))].copy()
 axpath=TABLES/'state_vectors/PAS_axes.primary.tsv.gz';axes=pd.read_csv(axpath,sep='\t',usecols=['species','study_id','setting','atlas_pas_id','primary_gene_id','chrom','strand','cut_site']);axes=axes[axes.setting.eq('supported50_p08')];lookup={tuple(r[k] for k in ['species','study_id','primary_gene_id','atlas_pas_id']):r for r in axes.to_dict('records')};assert len(lookup)==len(axes)
 eligpath=TABLES/'tandem/tandem_utr_eligibility.tsv.gz';elig=pd.read_csv(eligpath,sep='\t');elig_lookup={tuple(r[k] for k in ['species','study_id','gene_id']):r for r in elig.to_dict('records')}
 utrpath=TABLES/'tandem/annotation_terminal_utrs.tsv.gz';utrs=pd.read_csv(utrpath,sep='\t');utrby=dict(tuple(utrs.groupby(['species','gene_id'])))
 conf=json.loads(CONFIG.read_text());gtfs=[];models=[];inputs=[rosterpath,fgpath,effpath,axpath,eligpath,utrpath,CONFIG,HERE/'PROTOCOL.zh.md',HERE/'protocol.json']
 for sp in ('human','mouse'):
  path=ROOT/'notes/2026-10-03-polyapipe-rerun/statistics_package_v4r1/resources'/Path(conf['resources'][sp+'_ensembl_gtf']['path']).name
  assert sha(path)==conf['resources'][sp+'_ensembl_gtf']['sha256'];inputs.append(path);gtfs.append(path)
  models.append(parse_models(path,set(roster[sp+'_gene_id']),sp))
 model=pd.concat(models,ignore_index=True);save(model,'background_transcript_exon_CDS.tsv.gz');modelby=dict(tuple(model.groupby(['species','gene_id'])))
 empty_u=utrs.iloc[:0];empty_m=model.iloc[:0]
 events=[]
 for r in e.to_dict('records'):
  ids=json.loads(r['PAS_ids']);d=np.asarray(json.loads(r['delta_usage']),float);assert len(ids)==len(d)>=2
  assert abs(sum(d))<1e-8
  gi=min(range(len(ids)),key=lambda i:(-d[i],ids[i]));li=min(range(len(ids)),key=lambda i:(d[i],ids[i]))
  gain=lookup[(r['species'],r['study_id'],r['gene_id'],ids[gi])];loss=lookup[(r['species'],r['study_id'],r['gene_id'],ids[li])]
  assert gain['chrom']==loss['chrom'] and gain['strand']==loss['strand']
  cls,tids=classify(gain['cut_site'],loss['cut_site'],gain['chrom'],gain['strand'],utrby.get((r['species'],r['gene_id']),empty_u),modelby.get((r['species'],r['gene_id']),empty_m)) if d[gi]>0 and d[li]<0 else ('unknown','')
  eligrow=elig_lookup[(r['species'],r['study_id'],r['gene_id'])]
  iscore=bool(roster.set_index(r['species']+'_gene_id').loc[r['gene_id'],'foreground_002'])
  events.append(dict(species=r['species'],study_id=r['study_id'],route=r['route'],gene_id=r['gene_id'],gene_symbol=r['gene_symbol'],foreground_002=iscore,gain_pas=ids[gi],loss_pas=ids[li],gain_cut=int(gain['cut_site']),loss_cut=int(loss['cut_site']),chrom=gain['chrom'],strand=gain['strand'],gain_delta=float(d[gi]),loss_delta=float(d[li]),delta_DWUI=r['delta_DWUI'],gain_is_proximal=(gain['cut_site']<loss['cut_site']) if gain['strand']=='+' else (gain['cut_site']>loss['cut_site']),driver_class=cls,driver_compatible_transcripts=tids,PAS_number=r['PAS_number'],full_primary_axis_tandem=eligrow['tandem_utr_compatible'],full_primary_axis_compatible_transcripts=eligrow['compatible_transcripts']))
 ev=pd.DataFrame(events);assert not ev.duplicated(['species','study_id','route','gene_id']).any();save(ev,'route_event_inventory.tsv.gz')
 ev.to_csv(HERE.parent/'event_inventory.tsv.gz',sep='\t',index=False,na_rep='NA',float_format='%.17g',compression={'method':'gzip','mtime':0})
 (HERE.parent/'event_contract.json').write_text(json.dumps({'protocol_sha256':sha(HERE/'PROTOCOL.zh.md'),'inventory_sha256':sha(HERE.parent/'event_inventory.tsv.gz'),'unit':'species/study/route/gene','coordinate_rule':'native interbase cut; retained base0=cut-1 on plus, cut on minus','gene_classes_are_structural_compatibility':True,'all_751_ortholog_pairs_included':True,'rows':len(ev),'read_only_for_other_agents':True},indent=2)+'\n')
 genes=[];queries=[]
 for (sp,gene),g in e.groupby(['species','gene_id'],sort=True):
  units=sorted(g.study_id.unique());events_g=ev[ev.species.eq(sp)&ev.gene_id.eq(gene)]
  cls=sorted(events_g.driver_class.unique());status='pending';context=dict(species=sp,gene_id=gene,gene_symbol=g.gene_symbol.iloc[0],foreground_002=bool(events_g.foreground_002.iloc[0]),study_pattern=';'.join(units),route_pattern=';'.join(sorted(g.route.unique())),driver_classes=';'.join(cls),aggregate_driver_class=cls[0] if len(cls)==1 else 'mixed',all_drivers_strict=events_g.driver_class.eq('strict_tandem_UTR').all(),all_drivers_gain_proximal=events_g.gain_is_proximal.all(),n_routes=len(g),n_studies=len(units),PAS_number=float(g.PAS_number.median()),baseline=float(g.reference_DWUI.median()),mean_GEX=float(g.mean_GEX.median()),APA_coverage=float(g.APA_coverage.median()))
  z=[elig_lookup[(sp,u,gene)] for u in units];context['full_primary_axis_tandem']=all(r['tandem_utr_compatible'] for r in z)
  if not context['full_primary_axis_tandem']:status='non_tandem_full_primary_axis'
  elif not context['all_drivers_strict']:status='non_tandem_driver'
  else:
   common=set.intersection(*(set(str(r['compatible_transcripts']).split(';')) for r in z))
   if not common:status='no_common_transcript_across_studies'
   else:
    sites=[lookup[(sp,r['study_id'],gene,id_)] for r in g.to_dict('records') for id_ in json.loads(r['PAS_ids'])]
    if len(set(r['chrom'] for r in sites))!=1 or len(set(r['strand'] for r in sites))!=1:status='multichromosome_or_strand'
    else:
     strand=sites[0]['strand'];chrom=sites[0]['chrom'];ordered=sorted(set(int(r['cut_site']) for r in sites),reverse=strand=='-');pc,dc=ordered[0],ordered[-1];lo,hi=sorted((pc,dc));u=utrby.get((sp,gene),empty_u)
     fit=u[u.transcript_id.isin(common)&u.strand.eq(strand)&u.chrom.astype(str).str.removeprefix('chr').eq(str(chrom).removeprefix('chr'))&u.utr_start.sub(1).le(lo)&u.utr_end.ge(hi)].sort_values('transcript_id')
     if lo==hi:status='zero_length_interval'
     elif fit.empty:status='full_interval_outside_common_terminal_UTR'
     else:
      status='sequence_pending';context.update(chrom=chrom,strand=strand,proximal_cut=pc,distal_cut=dc,start0=lo,end0=hi,aUTR_length=hi-lo,chosen_transcript=fit.transcript_id.iloc[0],compatible_transcripts=';'.join(fit.transcript_id),native_PAS_union=';'.join(map(str,ordered)))
      queries.append(dict(species=sp,gene_id=gene,chrom=chrom,strand=strand,start0=lo,end0=hi))
  context['status']=status;genes.append(context)
 inventory=pd.DataFrame(genes);save(inventory,'gene_interval_inventory.before_sequence.tsv.gz');save(pd.DataFrame(queries),'sequence_queries.tsv')
 save(ev[ev.foreground_002],'core_route_events.tsv');save(inventory[inventory.foreground_002],'core_event_summary.tsv')
 (HERE/'input_hashes.json').write_text(json.dumps({str(p.relative_to(ROOT)):sha(p) for p in inputs},indent=2)+'\n')
 print('route events',len(ev),'core',int(ev.foreground_002.sum()));print(inventory.groupby(['species','status']).size().to_string())
 print('core classes',inventory[inventory.foreground_002].groupby(['species','aggregate_driver_class']).size().to_string())
if __name__=='__main__':main()
