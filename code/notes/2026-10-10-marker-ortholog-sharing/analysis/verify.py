"""Validate hashes, fixed calibration and all archived aggregate margins."""
import hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];OLD=ROOT/'notes/2026-10-10-ortholog-threshold-scan'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:pd.read_csv(p,sep='\t',float_precision='round_trip')
checks=0
def check(value,why):
 global checks
 assert value,why;checks+=1

def main():
 m=json.loads((HERE/'analysis_manifest.json').read_text())
 for p,h in m['input_sha256'].items():check(sha(ROOT/p)==h,p)
 for n,h in m['output_sha256'].items():check(sha(HERE/n)==h,n)
 f=read(HERE/'ortholog_3x3_threshold_scan.tsv');members=read(HERE/'measured_ortholog_members.tsv.gz');ev=read(HERE/'study_events.tsv.gz')
 check(len(f)==972,'3 defs*2 contrasts*9 thresholds*2 models*9 cells');check(len(f[f.tested])==432,'3 duplicated family formal slots')
 check(ev.loc[ev.event.eq(-1),'effect_max'].le(-ev.loc[ev.event.eq(-1),'threshold']).all(),'negative every series threshold')
 check(ev.loc[ev.event.eq(1),'effect_min'].ge(ev.loc[ev.event.eq(1),'threshold']).all(),'positive every series threshold')
 check(ev.loc[ev.event.eq(-1),'can_negative'].all() and ev.loc[ev.event.eq(1),'can_positive'].all(),'direction capacity')
 old=read(OLD/'ortholog_3x3_threshold_scan.tsv');old=old[old.contrast.eq('P_S')].copy()
 keys=['contrast','threshold','model','human_direction','mouse_direction'];columns=[c for c in old if not c.startswith('q_') and c not in keys]
 old=old.set_index(keys).sort_index()
 for d,part in f[f.contrast.eq('P_S')].groupby('definition'):
  part=part.set_index(keys).sort_index();pd.testing.assert_frame_equal(old[columns],part[columns],check_exact=True);checks_before=1;check(True,'PS exact '+d)
 archive=np.load(HERE/'null_draws.npz');oldnull=np.load(OLD/'null_draws.npz')
 with archive,oldnull:
  check(len(archive.files)==72,'54 fresh PM +18 immutable PS contexts')
  for k in archive.files:
   a=archive[k]
   if k.split('|')[1]=='P_S':check(np.array_equal(a,oldnull[k]),'immutable PS draws '+k);continue
   _,co,t,model,d=k.split('|');t=float(t);part=members[members.definition.eq(d)&members.threshold.eq(t)]
   hc=part.human_direction.value_counts().reindex(['proximal','distal','nonrecurrent'],fill_value=0).to_numpy();mc=part.mouse_direction.value_counts().reindex(['proximal','distal','nonrecurrent'],fill_value=0).to_numpy()
   check(a.shape==(100000,3,3),'shape '+k);check(np.array_equal(a.sum(2),np.broadcast_to(hc,(100000,3))),'every human margin '+k);check(np.array_equal(a.sum(1),np.broadcast_to(mc,(100000,3))),'every mouse margin '+k)
   out=f[f.contrast.eq(co)&f.threshold.eq(t)&f.model.eq(model)&f.definition.eq(d)]
   for i,hd in enumerate(['proximal','distal','nonrecurrent']):
    for j,md in enumerate(['proximal','distal','nonrecurrent']):
     r=out[out.human_direction.eq(hd)&out.mouse_direction.eq(md)].iloc[0];v=a[:,i,j]
     check(r.observed==int(((part.human_direction==hd)&(part.mouse_direction==md)).sum()),'observed '+k)
     check(r.raw_p_enrichment==(1+int((v>=r.observed).sum()))/100001 if r.tested or np.isfinite(r.raw_p_enrichment) else True,'raw p '+k)
     check(abs(float(v.mean())-r.null_mean)<1e-12,'null mean '+k)
     if not r.inferentially_mobile and r.tested:check(r.p_for_correction==1 and pd.isna(r.p_enrichment),'NE kept slot '+k)
 core=read(HERE/'original_core27_flow.tsv');check(len(core)==729,'all original27 opportunities retained');check(core.groupby(['definition','threshold']).size().eq(27).all(),'27 every context')
 check(core.loc[~core.measurable,'human_direction'].isna().all() and core.loc[~core.measurable,'mouse_direction'].isna().all(),'coverage loss not neutral')
 result={'status':'SELF_NUMERIC_QA_PASS_AWAITING_INDEPENDENT_QA','checks_passed':checks,'checks_failed':0,'checks_skipped':0,'fixed_PS_draw_contexts':18,'fresh_PM_draw_contexts':54,'aggregate_draws_verified':7200000,'script_sha256':sha(Path(__file__))}
 (HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
