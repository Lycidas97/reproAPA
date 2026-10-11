"""Independent post-seal site vector/correspondence/support audit; no site imports."""
from pathlib import Path
from collections import defaultdict
import hashlib,json
import numpy as np
import pandas as pd

B=Path('notes/2026-10-10-conserved-response');S=B/'sites';O=B/'review';T=Path('notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables')
CHECKS=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check(name,ok):
    CHECKS.append({'check':name,'pass':bool(ok)})
    if not ok:raise AssertionError(name)
def read(p,**kw):return pd.read_csv(p,sep='\t',float_precision='round_trip',low_memory=False,**kw)
def merged(ids,usage,mapping):
    out=defaultdict(float)
    for i,v in zip(ids,usage):out[mapping.get(i,'UNMAPPED:'+i)]+=v
    return dict(out)
def close_dict(a,b):return set(a)==set(b) and all(abs(a[k]-b[k])<2e-12 for k in a)

def run():
    # Parent must announce the site delivery sealed before this file is executed.
    manifest=S/'event_sequence_sealed_manifest.json';check('site manifest exists; no provisional results read',manifest.exists())
    check('announced event/sequence seal sha',sha(manifest)=='c7a4de318f0a4e500d1738787f21ca240031c401270d05c94696e7e76abc3645')
    for name,digest in json.loads(manifest.read_text())['files'].items():check('sealed event/sequence output '+name,sha(S/name)==digest)
    contexts=read(S/'site_event_contexts.tsv.gz');contexts=contexts[contexts.source_population.isin(['main48','hN79'])&contexts.gene_threshold.eq(.02)].copy()
    traces=read(S/'site_usage_traces.tsv.gz');traces=traces[traces.source_population.isin(['main48','hN79'])&traces.gene_threshold.eq(.02)].copy()
    routes=read(T/'endpoints/route_effects.tsv.gz');routes=routes[routes.setting.eq('supported50_p08')&((routes.contrast.eq('P_S')&routes.design.eq('true_source'))|(routes.contrast.eq('P_M')&routes.design.eq('same_library')))];routes=routes.set_index(['species','study_id','gene_id','contrast'],drop=False).sort_index()
    naive=read(T/'naive_state/library_pair_effects.tsv.gz');naive=naive[naive.version.eq('v2_primary')&naive.setting.eq('supported50_p08')&naive.node_a.eq('N')&naive.node_b.eq('P')&naive.library_a.str.endswith(('6','8'))&naive.library_b.str.endswith('7')]
    mapping=read(T/'cis/detected_to_annotation.tsv.gz');mapping=mapping.set_index('native_pas_id').id.to_dict()
    catalog=read(T/'cis/annotation_catalog.tsv.gz').set_index('id');assigned=read(T/'cis/reciprocal_one_to_one.tsv.gz');candidates=read(T/'cis/reciprocal_PAS_candidates.tsv.gz')
    # Replay assignment using an independent Python tuple ordering across all candidates.
    items=sorted(candidates.to_dict('records'),key=lambda x:(not x['reciprocal25'],not x['reciprocal50'],x['forward_distance'],x['backward_distance'],x['human_id'],x['mouse_id']))
    hh=set();mm=set();own=[]
    for x in items:
        if x['human_id'] in hh or x['mouse_id'] in mm:continue
        hh.add(x['human_id']);mm.add(x['mouse_id']);own.append((x['human_id'],x['mouse_id']))
    check('global frozen reciprocal greedy assignment',own==list(zip(assigned.human_id,assigned.mouse_id)))
    check('assigned unique both sides',assigned.human_id.is_unique and assigned.mouse_id.is_unique)
    # Feature metadata only; counts arrays are not accessed.
    npz=Path('notes/2026-10-10-marker-reference-sensitivity/inputs.npz')
    check('actual primary feature source sha',sha(npz)=='6652eb94bc70ad1e361992d2991a4e24df9c7d0fe83e694db155ce156bbaccb6')
    features={};cuts={};genes={}
    with np.load(npz,allow_pickle=False) as archive:
        for study in sorted(set(traces.study_id)):
            values=json.loads(str(archive[study+'_features_json']))
            keep=[x for x in values if str(x.get('keep_for_apa_primary')).lower()=='true' and x.get('primary_gene_id')]
            features[study]=keep;cuts[study]={x['atlas_pas_id']:x for x in keep}
            g=defaultdict(list)
            for x in keep:g[x['primary_gene_id']].append(x)
            genes[study]=g
    terminals=read(T/'tandem/annotation_terminal_utrs.tsv.gz')
    profiles={};max_vector_error=max_native_error=0.
    # Every primary context trace is replayed against actual original endpoint vectors.
    for row in traces.itertuples(index=False):
        gid=row.human_gene_id if row.species=='human' else row.mouse_gene_id
        if row.source_population=='hN79' and row.species=='human':
            r=naive[naive.gene_id.eq(gid)&naive.library_a.eq(row.series)];check('actual hN source pair unique',len(r)==1);r=r.iloc[0]
            ids=json.loads(r.PAS_ids);target=json.loads(r.usage_a);reference=json.loads(r.usage_b);raw_delta=r.delta_DWUI
        else:
            c='P_M' if row.source_population=='hN79' else row.contrast
            r=routes.loc[(row.species,row.study_id,gid,c)]
            if isinstance(r,pd.Series):r=r.to_frame().T
            if row.study_id=='GSE178325':r=r[r.route.str.endswith(row.series.removeprefix('178325_'))]
            elif row.study_id in ['GSE103221','GSE115943']:r=r[r.route.eq(row.series)]
            check('actual P-S source/P-M same-library source unique',len(r)==1);r=r.iloc[0]
            ids=json.loads(r.PAS_ids);target=json.loads(r.target_usage);reference=json.loads(r.reference_usage);raw_delta=r.delta_DWUI
        check('native vectors finite and sum1',np.isfinite(target).all() and np.isfinite(reference).all() and abs(sum(target)-1)<2e-12 and abs(sum(reference)-1)<2e-12)
        a,b=merged(ids,target,mapping),merged(ids,reference,mapping);ta,tb=json.loads(row.target),json.loads(row.reference)
        check('actual native annotation-merged target',close_dict(a,ta));check('actual native annotation-merged reference',close_dict(b,tb))
        max_vector_error=max(max_vector_error,max(abs(a[k]-ta[k]) for k in a),max(abs(b[k]-tb[k]) for k in b))
        check('original delta source exact',abs(row.original_delta-raw_delta)<2e-12)
        f=genes[row.study_id][gid];span=max(x['cut_site'] for x in f)-min(x['cut_site'] for x in f)
        check('full primary native span positive',span>0)
        weight=np.asarray([cuts[row.study_id][i]['cut_site']*(1 if cuts[row.study_id][i]['strand']=='+' else -1)/span for i in ids])
        replay=float((np.asarray(target)-np.asarray(reference))@weight);max_native_error=max(max_native_error,abs(replay-raw_delta));check('native original DWUI from actual cuts/full axis',abs(replay-raw_delta)<1e-9)
        key=(row.source_population,row.human_gene_id,row.contrast,row.tolerance,row.species)
        profiles.setdefault(key,[]).append({'study':row.study_id,'series':row.series,'target':a,'reference':b})
    reciprocal25=reciprocal50=0
    for row in contexts.to_dict('records'):
        hg,mg=row['human_gene_id'],row['mouse_gene_id'];tol=int(row['tolerance']);sp_profiles={};deltas={};driver={};strict={};supports={}
        for sp,gid in [('human',hg),('mouse',mg)]:
            pp=profiles[(row['source_population'],hg,row['contrast'],tol,sp)];sp_profiles[sp]=pp
            keys=sorted(set().union(*(set(x['target'])|set(x['reference']) for x in pp)));studies=defaultdict(list)
            for p in pp:studies[p['study']].append({k:p['target'].get(k,0)-p['reference'].get(k,0) for k in keys})
            dm={s:{k:sum(x[k] for x in values)/len(values) for k in keys} for s,values in studies.items()}
            delta={k:sum(x[k] for x in dm.values())/len(dm) for k in keys};gain=min(keys,key=lambda k:(-delta[k],k));loss=min(keys,key=lambda k:(delta[k],k));deltas[sp]=delta;driver[sp]=(gain,loss)
            check('complete-vector maximum gain/loss lexical ties',gain==row[sp+'_gain_id'] and loss==row[sp+'_loss_id'])
            check('equal-study mean driver deltas',abs(delta[gain]-row[sp+'_gain_delta'])<2e-12 and abs(delta[loss]-row[sp+'_loss_delta'])<2e-12)
            all_count=sum(all(x[gain]>=.02 and x[loss]<=-.02 for x in vals) for vals in studies.values());mean_count=sum(x[gain]>=.02 and x[loss]<=-.02 for x in dm.values());supports[sp]=all_count
            check('primary all-series study support',all_count==row[sp+'_gain_support_studies']);check('separate mean-study supplemental support',mean_count==row[sp+'_gain_support_studies_mean'])
            common=None
            for study in studies:
                f=genes[study][gid];chrset={str(x['chrom']).removeprefix('chr') for x in f};strandset={x['strand'] for x in f};valid=set()
                if len(f)>=2 and len(chrset)==1 and len(strandset)==1:
                    lo=min(x['cut_site']+(x['strand']=='-') for x in f);hi=max(x['cut_site']+(x['strand']=='-') for x in f)
                    tu=terminals[terminals.species.eq(sp)&terminals.gene_id.eq(gid)]
                    for u in tu.itertuples():
                        if str(u.chrom).removeprefix('chr') in chrset and u.strand in strandset and u.utr_start<=lo and u.utr_end>=hi:valid.add(u.transcript_id)
                common=valid if common is None else common&valid
            strict[sp]=bool(common)
            check('full actual primary axes common terminal UTR',strict[sp]==row[sp+'_strict_UTR'])
            if strict[sp]:check('all common transcripts literal',set(row[sp+'_common_transcripts'].split(';'))==common)
        pairs=assigned[assigned.human_gene_id.eq(hg)&assigned['reciprocal'+str(tol)]&assigned.human_id.isin(deltas['human'])&assigned.mouse_id.isin(deltas['mouse'])].copy()
        hpairs={r.human_id:r.mouse_id for r in pairs.itertuples()};mpairs={m:h for h,m in hpairs.items()}
        check('all mapped PAS count',len(pairs)==row['matched_PAS'])
        ordered=sorted(hpairs,key=lambda k:catalog.loc[k,'cut_site']*(1 if catalog.loc[k,'strand']=='+' else -1));hpos=[catalog.loc[k,'cut_site']*(1 if catalog.loc[k,'strand']=='+' else -1) for k in ordered];mpos=[catalog.loc[hpairs[k],'cut_site']*(1 if catalog.loc[hpairs[k],'strand']=='+' else -1) for k in ordered]
        order=len(pairs)>=2 and all(a<b for a,b in zip(hpos,hpos[1:])) and all(a<b for a,b in zip(mpos,mpos[1:]));check('reciprocal transcript order',order==row['order_eligible'])
        masses={};thresholds={};partners={}
        for sp,mp in [('human',hpairs),('mouse',mpairs)]:
            gain,loss=driver[sp];partners[sp]=gain in mp and loss in mp and gain!=loss;thresholds[sp]=gain!=loss and deltas[sp][gain]>=.02 and deltas[sp][loss]<=-.02
            mass=min(sum(p[node].get(k,0) for k in mp) for p in sp_profiles[sp] for node in ['target','reference']);masses[sp]=mass
            check('each actual study/state matched mass',abs(mass-row[sp+'_min_retained_mass'])<2e-12)
            check('full maximum driver partners/threshold',partners[sp]==row[sp+'_drivers_have_partners'] and thresholds[sp]==row[sp+'_driver_threshold'])
        massok=order and min(masses.values())>=.8;judge=all(partners.values()) and all(thresholds.values());hgain,hloss=driver['human'];mgain,mloss=driver['mouse'];double=judge and hpairs.get(hgain)==mgain and hpairs.get(hloss)==mloss;strictdouble=all(strict.values()) and massok and double
        check('maximum gain-to-gain loss-to-loss correspondence',double==row['double_driver_correspondence'])
        check('primary all-series two-independent-study flag',row['two_study_exact_event_both']==(strictdouble and supports['human']>=2 and supports['mouse']>=2))
        # Literal scope in the addendum: parent case only shared contrast; parent control only recurrent side.
        expected=[]
        if row['parent_case'] and row['context_role']=='shared_recurrent' and strictdouble:expected=[(hgain,mgain),(hloss,mloss)]
        if not row['parent_case'] and row['context_role']=='single_species_recurrent' and all(strict.values()) and massok:
            side=row['recurrent_side'];gain,loss=driver[side]
            if thresholds[side] and partners[side]:expected=[(gain,hpairs[gain]),(loss,hpairs[loss])] if side=='human' else [(mpairs[gain],gain),(mpairs[loss],loss)]
        actual=[tuple(v) for v in json.loads(row['sequence_driver_pairs'])];check('case/control actual contrast and recurrent-side driver scope',set(expected)==set(actual) and bool(expected)==row['sequence_candidate'])
        if tol==25:reciprocal25+=int(strictdouble)
        else:reciprocal50+=int(strictdouble)
    result={'status':'PASS','scope':'all primary main48/hN79 .02 contexts, all actual vectors, full gene axes, reciprocal drivers and both support definitions; sequence p/q audit delegated to root','primary_context_rows':len(contexts),'primary_trace_rows':len(traces),'global_candidate_pairs':len(candidates),'global_assigned_pairs':len(assigned),'max_original_native_DWUI_replay_error':max_native_error,'max_annotation_vector_error':max_vector_error,'strict_double_contexts_25':reciprocal25,'strict_double_contexts_50':reciprocal50,'checks':CHECKS,'site_manifest_sha256':sha(manifest),'studies_vote_clarification_sha256':sha(B/'SITE-STUDY-VOTE-CLARIFICATION.zh.md'),'site_contrast_addendum_sha256':sha(B/'SITE-CONTRAST-ADDENDUM.zh.md'),'actual_feature_npz_sha256':sha(npz),'audited_output_sha256':{p.name:sha(p) for p in [S/'site_event_contexts.tsv.gz',S/'site_usage_traces.tsv.gz']}}
    (O/'site_validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['checks','audited_output_sha256']}))

if __name__=='__main__':run()
