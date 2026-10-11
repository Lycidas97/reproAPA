"""Apply frozen Figure 5 annotation and inclusion rules to the expanded pool.

This stage evaluates geometry and inclusion, before recognition-site scanning.
P--M uses common physical libraries; P--S retains the Figure 5 observed source.
"""
from pathlib import Path
import copy
import hashlib
import json
import importlib.util
import numpy as np
import pandas as pd
from prepare_union import HERE, F5, MOUSE, ROUTES

OUT = HERE/'annotation'; OUT.mkdir(exist_ok=True)
SOURCES = ['polyadb_main','polyasite_broad','gencode_polya']

def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m

legacy = module('migration', F5/'annotation-pas-global-v1/migrate.py')
choose = module('tracks', F5/'site-recurrence-v1/figure5-publication-v14r4/prepare_annotated_tracks.py').choose
paired = module('paired', HERE.parents[1]/'notes/2026-09-29-core-sensitivity/figure5_paired_pm.py')

def effects(model, contrast, cutoff):
    if contrast == 'MP':
        s = paired.paired_state_stats(model, cutoff)
        return s['delta_usage'], s['m_counts'], s['p_counts'], len(s['libraries']), len(s['libraries'])
    s = legacy.state_stats(model)
    a, b = s[cutoff,contrast[0]], s[cutoff,contrast[1]]
    return b[0]-a[0], a[2], b[2], a[1], b[1]

def main():
    refpath = F5/'annotation-pas-global-v1/human_reference_source_records.tsv'
    ref = pd.read_csv(refpath,sep='\t').drop_duplicates(['source_family','target_chrom','target_strand','target_cut_site','catalog_pas_id'])
    keys = ['source_family','target_chrom','target_strand','target_cut_site']
    assert not ref.duplicated(keys).any()
    lookup = ref.set_index(keys).catalog_pas_id.to_dict()
    ledger, summary, regions, missing, checks = [], [], [], [], 0
    hashes = {str(refpath): hashlib.sha256(refpath.read_bytes()).hexdigest()}
    for sp, routes in ROUTES.items():
        for route, study in routes:
            match = MOUSE/('pas-annotation-support-v1/studies' if sp=='mouse' else 'apa-majority-shared-v1/human_annotation')/study/'primary_pas_reference_matches.tsv.gz'
            path = HERE/'native'/sp/route/'models.json'
            for p in [path,match]: hashes[str(p)] = hashlib.sha256(p.read_bytes()).hexdigest()
            table = pd.read_csv(match,sep='\t',low_memory=False).set_index('atlas_pas_id')
            models = json.loads(path.read_text()); transformed = {}
            for gid, raw in models.items():
                selected = []
                for pas in raw['pas']:
                    row = table.loc[pas['atlas_pas_id']]
                    assert row.primary_gene_id == gid and row.cut_site == pas['cut_site']
                    candidates = []
                    for source in SOURCES:
                        distance = row[source+'_distance_nt']
                        if pd.isna(distance) or not 0 <= distance <= 50: continue
                        key = (source,pas['chrom'],pas['strand'],int(row[source+'_nearest_cut']))
                        if sp=='human' and key not in lookup:
                            missing.append(key); continue
                        identity = row[source+'_catalog_pas_id'] if sp=='mouse' else lookup[key]
                        cut = int(identity.split(':')[5])
                        candidates.append(dict(annotation_id=identity,annotation_cut=cut,source=source,distance_nt=int(distance)))
                    chosen = choose(candidates)
                    assert chosen is not None, (sp,route,gid,pas)
                    selected.append(chosen[0])
                    ledger.append(dict(species=sp,route=route,gene_id=gid,gene_symbol=raw['gene_symbol'],raw_pas=pas['atlas_pas_id'],raw_cut=pas['cut_site'],**chosen[0]))
                axis = sorted({(c['annotation_id'],c['annotation_cut']) for c in selected},key=lambda x:(x[1],x[0]),reverse=raw['family']['strand']=='-')
                membership = np.array([[int(c['annotation_id']==aid) for aid,cut in axis] for c in selected])
                assert (membership.sum(1)==1).all()
                model = copy.deepcopy(raw)
                for profile in model['profiles'].values():
                    for obj in [profile]+profile.get('libraries',[]):
                        if 'counts' not in obj: continue
                        old = np.asarray(obj['counts']); new = old @ membership
                        assert np.array_equal(old,old.astype(np.int64)) and old.sum()==new.sum() and old.sum()==obj['total']
                        obj['counts'] = new.tolist()
                        obj['usage'] = (new/new.sum()).tolist() if new.sum() else [0.]*len(new)
                        checks += 1
                model['pas'] = [dict(atlas_pas_id=aid,annotation_id=aid,cut_site=cut,member_cut_site_min=cut,member_cut_site_max=cut,chrom=raw['family']['chrom'],strand=raw['family']['strand']) for aid,cut in axis]
                model['family']['pas_ids'] = [a[0] for a in axis]
                model['native_eligible_contrasts'] = raw['eligible_contrasts']
                if len(axis)<2: model['eligible_contrasts'] = []
                transformed[gid] = model
                for contrast in ['SP','MP']:
                    r20 = effects(raw,contrast,20); a20 = effects(model,contrast,20); a50 = effects(model,contrast,50)
                    valid = np.isfinite(a20[0]).all() and np.isfinite(a50[0]).all()
                    base = dict(species=sp,route=route,gene_id=gid,gene_symbol=raw['gene_symbol'],contrast=contrast,native_pas=len(selected),annotated_pas=len(axis),eligible=contrast in model['eligible_contrasts'],n_ref20=a20[3],n_target20=a20[4],n_ref50=a50[3],n_target50=a50[4],native_MPRO20=float(np.ptp(r20[0])),annotated_MPRO20=float(np.ptp(a20[0])),annotated_MPRO50=float(np.ptp(a50[0])))
                    summary.append(base)
                    if not valid or not base['eligible']: continue
                    fam = model['family']; cuts=np.array([a[1] for a in axis])
                    boundaries=cuts if fam['strand']=='+' else cuts-1
                    edges=sorted({fam['common_start0'],fam['common_end0']}|{int(v) for v in boundaries if fam['common_start0']<v<fam['common_end0']})
                    for start,end in zip(edges[:-1],edges[1:]):
                        mask = legacy.mask_for(start,end,fam['strand'],cuts)
                        if mask.min()==mask.max(): continue
                        e20, e50 = float(a20[0]@mask),float(a50[0]@mask)
                        observed=min(float(c@k) for c in a50[1:3] for k in [mask,1-mask])>=5
                        regions.append(dict(**base,chrom=fam['chrom'],strand=fam['strand'],start0=start,end0=end,length=end-start,mask=json.dumps(mask.tolist()),effect20=e20,effect50=e50,observed_both5=observed,strict_change=observed and abs(e20)>=.1 and abs(e50)>=.1 and e20*e50>0))
            dest=OUT/sp/route; dest.mkdir(parents=True,exist_ok=True)
            (dest/'models.json').write_text(json.dumps(transformed))
            print(sp,route,len(transformed),flush=True)
    assert not missing, ('Incomplete human reference lookup',len(missing))
    for name,data in [('pas_mapping.tsv.gz',ledger),('model_effects.tsv',summary),('interval_effects.tsv.gz',regions)]:
        pd.DataFrame(data).to_csv(OUT/name,sep='\t',index=False)
    (OUT/'QA.json').write_text(json.dumps(dict(status='PASS',count_conservation_checks=checks,missing_reference_keys=len(missing),input_sha256=hashes,scope='Candidate geometry and interval inclusion; no new recognition-site calls'),indent=2)+'\n')

if __name__=='__main__': main()
