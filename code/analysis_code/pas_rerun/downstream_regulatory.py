"""Rebuild terminal models, annotation identities and sequence sites for new PASs."""
import argparse
from collections import Counter,defaultdict
import copy
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess

import numpy as np
import pandas as pd
from scipy import sparse

from downstream_runtime import Context,Study

F5='analysis_code/release_2026-09-24/sources/F5/'
ROUTES={'mouse':[('GSE103221','CIP'),('GSE115943','WOT')],
        'human':[('GSE118258','GSE118258'),('GSE150311','GSE150311')]}
SOURCES=['polyadb_main','polyasite_broad','gencode_polya']


def retained_base(cut,strand):
    """Adapt interbase coordinates to the legacy structure functions explicitly."""
    return int(cut)+(strand=='-')


def numerical(ctx):
    structure=ctx.functions('notes/2026-10-01-mpro-biology/legacy_structure_consensus.py',
        ['stable_id','contains_cut','resolve_families','genomic_retention','classify_site'],{'defaultdict':defaultdict})
    migration=ctx.functions(F5+'annotation-pas-global-v1/migrate.py',['mask_for','state_stats'])
    paired=ctx.functions('notes/2026-09-29-core-sensitivity/figure5_paired_pm.py',['qualified_libraries','paired_state_stats'])
    return structure,migration,paired


def effect(model,contrast,threshold,migration,paired,common=True):
    if contrast=='MP' and common:
        p=paired.paired_state_stats(model,threshold)
        return p['delta_usage'],p['m_counts'],p['p_counts'],len(p['libraries']),len(p['libraries'])
    s=migration.state_stats(model);a,b=s[threshold,contrast[0]],s[threshold,contrast[1]]
    return b[0]-a[0],a[2],b[2],a[1],b[1]


def observed_regions(regions,common=True):
    """Match the observation gate to the estimator, retaining the raw audit columns."""
    regions=regions.copy()
    if common:
        paired=regions.contrast.eq('MP')
        regions.loc[paired,'observed_both5']=regions.loc[paired,'paired_observed_both5']
    else:
        for cutoff in [20,50]:regions[f'effect{cutoff}']=regions[f'separate_effect{cutoff}']
    regions['strict_change']=(regions.observed_both5&regions.effect20.abs().ge(.1)&
        regions.effect50.abs().ge(.1)&regions.effect20.mul(regions.effect50).gt(0))
    return regions


def dump(obj,path):
    def clean(v):
        if isinstance(v,dict):return {str(k):clean(x) for k,x in v.items()}
        if isinstance(v,(list,tuple)):return [clean(x) for x in v]
        if isinstance(v,np.ndarray):return clean(v.tolist())
        if isinstance(v,np.generic):return clean(v.item())
        if isinstance(v,float) and not np.isfinite(v):return None
        return v
    path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(clean(obj),allow_nan=False)+'\n')


def reference_lookup(ctx,sp):
    ref=ctx.read(ctx.resource('external_pas_source_records'),usecols=['species','source_family','target_chrom','target_strand','target_cut_site','catalog_pas_id'])
    key=['source_family','target_chrom','target_strand','target_cut_site']
    ref=ref[ref.species.eq(sp)].drop_duplicates(key+['catalog_pas_id'])
    assert not ref.duplicated(key).any(),'Original source cut maps to multiple catalog identities'
    return ref.set_index(key).catalog_pas_id.to_dict()


def strict_interval_pas(features,strict_pas):
    """Check QC on the primary axis; unavailable excluded PASs stay excluded."""
    keep=features.keep_for_apa_primary
    if keep.isna().any() or not keep.isin([True,False]).all():
        raise ValueError('Unknown/nonboolean PAS QC: keep_for_apa_primary')
    primary=features.loc[keep.astype(bool)].copy()
    for col in ['keep_for_apa_primary','transcript_model_conflict','multi_gene_ambiguous','internal_priming_flag']:
        if primary[col].isna().any() or not primary[col].isin([True,False]).all():
            raise ValueError('Unknown/nonboolean PAS QC: '+col)
        # A nullable object column retains object dtype after excluded rows leave.
        # Bitwise negation must act on booleans, not Python integer complements.
        primary[col]=primary[col].astype(bool)
    selected=strict_pas(primary)
    mask=pd.Series(False,index=features.index,dtype=bool)
    mask.loc[primary.index]=selected.astype(bool)
    return mask


def intervals(ctx):
    claims=ctx.read(ctx.analysis/'union_pool/figure5_union_claims.tsv')
    structure,migration,paired=numerical(ctx)
    prepare=ctx.functions('notes/2026-10-01-mpro-biology/legacy_prepare_pilot.py',['strict_pas','models_for_gene'])
    chooser=ctx.functions(F5+'site-recurrence-v1/figure5-publication-v14r4/prepare_annotated_tracks.py',['choose'])
    audit=[];mappings=[];regions=[];summaries=[];conservation=[]
    for sp,routes in ROUTES.items():
        sub=claims[claims.species.eq(sp)];wanted=set(sub.gene_id)
        utr=ctx.read(ctx.resource(sp+'_terminal_utrs'));catalog=ctx.read(ctx.resource(sp+'_catalog_utrs'))
        cg=dict(tuple(catalog[catalog.gene_id.isin(wanted)].groupby('gene_id')))
        models={gid:prepare.models_for_gene(g,cg.get(gid,catalog.iloc[:0])) for gid,g in utr[utr.gene_id.isin(wanted)].groupby('gene_id')}
        ref=reference_lookup(ctx,sp)
        for sid,route in routes:
            study=Study(ctx,sid);meta=study.meta[study.meta.route.eq(route)&study.meta.purity.ge(.8)]
            parts=[g for _,g in meta.groupby(['node','library_key']) if len(g)>=3 and g.n_cells.sum()>=30]
            meta=pd.concat(parts) if parts else meta.iloc[:0]
            f=study.f[strict_interval_pas(study.f,prepare.strict_pas)&study.f.primary_gene_id.isin(wanted)&study.f.matched_source_count_50nt.ge(1)].copy()
            fg=dict(tuple(f.groupby('primary_gene_id')));native={};annotated={};original=study.f.set_index('atlas_pas_id')
            for gid,gc in sub.groupby('gene_id',sort=True):
                gf=fg.get(gid,f.iloc[:0]).copy();gf['gene_id']=gid
                # Legacy functions use the retained nucleotide's one-based position.
                # All input/reference matching and exported genomic boundaries remain interbase.
                for col in ['member_cut_site_min','member_cut_site_max']:
                    if col not in gf:gf[col]=gf.cut_site
                for col in ['cut_site','member_cut_site_min','member_cut_site_max']:
                    gf[col]=gf[col]+gf.strand.eq('-').astype(int)
                families,inventory=structure.resolve_families(gf.to_dict('records'),models.get(gid,[]))
                active=[fam for fam in families if len(fam['pas_ids'])>=2]
                reason=('known_structure_review_Nopchap1' if gc.gene_symbol.iloc[0]=='Nopchap1' else
                    'no_supported_strict_pas' if gf.empty else 'no_terminal_model' if gid not in models else
                    'no_multi_pas_family' if not active else 'multiple_multi_pas_families' if len(active)>1 else '')
                if reason:
                    audit.extend(dict(species=sp,route=route,gene_id=gid,gene_symbol=gc.gene_symbol.iloc[0],contrast=c,status=reason) for c in gc.contrast)
                    continue
                fam=active[0];pas=sorted([p for p in inventory if p['family_id']==fam['family_id']],key=lambda p:p['member_cut_site_min'],reverse=fam['strand']=='-')
                indexes=[p['feature_index_1based']-1 for p in pas];gix=gf.feature_index_1based.to_numpy()-1;profiles={}
                family_counts=study.counts[:,indexes].toarray();gene_totals=np.asarray(study.counts[:,gix].sum(1)).ravel()
                for state,g in meta.groupby('node',sort=True):
                    libs=[]
                    for library,m in g.groupby('library_key',sort=True):
                        counts=family_counts[m.index].sum(0);total=int(counts.sum())
                        denom=int(gene_totals[m.index].sum())
                        libs.append(dict(library=library,total=total,coverage=total/denom if denom else 0.,counts=counts,
                            usage=counts/total if total else np.zeros(len(pas))))
                    counts=sum(l['counts'] for l in libs);total=int(counts.sum());denom=int(gene_totals[g.index].sum())
                    profiles[state]=dict(total=total,coverage=total/denom if denom else 0.,counts=counts,
                        usage=counts/total if total else np.zeros(len(pas)),libraries=libs)
                eligible=[]
                for contrast in gc.contrast:
                    status=('missing_state' if any(s not in profiles for s in contrast) else
                        'low_pooled_UMI' if min(profiles[s]['total'] for s in contrast)<50 else
                        'low_family_coverage' if min(profiles[s]['coverage'] for s in contrast)<.5 else 'ready')
                    if status=='ready':eligible.append(contrast)
                    audit.append(dict(species=sp,route=route,gene_id=gid,gene_symbol=gc.gene_symbol.iloc[0],contrast=contrast,status=status,
                        n_supported_pas=len(gf),n_multi_pas_families=len(active),inventory_statuses=json.dumps(dict(Counter(p['family_status'] for p in inventory)))))
                raw=dict(species=sp,route=route,gene_id=gid,gene_symbol=gc.gene_symbol.iloc[0],family=fam,pas=pas,profiles=profiles,eligible_contrasts=eligible,
                    model_cut_coordinate='1-based-retained-nucleotide; minus = interbase + 1')
                native[gid]=raw;selected=[]
                for pasrow in pas:
                    row=original.loc[pasrow['atlas_pas_id']];candidates=[]
                    for source in SOURCES:
                        distance=row[source+'_distance_nt']
                        if pd.isna(distance) or not 0<=distance<=50:continue
                        key=(source,row.chrom,row.strand,int(row[source+'_nearest_cut']))
                        identity=ref[key];bits=identity.split(':');assert bits[1]==sp and bits[3:5]==[row.chrom,row.strand]
                        candidates.append(dict(annotation_id=identity,annotation_cut=int(bits[5]),source=source,distance_nt=int(distance)))
                    chosen=chooser.choose(candidates);assert chosen is not None
                    selected.append(chosen[0]);mappings.append(dict(species=sp,route=route,gene_id=gid,gene_symbol=raw['gene_symbol'],
                        native_pas_id=row.name,native_cut_interbase=row.cut_site,**chosen[0],annotation_coordinate='0-based-interbase',
                        sources=';'.join(sorted({c['source'] for c in chosen}))))
                axis=sorted({(s['annotation_id'],s['annotation_cut']) for s in selected},key=lambda x:(x[1],x[0]),reverse=fam['strand']=='-')
                membership=np.array([[int(s['annotation_id']==a) for a,c in axis] for s in selected]);assert (membership.sum(1)==1).all()
                model=copy.deepcopy(raw)
                for state,profile in model['profiles'].items():
                    for obj in [profile]+profile['libraries']:
                        old=np.asarray(obj['counts']);new=old@membership;assert old.sum()==new.sum()==obj['total']
                        obj['counts']=new;obj['usage']=new/new.sum() if new.sum() else np.zeros(len(axis))
                        conservation.append(dict(species=sp,route=route,gene_id=gid,state=state,library=obj.get('library','pooled'),before=int(old.sum()),after=int(new.sum())))
                model['pas']=[dict(atlas_pas_id=a,annotation_id=a,cut_interbase=c,cut_site=retained_base(c,fam['strand']),
                    member_cut_site_min=retained_base(c,fam['strand']),member_cut_site_max=retained_base(c,fam['strand']),chrom=fam['chrom'],strand=fam['strand']) for a,c in axis]
                model['family']['pas_ids']=[a for a,c in axis];model['native_eligible_contrasts']=eligible
                if len(axis)<2:model['eligible_contrasts']=[]
                annotated[gid]=model;cuts=np.array([p['cut_site'] for p in model['pas']]);bounds=[c for a,c in axis]
                edges=sorted({fam['common_start0'],fam['common_end0']}|{c for c in bounds if fam['common_start0']<c<fam['common_end0']})
                for contrast in sorted(set(eligible)):
                    values={t:effect(model,contrast,t,migration,paired) for t in [20,50]}
                    separate={t:effect(model,contrast,t,migration,paired,False) for t in [20,50]}
                    summaries.append(dict(species=sp,route=route,gene_id=gid,contrast=contrast,native_pas=len(pas),annotated_pas=len(axis),
                        eligible=contrast in model['eligible_contrasts'],n_ref20=values[20][3],n_target20=values[20][4],n_ref50=values[50][3],n_target50=values[50][4]))
                    if len(axis)<2:continue
                    for start,end in zip(edges[:-1],edges[1:]):
                        mask=migration.mask_for(start,end,fam['strand'],cuts)
                        if mask.min()==mask.max():continue
                        support=lambda stats:all(min(float(c@mask),float(c@(1-mask)))>=5 for c in stats[50][1:3])
                        e={t:float(values[t][0]@mask) for t in [20,50]};observed=support(separate)
                        regions.append(dict(species=sp,route=route,gene_id=gid,gene_symbol=raw['gene_symbol'],contrast=contrast,
                            chrom=fam['chrom'],strand=fam['strand'],start0=start,end0=end,length=end-start,
                            region_id=f'{route}|{gid}|annotation:{start}-{end}',mask=json.dumps(mask.tolist()),eligible=True,effect20=e[20],effect50=e[50],
                            separate_effect20=float(separate[20][0]@mask),separate_effect50=float(separate[50][0]@mask),
                            observed_both5=observed,paired_observed_both5=support(values),
                            strict_change=observed and abs(e[20])>=.1 and abs(e[50])>=.1 and e[20]*e[50]>0))
            dump(native,ctx.out/'native'/sp/route/'models.json');dump(annotated,ctx.out/'annotation'/sp/route/'models.json')
            ctx.table(meta,f'native/{sp}/{route}/selected_metacells.tsv.gz')
            print(sp,route,len(native),len(annotated),flush=True)
    ctx.table(audit,'eligibility.tsv',columns=['species','route','gene_id','contrast','status'])
    ctx.table(mappings,'pas_mapping.tsv.gz',columns=['species','route','gene_id','native_pas_id','native_cut_interbase','annotation_id','annotation_cut'])
    ctx.table(conservation,'count_conservation.tsv.gz',columns=['species','route','gene_id','state','library','before','after'])
    ctx.table(summaries,'model_effects.tsv',columns=['species','route','gene_id','contrast','native_pas','annotated_pas','eligible'])
    ctx.table(regions,'all_coordinate_regions.tsv.gz',columns=['species','route','gene_id','contrast','chrom','strand','start0','end0','length','region_id','mask','effect20','effect50','observed_both5'])
    ctx.finish(dict(new_native_counts_and_candidates=True,all_GEX_qualified_group_counts=True,source_distance_identity_mapping=True,
        interbase_minus_boundary_converted=True,many_to_one_count_conservation=True,one_PAS_collapses_excluded=True,
       _PM_common_libraries=True,marginal_and_paired_support_separate=True))


def scan(ctx):
    structure,migration,paired=numerical(ctx)
    core=ctx.functions(F5+'full-cohort-v2/sources/execution-scripts/core.py',['ambiguous_overlaps','masks','change'])
    seeds=ctx.functions('notes/2026-10-01-figure5-union-update/resources/mouse_time_aware_integration_pilot_v1/regulatory_catalog_v1/scripts/build_catalog.py',
        ['rc','canonical_patterns','seed_hits'],{'defaultdict':defaultdict})
    bundle=json.loads(ctx.resource('regulatory_bundles').read_text());site_schema=None;coverage=[];scan_audit=[]
    for sp,routes in ROUTES.items():
        native={r:json.loads(ctx.register(ctx.analysis/f'intervals/native/{sp}/{r}/models.json').read_text()) for s,r in routes}
        wanted={m['region_id'] for gs in native.values() for g in gs.values() if g['eligible_contrasts'] for m in g['family']['models']}
        catalog=ctx.read(ctx.resource(sp+'_catalog_utrs')).set_index('region_id');regions=catalog.loc[sorted(wanted)].reset_index().to_dict('records')
        regionmap={r['region_id']:r for r in regions};seqs={};rid=None;pieces=[]
        with gzip.open(ctx.resource(sp+'_utr_fasta'),'rt') as handle:
            for line in handle:
                if line.startswith('>'):
                    if rid in wanted:seqs[rid]=''.join(pieces).upper()
                    rid=line[1:].split()[0];pieces=[]
                elif rid in wanted:pieces.append(line.strip())
            if rid in wanted:seqs[rid]=''.join(pieces).upper()
        assert set(seqs)==wanted
        out=ctx.out/sp;out.mkdir();fa=out/'query.fa'
        with fa.open('w') as handle:
            for r in regions:
                seq=seqs[r['region_id']];assert len(seq)==r['end0']-r['start0']
                handle.write('>'+r['region_id']+'\n'+seq+'\n')
        b=next(v for k,v in bundle.items() if '/'+sp+'_time_aware' in k)
        ann=ctx.read(ctx.resource(sp+'_mirna_annotations')).fillna('')
        mirnames={s:sorted({n for v in g.mirbase_exact_aliases for n in v.split(';') if n}) for s,g in ann.groupby('seed_dna')}
        provenance=ctx.resource(sp+'_motif_provenance')
        prov=pd.DataFrame(json.loads(provenance.read_text())) if provenance.suffix=='.json' else ctx.read(provenance)
        direct={mid:sorted(set(g.RBP_Name)) for mid,g in prov[prov.scan_status.eq('scan_ready')].groupby('Motif_ID')}
        transfer=defaultdict(set)
        for r in b.get('transfer_provenance',[]):transfer[r['Motif_ID']].add(r['mouse_symbol'])
        cached={rid:[] for rid in wanted};excluded=Counter()
        def add(rid,start,end,category,feature,kind):
            r=regionmap[rid];names=mirnames.get(feature,[]) if category=='miRNA' else direct.get(feature,[]) if category=='RBP' else sorted(transfer[feature])
            if not names:excluded['unmapped_names']+=1;return
            gs,ge=(r['start0']+start,r['start0']+end) if r['strand']=='+' else (r['end0']-end,r['end0']-start)
            cached[rid].append(dict(chrom=r['chrom'],strand=r['strand'],genomic_start0=gs,genomic_end0=ge,category=category,feature=feature,names=names,site_type=kind))
        fimo=ctx.resource('fimo');assert subprocess.check_output([str(fimo),'--version'],text=True).strip()=='5.5.2'
        for category,key in [('RBP',sp+'_direct_meme')]+([('RBP_transfer','mouse_transfer_meme')] if sp=='mouse' else []):
            dest=out/(category+'.fimo.tsv');argv=[str(fimo),'--text','--norc','--thresh','0.0001','--bfile',str(ctx.resource(sp+'_motif_background')),str(ctx.resource(key)),str(fa)]
            if wanted:
                with dest.open('w') as handle,(out/(category+'.stderr.log')).open('w') as log:subprocess.run(argv,stdout=handle,stderr=log,check=True)
                with dest.open() as handle:
                    for row in csv.DictReader(handle,delimiter='\t'):
                        if not row.get('motif_id') or row['motif_id'].startswith('#'):continue
                        start,end=int(row['start'])-1,int(row['stop']);rid=row['sequence_name']
                        assert row['strand']=='+' and seqs[rid][start:end]==row['matched_sequence'].upper()
                        add(rid,start,end,category,row['motif_id'],'PWM')
            else:dest.write_text('motif_id\tsequence_name\tstart\tstop\n')
            scan_audit.append(dict(species=sp,category=category,regions=len(wanted),argv=json.dumps(argv)))
        patterns=seeds.canonical_patterns(b['mirna_seeds']);matcher=re.compile('(?=('+'|'.join(sorted(patterns,key=lambda x:(-len(x),x)))+'))')
        for rid,seq in seqs.items():
            for start,length,seed,kind in seeds.seed_hits(seq,patterns,matcher):add(rid,start,start+length,'miRNA',seed,kind)
        rows=[];native_rows=[]
        for sid,route in routes:
            annotated=json.loads(ctx.register(ctx.analysis/f'intervals/annotation/{sp}/{route}/models.json').read_text())
            for gid,g in native[route].items():
                model=annotated[gid]
                if not g['eligible_contrasts']:continue
                hits={};spans=defaultdict(set)
                for m in g['family']['models']:
                    for h in cached[m['region_id']]:
                        if not structure.classify_site(h,g['family'],g['pas'])['site_status'].startswith('model_consensus'):
                            excluded['outside_consensus']+=1;continue
                        key=(h['category'],h['genomic_start0'],h['genomic_end0']);factor=h['feature'] if h['category']=='miRNA' else ';'.join(h['names'])
                        obj=hits.setdefault(key,dict(factors={},types=set()));obj['factors'][factor]=h['names'];obj['types'].add(h['site_type']);spans[key[0],factor].add(key[1:])
                ambiguous={k:core.ambiguous_overlaps(v) for k,v in spans.items()};nsite=Counter();fam=model['family']
                effects={(coordinate,c,t,common):effect(current,c,t,migration,paired,common)
                    for coordinate,current in [('native',g),('annotation',model)]
                    for c in current['eligible_contrasts'] for t in [20,50] for common in [False,True]}
                for (category,start,end),h in hits.items():
                    factors={k:v for k,v in h['factors'].items() if (start,end) not in ambiguous[category,k]}
                    if not factors:excluded['all_factor_links_overlap_ambiguous']+=1;continue
                    nsite[category]+=1
                    base=dict(species=sp,route=route,gene_id=gid,gene_symbol=g['gene_symbol'],family_id=fam['family_id'],category=category,chrom=fam['chrom'],strand=fam['strand'],start0=start,end0=end,
                        names=';'.join(sorted({n for v in factors.values() for n in v})),feature_groups=';'.join(sorted(factors)),site_types=';'.join(sorted(h['types'])))
                    for coordinate,current,destination in [('native',g,native_rows),('annotation',model,rows)]:
                        mask=core.masks(start,end,fam['strand'],current['pas']);exact=np.array_equal(mask[:,0],mask[:,1])
                        base_coordinate=dict(base,coordinate=coordinate,mask_json=json.dumps(mask.astype(int).tolist(),separators=(',',':')),
                            exact=exact,variable=bool(mask.min()!=mask.max()))
                        for co in current['eligible_contrasts']:
                            pa,pb=current['profiles'][co[0]],current['profiles'][co[1]]
                            pooled=core.change(pb['usage'],pa['usage'],mask)
                            row=dict(base_coordinate,contrast=co,pooled_lo=pooled[0],pooled_hi=pooled[1],reference_UMI=pa['total'],target_UMI=pb['total'])
                            for t in [20,50]:
                                for common,prefix in [(True,''),(False,'separate_')]:
                                    v=effects[coordinate,co,t,common];value=core.change(v[0],np.zeros(len(v[0])),mask)
                                    row.update({f'{prefix}eq{t}_lo':value[0],f'{prefix}eq{t}_hi':value[1],f'{prefix}eq{t}_n_reference':v[3],f'{prefix}eq{t}_n_target':v[4]})
                            destination.append(row)
                coverage.append(dict(species=sp,route=route,gene_id=gid,gene_symbol=g['gene_symbol'],**nsite))
        schema=['species','route','gene_id','gene_symbol','family_id','category','chrom','strand','start0','end0','names','feature_groups','site_types','mask_json','exact','variable','contrast','pooled_lo','pooled_hi','reference_UMI','target_UMI']+[f'{p}eq{t}_{x}' for p in ['','separate_'] for t in [20,50] for x in ['lo','hi','n_reference','n_target']]
        ctx.table(rows,f'full-cohort-v2/data/{sp}/site_effects.tsv.gz',columns=schema)
        ctx.table(native_rows,f'native/full-cohort-v2/data/{sp}/site_effects.tsv.gz',columns=schema)
        scan_audit.append(dict(species=sp,category='all',regions=len(wanted),sites=len(rows),exclusions=json.dumps(dict(excluded))))
    ctx.table(coverage,'site_coverage.tsv',columns=['species','route','gene_id','gene_symbol'])
    ctx.table(scan_audit,'scan_audit.tsv')
    ctx.finish(dict(new_union_regions_scanned=True,direct_and_transfer_rescanned=True,forward_FIMO_5_5_2_p1e4=True,
        exact_physical_sites_deduplicated=True,all_contrasts_recomputed=True,separate_and_paired_estimands_saved=True,
        native_and_annotation_coordinate_masks_recomputed=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);p.add_argument('--node',choices=['intervals','site_scan'],required=True)
    a=p.parse_args();{'intervals':intervals,'site_scan':scan}[a.node](Context(a.config,a.node,a.out))
