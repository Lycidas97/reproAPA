"""Six UB/read-end sensitivity rules on new native polyApipe windows and axes.

The primary analysis retains native polyApipe UR counts. These BAM measurements
use the historical CB/UB/read-end sensitivity definition; they are not asserted
to reproduce native caller counts or its read filters.
"""
import argparse
import bisect
from collections import Counter
from concurrent.futures import ThreadPoolExecutor,as_completed
import json
from pathlib import Path
import re
import time

import numpy as np
import pandas as pd
from scipy import sparse

from downstream_runtime import Context,Study
from polyapipe_prepare import digest,save
from polyapipe_effects import weights

MODES=['original','atlas_dedup','transferred_observed','transferred','primary_clipped','all_clipped']
SOURCE='notes/2026-09-26-endpoint-counting/common.py'


def canon(chrom):return 'chrM' if str(chrom) in ['MT','M','chrM'] else 'chr'+str(chrom).removeprefix('chr')


def gff(path):
    rows=[]
    for line in path.read_text().splitlines():
        if line.startswith('#') or not line.strip():continue
        c,_,_,start,end,_,strand,_,attrs=line.split('\t');a=dict(re.findall(r'(\w+)="([^";]*)"',attrs))
        rows.append(dict(site_id=a['peak'],chrom=canon(c),strand=strand,start=int(start)-1,end=int(end),
            endpoint=int(end)-1 if strand=='+' else int(start)-1))
    frame=pd.DataFrame(rows);assert frame.site_id.is_unique and frame.end.gt(frame.start).all();return frame


def make_windows(features,origins):
    x=origins.copy();x['left']=x.endpoint-x.start;x['right']=x.end-x.endpoint
    w=x.groupby('atlas_column_0based').agg(left=('left','median'),right=('right','median')).reset_index()
    w=features[['atlas_column_0based','atlas_pas_id','chrom','strand','cut_site','keep_for_apa_primary']].merge(w,on='atlas_column_0based',validate='one_to_one')
    w['endpoint']=w.cut_site-w.strand.eq('+').astype(int);w['start']=np.maximum(0,w.endpoint-np.floor(w.left)).astype(int)
    w['end']=(w.endpoint+np.ceil(w.right)).astype(int);assert w.end.gt(w.start).all()
    frames=[]
    for mode in ['transferred','primary_clipped','all_clipped']:
        a=w[w.keep_for_apa_primary.astype(str).str.lower().eq('true')].copy() if mode=='primary_clipped' else w.copy()
        a=a.sort_values(['chrom','strand','endpoint','atlas_pas_id']).reset_index(drop=True)
        if mode!='transferred':
            for _,ix in a.groupby(['chrom','strand']).groups.items():
                ix=np.asarray(ix);cuts=a.loc[ix,'endpoint'].to_numpy();assert (np.diff(cuts)>0).all()
                bounds=(cuts[:-1]+cuts[1:])//2
                a.loc[ix[1:],'start']=np.maximum(a.loc[ix[1:],'start'],bounds)
                a.loc[ix[:-1],'end']=np.minimum(a.loc[ix[:-1],'end'],bounds)
        frames.append(a.assign(mode=mode))
    return pd.concat(frames,ignore_index=True)


def recount_library(bam_path,binding,meta,features,samples,windows):
    import pysam
    barcode=dict(zip(binding.barcode.str.split('-').str[0],range(len(binding))));assert len(barcode)==len(binding)
    mci=pd.Index(meta.metacell_id).get_indexer(binding.metacell_id);assert (mci>=0).all()
    samples=samples.reset_index(drop=True);samples['site_index']=np.arange(len(samples));sample_targets=samples.target_index.to_numpy()
    tables={'sample':samples}|{m:windows[windows['mode'].eq(m)&windows.end.gt(windows.start)] for m in ['transferred','primary_clipped','all_clipped']}
    indexes={};intervals=[]
    for mode,frame in tables.items():
        for (chrom,strand),q in frame.groupby(['chrom','strand']):
            a=q.sort_values('start');starts=a.start.to_numpy();ends=a.end.to_numpy()
            indexes[mode,canon(chrom),strand]=(starts,ends,np.maximum.accumulate(ends),a.target_index.to_numpy(),a.site_index.to_numpy() if mode=='sample' else a.target_index.to_numpy())
        intervals.extend((canon(r.chrom),int(r.start),int(r.end)) for r in frame.itertuples())
    regions=[]
    for chrom,start,end in sorted(intervals):
        if regions and regions[-1][0]==chrom and start<=regions[-1][2]+100:regions[-1][2]=max(end,regions[-1][2])
        else:regions.append([chrom,start,end])
    shape=(len(meta),len(features));matrices={m:sparse.csr_matrix(shape,dtype=np.int64) for m in MODES if m!='transferred_observed'};audit=Counter()
    with pysam.AlignmentFile(bam_path,'rb') as bam:
        references={canon(c):c for c in bam.references};assert len(references)==len(bam.references)
        for chrom in sorted({r[0] for r in regions}):
            assert chrom in references,('missing_BAM_contig',chrom)
            molecules={m:set() for m in matrices}
            for _,start,end in (r for r in regions if r[0]==chrom):
                for read in bam.fetch(references[chrom],max(0,start),end):
                    audit['fetched_reads']+=1
                    if read.is_unmapped or read.is_duplicate or read.is_secondary or read.is_supplementary:continue
                    strand='-' if read.is_reverse else '+';ep=read.reference_start if read.is_reverse else read.reference_end-1
                    if not start<=ep<end:continue
                    if not read.has_tag('CB') or not read.has_tag('UB'):continue
                    cell=barcode.get(str(read.get_tag('CB')).split('-')[0]);umi=str(read.get_tag('UB'))
                    if cell is None:continue
                    audit['eligible_reads']+=1
                    for mode in tables:
                        index=indexes.get((mode,chrom,strand))
                        if index is None:continue
                        starts,ends,prefix,targets,sites=index;j=bisect.bisect_right(starts,ep)-1
                        while j>=0 and prefix[j]>ep:
                            if ends[j]>ep:
                                target=int(targets[j])
                                if mode=='sample':
                                    molecules['original'].add((int(sites[j]),cell,umi));molecules['atlas_dedup'].add((target,cell,umi))
                                else:molecules[mode].add((target,cell,umi))
                            j-=1
            assert molecules['all_clipped']<=molecules['primary_clipped']<=molecules['transferred']
            for mode,values in molecules.items():
                counted=Counter((int(mci[cell]),int(sample_targets[site]) if mode=='original' else site) for site,cell,umi in values)
                mat=sparse.csr_matrix((list(counted.values()),([k[0] for k in counted],[k[1] for k in counted])),shape=shape,dtype=np.int64)
                matrices[mode]+=mat
    assert (matrices['atlas_dedup']-matrices['original']).max()<=0
    called=np.zeros(len(features),bool);called[sample_targets]=True
    matrices['transferred_observed']=matrices['transferred'].multiply(called[None,:]).tocsr();matrices['transferred_observed'].eliminate_zeros()
    assert matrices['original'][:,~called].nnz==0
    return matrices,dict(audit)


def prepare(ctx,study,libraries,planned,native_root):
    study.load_gex();primary=study.primary
    f=primary[primary.matched_source_count_50nt.ge(1)&primary.primary_gene_id.isin(study.gex_index)].copy()
    f=f[f.groupby('primary_gene_id').atlas_pas_id.transform('size').ge(2)].copy().reset_index(drop=True);f['relative']=f.position
    mapping=ctx.read(study.root/'native_peak_to_atlas.tsv.gz');pieces=[]
    for lib in libraries.itertuples():
        path=ctx.register(native_root/'polyapipe'/lib.study_id/lib.library_id/'joint_polyA_peaks.gff');q=gff(path)
        old=mapping[mapping.library_id.eq(lib.library_id)][['site_id','atlas_column_0based','native_cut_site_1based']]
        q=q.merge(old,on='site_id',validate='one_to_one');assert len(q)==len(old)
        assert np.array_equal(q.endpoint+q.strand.eq('+').astype(int),q.native_cut_site_1based-q.strand.eq('-').astype(int))
        pieces.append(q.assign(library_id=lib.library_id))
    origins=pd.concat(pieces,ignore_index=True);windows=make_windows(study.f,origins)
    lookup=pd.Series(np.arange(len(f)),index=f.atlas_column_0based)
    for frame in [origins,windows]:frame['target_index']=frame.atlas_column_0based.map(lookup).fillna(-1).astype(int)
    windows=windows[windows.target_index.ge(0)].copy()
    for mode in MODES:
        bad=windows[windows['mode'].eq(mode)&windows.end.le(windows.start)].atlas_pas_id
        f[mode+'_available']=~f.atlas_pas_id.isin(bad)
    folder=study.id;ctx.table(f,folder+'/features.tsv.gz');ctx.table(windows,folder+'/windows.tsv.gz')
    ctx.table(origins[origins.target_index.ge(0)],folder+'/sample_windows.tsv.gz')
    jobs=[]
    for lib in libraries[libraries.library_key.isin(planned)].to_dict('records'):
        route,node,timepoint=planned[lib['library_key']]
        meta=study.meta[study.meta.library_key.eq(lib['library_key'])&study.meta.route.eq(route)&study.meta.node.eq(node)&study.meta.purity.ge(.8)].copy()
        assert len(meta)>=3 and meta.n_cells.sum()>=30 and meta.experimental_time.astype(str).eq(timepoint).all()
        binding=ctx.read(ctx.register(lib['members_path']));assert digest(lib['members_path'])==lib['members_sha256']
        binding=binding[binding.metacell_id.isin(meta.metacell_id)][['barcode','metacell_id']].copy();assert len(binding)==meta.n_cells.sum()
        bam=Path(lib['bam_path']);stat=bam.stat();assert stat.st_size==lib['bam_bytes'] and stat.st_mtime_ns==lib['bam_mtime_ns']
        assert digest(ctx.register(lib['bai_path']))==lib['bai_sha256']
        sample=origins[origins.library_id.eq(lib['library_id'])&origins.target_index.ge(0)].copy()
        prefix=folder+'/'+lib['library_id'];ctx.table(meta,prefix+'/metacells.tsv.gz');ctx.table(binding,prefix+'/cells.tsv.gz')
        jobs.append(dict(lib=lib,meta=meta,binding=binding,sample=sample,prefix=prefix))
    return f,windows,jobs


def estimate(ctx,study,features,jobs):
    math=ctx.functions(SOURCE,['library_stats','paired_stats']);study.load_gex();states={};genes=pd.Index(features.primary_gene_id.unique())
    for job in jobs:
        lib=job['lib'];meta=job['meta'];root=ctx.out/job['prefix'];matrices={mode:sparse.load_npz(ctx.register(root/(mode+'_mc.npz'))).tocsr() for mode in MODES}
        good=[]
        for mode,matrix in matrices.items():
            _,mask=math.library_stats(matrix,features,meta,study.gex_index,mode);good.append(mask)
        common=np.logical_and.reduce(good)
        gindex=study.gex_index.get_indexer(genes);assert (gindex>=0).all()
        ri=meta.index.to_numpy();gcounts=study.gex_counts[ri][:,gindex].toarray();gex=np.log2(1+gcounts/study.gex_total[ri,None]*1e6)
        gi=genes.get_indexer(features.primary_gene_id)
        projection=sparse.csr_matrix((np.ones(len(features),np.int64),(np.arange(len(features)),gi)),shape=(len(features),len(genes)))
        for view in ['reapplied','common_contributors']:
            for mode,matrix in matrices.items():
                stats,mask=math.library_stats(matrix,features,meta,study.gex_index,mode,common if view=='common_contributors' else None)
                counts=(matrix@projection).toarray();n=mask.sum(0)
                stats.update(mean_GEX=np.divide((np.expm1(gex*np.log(2))*mask).sum(0),n,out=np.zeros(len(n)),where=n>0),
                    APA_coverage=np.divide((counts*mask).sum(0),n,out=np.zeros(len(n)),where=n>0))
                states[view,mode,lib['library_key']]=stats
    pairs=[];design=study.pairs[study.pairs.contrast.eq('P_S')&study.pairs.design.eq('true_source')]
    available={j['lib']['library_key'] for j in jobs}
    design=design[design.library_key_a.isin(available)&design.library_key_b.isin(available)]
    for rec in design.to_dict('records'):
        for view in ['reapplied','common_contributors']:
            candidates={m:math.paired_stats(states[view,m,rec['library_key_a']],states[view,m,rec['library_key_b']],features)[1] for m in MODES}
            joint=np.logical_and.reduce(list(candidates.values()))
            for mode in MODES:
                mask=joint if view=='common_contributors' else candidates[mode]
                a,b=states[view,mode,rec['library_key_a']],states[view,mode,rec['library_key_b']]
                for i in np.flatnonzero(mask):
                    pairs.append(dict(rec,view=view,mode=mode,gene_id=genes[i],target_WPP=a['APA'][i],reference_WPP=b['APA'][i],delta=a['APA'][i]-b['APA'][i],
                        mean_GEX=(a['mean_GEX'][i]+b['mean_GEX'][i])/2,APA_coverage=(a['APA_coverage'][i]+b['APA_coverage'][i])/2,
                        n_mc_a=int(a['n_mc'][i]),n_mc_b=int(b['n_mc'][i]),n_cells_a=int(a['n_cells'][i]),n_cells_b=int(b['n_cells'][i])))
    cols=list(study.pairs.columns)+['view','mode','gene_id','target_WPP','reference_WPP','delta','mean_GEX','APA_coverage','n_mc_a','n_mc_b','n_cells_a','n_cells_b']
    p=pd.DataFrame(pairs,columns=cols);result=[]
    keys=['species','study_id','route','contrast','design','view','mode','gene_id']
    for key,q in p.groupby(keys):
        w=weights(q);row=dict(zip(keys,key))
        for value in ['target_WPP','reference_WPP','delta','mean_GEX','APA_coverage']:row[value]=float(w@q[value])
        row['gene_symbol']=features.loc[features.primary_gene_id.eq(row['gene_id']),'primary_gene_symbol'].iloc[0]
        row['PAS_number']=int(features.primary_gene_id.eq(row['gene_id']).sum());result.append(row)
    effects=pd.DataFrame(result,columns=keys+['target_WPP','reference_WPP','delta','mean_GEX','APA_coverage','gene_symbol','PAS_number'])
    common=effects[effects.view.eq('reapplied')].groupby(['route','gene_id'])['mode'].nunique();good=common[common.eq(6)].reset_index()[['route','gene_id']]
    fixed=effects[effects.view.eq('reapplied')].merge(good,on=['route','gene_id'],validate='many_to_one').assign(view='fixed_genes')
    fixedp=p[p.view.eq('reapplied')].merge(good,on=['route','gene_id'],validate='many_to_one').assign(view='fixed_genes')
    return pd.concat([effects,fixed],ignore_index=True),pd.concat([p,fixedp],ignore_index=True)


def run(ctx):
    assert (ctx.native/'annotation.complete.json').exists(),'Full new PAS annotation must finish before recount starts'
    libraries=ctx.read(ctx.resource('native_libraries'));protocol=json.loads(ctx.resource('legacy_endpoint_ledger').read_text())
    native_root=Path(ctx.config['native_caller_root']);groups={};total=0
    for rec in protocol['recount_inventory']:
        wanted=groups.setdefault(rec['study'],{})
        for role,libs in rec['libraries'].items():
            for lib,timepoint in libs:
                entry=(rec['route'],role,str(timepoint));assert lib not in wanted or wanted[lib]==entry;wanted[lib]=entry
    assert sum(map(len,groups.values()))==38
    summaries=[];all_effects=[];all_pairs=[]
    for sid,wanted in groups.items():
        study=Study(ctx,sid);lib=libraries[libraries.study_id.eq('WOT_GSE115943' if sid=='GSE115943' else sid)]
        assert set(wanted)<=set(lib.library_key)
        features,windows,jobs=prepare(ctx,study,lib,wanted,native_root)
        def execute(job):
            started=time.monotonic();l=job['lib'];matrices,audit=recount_library(l['bam_path'],job['binding'],job['meta'],features,job['sample'],windows)
            folder=ctx.out/job['prefix']
            for mode,mat in matrices.items():sparse.save_npz(folder/(mode+'_mc.npz'),mat)
            stat=Path(l['bam_path']).stat();assert stat.st_size==l['bam_bytes'] and stat.st_mtime_ns==l['bam_mtime_ns']
            record=dict(study_id=sid,library_key=l['library_key'],counts={m:int(x.sum()) for m,x in matrices.items()},
                seconds=time.monotonic()-started,bam_path=l['bam_path'],bam_bytes=stat.st_size,bam_mtime_ns=stat.st_mtime_ns,
                count_policy='CB/UB/end-point sensitivity; native UR remains primary',native_equivalence_asserted=False,**audit)
            (folder/'recount-audit.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record),flush=True);return record
        with ThreadPoolExecutor(max_workers=ctx.config.get('recount_workers',2)) as pool:
            futures=[pool.submit(execute,j) for j in jobs]
            for future in as_completed(futures):summaries.append(future.result());total+=1
        effects,pairs=estimate(ctx,study,features,jobs);all_effects.append(effects);all_pairs.append(pairs)
    assert total==38
    ctx.table(pd.concat(all_effects,ignore_index=True),'route_effects.tsv.gz');ctx.table(pd.concat(all_pairs,ignore_index=True),'pair_effects.tsv.gz')
    ctx.table([dict(r,counts=json.dumps(r['counts'])) for r in summaries],'library_recount_audit.tsv')
    ctx.finish(dict(all_38_prespecified_libraries=True,new_native_peak_windows=True,six_rules_recounted=True,
        fixed_genes_and_common_contributors_rebuilt=True,BAMs_read_only=True,count_subset_checks=True),
        primary_counts='native polyApipe UR',sensitivity_counts='CB/UB/read-end; original denotes sample-window sensitivity, not native caller counts',P_M_scope='not recounted')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);a=p.parse_args();run(Context(a.config,'recount',a.out))
