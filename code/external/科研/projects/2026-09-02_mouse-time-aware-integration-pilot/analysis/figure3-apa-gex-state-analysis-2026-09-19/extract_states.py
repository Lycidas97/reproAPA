from pathlib import Path
import json,gc,time
import numpy as np,pandas as pd
import anndata as ad
from scipy import sparse
ROOT=Path('/root/nfsdata/apadetect/analysis_projects'); M=ROOT/'mouse_time_aware_integration_pilot_v1'; H=ROOT/'human_time_aware_integration_v1'; O=ROOT/'figure3_apa_gex_state_20260919';O.mkdir(exist_ok=True)
REP=pd.read_csv(O/'library_sampling_units.tsv',sep='\t').set_index('library_key'); manifest={};checks=[]
def read(p,**kw):
 manifest[str(p)]={'bytes':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns};return pd.read_csv(p,sep='\t',**kw)
def save(d,p):d.to_csv(p,sep='\t',index=False,float_format='%.10g')
hg=ad.read_h5ad(H/'results_metacell_v1/release/human_metacell_all_gene_counts.h5ad')
for prof in sorted((M/'results_apa_majority_shared_v1/profiles').iterdir()):
 if not (prof/'selected_metacells.tsv.gz').exists():continue
 out=O/prof.name;out.mkdir(exist_ok=True)
 if (out/'complete.json').exists():continue
 m=read(prof/'selected_metacells.tsv.gz').sort_values('matrix_row').reset_index(drop=True);sp=m.species.iloc[0];src=(M/'results_metacell_pas_v7_v1/studies' if sp=='mouse' else H/'results_metacell_pas_v2/studies')/prof.name
 f=read(src/'features.tsv.gz');f=f[f.keep_for_apa_primary&f.primary_gene_id.notna()].copy();sup=(M/'results_pas_annotation_support_v1/studies' if sp=='mouse' else M/'results_apa_majority_shared_v1/human_annotation')/prof.name/'primary_pas_reference_matches.tsv.gz';f=f.merge(read(sup,usecols=['atlas_pas_id','matched_source_count_50nt']),on='atlas_pas_id',validate='one_to_one');v=f.groupby('primary_gene_id').agg(n=('cut_site','nunique'),c=('chrom','nunique'),s=('strand','nunique'));f=f[f.primary_gene_id.isin(v.query('n>=2 and c==1 and s==1').index)].copy();x=f.cut_site*np.where(f.strand.eq('+'),1,-1);lo=x.groupby(f.primary_gene_id).transform('min');hi=x.groupby(f.primary_gene_id).transform('max');f['position']=(x-lo)/(hi-lo)
 raw=sparse.load_npz(src/'metacell_by_pas_counts.npz').tocsr();ori=read(src/'metacells.tsv.gz');assert np.array_equal(ori.metacell_id.iloc[m.matrix_row],m.metacell_id);assert np.issubdtype(raw.data.dtype,np.integer)
 rows=[];grows=[];design=[]
 for lib,mm in m.groupby('library_key',sort=False):
  mm=mm.reset_index(drop=True);gp=M/'results_mouse_metacell_gex_completion_v1'/lib/'all_gene_counts.h5ad' if sp=='mouse' else H/'results_metacell_v1/release/human_metacell_all_gene_counts.h5ad';manifest[str(gp)]={'bytes':gp.stat().st_size,'mtime_ns':gp.stat().st_mtime_ns};g=ad.read_h5ad(gp) if sp=='mouse' else hg
  ids=pd.Index(g.var.gene_ids.astype(str) if 'gene_ids' in g.var else g.var_names.astype(str)).str.replace(r'\.[0-9]+$','',regex=True);assert ids.is_unique;sy=next((c for c in ['gene_symbols','gene_symbol','symbol','gene_name'] if c in g.var),None);symbols=g.var[sy].astype(str).to_numpy() if sy else g.var_names.astype(str).to_numpy();ix=g.obs_names.get_indexer(mm.metacell_id);assert (ix>=0).all();gx=g.X[ix].tocsr();total=np.asarray(gx.sum(1)).ravel();assert (total>0).all() and (gx.data==np.floor(gx.data)).all();assert np.array_equal(g.obs.iloc[ix].n_cells,mm.n_cells);gm=gx.toarray();gl=np.log2(1+gm/total[:,None]*1e6)
  b=dict(species=sp,study_id=mm.study_id.iloc[0],route=mm.route.iloc[0],library_key=lib,time=str(mm.experimental_time.iloc[0]),sample_unit=REP.loc[lib,'sampling_unit'],replicate_status=REP.loc[lib,'replicate_status']);checks.append(b|dict(exact_MC=True,integer_GEX=True,integer_PAS=True))
  fi=ids.get_indexer(f.primary_gene_id);fl=f[fi>=0].reset_index(drop=True);gg0=fi[fi>=0];cc0=raw[mm.matrix_row.to_numpy()][:,fl.feature_index_1based.to_numpy()-1]
  for setting,purity,supported in [('supported50_p08',.8,True),('supported50_p09',.9,True),('expanded_p08',.8,False)]:
   keep=fl.matched_source_count_50nt.to_numpy()>=1 if supported else np.ones(len(fl),bool);ff=fl[keep];gg=gg0[keep];cc=cc0[:,keep].tocsr();pos=ff.position.to_numpy();mp=sparse.csr_matrix((np.ones(len(gg)),(np.arange(len(gg)),gg)),shape=(len(gg),len(ids)));wm=sparse.csr_matrix((pos,(np.arange(len(gg)),gg)),shape=mp.shape);ct=(cc@mp).toarray();nu=(cc@wm).toarray();ap=np.divide(nu,ct,out=np.zeros_like(nu),where=ct>0);good=ct>=20
   features={k:np.flatnonzero(gg==k) for k in np.unique(gg)}
   for st in ['S','R','E','P','M']:
    ii=np.flatnonzero(mm.node.eq(st)&mm.purity.ge(purity))
    if len(ii)<3 or mm.iloc[ii].n_cells.sum()<30:continue
    ca=gm[ii].sum(0);cpm=ca/total[ii].sum()*1e6;allmean=gl[ii].mean(0)
    if setting=='supported50_p08':
     d=pd.DataFrame(dict(gene_id=ids,gene_symbol=symbols,count=ca,CPM=cpm,GEX=allmean))
     for k,vv in b.items():d[k]=vv
     d['state']=st;d['n_mc']=len(ii);d['n_cells']=mm.iloc[ii].n_cells.sum();grows.append(d[d['count']>0]);design.append(b|dict(state=st,n_mc=len(ii),n_cells=int(mm.iloc[ii].n_cells.sum())))
    goodi=good[ii];n=goodi.sum(0);cells=(goodi*mm.iloc[ii].n_cells.to_numpy()[:,None]).sum(0);valid=np.flatnonzero((n>=3)&(cells>=30));co=cc[ii].tocoo();q=good[ii[co.row],gg[co.col]];sc=np.bincount(co.col,weights=co.data*q,minlength=len(gg));nd=np.bincount(co.col,weights=q.astype(int),minlength=len(gg));ame=np.divide((ap[ii]*goodi).sum(0),n,out=np.zeros(len(ids)),where=n>0);gme=np.divide((gl[ii]*goodi).sum(0),n,out=np.zeros(len(ids)),where=n>0)
    for k in valid:
     ss=features.get(k,[])
     if len(ss)<2:continue
     rows.append(b|dict(setting=setting,state=st,gene_id=ids[k],gene_symbol=symbols[k],APA=ame[k],GEX=gme[k],GEX_all=allmean[k],GEX_CPM=cpm[k],GEX_count=int(ca[k]),n_mc=int(n[k]),n_cells=int(cells[k]),pas_counts=int(ct[ii,k][goodi[:,k]].sum()),site_counts=json.dumps(sc[ss].astype(int).tolist(),separators=(',',':')),site_mc=json.dumps(nd[ss].astype(int).tolist(),separators=(',',':'))))
  print('LIB',prof.name,lib,len(rows),flush=True);del gm,gl,ct,nu,ap,gx,cc0;gc.collect()
 save(pd.DataFrame(rows),out/'apa_gex_states.tsv.gz');save(pd.concat(grows,ignore_index=True),out/'full_gex_states.tsv.gz');save(pd.DataFrame(design),out/'sampling_design.tsv');(out/'complete.json').write_text(json.dumps(dict(APA_rows=len(rows),GEX_rows=sum(map(len,grows)))));save(pd.DataFrame(checks),O/'count_checks.tsv');(O/'manifest.json').write_text(json.dumps(manifest,indent=2));print('DONE',prof.name,flush=True);del raw,rows,grows;gc.collect()
(O/'complete.ok').write_text('complete\n')
