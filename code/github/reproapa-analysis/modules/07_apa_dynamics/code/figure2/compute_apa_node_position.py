"""Average native-study relative PAS position on coarse biological nodes."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

MAIN_NODES=['S','R','E','T','M','B','N','I','K']
ROUTES={'wot_common_to_2i':'WOT','gse103221_cip':'CIP','gse103221_osk':'OSK',
    'gse171127_d3_getm':'171127_D3_GETM','gse171127_d3_oskm':'171127_D3_OSKM',
    'gse171127_d6_pooled':'171127_late_pooled','gse171127_d12_pooled':'171127_late_pooled'}


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()


def run(root,out):
    out.mkdir(parents=True,exist_ok=True)
    meta=pd.read_csv(root/'results_apa_metacell_study_local_v1/inputs/metacell_coarse_annotations.tsv',sep='\t')
    all_nodes=[];all_libraries=[];audits=[];selections=[];structures=[]
    for source in sorted((root/'results_metacell_pas_v7_v1/studies').iterdir()):
        study=source.name
        hashes={n:sha(source/n) for n in ['features.tsv.gz','metacells.tsv.gz','metacell_by_pas_counts.npz']}
        m=meta[meta.pas_dataset_id.eq(study)].sort_values('pas_matrix_row_index_0based').reset_index(drop=True)
        original=pd.read_csv(source/'metacells.tsv.gz',sep='\t')
        assert np.array_equal(m.metacell_id,original.metacell_id)
        assert np.array_equal(m.pas_matrix_row_index_0based,np.arange(len(m)))
        m['display_route']=m.dataset.map(ROUTES)
        f=pd.read_csv(source/'features.tsv.gz',sep='\t',low_memory=False)
        f=f[f.keep_for_apa_primary & f.primary_gene_id.notna()].sort_values(['primary_gene_id','cut_site','atlas_pas_id']).reset_index(drop=True)
        genes=pd.Index(f.primary_gene_id.unique());gi=genes.get_indexer(f.primary_gene_id)
        starts=np.r_[0,np.flatnonzero(np.diff(gi))+1];ends=np.r_[starts[1:],len(f)]
        chrom_n=f.groupby('primary_gene_id',sort=False).chrom.nunique().to_numpy()
        strand_n=f.groupby('primary_gene_id',sort=False).strand.nunique().to_numpy()
        strand=f.strand.iloc[starts].to_numpy()
        pos=f.cut_site.to_numpy(float)*np.where(f.strand.eq('+'),1.,-1.)
        lo=np.minimum.reduceat(pos,starts);hi=np.maximum.reduceat(pos,starts)
        coordinate_ok=(chrom_n==1)&(strand_n==1)&np.isin(strand,['+','-'])&(hi>lo)&(ends-starts>=2)
        relative=np.divide(pos-lo[gi],(hi-lo)[gi],out=np.zeros(len(f)),where=(hi-lo)[gi]>0)
        assert ((relative>=0)&(relative<=1)).all()
        mapping=sparse.csr_matrix((np.ones(len(f)),(np.arange(len(f)),gi)),shape=(len(f),len(genes)))
        weights=sparse.csr_matrix((relative,(np.arange(len(f)),gi)),shape=mapping.shape)
        raw=sparse.load_npz(source/'metacell_by_pas_counts.npz').tocsr()
        c=raw[:,f.feature_index_1based.to_numpy(int)-1].tocsr();del raw
        totals=(c@mapping).toarray();numerator=(c@weights).toarray()
        index=np.divide(numerator,totals,out=np.full(totals.shape,np.nan),where=totals>=20)
        good=(totals>=20)&coordinate_ok
        assert np.isnan(index[totals<20]).all()
        assert np.nanmin(index)>=-1e-10 and np.nanmax(index)<=1+1e-10
        # Independent direct checks of weighted relative positions on raw PAS slices.
        checks=0
        for j in np.flatnonzero(coordinate_ok)[::max(1,int(coordinate_ok.sum()/40))]:
            eligible_mc=np.flatnonzero(good[:,j])
            for i in eligible_mc[:3]:
                v=c[i,starts[j]:ends[j]].toarray().ravel()
                direct=float(np.dot(v,relative[starts[j]:ends[j]])/v.sum())
                assert np.isclose(direct,index[i,j],atol=1e-12)
                checks+=1
        annotation=pd.DataFrame(dict(study=study,gene_id=genes,gene_symbol=f.primary_gene_symbol.iloc[starts].to_numpy(),
            n_primary_pas=ends-starts,coordinate_eligible=coordinate_ok))
        structures.append(annotation)
        valid_meta=m.coarse_label_purity.ge(.8)&m.dominant_coarse_state.isin(MAIN_NODES)&m.day.ge(0)
        assert m.loc[valid_meta,'display_route'].notna().all()
        for route,route_meta in m[valid_meta].groupby('display_route'):
            selected=[];groups=[]
            for (node,lib),gm in route_meta.groupby(['dominant_coarse_state','library_key']):
                if len(gm)<3 or gm.n_cells.sum()<30:continue
                ix=gm.index.to_numpy();selected.extend(ix);groups.append((node,lib,ix))
            selected=np.array(sorted(set(selected)),dtype=int)
            if not len(selected):continue
            # Multi-PAS support is evaluated over this route; switches may use one PAS per node.
            z=c[selected].tocoo();keep=good[selected][z.row,gi[z.col]]
            count=np.bincount(z.col[keep],weights=z.data[keep],minlength=len(f))
            detected=np.bincount(z.col[keep],minlength=len(f))
            support=np.add.reduceat(((count>=10)&(detected>=2)).astype(int),starts)
            route_ok=coordinate_ok&(support>=2)
            libs=[]
            for node,lib,ix in groups:
                ng=good[ix].sum(0);ok=(ng>=3)&route_ok
                mu=np.divide(np.where(good[ix],index[ix],0).sum(0),ng,out=np.full(len(genes),np.nan),where=ng>0)
                mu[~ok]=np.nan
                part=annotation.loc[ok,['study','gene_id','gene_symbol']].copy()
                part['route']=route;part['node']=node;part['library_key']=lib
                part['day']=float(m.day.iloc[ix[0]]);part['position_index']=mu[ok]
                part['qualified_metacells']=ng[ok]
                libs.append(part)
                sm=m.iloc[ix][['metacell_id','library_key','dataset','day','dominant_coarse_state',
                    'coarse_label_purity','n_cells','pas_matrix_row_index_0based']].copy()
                sm['route']=route;sm['study']=study;selections.append(sm)
            libtable=pd.concat(libs);all_libraries.append(libtable)
            for node in sorted({g[0] for g in groups},key=MAIN_NODES.index):
                values=libtable[libtable.node.eq(node)]
                ag=values.groupby('gene_id').agg(position_index=('position_index','mean'),
                    n_libraries=('library_key','nunique'),n_qualified_metacells=('qualified_metacells','sum'),
                    n_dates=('day','nunique'),day_min=('day','min'),day_max=('day','max')).reindex(genes)
                row=annotation.copy();row['route']=route;row['node']=node
                for col in ag:row[col]=ag[col].to_numpy()
                row['supported_pas_in_route']=support
                row['status']=np.select([~coordinate_ok,~route_ok,row.n_libraries.isna()],
                    ['unusable_relative_coordinate','fewer_than_two_supported_pas','insufficient_node_coverage'],default='available')
                row['statistical_filter']='none'
                all_nodes.append(row)
        assert hashes=={n:sha(source/n) for n in hashes}
        audits.append(dict(study=study,source_sha256=hashes,metacells=len(m),primary_pas=len(f),
            genes=len(genes),relative_coordinate_genes=int(coordinate_ok.sum()),direct_mc_recomputations=checks,
            original_matrix_unchanged=True))
        print(json.dumps(audits[-1]),flush=True)
    nodes=pd.concat(all_nodes,ignore_index=True)
    library=pd.concat(all_libraries,ignore_index=True)
    assert not nodes.duplicated(['study','route','node','gene_id']).any()
    assert not library.duplicated(['route','node','library_key','gene_id']).any()
    assert nodes.loc[nodes.status.ne('available'),'position_index'].isna().all()
    reference=nodes[nodes.node.eq('S')][['route','gene_id','position_index']].rename(columns={'position_index':'source_S_index'})
    nodes=nodes.merge(reference,on=['route','gene_id'],how='left',validate='many_to_one')
    nodes['delta_from_S']=nodes.position_index-nodes.source_S_index
    nodes.to_csv(out/'gene_node_positions.tsv.gz',sep='\t',index=False)
    library.to_csv(out/'gene_library_node_positions.tsv.gz',sep='\t',index=False)
    pd.concat(selections).drop_duplicates().to_csv(out/'selected_metacells.tsv.gz',sep='\t',index=False)
    pd.concat(structures).to_csv(out/'gene_coordinate_availability.tsv.gz',sep='\t',index=False)
    nodes.groupby(['route','node']).agg(genes=('position_index','count'),
        median_libraries=('n_libraries','median'),median_metacells=('n_qualified_metacells','median')).to_csv(out/'node_coverage.tsv',sep='\t')
    (out/'readback_audit.json').write_text(json.dumps(dict(status='PASS',source_audits=audits,
        node_rows=len(nodes),library_rows=len(library),coordinate_range=[0,1],
        negative_day_reference_excluded=True,absolute_lengths=False,statistical_tests=False),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();run(a.root,a.out)
