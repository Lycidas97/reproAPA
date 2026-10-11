import ast
import sys
import igraph as ig
import leidenalg
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from common import *


def run():
    out=RUN/'A07';x=obs().set_index('cell_id')
    parameters=pd.read_csv(ROOT/'inputs/S4-local-partitions.tsv',sep='\t').set_index('analysis_unit')
    scopes=json.loads((RUN/'runtime/scopes.json').read_text())
    records=[];source_records=[]
    for unit,row in parameters.iterrows():
        folder=out/unit;folder.mkdir(exist_ok=True)
        if unit=='wot_common_to_2i':
            wot=Path('/root/nfsdata/multica/worktrees/apadetect/DOC-54-wot-readiness/analyses/wot-unsupervised-first-topology-audit-v1')
            source=wot/'intermediate/knn_common_to_2i_program_blind_seed218_k50.npz'
            coords_path=wot/'unsupervised_cell_coordinates.tsv.gz'
            z=np.load(source)
            code=(RUN/'runtime/sources/wot_run_unsupervised.py').read_text()
            function=next(n for n in ast.parse(code).body if isinstance(n,ast.FunctionDef) and n.name=='connectivity')
            namespace={'np':np,'sparse':sparse}
            exec(compile(ast.Module(body=[function],type_ignores=[]),str(source),'exec'),namespace)
            graph=namespace['connectivity'](z['indices'],z['distances'],30)
            cells=pd.read_csv(coords_path,sep='\t');cells=cells[cells.scope.eq('common_to_2i')].copy()
            cells['cell_id']=cells.cell_id.str.replace('GSE115943_GSE122662|','GSE115943|',regex=False)
            selected=float(row.selected_resolution);grid=[.3,.45,.6]
            src=[source,coords_path]
        else:
            scope=Path(scopes[unit]);source=scope/'embedding/knn_graph.npz'
            graph=sparse.load_npz(source).tocsr()
            cells=pd.read_csv(scope/'prepared/cells.tsv.gz',sep='\t')
            frozen=pd.read_csv(scope/'frozen/cell_clusters.tsv.gz',sep='\t')
            assert np.array_equal(cells.cell_id,frozen.cell_id)
            selected=json.loads((scope/'frozen/metadata.json').read_text())['selected_resolution']
            grid=sorted(json.loads((scope/'sweep/metadata.json').read_text())['resolutions'])
            src=[source,scope/'prepared/cells.tsv.gz',scope/'frozen/metadata.json',scope/'sweep/metadata.json']
        assert graph.shape==(int(row.n_cells),int(row.n_cells)) and len(cells)==row.n_cells
        assert cells.cell_id.is_unique and cells.cell_id.isin(x.index).all()
        labels=x.reindex(cells.cell_id).state_key.to_numpy()
        assert len(np.unique(labels))==row.n_local_states
        variants=sorted(set([selected,max([g for g in grid if g<selected-1e-6],default=selected/2),
                             min([g for g in grid if g>selected+1e-6],default=selected*1.5)]))
        write_json(folder/'config.json',dict(unit=unit,baseline_resolution=selected,resolutions=variants,
            baseline_seed=int(row.seed),variant_seeds=[SEED,SEED+1,SEED+2],cells=len(cells),
            fixed_original_graph=True,WOT_unified_pipeline_in_this_stage=False))
        source_records.extend(dict(path=str(p),sha256=sha(p)) for p in src)
        write_json(out/'source_receipt.json',source_records)
        upper=sparse.triu(graph,k=1).tocoo()
        g=ig.Graph(n=graph.shape[0],edges=list(zip(upper.row.tolist(),upper.col.tolist())),directed=False)
        g.es['weight']=upper.data.astype(float).tolist()
        del upper,graph
        combinations=[(selected,int(row.seed),'baseline_replay')]+[(r,s,'sensitivity') for r in variants for s in [SEED,SEED+1,SEED+2]]
        for resolution,seed,role in combinations:
            progress('A07','RUNNING',stage='Leiden_resolution_and_seed',unit=unit,resolution=resolution,seed=seed,role=role)
            partition=leidenalg.find_partition(g,leidenalg.RBConfigurationVertexPartition,
                weights='weight',resolution_parameter=resolution,seed=seed,n_iterations=-1)
            alternative=np.asarray(partition.membership,np.int32)
            metrics=dict(unit=unit,resolution=resolution,seed=seed,role=role,
                original_states=len(np.unique(labels)),alternative_clusters=len(np.unique(alternative)),
                ARI=float(adjusted_rand_score(labels,alternative)),NMI=float(normalized_mutual_info_score(labels,alternative)))
            name=f'{role}_r{resolution:.8f}_seed{seed}'
            result=cells[['cell_id']].copy();result['original_state']=labels;result['alternative_cluster']=alternative
            result.to_csv(folder/(name+'.tsv.gz'),sep='\t',index=False)
            pd.crosstab(result.original_state,result.alternative_cluster).to_csv(folder/(name+'_overlap.tsv'),sep='\t')
            records.append(metrics);pd.DataFrame(records).to_csv(out/'partition_metrics.tsv',sep='\t',index=False)
            if role=='baseline_replay':
                # A changed installed Leiden implementation must remain visible, not hidden by choosing a seed.
                write_json(folder/'baseline_replay.json',dict(**metrics,exact_membership_reproduced=metrics['ARI']>1-1e-12))
        del g
    finish('A07','fixed_graph_partition_grid',runs=len(records),analysis_units=len(parameters),
        next_steps=['WOT full matched-universe uniform chooseR pipeline','A04 family mapping of alternatives','A05 APA reaggregation'],
        boundary='Original graphs fixed; seed/resolution sensitivity is separate from whole-pipeline validation.')


if __name__=='__main__': execute('A07',run)
