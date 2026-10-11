"""Bind frozen cell coordinates and verified membership to the candidate renderer."""
import json
from pathlib import Path
import pandas as pd

from polyapipe_prepare import digest,save

PROJECT=Path(__file__).resolve().parents[2]
OUT=PROJECT/'notes/2026-10-03-polyapipe-rerun/render_resources'


def main():
    OUT.mkdir(exist_ok=False)
    registry=pd.read_csv(PROJECT/'analysis_code/pas_rerun/frozen-membership.tsv',sep='\t')
    sources={};frames=[]
    for row in registry[registry.role.eq('membership')].itertuples():
        path=Path(row.path);assert digest(path)==row.sha256,path;sources[str(path)]=row.sha256
        frames.append(pd.read_csv(path,sep='\t',compression='gzip',usecols=['cell_id','metacell_id']))
    members=pd.concat(frames,ignore_index=True);assert len(members)==551201 and members.cell_id.is_unique
    embedding=PROJECT.parent/'科研/projects/2026-09-10_human-mouse-time-aware-umap/data/results_raw_cell_species_integrated_v3_faiss2/state_and_trajectory_display_v13/raw_cell_stateflow_groups.tsv.gz'
    sources[str(embedding)]=digest(embedding)
    cells=pd.read_csv(embedding,sep='\t',usecols=['cell_id','species_key','study_id','library_key','stateflow_code','state_group','UMAP1','UMAP2'])
    cells=cells.merge(members,on='cell_id',validate='one_to_one');assert len(cells)==551201 and cells.metacell_id.nunique()==11104
    # Verify the cell-state table is exactly the one used to freeze the framework.
    counts=pd.crosstab(cells.metacell_id,cells.state_group)
    frozen=pd.read_csv(PROJECT/'notes/2026-10-03-polyapipe-rerun/frozen_framework/metacell_states.tsv.gz',sep='\t').set_index('metacell_id')
    purity=counts.max(axis=1)/counts.sum(axis=1)
    import numpy as np
    np.testing.assert_allclose(purity.reindex(frozen.index),frozen.state_purity,rtol=0,atol=1e-15)
    assert counts.idxmax(axis=1).reindex(frozen.index).eq(frozen.state_group).all()
    save(cells,OUT/'frozen_cell_embedding.tsv.gz')
    backup=json.loads((PROJECT/'notes/2026-10-03-polyapipe-rerun/manuscript-backup.json').read_text())
    original=Path(backup['backup'])/'Fig/Figure.1.pdf'
    baseline=pd.read_csv(PROJECT/'analysis_code/pas_rerun/baseline-files.tsv',sep='\t').set_index('file')
    assert digest(original)==baseline.loc['Fig/Figure.1.pdf','sha256']
    base=json.loads((PROJECT/'notes/2026-10-03-polyapipe-rerun/statistics_package_v12/config.json').read_text())
    config=dict(analysis_root=str(PROJECT/'notes/2026-10-03-polyapipe-rerun/statistics_results'),
        run_id=base['run_id'],input_fingerprint=base['input_fingerprint'],resources={
            'frozen_cell_embedding':dict(path=str(OUT/'frozen_cell_embedding.tsv.gz'),sha256=digest(OUT/'frozen_cell_embedding.tsv.gz')),
            'original_figure1':dict(path=str(original),sha256=digest(original))})
    historical=PROJECT/'notes/2026-09-22-pas-support-layout/source_data/study_step_effects.tsv'
    config['resources']['historical_counting_steps']=dict(path=str(historical),sha256=digest(historical),evidence_scope='historical_original_pipeline')
    (OUT/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    (OUT/'audit.json').write_text(json.dumps(dict(cells=len(cells),metacells=cells.metacell_id.nunique(),
        frozen_purity_recovered=True,source_sha256=sources),indent=2)+'\n')
    print(json.dumps(dict(cells=len(cells),config=str(OUT/'config.json'))))


if __name__=='__main__':main()
