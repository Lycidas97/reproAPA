"""Freeze executable statistics and read-only resources for the running native intake."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

from polyapipe_prepare import digest

PROJECT=Path(__file__).resolve().parents[2]
PROJECTS=PROJECT.parent/'科研/projects'
REMOTE=Path('/root/nfsdata/apadetect/analysis_projects/reproapa_polyapipe_rerun_20261003_v1')


def main(a):
    a.out.mkdir(parents=True,exist_ok=False);code=a.out/'code';code.mkdir();resources=a.out/'resources';resources.mkdir();legacy=a.out/'legacy';legacy.mkdir()
    dest=REMOTE/a.out.name
    for path in Path(__file__).parent.glob('*.py'):
        if path.name=='downstream_common.py':continue
        shutil.copy2(path,code/path.name)
    source_paths=[
        'notes/2026-09-30-multimetric-statistics/inference/analyze.py',
        'notes/2026-09-30-multimetric-statistics/metric_inference/analyze.py',
        'notes/2026-09-30-multimetric-statistics/pool_inference/analyze.py',
        'notes/2026-09-30-multimetric-statistics/pair_count_tests/analyze.py',
        'notes/2026-09-30-multimetric-statistics/nulls/analyze.py',
        'notes/2026-10-02-null-model-controls/common.py','notes/2026-10-02-null-model-controls/recurrence.py',
        'notes/2026-10-02-null-model-controls/functional.py','notes/2026-10-02-reviewer-followup/trajectory_noise.py',
        'notes/2026-10-02-reviewer-followup/apa_followup.py','notes/2026-10-02-reviewer-followup/same_state.py',
        'notes/2026-09-14-integration-reliability/evaluate_reaggregated_apa_effects.py']
    source_paths+=['github/reproapa-analysis/revisions/20260930/external/2026-09-26_reproapa-go-analysis/scripts/'+n
        for n in ['run_go.py','run_go_nested.py','run_kegg_nested.py']]
    source_paths+=['../科研/projects/2026-09-02_mouse-time-aware-integration-pilot/scripts/'+n
        for n in ['test_figure2_apa_manifold_coherence.py','robustness_figure2_apa_knn.py']]
    overrides={}
    for relative in source_paths:
        path=PROJECT/relative;name=digest(path)[:12]+'_'+path.name;shutil.copy2(path,legacy/name)
        overrides[relative]=dict(path=str(dest/'legacy'/name),sha256=digest(path))
    go=PROJECTS/'2026-09-26_reproapa-go-analysis'
    mapping=dict(ortholog=PROJECTS/'2026-09-13_reproapa-figure4-real/data/figure4-program-v9/ortholog_universe.tsv',
        go_obo=go/'sources/go-basic.obo',go_human=go/'sources/human_go_ids.tsv',go_mouse=go/'sources/mouse_go_ids.tsv',
        expression_cycle_inputs=PROJECT/'notes/2026-09-26-state-comparability-audit/expression_inputs.npz',
        frozen_metacell_states=PROJECT/'notes/2026-10-03-polyapipe-rerun/frozen_framework/metacell_states.tsv.gz')
    integration=PROJECTS/'2026-09-10_human-mouse-time-aware-umap/data/results_integrated_species_scvi_v2'
    mapping.update(expression_graph_nodes=integration/'common_metacell_nodes.tsv.gz',
        expression_latent=integration/'species_balanced_scvi_latent.npy',expression_graph=integration/'local_expression_kernel.npz')
    for sp,org,title in [('human','hsa','Homo_sapiens'),('mouse','mmu','Mus_musculus')]:
        mapping['kegg_'+sp+'_gene_info']=go/'sources_kegg'/(title+'.gene_info.gz')
        for kind in ['ncbi_conv','gene_pathway','pathways']:mapping['kegg_'+org+'_'+kind]=go/'sources_kegg'/f'{org}_{kind}.tsv'
    resource_manifest={}
    for key,path in mapping.items():
        assert not (resources/path.name).exists();shutil.copy2(path,resources/path.name)
        resource_manifest[key]=dict(path=str(dest/'resources'/path.name),sha256=digest(path))
    jobs={}
    def job(node,script,depends=(),argument=True):
        jobs[node]=dict(script=str(dest/'code'/script),arguments=['--node',node] if argument else [],depends_on=list(depends))
    job('state_vectors','downstream_bundle.py',argument=False)
    job('endpoints','downstream_cohorts.py',['state_vectors']);job('metrics','downstream_cohorts.py',['endpoints'])
    for node in ['well_inference','pool_inference','count_inference']:job(node,'downstream_inference.py',['metrics'])
    job('recurrence','downstream_recurrence.py',['metrics'])
    for node in ['cross_sharing','recurrence_controls']:job(node,'downstream_recurrence.py',['recurrence'])
    job('profiles','downstream_profiles.py',['state_vectors'])
    for node in ['profile_controls','profile_gex']:job(node,'downstream_profiles.py',['profiles'])
    job('moran','downstream_moran.py',['profiles'],False)
    job('measurement','downstream_measurement.py',['recurrence'],False)
    job('functional','downstream_functional.py',['recurrence'],False)
    job('dwui_pool','downstream_cohorts.py',['state_vectors'])
    job('union_pool','downstream_cohorts.py',['dwui_pool','recurrence'])
    job('matched_pms','downstream_cohorts.py',['recurrence','well_inference'])
    intake=PROJECT/'notes/2026-10-03-polyapipe-rerun/run/intake.native.json'
    config=dict(project_root=str(dest/'project'),native_run=str(REMOTE/'run'),frozen_gex=str(REMOTE/'frozen_gex'),
        analysis_root=str(REMOTE/'statistics_v1'),run_id='polyapipe_native_20261003_v1',input_fingerprint=digest(intake),
        resources=resource_manifest,source_overrides=overrides,jobs=jobs,count_workers=4,
        code_sha256={str(dest/'code'/p.name):digest(p) for p in code.glob('*.py')})
    (a.out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    (a.out/'package.json').write_text(json.dumps(dict(scope='bound statistics; manuscript pipeline remains incomplete',
        nodes=list(jobs),files={str(p.relative_to(a.out)):digest(p) for p in a.out.rglob('*') if p.is_file()}),indent=2)+'\n')
    print(json.dumps(dict(package=str(a.out),remote_destination=str(dest),nodes=len(jobs),bytes=sum(p.stat().st_size for p in a.out.rglob('*') if p.is_file()))))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);main(p.parse_args())
