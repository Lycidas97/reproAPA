"""Bind one prepared PAS input to a fresh, complete manuscript statistics DAG."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys
from pipeline import digest,read,write,require,topological

HISTORICAL={'recount','consensus','caller_comparison','reference_comparison'}
EXTENSIONS={
 'supplement_statistics':('downstream_supplement_statistics.py',['panel_statistics','profiles','shape_ortholog','candidate_screen','matched_pms','recurrence','functional','cross_continuous','cross_overlap','site_statistics','examples','geometry','axis_extensions']),
 'direction_difference':('downstream_direction.py',['cross_sharing']),
 'external_mirna':('downstream_mirna.py',['motif_background','site_scan']),
 'shape_ortholog':('downstream_shape_ortholog.py',['profiles','cross_overlap']),
 'naive_state':('downstream_naive.py',['tandem','global']),
 'same_pas':('downstream_same_pas.py',['matched_pms','cis']),
 'candidate_screen':('downstream_candidate_screen.py',['matched_pms','recurrence','site_statistics','well_inference','state_vectors']),
 'axis_extensions':('downstream_axis_extensions.py',['cis','global','profiles','apa_gex','well_inference','pool_inference','intermediate','external_pas','longread','ago2']),
 'geometry':('downstream_geometry.py',['candidate_screen','state_vectors','intervals','profiles','global','cross_trace']),
 'panel_statistics':('downstream_panel_statistics.py',['candidate_screen','geometry','naive_state','shape_ortholog','same_pas','axis_extensions','measurement','profile_controls','site_statistics','weighted_retention','cis','cross_sharing_range10','figure_scores','shared_projection','rbp_expression','moran','functional','profile_gex','count_inference','pas_qc','metacell','motif_controls','original_pool_sites','coordinate_states','examples','external_mirna','direction_difference']),
}


def primary_reference(manifest,annotation):
    variants={manifest.get('reference_variant'),annotation.get('support_reference_variant')}-{None}
    require(variants=={'v2_primary'},'Current manuscript requires an unambiguous v2 integrated reference')


def make_jobs(environment,code,python):
    jobs={n:dict(script=str(code/Path(j['script']).name),arguments=j.get('arguments',[]),
                 depends_on=list(j['depends_on']),python=python) for n,j in environment['jobs'].items() if n not in HISTORICAL}
    jobs['core']=dict(script=str(code/'downstream_core.py'),arguments=[],depends_on=[],python=python)
    jobs['pas_qc']['script']=str(code/'downstream_catalog_qc.py')
    jobs['cross_sharing_range10']=dict(script=str(code/'complete_cross_sharing_range10.py'),arguments=[],depends_on=['recurrence'],python=python)
    for n,(script,deps) in EXTENSIONS.items():jobs[n]=dict(script=str(code/script),arguments=[],depends_on=deps,python=python)
    jobs['examples']['depends_on']=sorted(set(jobs['examples']['depends_on'])|{'candidate_screen'})
    for n,j in jobs.items():
        if n!='core':j['depends_on']=sorted(set(j['depends_on'])|{'core'})
    jobs['bundle']=dict(script=str(code/'export_statistics_bundle.py'),arguments=[],depends_on=sorted(jobs),python=python)
    topological(jobs)
    for j in jobs.values():require(Path(j['script']).is_file(),f"Missing implementation {j['script']}")
    return jobs


def prepare(environment_path,native,out,plan_only=False,python=None,core_workers=None):
    env=read(environment_path);source=Path(native).resolve();out=Path(out).resolve()
    core_workers=env.get('core_workers',2) if core_workers is None else core_workers
    require(isinstance(core_workers,int) and core_workers>0,'core_workers must be a positive integer')
    require(not out.exists(),'Use a new pipeline package directory; existing runs are immutable')
    require(not out.is_relative_to(source),'Pipeline output must be outside the prepared input')
    require(not any(p in out.parts for p in ['Fig','main','supplement','production-v15']),'Never write a candidate run into production')
    out.mkdir(parents=True);code=out/'code';code.mkdir()
    for p in Path(__file__).parent.glob('*.py'):
        if p.name.startswith(('render_','test_','prepare_statistics')):continue
        shutil.copy2(p,code/p.name)
    contracts=out/'contracts';shutil.copytree(Path(environment_path).resolve().parent/'contracts',contracts)
    jobs=make_jobs(env,code,python or env['python'])
    config=dict(project_root=env['project_root'],native_run=str(out/'statistics/core/native'),prepared_native_source=str(source),
        frozen_gex=env['frozen_gex'],analysis_root=str(out/'statistics'),reference_variant='v2_primary',
        count_workers=env.get('count_workers',4),core_workers=core_workers,example_candidates_per_role=3,
        resources={k:v for k,v in env['resources'].items() if not k.endswith('_resume')},
        source_overrides=env['source_overrides'],jobs=jobs,contracts_root=str(contracts),
        input_files={},code_sha256={str(p):digest(p) for p in sorted(code.glob('*.py'))},
        policy=dict(analysis_routes='parallel_primary_sources',reference='v2_only_v3_excluded_by_user',recount=False,production_replacement=False,
                    statistics_and_rendering='separate_commands',historical_audits='not_current_results'),
        run_id='UNBOUND_PLAN',input_fingerprint='UNBOUND_PLAN',prepared=False)
    for p in sorted(contracts.rglob('*')):
        if p.is_file():config['code_sha256'][str(p)]=digest(p)
    if not plan_only:
        manifest=read(source/'intake.native.json')
        require((source/'annotation.complete.json').is_file(),'Full annotation receipt is required')
        primary_reference(manifest,read(source/'annotation.complete.json'))
        require(len(manifest['studies'])==7 and sum(r['libraries'] for r in manifest['studies'])==161,'All 7 studies / 161 frozen libraries required')
        require(manifest.get('status')!='prepared_pilot','Pilot cannot be a full manuscript input')
        inputs={}
        for name in ['intake.native.json','annotation.complete.json','frozen_native_membership_audit.json','frozen_native_membership_audit.tsv']:
            p=source/name;inputs[str(p)]=digest(p)
        # Every route freezes its full matrix axis plus its own eligibility flags.
        for row in manifest['studies']:
            sid=row['study_id']
            for name in ['features.tsv.gz','metacells.tsv.gz','selected_metacells.tsv.gz','sampling_design.tsv','endpoint_library_design.tsv','metacell_by_pas_counts.npz']:
                p=source/'studies'/sid/name;h=digest(p);inputs[str(p)]=h
                if 'outputs' in manifest:require(h==manifest['outputs'][sid][name],f'Prepared input changed: {p}')
            for name in ['metacells.tsv.gz','genes.tsv.gz','metacell_by_gene_counts.npz']:
                p=Path(env['frozen_gex'])/sid/name;inputs[str(p)]=digest(p)
        for key,item in {**config['resources'],**config['source_overrides']}.items():
            p=Path(item['path']);h=digest(p);require(h==item['sha256'],f'Pinned resource changed: {key}')
            inputs[str(p)]=h
        inputs[str(Path(environment_path).resolve())]=digest(environment_path)
        config.update(prepared=True,run_id=manifest['run_id'],analysis_route=manifest.get('analysis_route','base'),
            caller=manifest['caller'],input_files=inputs,input_fingerprint=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest())
    write(out/'pipeline.config.json',config)
    write(out/'plan.json',dict(prepared=config['prepared'],order=topological(jobs),statistics_nodes=len(jobs),
        historical_only=sorted(HISTORICAL),production_replacement=False,rendering_performed=False))
    return config


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--environment',type=Path,required=True)
    p.add_argument('--native',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--plan-only',action='store_true');p.add_argument('--python')
    a=p.parse_args();c=prepare(a.environment,a.native,a.out,a.plan_only,a.python)
    print(json.dumps(dict(config=str(a.out/'pipeline.config.json'),nodes=len(c['jobs']),prepared=c['prepared'])))
