"""Freeze additional modules, retaining the first running package unchanged."""
import argparse
import json
from pathlib import Path
import shutil

from polyapipe_prepare import digest
from prepare_statistics_package import PROJECT,PROJECTS,REMOTE

F5='analysis_code/release_2026-09-24/sources/F5/'


def main(args):
    out=args.out;out.mkdir(parents=True,exist_ok=False);dest=REMOTE/out.name
    for name in ['code','legacy','resources']:(out/name).mkdir()
    config=json.loads(args.base.read_text())
    for source in Path(__file__).parent.glob('*.py'):
        if source.name!='downstream_common.py':shutil.copy2(source,out/'code'/source.name)
    sources=['notes/2026-10-01-mpro-biology/legacy_structure_consensus.py','notes/2026-10-01-mpro-biology/legacy_prepare_pilot.py',
        'notes/2026-09-29-core-sensitivity/figure5_paired_pm.py','notes/2026-09-29-core-sensitivity/figure4_matched_overlap.py',
        'notes/2026-10-02-null-model-controls/cross_species.py','notes/2026-10-02-null-model-controls/motif.py',
        'notes/2026-10-02-reviewer-followup/motif_background.py','notes/2026-10-02-null-result-diagnostics/motif_sensitivity.py',
        'notes/2026-10-02-rbp-expression/analyze.py',
        'notes/2026-10-01-figure5-union-update/resources/mouse_time_aware_integration_pilot_v1/regulatory_catalog_v1/scripts/build_catalog.py',
        F5+'annotation-pas-global-v1/migrate.py',F5+'site-recurrence-v1/figure5-publication-v14r4/prepare_annotated_tracks.py',
        F5+'full-cohort-v2/sources/execution-scripts/core.py',F5+'site-recurrence-v1/analyze.py',F5+'site-recurrence-v1/consistency.py']
    for relative in sources:
        path=PROJECT/relative;name=digest(path)[:12]+'_'+path.name;shutil.copy2(path,out/'legacy'/name)
        config['source_overrides'][relative]=dict(path=str(dest/'legacy'/name),sha256=digest(path))
    union=PROJECT/'notes/2026-10-01-mpro-biology';res=PROJECT/'notes/2026-10-01-figure5-union-update/resources'
    mouse=res/'mouse_time_aware_integration_pilot_v1';human=res/'human_time_aware_integration_v1'
    mapping=dict(regulatory_bundles=union/'remote_resource_paths.json',
        mouse_terminal_utrs=PROJECTS/'2026-09-02_mouse-time-aware-integration-pilot/data/tandem-apa-v1/all_reference_terminal_utr_candidates.tsv.gz',
        human_terminal_utrs=union/'human_terminal_utrs.tsv.gz',mouse_catalog_utrs=union/'mouse_catalog_utrs.tsv.gz',human_catalog_utrs=union/'human_catalog_utrs.tsv.gz',
        mouse_utr_fasta=res/'mouse_sequences.fa.gz',human_utr_fasta=res/'human_sequences.fa.gz',
        mouse_mirna_annotations=union/'mouse_mirna_annotations.tsv',human_mirna_annotations=res/'regulatory_resource_revision_v1/resources/human_mirna_annotations.tsv',
        mouse_motif_provenance=union/'mouse_reference_motif_provenance.tsv.gz',human_motif_provenance=res/'human_motif_provenance.json',
        human_direct_meme=human/'regulatory_query_v1/resources/human_direct_scannable.meme',mouse_transfer_meme=mouse/'regulatory_query_v2/resources/transferred_scannable.meme',
        mouse_motif_background=mouse/'regulatory_catalog_v1/background.txt',human_motif_background=res/'human_background.txt')
    for key,path in mapping.items():
        name=key+'__'+path.name;shutil.copy2(path,out/'resources'/name)
        config['resources'][key]=dict(path=str(dest/'resources'/name),sha256=digest(path))
    config['resources'].update(json.loads(args.remote_resources.read_text()))
    jobs={}
    def job(node,script,deps,argument=True):jobs[node]=dict(script=str(dest/'code'/script),arguments=['--node',node] if argument else [],depends_on=deps)
    job('cross_overlap','downstream_cross.py',['state_vectors']);job('cross_controls','downstream_cross.py',['cross_overlap'])
    job('intervals','downstream_regulatory.py',['union_pool']);job('site_scan','downstream_regulatory.py',['intervals'])
    job('site_statistics','downstream_sites.py',['site_scan']);job('original_pool_sites','downstream_sites.py',['site_scan','dwui_pool'])
    job('motif_background','downstream_motifs.py',['site_statistics']);job('motif_controls','downstream_motifs.py',['motif_background'])
    job('weighted_retention','downstream_motifs.py',['motif_controls']);job('rbp_expression','downstream_rbp.py',['site_statistics'],False)
    config['jobs']=jobs;config['upstream_dispatch_status']=[str(Path(config['analysis_root'])/'dispatch-status.json')]
    config['code_sha256']={str(dest/'code'/p.name):digest(p) for p in (out/'code').glob('*.py')}
    (out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    (out/'package.json').write_text(json.dumps(dict(scope='additional statistical modules; full manuscript rerun remains incomplete',nodes=list(jobs),
        files={str(p.relative_to(out)):digest(p) for p in out.rglob('*') if p.is_file()}),indent=2)+'\n')
    print(json.dumps(dict(package=str(out),remote_destination=str(dest),nodes=len(jobs),bytes=sum(p.stat().st_size for p in out.rglob('*') if p.is_file()))))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--remote-resources',type=Path,required=True);main(p.parse_args())
