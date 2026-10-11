"""Supplement tests preserve denominators, site identity and source isolation."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from downstream_supplement_statistics import dual_points,relationship_tables,select_partial_sites


class SupplementStatisticsTests(unittest.TestCase):
    def test_dual_display_keeps_rejected_measurable_candidates(self):
        rec=pd.DataFrame([dict(species='mouse',gene_id=g,metric='WPP',threshold=.02,series_matching=True,contrast='P_M',direction='positive') for g in ['yes','no','missing']])
        primary=pd.DataFrame([dict(species='mouse',setting='supported50_p08',gene_id=g,route=route,DWUI_PS=.03,DWUI_PM=.04 if g=='yes' else -.04)
            for g in ['yes','no','missing'] for route in (['CIP'] if g=='missing' else ['CIP','WOT'])])
        points,audit=dual_points(primary,rec)
        self.assertEqual(set(points.gene_id),{'yes','no'})
        self.assertEqual(points.groupby('gene_id').dual_qualified_gene.first().to_dict(),{'no':False,'yes':True})
        self.assertFalse(audit.set_index('gene_id').loc['missing','matched_measurable'])

    def test_relationship_union_is_not_sum_of_contrasts_or_gene_symbols(self):
        rows=[dict(species='mouse',category='miRNA',gene_id=g,factor=f,contrast=c,exact_partial_recurrence=True)
            for g,f,c in [('g1','seed','SP'),('g1','seed','MP'),('g2','seed','MP'),('g1','other','MP')]]
        members,counts=relationship_tables(pd.DataFrame(rows))
        q=counts.query("species=='mouse' and category=='miRNA'").iloc[0]
        self.assertEqual((q.relationships,q.genes),(3,2))
        self.assertEqual(members.query("gene_id=='g1' and factor=='seed'").contrasts.iloc[0],'MP;SP')
        self.assertEqual(counts.query("species=='human'").relationships.sum(),0)

    def fixture(self):
        partial=pd.DataFrame([dict(species='mouse',gene_id='new',gene_symbol_a='NewGene',factor='newseed',category='miRNA',contrast='MP',
            exact_partial_recurrence=True,shared_recurrent_changed_sites='chr1:10-17',shared_retained_sites='chr1:20-27')])
        exact=[];features=[]
        for role,start,mask in [('changed',10,[0,1]),('retained',20,[1,1])]:
            exact.append(dict(species='mouse',gene_id='new',factor='newseed',category='miRNA',contrast='MP',chrom='chr1',strand='+',start0=start,end0=start+7,
                recurrent_strict=role=='changed',retained_a=role=='retained',retained_b=role=='retained',route_a='CIP',route_b='WOT',
                mask_a=json.dumps([[x,x] for x in mask]),mask_b=json.dumps([[x,x] for x in mask])))
            for route in ['CIP','WOT']:features.append(dict(species='mouse',gene_id='new',route=route,kind='site',label='newseed',start0=start,end0=start+7,
                feature_id=f'{route}|{start}',mask=json.dumps(mask)))
        selection=pd.DataFrame([dict(role='partial_retention',species='mouse',rank=1,gene_symbol='NewGene')])
        return partial,pd.DataFrame(exact),pd.DataFrame(features),selection

    def test_site_pair_requires_both_routes_and_keeps_exact_identity(self):
        partial,exact,features,selection=self.fixture()
        sites,status=select_partial_sites(partial,exact,features,selection)
        self.assertEqual(len(sites),4);self.assertEqual(set(sites.gene_id),{'new'})
        self.assertEqual(set(sites.site_role),{'changed','retained'})
        missing=features[~(features.route.eq('WOT')&features.start0.eq(20))]
        sites,status=select_partial_sites(partial,exact,missing,selection)
        self.assertTrue(sites.empty)
        self.assertIn('no_current_exported_exact_site_pair',status.status.tolist())

    def test_site_pair_rejects_a_mask_that_uses_a_different_PAS_order(self):
        partial,exact,features,selection=self.fixture();features.loc[0,'mask']='[1, 0]'
        with self.assertRaises(AssertionError):select_partial_sites(partial,exact,features,selection)

    def test_empty_partial_cohort_does_not_restore_an_old_example(self):
        partial,exact,features,selection=self.fixture()
        selected,status=select_partial_sites(partial.iloc[:0],exact,features,selection.iloc[:0])
        self.assertTrue(selected.empty);self.assertEqual(len(status),2)

    def test_late_selected_site_can_be_exported_from_current_models(self):
        partial,exact,features,selection=self.fixture()
        selected,status=select_partial_sites(partial,exact,features.iloc[:0],selection,allow_new_features=True)
        self.assertEqual(len(selected),4)
        self.assertEqual(set(selected.factor),{'newseed'})
        self.assertTrue(selected.feature_id.str.contains('new').all())

    def test_completed_run_extension_preserves_main_contract_and_original_bytes(self):
        import shutil
        from test_result_records import bundle
        from pipeline import digest,read,write
        from result_records import ResultBundle,export_records
        from complete_supplement_bundle import complete
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source,identity,_=bundle(root,'base',-.4)
            manifest=read(source/'table-manifest.json');manifest['node_receipts']={'global':'fixture'}
            write(source/'table-manifest.json',manifest)
            panels=read(source/'plot-contract.json')['panels']
            export_records(source,identity,panels,manifest['tables'],source/'contracts')
            marker=read(source/'bundle-complete.json');marker.update(table_manifest_sha256=digest(source/'table-manifest.json'),result_records_sha256=digest(source/'result-records.json'))
            write(source/'bundle-complete.json',marker)
            (source/'contracts/required-statistics-tables.tsv').write_text('path\trequired_columns\nglobal/effects.tsv\troute;estimate;n\n')
            shutil.copytree(source/'contracts',source/'manuscript-contracts')
            config=root/'config.json';write(config,dict(identity,project_root=str(root),native_run='UNUSED',frozen_gex='UNUSED'))
            contracts=root/'new';contracts.mkdir()
            (contracts/'supplement-producer-bindings.tsv').write_text('panel\tproducers\nD11D\tsupplement_statistics\n')
            (contracts/'supplement-statistics-tables.tsv').write_text('path\trequired_columns\nsupplement_statistics/points.tsv\tgene_id;qualified\n')
            (contracts/'supplement-data-interface.tsv').write_text('component_id\tpanel\ttable\tread_columns\tselection\tunit_or_uncertainty\nD11D:supplement_statistics/points.tsv\tD11D\tsupplement_statistics/points.tsv\tgene_id;qualified\tall measurable\tgenes\n')
            # A concurrent main-interface edit must not leak into the completed source run.
            (contracts/'panel-data-interface.tsv').write_text('INVALID MAIN CONTRACT MUST NOT BE READ')
            before={p:digest(p) for p in source.rglob('*') if p.is_file()}
            def compute(ctx):
                ctx.table([dict(gene_id='fresh',qualified=False)],'points.tsv');ctx.finish(dict(current_input=True))
            with patch('complete_supplement_bundle.supplement',side_effect=compute),patch('complete_supplement_bundle.panel_contracts',return_value=([
                dict(panel='F1D',nodes=['global']),dict(panel='D11D',nodes=['supplement_statistics'])],[])):
                result=complete(source,config,contracts,root/'extended')
            reader=ResultBundle(result,identity)
            table=reader.table('D11D','D11D:supplement_statistics/points.tsv')
            self.assertIn('fresh',table.read_text());self.assertTrue(all(digest(p)==h for p,h in before.items()))
            self.assertEqual(digest(result/'manuscript-contracts/panel-data-interface.tsv'),digest(source/'manuscript-contracts/panel-data-interface.tsv'))
            write(config,dict(read(config),input_fingerprint='different-source'))
            with self.assertRaisesRegex(ValueError,'different PAS input'):complete(source,config,contracts,root/'invalid')
            self.assertFalse((root/'invalid').exists())

    def test_gap_tables_are_gated_and_bound_to_the_new_DAG(self):
        import csv,sys
        from prepare_pipeline import make_jobs
        from export_statistics_bundle import panel_contracts,require_numerical_products
        root=Path(__file__).resolve().parents[2];folder=root/'notes/2026-10-05-apadetect-analysis-intake'
        jobs=make_jobs(json.loads((folder/'pipeline.environment.json').read_text()),root/'analysis_code/pas_rerun',sys.executable)
        self.assertIn('supplement_statistics',jobs['bundle']['depends_on'])
        panels,_=panel_contracts(folder/'contracts',jobs);panels={p['panel']:p for p in panels}
        read=lambda name:list(csv.DictReader((folder/'contracts'/name).open(),delimiter='\t'))
        rules=read('required-statistics-tables.tsv')+read('supplement-statistics-tables.tsv')
        schemas={}
        for rule in rules:schemas.setdefault(rule['path'],set()).update(rule['required_columns'].split(';'))
        for r in read('supplement-data-interface.tsv'):
            self.assertIn(r['table'].split('/')[0],panels[r['panel']]['nodes'])
            self.assertTrue(set(r['read_columns'].split(';'))<=schemas[r['table']])
        rows=[dict(path='tables/'+path,columns=sorted(cols)) for path,cols in schemas.items()]
        require_numerical_products(folder/'contracts',rows)
        with self.assertRaisesRegex(ValueError,'Missing required statistical product'):
            require_numerical_products(folder/'contracts',[r for r in rows if not r['path'].endswith('dual-reference-all-measurable.tsv.gz')])

    def test_completed_old_count_waits_for_repair_without_a_sticky_failure(self):
        from continue_supplement_series import run
        from pipeline import write,read
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);series=root/'series';contracts=root/'contracts';contracts.mkdir()
            bundle=root/'bundle';config=root/'config.json';analysis=root/'statistics'
            write(series/'source-runs.json',dict(routes=[dict(source_id='base',status='statistics_complete',bundle=str(bundle),config=str(config))]))
            write(config,dict(analysis_root=str(analysis)))
            write(analysis/'pipeline-status.json',dict(phase='statistics_complete'))
            write(bundle/'bundle-complete.json',dict(status='complete'))
            write(bundle/'tables/count_inference/audit.json',dict(checks=dict(fresh_comparison_BH=True)))
            with patch('continue_supplement_series.subprocess.run') as execute:
                self.assertEqual(run(series,contracts,root/'continuation',once=True),0)
                execute.assert_not_called()
            state=read(root/'continuation/status.json')
            self.assertEqual(state['routes']['base']['status'],'waiting_for_corrected_count_inference')
            self.assertEqual(state['status'],'waiting_for_source_bundles')

    def test_continuation_waits_for_a_complete_bundle_and_preserves_failed_status(self):
        from continue_supplement_series import run
        from pipeline import write,read
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);series=root/'series';contracts=root/'contracts';contracts.mkdir()
            (contracts/'supplement-statistics-tables.tsv').write_text('fixed-contract')
            rows=[dict(source_id=k,status=s,bundle=str(series/k/'bundle'),config=str(series/k/'config.json')) for k,s in [('base','running'),('asnr','failed')]]
            write(series/'source-runs.json',dict(status='running',routes=rows))
            with patch('continue_supplement_series.subprocess.run') as execute:
                self.assertEqual(run(series,contracts,root/'continuation',once=True),0)
                execute.assert_not_called()
            state=read(root/'continuation/status.json')
            self.assertEqual(state['routes']['base']['status'],'waiting_for_source_bundle')
            self.assertEqual(state['routes']['asnr']['status'],'upstream_failed')
            (contracts/'supplement-statistics-tables.tsv').write_text('changed-contract')
            with self.assertRaisesRegex(ValueError,'contracts changed'):run(series,contracts,root/'continuation',once=True)


if __name__=='__main__':unittest.main()
