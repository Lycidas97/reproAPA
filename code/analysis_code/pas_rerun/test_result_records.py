"""Intent: input switching cannot mix sources, change estimators or silently use old plots."""
import csv
import json
from pathlib import Path
import sys
import tempfile
import threading
import unittest
from unittest.mock import patch
from pipeline import digest,write,read
from result_records import export_records,ResultBundle,collect_sources,source_identity
from render_sources import plan,execute as render
from run_sources import execute_sources

ROOT=Path(__file__).resolve().parents[2]


def bundle(root,name,value,caller='APADETECT',route='base',builder=None,ready=True):
    out=root/name;out.mkdir();contracts=out/'contracts';contracts.mkdir()
    identity=source_identity(dict(run_id=name,input_fingerprint='input-'+name,caller=caller,
                                 analysis_route=route,reference_variant='v2_primary'))
    table=out/'tables/global/effects.tsv';table.parent.mkdir(parents=True)
    table.write_text(f'route\testimate\tn\nCIP\t{value}\t17\n')
    columns=['route','estimate','n'];tables=[dict(path='tables/global/effects.tsv',sha256=digest(table),columns=columns,node='global')]
    panels=[dict(panel='F1D',tables=[tables[0]['path']])]
    write(out/'plot-contract.json',dict(identity,panels=panels));write(out/'table-manifest.json',dict(tables=tables))
    (contracts/'panel-data-interface.tsv').write_text('component_id\tpanel\ttable\tread_columns\tselection\tunit_or_uncertainty\n'
        'F1D:global/effects.tsv\tF1D\tglobal/effects.tsv\troute;estimate;n\tCIP\tgenes\n')
    write(contracts/'figure-renderers.json',dict(figures=[] if builder is None else [dict(figure='F1',script=builder.name,
        sha256=digest(builder),panels=['F1D'],interface='pas-panel-result-records/1' if ready else 'legacy_direct_tables',
        status='ready' if ready else 'requires_records_adapter')]))
    record=export_records(out,identity,panels,tables,contracts)
    write(out/'bundle-complete.json',dict(status='statistics_bundle_complete',run_id=name,
        plot_contract_sha256=digest(out/'plot-contract.json'),table_manifest_sha256=digest(out/'table-manifest.json'),
        result_records_sha256=digest(out/'result-records.json')))
    return out,identity,record


class ResultRecordTests(unittest.TestCase):
    def test_same_schema_and_component_keys_different_input_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);a,ia,ra=bundle(root,'base',-.4);b,ib,rb=bundle(root,'asnr',.2,route='reject_original_point')
            c,ic,rc=bundle(root,'polyapipe',.8,caller='polyApipe')
            self.assertEqual(len({r['statistics_schema_id'] for r in [ra,rb,rc]}),1)
            self.assertEqual(len({r['components'][0]['component_id'] for r in [ra,rb,rc]}),1)
            self.assertEqual(len({r['components'][0]['record_id'] for r in [ra,rb,rc]}),3)
            for path,identity,value in [(a,ia,'-0.4'),(b,ib,'0.2'),(c,ic,'0.8')]:
                reader=ResultBundle(path,identity);p=reader.table('F1D','F1D:global/effects.tsv')
                with p.open() as f:row=next(csv.DictReader(f,delimiter='\t'))
                self.assertEqual(row['estimate'],value);self.assertEqual(row['route'],'CIP')
                self.assertFalse(reader.provenance()['statistics_recomputed'])

    def test_wrong_source_changed_bytes_and_wrong_panel_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);a,ia,_=bundle(root,'a',-.4);_,ib,_=bundle(root,'b',.2)
            with self.assertRaisesRegex(ValueError,'Requested source'):ResultBundle(a,ib)
            reader=ResultBundle(a,ia)
            with self.assertRaisesRegex(ValueError,'another panel'):reader.table('F2A','F1D:global/effects.tsv')
            (a/'tables/global/effects.tsv').write_text('route\testimate\tn\nCIP\t999\t17\n')
            with self.assertRaisesRegex(ValueError,'changed'):reader.table('F1D','F1D:global/effects.tsv')

    def test_legacy_builder_preflight_cannot_generate_old_figure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);script=root/'old.py';script.write_text('raise RuntimeError("must never execute")\n')
            a,_,_=bundle(root,'a',1.,builder=script,ready=False)
            series=root/'series.json';collect_sources([dict(source_id='base',bundle=a)],series)
            p=plan(series,root,root/'candidate')
            self.assertEqual(p['status'],'blocked_by_unadapted_builders')
            with self.assertRaisesRegex(ValueError,'adapters are incomplete'):render(p,root/'candidate',1)
            self.assertFalse((root/'candidate/base/F1').exists())

    def test_two_sources_dispatch_same_builder_to_isolated_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);script=root/'figure.py'
            script.write_text('''import argparse,json,csv,sys
from pathlib import Path
sys.path.insert(0, CODE)
from result_records import ResultBundle
from pipeline import write,digest
p=argparse.ArgumentParser();p.add_argument('--records');p.add_argument('--out');a=p.parse_args()
r=ResultBundle(Path(a.records).parent);f=r.table('F1D','F1D:global/effects.tsv')
with f.open() as stream:value=next(csv.DictReader(stream,delimiter='\\t'))['estimate']
out=Path(a.out)
artifacts=[]
for name in ['fixture.pdf','fixture.png']:
 file=out/name;file.write_text('SYNTHETIC DISPATCH TEST ONLY '+value);artifacts.append(dict(path=name,sha256=digest(file)))
receipt=r.provenance();receipt.update(figure='F1',panels=['F1D'],artifacts=artifacts)
write(out/'render-receipt.json',receipt)
'''.replace('CODE',repr(str(ROOT/'analysis_code/pas_rerun'))))
            a,_,_=bundle(root,'a',-1.,builder=script);b,_,_=bundle(root,'b',2.,route='ASNR',builder=script)
            series=root/'series.json';collect_sources([dict(source_id='base',bundle=a),dict(source_id='asnr',bundle=b)],series)
            before={str(p):digest(p) for path in [a,b] for p in path.rglob('*') if p.is_file()}
            p=plan(series,root,root/'candidate');render(p,root/'candidate',2)
            self.assertIn('-1.0',(root/'candidate/base/F1/fixture.pdf').read_text())
            self.assertIn('2.0',(root/'candidate/asnr/F1/fixture.pdf').read_text())
            self.assertTrue(all(digest(Path(p))==h for p,h in before.items()))
            self.assertEqual(read(root/'candidate/render-status.json')['status'],'candidate_render_complete')

    def test_one_failed_source_does_not_cancel_independent_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);barrier=threading.Barrier(2);called=[]
            def execute(config,workers):
                barrier.wait(timeout=5);called.append(config)
                if config.endswith('bad.json'):raise ValueError('bad input fixture')
            runs=[dict(source_id=n,route=n,config=str(root/(n+'.json')),bundle=str(root/n),identity={}) for n in ['good','bad']]
            with patch('run_sources.execute',execute),patch('run_sources.collect_sources') as collect:
                with self.assertRaisesRegex(ValueError,'Failed sources'):execute_sources(runs,root,1,2)
                self.assertEqual(collect.call_args.args[0][0]['source_id'],'good')
            self.assertEqual(len(called),2)
            state=read(root/'source-runs.json');self.assertEqual(state['status'],'partial_failure')
            self.assertEqual({r['source_id']:r['status'] for r in state['routes']},{'good':'statistics_complete','bad':'failed'})

    def test_same_bundle_cannot_be_disguised_as_another_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);a,_,_=bundle(root,'a',1.)
            with self.assertRaisesRegex(ValueError,'same source bundle'):
                collect_sources([dict(source_id='a',bundle=a),dict(source_id='b',bundle=a)],root/'series.json')

    def test_result_index_failure_is_not_reported_as_complete(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            runs=[dict(source_id='a',config=str(root/'a.json'),bundle=str(root/'a'),identity={})]
            with patch('run_sources.execute'),patch('run_sources.collect_sources',side_effect=ValueError('incompatible schema')):
                with self.assertRaisesRegex(ValueError,'incompatible schema'):execute_sources(runs,root,1)
            self.assertEqual(read(root/'source-runs.json')['status'],'result_index_failed')

    def test_ready_source_starts_statistics_without_waiting_for_other_preparations(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);started=threading.Event();events=[]
            runs=[dict(source_id=n,config=str(root/(n+'.json')),bundle=str(root/n),identity={}) for n in ['ready','slow']]
            def prepare(row):
                if row['source_id']=='slow':
                    self.assertTrue(started.wait(timeout=5),'A ready source must not wait for every preparation')
                events.append('prepared:'+row['source_id']);return row
            def execute(config,workers):
                if config.endswith('ready.json'):started.set()
            with patch('run_sources.execute',execute),patch('run_sources.collect_sources'):
                execute_sources(runs,root,1,2,prepare)
            self.assertEqual(events,['prepared:ready','prepared:slow'])

    def test_preparation_failure_isolated_and_prepare_only_never_computes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            runs=[dict(source_id=n,config=str(root/(n+'.json'))) for n in ['good','bad']]
            def prepare(row):
                if row['source_id']=='bad':raise ValueError('invalid source')
                return row
            with patch('run_sources.execute') as execute,patch('run_sources.collect_sources') as collect:
                with self.assertRaisesRegex(ValueError,'Failed sources'):execute_sources(runs,root,1,2,prepare,True)
                execute.assert_not_called();collect.assert_not_called()
            state=read(root/'source-runs.json')
            self.assertEqual({r['source_id']:r['status'] for r in state['routes']},{'good':'prepared','bad':'failed'})

    def test_interface_upgrade_preserves_computed_values_and_original_bundle(self):
        from upgrade_result_interface import upgrade
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);a,identity,_=bundle(root,'v8',-.4)
            manifest=read(a/'table-manifest.json');manifest['node_receipts']={'global':'verified-fixture'}
            write(a/'table-manifest.json',manifest)
            marker=read(a/'bundle-complete.json');marker['table_manifest_sha256']=digest(a/'table-manifest.json')
            marker.pop('result_records_sha256');write(a/'bundle-complete.json',marker)
            (a/'result-records.json').unlink();(a/'result-records.tsv').unlink()
            (a/'contracts/required-statistics-tables.tsv').write_text('path\trequired_columns\nglobal/effects.tsv\troute;estimate;n\n')
            before={str(p):digest(p) for p in a.rglob('*') if p.is_file()}
            panels=[dict(panel='F1D',nodes=['global'])]
            with patch('upgrade_result_interface.panel_contracts',return_value=(panels,[])):
                upgraded=upgrade(a,a/'contracts',root/'upgraded')
            reader=ResultBundle(upgraded,identity)
            table=reader.table('F1D','F1D:global/effects.tsv')
            self.assertEqual(digest(table),digest(a/'tables/global/effects.tsv'))
            self.assertTrue(all(digest(Path(p))==h for p,h in before.items()))
            self.assertTrue(read(upgraded/'bundle-complete.json')['interface_only_upgrade'])
            self.assertFalse(read(upgraded/'interface-parent.json')['source_statistics_recomputed'])
            # Missing numerical products must stop before creating a partial interface package.
            (a/'contracts/required-statistics-tables.tsv').write_text('path\trequired_columns\nmissing/table.tsv\testimate\n')
            with self.assertRaisesRegex(ValueError,'Missing required statistical product'):
                upgrade(a,a/'contracts',root/'missing')
            self.assertFalse((root/'missing').exists())

    def test_current_main_component_contract_is_fully_bound(self):
        from manuscript_contract import MAIN
        from export_statistics_bundle import require_numerical_products
        folder=ROOT/'notes/2026-10-05-apadetect-analysis-intake/contracts'
        with (folder/'panel-data-interface.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        self.assertEqual(len(rows),82);self.assertEqual(len({r['panel'] for r in rows}),27)   # +F3C recurrence_members (records adapter)
        for r in rows:self.assertTrue(r['table'].split('/')[0] in MAIN[r['panel']] or r['table'].startswith('frozen/'),r)
        with (folder/'required-statistics-tables.tsv').open() as f:required=list(csv.DictReader(f,delimiter='\t'))
        self.assertEqual(len(required),278)
        by_path={r['path']:set(r['required_columns'].split(';')) for r in required}
        for r in rows:
            self.assertIn(r['table'],by_path)
            self.assertTrue(set(r['read_columns'].split(';'))<=by_path[r['table']],r['table'])


if __name__=='__main__':unittest.main()
