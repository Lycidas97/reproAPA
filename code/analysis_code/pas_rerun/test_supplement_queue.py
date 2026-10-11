"""Queue handovers must not lose routes, reuse stale paths or accept early bundles."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch,Mock

from pipeline import read,write
from continue_supplement_series import run,source_rows,source_bindings


class SupplementQueueTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name).resolve();self.series=self.root/'families';self.out=self.root/'completion'
        self.contracts=self.root/'contracts';self.contracts.mkdir()
        (self.contracts/'supplement-statistics-tables.tsv').write_text('fixed contract')
        groups={f:[dict(route=k,family=f) for k in keys] for f,keys in [
            ('uniform_count_base',['base']),('ASNR',['asnr_old']+[f'asnr_{i}' for i in range(5)]),
            ('ASNR_plus_terminal',[f'terminal_{i}' for i in range(18)])]}
        adopted={k:dict(config=str(self.root/'old/runs'/k/'pipeline.config.json'),
                       bundle=str(self.root/'old/runs'/k/'statistics/bundle')) for k in ['base','asnr_old']}
        write(self.series/'family-plan.json',dict(groups=groups,adoption=dict(routes=adopted)))
        for f,items in groups.items():
            write(self.series/'family-status'/f'{f}.json',dict(routes=[dict(route=r['route'],status='running') for r in items]))

    def change_status(self,family,route,status,**extra):
        p=self.series/'family-status'/f'{family}.json';data=read(p)
        next(r for r in data['routes'] if r['route']==route).update(status=status,**extra);write(p,data)

    def bundle_marker(self,key,phase='statistics_complete'):
        row=next(r for r in source_rows(self.series) if r['source_id']==key)
        bundle=Path(row['bundle']);analysis=bundle.parent
        write(bundle/'bundle-complete.json',{})
        write(row['config'],dict(analysis_root=str(analysis),analysis_route=key))
        write(analysis/'pipeline-status.json',dict(phase=phase))
        return row

    def test_family_plan_covers_all_25_and_overrides_stale_legacy_queue(self):
        write(self.series/'source-runs.json',dict(routes=[dict(source_id='stale')]))
        rows=source_rows(self.series);bindings=source_bindings(rows)
        self.assertEqual(len(bindings),25)
        self.assertEqual(sum(r['adopted'] for r in rows),2)
        for r in rows:
            expected=self.root/'old' if r['adopted'] else self.series
            self.assertEqual(Path(r['bundle']),expected/'runs'/r['source_id']/'statistics/bundle')

    def test_failure_is_read_from_new_family_without_spawning_science(self):
        self.change_status('ASNR','asnr_0','failed',error='current queue failure')
        with patch('continue_supplement_series.subprocess.run') as execute:
            self.assertEqual(run(self.series,self.contracts,self.out,once=True),0);execute.assert_not_called()
        state=read(self.out/'status.json')
        self.assertEqual(len(state['routes']),25)
        self.assertEqual(state['routes']['asnr_0']['status'],'upstream_failed')
        self.assertEqual(state['routes']['asnr_0']['error'],'current queue failure')
        self.assertEqual(state['routes']['base']['status'],'waiting_for_source_bundle')

    def test_early_bundle_waits_for_queue_acceptance_including_adopted_route(self):
        for key in ['base','asnr_0']:self.bundle_marker(key)
        with patch('continue_supplement_series.subprocess.run') as execute:
            run(self.series,self.contracts,self.out,once=True);execute.assert_not_called()

    def test_queue_acceptance_also_requires_final_dag_acceptance(self):
        self.bundle_marker('base',phase='running')
        self.change_status('uniform_count_base','base','statistics_complete')
        with patch('continue_supplement_series.subprocess.run') as execute:
            with self.assertRaisesRegex(ValueError,'final DAG acceptance'):
                run(self.series,self.contracts,self.out,once=True)
            execute.assert_not_called()

    def test_accepted_new_route_uses_current_paths_once_only(self):
        row=self.bundle_marker('asnr_0');self.change_status('ASNR','asnr_0','statistics_complete',bundle=row['bundle'])
        with patch('continue_supplement_series.subprocess.run',return_value=Mock(returncode=0)) as execute, \
             patch('continue_supplement_series.ResultBundle',return_value=Mock(identity={'source_key':'fixture'})):
            run(self.series,self.contracts,self.out,once=True)
            run(self.series,self.contracts,self.out,once=True)
            execute.assert_called_once();argv=execute.call_args.args[0]
            self.assertEqual(argv[argv.index('--bundle')+1],row['bundle'])
            self.assertEqual(argv[argv.index('--config')+1],row['config'])
        self.assertEqual(read(self.out/'status.json')['routes']['asnr_0']['status'],'complete')

    def test_changed_adoption_cannot_silently_rebind_on_restart(self):
        run(self.series,self.contracts,self.out,once=True)
        p=self.series/'family-plan.json';data=read(p)
        data['adoption']['routes']['base']['bundle']=str(self.root/'unrelated/bundle');write(p,data)
        with self.assertRaisesRegex(ValueError,'source bindings changed'):
            run(self.series,self.contracts,self.out,once=True)

    def test_incomplete_or_duplicate_family_rows_cannot_drop_sources(self):
        p=self.series/'family-status/ASNR.json';data=read(p);original=list(data['routes'])
        for bad in [original[:-1],original[:-1]+[original[0]]]:
            write(p,dict(routes=bad))
            with self.assertRaisesRegex(ValueError,'membership changed'):source_rows(self.series)

    def test_misbound_completed_family_cannot_supply_a_different_bundle(self):
        self.change_status('ASNR','asnr_0','statistics_complete',bundle=str(self.root/'old/stale/bundle'))
        with self.assertRaisesRegex(ValueError,'bundle binding changed'):source_rows(self.series)

    def test_common_main_status_does_not_accept_core_or_hide_failure(self):
        series=self.root/'common-series'
        row=dict(source_id='polyapipe_common',config=str(series/'config.json'),bundle=str(series/'bundle'))
        for status,expected in [('CORE_AND_MAIN_STATISTICS_RUNNING','CORE_AND_MAIN_STATISTICS_RUNNING'),
                                ('MAIN_STATISTICS_COMPLETE','statistics_complete'),('FAILED','failed')]:
            write(series/'source-runs.json',dict(routes=[dict(row,status=status)]))
            self.assertEqual(source_rows(series)[0]['status'],expected)
        with patch('continue_supplement_series.subprocess.run') as execute:
            self.assertEqual(run(series,self.contracts,self.out,once=True),1);execute.assert_not_called()
        self.assertEqual(read(self.out/'status.json')['routes']['polyapipe_common']['status'],'upstream_failed')

    def test_common_main_complete_dispatches_only_after_own_dag_acceptance(self):
        series=self.root/'common-series';bundle=series/'statistics/bundle';config=series/'pipeline.config.json'
        write(bundle/'bundle-complete.json',{})
        write(config,dict(analysis_root=str(bundle.parent),analysis_route='polyapipe_common'))
        write(bundle.parent/'pipeline-status.json',dict(phase='statistics_complete'))
        row=dict(source_id='polyapipe_common',config=str(config),bundle=str(bundle),status='CORE_AND_MAIN_STATISTICS_RUNNING')
        write(series/'source-runs.json',dict(routes=[row]))
        with patch('continue_supplement_series.subprocess.run',return_value=Mock(returncode=1)) as execute:
            self.assertEqual(run(series,self.contracts,self.out,once=True),0);execute.assert_not_called()
            write(series/'source-runs.json',dict(routes=[dict(row,status='MAIN_STATISTICS_COMPLETE')]))
            self.assertEqual(run(series,self.contracts,self.out,once=True),1);execute.assert_called_once()
            argv=execute.call_args.args[0]
            self.assertEqual(argv[argv.index('--bundle')+1],str(bundle))
            self.assertEqual(argv[argv.index('--config')+1],str(config))


if __name__=='__main__':unittest.main()
