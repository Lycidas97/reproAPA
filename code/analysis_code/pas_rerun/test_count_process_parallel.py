"""Process parallelism must preserve each conditional test and ordered result."""
import json
import os
from pathlib import Path
import tempfile
import unittest
import numpy as np
from scipy.stats import fisher_exact,random_table
from downstream_runtime import Context
from downstream_inference import BASE,count_results


class CountProcessTests(unittest.TestCase):
    def context(self,root):
        if os.environ.get('PAS_TEST_CONFIG'):
            return Context(os.environ['PAS_TEST_CONFIG'],'count_inference',root/'count')
        project=Path(__file__).resolve().parents[2]
        cfg=root/'config.json'
        cfg.write_text(json.dumps(dict(project_root=str(project),native_run=str(root/'native'),
            frozen_gex=str(root/'gex'),analysis_root=str(root/'analysis'),run_id='test-process',
            input_fingerprint='test-input',source_overrides={})))
        return Context(cfg,'count_inference',root/'count')

    def test_one_and_multiple_processes_match_exact_and_MC_results_bit_for_bit(self):
        with tempfile.TemporaryDirectory() as tmp:
            ctx=self.context(Path(tmp))
            from conditional_count_test import fisher as f
            # Multiple chunks exercise more than one worker; zero columns are test-only removals.
            jobs=[(i,np.array([12,0,4]) if i%3==0 else np.array([12,3,4]),
                   np.array([3,0,14]) if i%3==0 else np.array([3,11,14]),12300+i) for i in range(40)]
            expected=[f(job) for job in jobs]
            self.assertEqual(list(count_results(ctx,jobs,1)),expected)
            self.assertEqual(list(count_results(ctx,jobs,3)),expected)
            self.assertEqual({r[2] for r in expected},{'Fisher_exact_2x2','Fisher_MC_2xn'})


if __name__=='__main__':unittest.main()
