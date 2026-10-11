"""A failed node must be visible while an independent sibling is still finishing."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import sys
import tempfile
import time
import unittest
from pipeline import execute,read,write,digest


class FailureStatusTests(unittest.TestCase):
    def test_failure_is_published_before_waiting_for_running_sibling(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);script=root/'node.py';release=root/'release';started=root/'started'
            script.write_text("""import argparse,pathlib,time,sys
p=argparse.ArgumentParser();p.add_argument('--config');p.add_argument('--out');a=p.parse_args()
root=pathlib.Path(a.config).parent
if pathlib.Path(a.out).name=='slow':
 (root/'started').touch()
 while not (root/'release').exists():time.sleep(.02)
else:
 while not (root/'started').exists():time.sleep(.02)
sys.exit(1)
""")
            cfg=root/'config.json';write(cfg,dict(prepared=True,run_id='test',input_fingerprint='input',
                analysis_root=str(root/'run'),jobs={n:dict(script=str(script),arguments=[],depends_on=[],python=sys.executable) for n in ['slow','fail']},
                code_sha256={str(script):digest(script)},input_files={},resources={}))
            with ThreadPoolExecutor(max_workers=1) as pool:
                future=pool.submit(execute,cfg,2)
                try:
                    deadline=time.monotonic()+10;state={}
                    while time.monotonic()<deadline:
                        p=root/'run/pipeline-status.json'
                        if p.exists():state=read(p)
                        if state.get('phase')=='failed':break
                        time.sleep(.02)
                    self.assertEqual(state.get('phase'),'failed')
                    self.assertEqual(state['failed_node'],'fail')
                    self.assertEqual(state['draining_nodes'],['slow'])
                    self.assertFalse(future.done())
                finally:release.touch()
                with self.assertRaises(Exception):future.result(timeout=10)


if __name__=='__main__':unittest.main()
