"""Run the existing Figure 5 terminal-family rules against local count matrices."""
from pathlib import Path
from collections import Counter
from types import SimpleNamespace
import ast
import json
import hashlib
import importlib.util
import numpy as np
import pandas as pd
from scipy import sparse

HERE = Path(__file__).resolve().parent
PROJECTS = Path('/Users/lycidas/work/科研/projects')
F5 = PROJECTS/'2026-09-14_figure5-regulatory-results'
MOUSE = PROJECTS/'2026-09-02_mouse-time-aware-integration-pilot/data'
HUMAN = PROJECTS/'2026-09-09_human-time-aware-integration/data'
ROUTES = {'mouse': [('CIP','gse103221_pluripotency_reprogramming_core'),
                    ('WOT','gse115943_pluripotency_reprogramming_core')],
          'human': [('GSE118258','gse118258_human_pluripotency_reprogramming'),
                    ('GSE150311','gse150311_human_pluripotency_reprogramming')]}

def functions(path, names, namespace):
    """Load only frozen pure functions, excluding original remote execution drivers."""
    source = path.read_text()
    tree = ast.parse(source)
    body = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names]
    assert {n.name for n in body} == set(names)
    # Only filesystem aliases differ from the original prepare function.
    code = '\n\n'.join(ast.get_source_segment(source, n) for n in body)
    code = code.replace('results_metacell_pas_v7_v1', 'metacell-pas-v7-v1')
    code = code.replace('results_pas_annotation_support_v1', 'pas-annotation-support-v1')
    exec(compile(code, str(path), 'exec'), namespace)

def main():
    ns = dict(np=np, pd=pd, hashlib=hashlib, sparse=sparse, Counter=Counter, json=json)
    functions(HERE/'legacy_prepare_pilot.py', ['sha','strict_pas','models_for_gene'], ns)
    spec = importlib.util.spec_from_file_location('structure', HERE/'legacy_structure_consensus.py')
    structure = importlib.util.module_from_spec(spec); spec.loader.exec_module(structure)
    legacy = SimpleNamespace(MOUSE=MOUSE, HUMAN=HUMAN, PROFILES=MOUSE/'apa-majority-shared-v1',
        sha=ns['sha'], strict_pas=ns['strict_pas'], resolve_families=structure.resolve_families,
        read=lambda p: pd.read_csv(p, sep='\t', low_memory=False))
    ns['legacy'] = legacy
    functions(F5/'full-cohort-v2/run_full.py', ['dump','prepare'], ns)
    claims = pd.read_csv(HERE/'figure5_union_claims.tsv', sep='\t')
    for sp in ['mouse','human']:
        out = HERE/'native'/sp
        if (out/'complete.json').exists():
            print(sp, 'already complete', flush=True); continue
        out.mkdir(parents=True, exist_ok=False)
        sub = claims[claims.species.eq(sp)]
        utr = MOUSE/'tandem-apa-v1/all_reference_terminal_utr_candidates.tsv.gz' if sp=='mouse' else HERE/'human_terminal_utrs.tsv.gz'
        cpath = HERE/f'{sp}_catalog_utrs.tsv.gz'
        hashes = {str(p): ns['sha'](p) for p in [utr,cpath,HERE/'figure5_union_claims.tsv',F5/'full-cohort-v2/run_full.py',HERE/'legacy_prepare_pilot.py',HERE/'legacy_structure_consensus.py']}
        u = legacy.read(utr); c = legacy.read(cpath)
        cg = dict(tuple(c[c.gene_id.isin(sub.gene_id)].groupby('gene_id')))
        models = {gid: ns['models_for_gene'](g, cg.get(gid,c.iloc[:0]))
                  for gid,g in u[u.gene_id.isin(sub.gene_id)].groupby('gene_id')}
        audit, replay = [], []
        for route, study in ROUTES[sp]:
            genes, rows = ns['prepare'](sp,route,study,sub,models,out,hashes)
            audit.extend(rows)
            genes = json.loads((out/route/'models.json').read_text())
            old = json.loads((F5/f'full-cohort-v2/data/{sp}/{route}/models.json').read_text())
            for gid, model in old.items():
                new = genes[gid]
                assert model['family'] == new['family'] and model['pas'] == new['pas']
                for state, profile in model['profiles'].items():
                    assert profile == new['profiles'][state], (sp, route, gid, state)
                    replay.append((route,gid,state))
        pd.DataFrame(audit).to_csv(out/'eligibility.tsv',sep='\t',index=False)
        (out/'complete.json').write_text(json.dumps(dict(status='PASS',input_sha256=hashes,
            candidate_genes=sub.gene_id.nunique(),original_profiles_reproduced=len(replay),
            scope='Terminal models and raw library counts; recognition sites not yet recomputed'),indent=2)+'\n')
        print(sp, 'original profiles reproduced',len(replay),flush=True)

if __name__ == '__main__':
    main()
