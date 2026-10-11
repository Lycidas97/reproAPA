"""Freeze motif candidates and collect their raw expression without remote writes."""
from pathlib import Path
import hashlib
import json
import subprocess
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
M = Path('/Users/lycidas/work/科研/projects/2026-09-02_mouse-time-aware-integration-pilot')
STUDIES = {
    'GSE103221': ['mouse', 'gse103221_pluripotency_reprogramming_core', 'CIP'],
    'GSE115943': ['mouse', 'gse115943_pluripotency_reprogramming_core', 'WOT'],
    'GSE118258': ['human', 'gse118258_human_pluripotency_reprogramming', 'GSE118258'],
    'GSE150311': ['human', 'gse150311_human_pluripotency_reprogramming', 'GSE150311'],
}


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    path = ROOT / 'notes/2026-10-01-figure5-union-update/analysis-paired-support/exact_site_cross_study.tsv.gz'
    d = pd.read_csv(path, sep='\t')
    d = d[d.category.eq('RBP') & d.recurrent_strict].copy()
    assert d.provenance_a.eq(d.provenance_b).all()
    assert (d.effect20_a * d.effect20_b > 0).all()
    d['evidence'] = d.provenance_a.map({'RBP': 'direct_species', 'RBP_transfer': 'transferred_motif'})
    assert d.evidence.notna().all()
    d = d.drop_duplicates(['species', 'factor', 'gene_id', 'chrom', 'strand', 'start0', 'end0', 'contrast'])
    d.to_csv(HERE / 'candidate_sites.tsv.gz', sep='\t', index=False, compression={'method': 'gzip', 'mtime': 0})
    c = d.groupby(['species', 'factor', 'evidence']).agg(target_genes=('gene_id', 'nunique'),
        site_contrast_links=('factor', 'size')).reset_index()
    assert not c.duplicated(['species', 'factor']).any()
    c.to_csv(HERE / 'candidates.tsv', sep='\t', index=False)
    config = dict(studies=STUDIES, factors={sp: sorted(c[c.species.eq(sp)].factor) for sp in ['mouse', 'human']},
        metadata_sha256={study: digest(M / 'data/apa-majority-shared-v1/profiles' / spec[1] / 'selected_metacells.tsv.gz') for study, spec in STUDIES.items()})
    (HERE / 'config.json').write_text(json.dumps(config, indent=2) + '\n')
    (HERE / 'candidate_source.json').write_text(json.dumps(dict(path=str(path), sha256=digest(path)), indent=2) + '\n')
    code = 'CONFIG = ' + repr(config) + '\n' + (HERE / 'extract_remote.py').read_text()
    command = 'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /root/miniconda3/envs/scanpy/bin/python -'
    partial = HERE / 'expression_inputs.partial.npz'
    with partial.open('wb') as output, (HERE / 'extract.log').open('wb') as log:
        subprocess.run(['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=20', '213docker', command],
            input=code.encode(), stdout=output, stderr=log, check=True)
    partial.rename(HERE / 'expression_inputs.npz')
    print(c.groupby(['species', 'evidence']).size().to_string())
    print((HERE / 'extract.log').read_text())


if __name__ == '__main__':
    main()
