"""Small read-only reciprocal-chain batch and frozen one-base PAS calibration."""
import hashlib
import json
import subprocess

import pandas as pd

from analyze_sites import HERE, ROOT, TABLES, read, save, INPUTS


def main():
    queries = read(HERE / 'seed_mapping_queries.tsv')
    config_path = ROOT / 'notes/2026-10-06-dual-main-plot-data/provenance/polyapipe-pipeline.config.json'
    config = json.loads(config_path.read_text())
    calibration = []
    for species in ['human', 'mouse']:
        source = read(TABLES / f'cis/liftover/{species}_reciprocal.audit.tsv.gz')
        for (_, _), group in source.groupby(['mapping_status', 'strand']):
            for row in group.sort_values('id').head(4).itertuples():
                start = int(row.cut_site) - int(row.strand == '+')
                calibration.append(dict(query_id='CAL:' + row.id, species=species, gene_id=row.id.split('|')[1],
                                        seed='CAL', chrom=row.chrom, strand=row.strand, start0=start, end0=start + 1,
                                        expected_status=row.mapping_status, expected_chrom=row.mapped_chrom,
                                        expected_strand=row.mapped_strand, expected_cut=row.mapped_cut))
    calibration = pd.DataFrame(calibration)
    save(calibration, 'chain_calibration_queries.tsv')
    all_queries = pd.concat([queries, calibration[queries.columns]], ignore_index=True).to_dict('records')
    code = 'REMOTE_RESOURCES=' + repr(config['resources']) + '\nREMOTE_QUERIES=' + repr(all_queries) + '\n' + (HERE / 'chain_mapper.py').read_text()
    command = 'PYTHONDONTWRITEBYTECODE=1 /root/miniconda3/envs/apasim_analysis/bin/python -'
    with (HERE / 'seed_chain_maps.partial.jsonl').open('wb') as output, (HERE / 'seed_chain_remote.log').open('wb') as log:
        subprocess.run(['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=20', '213docker', command],
                       input=code.encode(), stdout=output, stderr=log, check=True)
    results = [json.loads(line) for line in (HERE / 'seed_chain_maps.partial.jsonl').read_text().splitlines()]
    lookup = {r['query_id']: r for r in results}
    assert len(results) == len(all_queries) == len(lookup)
    validations = []
    for row in calibration.itertuples():
        result = lookup[row.query_id]
        assert result['status'] == row.expected_status, (row.query_id, result, row.expected_status)
        if row.expected_status == 'unique':
            hit = result['hits'][0]
            cut = hit['end0'] if hit['strand'] == '+' else hit['start0']
            assert (hit['chrom'], hit['strand'], cut) == (row.expected_chrom, row.expected_strand, row.expected_cut)
        validations.append(dict(query_id=row.query_id, actual_status=result['status'], expected_status=row.expected_status, PASS=True))
    save(pd.DataFrame(validations), 'chain_calibration_validation.tsv')
    (HERE / 'seed_chain_maps.partial.jsonl').rename(HERE / 'seed_chain_maps.jsonl')
    report = dict(status='PASS', seed_queries=len(queries), calibration_queries=len(calibration),
                  calibrated_point_mapping_exact=True, resources=[json.loads(line) for line in (HERE / 'seed_chain_remote.log').read_text().splitlines() if line.startswith('{')],
                  config_sha256=hashlib.sha256(config_path.read_bytes()).hexdigest(),
                  query_sha256=hashlib.sha256((HERE / 'seed_mapping_queries.tsv').read_bytes()).hexdigest(),
                  mapper_sha256=hashlib.sha256((HERE / 'chain_mapper.py').read_bytes()).hexdigest())
    (HERE / 'seed_chain_validation.json').write_text(json.dumps(report, indent=2) + '\n')
    (HERE / 'chain_input_hashes.json').write_text(json.dumps(INPUTS, indent=2) + '\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
