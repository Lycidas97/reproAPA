"""Build the scan-ready CisBP-RNA motif-to-source-to-target ledger."""
import csv
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESOURCE = ROOT / 'projects/2026-05-25_cell-reprogramming-scrna-datasets/analysis'

LAYERS = [
    ('human_direct', RESOURCE / 'human-regulatory-query-v1/resources/human_direct_motif_audit.tsv',
     RESOURCE / 'human-regulatory-query-v1/resources/human_direct_scannable.meme', 'scan_ready',
     'Homo_sapiens', 'CisBP-RNA2_Homo_sapiens.zip'),
    ('mouse_direct', RESOURCE / 'mouse-regulatory-catalog-v1/outputs/mouse_direct_motif_provenance.tsv',
     None, 'scan_ready', 'Mus_musculus', 'CisBP-RNA2_Mus_musculus.zip'),
    ('human_to_mouse_transfer', RESOURCE / 'mouse-regulatory-query-v2/resources/human_to_mouse_motif_audit.tsv',
     RESOURCE / 'mouse-regulatory-query-v2/resources/transferred_scannable.meme', 'ready',
     'Mus_musculus', 'CisBP-RNA2_Homo_sapiens.zip'),
]

FIELDS = [
    'layer', 'database', 'database_build', 'source_archive', 'motif_id', 'cisbp_rbp_id',
    'source_species', 'source_gene_id', 'source_rbp', 'source_evidence', 'motif_source_id',
    'original_assay_id', 'original_assay_type', 'original_source', 'original_author',
    'original_year', 'original_pmid', 'original_source_version', 'target_species',
    'target_gene_id', 'target_rbp', 'transfer_gate', 'ortholog_ensembl_release',
    'rbd_evidence', 'scan_status', 'audit_file',
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def meme_ids(path):
    return {line.split()[1] for line in path.read_text().splitlines() if line.startswith('MOTIF ')}


def main():
    all_rows, audits = [], []
    hashes = {str(Path(__file__)): sha(Path(__file__))}
    for layer, audit, meme, status, target_species, archive in LAYERS:
        with audit.open(newline='') as handle:
            source = list(csv.DictReader(handle, delimiter='\t'))
        hashes[str(audit)] = sha(audit)
        ready = [r for r in source if r['scan_status'] == status]
        ids = {r['Motif_ID'] for r in ready}
        assert len(ready) == len(ids), (layer, 'motif IDs not unique')
        assert all(mid.endswith('_2.00') for mid in ids)
        if meme:
            hashes[str(meme)] = sha(meme)
            assert ids == meme_ids(meme), (layer, len(ids - meme_ids(meme)), len(meme_ids(meme) - ids))
        for r in ready:
            assert r['RBP_Status'] == 'D' and r['RBP_Species'] in ('Homo_sapiens', 'Mus_musculus')
            transfer = layer == 'human_to_mouse_transfer'
            if transfer:
                assert r['RBP_Species'] == 'Homo_sapiens'
                assert r['mouse_gene_id'].startswith('ENSMUSG') and r['mouse_symbol']
                assert r['gate_status'] == 'one2one_plus_database_RBD_transfer'
                assert r['ensembl_release'] == '115'
                assert r['rbd_identity_recomputed'].lower() == 'false'
            else:
                assert r['RBP_Species'] == target_species
            gene_id = r.get('gene_id') or r.get('DBID')
            all_rows.append({
                'layer': layer, 'database': 'CisBP-RNA', 'database_build': '2.00',
                'source_archive': archive, 'motif_id': r['Motif_ID'],
                'cisbp_rbp_id': r['RBP_ID'], 'source_species': r['RBP_Species'],
                'source_gene_id': gene_id, 'source_rbp': r['RBP_Name'],
                'source_evidence': 'source-RBP direct experimental motif',
                'motif_source_id': r['MSource_ID'],
                'original_assay_id': r.get('assay_id') or r.get('DBID.1'),
                'original_assay_type': r['Motif_Type'],
                'original_source': r['MSource_Identifier'],
                'original_author': r['MSource_Author'],
                'original_year': r['MSource_Year'],
                'original_pmid': r['PMID'],
                'original_source_version': r['MSource_Version'],
                'target_species': target_species,
                'target_gene_id': r['mouse_gene_id'] if transfer else gene_id,
                'target_rbp': r['mouse_symbol'] if transfer else r['RBP_Name'],
                'transfer_gate': r['gate_status'] if transfer else 'none',
                'ortholog_ensembl_release': r['ensembl_release'] if transfer else '',
                'rbd_evidence': r['rbd_evidence'] if transfer else '',
                'scan_status': r['scan_status'], 'audit_file': str(audit),
            })
        audits.append({'layer': layer, 'audited_rows': len(source),
                       'scan_ready_motifs': len(ready),
                       'scan_ready_factors': len({r['RBP_Name'] for r in ready}),
                       'meme_id_match': 'PASS' if meme else 'not_local'})
    out = HERE / 'motif_source_transfer.tsv'
    with out.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter='\t')
        writer.writeheader()
        writer.writerows(sorted(all_rows, key=lambda r: (r['layer'], r['source_rbp'], r['motif_id'])))
    assert len(all_rows) == 468
    record = {'status': 'PASS', 'database': 'CisBP-RNA', 'build': '2.00',
              'scannable_motif_records': len(all_rows), 'layers': audits,
              'input_sha256': hashes, 'output_sha256': sha(out)}
    (HERE / 'verification.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps({k: v for k, v in record.items() if k != 'input_sha256'}, indent=2))


if __name__ == '__main__':
    main()
