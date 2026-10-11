"""Assemble completed annotation chunks in study order with bounded memory.

Chunking preserves each study's coordinate-table order. Check that invariant
for every row instead of collecting all studies into one large DataFrame.
"""
import argparse
from contextlib import ExitStack
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import time


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def progress(phase,**extra):
    print(json.dumps(dict(phase=phase,time_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),**extra)),flush=True)


def finalize(root, *, verify_lock=True, backend='bedtools', resumed_legacy_chunks=True):
    root=Path(root)
    marker=root/'annotation.complete.json'
    assert not marker.exists(), 'Annotation already finalized'
    manifest=json.loads((root/'intake.native.json').read_text())
    expected={x['study_id']:x['atlas_pas'] for x in manifest['studies']}
    lock=root/('annotation_inputs_bedtools.json' if backend=='bedtools' else 'annotation_inputs.json')
    if verify_lock:
        progress('verifying_annotation_inputs')
        for path,value in json.loads(lock.read_text()).items():
            assert sha(path)==value, f'Annotation input changed: {path}'
    chunks=sorted((root/'annotation_chunks').glob('*/*/complete.json'))
    assert chunks, 'No completed chunks'
    assert len(chunks)==len(list((root/'annotation_chunks').glob('*/*/features.tsv.gz'))), 'Incomplete chunk output exists'
    for study in expected:
        assert not (root/'studies'/study/'features.tsv.gz').exists(), f'Final table already exists: {study}'
    counts={study:0 for study in expected};writers={};readers={};headers={};temp={}
    annotation_header=None;chunk_receipts=[]
    progress('assembling_study_tables',chunks=len(chunks),expected_rows=sum(expected.values()))
    with ExitStack() as stack:
        for study in expected:
            folder=root/'studies'/study
            source=stack.enter_context(gzip.open(folder/'coordinate_features.tsv.gz','rt',newline=''))
            readers[study]=csv.DictReader(source,delimiter='\t')
            headers[study]=readers[study].fieldnames
            assert 'atlas_pas_id' in headers[study]
            temp[study]=folder/'features.assembling.tsv.gz'
            handle=stack.enter_context(temp[study].open('wb'))
            compressed=stack.enter_context(gzip.GzipFile(filename='',mode='wb',fileobj=handle,mtime=0))
            text=stack.enter_context(io.TextIOWrapper(compressed,encoding='utf-8',newline=''))
            writers[study]=csv.writer(text,delimiter='\t',lineterminator='\n')
        for i,path in enumerate(chunks,1):
            receipt=json.loads(path.read_text());feature=path.parent/'features.tsv.gz'
            assert sha(feature)==receipt['features_sha256'], f'Chunk checksum changed: {feature}'
            assert all((path.parent/n).is_file() for n in ['sequence_qc.tsv.gz','annotation_candidates.tsv.gz'])
            assert path.parent.parent.name==receipt['species'] and int(path.parent.name)==receipt['block']
            chunk_receipts.append(dict(path=str(path),sha256=sha(path),features_sha256=receipt['features_sha256']))
            n=0
            with gzip.open(feature,'rt',newline='') as f:
                data=csv.DictReader(f,delimiter='\t')
                if annotation_header is None:
                    annotation_header=data.fieldnames
                    extras={s:[c for c in annotation_header if c not in headers[s]] for s in expected}
                    for s in expected:writers[s].writerow(headers[s]+extras[s])
                assert data.fieldnames==annotation_header, f'Chunk schema changed: {feature}'
                for row in data:
                    pas=row['atlas_pas_id'];prefix,study,_=pas.split(':',2)
                    assert prefix=='polyApipe' and study in expected, f'Unexpected study PAS: {pas}'
                    raw=next(readers[study],None)
                    assert raw is not None and raw['atlas_pas_id']==pas, f'Missing, duplicate or reordered PAS: {pas}'
                    assert row['species']==receipt['species']
                    assert row['primary_assignment_status'], f'Missing assignment: {pas}'
                    counts[study]+=1;n+=1
                    assert int(raw['feature_index_1based'])==counts[study], f'Matrix axis changed: {pas}'
                    writers[study].writerow([raw[c] for c in headers[study]]+[row[c] for c in extras[study]])
            assert n==receipt['rows'], f'Chunk row count changed: {feature}'
            if i%100==0 or i==len(chunks):progress('assembling_study_tables',chunks_done=i,rows_done=sum(counts.values()),studies=counts)
        assert counts==expected, f'Study row counts differ: {counts} != {expected}'
        for study,reader in readers.items():
            assert next(reader,None) is None, f'Unannotated coordinate rows: {study}'
    # Publish only after every input, row, matrix-axis position and study count passed.
    outputs={}
    for study,path in temp.items():
        target=path.with_name('features.tsv.gz');path.replace(target)
        outputs[study]=dict(rows=counts[study],path=str(target),sha256=sha(target))
    record=dict(status='complete',rows=sum(counts.values()),chunks=len(chunks),backend=backend,
        resumed_legacy_chunks=resumed_legacy_chunks,assembly='streaming_study_order',time_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
        checks=dict(all_sites_annotated=True,external_original_source_coordinates=True,chunk_checksums=True,
                    study_matrix_axis_preserved=True,study_counts_match=True),
        assembly_code_sha256=sha(__file__),annotation_input_lock_sha256=sha(lock),studies=outputs)
    audit=root/'annotation_assembly_inputs.json'
    audit.write_text(json.dumps(chunk_receipts,indent=2)+'\n')
    record['assembly_input_audit_sha256']=sha(audit)
    temporary=marker.with_name('annotation.complete.tmp.json')
    temporary.write_text(json.dumps(record,indent=2)+'\n');temporary.replace(marker)
    progress('annotation_complete',rows=record['rows'],studies=counts)
    return record


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True)
    a=p.parse_args()
    try:finalize(a.run.resolve())
    except Exception as error:
        progress('annotation_assembly_failed',error=repr(error))
        raise
