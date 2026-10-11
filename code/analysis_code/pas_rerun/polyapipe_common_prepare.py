"""Prepare accepted, study-specific common GFF counts without reclustering PAS."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np
import pandas as pd
from scipy import sparse
from polyapipe_prepare import digest, save


def write(path, value):
    path = Path(path)
    tmp = path.with_suffix('.tmp.json')
    tmp.write_text(json.dumps(value, indent=2) + '\n'); tmp.replace(path)


def checked(record):
    path = Path(record['path'])
    assert digest(path) == record['sha256'], f'Changed accepted input: {path}'
    return path


def catalogue(record):
    rows = []
    with checked(record['catalogue']).open() as handle:
        for line in handle:
            if line.startswith('#') or not line.strip(): continue
            chrom, _, kind, start, end, _, strand, _, attributes = line.rstrip().split('\t')
            attrs = dict(re.findall(r'(\w+)="([^"]*)"', attributes))
            assert kind == 'polyAends' and strand in ['+', '-'] and int(end) >= int(start) >= 1
            assert attrs['misprime'] in ['True', 'False']
            rows.append(dict(chrom=chrom, strand=strand, native_peak_id=attrs['peak'],
                native_window_start1=int(start), native_window_end1=int(end),
                native_misprime=attrs['misprime'] == 'True', native_peakdepth=int(attrs['peakdepth']),
                cut_site=int(end) if strand == '+' else int(start)-1))
    f = pd.DataFrame(rows).sort_values(['chrom', 'strand', 'cut_site', 'native_peak_id']).reset_index(drop=True)
    assert len(f) == record['native_pas_count'] and f.native_peak_id.is_unique
    f['atlas_column_0based'] = np.arange(len(f)); f['feature_index_1based'] = np.arange(1, len(f)+1)
    f['atlas_pas_id'] = 'polyApipe:' + record['study_id'] + ':' + f.native_peak_id
    f['species'] = record['species']
    f['genome_build'] = 'GRCh38_2024_A' if record['species'] == 'human' else 'GRCm39_2024_A'
    f['member_cut_site_min'] = f.cut_site; f['member_cut_site_max'] = f.cut_site
    f['cluster_width_bp'] = 0; f['native_peaks'] = 1
    f['annotation_cut_site_0based_base'] = f.cut_site - f.strand.eq('+').astype(int)
    return f


def member_axis(members, meta, library):
    assert members.cell_id.is_unique and members.barcode.is_unique
    assert len(members) == library['cells'] and members.metacell_id.nunique() == library['metacells']
    cb = members.barcode.str.split('-').str[0]
    assert cb.is_unique, 'Library-local normalized barcode collision'
    assert set(members.metacell_id) == set(meta.metacell_id)
    sizes = members.groupby('metacell_id').size().reindex(meta.metacell_id)
    np.testing.assert_array_equal(sizes.to_numpy(), meta.n_cells.to_numpy())
    cell_study = 'GSE115943' if library['study_id'] == 'WOT_GSE115943' else library['study_id']
    cell_library = library.get('library_key', library['library_id'])
    cell_library = cell_library.split('::',1)[1] if '::' in cell_library else library['library_id']
    assert members.cell_id.str.startswith(cell_study+'|'+cell_library+'|').all()
    index = pd.Index(meta.metacell_id).get_indexer(members.metacell_id)
    assert (index >= 0).all()
    return pd.DataFrame(dict(barcode=cb, row=index))


def library_matrix(task):
    library, dest, frozen, executable = task
    dest, frozen = Path(dest), Path(frozen)
    target = dest/'libraries'/library['library_id']; target.mkdir(parents=True, exist_ok=False)
    members = pd.read_csv(checked(library['frozen_members']), sep='\t')
    barcodes = checked(library['frozen_barcodes']).read_text().splitlines()
    # The frozen whitelist retains -1; only the count file uses normalized CB.
    assert set(barcodes) == set(members.barcode) and len(barcodes) == len(members)
    meta = pd.read_csv(frozen/'metacells.tsv.gz', sep='\t')
    meta = meta[meta.library_key.eq(library['library_key'])].copy()
    mapping = member_axis(members, meta, library)
    mapping.to_csv(target/'CB-map.tsv', sep='\t', header=False, index=False)
    ncols = int((dest/'PAS-columns.txt').read_text())
    counts = checked(library['counts']); original_stat = counts.stat()
    binary, audit = target/'matrix.bin', target/'aggregate.json'
    subprocess.run([str(executable), str(counts), str(dest/'PAS-map.tsv'), str(target/'CB-map.tsv'),
        str(len(meta)), str(ncols), str(binary), str(audit)], check=True)
    assert (counts.stat().st_size, counts.stat().st_mtime_ns) == (original_stat.st_size, original_stat.st_mtime_ns)
    a = json.loads(audit.read_text())
    assert a['records'] == library['nonzero_cell_pas_records'] and a['umi'] == library['umi_total']
    assert a['observed_cells'] == library['observed_count_cells']
    data = np.fromfile(binary, dtype=[('row','<u4'),('column','<u4'),('count','<u8')])
    assert len(data) == a['nnz'] and (data['count'] <= np.iinfo(np.int64).max).all()
    matrix = sparse.csr_matrix((data['count'].astype(np.int64), (data['row'], data['column'])), shape=(len(meta), ncols))
    assert int(matrix.sum()) == a['umi'] and matrix.nnz == len(data)
    sparse.save_npz(target/'metacell_by_pas_counts.npz', matrix)
    save(meta[['metacell_id', 'n_cells']], target/'metacells.tsv.gz')
    binary.unlink()
    a.update(library_id=library['library_id'], library_key=library['library_key'],
        integer_molecules=a['umi'], integer_count_conservation=True, cells=len(members), metacells=len(meta),
        inputs={r['path']:r['sha256'] for r in [library['counts'],library['frozen_members'],library['frozen_barcodes']]})
    write(target/'complete.json', a)
    print(json.dumps(dict(library=library['library_id'], phase='aggregated', UMI=a['umi'])), flush=True)
    return a


def prepare(handoff, frozen, out, executable, workers=24):
    handoff, frozen, out = Path(handoff), Path(frozen), Path(out)
    h = json.loads(handoff.read_text()); assert h['source_id'] == 'polyapipe_common' and h['reference_variant'] == 'v2_primary'
    assert h['library_count'] == 161 and len(h['studies']) == 7
    out.mkdir(parents=True, exist_ok=False)
    inputs = {str(handoff):digest(handoff), str(executable):digest(executable), str(Path(__file__).resolve()):digest(__file__)}
    summaries, all_audits, tasks = [], [], []
    for study in h['studies']:
        sid = study['study_id']; dest = out/'studies'/sid; dest.mkdir(parents=True)
        f = catalogue(study); save(f, dest/'coordinate_features.tsv.gz')
        inputs[study['catalogue']['path']] = study['catalogue']['sha256']
        f[['native_peak_id','atlas_column_0based']].to_csv(dest/'PAS-map.tsv', sep='\t', header=False, index=False)
        (dest/'PAS-columns.txt').write_text(str(len(f)))
        tasks.extend((lib,str(dest),str(frozen/'studies'/sid),str(executable)) for lib in h['libraries'] if lib['study_id'] == sid)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        all_audits = list(pool.map(library_matrix, tasks, chunksize=1))
    for study in h['studies']:
        sid = study['study_id']; dest = out/'studies'/sid; source = frozen/'studies'/sid
        f = pd.read_csv(dest/'coordinate_features.tsv.gz', sep='\t')
        audits = [a for a in all_audits if any(l['library_id']==a['library_id'] and l['study_id']==sid for l in h['libraries'])]
        meta = pd.read_csv(source/'metacells.tsv.gz', sep='\t')
        blocks, indexes, totals = [], [], []
        for a in audits:
            root = dest/'libraries'/a['library_id']; c = sparse.load_npz(root/'metacell_by_pas_counts.npz')
            blocks.append(c); indexes.extend(pd.read_csv(root/'metacells.tsv.gz', sep='\t').metacell_id)
            totals.append(sparse.csr_matrix(c.sum(0)))
            inputs.update(a['inputs'])
        assert len(set(indexes)) == len(indexes) == study['metacells']
        order = pd.Index(indexes).get_indexer(meta.metacell_id); assert (order >= 0).all()
        counts = sparse.vstack(blocks, format='csr')[order].tocsr()
        f['current_target_umi'] = np.asarray(counts.sum(0)).ravel()
        library_counts = sparse.vstack(totals, format='csr')
        f['n_libraries'] = np.asarray((library_counts > 0).sum(0)).ravel()
        f['atlas_status'] = np.where(f.n_libraries >= 2, 'core', 'single_library_or_unobserved')
        save(f, dest/'coordinate_features.tsv.gz'); sparse.save_npz(dest/'metacell_by_pas_counts.npz', counts)
        sparse.save_npz(dest/'library_by_pas_counts.npz', library_counts)
        save(pd.DataFrame(audits).drop(columns='inputs'), dest/'native_library_audit.tsv')
        for name in ['metacells.tsv.gz','all_metacells.tsv.gz','selected_metacells.tsv.gz','sampling_design.tsv','endpoint_library_design.tsv']:
            inputs[str(source/name)] = digest(source/name)
            if name in ['metacells.tsv.gz','all_metacells.tsv.gz']:
                m = pd.read_csv(source/name, sep='\t')
                stale = [c for c in m if 'pas' in c.lower() or c.startswith('genes_with_ge2') or c=='primary_detected_genes']
                save(m.drop(columns=stale), dest/name)
            else: shutil.copy2(source/name, dest/name)
        selected = pd.read_csv(dest/'selected_metacells.tsv.gz', sep='\t')
        np.testing.assert_array_equal(meta.iloc[selected.matrix_row].metacell_id, selected.metacell_id)
        assert int(meta.n_cells.sum()) == study['cells'] and int(counts.sum()) == sum(a['umi'] for a in audits)
        summaries.append(dict(study_id=sid,species=study['species'],libraries=len(audits),cells=study['cells'],
            metacells=len(meta),atlas_pas=len(f),molecules=int(counts.sum()),count_conservation=True))
        for a in audits:
            (dest/'libraries'/a['library_id']/'metacell_by_pas_counts.npz').unlink()
    assert sum(s['molecules'] for s in summaries) == h['umi_total']
    assert sum(s['cells'] for s in summaries) == 551201 and sum(s['metacells'] for s in summaries) == 11104
    save(pd.DataFrame(all_audits).drop(columns='inputs'),out/'frozen_native_membership_audit.tsv')
    write(out/'frozen_native_membership_audit.json',dict(status='PASS',libraries=161,cells=551201,metacells=11104,
        checks=dict(exact_frozen_barcode_membership=True,integer_counts_conserved=True,full_native_PAS_identity=True,zero_rows_preserved=True)))
    write(out/'intake.native.json',dict(status='native_atlas_prepared_not_full_rerun',caller='polyApipe',
        source_id='polyapipe_common',analysis_route='polyapipe_common',reference_variant='v2_primary',
        run_id=h['proposed_run_id'],count_input_fingerprint=h['count_input_fingerprint'],
        cut_coordinate_system='0-based-interbase',catalogue_scope='separate common GFF per study; no reclustering',
        inputs=inputs,studies=summaries))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ['handoff','frozen','out','executable']: p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--workers',type=int,default=24); a = p.parse_args()
    prepare(a.handoff,a.frozen,a.out,a.executable,a.workers)
