"""Freeze an outcome-independent mouse CIP pilot from existing count/state axes."""
import hashlib
import json
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

BASE = Path('/root/nfsdata/apadetect/analysis_projects')
ROOT = BASE/'mouse_time_aware_integration_pilot_v1'
OUT = Path(__file__).resolve().parents[1]
STUDY = 'gse103221_pluripotency_reprogramming_core'
STATES = ['S', 'R', 'E', 'T', 'M']
sys.path.insert(0, str(BASE/'pas_regulatory_comparison_interface_v1/scripts'))
from reference_model import included_lengths


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def write(frame, name):
    frame.to_csv(OUT/'inputs'/name, sep='\t', index=False)


def strict_pas(f):
    required = ['keep_for_apa_primary', 'transcript_model_conflict', 'multi_gene_ambiguous', 'internal_priming_flag']
    for col in required:
        if f[col].isna().any() or not f[col].isin([True, False]).all():
            raise ValueError('Unknown/nonboolean PAS QC: '+col)
    return (f.keep_for_apa_primary & ~f.transcript_model_conflict & ~f.multi_gene_ambiguous
            & ~f.internal_priming_flag & f.sequence_qc_status.eq('computed') & f.primary_gene_id.notna())


def models_for_gene(u, catalog):
    u = u.copy()
    u['anchor'] = np.where(u.strand.eq('+'), u.utr_start_1based, u.utr_end_1based)
    u['exon_anchor'] = np.where(u.strand.eq('+'), u.exon_start_1based, u.exon_end_1based)
    keys = ['gene_id', 'chrom', 'strand', 'stop_start', 'stop_end', 'anchor', 'exon_anchor']
    models = []
    for key, group in u.groupby(keys, sort=True):
        best = group.sort_values(['utr_length_nt', 'transcript_id'], ascending=[False, True]).iloc[0]
        start, end = int(best.utr_start_1based)-1, int(best.utr_end_1based)
        hits = catalog[catalog.annotated_start0.eq(start) & catalog.annotated_end0.eq(end)
                       & catalog.chrom.eq(best.chrom) & catalog.strand.eq(best.strand)]
        sid = 'STR_'+hashlib.sha256('|'.join(map(str,key)).encode()).hexdigest()[:16]
        models.append(dict(structure_id=sid, gene_id=best.gene_id, chrom=best.chrom, strand=best.strand,
            start0=start, end0=end, annotated_length=end-start,
            stop_start=int(best.stop_start), stop_end=int(best.stop_end), exon_anchor=int(best.exon_anchor),
            transcript_ids=';'.join(sorted(group.transcript_id.unique())),
            region_id=hits.region_id.iloc[0] if len(hits)==1 else '',
            catalog_status='matched' if len(hits)==1 else 'missing_or_ambiguous_catalog_region'))
    return models


def assign_pas(pas, models):
    rows = []
    for p in pas.to_dict('records'):
        matches = [m for m in models if p['chrom']==m['chrom'] and p['strand']==m['strand']
            and m['start0'] < p['member_cut_site_min'] <= p['member_cut_site_max'] <= m['end0']]
        p.update(gene_id=p['primary_gene_id'], structure_id='', region_id='', included_min=0, included_max=0,
                 structure_status='no_annotated_terminal_structure', n_compatible_structures=len(matches))
        if len(matches)>1:
            p['structure_status'] = 'ambiguous_terminal_structure'
        elif len(matches)==1:
            m = matches[0]
            p['structure_status'] = m['catalog_status']
            if m['catalog_status']=='matched':
                lo,hi = included_lengths(m['start0'],m['end0'],p['strand'],p['member_cut_site_min'],p['member_cut_site_max'])
                p.update(structure_id=m['structure_id'],region_id=m['region_id'],included_min=lo,included_max=hi)
        rows.append(p)
    # Cluster intervals, not representative cut points, determine whether ordering is safe.
    for sid in {r['structure_id'] for r in rows} - {''}:
        members = sorted([r for r in rows if r['structure_id']==sid], key=lambda r:r['included_min'])
        if any(a['included_max']>=b['included_min'] for a,b in zip(members,members[1:])):
            for r in members:
                r.update(structure_id='', region_id='', structure_status='overlapping_pas_clusters')
    return rows


def main():
    (OUT/'inputs').mkdir(exist_ok=False)
    source = ROOT/'results_metacell_pas_v7_v1/studies'/STUDY
    paths = dict(features=source/'features.tsv.gz', matrix=source/'metacell_by_pas_counts.npz',
        metacells=source/'metacells.tsv.gz', annotations=ROOT/'results_apa_metacell_study_local_v1/inputs/metacell_coarse_annotations.tsv',
        utrs=ROOT/'results_tandem_apa_v1/all_reference_terminal_utr_candidates.tsv.gz',
        catalog=ROOT/'regulatory_catalog_v1/catalog.sqlite', bundle=ROOT/'regulatory_query_v3/resources/bundle.json',
        reference_math=BASE/'pas_regulatory_comparison_interface_v1/scripts/reference_model.py')
    manifest = {k:dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size) for k,p in paths.items()}
    m = pd.read_csv(paths['annotations'],sep='\t')
    m = m[m.pas_dataset_id.eq(STUDY)].sort_values('pas_matrix_row_index_0based').reset_index(drop=True)
    axis = pd.read_csv(paths['metacells'],sep='\t')
    assert m.metacell_id.is_unique and np.array_equal(m.metacell_id,axis.metacell_id)
    assert np.array_equal(m.pas_matrix_row_index_0based,np.arange(len(m)))
    assert np.array_equal(m.library_key,axis.library_key)
    f = pd.read_csv(paths['features'],sep='\t',low_memory=False)
    assert f.atlas_pas_id.is_unique and np.array_equal(f.feature_index_1based,np.arange(1,len(f)+1))
    assert set(f.species)=={'mouse'} and set(f.genome_build)=={'GRCm39_2024_A'}
    raw = sparse.load_npz(paths['matrix']).tocsr()
    assert raw.shape==(len(m),len(f))
    assert np.isfinite(raw.data).all() and (raw.data>=0).all() and np.equal(raw.data,np.floor(raw.data)).all()
    valid = (m.dataset.eq('gse103221_cip') & m.coarse_label_purity.ge(.8)
             & m.dominant_coarse_state.isin(STATES) & m.day.ge(0) & ~m.coarse_label_tie)
    selected = []
    group_audit = []
    for (node,lib),gm in m[valid].groupby(['dominant_coarse_state','library_key']):
        ok = len(gm)>=3 and gm.n_cells.sum()>=30
        group_audit.append(dict(state=node,library_key=lib,n_metacells=len(gm),n_cells=int(gm.n_cells.sum()),
            days=';'.join(map(str,sorted(gm.day.unique()))),included=ok))
        if ok:
            selected.extend(gm.index)
    selected = sorted(selected)
    sm = m.iloc[selected].copy()
    write(sm,'selected_metacells.tsv.gz')
    write(pd.DataFrame(group_audit),'group_library_audit.tsv')
    write(sm.groupby('dominant_coarse_state').agg(n_metacells=('metacell_id','size'),n_cells=('n_cells','sum'),
        n_libraries=('library_key','nunique'),label=('dominant_coarse_label','first')).reset_index(),'groups.tsv')
    counts = {s:np.asarray(raw[m.index[m.index.isin(selected)&m.dominant_coarse_state.eq(s)]].sum(axis=0)).ravel().astype(np.int64) for s in STATES}
    assert all(v.sum()>0 for v in counts.values())
    keep = strict_pas(f)
    eligible = f[keep].copy()
    u = pd.read_csv(paths['utrs'],sep='\t')
    db = sqlite3.connect(f"file:{paths['catalog']}?mode=ro",uri=True)
    catalog = pd.read_sql_query('SELECT u.*,r.chrom,r.strand,r.start0,r.end0 FROM utrs u JOIN regions r USING(region_id)',db)
    u_groups = dict(tuple(u.groupby('gene_id')))
    c_groups = dict(tuple(catalog.groupby('gene_id')))
    candidates, all_rows, all_models = [], {}, {}
    for gene,gf in eligible.groupby('primary_gene_id',sort=True):
        if len(gf)<2 or gene not in u_groups:
            continue
        models = models_for_gene(u_groups[gene],c_groups.get(gene,catalog.iloc[:0]))
        rows = assign_pas(gf,models)
        all_rows[gene],all_models[gene] = rows,models
        gene_ix = gf.feature_index_1based.to_numpy(int)-1
        for model in models:
            pr = [r for r in rows if r['structure_id']==model['structure_id']]
            if len(pr)<2:
                continue
            ix = [r['feature_index_1based']-1 for r in pr]
            totals = [int(counts[s][ix].sum()) for s in STATES]
            coverage = [counts[s][ix].sum()/counts[s][gene_ix].sum() if counts[s][gene_ix].sum() else 0 for s in STATES]
            candidates.append(dict(gene_id=gene,gene_symbol=gf.primary_gene_symbol.iloc[0],strand=model['strand'],
                structure_id=model['structure_id'],n_pas=len(pr),n_gene_pas=len(gf),min_state_umi=min(totals),
                min_coverage=min(coverage),all_gene_pas_resolved=len(pr)==len(gf)))
    candidate = pd.DataFrame(candidates)
    write(candidate,'candidate_structure_inventory.tsv.gz')
    ready = candidate[(candidate.n_pas>=3)&(candidate.min_state_umi>=50)&(candidate.min_coverage>=.5)]
    chosen = []
    for strand in ['+','-']:
        ranking = ready[ready.strand.eq(strand)].sort_values(['min_state_umi','gene_id'],ascending=[False,True]).drop_duplicates('gene_id')
        chosen.extend(ranking.head(4).gene_id)
    # Fixed structural probes are added before any regulatory query, regardless of support/effect.
    for symbol in ['Scml2','Dnmt3a','Hmgn2']:
        found = eligible.loc[eligible.primary_gene_symbol.eq(symbol),'primary_gene_id'].unique()
        if len(found)==1 and found[0] in all_rows:
            chosen.append(found[0])
    chosen = list(dict.fromkeys(chosen))
    assert len(chosen)>=8
    genes = pd.DataFrame([dict(query_id=g,gene_id=g,gene_symbol=eligible.loc[eligible.primary_gene_id.eq(g),'primary_gene_symbol'].iloc[0],
        selection='count_support_strand_stratified' if g in ready.sort_values('min_state_umi').gene_id.values and g in chosen[:8] else 'predeclared_structural_probe') for g in chosen])
    write(genes,'genes.tsv')
    rows = [r for g in chosen for r in all_rows[g]]
    models = [r for g in chosen for r in all_models[g]]
    write(pd.DataFrame(rows),'gene_pas.tsv')
    write(pd.DataFrame(models),'structures.tsv')
    write(f[f.primary_gene_id.isin(chosen)].assign(pilot_qc_pass=keep[f.primary_gene_id.isin(chosen)]),'all_gene_pas_qc.tsv')
    cr,lr = [],[]
    for r in rows:
        j = r['feature_index_1based']-1
        for s in STATES:
            cr.append(dict(gene_id=r['gene_id'],atlas_pas_id=r['atlas_pas_id'],state=s,count=int(counts[s][j])))
        for (s,lib),gm in sm.groupby(['dominant_coarse_state','library_key']):
            lr.append(dict(gene_id=r['gene_id'],atlas_pas_id=r['atlas_pas_id'],state=s,library_key=lib,
                count=int(raw[gm.index,j].sum())))
    write(pd.DataFrame(cr),'pas_state_counts.tsv')
    write(pd.DataFrame(lr),'pas_library_state_counts.tsv.gz')
    query=[]
    for model in models:
        if sum(r['structure_id']==model['structure_id'] for r in rows)>=2:
            query.append(dict(query_id=model['structure_id'],gene_id=model['gene_id'],study_id=STUDY,region_id=model['region_id']))
    write(pd.DataFrame(query),'resource_query.tsv')
    # Export a small independent re-readable source matrix, not another normalization.
    cols = [r['feature_index_1based']-1 for r in rows]
    sparse.save_npz(OUT/'inputs/selected_metacell_by_gene_pas_counts.npz',raw[selected][:,cols])
    manifest.update(configuration=dict(study=STUDY,route='CIP',states=STATES,min_event_umi=50,min_coverage=.5,
        effect_tolerance=.1,selection='top4_by_min_state_event_UMI_per_strand_plus_named_probes_no_effect_selection',
        contrasts=[['S','R'],['R','E'],['E','T'],['R','M'],['S','T']],statistical_scope='pooled_descriptive'),
        audit=dict(matrix_shape=list(raw.shape),source_metacells=len(m),selected_metacells=len(selected),
            strict_pas=int(keep.sum()),candidate_structures=len(candidate),ready_ge3pas_structures=len(ready),
            selected_genes=len(chosen),query_structures=len(query),selected_pas=len(rows)))
    for k,p in paths.items():
        assert sha(p)==manifest[k]['sha256'], 'Input changed during preparation: '+k
    (OUT/'inputs/source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest['audit']),flush=True)
    print(genes.to_string(index=False))


if __name__=='__main__':
    main()
