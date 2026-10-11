"""Candidate RBP expression and predicted target-site retention, by study."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
KEY = ['species', 'study_id', 'route', 'factor']


def save(frame, name):
    frame.to_csv(HERE / name, sep='\t', index=False, na_rep='NA', float_format='%.12g',
        compression={'method': 'gzip', 'mtime': 0} if name.endswith('.gz') else None)


def library_rows(meta, counts, totals, factors, study, species):
    rows = []
    for (route, library, state), indices in meta.groupby(['route', 'library_key', 'node']).groups.items():
        ix = np.array(list(indices)); part = meta.loc[ix]
        assert part.experimental_time.nunique(dropna=False) == 1
        for j, factor in enumerate(factors):
            c = counts[ix, j]; available = (c >= 0).all()
            cpm = c / totals[ix] * 1e6 if available else np.full(len(ix), np.nan)
            rows.append(dict(species=species, study_id=study, route=route, factor=factor,
                library_key=library, state=state, time=str(part.experimental_time.iloc[0]),
                n_metacells=len(ix), n_cells=int(part.n_cells.sum()),
                qualified=len(ix) >= 3 and part.n_cells.sum() >= 30,
                mapped=available, count=int(c.sum()) if available else np.nan,
                total_gex=int(totals[ix].sum()),
                pooled_CPM=float(c.sum()/totals[ix].sum()*1e6) if available else np.nan,
                mean_log2CPM1=float(np.log2(cpm+1).mean()), mean_CPM=float(cpm.mean()),
                detected_metacells=int((c > 0).sum()) if available else np.nan,
                detected_fraction=float((c > 0).mean()) if available else np.nan))
    return pd.DataFrame(rows)


def state_rows(libraries):
    rows = []
    for key, all_rows in libraries.groupby(KEY):
        for state in 'SREPM':
            raw = all_rows[all_rows.state.eq(state)]
            q = raw[raw.qualified & raw.mapped]
            rec = dict(zip(KEY, key), state=state, n_libraries=len(q),
                n_qualified_libraries=int(raw.qualified.sum()), n_metacells=int(q.n_metacells.sum()),
                n_cells=int(q.n_cells.sum()), n_times=q.time.nunique(),
                status='available' if len(q) else 'absent_state' if raw.empty else 'unmapped_or_undercovered')
            for column in ['mean_log2CPM1', 'pooled_CPM', 'detected_fraction']:
                rec[column] = q[column].mean()
            rec.update(log_min=q.mean_log2CPM1.min(), log_max=q.mean_log2CPM1.max(),
                time_equal_log2CPM1=q.groupby('time').mean_log2CPM1.mean().mean(),
                expressed_any_2CPM=bool(q.pooled_CPM.ge(2).any()))
            rows.append(rec)
    return pd.DataFrame(rows)


def contrast_rows(libraries, states):
    rows, pair_rows = [], []
    for key, raw in libraries.groupby(KEY):
        q = raw[raw.qualified & raw.mapped]
        s = states.set_index(KEY).loc[[key]].set_index('state')
        for contrast, reference in [('SP', 'S'), ('MP', 'M')]:
            rec = dict(zip(KEY, key), contrast=contrast, target='P', reference=reference,
                n_reference_libraries=int(s.loc[reference, 'n_libraries']),
                n_target_libraries=int(s.loc['P', 'n_libraries']), n_pairs=0,
                delta_log2CPM1=np.nan, time_equal_delta=np.nan,
                target_log2CPM1=np.nan, reference_log2CPM1=np.nan,
                target_CPM=np.nan, reference_CPM=np.nan,
                expression_in_either_state=False, expression_in_both_states=False,
                status='unavailable')
            if contrast == 'SP' and s.loc[['S', 'P'], 'n_libraries'].gt(0).all():
                rec.update(status='descriptive_unpaired',
                    delta_log2CPM1=s.loc['P', 'mean_log2CPM1']-s.loc['S', 'mean_log2CPM1'],
                    time_equal_delta=s.loc['P', 'time_equal_log2CPM1']-s.loc['S', 'time_equal_log2CPM1'])
                target, ref = q[q.state.eq('P')], q[q.state.eq('S')]
            elif contrast == 'MP':
                pairs = q[q.state.eq('P')].merge(q[q.state.eq('M')], on='library_key', suffixes=('_P', '_M'), validate='one_to_one')
                if pairs.empty:
                    rows.append(rec); continue
                pairs['delta'] = pairs.mean_log2CPM1_P-pairs.mean_log2CPM1_M
                rec.update(status='descriptive_paired', n_pairs=len(pairs),
                    delta_log2CPM1=pairs.delta.mean(),
                    time_equal_delta=pairs.groupby('time_P').delta.mean().mean(),
                    n_reference_libraries=len(pairs), n_target_libraries=len(pairs))
                for r in pairs.itertuples():
                    pair_rows.append(dict(zip(KEY, key), library_key=r.library_key, time=r.time_P,
                        target_log2CPM1=r.mean_log2CPM1_P, reference_log2CPM1=r.mean_log2CPM1_M,
                        delta_log2CPM1=r.delta, target_CPM=r.pooled_CPM_P, reference_CPM=r.pooled_CPM_M))
                target = q[q.state.eq('P') & q.library_key.isin(pairs.library_key)]
                ref = q[q.state.eq('M') & q.library_key.isin(pairs.library_key)]
            else:
                rows.append(rec); continue
            rec.update(target_log2CPM1=target.mean_log2CPM1.mean(), reference_log2CPM1=ref.mean_log2CPM1.mean(),
                target_CPM=target.pooled_CPM.mean(), reference_CPM=ref.pooled_CPM.mean(),
                expression_in_either_state=bool(target.pooled_CPM.ge(2).any() or ref.pooled_CPM.ge(2).any()),
                expression_in_both_states=bool(target.pooled_CPM.ge(2).any() and ref.pooled_CPM.ge(2).any()))
            rows.append(rec)
    return pd.DataFrame(rows), pd.DataFrame(pair_rows)


def site_rows(sites):
    rows = []
    for key, q in sites.groupby(['species', 'factor', 'evidence', 'contrast']):
        for suffix in ['a', 'b']:
            for threshold in [20, 50]:
                e = q[f'effect{threshold}_{suffix}']
                rows.append(dict(zip(['species', 'factor', 'evidence', 'contrast'], key),
                    route=q[f'route_{suffix}'].iloc[0], threshold=threshold,
                    n_target_genes=q.gene_id.nunique(), n_sites=len(q),
                    n_sites_increased=int(e.gt(0).sum()), n_sites_decreased=int(e.lt(0).sum()),
                    weighted_site_change=e.sum(), mean_site_inclusion_delta=e.mean(),
                    mean_gene_inclusion_delta=q.assign(effect=e).groupby('gene_id').effect.mean().mean()))
    return pd.DataFrame(rows)


def main():
    archive = np.load(HERE / 'expression_inputs.npz', allow_pickle=False)
    config = json.loads((HERE / 'config.json').read_text())
    candidates = pd.read_csv(HERE / 'candidates.tsv', sep='\t')
    libraries = []
    for study, (species, folder, route) in config['studies'].items():
        meta = pd.DataFrame(json.loads(str(archive[study+'_meta_json'])))
        counts, total = archive[study+'_counts'], archive[study+'_total_gex']
        assert (counts >= -1).all() and (total > 0).all()
        libraries.append(library_rows(meta, counts, total, archive[study+'_factors'], study, species))
    libraries = pd.concat(libraries, ignore_index=True)
    libraries = libraries.merge(candidates, on=['species', 'factor'], validate='many_to_one')
    states = state_rows(libraries)
    contrasts, pairs = contrast_rows(libraries, states)
    for f in [states, contrasts, pairs]:
        f['evidence'] = f.set_index(['species', 'factor']).index.map(candidates.set_index(['species', 'factor']).evidence)
    site = site_rows(pd.read_csv(HERE / 'candidate_sites.tsv.gz', sep='\t'))
    linked = contrasts.merge(site[site.threshold.eq(20)].drop(columns=['threshold', 'evidence']),
        on=['species', 'factor', 'route', 'contrast'], how='left', validate='one_to_one')
    cross = []
    for key, q in linked.groupby(['species', 'factor', 'evidence', 'contrast']):
        q = q.sort_values('study_id'); assert len(q) == 2
        a, b = q.iloc[0], q.iloc[1]
        da, db = a.delta_log2CPM1, b.delta_log2CPM1
        status = ('unavailable' if not np.isfinite(da+db) else 'both_increased' if min(da,db)>1e-12
            else 'both_decreased' if max(da,db)<-1e-12 else 'zero_or_tied' if min(abs(da),abs(db))<=1e-12 else 'opposite')
        cross.append(dict(zip(['species', 'factor', 'evidence', 'contrast'], key),
            study_a=a.study_id, study_b=b.study_id, effect_a=da, effect_b=db,
            time_effect_a=a.time_equal_delta, time_effect_b=b.time_equal_delta,
            direction=status, expressed_in_both_studies=bool(a.expression_in_either_state and b.expression_in_either_state),
            supported_both_endpoints_both_studies=bool(a.expression_in_both_states and b.expression_in_both_states),
            has_recurrent_sites=bool(pd.notna(a.n_sites) and a.n_sites > 0),
            n_target_genes=a.n_target_genes, n_sites=a.n_sites,
            sign_stable_time_weighting=bool(np.isfinite(da+db+a.time_equal_delta+b.time_equal_delta)
                and np.sign(da)==np.sign(a.time_equal_delta) and np.sign(db)==np.sign(b.time_equal_delta))))
    cross = pd.DataFrame(cross)
    mapping = pd.DataFrame(json.loads(str(archive['mapping_json'])))
    for f, name in [(libraries, 'library_expression.tsv.gz'), (states, 'state_expression.tsv'),
        (contrasts, 'expression_contrasts.tsv'), (pairs, 'paired_library_effects.tsv'),
        (site, 'site_retention.tsv'), (linked, 'expression_site_links.tsv'),
        (cross, 'cross_study_summary.tsv'), (mapping, 'gene_mapping.tsv')]:
        save(f, name)
    (HERE / 'remote_input_sha256.json').write_text(str(archive['provenance_json'])+'\n')
    summary = dict(candidates=candidates.groupby(['species','evidence']).size().to_string(),
        mapping_status=mapping.groupby(['species','status']).factor.nunique().to_string(),
        cross_study=cross[cross.has_recurrent_sites].groupby(['species','evidence','contrast','direction']).size().to_string(),
        inputs_sha256={name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['expression_inputs.npz','candidate_sites.tsv.gz','candidates.tsv','config.json']})
    (HERE/'analysis_audit.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    main()
