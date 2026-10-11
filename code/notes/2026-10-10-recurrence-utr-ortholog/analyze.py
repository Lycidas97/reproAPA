"""Tandem recurrence and fixed-margin ortholog direction correspondence."""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TABLES = ROOT / 'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables'
sys.path.insert(0, str(ROOT / 'notes/2026-10-10-directional-recurrence'))
from analyze import observed_members, prepare_units, simulate, seed, bh, digest

LABELS = ['proximal', 'distal', 'nonrecurrent']


def save(frame, name):
    compression = {'method': 'gzip', 'mtime': 0} if name.endswith('.gz') else None
    frame.to_csv(HERE / name, sep='\t', index=False, na_rep='NA',
                 float_format='%.17g', compression=compression)


def quantiles(values, n):
    assert values.notna().all()
    return pd.qcut(values, n, labels=False, duplicates='drop').fillna(0).astype(int)


def eligible_tandem(frame, eligibility):
    joined = frame.merge(eligibility[['species', 'study_id', 'gene_id', 'tandem_utr_compatible']],
                         left_on=['species', 'unit', 'gene_id'],
                         right_on=['species', 'study_id', 'gene_id'], validate='many_to_one')
    assert len(joined) == len(frame) and joined.tandem_utr_compatible.notna().all()
    gates = joined.groupby('gene_id').tandem_utr_compatible.all()
    return set(gates.index[gates])


def common_pairs(frame, orthologs):
    measured = {}
    for species, part in frame.groupby('species'):
        counts = part.groupby('gene_id').unit.nunique()
        measured[species] = set(counts.index[counts.ge(2)])
    return orthologs[orthologs.human_gene_id.isin(measured['human']) &
                     orthologs.mouse_gene_id.isin(measured['mouse'])].copy().reset_index(drop=True)


def gene_features(frame, ids):
    rows = []
    for gene, part in frame[frame.gene_id.isin(ids)].groupby('gene_id', sort=True):
        part = part.sort_values('unit')
        rows.append(dict(gene_id=gene, measured_pattern=';'.join(
            part.unit + ':' + part.series_pattern),
            proximal_capacity=int(part.can_negative.sum()), distal_capacity=int(part.can_positive.sum()),
            PAS_number=float(part.PAS_number.median()), mean_GEX=float(part.mean_GEX.median()),
            APA_coverage=float(part.APA_coverage.median()), baseline=float(part.baseline.median())))
    out = pd.DataFrame(rows).set_index('gene_id').reindex(ids)
    assert out.notna().all().all()
    for name, n in [('PAS_number', 3), ('mean_GEX', 3), ('APA_coverage', 3), ('baseline', 5)]:
        out[name + '_bin'] = quantiles(out[name], n)
    return out.reset_index()


def class_labels(frame, ids):
    roster = observed_members(frame).set_index('gene_id').direction
    return roster.reindex(ids).map({'proximal': 0, 'distal': 1}).fillna(2).astype(np.int8).to_numpy()


def partitions(features, model):
    cols = ['measured_pattern', 'proximal_capacity', 'distal_capacity', 'PAS_number_bin']
    if model == 'full_assay':
        cols += ['mean_GEX_bin', 'APA_coverage_bin', 'baseline_bin']
    return [np.asarray(ii, int) for ii in features.groupby(cols, sort=True).indices.values()]


def class_probabilities(labels, cells):
    p = np.zeros((len(labels), 3))
    for ii in cells:
        p[ii] = np.bincount(labels[ii], minlength=3) / len(ii)
    return p


def draw_classes(labels, cells, size, random):
    out = np.broadcast_to(labels, (size, len(labels))).copy()
    for ii in cells:
        values = labels[ii]
        if np.unique(values).size > 1:
            out[:, ii] = random.permuted(np.broadcast_to(values, (size, len(ii))), axis=1)
        for label in range(3):
            assert (out[:, ii] == label).sum(1).tolist() == [int((values == label).sum())] * size
    return out


def cross_draws(h, m, hg, mg, random, n=10000):
    result = np.zeros((n, 3, 3), int)
    for start in range(0, n, 250):
        size = min(250, n-start)
        hd = draw_classes(h, hg, size, random)
        md = draw_classes(m, mg, size, random)
        for a in range(3):
            for b in range(3):
                result[start:start+size, a, b] = ((hd == a) & (md == b)).sum(1)
        assert (result[start:start+size].sum((1, 2)) == len(h)).all()
    return result


def null_row(observed, values, expected):
    sd = float(values.std(ddof=1))
    assert abs(values.mean()-expected) <= 6*sd/np.sqrt(len(values)) + 1e-9
    tail = int((values >= observed).sum())
    return dict(observed=int(observed), analytic_expected=float(expected), null_mean=float(values.mean()),
        null_sd=sd, null_p025=float(np.quantile(values, .025)), null_p975=float(np.quantile(values, .975)),
        excess=float(observed-expected), fold=float(observed/expected) if expected else np.nan,
        p_enrichment=(1+tail)/(1+len(values)), tail_count=tail, permutations=len(values),
        inferentially_mobile=sd > 0)


def main():
    paths = [HERE/'protocol.json', TABLES/'recurrence_controls/recurrence_frozen_covariates.tsv.gz',
             TABLES/'tandem/tandem_utr_eligibility.tsv.gz', TABLES/'tandem/tandem_utr_effects.tsv.gz',
             TABLES/'cross_overlap/ortholog_matched_features.tsv.gz',
             ROOT/'notes/2026-10-10-directional-recurrence/observed_members.tsv.gz',
             ROOT/'notes/2026-10-10-directional-recurrence/analyze.py',
             ROOT/'notes/2026-10-02-null-model-controls/recurrence.py',
             ROOT/'notes/2026-10-02-null-model-controls/common.py',
             ROOT/'analysis_code/pas_rerun/downstream_runtime.py']
    hashes = {str(p.relative_to(ROOT)): digest(p) for p in paths}
    f = pd.read_csv(paths[1], sep='\t', float_precision='round_trip')
    eligibility = pd.read_csv(paths[2], sep='\t')
    effects = pd.read_csv(paths[3], sep='\t', float_precision='round_trip')
    effects = effects[effects.setting.eq('supported50_p08')].copy()
    assert np.allclose(effects.delta_nt, effects.delta*effects.native_axis_span_nt, atol=1e-9, rtol=0)
    orth = pd.read_csv(paths[4], sep='\t')[['human_gene_id', 'mouse_gene_id', 'human_gene_symbol',
                                          'mouse_gene_symbol', 'is_high_confidence', 'homology_type']]
    assert len(orth) == 16092 and orth.is_high_confidence.eq(1).all()
    assert orth.homology_type.eq('ortholog_one2one').all()
    assert orth.human_gene_id.is_unique and orth.mouse_gene_id.is_unique
    original = pd.read_csv(paths[5], sep='\t')
    tandem, length, subset_members, margins, inventories, archives = [], [], [], [], [], {}
    for (sp, co, threshold), part in f.groupby(['species', 'contrast', 'threshold'], sort=True):
        context = dict(species=sp, contrast=co, threshold=threshold)
        observed = observed_members(part)
        old = original[original.species.eq(sp)&original.contrast.eq(co)&original.threshold.eq(threshold)]
        assert set(zip(observed.gene_id, observed.direction)) == set(zip(old.gene_id, old.direction))
        ids = eligible_tandem(part, eligibility)
        sub = part[part.gene_id.isin(ids)].copy()
        roster = observed_members(sub)
        genes, units, expected, cells = prepare_units(sub)
        key = '|'.join(map(str, ('tandem', sp, co, threshold)))
        run_seed = seed('utr_direction', sp, co, threshold)
        null = simulate(units, len(genes), 10000, np.random.default_rng(run_seed))
        archives[key] = null
        denominator = int(sub.groupby('gene_id').unit.nunique().ge(2).sum())
        inventories.append(context | dict(all_measured_genes=part.gene_id.nunique(),
            recurrence_measurable_genes=int(part.groupby('gene_id').unit.nunique().ge(2).sum()),
            tandem_inventory_genes=len(ids), tandem_recurrence_measurable_genes=denominator,
            original_recurrent_genes=len(observed), tandem_recurrent_genes=len(roster)))
        subset_members.append(roster.assign(**context))
        margins.append(cells.assign(**context))
        for i, direction in enumerate(LABELS[:2]):
            count = int(roster.direction.eq(direction).sum())
            tandem.append(context | dict(direction=direction, measurable_genes=denominator, seed=run_seed,
                                           **null_row(count, null[:, i], expected[direction])))
        design = 'true_source' if co == 'P_S' else 'same_library'
        z = effects[effects.species.eq(sp)&effects.contrast.eq(co)&effects.design.eq(design)&
                    effects.study_id.isin(part.unit.unique())&effects.gene_id.isin(roster.gene_id)].merge(
            roster[['gene_id', 'direction']], on='gene_id', validate='many_to_one')
        contributors = part[['unit', 'gene_id', 'series_pattern', 'effect_mean']].rename(columns={'unit': 'study_id'})
        z = z.merge(contributors, on=['study_id', 'gene_id'], validate='many_to_one')
        z = z[[r in pattern.split(';') for r, pattern in zip(z.route, z.series_pattern)]]
        recovered = z.groupby(['study_id', 'gene_id']).agg(delta=('delta', 'mean'), original=('effect_mean', 'first'))
        assert np.allclose(recovered.delta, recovered.original, atol=1e-12, rtol=0)
        assert z.tandem_utr_compatible.all() and np.sign(z.delta_nt).eq(np.sign(z.delta)).all()
        for (sid, route, direction), q in z.groupby(['study_id', 'route', 'direction']):
            length.append(context | dict(study_id=sid, route=route, direction=direction, level='route',
                n_genes=len(q), mean_DWUI=q.delta.mean(), median_DWUI=q.delta.median(),
                mean_delta_nt=q.delta_nt.mean(), median_delta_nt=q.delta_nt.median(),
                n_opposing=int((q.delta*np.where(direction=='distal', 1, -1) < 0).sum())))
        study = z.groupby(['direction', 'gene_id', 'study_id'])[['delta', 'delta_nt']].mean()
        gene = study.groupby(['direction', 'gene_id']).mean().reset_index()
        for direction, q in gene.groupby('direction'):
            length.append(context | dict(study_id='study_equal', route='gene_equal', direction=direction,
                level='gene_equal_study_equal', n_genes=len(q), mean_DWUI=q.delta.mean(),
                median_DWUI=q.delta.median(), mean_delta_nt=q.delta_nt.mean(), median_delta_nt=q.delta_nt.median(),
                n_opposing=int((q.delta*np.where(direction=='distal', 1, -1) < 0).sum())))
        print(key, 'recurrent', len(roster), 'eligible', denominator, flush=True)
    t = pd.DataFrame(tandem)
    t['q_eight_per_threshold'] = t.groupby('threshold').p_enrichment.transform(bh)
    t['q_joint_16'] = bh(t.p_enrichment)
    save(t, 'tandem_directional_recurrence.tsv')
    save(pd.DataFrame(length), 'tandem_recurrent_length.tsv')
    save(pd.DataFrame(inventories), 'tandem_inventory.tsv')
    save(pd.concat(subset_members, ignore_index=True), 'tandem_members.tsv.gz')
    save(pd.concat(margins, ignore_index=True), 'tandem_stratum_margins.tsv.gz')
    cross, paired, diagnostics, features = [], [], [], []
    for (co, threshold), part in f.groupby(['contrast', 'threshold'], sort=True):
        pairs = common_pairs(part, orth)
        context = dict(contrast=co, threshold=threshold, measured_orthologs=len(pairs))
        labels, feats = {}, {}
        for sp in ['human', 'mouse']:
            ids = pairs[sp+'_gene_id'].tolist()
            sf = part[part.species.eq(sp)]
            feats[sp] = gene_features(sf, ids)
            labels[sp] = class_labels(sf, ids)
            features.append(feats[sp].assign(species=sp, contrast=co, threshold=threshold))
            pairs[sp+'_direction'] = np.asarray(LABELS)[labels[sp]]
        paired.append(pairs.assign(contrast=co, threshold=threshold))
        h, m = labels['human'], labels['mouse']
        for model in ['measurement_PAS_capacity', 'full_assay']:
            cells = {sp: partitions(feats[sp], model) for sp in ['human', 'mouse']}
            probabilities = {sp: class_probabilities(labels[sp], cells[sp]) for sp in ['human', 'mouse']}
            analytic = probabilities['human'].T @ probabilities['mouse']
            run_seed = seed('ortholog_direction', co, threshold, model)
            null = cross_draws(h, m, cells['human'], cells['mouse'], np.random.default_rng(run_seed))
            archives['|'.join(map(str, ('ortholog', co, threshold, model)))] = null
            for sp in ['human', 'mouse']:
                for j, ii in enumerate(cells[sp]):
                    counts = np.bincount(labels[sp][ii], minlength=3)
                    diagnostics.append(context | dict(model=model, species=sp, stratum=j, n=len(ii),
                        proximal=int(counts[0]), distal=int(counts[1]), nonrecurrent=int(counts[2]),
                        expected_changed_fraction=1-float((counts/len(ii))@(counts/len(ii)))))
            for a, hd in enumerate(LABELS):
                for b, md in enumerate(LABELS):
                    row = context | dict(model=model, human_direction=hd, mouse_direction=md, seed=run_seed,
                        **null_row(int(((h==a)&(m==b)).sum()), null[:, a, b], analytic[a, b]))
                    row['tested'] = a < 2 and b < 2
                    cross.append(row)
            print('ortholog', co, threshold, model, 'pairs', len(pairs), flush=True)
    c = pd.DataFrame(cross)
    c.loc[~c.tested, 'p_enrichment'] = np.nan
    c['q_eight_per_threshold_model'] = np.nan
    c['q_joint_16_per_model'] = np.nan
    for _, ii in c[c.tested].groupby(['threshold', 'model']).groups.items():
        c.loc[ii, 'q_eight_per_threshold_model'] = bh(c.loc[ii, 'p_enrichment'])
    for _, ii in c[c.tested].groupby('model').groups.items():
        c.loc[ii, 'q_joint_16_per_model'] = bh(c.loc[ii, 'p_enrichment'])
    save(c, 'ortholog_3x3.tsv')
    save(pd.concat(paired, ignore_index=True), 'measured_ortholog_members.tsv.gz')
    save(pd.concat(features, ignore_index=True), 'ortholog_assay_features.tsv.gz')
    save(pd.DataFrame(diagnostics), 'ortholog_stratum_diagnostics.tsv.gz')
    np.savez_compressed(HERE/'null_draws.npz', **archives)
    assert all(digest(ROOT/p) == value for p, value in hashes.items())
    output_names = ['tandem_directional_recurrence.tsv', 'tandem_recurrent_length.tsv', 'tandem_inventory.tsv',
        'tandem_members.tsv.gz', 'tandem_stratum_margins.tsv.gz', 'ortholog_3x3.tsv',
        'measured_ortholog_members.tsv.gz', 'ortholog_assay_features.tsv.gz',
        'ortholog_stratum_diagnostics.tsv.gz', 'null_draws.npz']
    (HERE/'manifest.json').write_text(json.dumps(dict(status='PASS', input_sha256=hashes,
        output_sha256={name: digest(HERE/name) for name in output_names}, script_sha256=digest(Path(__file__)),
        python_executable=sys.executable, python=sys.version, numpy=np.__version__, pandas=pd.__version__,
        original_directional_rosters_verified=True, analytic_null_means_verified=True,
        all_draws_preserve_stratum_margins=True, full_ortholog_identity_count=len(orth)), indent=2)+'\n')


if __name__ == '__main__':
    main()
