"""Adapt one frozen analysis result set to coordinate-agnostic figure tables."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

PHY = ['species', 'contrast', 'category', 'gene_id', 'chrom', 'strand', 'start0', 'end0']
PAIR = ['species', 'contrast', 'gene_id', 'region_id_a', 'region_id_b']
ROUTES = {'mouse': ('CIP', 'WOT'), 'human': ('GSE118258', 'GSE150311')}


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def classify_availability(a, b, factors_a, factors_b, observed_a, observed_b):
    if not (a and b):
        return 'single_annotation'
    if not set(factors_a).intersection(factors_b):
        return 'identity_unmatched'
    return 'both_supported' if observed_a and observed_b else 'not_jointly_evaluable'


def first_observed(values, cutoff):
    # Missing intermediate states must never be promoted to a late onset.
    a, b = values[20], values[50]
    if not np.isfinite(np.r_[a, b]).all():
        return 'incomplete'
    da, db = a - a[0], b - b[0]
    if da[-1] * db[-1] <= 0 or min(abs(da[-1]), abs(db[-1])) < cutoff:
        return 'below'
    for i, label in enumerate('REP', 1):
        if da[i] * da[-1] > 0 and db[i] * db[-1] > 0 and min(abs(da[i]), abs(db[i])) >= cutoff:
            return label
    raise AssertionError('P must pass if the endpoint qualified')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--freeze', action='store_true', help='Create a new input lock, never overwrite one.')
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text())
    assert cfg['schema_version'] == 1
    assert cfg['coordinate_basis'] in ['original_pas', 'matched_annotation_pas']
    paths = {k: (args.config.parent / v).resolve() for k, v in cfg['sources'].items()}
    stamp = {k: {'path': str(p), 'sha256': digest(p)} for k, p in paths.items()}
    args.out.mkdir(parents=True, exist_ok=True)
    lockpath = args.out / 'input_lock.json'
    lock = {'schema_version': 1, 'result_set_id': cfg['result_set_id'],
            'coordinate_basis': cfg['coordinate_basis'], 'sources': stamp}
    if args.freeze:
        if lockpath.exists():
            raise FileExistsError('Use a new output directory for a new result set.')
        lockpath.write_text(json.dumps(lock, indent=2))
    assert json.loads(lockpath.read_text()) == lock, 'Inputs changed; do not mix runs. Freeze into a new output directory.'
    tables = args.out / 'tables'
    tables.mkdir(exist_ok=True)
    qa = {'result_set_id': cfg['result_set_id'], 'coordinate_basis': cfg['coordinate_basis'], 'checks': {}}

    def read(key, **kwargs):
        return pd.read_csv(paths[key], sep='\t', **kwargs)

    def save(name, frame):
        frame.to_csv(tables / (name + '.tsv.gz'), sep='\t', index=False)

    # Site-level availability: outer join physical coordinates, not an inner-join denominator.
    inv = read('inventory', usecols=PHY + ['route', 'factor', 'observed', 'quantified'])
    flags = inv.groupby(PHY + ['route'], observed=True)[['observed', 'quantified']].nunique()
    assert flags.to_numpy().max() == 1, 'Same physical site has conflicting observation flags.'
    physical = inv.groupby(PHY + ['route'], observed=True).agg(
        factors=('factor', frozenset), observed=('observed', 'first'), quantified=('quantified', 'first')).reset_index()
    availability = []
    for sp, (ra, rb) in ROUTES.items():
        a = physical[physical.route.eq(ra)].drop(columns='route')
        b = physical[physical.route.eq(rb)].drop(columns='route')
        z = a.merge(b, on=PHY, how='outer', suffixes=('_a', '_b'), indicator=True, validate='one_to_one')
        labels = []
        for _, r in z.iterrows():
            both = r['_merge'] == 'both'
            labels.append(classify_availability(r['_merge'] != 'right_only', r['_merge'] != 'left_only',
                r.factors_a if both else (), r.factors_b if both else (),
                bool(r.observed_a and r.quantified_a) if both else False,
                bool(r.observed_b and r.quantified_b) if both else False))
        z['availability'] = labels
        availability.append(z.groupby(['species', 'contrast', 'category', 'availability']).size().reset_index(name='n'))
        qa['checks'][sp + '_outer_physical_sites'] = len(z)
    save('availability', pd.concat(availability))
    del inv, physical

    exact = read('exact')
    common = exact[exact.both_observed].copy()
    stateflags = common.groupby(PHY)[['both_changed', 'same_direction']].nunique()
    assert stateflags.to_numpy().max() == 1
    phy = common.drop_duplicates(PHY).copy()
    phy['fate'] = np.where(~phy.both_changed, 'below', np.where(phy.same_direction, 'concordant', 'discordant'))
    save('supported_fate', phy.groupby(['species', 'contrast', 'category', 'fate']).size().reset_index(name='n'))
    availability_df = pd.concat(availability)
    supported = availability_df[availability_df.availability.eq('both_supported')].set_index(['species','contrast','category']).n
    assert supported.to_dict() == phy.groupby(['species','contrast','category']).size().to_dict()
    unitrows = []
    for key, q in common[common.both_changed].groupby(['species','contrast','category']):
        p = q.drop_duplicates(PHY)
        r = q.drop_duplicates(PAIR)
        unitrows.append(dict(zip(['species','contrast','category'], key)) | dict(
            site_n=len(p), region_n=len(r), site_c=p.same_direction.mean(), region_c=r.same_direction.mean()))
    units = pd.DataFrame(unitrows)
    save('units', units)
    stats = read('statistics')
    primary = stats[stats.unit.eq('region_pair')]
    for r in units.itertuples():
        ref = primary[primary.species.eq(r.species)&primary.contrast.eq(r.contrast)&primary.category.eq(r.category)].iloc[0]
        assert r.region_n == ref.n and np.isclose(r.region_c, ref.cobs)
    save('statistics', stats)
    gf = read('classification')
    save('classification', gf.groupby(['species','route','contrast','category','classification']).size().reset_index(name='n'))
    partial = read('partial')
    for r in partial.itertuples():
        q = gf[gf.route.eq(r.route)&gf.contrast.eq(r.contrast)&gf.category.eq(r.category)]
        q = q[q.n_evaluable.ge(2)&q.n_changed.gt(0)]
        assert len(q)==r.changed_multisite_pairs and int(q.has_partial.sum())==r.partial_including_mixed
    save('partial', partial)
    del gf
    cp = read('crosspartial')
    cp = cp[cp.exact_partial_recurrence].copy()
    save('crosspartial', cp)
    save('exact_partial', cp.groupby(['species','contrast','category']).agg(pairs=('factor','size'),genes=('gene_id','nunique')).reset_index())
    unique_cp = cp.drop_duplicates(['species','category','gene_id','factor'])
    save('gene_contribution', unique_cp.groupby(['species','category','gene_id','gene_symbol_a']).size().reset_index(name='pairs'))

    # Use exact main-figure row/mask identities; never join state curves by gene alone.
    rows, hv = read('heat_rows'), read('heat_values')
    recurrent = exact[exact.recurrent_strict].copy()
    annotation = recurrent.groupby(PAIR).category.agg(lambda x:';'.join(sorted(set(x)))).reset_index(name='contrast_categories')
    heat = rows.assign(contrast=rows.contrasts.str.split(';')).explode('contrast').merge(annotation, on=PAIR, validate='one_to_one')
    assert len(heat) == sum(rows.contrasts.str.split(';').str.len())
    effects = recurrent.drop_duplicates(PAIR)[PAIR+['effect20_a','effect50_a','effect20_b','effect50_b']]
    heat = heat.merge(effects, on=PAIR, validate='one_to_one')
    heat['category_order']=heat.contrast_categories.map({'miRNA':0,'RBP':1,'RBP;miRNA':2})
    heat=heat.sort_values(['species','contrast','category_order','effect20_a','gene_id','row_id'])
    heat['display_row']=heat.groupby(['species','contrast']).cumcount()+1
    for r in heat.itertuples():
        ref = 'S' if r.contrast=='SP' else 'M'
        for side in ['a','b']:
            for threshold in [20,50]:
                q=hv[hv.row_id.eq(r.row_id)&hv.side.eq(side)&hv.threshold.eq(threshold)].set_index('state')
                delta=q.loc['P','inclusion']-q.loc[ref,'inclusion']
                assert np.isclose(delta,getattr(r,f'effect{threshold}_{side}'))
    save('heat_rows',heat)
    save('heat_values',hv)
    save('state_inventory',read('state_inventory'))
    coverage=hv[hv.threshold.eq(20)].groupby(['species','route','state']).agg(
        available=('inclusion','count'),total=('row_id','size'),min_libraries=('n_libraries','min'),max_libraries=('n_libraries','max')).reset_index()
    save('coverage',coverage)
    thresholdrows=[]
    for (row_id,route),g in hv.groupby(['row_id','route']):
        v={t:g[g.threshold.eq(t)].set_index('state').reindex(list('SREP')).inclusion.to_numpy() for t in [20,50]}
        for cutoff in [.05,.1,.15,.2]:
            thresholdrows.append(dict(row_id=row_id,route=route,species=g.species.iloc[0],gene=g.gene.iloc[0],
                cutoff=cutoff,first_observed=first_observed(v,cutoff)))
    save('thresholds',pd.DataFrame(thresholdrows))
    qa['checks']['heat_rows_by_comparison']=heat.groupby(['species','contrast']).size().to_dict().__str__()

    cases = read('case_profiles')[['gene','route','pattern_id']].drop_duplicates()
    assert not cases.duplicated(['gene','route']).any()
    lib=read('library_profiles').merge(cases,on='pattern_id',validate='many_to_one')
    save('case_libraries',lib)
    means=lib.groupby(['gene','route','threshold','state']).agg(inclusion=('inclusion','mean'),libraries=('library','nunique')).reset_index()
    ref=read('case_profiles')
    test=means[means.threshold.eq(20)].merge(ref,on=['gene','route','state','threshold'],suffixes=('_new','_ref'),validate='one_to_one')
    assert len(test)==len(ref) and np.allclose(test.inclusion_new,test.inclusion_ref)
    save('case_means',means)
    ct=[]
    for (gene,route),g in means.groupby(['gene','route']):
        v={t:g[g.threshold.eq(t)].set_index('state').reindex(list('SREP')).inclusion.to_numpy() for t in [20,50]}
        for cutoff in [.05,.1,.15,.2]:
            ct.append(dict(gene=gene,route=route,cutoff=cutoff,first_observed=first_observed(v,cutoff)))
    save('case_thresholds',pd.DataFrame(ct))
    save('case_sensitivity',read('sensitivity').merge(cases,on='pattern_id',validate='many_to_one'))
    occupancy=read('occupancy')
    for r in occupancy.groupby(['gene','route','site']).included:
        assert set(r[1]) <= {0,1}
    save('occupancy',occupancy)
    save('case_sites',read('case_sites'))
    expr=read('expression')
    save('expression',expr[expr.gene_symbol.str.upper().eq('PCBP1')])
    save('candidates',read('candidates'))
    save('glyphs',read('glyphs'))
    save('datasets',read('datasets'))
    for key in ['eclip','sequence']:
        (tables/(key+'.json')).write_text(paths[key].read_text())
    # The figure contract describes provenance; it does not relabel old calculations as new ones.
    manifest={'schema_version':1,'result_set_id':cfg['result_set_id'],'coordinate_basis':cfg['coordinate_basis'],
              'files':{p.name:digest(p) for p in sorted(tables.iterdir())},
              'input_lock_sha256':digest(lockpath),'adapter_sha256':digest(Path(__file__))}
    (args.out/'bundle.json').write_text(json.dumps(manifest,indent=2))
    qa['status']='PASS'
    qa['limitations']=['No Figure5-specific state-membership exclusion rerun.',
                      'Threshold-specific gene-label permutation tests were not added.',
                      'PAS-coordinate version is producer-declared; source hashes and cross-table identities are verified.']
    (args.out/'data_QA.json').write_text(json.dumps(qa,indent=2))
    print(json.dumps(qa,indent=2))


if __name__=='__main__':
    main()
