"""Independent labels, covariate bins, analytical null and seeded-prefix replay.

Uses the review's independently rebuilt study covariates. No production helper
is imported. Full null tails/BH are independently audited by root_verify.py.
"""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASE = HERE.parent
PROD = BASE / 'analysis'
OLD = ROOT / 'notes/2026-10-10-ortholog-threshold-scan'
LABELS = ['proximal', 'distal', 'nonrecurrent']
PREFIX = 1000  # Four complete legacy batches; chosen before looking at draws.


def read(path):
    return pd.read_csv(path, sep='\t', float_precision='round_trip')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def close(actual, expected, name, atol=2e-10):
    np.testing.assert_allclose(actual, expected, atol=atol, rtol=0,
                               equal_nan=True, err_msg=name)


def summarize(source, gene_ids, threshold):
    rows, direction, support = [], [], []
    for gene in gene_ids:
        part = source[source.gene_id.eq(gene)].sort_values('unit')
        assert part.unit.is_unique and len(part) >= 2
        hits_proximal = part.effect_max.le(-threshold)
        hits_distal = part.effect_min.ge(threshold)
        assert not (hits_proximal & hits_distal).any()
        p, d = int(hits_proximal.sum()), int(hits_distal.sum())
        assert not (p >= 2 and d >= 2)
        direction.append(0 if p >= 2 else 1 if d >= 2 else 2)
        support.append((p, d, len(part), min(p, d) > 0))
        rows.append(dict(gene_id=gene,
            measured_pattern=';'.join(part.unit + ':' + part.series_pattern),
            proximal_capacity=int(part.capacity_negative.ge(threshold-1e-9).sum()),
            distal_capacity=int(part.capacity_positive.ge(threshold-1e-9).sum()),
            **{col:float(part[col].median()) for col in
               ['PAS_number', 'mean_GEX', 'APA_coverage', 'baseline']}))
    features = pd.DataFrame(rows)
    for col, bins in [('PAS_number',3), ('mean_GEX',3), ('APA_coverage',3), ('baseline',5)]:
        features[col+'_bin'] = pd.qcut(features[col], bins, labels=False,
                                      duplicates='drop').fillna(0).astype(int)
    return features, np.array(direction, dtype=np.int8), support


def strata(features, model):
    columns = ['measured_pattern','proximal_capacity','distal_capacity','PAS_number_bin']
    if model == 'full_assay':
        columns += ['mean_GEX_bin','APA_coverage_bin','baseline_bin']
    return [np.array(indices, dtype=int) for indices in
            features.groupby(columns, sort=True).indices.values()]


def probabilities_and_mobility(features, labels, groups):
    probabilities = np.zeros((len(labels),3))
    diagnostics = []
    for j, indices in enumerate(groups):
        counts = np.bincount(labels[indices], minlength=3)
        probabilities[indices] = counts / len(indices)
        diagnostics.append(dict(stratum=j,n=len(indices),proximal=int(counts[0]),
            distal=int(counts[1]),nonrecurrent=int(counts[2]),
            singleton=len(indices)==1,homogeneous=np.count_nonzero(counts)==1,
            expected_changed_fraction=1-float(np.sum((counts/len(indices))**2)),
            gene_ids=';'.join(features.iloc[indices].gene_id)))
    diagnostics = pd.DataFrame(diagnostics)
    mobility = dict(measured_genes=len(labels),n_strata=len(groups),
        n_singleton_strata=int(diagnostics.singleton.sum()),
        largest_stratum=int(diagnostics.n.max()),
        singleton_gene_fraction=float(diagnostics.loc[diagnostics.singleton,'n'].sum()/len(labels)),
        frozen_label_fraction=float(diagnostics.loc[diagnostics.homogeneous,'n'].sum()/len(labels)),
        expected_changed_fraction=float(np.sum(diagnostics.n*diagnostics.expected_changed_fraction)/len(labels)))
    return probabilities, diagnostics, mobility


def draw_species(labels, groups, random, n):
    draws = np.tile(labels, (n,1))
    for indices in groups:
        fixed = labels[indices]
        if len(np.unique(fixed)) > 1:
            draws[:,indices] = random.permuted(np.tile(fixed,(n,1)), axis=1)
        for label in range(3):
            assert np.all(np.sum(draws[:,indices]==label,axis=1)==np.sum(fixed==label))
    return draws


def replay(labels, groups, seed):
    random = np.random.default_rng(seed)
    out = np.zeros((PREFIX,3,3),dtype=int)
    for begin in range(0,PREFIX,250):
        n = min(250,PREFIX-begin)
        human = draw_species(labels['human'],groups['human'],random,n)
        mouse = draw_species(labels['mouse'],groups['mouse'],random,n)
        for i in range(3):
            for j in range(3):
                out[begin:begin+n,i,j] = np.sum((human==i)&(mouse==j),axis=1)
        assert np.all(out[begin:begin+n].sum(axis=(1,2))==len(labels['human']))
    return out


def main():
    protocol = json.loads((BASE/'protocol.json').read_text())
    manifest = json.loads((PROD/'analysis_manifest.json').read_text())
    for path, expected in manifest['input_sha256'].items():
        assert digest(ROOT/path)==expected, path
    for path, expected in manifest['output_sha256'].items():
        assert digest(PROD/path)==expected, path
    original_hashes = json.loads((HERE/'preflight_manifest.json').read_text())['input_sha256']
    for path, expected in original_hashes.items():
        assert digest(ROOT/path)==expected, path
    studies = read(HERE/'independent_study_covariates.tsv.gz')
    pairs = read(PROD/'fixed_background.tsv')
    assert len(pairs)==728 and pairs.human_gene_id.is_monotonic_increasing
    production_features = read(PROD/'ortholog_assay_features.tsv.gz')
    production_members = read(PROD/'measured_ortholog_members.tsv.gz')
    scan = read(PROD/'ortholog_3x3_threshold_scan.tsv')
    diagnostics = read(PROD/'ortholog_stratum_diagnostics.tsv.gz')
    mobility = read(PROD/'mobility_summary.tsv')
    checked = []
    with np.load(PROD/'null_draws.npz') as archives:
        for definition, sources in protocol['definitions'].items():
            source = {sp:studies[studies.species.eq(sp)&studies.source_contrast.eq(sources[sp])]
                      for sp in ['human','mouse']}
            for threshold in protocol['thresholds']:
                features, labels = {}, {}
                member = production_members[production_members.definition.eq(definition)&
                                             production_members.threshold.eq(threshold)]
                assert np.array_equal(member.human_gene_id,pairs.human_gene_id)
                for sp in ['human','mouse']:
                    features[sp],labels[sp],support = summarize(source[sp],pairs[sp+'_gene_id'],threshold)
                    actual = production_features[production_features.definition.eq(definition)&
                        production_features.threshold.eq(threshold)&production_features.species.eq(sp)]
                    assert len(actual)==728
                    assert np.array_equal(actual.gene_id,features[sp].gene_id)
                    for col in features[sp].columns:
                        if col in ['gene_id','measured_pattern']:
                            assert np.array_equal(actual[col],features[sp][col]), col
                        else:
                            close(actual[col],features[sp][col],(definition,threshold,sp,col))
                    assert np.array_equal(member[sp+'_direction'],np.array(LABELS)[labels[sp]])
                    for j,name in enumerate(['proximal_support','distal_support','measured_studies','opposing_study_allowed']):
                        assert np.array_equal(member[sp+'_'+name],np.array(support)[:,j])
                for model in protocol['models']:
                    groups = {sp:strata(features[sp],model) for sp in ['human','mouse']}
                    probabilities = {}
                    for sp in ['human','mouse']:
                        probabilities[sp],ind_diag,ind_mobile = probabilities_and_mobility(features[sp],labels[sp],groups[sp])
                        current = diagnostics[diagnostics.definition.eq(definition)&diagnostics.threshold.eq(threshold)&
                                              diagnostics.model.eq(model)&diagnostics.species.eq(sp)]
                        assert len(current)==len(ind_diag)
                        for col in ind_diag:
                            if col in ['gene_ids']:
                                assert np.array_equal(current[col],ind_diag[col]), col
                            else:
                                close(current[col],ind_diag[col],(definition,threshold,model,sp,col))
                        current_mobile = mobility[mobility.definition.eq(definition)&mobility.threshold.eq(threshold)&
                                                  mobility.model.eq(model)&mobility.species.eq(sp)]
                        assert len(current_mobile)==1
                        for col,value in ind_mobile.items():
                            close(current_mobile[col].iloc[0],value,(definition,threshold,model,sp,col))
                    expected = probabilities['human'].T @ probabilities['mouse']
                    observed = np.zeros((3,3),int)
                    for i in range(3):
                        for j in range(3):
                            observed[i,j] = np.sum((labels['human']==i)&(labels['mouse']==j))
                    current = scan[scan.contrast.eq('P_M')&scan.definition.eq(definition)&
                                   scan.threshold.eq(threshold)&scan.model.eq(model)]
                    assert len(current)==9
                    seed = 20261003 + int(hashlib.sha256('|'.join(map(str,
                        ['ortholog_direction','P_M',threshold,model])).encode()).hexdigest()[:12],16)
                    assert np.all(current.seed.eq(seed))
                    for i,hd in enumerate(LABELS):
                        for j,md in enumerate(LABELS):
                            row = current[current.human_direction.eq(hd)&current.mouse_direction.eq(md)].iloc[0]
                            assert row.observed==observed[i,j]
                            close(row.analytic_expected,expected[i,j],(definition,threshold,model,hd,md))
                    key = '|'.join(map(str,['ortholog','P_M',threshold,model,definition]))
                    independent = replay(labels,groups,seed)
                    assert archives[key].shape==(100000,3,3)
                    assert np.array_equal(independent,archives[key][:PREFIX]), key
                    checked.append(dict(definition=definition,threshold=threshold,model=model,
                                        observed=observed.tolist(),expected=expected.tolist(),
                                        replayed_prefix=PREFIX,seed=seed))
        with np.load(OLD/'null_draws.npz') as old:
            ps_keys = [key for key in old.files if key.split('|')[1]=='P_S']
            assert len(ps_keys)==18
            for key in ps_keys:
                assert np.array_equal(old[key],archives[key]), key
    report = dict(status='PASS',independently_rebuilt_actual_covariate_inputs=True,
        production_helpers_imported=False,contexts_checked=len(checked),
        feature_rows_checked=len(production_features),member_rows_checked=len(production_members),
        analytical_cells_checked=len(checked)*9,seeded_prefix_draws_per_context=PREFIX,
        seeded_prefix_draws_total=PREFIX*len(checked),all_prefix_stratum_margins_preserved=True,
        full_archive_tail_BH_verified_by='root/root_verify.py',full_PS_archive_contexts_unchanged=18,
        frozen_legacy_inputs_unchanged=len(original_hashes),
        input_sha256={str(path.relative_to(ROOT)):digest(path) for path in
            [BASE/'PROTOCOL.zh.md',BASE/'protocol.json',HERE/'independent_study_covariates.tsv.gz',
             PROD/'analysis_manifest.json',PROD/'null_draws.npz']},
        python=sys.version,numpy=np.__version__,pandas=pd.__version__)
    (HERE/'null_design_QA.json').write_text(json.dumps(report,indent=2)+'\n')
    (HERE/'independent_null_contexts.json').write_text(json.dumps(checked,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
