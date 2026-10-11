"""Fixed-candidate route effects, study votes, LOSO and original-family profile sensitivity."""
from pathlib import Path
import os
for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[key] = '2'
os.environ['CUDA_VISIBLE_DEVICES'] = ''
import argparse
import fcntl
import hashlib
import json
import resource
import subprocess
import time
import traceback
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

BASE = Path('/root/nfsdata/apadetect/analysis_projects')
ROOT = BASE / 'reproapa_integration_validation_20260913'
RUN = ROOT / 'execution_v1'
PROFILES = RUN / 'A05/reaggregation_profiles_v1'
OUT = RUN / 'A05/reaggregated_effects_v1'
HUMAN = {'GSE118258': ['GSE118258'], 'GSE150311': ['GSE150311'], 'GSE178325': ['178325_0618'], 'GSE221739': ['GSE221739']}
BOTH = {'GSE115943': ['WOT'], 'GSE103221': ['CIP', 'OSK']}
CIP = {'GSE115943': ['WOT'], 'GSE103221': ['CIP']}
CONTRASTS = {'human': {'SP': HUMAN, 'SM': HUMAN, 'MP': HUMAN},
             'mouse': {'SR': BOTH, 'RE': BOTH, 'EP': CIP, 'SP': CIP, 'MP': CIP,
                       'SM': {**BOTH, 'GSE171127': ['171127_GETM', '171127_OSKM']}}}
STUDIES = {'human': sorted(HUMAN), 'mouse': ['GSE103221', 'GSE115943', 'GSE171127']}
AMPLITUDES = [.005, .01, .02, .03]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024**2), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return pd.read_csv(path, sep='\t', low_memory=False)


def dump(path, record):
    temp = path.with_suffix('.tmp'); temp.write_text(json.dumps(record, indent=2) + '\n'); temp.replace(path)


def save(frame, path):
    frame.to_csv(path, sep='\t', index=False, compression={'method': 'gzip', 'mtime': 0} if str(path).endswith('.gz') else None)


def status(phase, **extra):
    row = dict(phase=phase, pid=os.getpid(), updated_unix=time.time(), **extra)
    dump(OUT / 'status.json', row); print(json.dumps(row), flush=True)


def wait_profiles(paths, wait):
    while True:
        missing = [p for p in paths if not (p / 'complete.json').exists()]
        if not missing:
            break
        assert wait, ('Required profiles incomplete', [str(p) for p in missing])
        state = json.loads((PROFILES / 'status.json').read_text())
        live = subprocess.run(['ps', '-p', str(state['pid']), '-o', 'stat='], capture_output=True, text=True)
        assert live.returncode == 0 and live.stdout.strip() and not live.stdout.strip().startswith('Z'), state
        assert state['phase'] not in ['FAILED', 'REQUESTED_PROFILES_COMPLETE_REVIEW_PENDING'], state
        status('WAITING_FOR_PROFILES', first_missing=str(missing[0]), missing_studies=len(missing), producer_pid=state['pid'])
        time.sleep(30)
    sources = []
    for path in paths:
        record = json.loads((path / 'complete.json').read_text())
        for item in record['sources'] + record['outputs']:
            assert sha(item['path']) == item['sha256'], item['path']
        sources += [path / 'complete.json', path / 'gene_node_positions.tsv.gz']
    return sources


def study_votes(delta, amplitude):
    complete = np.isfinite(delta).all(axis=0)
    positive = complete & (delta >= amplitude).all(axis=0)
    negative = complete & (delta <= -amplitude).all(axis=0)
    conflict = complete & (delta >= amplitude).any(axis=0) & (delta <= -amplitude).any(axis=0)
    state = np.select([~complete, conflict, positive, negative],
                      ['not_evaluable', 'route_conflict', 'positive', 'negative'], default='weak')
    return positive.astype(int) - negative.astype(int), state


def majority(signs):
    need = signs.shape[0] // 2 + 1
    pos = (signs == 1).sum(axis=0); neg = (signs == -1).sum(axis=0)
    return np.where(pos >= need, 1, np.where(neg >= need, -1, 0)), pos, neg, need


def comparison(reference, current, complete_studies, need, assessable):
    return np.select([np.full(len(reference), not assessable), complete_studies < need,
         (reference != 0) & (current == reference), (reference != 0) & (current == -reference),
         (reference != 0) & (current == 0), (reference == 0) & (current != 0)],
        ['fewer_than_two_studies', 'insufficient_study_coverage', 'same_supported_direction',
         'opposite_supported_direction', 'no_longer_supported', 'new_support_on_fixed_candidate'],
        default='no_supported_direction')


def mean_available(values):
    count = np.isfinite(values).sum(axis=0)
    return np.divide(np.nansum(values, axis=0), count, out=np.full(values.shape[1], np.nan), where=count > 0)


def route_effects(data, members, species):
    genes = pd.Index(members.gene_id)
    frames = []
    for contrast, studies in CONTRASTS[species].items():
        a, b = contrast
        for study, routes in studies.items():
            for route in routes:
                part = data[data.route.eq(route)]
                assert part.study_id.eq(study).all()
                indexed = part.set_index(['gene_id', 'node'])
                assert indexed.index.is_unique
                left = indexed.reindex(pd.MultiIndex.from_product([genes, [a]]))
                right = indexed.reindex(pd.MultiIndex.from_product([genes, [b]]))
                frame = pd.DataFrame(dict(gene_id=genes, contrast=contrast, study_id=study, route=route,
                    from_node=a, to_node=b, from_position=left.position_index.to_numpy(),
                    to_position=right.position_index.to_numpy()))
                frame['delta'] = frame.to_position - frame.from_position
                for label, source in [('from', left), ('to', right)]:
                    frame[label + '_status'] = source.status.fillna('missing_node_row').to_numpy()
                    frame[label + '_libraries'] = source.n_libraries.to_numpy()
                    frame[label + '_units'] = source.qualified_metacells.to_numpy()
                frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def matrices(effects, genes, species, contrast):
    views = {}
    for study, routes in CONTRASTS[species][contrast].items():
        views[study] = np.stack([effects[effects.contrast.eq(contrast) & effects.route.eq(route)].set_index('gene_id').delta.reindex(genes).to_numpy() for route in routes])
    return views


def group_patterns(effects, original, members, species, variant):
    routes = list(dict.fromkeys(r for rr in CONTRASTS[species]['SP'].values() for r in rr))
    curves = []; patterns = []
    for level in ['family', 'subtype']:
        for group, entire in members.groupby(level, sort=True):
            family = entire.family.unique(); assert len(family) == 1
            expected = -1 if family[0] == 'C1' else 1
            for subset in ['all_original_members', 'original_stable_members']:
                ids = entire.gene_id if subset == 'all_original_members' else entire.loc[entire.stable_member, 'gene_id']
                route_means = {}; original_means = {}; route_coverage = {}
                for route in routes:
                    got = effects[effects.contrast.eq('SP') & effects.route.eq(route)].set_index('gene_id').delta.reindex(ids)
                    old = original[original.contrast.eq('SP') & original.route.eq(route)].set_index('gene_id').delta.reindex(ids)
                    paired = got.notna() & old.notna()
                    route_means[route] = float(got.mean()); original_means[route] = float(old.mean())
                    route_coverage[route] = int(got.notna().sum()) / len(ids) if len(ids) else np.nan
                    study = next(s for s, rr in CONTRASTS[species]['SP'].items() if route in rr)
                    curves.append(dict(variant=variant, species=species, group_level=level, group=group, family=family[0],
                         subset=subset, study_id=study, route=route, fixed_group_candidates=len(ids),
                         available_genes=int(got.notna().sum()), original_available_genes=int(old.notna().sum()),
                         paired_genes=int(paired.sum()), mean_SP=route_means[route], original_mean_SP=original_means[route],
                         paired_mean_SP_change=float((got[paired] - old[paired]).mean()),
                         available_fraction=route_coverage[route]))
                for omitted in ['NONE', *STUDIES[species]]:
                    required = {s: rr for s, rr in CONTRASTS[species]['SP'].items() if s != omitted}
                    study_mean = np.array([np.mean([route_means[r] for r in rr]) for rr in required.values()])
                    old_mean = np.array([np.mean([original_means[r] for r in rr]) for rr in required.values()])
                    fractions = [route_coverage[r] for rr in required.values() for r in rr]
                    for amp in AMPLITUDES:
                        state = ('fewer_than_two_studies' if len(required) < 2 else 'missing_study_pattern' if not np.isfinite(study_mean).all() else
                                 'same_family_pattern' if (expected * study_mean >= amp).all() else
                                 'opposite_family_pattern' if (-expected * study_mean >= amp).all() else 'mixed_or_weak_pattern')
                        patterns.append(dict(variant=variant, species=species, group_level=level, group=group, family=family[0],
                             subset=subset, excluded_study=omitted, amplitude=amp, required_studies=len(required),
                             evaluable_studies=int(np.isfinite(study_mean).sum()), fixed_group_candidates=len(ids),
                             minimum_route_candidate_fraction=float(np.min(fractions)) if fractions else np.nan,
                             mean_study_SP=float(np.mean(study_mean)), original_mean_study_SP=float(np.mean(old_mean)),
                             status=state, frozen_group_membership=True))
    return pd.DataFrame(curves), pd.DataFrame(patterns)


def calculate(species, variant, paths, members, original, baseline_votes, baseline_contrasts, fixed_inputs):
    folder = OUT / species / variant; folder.mkdir(parents=True, exist_ok=True)
    if (folder / 'complete.json').exists():
        done = json.loads((folder / 'complete.json').read_text())
        for item in done['sources'] + done['outputs']:
            assert sha(item['path']) == item['sha256'], item['path']
        return
    began = time.time(); status('RUNNING', species=species, variant=variant, stage='fixed_candidate_effects')
    data = pd.concat([read(p / 'gene_node_positions.tsv.gz') for p in paths], ignore_index=True)
    assert data.species.eq(species).all() and data.layer.eq('supported50').all()
    genes = pd.Index(members.gene_id); effects = route_effects(data, members, species)
    sources = fixed_inputs + [dict(path=str(p / name), sha256=sha(p / name)) for p in paths for name in ['complete.json', 'gene_node_positions.tsv.gz']]
    vote_rows = []; evidence = []; summary = []
    for contrast, studies in CONTRASTS[species].items():
        current = matrices(effects, genes, species, contrast)
        reference = matrices(original, genes, species, contrast)
        names = list(studies)
        current_mean = np.stack([current[s].mean(axis=0) for s in names])
        original_mean = np.stack([reference[s].mean(axis=0) for s in names])
        for amp in AMPLITUDES:
            vs = []; bs = []; complete = []
            for study in names:
                v, state = study_votes(current[study], amp); b, _ = study_votes(reference[study], amp)
                vs.append(v); bs.append(b); complete.append(np.isfinite(current[study]).all(axis=0))
                vote_rows.append(pd.DataFrame(dict(species=species, variant=variant, gene_id=genes, contrast=contrast,
                     study_id=study, amplitude=amp, study_vote=v, original_study_vote=b, status=state,
                     n_required_routes=len(studies[study]), study_mean_delta=current[study].mean(axis=0),
                     original_study_mean_delta=reference[study].mean(axis=0))))
            vs = np.stack(vs); bs = np.stack(bs); complete = np.stack(complete)
            full_baseline, _, _, _ = majority(bs)
            for omitted in ['NONE', *STUDIES[species]]:
                ix = [i for i, s in enumerate(names) if s != omitted]
                sign, pos, neg, need = majority(vs[ix]); fold_baseline, _, _, _ = majority(bs[ix])
                ncomplete = complete[ix].sum(axis=0); assessable = len(ix) >= 2
                common = np.isfinite(current_mean[ix]) & np.isfinite(original_mean[ix])
                count = common.sum(axis=0)
                cur_paired = mean_available(np.where(common, current_mean[ix], np.nan))
                old_paired = mean_available(np.where(common, original_mean[ix], np.nan))
                change = cur_paired - old_paired
                state = comparison(full_baseline, sign, ncomplete, need, assessable)
                within = comparison(fold_baseline, sign, ncomplete, need, assessable)
                frame = pd.DataFrame(dict(species=species, variant=variant, gene_id=genes, contrast=contrast,
                     amplitude=amp, excluded_study=omitted, excluded_study_in_contrast=omitted in names,
                     original_full_direction=full_baseline, original_fold_direction=fold_baseline,
                     recomputed_direction=sign.astype(float) if assessable else np.full(len(genes), np.nan),
                     positive_studies=pos, negative_studies=neg, n_required_studies=len(ix), required_support=need,
                     cross_study_assessable=assessable, complete_studies=ncomplete,
                     baseline_full_mean_delta=mean_available(original_mean),
                     current_mean_delta=mean_available(current_mean[ix]), original_fold_mean_delta=mean_available(original_mean[ix]),
                     paired_studies=count, paired_current_mean_delta=cur_paired, paired_original_mean_delta=old_paired,
                     paired_mean_delta_change=change, status=state, filter_status_within_fold=within,
                     fixed_candidate_denominator=len(genes)))
                evidence.append(frame)
                usable = (count >= 2) & np.isfinite(change)
                rank = float(spearmanr(cur_paired[usable], old_paired[usable]).statistic) if usable.sum() >= 3 and np.unique(cur_paired[usable]).size > 1 and np.unique(old_paired[usable]).size > 1 else np.nan
                for state_name, number in pd.Series(state).value_counts().sort_index().items():
                    summary.append(dict(species=species, variant=variant, contrast=contrast, amplitude=amp,
                         excluded_study=omitted, status=state_name, genes=int(number), fixed_candidate_denominator=len(genes),
                         fraction=int(number) / len(genes), paired_two_study_genes=int(usable.sum()),
                         paired_effect_spearman=rank, median_absolute_paired_effect_change=float(np.median(np.abs(change[usable]))) if usable.any() else np.nan))
    votes = pd.concat(vote_rows, ignore_index=True); claims = pd.concat(evidence, ignore_index=True)
    curves, patterns = group_patterns(effects, original, members, species, variant)
    replay = None
    if variant == 'original_all_states':
        keys = ['gene_id', 'contrast', 'study_id', 'route']
        a = effects.set_index(keys).sort_index(); b = original.set_index(keys).reindex(a.index)
        np.testing.assert_allclose(a[['from_position', 'to_position', 'delta']], b[['from_position', 'to_position', 'delta']], rtol=1e-12, atol=1e-12, equal_nan=True)
        keys = ['gene_id', 'contrast', 'study_id']
        a = votes[votes.amplitude.eq(.01)].set_index(keys).sort_index(); b = baseline_votes.set_index(keys).reindex(a.index)
        np.testing.assert_array_equal(a.study_vote, b.study_vote); assert a.status.equals(b.status)
        keys = ['gene_id', 'contrast']
        a = claims[claims.amplitude.eq(.01) & claims.excluded_study.eq('NONE')].set_index(keys).sort_index()
        b = baseline_contrasts.set_index(keys).reindex(a.index)
        np.testing.assert_array_equal(a.recomputed_direction, b.direction)
        old_group = read(ROOT / 'prior_evidence/APA' / species / 'family_assignment_evidence.tsv')
        a = curves[curves.group_level.eq('subtype') & curves.subset.eq('all_original_members')].set_index(['group', 'study_id', 'route']).sort_index()
        b = old_group.set_index(['subtype', 'study_id', 'route']).sort_index()
        assert list(a.index) == list(b.index)
        np.testing.assert_allclose(a.mean_SP, b.mean_SP, rtol=1e-12, atol=1e-12)
        replay = dict(status='ORIGINAL_ROUTE_VOTE_DIRECTION_AND_GROUP_PATTERN_REPLAY_PASS', route_rows=len(effects),
                      study_votes=len(baseline_votes), contrast_directions=len(baseline_contrasts), subtype_route_means=len(b))
    files = []
    for name, frame in [('route_effects', effects), ('study_votes', votes), ('candidate_claims', claims),
                         ('fixed_group_route_SP', curves), ('fixed_group_pattern_status', patterns), ('candidate_summary', pd.DataFrame(summary))]:
        path = folder / (name + '.tsv.gz'); save(frame, path); files.append(path)
    for item in sources:
        assert sha(item['path']) == item['sha256'], item['path']
    dump(folder / 'complete.json', dict(status='EFFECTS_COMPLETE_INDEPENDENT_REVIEW_PENDING', completed_unix=time.time(),
         seconds=time.time() - began, species=species, variant=variant, fixed_candidates=len(genes),
         route_effect_rows=len(effects), vote_rows=len(votes), candidate_claim_rows=len(claims),
         baseline_replay=replay, sources=sources, outputs=[dict(path=str(p), sha256=sha(p)) for p in files],
         reclustering_membership_evaluated=False, ortholog_pairs_evaluated=False, original_state_local_effects_evaluated=False))
    status('VARIANT_COMPLETE', species=species, variant=variant, seconds=time.time() - began)


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--species', nargs='+', choices=['human', 'mouse'], default=['human', 'mouse'])
    parser.add_argument('--variants', nargs='+'); parser.add_argument('--wait', action='store_true'); args = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_AS, (64 * 1024**3, 64 * 1024**3))
    OUT.mkdir(parents=True, exist_ok=True)
    lock = (OUT / 'run.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    original_inputs = json.loads((RUN / 'A05/input_receipt.json').read_text())
    for item in original_inputs:
        assert sha(item['path']) == item['sha256'], item['path']
    manifest = read(RUN / 'A05/candidate_cell_counts_v1/libraries.tsv')
    available = read(RUN / 'A05/reaggregation_units_v1/variant_coverage.tsv')
    allowed = available.loc[available.node_scope.eq('coarse'), 'variant'].tolist()
    variants = args.variants if args.variants else allowed
    assert set(variants).issubset(allowed) and len(variants) == len(set(variants))
    fixed = [Path(__file__), OUT / 'A05-effects-design.zh.md', RUN / 'A05/input_receipt.json',
             RUN / 'A05/reaggregation_units_v1/verification.json',
             RUN / 'A05/candidate_cell_counts_v1/libraries.tsv',
             RUN / 'A05/reaggregation_units_v1/variant_coverage.tsv']
    fixed += [ROOT / 'prior_evidence/APA' / species / 'family_assignment_evidence.tsv'
              for species in ['human', 'mouse']]
    inputs = [dict(path=str(p), sha256=sha(p)) for p in fixed] + original_inputs
    for species in args.species:
        members = read(ROOT / 'prior_evidence/APA' / species / 'aligned_gene_membership.tsv')
        assert members.gene_id.is_unique and len(members) == (1194 if species == 'human' else 1569)
        baseline = BASE / 'mouse_time_aware_integration_pilot_v1/results_apa_shared_local_v1' / (species + '_supported50')
        original = read(baseline / 'all_route_effects.tsv.gz'); original = original[original.gene_id.isin(members.gene_id)]
        votes = read(baseline / 'all_study_votes.tsv.gz'); votes = votes[votes.gene_id.isin(members.gene_id)]
        contrasts = read(baseline / 'all_gene_contrast_evidence.tsv.gz'); contrasts = contrasts[contrasts.gene_id.isin(members.gene_id)]
        assert set(contrasts.contrast) == set(CONTRASTS[species])
        studies = sorted(manifest.loc[manifest.species.eq(species), 'pas_dataset_id'].unique())
        for variant in variants:
            paths = [PROFILES / study / variant for study in studies]
            wait_profiles(paths, args.wait)
            if variant != 'original_all_states':
                done = json.loads((OUT / species / 'original_all_states/complete.json').read_text())
                assert done['baseline_replay']['status'] == 'ORIGINAL_ROUTE_VOTE_DIRECTION_AND_GROUP_PATTERN_REPLAY_PASS'
            calculate(species, variant, paths, members, original, votes, contrasts, inputs)
    status('REQUESTED_EFFECTS_COMPLETE_REVIEW_PENDING', species=args.species, variants=len(variants))


if __name__ == '__main__':
    try:
        main()
    except BaseException as exc:
        OUT.mkdir(parents=True, exist_ok=True); status('FAILED', error=str(exc), traceback=traceback.format_exc()); raise
