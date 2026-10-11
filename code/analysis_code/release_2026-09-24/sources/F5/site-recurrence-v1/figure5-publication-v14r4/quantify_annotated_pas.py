"""Sum raw counts by annotated PAS, then recompute site inclusion from annotation cuts."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent.parent


def main():
    mapping = pd.read_csv(ROOT / 'panel_EF_pas_annotation_mapping.tsv', sep='\t')
    sites = pd.read_csv(ROOT / 'panel_EF_sites.tsv', sep='\t')
    assert mapping.status.eq('matched').all(), 'Unmatched PAS require an explicit coverage policy.'
    summaries, libraries, counts_out, masks_out, effects, changes = [], [], [], [], [], []
    source_hashes = {str(ROOT/name): hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                     for name in ['panel_EF_pas_annotation_mapping.tsv', 'panel_EF_sites.tsv']}
    checked, max_oracle_error = 0, 0.
    for route in ['CIP', 'WOT']:
        path = PROJECT / f'full-cohort-v2/data/mouse/{route}/models.json'
        source_hashes[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        models = json.loads(path.read_text())
        for (gid, gene), m in mapping[mapping.route.eq(route)].groupby(['gene_id', 'gene']):
            model = models[gid]
            raw_pas = model['pas']
            m = m.set_index('atlas_pas_id').loc[[p['atlas_pas_id'] for p in raw_pas]]
            assert np.array_equal(m.raw_cut, [p['cut_site'] for p in raw_pas])
            ann = m[['annotation_id', 'annotation_cut']].drop_duplicates().sort_values('annotation_cut', ascending=False)
            assert not ann.annotation_id.duplicated().any()
            ann_ids = ann.annotation_id.tolist()
            membership = np.array([[int(v == a) for a in ann_ids] for v in m.annotation_id])
            assert (membership.sum(axis=1) == 1).all()
            for site in sites[sites.gene.eq(gene)].itertuples():
                assert site.strand == '-'
                mask = int(site.start0) >= ann.annotation_cut.to_numpy(int)-1
                for aid, cut, included in zip(ann_ids, ann.annotation_cut, mask):
                    masks_out.append(dict(gene=gene, route=route, site=site.site, annotation_id=aid,
                                          annotation_cut=int(cut), included=bool(included)))
                for threshold in [20, 50]:
                    for state in 'SREPM':
                        values, raw_values = [], []
                        for lib in model['profiles'].get(state, {}).get('libraries', []):
                            raw_counts = np.array(lib['counts'], dtype=np.int64)
                            assert np.isclose(raw_counts.sum(), lib['total'])
                            counts = raw_counts @ membership
                            assert counts.sum() == raw_counts.sum()
                            # No PAS are dropped here, so coverage and eligibility remain unchanged.
                            if counts.sum() < threshold or lib['coverage'] < .5:
                                continue
                            assert np.allclose(raw_counts/raw_counts.sum(), lib['usage'])
                            value = float(counts @ mask / counts.sum())
                            direct = sum(int(c) for c, cut in zip(raw_counts, m.annotation_cut)
                                         if int(site.start0) >= int(cut)-1) / raw_counts.sum()
                            assert np.isclose(value, direct, atol=1e-12)
                            max_oracle_error = max(max_oracle_error, abs(value-direct))
                            raw_mask_lo = np.array([int(site.start0) >= p['member_cut_site_max']-1 for p in raw_pas])
                            raw_mask_hi = np.array([int(site.start0) >= p['member_cut_site_min']-1 for p in raw_pas])
                            assert np.array_equal(raw_mask_lo, raw_mask_hi), 'Original site has ambiguous PAS bounds.'
                            raw_value = float(raw_counts @ raw_mask_lo / raw_counts.sum())
                            values.append(value)
                            raw_values.append(raw_value)
                            libraries.append(dict(gene=gene, route=route, site=site.site, threshold=threshold,
                                state=state, library=lib['library'], total_family_UMI=int(counts.sum()),
                                coverage=lib['coverage'], inclusion=value, raw_inclusion=raw_value))
                            checked += 1
                            if threshold == 20 and site.site == sites[sites.gene.eq(gene)].site.min():
                                for aid, cut, count in zip(ann_ids, ann.annotation_cut, counts):
                                    counts_out.append(dict(gene=gene, route=route, state=state,
                                        library=lib['library'], annotation_id=aid, annotation_cut=int(cut),
                                        count=int(count), usage=float(count/counts.sum())))
                        assert values, (gene, route, site.site, threshold, state)
                        summaries.append(dict(gene=gene, route=route, gene_id=gid, site=site.site,
                            threshold=threshold, state=state, n_libraries=len(values), inclusion=np.mean(values),
                            minimum=min(values), maximum=max(values), raw_inclusion=np.mean(raw_values)))
    profiles = pd.DataFrame(summaries)
    for (gene, route, site, threshold), g in profiles.groupby(['gene', 'route', 'site', 'threshold']):
        g = g.set_index('state')
        for contrast, ref in [('SP', 'S'), ('MP', 'M'), ('EP', 'E')]:
            value = float((g.loc['P', 'inclusion']-g.loc[ref, 'inclusion'])*100)
            raw_value = float((g.loc['P', 'raw_inclusion']-g.loc[ref, 'raw_inclusion'])*100)
            effects.append(dict(gene=gene, route=route, site=site, threshold=threshold,
                contrast=contrast, effect_pp=value, raw_effect_pp=raw_value, difference_pp=value-raw_value))
    old = pd.read_csv(ROOT / 'panel_EF_profiles.tsv', sep='\t')
    for row in old.itertuples():
        site = 1 if row.gene == 'Hmgn2' else 2
        x = profiles[profiles.gene.eq(row.gene)&profiles.route.eq(row.route)&profiles.site.eq(site)&
                     profiles.threshold.eq(20)&profiles.state.eq(row.state)].iloc[0]
        assert np.isclose(x.raw_inclusion, row.inclusion), (row.gene, row.route, row.state)
        changes.append(dict(gene=row.gene, route=row.route, state=row.state,
            old_inclusion=row.inclusion, new_inclusion=x.inclusion, difference_pp=(x.inclusion-row.inclusion)*100))
    for filename, data in [('annotation_site_profiles.tsv', summaries), ('annotation_site_effects.tsv', effects),
                           ('annotation_library_profiles.tsv', libraries), ('annotation_pas_counts.tsv', counts_out),
                           ('annotation_site_pas_masks.tsv', masks_out), ('annotation_curve_changes.tsv', changes)]:
        pd.DataFrame(data).to_csv(ROOT / filename, sep='\t', index=False)
    effect_frame = pd.DataFrame(effects)
    changed = [r for r in changes if abs(r['difference_pp']) > 1e-10]
    assert {(r['gene'], r['route'], r['state']) for r in changed} == {
        ('Hmgn2', 'CIP', 'R'), ('Hmgn2', 'CIP', 'E')}, 'Frozen annotation-boundary regression changed.'
    assert all(r['difference_pp'] < 0 for r in changed)
    qa = dict(status='PASS', quantification='annotation-PAS counts summed before usage and site-inclusion calculation',
        counts_conserved=True, no_unmatched_pas=True, coverage_and_library_eligibility_unchanged=True,
        checked_library_site_threshold_rows=checked, max_direct_oracle_error=max_oracle_error,
        original_profiles_reconstructed=True, thresholds=[20,50], input_sha256=source_hashes,
        changed_state_rows=changed,
        effects20=effect_frame[effect_frame.threshold.eq(20)&effect_frame.site.isin([1,2])].to_dict('records'))
    (ROOT / 'annotation_quantification_QA.json').write_text(json.dumps(qa, indent=2))
    print(json.dumps({k:v for k,v in qa.items() if k not in ['input_sha256','effects20']}, indent=2))


if __name__ == '__main__':
    main()
