"""Conservative annotation-compatible tandem-UTR subset; no isoform inference."""
from collections import defaultdict
import gzip
import json
import re
import numpy as np
import pandas as pd
from analyze import HERE, PROJECTS, M, H, PROFILES, read, save, sha, INPUTS, metrics, orthologs


def terminal_utrs(path, genes, species):
    INPUTS[str(path)] = sha(path)
    transcripts = {}
    with gzip.open(path, 'rt') as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            a = line.rstrip().split('\t')
            if a[2] not in ['exon', 'stop_codon']:
                continue
            match = re.search(r'gene_id "([^"]+)"', a[8])
            if not match or match[1] not in genes:
                continue
            attr = dict(re.findall(r'(\w+) "([^"]*)"', a[8]))
            tid = attr.get('transcript_id')
            if not tid:
                continue
            start, end = int(a[3]), int(a[4])
            if tid not in transcripts:
                transcripts[tid] = dict(species=species, gene_id=match[1], transcript_id=tid,
                    chrom=a[0].removeprefix('chr'), strand=a[6], exons=[], stops=[])
            rec = transcripts[tid]
            if a[2] == 'exon':
                # Only the terminal exon is needed; do not retain whole gene models.
                if not rec['exons'] or (a[6] == '+' and end > rec['exons'][0][1]) or (a[6] == '-' and start < rec['exons'][0][0]):
                    rec['exons'] = [(start, end)]
            else:
                rec['stops'].append((start, end))
    rows = []
    for rec in transcripts.values():
        if not rec['exons'] or not rec['stops']:
            continue
        left, right = rec['exons'][0]
        if not all(left <= a <= b <= right for a, b in rec['stops']):
            continue
        low = max(b for a, b in rec['stops']) + 1 if rec['strand'] == '+' else left
        high = right if rec['strand'] == '+' else min(a for a, b in rec['stops']) - 1
        if low <= high:
            rows.append({k: v for k, v in rec.items() if k not in ['exons', 'stops']} | dict(utr_start=low, utr_end=high))
    return pd.DataFrame(rows)


def main():
    effects = read(HERE / 'qualified_route_effects.tsv.gz')
    root = PROJECTS / '2026-09-13_human-mouse-ortholog-mapping/sources/web/ensembl116/genes'
    utrs = []
    for species, name in [('human', 'homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz'),
                          ('mouse', 'mus_musculus/Mus_musculus.GRCm39.116.gtf.gz')]:
        genes = set(effects.loc[effects.species.eq(species), 'gene_id'])
        utrs.append(terminal_utrs(root / name, genes, species))
        print('Terminal UTRs parsed', species, len(utrs[-1]), flush=True)
    utrs = pd.concat(utrs, ignore_index=True)
    save(utrs, 'annotation_terminal_utrs.tsv.gz')
    by_gene = dict(tuple(utrs.groupby('gene_id')))
    memberships = []
    folders = {}
    for path in sorted(PROFILES.glob('*/selected_metacells.tsv.gz')):
        m = read(path)
        study, species = m.study_id.iloc[0], m.species.iloc[0]
        if study not in set(effects.study_id):
            continue
        folder = path.parent.name
        folders[study] = folder
        source = (M / 'data/metacell-pas-v7-v1/studies' if species == 'mouse' else H / 'data/results_metacell_pas_v2/studies') / folder
        f = read(source / 'features.tsv.gz', usecols=['atlas_pas_id', 'cut_site', 'chrom', 'strand', 'primary_gene_id', 'keep_for_apa_primary'])
        f = f[f.keep_for_apa_primary & f.primary_gene_id.isin(effects.loc[effects.study_id.eq(study), 'gene_id'])].copy()
        f['chrom'] = f.chrom.astype(str).str.removeprefix('chr')
        f['signed_position'] = f.cut_site * np.where(f.strand.eq('+'), 1, -1)
        span = f.groupby('primary_gene_id').signed_position.agg(lambda x: x.max() - x.min())
        support = (M / 'data/pas-annotation-support-v1/studies' if species == 'mouse' else M / 'data/apa-majority-shared-v1/human_annotation') / folder / 'primary_pas_reference_matches.tsv.gz'
        f = f.merge(read(support, usecols=['atlas_pas_id', 'matched_source_count_50nt']), on='atlas_pas_id', validate='one_to_one')
        for layer, part in [('supported50', f[f.matched_source_count_50nt.ge(1)]), ('expanded_primary', f)]:
            for gene, q in part.groupby('primary_gene_id'):
                u = by_gene.get(gene)
                fits = pd.DataFrame()
                if u is not None and len(q) >= 2 and q.chrom.nunique() == 1 and q.strand.nunique() == 1:
                    fits = u[u.chrom.eq(q.chrom.iloc[0]) & u.strand.eq(q.strand.iloc[0]) &
                             u.utr_start.le(q.cut_site.min()) & u.utr_end.ge(q.cut_site.max())]
                memberships.append(dict(species=species, study_id=study, gene_id=gene, layer=layer,
                    retained_pas=len(q), native_axis_span_nt=span.loc[gene], tandem_utr_compatible=not fits.empty,
                    compatible_transcripts=';'.join(fits.transcript_id) if not fits.empty else '',
                    native_retained_min=q.cut_site.min(), native_retained_max=q.cut_site.max()))
        print('PAS annotation checked', study, flush=True)
    membership = pd.DataFrame(memberships)
    save(membership, 'tandem_utr_eligibility.tsv.gz')
    effects['layer'] = np.where(effects.setting.eq('expanded_p08'), 'expanded_primary', 'supported50')
    x = effects.merge(membership, on=['species', 'study_id', 'gene_id', 'layer'], validate='many_to_one')
    assert len(x) == len(effects)
    x = x[x.tandem_utr_compatible].copy()
    x['delta_nt'] = x.delta * x.native_axis_span_nt
    assert ((np.sign(x.delta_nt) == np.sign(x.delta)) | x.delta.eq(0)).all()
    save(x, 'tandem_utr_effects.tsv.gz')
    keys = ['species', 'study_id', 'route', 'contrast', 'design']
    rows = []
    for key, g in x.groupby(keys):
        fixed = set.intersection(*(set(q.gene_id) for _, q in g.groupby('setting')))
        for setting, q in g.groupby('setting'):
            for background, ids in [('available', set(q.gene_id)), ('fixed_settings_within_route', fixed)]:
                z = q[q.gene_id.isin(ids)]
                if z.empty:
                    continue
                rows.append(dict(zip(keys, key), setting=setting, background=background,
                    **metrics(z.delta), mean_delta_nt=z.delta_nt.mean(), median_delta_nt=z.delta_nt.median()))
    save(pd.DataFrame(rows), 'tandem_utr_summary.tsv')
    common_rows, sets = [], {}
    eligible = x[x.design.isin(['true_source', 'same_library'])]
    for (sp, con), q in eligible.groupby(['species', 'contrast']):
        sets[sp, con] = set.intersection(*(set(z.gene_id) for _, z in q.groupby(['route', 'setting'])))
    orth = orthologs()
    for con in ['P_S', 'P_M']:
        pairs = orth[orth.human_gene_id.isin(sets.get(('human', con), set())) & orth.mouse_gene_id.isin(sets.get(('mouse', con), set()))]
        for sp in ['human', 'mouse']:
            for background, ids in [('common_within_species_all_settings', sets.get((sp, con), set())),
                                    ('common_one_to_one_ortholog_all_settings', set(pairs[sp + '_gene_id']))]:
                for key, q in eligible[eligible.species.eq(sp) & eligible.contrast.eq(con)].groupby(keys + ['setting']):
                    z = q[q.gene_id.isin(ids)]
                    if z.empty:
                        continue
                    common_rows.append(dict(zip(keys + ['setting'], key), background=background, **metrics(z.delta),
                                            mean_delta_nt=z.delta_nt.mean(), median_delta_nt=z.delta_nt.median()))
    save(pd.DataFrame(common_rows), 'tandem_utr_common_background_summary.tsv')
    (HERE / 'tandem_utr_input_sha256.json').write_text(json.dumps(INPUTS, indent=2) + '\n')
    print(pd.DataFrame(rows).query("setting == 'supported50_p08' and background == 'available'")[keys + ['n', 'mean', 'median', 'mean_delta_nt', 'median_delta_nt']].to_string(index=False))


if __name__ == '__main__':
    main()
