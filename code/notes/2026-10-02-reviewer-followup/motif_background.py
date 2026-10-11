"""Same-gene matched stable UTR windows, preserving the frozen hit inventory."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from apa_followup import HERE, ROOT, read, save, rng_for, nullmod

U = ROOT / 'notes/2026-10-01-figure5-union-update'
MODELS = ROOT / 'notes/2026-10-01-mpro-biology/annotation'
ROUTES = {'mouse': ['CIP', 'WOT'], 'human': ['GSE118258', 'GSE150311']}


def fasta(path):
    out, key, pieces = {}, None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            if key is not None:
                out[key] = ''.join(pieces)
            key, pieces = line[1:].split()[0], []
        else:
            pieces.append(line.strip())
    if key is not None:
        out[key] = ''.join(pieces)
    return out


def runs(parts, label):
    blocks = []
    for a, b, kind in parts:
        if kind == label:
            if blocks and blocks[-1][1] == a:
                blocks[-1] = (blocks[-1][0], b)
            else:
                blocks.append((a, b))
    return blocks


def main():
    regions = read(U / 'all_coordinate_regions.tsv.gz')
    pairs, failures, coverage = [], [], []
    for species, routes in ROUTES.items():
        models = {r: json.loads((MODELS / species / r / 'models.json').read_text()) for r in routes}
        seqs = fasta(U / species / 'query.fa')
        hits = read(U / f'full-cohort-v2/data/{species}/site_effects.tsv.gz',
                    usecols=['gene_id', 'category', 'chrom', 'strand', 'start0', 'end0'])
        hits = hits[hits.category.isin(['miRNA', 'RBP'])].drop_duplicates()
        hitgroups = {g: b for g, b in hits.groupby('gene_id')}
        regs = regions[regions.species.eq(species)]
        for (contrast, gene), group in regs.groupby(['contrast', 'gene_id']):
            if set(group.route) != set(routes):
                continue
            families = [models[r][gene]['family'] for r in routes]
            if len({(f['chrom'], f['strand']) for f in families}) != 1:
                continue
            lo, hi = max(f['common_start0'] for f in families), min(f['common_end0'] for f in families)
            if hi <= lo:
                continue
            f = families[0]
            available = [m for m in f['models'] if m['region_id'] in seqs and m['start0'] <= lo and m['end0'] >= hi]
            if not available:
                continue
            model = available[0]
            seq = seqs[model['region_id']]
            # Convert transcription-oriented sequence back to genomic ascending order for GC only.
            if f['strand'] == '-':
                seq = seq[::-1]
            seq = seq[lo - model['start0']:hi - model['start0']].upper()
            assert len(seq) == hi - lo
            gc = np.r_[0, np.cumsum(np.array(list(seq)) == 'G') + np.cumsum(np.array(list(seq)) == 'C')]
            cuts = sorted({lo, hi} | {max(lo, min(hi, int(x))) for x in np.r_[group.start0, group.end0]})
            parts = []
            for a, b in zip(cuts[:-1], cuts[1:]):
                effects, changed, stable = [], True, True
                for route in routes:
                    block = group[group.route.eq(route) & group.start0.le(a) & group.end0.ge(b)]
                    if block.empty:
                        # A region upstream of every retained PAS is structurally constitutive.
                        model_r = models[route][gene]
                        pos = np.array([p['cut_site'] - 1 for p in model_r['pas']])
                        constitutive = b <= pos.min() if f['strand'] == '+' else a >= pos.max()
                        changed = False
                        stable &= constitutive
                        effects.extend([0., 0.] if constitutive else [np.nan, np.nan])
                    else:
                        assert len(block) == 1
                        row = block.iloc[0]
                        effects.extend([row.effect20, row.effect50])
                        changed &= bool(row.observed_both5 and row.strict_change and row.eligible)
                        stable &= bool(row.eligible and np.isfinite(row.effect20 + row.effect50) and max(abs(row.effect20), abs(row.effect50)) < .05)
                changed &= np.isfinite(effects).all() and (min(effects) >= .1 or max(effects) <= -.1)
                parts.append((a, b, 'changed' if changed else 'stable' if stable else 'other'))
            stable_blocks = runs(parts, 'stable')
            changes = [(a, b) for a, b, kind in parts if kind == 'changed']
            coverage.append(dict(species=species, contrast=contrast, gene_id=gene, changed_intervals=len(changes), stable_bp=sum(b-a for a,b in stable_blocks)))
            hg = hitgroups.get(gene, hits.iloc[:0])
            for a, b in changes:
                length = b - a
                center = (a + b) / 2
                fraction = (gc[b-lo] - gc[a-lo]) / length
                best = None
                for c, d in stable_blocks:
                    if d - c < length:
                        continue
                    starts = np.arange(c, d - length + 1, dtype=int)
                    frac = (gc[starts + length-lo] - gc[starts-lo]) / length
                    distance = abs(starts + length / 2 - center) / (hi - lo)
                    good = (abs(frac-fraction) <= .05) & (distance <= .2)
                    if not good.any():
                        continue
                    starts, frac, distance = starts[good], frac[good], distance[good]
                    score = abs(frac - fraction) / .05 + distance / .2
                    j = int(np.argmin(score))
                    candidate = (float(score[j]), int(starts[j]), float(frac[j]), float(distance[j]))
                    if best is None or candidate < best:
                        best = candidate
                if best is None:
                    failures.append(dict(species=species, contrast=contrast, gene_id=gene, start0=a, end0=b, status='no_length_GC_position_matched_stable_window'))
                    continue
                _, c, control_gc, distance = best
                for category in ['miRNA', 'RBP']:
                    h = hg[hg.category.eq(category)]
                    n_change = int((h.start0.ge(a) & h.end0.le(b)).sum())
                    n_control = int((h.start0.ge(c) & h.end0.le(c + length)).sum())
                    pairs.append(dict(species=species, contrast=contrast, gene_id=gene, category=category,
                        start0=a, end0=b, control_start0=c, control_end0=c+length, length=length,
                        GC_change=fraction, GC_control=control_gc, relative_position_distance=distance,
                        changed_sites=n_change, control_sites=n_control,
                        changed_density=1000*n_change/length, control_density=1000*n_control/length,
                        density_difference=1000*(n_change-n_control)/length))
    out = HERE / 'motif-background'
    save(pd.DataFrame(pairs), out / 'matched_windows.tsv.gz')
    save(pd.DataFrame(failures), out / 'unmatched_windows.tsv')
    save(pd.DataFrame(coverage), out / 'coverage.tsv')
    frame = pd.DataFrame(pairs)
    genes = frame.groupby(['species','contrast','category','gene_id'])[['density_difference','changed_density','control_density']].mean().reset_index()
    save(genes, out / 'gene_density.tsv')
    summaries = []
    for species in ROUTES:
        for contrast in ['SP','MP']:
            for category in ['miRNA','RBP']:
                g = genes[genes.species.eq(species) & genes.contrast.eq(contrast) & genes.category.eq(category)]
                x = g.density_difference.to_numpy()
                row = dict(species=species, contrast=contrast, category=category, n_genes=len(x),
                           mean_difference=float(x.mean()) if len(x) else np.nan,
                           mean_changed_density=float(g.changed_density.mean()), mean_control_density=float(g.control_density.mean()))
                if len(x) >= 10:
                    rng = rng_for('motif',species,contrast,category)
                    boot = x[rng.integers(len(x),size=(5000,len(x)))].mean(1)
                    null = (rng.choice([-1,1],size=(10000,len(x))) * x).mean(1)
                    row.update(lo=float(np.quantile(boot,.025)),hi=float(np.quantile(boot,.975)),
                               p=(1+int((abs(null)>=abs(x.mean())-1e-12).sum()))/10001,status='matched_gene_permutation')
                else:
                    row.update(lo=np.nan,hi=np.nan,p=1.,status='fewer_than_10_genes_descriptive_only')
                summaries.append(row)
    summary = pd.DataFrame(summaries)
    summary['q_8_tests'] = nullmod.bh(summary.p)
    save(summary,out / 'density_summary.tsv')
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
