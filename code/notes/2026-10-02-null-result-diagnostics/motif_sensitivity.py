"""Post-hoc full matching grid; never choose controls using their motif counts."""
from pathlib import Path
import sys, json, hashlib, itertools
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'notes/2026-10-02-reviewer-followup'))
from apa_followup import read, save, rng_for, nullmod, INPUTS
from motif_background import fasta, runs, U, MODELS, ROUTES

GRID = [(g, p) for g in [.05, .10, .20, np.inf] for p in [.20, .40, np.inf]]


def name(g, p):
    return f'GC{g:g}_POS{p:g}'


def choose(candidates, g, p):
    good = candidates[(candidates[:, 2] <= g) & (candidates[:, 3] <= p)]
    if not len(good):
        return None
    return good[np.lexsort((good[:, 0], good[:, 2] / .05 + good[:, 3] / .2))[0]]


def build():
    regions = read(U / 'all_coordinate_regions.tsv.gz')
    pairs, attrition, inventory = [], [], []
    for species, routes in ROUTES.items():
        models = {}
        for r in routes:
            path = MODELS / species / r / 'models.json'
            INPUTS[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
            models[r] = json.loads(path.read_text())
        path = U / species / 'query.fa'
        INPUTS[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        seqs = fasta(path)
        hits = read(U / f'full-cohort-v2/data/{species}/site_effects.tsv.gz',
                    usecols=['gene_id', 'category', 'chrom', 'strand', 'start0', 'end0'])
        hits = hits[hits.category.isin(['miRNA', 'RBP'])].drop_duplicates()
        for category, h in hits.groupby('category'):
            inventory.append(dict(species=species, category=category, physical_sites=len(h), genes=h.gene_id.nunique()))
        hitgroups = dict(tuple(hits.groupby('gene_id')))
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
            if f['strand'] == '-':
                seq = seq[::-1]
            seq = seq[lo-model['start0']:hi-model['start0']].upper()
            assert len(seq) == hi-lo
            gc = np.r_[0, np.cumsum(np.array(list(seq)) == 'G') + np.cumsum(np.array(list(seq)) == 'C')]
            cuts = sorted({lo, hi} | {max(lo, min(hi, int(x))) for x in np.r_[group.start0, group.end0]})
            parts = []
            for a, b in zip(cuts[:-1], cuts[1:]):
                effects, changed, stable = [], True, True
                for route in routes:
                    block = group[group.route.eq(route) & group.start0.le(a) & group.end0.ge(b)]
                    if block.empty:
                        pos = np.array([p['cut_site']-1 for p in models[route][gene]['pas']])
                        constitutive = b <= pos.min() if f['strand'] == '+' else a >= pos.max()
                        changed = False
                        stable &= constitutive
                        effects.extend([0., 0.] if constitutive else [np.nan, np.nan])
                    else:
                        assert len(block) == 1
                        row = block.iloc[0]
                        effects.extend([row.effect20, row.effect50])
                        changed &= bool(row.observed_both5 and row.strict_change and row.eligible)
                        stable &= bool(row.eligible and np.isfinite(row.effect20+row.effect50) and max(abs(row.effect20), abs(row.effect50)) < .05)
                changed &= np.isfinite(effects).all() and (min(effects) >= .1 or max(effects) <= -.1)
                parts.append((a, b, 'changed' if changed else 'stable' if stable else 'other'))
            blocks = runs(parts, 'stable')
            hg = hitgroups.get(gene, hits.iloc[:0])
            for a, b, kind in parts:
                if kind != 'changed':
                    continue
                length = b-a
                fraction = (gc[b-lo]-gc[a-lo])/length
                starts = np.concatenate([np.arange(c, d-length+1) for c,d in blocks if d-c >= length]) if any(d-c >= length for c,d in blocks) else np.array([], dtype=int)
                fractions = (gc[starts+length-lo]-gc[starts-lo])/length
                distances = abs(starts+length/2-(a+b)/2)/(hi-lo)
                candidates = np.column_stack([starts, fractions, abs(fractions-fraction), distances])
                base = dict(species=species, contrast=contrast, gene_id=gene, start0=a, end0=b, length=length)
                selected = choose(candidates, .05, .20)
                status = ('no_stable_sequence' if not blocks else 'no_equal_length_window' if not len(starts)
                          else 'no_GC_match' if not (candidates[:,2] <= .05).any()
                          else 'no_joint_position_match' if selected is None else 'matched')
                attrition.append(dict(**base, status=status, stable_bp=sum(d-c for c,d in blocks),
                    longest_stable=max([d-c for c,d in blocks], default=0), equal_length_windows=len(starts),
                    gc_alone_windows=int((candidates[:,2] <= .05).sum()),
                    position_alone_windows=int((candidates[:,3] <= .2).sum()),
                    joint_windows=int(((candidates[:,2] <= .05)&(candidates[:,3] <= .2)).sum())))
                for g,p in GRID:
                    control = choose(candidates, g, p)
                    if control is None:
                        continue
                    c, control_gc, gc_diff, distance = control
                    c = int(c)
                    for category in ['miRNA', 'RBP']:
                        h = hg[hg.category.eq(category)]
                        nc = int((h.start0.ge(a)&h.end0.le(b)).sum())
                        ns = int((h.start0.ge(c)&h.end0.le(c+length)).sum())
                        pairs.append(dict(**base, setting=name(g,p), category=category,
                            baseline_window=selected is not None, control_start0=c, control_end0=c+length,
                            GC_change=fraction, GC_control=control_gc, GC_difference=gc_diff, relative_position_distance=distance,
                            changed_sites=nc, control_sites=ns, changed_density=1000*nc/length,
                            control_density=1000*ns/length, density_difference=1000*(nc-ns)/length))
    save(pd.DataFrame(attrition), HERE/'motif_attrition.tsv')
    save(pd.DataFrame(inventory), HERE/'motif_inventory.tsv')
    out = pd.DataFrame(pairs)
    save(out, HERE/'motif_windows.tsv.gz')
    return out


def summarize(frame):
    out, gout = [], []
    for cohort in ['dynamic', 'fixed_baseline_windows']:
        f = frame if cohort == 'dynamic' else frame[frame.baseline_window]
        genes = f.groupby(['setting','species','contrast','category','gene_id'])[['density_difference','changed_density','control_density']].mean().reset_index()
        gout.append(genes.assign(cohort=cohort))
        for (g,p), species, contrast, category in itertools.product(GRID, ROUTES, ['SP','MP'], ['miRNA','RBP']):
            b = genes[genes.setting.eq(name(g,p))&genes.species.eq(species)&genes.contrast.eq(contrast)&genes.category.eq(category)]
            x = b.density_difference.to_numpy()
            w = f[f.setting.eq(name(g,p))&f.species.eq(species)&f.contrast.eq(contrast)&f.category.eq(category)]
            nz = x[abs(x)>1e-12]
            mean = x.mean() if len(x) else np.nan
            row = dict(cohort=cohort, setting=name(g,p), species=species, contrast=contrast, category=category,
                n_genes=len(x), n_windows=len(w), n_nonzero=len(nz), n_positive=int((x>1e-12).sum()),
                n_negative=int((x < -1e-12).sum()), both_zero_genes=int(((b.changed_density==0)&(b.control_density==0)).sum()),
                mean_difference=mean, mean_changed_density=b.changed_density.mean(), mean_control_density=b.control_density.mean(),
                median_length=w.length.median(), median_GC_difference=w.GC_difference.median(), median_position_distance=w.relative_position_distance.median())
            if len(x) >= 10:
                # Same stream in every setting permits an exact replay of the primary baseline.
                rng = rng_for('motif', species, contrast, category)
                boot = x[rng.integers(len(x),size=(5000,len(x)))].mean(1)
                null = (rng.choice([-1,1],size=(10000,len(x)))*x).mean(1)
                p_mc = (1+int((abs(null)>=abs(mean)-1e-12).sum()))/10001
                if len(nz) <= 20:
                    sums = np.array([0.])
                    for v in nz:
                        sums = np.r_[sums+v,sums-v]
                    p_value = float((abs(sums)>=abs(x.sum())-1e-10).mean())
                    method = 'exact_nonzero_signs'
                else:
                    p_value, method = p_mc, '10000_sign_flips'
                row.update(lo=np.quantile(boot,.025), hi=np.quantile(boot,.975), p=p_value, p_mc=p_mc,
                           method=method, sign_p_floor=2**(1-len(nz)) if len(nz) else 1.)
            else:
                row.update(lo=np.nan,hi=np.nan,p=1.,p_mc=1.,method='fewer_than_10_genes_descriptive_only',sign_p_floor=np.nan)
            if len(x)>1:
                loo=(x.sum()-x)/(len(x)-1)
                row.update(loo_min=loo.min(),loo_max=loo.max())
            out.append(row)
    summary=pd.DataFrame(out)
    summary['q_96']=summary.groupby('cohort').p.transform(lambda v: nullmod.bh(v))
    summary['q_setting_8']=summary.groupby(['cohort','setting']).p.transform(lambda v: nullmod.bh(v))
    save(pd.concat(gout),HERE/'motif_genes.tsv.gz')
    save(summary,HERE/'motif_summary.tsv')
    old=read(ROOT/'notes/2026-10-02-reviewer-followup/motif-background/density_summary.tsv')
    base=summary[summary.cohort.eq('dynamic')&summary.setting.eq(name(.05,.2))].merge(old,on=['species','contrast','category'],suffixes=('_new','_old'))
    assert len(base)==8
    for col in ['n_genes','mean_difference','mean_changed_density','mean_control_density','lo','hi']:
        assert np.allclose(base[col+'_new'],base[col+'_old'],equal_nan=True),col
    assert np.allclose(base.p_mc,base.p_old)
    save(base,HERE/'motif_baseline_replay.tsv')
    (HERE/'motif_input_sha256.json').write_text(json.dumps(INPUTS,indent=2)+'\n')
    print(summary[summary.cohort.eq('dynamic')&summary.setting.isin([name(.05,.2),name(np.inf,np.inf)])].to_string(index=False),flush=True)


if __name__=='__main__':
    summarize(build())
