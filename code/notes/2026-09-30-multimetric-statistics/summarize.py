"""Summarize frozen APA results and render the supplementary figures."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import render_figures

HERE = Path(__file__).resolve().parent
NAMES = {'WPP': 'DWUI', 'RWUI': 'RWUI', 'MPRO': 'MPRO'}


def save(frame, name):
    frame.to_csv(HERE / name, sep='\t', index=False, na_rep='NA', float_format='%.17g')


def tables(routes):
    summaries, concordance, hidden = [], [], []
    for (species, study, route, contrast), f in routes.groupby(['species', 'study_id', 'route', 'contrast']):
        context = dict(species=species, study_id=study, route=route, contrast=contrast)
        for subset, g in [('all', f), ('PAS_ge3', f[f.PAS_number.ge(3)])]:
            for metric in ['WPP', 'RWUI', 'PDUI', 'PPUI']:
                d = g['delta_' + metric]
                summaries.append(context | dict(subset=subset, metric=NAMES.get(metric, metric), n=len(g),
                    mean=d.mean(), median=d.median(), negative_002=int((d <= -.02).sum()), positive_002=int((d >= .02).sum()),
                    negative_005=int((d <= -.05).sum()), positive_005=int((d >= .05).sum())))
            for metric, direction in [('RWUI', 1), ('PDUI', 1), ('PPUI', -1)]:
                a, b = g.delta_WPP, g['delta_' + metric] * direction
                strong = a.abs().ge(.02) & b.abs().ge(.02)
                concordance.append(context | dict(subset=subset, metric=metric, n=len(g),
                    rho=spearmanr(a, b).statistic, both_strong_n=int(strong.sum()),
                    same_direction_both_strong=int(((a * b > 0) & strong).sum()),
                    same_direction_fraction=float((a[strong] * b[strong] > 0).mean()) if strong.any() else np.nan))
            for threshold in [.10, .20]:
                yes = g.delta_WPP.abs().lt(.02) & g.MPRO.ge(threshold)
                hidden.append(context | dict(subset=subset, MPRO_threshold=threshold, n=len(g),
                    hidden_count=int(yes.sum()), hidden_fraction=float(yes.mean())))
    save(pd.DataFrame(summaries), 'global_metric_summary.tsv')
    save(pd.DataFrame(concordance), 'metric_concordance.tsv')
    save(pd.DataFrame(hidden), 'hidden_redistribution_summary.tsv')
    candidates = routes[routes.PAS_number.ge(3) & routes.PAS_number.le(6) &
                        routes.MPRO.ge(.10) & routes.delta_WPP.abs().lt(.02)].sort_values('MPRO', ascending=False)
    axes = pd.read_csv(HERE / 'frozen_PAS_axes.tsv.gz', sep='\t')
    suitable = []
    for row in candidates.itertuples():
        a, b = np.array(json.loads(row.target_usage)), np.array(json.loads(row.reference_usage))
        if (a >= .10).sum() < 2 or (b >= .10).sum() < 2 or row.APA_coverage < 50:
            continue
        ix = axes[axes.study_id.eq(row.study_id)].set_index('atlas_pas_id').loc[json.loads(row.PAS_ids)]
        changed = np.abs(a-b) >= .15
        if changed.sum() < 3 or np.ptp(ix.position.to_numpy()[changed]) < .35:
            continue
        suitable.append(row.Index)
    examples = candidates.loc[suitable].drop_duplicates('gene_id').head(2)
    assert len(examples) == 2
    save(examples, 'composition_examples.tsv')
    return pd.DataFrame(summaries), pd.DataFrame(concordance), pd.DataFrame(hidden), examples



def main():
    tables(pd.read_csv(HERE / 'route_effects.tsv.gz', sep='\t'))
    render_figures.main()


if __name__ == '__main__':
    main()
