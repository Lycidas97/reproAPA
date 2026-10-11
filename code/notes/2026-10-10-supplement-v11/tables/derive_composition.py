"""Route-mean composition by sign (Figure 2c values) as a table: the same computation as p_s_effects() in
notes/2026-10-10-divergence-main-figures/compose.py, with its asserts against the frozen global_direction_summary."""
import hashlib, json
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[3]
B = 'notes/2026-10-06-supplement-polyapipe/contract-r2/bundle/tables/global/'
reads = {}
def read(rel):
    p = ROOT/rel; reads[rel] = hashlib.sha256(p.read_bytes()).hexdigest(); return pd.read_csv(p, sep='\t')
q = read(B+'qualified_date_effects.tsv.gz'); q = q[(q.setting == 'supported50_p08') & (q.contrast == 'P_S')]
g = q.groupby(['route', 'gene_id']).APA_delta.mean().reset_index()
b = read(B+'background_gene_membership.tsv.gz')
b = b[(b.setting == 'supported50_p08') & (b.contrast == 'P_S') & (b.background == 'common_one_to_one_ortholog_all_settings')]
s = read(B+'global_direction_summary.tsv')
s = s[(s.setting == 'supported50_p08') & (s.contrast == 'P_S') & (s.background == 'common_one_to_one_ortholog_all_settings')]
NAME = {'GSE118258': ('human', 'Xing2020'), '178325_0618': ('human', 'Guan2022–0618'), '178325_0809': ('human', 'Guan2022–0809'),
        'CIP': ('mouse', 'Guo2019–CIP'), 'WOT': ('mouse', 'Schiebinger2019–2i')}
PCT = [10, 25, 50, 75, 90]; rows = []; R = {}
for r, (sp, nm) in NAME.items():
    x = g[(g.route == r) & g.gene_id.isin(b[b.route == r].gene_id)].APA_delta.values
    row = s[s.route == r].iloc[0]
    assert len(x) == row.n and abs(x.mean()-row['mean']) < 1e-9 and (x <= -0.02).sum() == row.proximal_002 and (x >= 0.02).sum() == row.distal_002
    pr, di = x[x <= -0.02], x[x >= 0.02]
    v = dict(n_genes=len(x), mean_dDWUI=x.mean(), distal_arm=x[x > 0].sum()/len(x), proximal_arm=x[x < 0].sum()/len(x),
             n_proximal_responders=len(pr), n_distal_responders=len(di), fraction_proximal_responders=len(pr)/len(x), fraction_distal_responders=len(di)/len(x),
             mean_dDWUI_proximal_responders=pr.mean(), mean_dDWUI_distal_responders=di.mean(),
             **{f'proximal_responders_p{p}': v_ for p, v_ in zip(PCT, np.percentile(pr, PCT))}, **{f'distal_responders_p{p}': v_ for p, v_ in zip(PCT, np.percentile(di, PCT))})
    assert abs(v['distal_arm']-row.positive_contribution) < 1e-9 and abs(v['proximal_arm']-row.negative_contribution) < 1e-9
    R[r] = v; rows.append(dict(level='route', species=sp, study=nm, route=r, **v))
avg = lambda keys, f: {k: f(k) for k in keys}
K = [k for k in R['CIP'] if k not in ('n_genes', 'n_proximal_responders', 'n_distal_responders')]
hum = avg(K, lambda k: (R['GSE118258'][k]+(R['178325_0618'][k]+R['178325_0809'][k])/2)/2)
mou = avg(K, lambda k: (R['CIP'][k]+R['WOT'][k])/2)
rows.append(dict(level='species_study_equal', species='human', study='Xing2020, Guan2022 (series averaged)', route='', **hum))
rows.append(dict(level='species_study_equal', species='mouse', study='Guo2019–CIP, Schiebinger2019–2i', route='', **mou))
dpos, dneg = hum['distal_arm']-mou['distal_arm'], hum['proximal_arm']-mou['proximal_arm']; tot = dpos+dneg
rows.append(dict(level='human_minus_mouse', species='', study='', route='', mean_dDWUI=tot, distal_arm=dpos, proximal_arm=dneg,
                 distal_arm_share=dpos/tot, proximal_arm_share=dneg/tot))
out = pd.DataFrame(rows); out.to_csv(Path(__file__).parent/'derived/route_mean_composition.tsv', sep='\t', index=False)
json.dump(dict(reads=reads, rows=len(out), note='mirrors compose.py p_s_effects(); asserts against global_direction_summary passed'),
          open(Path(__file__).parent/'derived/route_mean_composition.receipt.json', 'w'), indent=1)
print(out[['level', 'species', 'route', 'mean_dDWUI', 'distal_arm', 'proximal_arm', 'fraction_proximal_responders']].to_string())
print('shares', dpos/tot, dneg/tot)
