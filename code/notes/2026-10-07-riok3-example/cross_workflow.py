"""Riok3 in both PAS workflows: one summary table plus the per-library inclusion rows of each (2026-10-07)."""
import json, hashlib
from pathlib import Path
import pandas as pd
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
P50 = ROOT/'notes/2026-10-06-dual-main-plot-data/apadetect-p50/main/tables'
D = HERE/'tables'; G = 'ENSMUSG00000024404'; IV = (12288586, 12289475)
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
inputs = {str(p.relative_to(ROOT)): sha(p) for p in [P50/'geometry/axis-decomposition.tsv', P50/'examples/inclusion_paired_libraries.tsv.gz', P50/'panel_statistics/inclusion-paired-means.tsv', P50/'candidate_screen/recurrent_region_pairs.tsv.gz']}
ax_p = pd.read_csv(D/'axis-decomposition.tsv', sep='\t'); ax_a = pd.read_csv(P50/'geometry/axis-decomposition.tsv', sep='\t'); ax_a = ax_a[ax_a.gene_id.eq(G)]
ip_p = pd.read_csv(D/'inclusion_paired_libraries.tsv', sep='\t'); ip_p = ip_p[ip_p.start0.eq(IV[0]) & ip_p.end0.eq(IV[1])]
ip_a = pd.read_csv(P50/'examples/inclusion_paired_libraries.tsv.gz', sep='\t'); ip_a = ip_a[ip_a.gene_id.eq(G) & ip_a.kind.eq('interval')]
# p50 regions containing the shared interval (CIP region equals it; WOT region extends downstream; the complete-retention mask is the same)
rr = pd.read_csv(P50/'candidate_screen/recurrent_region_pairs.tsv.gz', sep='\t'); rr = rr[rr.gene_id.eq(G) & rr.contrast.eq('MP')]
src = {}
for r in rr.itertuples():
    for side in 'ab':
        src[getattr(r, 'route_'+side)] = (int(getattr(r, 'region_start0_'+side)), int(getattr(r, 'region_end0_'+side)))
ip_a = pd.concat([ip_a[ip_a.route.eq(rt) & ip_a.start0.eq(s) & ip_a.end0.eq(e)] for rt, (s, e) in src.items()])
rows = []
for wf, ax, ip in [('polyApipe', ax_p, ip_p), ('APADETECT p50', ax_a, ip_a)]:
    for a in ax.sort_values('route').itertuples():
        row = dict(workflow=wf, route=a.route, study_id=a.study_id, native_pas=int(a.native_pas), native_pairs=int(a.n_native_pairs),
                   native_delta_DWUI=a.native_delta_DWUI, native_MPRO=a.native_MPRO, terminal_native_pas=int(a.terminal_native_pas),
                   terminal_shift_nt=a.terminal_shift_nt, within_family_DWUI=a.within_family_DWUI, family_mass_DWUI=a.family_mass_DWUI,
                   nonterminal_DWUI=a.nonterminal_DWUI)
        for c in (20, 50):
            q = ip[ip.route.eq(a.route) & ip.cutoff.eq(c)]
            row[f'interval_pairs_{c}'] = len(q); row[f'delta_I_pp_{c}'] = q.delta.mean()*100
            row[f'libraries_{c}'] = ';'.join(sorted(q.library_key))
        row['interval_source_region'] = f'{int(ip.start0.iloc[0]) if wf=="polyApipe" else src[a.route][0]}-{IV[1] if wf=="polyApipe" else src[a.route][1]}'
        rows.append(row)
S = pd.DataFrame(rows)
for c in (20, 50):
    for rt in ('CIP', 'WOT'):
        a, b = S[(S.workflow == 'polyApipe') & (S.route == rt)][f'libraries_{c}'].iloc[0], S[(S.workflow != 'polyApipe') & (S.route == rt)][f'libraries_{c}'].iloc[0]
        print(rt, c, 'same paired libraries in both workflows:', a == b)
S['shared_interval'] = 'chr18:[12288586,12289475) +, 889 nt'
S.to_csv(HERE/'tables/riok3-cross-workflow-summary.tsv', sep='\t', index=False)
ip_a.to_csv(HERE/'tables/inclusion_paired_libraries.p50.tsv', sep='\t', index=False)
json.dump(inputs, open(HERE/'cross-workflow-inputs.json', 'w'), indent=1)
print(S.drop(columns=[c for c in S if c.startswith('libraries_')]).T)
