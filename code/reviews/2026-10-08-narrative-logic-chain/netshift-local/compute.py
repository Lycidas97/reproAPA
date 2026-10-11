"""Descriptive: gene-wide DWUI vs local strict recurrent interval changes (polyApipe contract-r2).
For each species x contrast, genes with >=1 strict recurrent interval pair (both routes).
Local direction = sign of interval inclusion change (dI<0: interval lost -> expected proximal, dDWUI<0).
Classify each gene by its gene-wide dDWUI in the same two routes and contrast."""
import csv, gzip, collections, sys, os
T = sys.argv[1]; OUT = sys.argv[2]
def rd(p):
    f = gzip.open(p,'rt') if p.endswith('.gz') else open(p)
    return csv.DictReader(f, delimiter='\t')
# gene-wide effects
eff = {}
# endpoint estimates used in the main text (supported50_p08); P_M paired, P_S documented source
cmap = {'P_M': 'MP', 'P_S': 'SP'}
for r in rd(f'{T}/endpoints/route_effects.tsv.gz'):
    if r['setting'] != 'supported50_p08' or r['contrast'] not in cmap: continue
    if r['delta_DWUI'] in ('', 'nan', 'NA'): continue
    eff[(r['species'], r['route'], r['gene_id'], cmap[r['contrast']])] = float(r['delta_DWUI'])
# admission path
adm = collections.defaultdict(set)
for r in rd(f'{T}/union_pool/figure5_union_claims.tsv'):
    adm[(r['species'], r['gene_id'], r['contrast'])] |= set(r['admission'].split(';'))
# strict recurrent interval pairs, dedup by region pair
pairs = {}
for r in rd(f'{T}/candidate_screen/recurrent_region_pairs.tsv.gz'):
    if r['recurrent_strict'] != 'True': continue
    k = (r['species'], r['contrast'], r['gene_id'], r['region_id_a'], r['region_id_b'])
    if k in pairs: continue
    ea = [float(r[x]) for x in ('effect20_a','effect50_a') if r[x] not in ('','nan')]
    eb = [float(r[x]) for x in ('effect20_b','effect50_b') if r[x] not in ('','nan')]
    pairs[k] = dict(route_a=r['route_a'], route_b=r['route_b'], sym=r['gene_symbol_a'],
                    dI=(sum(ea)/len(ea) + sum(eb)/len(eb))/2)
genes = collections.defaultdict(list)
for (sp, c, g, *_), v in pairs.items(): genes[(sp, c, g)].append(v)
rows = []
for (sp, c, g), vs in genes.items():
    signs = {1 if v['dI'] > 0 else -1 for v in vs}
    ra, rb = vs[0]['route_a'], vs[0]['route_b']
    da, db = eff.get((sp, ra, g, c)), eff.get((sp, rb, g, c))
    if len(signs) > 1: cls = 'mixed_local_direction'
    elif da is None or db is None: cls = 'gene_wide_unavailable'
    else:
        s = signs.pop()
        ok = [d*s > 0 for d in (da, db)]
        big = [abs(d) >= 0.02 for d in (da, db)]
        if all(ok) and all(big): cls = 'captured_same_sign_ge0.02_both'
        elif all(ok): cls = 'same_sign_but_lt0.02_in_ge1'
        elif any(not o and abs(d) >= 0.01 for o, d in zip(ok, (da, db))): cls = 'opposite_sign_ge0.01_in_ge1'
        else: cls = 'opposite_sign_lt0.01_only'
    a = adm.get((sp, g, c), set())
    path = 'MPRO_only' if a and all(x.startswith('MPRO') for x in a) else ('DWUI_path' if a else 'unknown')
    rows.append(dict(species=sp, contrast=c, gene_id=g, gene_symbol=vs[0]['sym'], n_strict_pairs=len(vs),
                     local_direction='mixed' if cls=='mixed_local_direction' else ('loss' if vs[0]['dI']<0 else 'gain'),
                     mean_dI=sum(v['dI'] for v in vs)/len(vs), route_a=ra, route_b=rb,
                     dDWUI_a=da, dDWUI_b=db, admission=path, cls=cls))
os.makedirs(OUT, exist_ok=True)
with open(f'{OUT}/gene_classes.tsv','w') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t'); w.writeheader(); w.writerows(rows)
order = ['captured_same_sign_ge0.02_both','same_sign_but_lt0.02_in_ge1','opposite_sign_lt0.01_only',
         'opposite_sign_ge0.01_in_ge1','mixed_local_direction','gene_wide_unavailable']
summ = collections.Counter((r['species'], r['contrast'], r['admission'], r['cls']) for r in rows)
with open(f'{OUT}/summary.tsv','w') as f:
    f.write('species\tcontrast\tadmission\t' + '\t'.join(order) + '\ttotal\n')
    for sp in ('mouse','human'):
        for c in ('MP','SP'):
            for ad in ('all','DWUI_path','MPRO_only'):
                cnt = [sum(v for (s,cc,a,k),v in summ.items() if s==sp and cc==c and k==o and (ad=='all' or a==ad)) for o in order]
                if sum(cnt): f.write(f'{sp}\t{c}\t{ad}\t' + '\t'.join(map(str,cnt)) + f'\t{sum(cnt)}\n')
print(open(f'{OUT}/summary.tsv').read())
for r in rows:
    if r['gene_symbol'] in ('Riok3','Trip12','Fxr1','Eif4ebp2'): print(r)
