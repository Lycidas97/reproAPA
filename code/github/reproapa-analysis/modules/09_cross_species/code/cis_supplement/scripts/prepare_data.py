from pathlib import Path
import shutil,hashlib,json
import pandas as pd
p=Path(__file__).resolve().parents[1]
u=p/'sources/upstream'
manifest=json.loads((p/'sources/manifest.json').read_text())
for item in manifest:
    snapshot=p/item['snapshot']
    assert hashlib.sha256(snapshot.read_bytes()).hexdigest()==item['sha256'],snapshot
mi=pd.read_csv(u/'data/figure4-cis-architecture-v27/annot/model_input_with_seq_annot.tsv.gz',sep='\t')
primary=mi[(mi.layer=='supported50')&(mi.tau==.02)].copy()
keys=['layer','contrast','discovery','tau','human_gene_id']
assert not mi.duplicated(keys).any()
for c in [x[:-8] for x in mi if x.endswith('_missing')]:
 primary.loc[primary[c+'_missing']==1,c]=float('nan')
primary.to_csv(p/'data/primary_observed_features.tsv.gz',sep='\t',index=False)
match=pd.read_csv(u/'data/figure4-recurrence-validation-v26/matched_comparators.tsv.gz',sep='\t')
match=match[(match.layer=='supported50')&(match.tau==.02)].copy()
features=['A_g','driver_conserved','frac_conserved_discovery_pas','phastcons_driver_mean','motif_rank_driver_mean']
match=match.merge(primary[keys+features],on=keys,validate='many_to_one')
match['class']=match.apply(lambda r:'Matched background' if r.role=='control' else ('Shared' if r.validation_same_direction else 'Recurrent-only'),axis=1)
match.to_csv(p/'data/descriptive_matched_features.tsv.gz',sep='\t',index=False)
counts=match.groupby(['contrast','discovery','class']).agg(unique_genes=('human_gene_id','nunique'),records=('human_gene_id','size'),weight_sum=('match_weight','sum')).reset_index()
counts.to_csv(p/'analysis/descriptive_class_counts.tsv',sep='\t',index=False)
print(counts.to_string(index=False))
