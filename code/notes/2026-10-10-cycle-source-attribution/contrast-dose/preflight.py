"""Coverage/design only: no cycle scores or APA effect columns read."""
from pathlib import Path
import json
import pandas as pd
H=Path(__file__).resolve().parent
R=H.parents[2]
B=R/'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables'
cols=['species','study_id','route','contrast','design','setting','gene_id','library_key_a','library_key_b','time_a','time_b','sample_unit_a','sample_unit_b','metacell_ids_a','metacell_ids_b','n_mc_a','n_mc_b','weight']
p=pd.read_csv(B/'state_vectors/pair_effects.primary.tsv.gz',sep='\t',usecols=cols)
p=p[p.setting.eq('supported50_p08')].copy()
s=pd.read_csv(R/'notes/2026-10-10-proliferation-apa/metacell_cycle_scores.tsv.gz',sep='\t',usecols=['study_id','metacell_id','route','library_key','node','purity','n_cells','experimental_time'])
assert not s.duplicated(['study_id','metacell_id']).any()
lookup={(r.study_id,r.metacell_id):r for r in s.itertuples()}
audit=[]
for side,state in [('a','P'),('b',None)]:
 for row in p.itertuples():
  ids=json.loads(getattr(row,'metacell_ids_'+side))
  for mid in ids:
   k=(row.study_id,mid)
   if k not in lookup:
    audit.append(dict(study_id=row.study_id,metacell_id=mid,problem='score_missing'))
    continue
   rec=lookup[k]
   expected=state or row.contrast[-1]
   if rec.route!=row.route or rec.library_key!=getattr(row,'library_key_'+side) or rec.node!=expected or rec.purity<.8:
    audit.append(dict(study_id=row.study_id,metacell_id=mid,problem='identity_state_library_mismatch'))
assert not audit,audit[:3]
rows=[]
for key,part in p.groupby(['species','study_id','route','contrast','design'],sort=True):
 u=part.drop_duplicates(['library_key_a','library_key_b'])
 rows.append(dict(zip(['species','study_id','route','contrast','design'],key))|dict(genes=part.gene_id.nunique(),pair_rows=len(part),target_libraries=u.library_key_a.nunique(),reference_libraries=u.library_key_b.nunique(),physical_library_pairs=len(u),same_library_pairs=int((u.library_key_a==u.library_key_b).sum()),same_time_pairs=int((u.time_a==u.time_b).sum()),all_members_cycle_covered=True))
pd.DataFrame(rows).to_csv(H/'endpoint_coverage.tsv',sep='\t',index=False)
p[['species','study_id','route','contrast','design','library_key_a','library_key_b','time_a','time_b','sample_unit_a','sample_unit_b']].drop_duplicates().to_csv(H/'endpoint_library_pairs.tsv',sep='\t',index=False)
t=pd.read_csv(B/'matched_pms/weighted_triplets.tsv.gz',sep='\t',usecols=['species','study_id','route','setting','gene_id','library_key_a','library_key_b','triplet_weight','metacell_ids_a','metacell_ids_b'])
t=t[t.setting.eq('supported50_p08')]
q=p[p.contrast.eq('P_M')].set_index(['study_id','route','gene_id','library_key_a'])
assert q.index.is_unique
joined=t.merge(q[['metacell_ids_b','library_key_b']].rename(columns={'metacell_ids_b':'M_metacell_ids','library_key_b':'M_library_key'}),left_on=['study_id','route','gene_id','library_key_a'],right_index=True,how='left',validate='many_to_one')
assert joined.M_metacell_ids.notna().all()
assert joined.library_key_a.eq(joined.M_library_key).all()
matched=[]
for key,part in joined.groupby(['species','study_id','route']):
 matched.append(dict(zip(['species','study_id','route'],key))|dict(genes=part.gene_id.nunique(),triplet_rows=len(part),target_libraries=part.library_key_a.nunique(),source_libraries=part.library_key_b.nunique(),all_PM_same_library=True,all_members_cycle_covered=True))
pd.DataFrame(matched).to_csv(H/'matched_triplet_coverage.tsv',sep='\t',index=False)
(H/'preflight_validation.json').write_text(json.dumps({'status':'PASS','primary_pair_rows':len(p),'cycle_metacells':len(s),'all_pair_members_score_covered':True,'matched_triplets':len(t),'matched_M_members_recovered_from_frozen_PM_pairs':True,'new_cycle_or_APA_effects_read_or_computed':False},indent=2)+'\n')
print(pd.DataFrame(rows).to_string(index=False))
print(pd.DataFrame(matched).to_string(index=False))
