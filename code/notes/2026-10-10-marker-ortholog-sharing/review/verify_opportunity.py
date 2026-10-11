"""Independent reconstruction of the frozen joint bridge identity background."""
from pathlib import Path
import hashlib,json

import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
TASK=HERE.parent
MARKER=ROOT/'notes/2026-10-10-marker-reference-sensitivity'
OLD=ROOT/'notes/2026-10-10-ortholog-threshold-scan'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(x,name):x.to_csv(HERE/name,sep='\t',index=False,na_rep='NA')

def main():
 protocol=json.loads((TASK/'protocol.json').read_text())
 assert protocol['background_n']==728
 cols=['study_id','route','gene_id','contrast','background','n_libraries',
       'contributor_libraries','weights','PAS_ids','version','setting']
 q=pd.read_csv(MARKER/'unique_support_bridge_route_effects.tsv.gz',sep='\t',usecols=cols)
 q=q[q.version.eq(protocol['marker_version']) & q.setting.eq(protocol['setting']) & q.background.eq(protocol['bridge_background'])]
 key=['study_id','route','gene_id']
 a=q[q.contrast.eq('D00')].set_index(key).sort_index()
 b=q[q.contrast.eq('D01')].set_index(key).sort_index()
 assert a.index.equals(b.index) and a.index.is_unique
 assert a[['n_libraries','contributor_libraries','weights','PAS_ids']].equals(b[['n_libraries','contributor_libraries','weights','PAS_ids']])
 opp=a.reset_index()
 allowed={('GSE103221','CIP'),('GSE115943','WOT'),('GSE118258','GSE118258'),
          ('GSE150311','GSE150311'),('GSE221739','GSE221739')}
 assert set(zip(opp.study_id,opp.route))==allowed
 species={x:'mouse' if x in ['GSE103221','GSE115943'] else 'human' for x,_ in allowed}
 opp['species']=opp.study_id.map(species)
 measured={s:set(g.groupby('gene_id').study_id.nunique().loc[lambda x:x>=2].index) for s,g in opp.groupby('species')}
 old=pd.read_csv(OLD/'fixed_ortholog_background.tsv',sep='\t')
 old=old[old.contrast.eq('P_M')].copy().reset_index(drop=True)
 assert len(old)==751 and old.human_gene_id.is_unique and old.mouse_gene_id.is_unique
 old['human_bridge_ge2']=old.human_gene_id.isin(measured['human'])
 old['mouse_bridge_ge2']=old.mouse_gene_id.isin(measured['mouse'])
 old['joint_eligible']=old.human_bridge_ge2 & old.mouse_bridge_ge2
 assert int(old.human_bridge_ge2.sum())==735 and int(old.joint_eligible.sum())==728
 joint=old[old.joint_eligible].copy().reset_index(drop=True)
 save(old,'independent_eligibility751.tsv');save(joint,'independent_background728.tsv')
 previous=pd.read_csv(OLD/'measured_ortholog_members.tsv.gz',sep='\t')
 previous=previous[previous.contrast.eq('P_M') & previous.threshold.eq(.02) & previous.human_direction.eq('proximal') & previous.mouse_direction.eq('proximal')]
 assert len(previous)==27
 audit=previous[['human_gene_id','mouse_gene_id']].merge(old,on=['human_gene_id','mouse_gene_id'],how='left',validate='one_to_one')
 assert len(audit)==27 and audit.joint_eligible.sum()==25
 for sp in ['human','mouse']:
  counts=opp[opp.species.eq(sp)].groupby('gene_id').study_id.nunique()
  audit[sp+'_bridge_studies']=audit[sp+'_gene_id'].map(counts).fillna(0).astype(int)
 audit['qualification_status']=audit.joint_eligible.map({True:'measurable',False:'not_estimable'})
 save(audit,'independent_original27_qualification.tsv')
 hashes={str(p.relative_to(ROOT)):sha(p) for p in [TASK/'PROTOCOL.zh.md',TASK/'protocol.json',
         MARKER/'unique_support_bridge_route_effects.tsv.gz',OLD/'fixed_ortholog_background.tsv',OLD/'measured_ortholog_members.tsv.gz']}
 out={'status':'PASS','background751':751,'human_only_eligible':735,'joint_background':728,
      'original27_measurable':25,'original27_not_estimable':2,
      'D00_D01_contributor_PAS_weight_identity_pass':True,'new_direction_labels_computed':False,
      'input_sha256':hashes,'not_estimable_rows':audit[~audit.joint_eligible][['human_gene_id','mouse_gene_id','human_bridge_studies','mouse_bridge_studies']].to_dict('records')}
 (HERE/'opportunity_QA.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(out,indent=2))

if __name__=='__main__':main()
