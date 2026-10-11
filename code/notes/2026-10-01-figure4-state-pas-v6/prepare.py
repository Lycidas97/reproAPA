"""Freeze Figure 4 inputs and verify each PAS illustration against its native DWUI."""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
DATA=HERE/'data'
DATA.mkdir(exist_ok=True)
manifest={}
def read(source,name):
    p=ROOT/source
    manifest[str(source)]=hashlib.sha256(p.read_bytes()).hexdigest()
    d=pd.read_csv(p,sep='\t')
    d.to_csv(DATA/name,sep='\t',index=False)
    return d

pair=read('table/Supplementary-Source-Data-20260928/opposite-background-sharing/sharing_null_summary.tsv','pair_tests.tsv')
read('table/Supplementary-Source-Data-20260928/opposite-background-sharing/pair_background_global_directions.tsv','pair_backgrounds.tsv')
read('notes/2026-09-29-core-sensitivity/figure4-matched-overlap/overlap_null_results.tsv','set_tests.tsv')
values=read('Fig/main-figure-candidates/20260927-global-background-link/Figure.4.source-data.tsv','member_values.tsv')
flags=read('Fig/main-figure-candidates/20260927-global-background-link/Figure.4.member-outcomes.tsv','member_flags.tsv')
assert len(flags)==18
p=ROOT/'notes/2026-09-30-multimetric-statistics/route_effects.tsv.gz'
manifest[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
d=pd.read_csv(p,sep='\t')
d=d[d.gene_symbol.str.upper().isin(['SZRD1','RABIF'])&d.route.isin(['GSE118258','CIP','WOT'])].copy()
d['display_gene']=d.gene_symbol.str.upper()
assert len(d)==12
pas=[];checks=[];common_p=[]
for (gene,route),g in d.groupby(['display_gene','route']):
    ps=g[g.contrast.eq('P_S')].iloc[0];pm=g[g.contrast.eq('P_M')].iloc[0]
    assert ps.design=='true_source'
    assert ps.PAS_ids==pm.PAS_ids
    same=bool(np.allclose(json.loads(ps.target_usage),json.loads(pm.target_usage),atol=1e-12))
    common_p.append(dict(gene=gene,route=route,same_P_vector=same))
    if route in ['GSE118258','CIP']:
        assert same
for r in d.itertuples():
    ids=json.loads(r.PAS_ids)
    coords=np.array([int(i.split(':')[3])for i in ids])
    strand=ids[0].split(':')[4]
    dist=coords-coords.min() if strand=='+' else coords.max()-coords
    target=np.array(json.loads(r.target_usage));reference=np.array(json.loads(r.reference_usage))
    np.testing.assert_allclose([target.sum(),reference.sum()],[1,1],atol=1e-12)
    effect=(target-reference)@dist/dist.max()
    np.testing.assert_allclose(effect,r.delta_WPP,atol=1e-12)
    checks.append(dict(gene=r.display_gene,route=r.route,contrast=r.contrast,reported=r.delta_WPP,reconstructed=effect,span_nt=int(dist.max()),n_pas=len(ids)))
    for rank,i in enumerate(np.argsort(dist),1):
        pas.append(dict(gene=r.display_gene,route=r.route,contrast=r.contrast,pas_id=ids[i],pas_rank=rank,chrom=ids[i].split(':')[2],coordinate=int(coords[i]),strand=strand,distance_nt=int(dist[i]),P=target[i],reference=reference[i],reference_state='S' if r.contrast=='P_S' else 'M'))
d.to_csv(DATA/'case_effects.tsv',sep='\t',index=False)
pd.DataFrame(pas).to_csv(DATA/'case_PAS_usage.tsv',sep='\t',index=False)
pd.DataFrame(checks).to_csv(DATA/'case_representation_checks.tsv',sep='\t',index=False)
pd.DataFrame(common_p).to_csv(DATA/'case_common_P_checks.tsv',sep='\t',index=False)
# Confirm SZRD1 examples agree with the matrix carried by the original main figure.
for route in ['GSE118258','CIP','WOT']:
    old=values[values.kind.eq('gene')&values.human_symbol.eq('SZRD1')&values.route.eq(route)].delta.iloc[0]
    new=d[d.display_gene.eq('SZRD1')&d.route.eq(route)&d.contrast.eq('P_S')].delta_WPP.iloc[0]
    np.testing.assert_allclose(old,new,atol=1e-9)
protected={}
for p in [ROOT/'Fig/Figure.4.pdf',ROOT/'main/results/cross-species-apa-recurrence-is-concentrated-in-selective-ortholog-state-relationships.tex']:
    protected[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
(HERE/'source-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
(HERE/'protected-files.json').write_text(json.dumps(protected,indent=2)+'\n')
(HERE/'data-verification.json').write_text(json.dumps({'status':'PASS','endpoint_rows':len(d),'all_PAS_vectors_sum_to_one':True,'all_native_DWUI_reconstructed':True,'same_P_vector_across_references_in_displayed_case_routes':True,'displayed_case_routes':['GSE118258','CIP'],'WOT_not_used_for_three_state_profiles':'P contributor summaries can differ across P-S and P-M','SZRD1_matches_original_main_matrix':True,'inference_rerun':False},indent=2)+'\n')
print('Frozen sources and 12 native endpoint/PAS checks: PASS')
