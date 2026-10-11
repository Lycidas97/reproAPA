"""Audit descriptive coarse states and aggregate existing OT within each time pair."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

P = Path(__file__).resolve().parents[1]
D = P/'data/state-flow-v1'
LABELS = {'S':'成纤维来源 / 来源样','R':'来源退出与重塑','E':'上皮化 / 混合过渡',
          'P':'多能性程序获得','N':'naive 相关候选','M':'间质 / ECM 保留',
          'K':'角质细胞样候选','X':'内胚层相关分泌程序','U':'混合程序 / 未明确集合'}
# Explicit review decisions; original source annotations are retained below.
GROUPS = [
 ('GSE118258','S','001 002 003','FN1;GREM1;THBS1;BGN','成纤维来源程序；周期亚型保留'),
 ('GSE118258','R','011','CD248;LOX;PRRX2;CAV1','来源相关早期重塑；过程定位候选，不是新身份'),
 ('GSE118258','E','004 006','KRT8;SPINT2;CD24;MIR302CHG','上皮相关及早期目标程序；原工作注释 unresolved，本次仅作宽过程候选'),
 ('GSE118258','E','005','H19;MFAP5;COL3A1;SOX4','来源与目标相关混合群；候选纳入过渡窗口，不称纯 MET'),
 ('GSE118258','P','010','TDGF1;ESRG;GDF3;NANOG;LIN28A','多能性支持成套增加；naive/primed 混合信号，不确证亚型'),
 ('GSE118258','M','007 009 012','FN1;COL6A3;MFAP5;LUM','来源/ECM 持续或晚期保留，不等同功能性失败'),
 ('GSE118258','K','008','KRT6A;KRT17;LY6D;SLURP2','角质程序高；原上皮化候选细化为角质样，保留候选'),
 ('GSE150311','S','001 002 003 004 007 008 009','COL1A1;COL6A2;DCN;PRRX1','成纤维样集合，含 IFN 高亚型；混合 capture 无时间箭头'),
 ('GSE150311','K','006','SLURP1;SLURP2;LY6D;KLF4','角质程序候选；原 process_dominant 不作为具体身份'),
 ('GSE150311','E','010 013','KRT8;KRT18;MIR302CHG;COL1A1','上皮或来源/目标混合窗口；无确定时间方向'),
 ('GSE150311','M','011 012 015','COL1A2;COL6A2;DCN;S100A4','保留间质成分；C015 原 unresolved，仅作宽家族候选'),
 ('GSE150311','N','005 017 018 020 021 023','DNMT3L;DPPA3;NLRP7;ALPG','naive 相关成套程序；跨培养条件亚型不等价，非鼠 naive 对应检验'),
 ('GSE150311','P','016 024 025 026 027','LIN28A;TDGF1;DPPA4;ESRG','多能性程序家族；C016 原 unresolved、C024 原 naive 候选仍保留亚型不确定'),
 ('GSE150311','U','014 019 022','NEFL;TUBB3;HAND1;SOX11','神经或 HAND1 上皮混合成分；不确证神经/滋养层命运'),
 ('GSE178325','S','001 002 003 013','SERPINE1;CAV1;PRRX1;FN1','来源/Stage I day 0 成纤维样；旧 process_dominant 不足以概括身份'),
 ('GSE178325','R','004 005 006 007 008','MMP3;KRT19;RARRES1;RBP1;APCDD1','早中期间质重塑与上皮成分变化；共同过程窗口而非单一可塑性身份'),
 ('GSE178325','X','009','TTR;APOA1;AFP;FABP1;GATA6','内胚层相关分泌程序；与旧 naive 候选冲突，拒绝并入 naive 或成熟多能性；不确证 XEN'),
 ('GSE178325','M','010 011','POSTN;IGFBP5;COL3A1;MGP','晚期间质/ECM 与上皮混合分支；不等同单一失败终点'),
 ('GSE178325','P','012 014 015 016 017 018','DPPA4;L1TD1;ESRG;DNMT3B;TDGF1','培养终点及成熟参考的目标程序；没有现有时间 OT 入边，参考与路线终点分列'),
 ('GSE221739','S','003','DCN;LUM;COL3A1;TWIST2','观测 D3 起的来源样；本项目缺 D0，不能作为已观测的未诱导来源'),
 ('GSE221739','M','001','COL5A1;COL6A3;FN1;ACTA2;THBS1','ECM/收缩性间质程序；细胞周期和稳定性标记保留'),
 ('GSE221739','P','002 004','DPPA4;DNMT3B;ESRG;L1TD1;TDGF1','多能性相关程序强，C002/C004 内部亚型不拆成必经阶段；不靠外源 POU5F1/SOX2/NANOG/LIN28 判定'),
 ('GSE221739','U','005','RPL23A;RPL7;RPL11;GAPDH','低 RNA/核糖体主导混合群；不以稀疏目标基因或阶段强指定身份'),
]


def main():
    primary_path=P/'data/results_v1/prepared/primary_cells.tsv.gz'
    flow_path=P/'data/results_temporal_v1/layout/state_condition_flows.tsv.gz'
    inputs=[primary_path,flow_path,P/'data/results_v1/prepared/temporal_edge_plan.tsv',D/'state_gene_expression.tsv']
    hashes={str(f.relative_to(P)):hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs}
    cells=pd.read_csv(primary_path,sep='\t',low_memory=False)
    assert cells.cell_id.is_unique and len(cells)==215114
    rows=[]
    for study,coarse,clusters,markers,reason in GROUPS:
        for cluster in clusters.split():
            rows.append(dict(study_id=study,cluster_id='C'+cluster,coarse_state=coarse,
                coarse_label=LABELS[coarse],review_marker_set=markers,review_qualification=reason,
                annotation_level='descriptive_candidate',cross_study_equivalence_validated=False))
    mapping=pd.DataFrame(rows)
    counts=cells.groupby(['study_id','cluster_id']).size().rename('n_cells').reset_index()
    mapping=mapping.merge(counts,on=['study_id','cluster_id'],validate='one_to_one')
    assert len(mapping)==62 and mapping.n_cells.sum()==215114
    mapping['state_id']=mapping.study_id+'::studywide::'+mapping.cluster_id
    sources=[]; binding=[]
    for study,g in cells.groupby('study_id'):
        a=D/'source_annotations'/study
        if (a/'cell_annotation.tsv.gz').exists():
            original=pd.read_csv(a/'cell_annotation.tsv.gz',sep='\t')
            original=original[original.qc_tier.eq('primary_pass')]
            linked=g.merge(original[['cell_id','local_cluster_id']],on='cell_id',validate='one_to_one')
            assert len(g)==len(original)==len(linked) and linked.cluster_id.eq(linked.local_cluster_id).all()
            work=pd.read_csv(a/'cluster_working_annotation.tsv',sep='\t').rename(columns={'local_state_id':'cluster_id'})
            sources.append(work[['study_id','cluster_id','working_semantic_annotation','confirmatory_annotation',
                                 'contributing_programs','manual_hold','prohibited_claim']])
            binding.append(dict(study_id=study,cells=len(g),clusters=g.cluster_id.nunique(),exact_binding=True,
                                source='existing_working_annotation'))
        else:
            binding.append(dict(study_id=study,cells=len(g),clusters=g.cluster_id.nunique(),exact_binding=True,
                                source='current_frozen_QC_marker_and_full_GEX'))
    mapping=mapping.merge(pd.concat(sources),on=['study_id','cluster_id'],how='left',validate='one_to_one')
    mapping['original_annotation_preserved']=True
    composition=cells.groupby(['original_state_id','library_key','experimental_time','condition'],dropna=False).size().rename('cells').reset_index()
    for state,g in composition.groupby('original_state_id'):
        j=mapping.state_id.eq(state)
        t=g.groupby('experimental_time').cells.sum().sort_values(ascending=False)
        mapping.loc[j,'time_composition']=';'.join(f'{k}:{v}' for k,v in t.items())
        mapping.loc[j,'library_count']=g.library_key.nunique()
    flow=pd.read_csv(flow_path,sep='\t')
    lookup=mapping.set_index('state_id').coarse_state
    for side in ['source','target']:
        flow[side+'_state']=flow[side+'_state_condition'].str.split('|').str[0]
        flow[side+'_condition']=flow[side+'_state_condition'].str.split('|').str[1]
        flow['coarse_'+side]=flow[side+'_state'].map(lookup)
        assert flow['coarse_'+side].notna().all()
    assert flow.pair.nunique()==34
    assert np.allclose(flow.groupby('pair').mass_fraction.sum(),1)
    assert not flow.source_state.str.startswith('GSE150311').any()
    keys=['pair','route_id','time_from','time_to','relation','coarse_source','coarse_target']
    coarse=flow.groupby(keys,dropna=False).mass_fraction.sum().reset_index()
    coarse['source_fraction']=coarse.mass_fraction/coarse.groupby(['pair','coarse_source']).mass_fraction.transform('sum')
    coarse['target_fraction']=coarse.mass_fraction/coarse.groupby(['pair','coarse_target']).mass_fraction.transform('sum')
    for side in ['source','target']:
        # A local state-condition count repeats across opposite targets; count it once.
        n=flow.drop_duplicates(['pair',side+'_state_condition']).groupby(['pair','coarse_'+side])[side+'_cells'].sum().rename(side+'_cells').reset_index()
        coarse=coarse.merge(n,on=['pair','coarse_'+side],validate='many_to_one')
    coarse['internal_to_coarse']=coarse.coarse_source.eq(coarse.coarse_target)
    coarse['edge_id']=[f'H{i+1:04d}' for i in range(len(coarse))]
    coarse['display_major']=((~coarse.internal_to_coarse)&(coarse.mass_fraction>=.05)&
        (coarse.source_fraction>=.15)&(coarse.source_cells>=50)&(coarse.target_cells>=50))
    coarse['display_rule']='cross_coarse; pair_mass>=0.05; source_fraction>=0.15; both_end_counts>=50'
    # Keep condition-specific source/target bins, especially D16 positive/negative sorting.
    ckeys=keys+['source_condition','target_condition']
    conditions=flow.groupby(ckeys,dropna=False).mass_fraction.sum().reset_index()
    conditions['source_condition_fraction']=conditions.mass_fraction/conditions.groupby(['pair','coarse_source','source_condition']).mass_fraction.transform('sum')
    conditions['target_condition_fraction']=conditions.mass_fraction/conditions.groupby(['pair','coarse_target','target_condition']).mass_fraction.transform('sum')
    summary=mapping.groupby(['coarse_state','coarse_label'],sort=False).agg(original_states=('state_id','size'),cells=('n_cells','sum'),studies=('study_id',lambda x:';'.join(sorted(set(x))))).reset_index()
    for name,table in [('local_to_coarse_state',mapping),('annotation_binding_checks',pd.DataFrame(binding)),
                       ('state_library_time_composition',composition),('all_local_state_flows',flow),
                       ('all_coarse_time_flows',coarse),('condition_preserving_coarse_flows',conditions),
                       ('major_coarse_time_flows',coarse[coarse.display_major]),('coarse_state_summary',summary)]:
        table.to_csv(D/(name+'.tsv'),sep='\t',index=False)
    assert np.allclose(coarse.groupby('pair').mass_fraction.sum(),1)
    assert np.allclose(conditions.groupby('pair').mass_fraction.sum(),1)
    review=dict(status='DATA_PASS_PENDING_FIGURE_REVIEW',cells=len(cells),original_states=len(mapping),
        exact_existing_annotation_cells=sum(b['cells'] for b in binding if b['source']=='existing_working_annotation'),
        original_time_pairs=34,all_coarse_time_rows=len(coarse),internal_rows=int(coarse.internal_to_coarse.sum()),
        major_external_rows=int(coarse.display_major.sum()),source_hashes=hashes,
        no_cross_pair_probability_sum=True,no_new_OT_or_embedding_fit=True,
        no_confirmatory_identity_upgrade=True,no_cross_study_equivalence_test=True)
    (D/'data_review.json').write_text(json.dumps(review,ensure_ascii=False,indent=2)+'\n')
    print(summary.to_string(index=False))
    print(coarse[coarse.display_major][['edge_id','route_id','time_from','time_to','coarse_source','coarse_target','mass_fraction','source_fraction']].round(3).to_string(index=False))


if __name__=='__main__':main()
