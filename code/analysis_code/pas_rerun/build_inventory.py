#!/usr/bin/env python3
"""Build the dated PAS rerun inventory from the active split manuscript only."""
import csv
import hashlib
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MOD = 'github/reproapa-analysis/modules/'
MOUSE = '../科研/projects/2026-09-02_mouse-time-aware-integration-pilot/scripts/'
# The separately running SCAPE launcher is not a native-rerun input. Its original
# inventory hash remains in baseline-reconciliation.json and legacy-entrypoints.
SEPARATE_CALLER_LAUNCHERS = {'notes/2026-10-02-scape-polyapipe-comparison/scripts/run_scape_atlas.py'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_tsv(name, rows):
    with (HERE / name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)


# IDs are stable analysis identities, independent of changing figure numbers.
# Programs below are migration sources, NOT executable adapters to new inputs.
DEFINITIONS = [
('pas_qc', '', 'PAS目录、支持率、序列信号、读端轮廓和过滤审计', 'catalog;matrices;calls;genome;references;read_evidence', 'supported50与expanded_primary；按site/UMI分别统计；保留过滤失败和投影零', [MOD+'05_pas_catalog/code/reliability']),
('recount', 'pas_qc', '原流程六规则计数的历史审计', 'frozen_historical_recount;scope_declaration', '用户选择仅保留历史审计；不做新的BAM重计数；所有输出明确标记非本轮polyApipe证据', ['notes/2026-09-26-endpoint-counting']),
('metacell', 'pas_qc', '按固定成员重新聚合PAS及计数守恒', 'matrices;cells;memberships;gex', '成员不变时复用SEACells，重建PAS覆盖；不得复用旧PAS-qualified成员', [MOD+'06_metacell_wpp/code/human/aggregate_human_metacell_pas.py', MOD+'06_metacell_wpp/code/mouse/aggregate_metacell_pas_v7.py']),
('state_vectors', 'metacell', '原生PAS轴、完整比例向量和状态均值', 'catalog;metacell_counts;state_registry', 'DWUI轴由完整primary目录定义；supported50筛选后不重新归一坐标；保留零计数PAS', [MOD+'07_apa_dynamics/code/figure2/compute_apa_node_position.py', 'notes/2026-09-30-multimetric-statistics/extract.py']),
('endpoints', 'state_vectors', '重建全部合格library pair、P–S/P–M及权重', 'state_vectors;library_design;contrast_design', 'observed-S/documented-source/paired P–M分开；按gene重算资格与日期/单元层级权重', ['notes/2026-09-29-core-sensitivity/export_endpoint_route_wpp.py', 'notes/2026-09-26-state-comparability-audit']),
('global', 'endpoints recount', '总体方向、固定背景、来源和覆盖敏感性', 'endpoint_effects;recount_effects;gex', 'available/intersection/common-contributor各自重算；purity .8/.9、≥3/5 MC、supported/expanded', ['notes/2026-09-26-global-direction-audit/analyze.py', 'notes/2026-09-29-core-sensitivity/summarize_route_wpp.py']),
('tandem', 'global', 'tandem 3′UTR与nt尺度位置变化', 'catalog;coding_transcripts;global_effects', '所有primary PAS须落在同一含stop codon终末外显子3′UTR；ΔDWUI×native span', ['notes/2026-09-26-global-direction-audit/tandem_utr.py']),
('metrics', 'endpoints', 'DWUI/RWUI/PDUI/PPUI/MPRO', 'qualified_pairs;complete_pas_vectors', 'MPRO=range(加权后Δp)，非先算range再平均；RWUI秩i/n；0.02/0.05与0.10/0.20', ['notes/2026-09-30-multimetric-statistics/extract.py', 'notes/2026-09-30-multimetric-statistics/summarize.py']),
('profiles', 'state_vectors', '重选轨迹候选、k=2–6评估、k=3分类、L/D/T', 'full_state_profiles;paired_studies;gex', '重新应用range≥.01和paired r≥.6；540/897为旧结果，不作新断言；类标按新centroid排序', [MOD+'07_apa_dynamics/code/figure2/cluster_apa_human_mouse_joint_v1.py', MOD+'07_apa_dynamics/code/figure2/analyze_figure2_dynamic_architecture_v1.py']),
('profile_controls', 'profiles', '轨迹扰动、深度、选样感知噪声模型', 'full_eligible_profiles;metacell_counts;profile_membership', '80%子采样24次/库权重32次/seed/layer/leave-route；20/40/80 UMI；500全背景单调模拟重选与聚类', [MOD+'07_apa_dynamics/code/figure2/prepare_figure2_core_supplement_v1.py', 'notes/2026-10-02-reviewer-followup/trajectory_noise.py']),
('moran', 'profiles', '表达图上的APA组织及coverage校正', 'gex_graph;class_scores;metacell_coverage', 'library-state residual置换；k=10/15/20/30；70/80/90%固定gene panel；各检验族分开', [MOUSE+'test_figure2_apa_manifold_coherence.py', MOUSE+'robustness_figure2_apa_manifold.py', MOUSE+'robustness_figure2_apa_knn.py', MOUSE+'project_k3_group_wpp_umap.py']),
('profile_gex', 'profiles', '轨迹APA–GEX相关与相邻转变同向率', 'state_profiles;matched_gex', '相同状态/APA合格MC重算GEX；zero变化保留分母', [MOD+'07_apa_dynamics/code/figure2/analyze_figure2_cde_dynamics_v1.py']),
('recurrence', 'metrics', '物种内多指标复现、LOSO与study票', 'metric_effects;eligibility;matching_covariates', '每GEO一次票且all-evaluable-series一致；10k置换；8/metric与24跨metric BH；全背景含non-event', [MOD+'08_recurrence/code/reference_wpp_null/analyze.py', 'notes/2026-09-30-multimetric-statistics/nulls/analyze.py']),
('consensus', 'recount', '原流程六规则复现与联合共识的历史审计', 'frozen_historical_consensus;scope_declaration', '保留原始固定基因/共同贡献者定义及零分布，明确标记historical，不纳入本轮native显著性族', ['notes/2026-09-29-core-sensitivity/figure3_recount_recurrence.py', 'notes/2026-09-29-core-sensitivity/figure3_joint_consensus_null.py']),
('measurement', 'recurrence', 'bootstrap、same-state、等深度、cycle、common-axis控制', 'all_qualified_pairs;metacell_vectors;cycle_scores', '500 bootstrap共享source draw；100 halves；20/40UMI各20次；exact series identity；112 planned positions含不可用P=1', ['notes/2026-10-02-reviewer-followup/apa_followup.py', 'notes/2026-10-02-reviewer-followup/same_state.py']),
('recurrence_controls', 'recurrence', '嵌套匹配与复现统计功效', 'full_event_labels;reference_capacity;covariates', 'exact study/series→expression/coverage→baseline→PAS；PAS-first路径；保持margin信号注入', ['notes/2026-10-02-null-model-controls/recurrence.py']),
('well_inference', 'metrics', 'Schiebinger独立well的P–M五指标与P–S推断', 'well_design;qualified_vectors', '重选≥4日期/各≥2well；同步sign作用完整向量；BH=5×新n，不能固定1340；P–S来源方差只进入一次', ['notes/2026-09-30-multimetric-statistics/inference/analyze.py', 'notes/2026-09-30-multimetric-statistics/metric_inference/analyze.py']),
('pool_inference', 'metrics', 'Panariello独立pool推断', 'pool_design;qualified_vectors', 'technical captures先合并；≥4pool/≥2日期各2pool；全部sign；报告最小可达P/q', ['notes/2026-09-30-multimetric-statistics/pool_inference/analyze.py']),
('count_inference', 'endpoints', '条件PAS计数组成推断', 'integer_pair_counts', '2×n Fisher及131071抽样；剔双零列仅限test；按comparison重建BH族与MC误差', ['notes/2026-09-30-multimetric-statistics/pair_count_tests/analyze.py']),
('matched_pms', 'recurrence well_inference', 'DWUI-only/shared/MPRO-only及共同P/M/S', 'new_recurrent_sets;complete_triplets;vectors', '新成员主分析与旧成员migration分开；相同贡献者权重；验证向量P−S=(P−M)+(M−S)；MPRO本身不满足加法', ['notes/2026-10-01-figure3-matched-metrics/extract_matched.py', 'notes/2026-10-01-figure3-matched-metrics/summarize_states.py', 'notes/2026-10-01-dwui-only-followup/analyze.py', 'notes/2026-09-27-source-state-connection/analyze.py']),
('apa_gex', 'recurrence', 'APA相同/GEX相反、双方向零分布', 'recurrent_effects;matched_gex;covariates', 'response-identity与direction-only分开；GEX .25；取两shuffle更大P；5bin/3bin；8-test族', [MOD+'08_recurrence/code/common_reference_20260925/analysis/figure3-recurrence-common-reference-2026-09-25/analyze.py']),
('intermediate', 'recurrence', 'R/E中间状态回溯与成员敏感性', 'endpoint_members;state_vectors;expression_support', '重建admission与R/E各自分母；新主队列/本轮固定队列敏感性；missing非零', [MOD+'08_recurrence/code/intermediate_states', 'notes/2026-09-26-state-comparability-audit/analyze.py']),
('functional', 'recurrence', 'GO/KEGG嵌套背景、精确尾概率与功效', 'G0_expression;G1_testable;G2_recurrent;annotations', '重建G0/G1/G2、strata与全term族；精确hypergeometric卷积优先；28 analyses/database校正；版本锁定pandas', ['github/reproapa-analysis/revisions/20260930/external/2026-09-26_reproapa-go-analysis/scripts', 'notes/2026-10-02-null-result-diagnostics/functional_exact.py', 'notes/2026-10-02-null-model-controls/functional.py']),
('external_pas', 'consensus recurrence', 'GSE99922独立PAS-seq比较', 'external_pas_seq;new_and_old_candidate_rosters;reference_catalog', '2MEF+2untransfected iPSC；native/common vectors；outcome-blind唯一matching；旧roster另报失效', ['notes/2026-10-02-gse99922-pilot/expanded_measurements.py', 'notes/2026-10-02-gse99922-pilot/vector_analysis.py', 'notes/2026-10-02-gse99922-pilot/compare_original_vectors.py']),
('cross_continuous', 'recurrence', '跨物种连续效应、匹配状态、复现集外推', 'orthologs;state_effects;covariates', 'one-to-one背景；ortholog block；复用control按unique ortholog处理；matched12-test族', [MOD+'09_cross_species/code/figure4/prepare_two_contrasts_v12.py', MOD+'09_cross_species/code/figure4/state_specific_recurrence_v25.py', MOD+'09_cross_species/code/figure4/audit_figure4c_reused_controls_v53.py']),
('cross_overlap', 'recurrence', 'observed-S P–S/paired P–M近端共享集', 'orthologs;state_effects;native_metrics', '重新筛选成员；旧18/10/26仅migration；16,092为旧映射背景记录；20k node-reference置换', ['notes/2026-09-29-core-sensitivity/figure4_matched_overlap.py', MOD+'09_cross_species/code/figure4/ortholog_pair_permutation_v24.py']),
('cross_sharing', 'recurrence', 'documented-source四study pair多指标共享', 'orthologs;documented_source_metrics;covariates', 'Q5/R10；DWUI/RWUI/MPRO主24族；DWUI8探索/24含方向族/全56敏感性须独立ID', ['notes/2026-09-30-multimetric-statistics/nulls/analyze.py']),
('cross_controls', 'cross_overlap', '跨物种嵌套条件零模型、校准与功效', 'new_fixed_cohort;labels;capacity;covariates', '在本轮新cohort内固定；20k null/3000 calibration；OR1/2/4/8；2 contrasts/model BH', ['notes/2026-10-02-null-model-controls/cross_species.py']),
('cross_trace', 'cross_overlap global', '共享成员跨背景/来源/setting追踪', 'new_shared_sets;global_backgrounds;state_effects', '重选结果决定行数；新set和旧18/10迁移表各自完整；不强行保留不支持的实例', ['notes/2026-10-01-figure4-state-pas-v6/prepare.py', MOD+'09_cross_species/code/figure4/prepare_compact_tables_v54.py']),
('cis', 'cross_continuous cross_overlap', 'PAS同源、driver特征、Firth模型与CV', 'orthologs;pas_drivers;genomes;chains;phastcons;motifs', '重新定义driver和特征；32 configurations/4999bootstrap；fold内impute/scale；10×5CV与2000gene uncertainty', [MOD+'09_cross_species/code/cis_supplement/scripts/prepare_data.py', MOD+'09_cross_species/code/cis_supplement/scripts/recompute_models.py', MOD+'09_cross_species/code/cis_supplement/upstream_scripts/cis-architecture-v27/03b_driver_annot.py', 'notes/2026-10-02-null-result-diagnostics/other_diagnostics.py']),
('dwui_pool', 'state_vectors', 'Figure5独立DWUI .01固定多数候选池', 'coarse_state_registry;full_state_effects', '按预定study数3/4或2/2投票；missing/conflicting abstain；与Figure3 .02复现分开', [MOD+'10_regulatory_sites/code/regulatory_sites/candidate-pool-v1/build_pool.py']),
('union_pool', 'dwui_pool recurrence', 'DWUI–MPRO并集和加入来源', 'dwui_pool;MPRO_0p10_recurrence', '按contrast追踪basis；2149/2343不是固定输入', ['notes/2026-10-01-mpro-biology/prepare_union.py']),
('intervals', 'union_pool state_vectors', 'terminal-family映射、原子区间和inclusion', 'catalog;terminal_annotations;new_candidate_pools;counts', 'native→annotated identity ≤50nt；family合并≥2PAS；20/50UMI和coverage .5；P–M common library；included/excluded5', ['notes/2026-10-01-mpro-biology/annotate_union.py', MOD+'10_regulatory_sites/code/regulatory_sites/annotation-pas-global-v1/migrate.py']),
('site_scan', 'intervals', 'miRNA seed、direct RBP和transfer扫描', 'reference_UTRs;seed_sets;motif_banks;genome', '完整新并集重新scan；FIMO P<1e-4 forward-only；exact物理site去重；direct/transfer分开', ['notes/2026-10-01-figure5-union-update/scan_union.py']),
('site_statistics', 'site_scan', 'union位点/region复现、一致性和部分保留', 'interval_effects;site_hits;family_models', '新gene/region分母；5000/2000 gene bootstrap；19999 permutation；小样本描述；common与separate规则不混合', ['notes/2026-10-01-figure5-union-update/analyze_union.py', MOD+'10_regulatory_sites/code/regulatory_sites/site-recurrence-v1']),
('original_pool_sites', 'dwui_pool site_scan', 'DWUI-only池的S10–S15配套统计', 'new_dwui_pool;site_hits;separate_state_inclusion', '重跑原设计的DWUI-only池；另存历史固定roster migration；不得用union覆盖Tables S10–S15', [MOD+'10_regulatory_sites/code/regulatory_sites/full-cohort-v2/run_full.py', MOD+'10_regulatory_sites/code/regulatory_sites/annotation-pas-state-revision-v1']),
('coordinate_states', 'original_pool_sites', 'native/annotated坐标及R/E成员审计', 'source_to_catalog;family_counts;membership_variants', '记录距离/多对一/onePAS collapse；C012/L00 exclusions；first-crossing5/10/15/20pp；不可测状态不补', ['notes/2026-09-29-core-sensitivity/figure5_coordinate_audit.py', 'notes/2026-09-29-core-sensitivity/figure5_paired_pm.py']),
('motif_background', 'site_statistics', '同gene稳定窗口匹配和density', 'changed_intervals;stable_intervals;sequence;site_hits', '等长；GC≤.05、position≤.20；5000bootstrap/10k sign；8族；<10gene保留P1', ['notes/2026-10-02-reviewer-followup/motif_background.py']),
('motif_controls', 'motif_background', '匹配网格、跨gene背景与检测功效', 'matched_windows;site_hits;all_changed_intervals', '12settings/96族；transfer24；6background48；LOO；0/2/5/10/20sites/kb；full/subset injection', ['notes/2026-10-02-null-result-diagnostics/motif_sensitivity.py', 'notes/2026-10-02-null-result-diagnostics/transfer_sensitivity.py', 'notes/2026-10-02-null-model-controls/motif.py']),
('weighted_retention', 'motif_controls', '按ΔI加权的预测位点净保留', 'all_strict_changed_intervals;complete_site_hits', '匹配与未匹配均保留；W=sum(n_i ΔI_i)；20/50层；gene正/负/零分开', ['notes/2026-10-02-null-model-controls/summarize.py']),
('rbp_expression', 'site_statistics', '重选RBP候选并提取表达上下文', 'recurrent_factor_sites;all_gene_gex;memberships', '候选随PAS更新；所有鼠direct及human top6按target数预选；CPM≥2仅flag；transfer另列', ['notes/2026-10-02-rbp-expression/collect.py', 'notes/2026-10-02-rbp-expression/extract_remote.py', 'notes/2026-10-02-rbp-expression/analyze.py']),
('examples', 'matched_pms cross_trace site_statistics well_inference', '全部主图实例PAS曲线、MC分布与资格', 'new_candidate_sets;native_vectors;family_inclusion', 'Slc38a2/Nudt21/SZRD1/RABIF/Rtf2/AP1S2/Hmgn2/Eif4ebp2/Rad23b逐一复核；失败保留reason，不套旧数', ['notes/2026-10-01-figure4-state-pas-v6/prepare.py', 'notes/2026-10-01-figure-revision-drafts/extract_case_metacells.py']),
('longread', 'examples', '外部长读长terminal结构重新相交', 'current_features;GSE282319_gtf_tpm_sqanti;chain;genome', '资源可复用但新区间必须重新相交；唯一full-length reciprocal同链；原3基因另存追踪', ['notes/2026-10-02-GSE282319-long-read/analyze.py', 'notes/2026-10-02-GSE282319-long-read/supplement.py']),
('ago2', 'examples', 'Rad23b区间/seed的AGO2 overlap', 'current_rad23b_interval_seed;GSE25310_tracks;chains;genomes', '不锁旧415nt/145bp；feature消失要明确报告；same-strand≥1bp；WT1技术库同一生物单元', ['notes/2026-10-02-rad23b-ago2-gse25310/check_overlap.py']),
('caller_comparison', 'recurrence examples', '预留SCAPE-APA/polyApipe全七研究比较', 'alternative_caller_catalog_counts;frozen_legacy_effects;corrected_effects', 'native和common PAS/MC/pairs两层；25/50/100nt双分母；not-run≠not-measurable；旧list迁移与新list分开', ['notes/2026-10-02-scape-polyapipe-comparison', 'notes/2026-10-03-polyapipe-asnr-junction-postfilter/analyze.py']),
]

FIG_NODES = {
 'Figure 1': 'global state_vectors consensus', 'Figure 2': 'profiles profile_controls moran profile_gex',
 'Figure 3': 'recurrence consensus measurement recurrence_controls well_inference matched_pms examples',
 'Figure 4': 'cross_continuous cross_overlap cross_sharing cross_trace examples',
 'Figure 5': 'site_statistics weighted_retention rbp_expression examples',
 'Figure S1': '', 'Figure S2': '', 'Figure S3': '',
 'Figure S4': 'pas_qc recount profile_controls', 'Figure S5': 'caller_comparison',
 'Figure S6': 'global tandem metrics', 'Figure S7': 'profiles profile_controls',
 'Figure S8': 'moran profile_gex', 'Figure S9': 'measurement recurrence_controls recurrence',
 'Figure S10': 'well_inference pool_inference count_inference', 'Figure S11': 'matched_pms examples',
 'Figure S12': 'apa_gex examples', 'Figure S13': 'intermediate', 'Figure S14': 'functional',
 'Figure S15': 'cross_continuous', 'Figure S16': 'cross_sharing cross_overlap cross_controls',
 'Figure S17': 'cross_trace', 'Figure S18': 'cis',
 'Figure S19': 'weighted_retention union_pool site_statistics', 'Figure S20': 'motif_background motif_controls',
 'Figure S21': 'rbp_expression', 'Figure S22': 'coordinate_states site_statistics longread',
 'Figure S23': 'examples longread ago2',
}


def source_files():
    out = []
    for doc in ['main', 'supplement']:
        manifest = ROOT / doc / '.texsplit.json'
        m = json.loads(manifest.read_text())
        out.append(manifest)
        for item in m['backbone'] + m['sections']:
            out.append(ROOT / doc / item['file'])
            out.extend(ROOT / doc / s['file'] for s in item.get('subsections', []))
    return out


def main():
    nodes, programs = [], []
    for node, deps, title, inputs, rules, sources in DEFINITIONS:
        files = set()
        for name in sources:
            p = ROOT / name
            if not p.exists():
                raise ValueError(f'Missing registered migration source: {name}')
            files.update([p] if p.is_file() else (q for q in p.rglob('*') if q.suffix in {'.py', '.sh', '.R'} and not any(x in q.parts for x in ['before', '__pycache__', 'data', 'original'])))
        nodes.append(dict(id=node, depends_on=deps.split(), title=title, inputs=inputs.split(';'), rules=rules,
                          sources=sources, adapter_status='implemented_fixture_tested_runtime_unbound' if node=='metacell' else 'unbound', products=['result.json']))
        for p in sorted(files):
            t = p.read_text(errors='replace')
            programs.append(dict(node=node, file=str(p.relative_to(ROOT)), sha256=sha(p),
                                 absolute_path_literals=len(re.findall(r'[\'\"]/(?:Users|home|data|mnt|work)/', t)),
                                 assertion_lines=';'.join(str(i) for i, line in enumerate(t.splitlines(), 1) if re.search(r'\bassert\b', line)),
                                 argparse='yes' if 'argparse' in t else 'no', status='migration_source_not_adapter'))
    all_ids = [n['id'] for n in nodes]
    nodes += [dict(id='render', depends_on=all_ids, title='从本轮panel-data重绘所有受影响主补图', inputs=['panel_data'], rules='不得裁切旧PAS图层；同一run_id；源数据/样本量/Pq从表读取；空结果有效但须明确显示', sources=['notes/2026-10-02-main-figure-revision', 'notes/2026-10-03-supplement-claude-dense'], adapter_status='unbound', products=['result.json']),
              dict(id='tables', depends_on=all_ids, title='重新生成补表S8–S30及完整source-data', inputs=['numeric_results'], rules='S1–S7条件复用；S17规则+新n；S10–S15 DWUI-only，S19 union；S29全精度记录；S30不得假完成', sources=['supplement/supplementary-tables/main.tex', 'table/README.md'], adapter_status='unbound', products=['result.json']),
              dict(id='text', depends_on=['render', 'tables'], title='摘要/正文/讨论/结论/图注/补文同步', inputs=['claim_ledger', 'numeric_results'], rules='逐项处理text-targets.tsv；结果改变则重新写结论，禁止全局数字替换', sources=['main', 'supplement'], adapter_status='unbound', products=['result.json']),
              dict(id='release', depends_on=['text'], title='组装候选稿、结构/编译/视觉检查和结果清单', inputs=['graphics', 'tables', 'split_sources'], rules='candidate目录完整main/supplement/Fig/table；make check+make pdf；页图复核；receipt齐全方能发布', sources=['scripts/tex_sections.py', 'Makefile'], adapter_status='unbound', products=['result.json'])]
    (HERE / 'workflow.json').write_text(json.dumps(dict(schema_version=1, baseline_date='2026-10-03', scope='current manuscript; PAS replacement; expression framework frozen unless explicitly changed', nodes=nodes), ensure_ascii=False, indent=2)+'\n')
    write_tsv('analyses.tsv', [dict(id=n['id'], depends_on=';'.join(n['depends_on']), analysis=n['title'], inputs=';'.join(n['inputs']), rules=n['rules'], migration_sources=';'.join(n['sources']), adapter_status=n['adapter_status']) for n in nodes])
    write_tsv('legacy-entrypoints.tsv', programs)
    figures = list(csv.DictReader((ROOT/'Fig/figure-manifest.tsv').open(), delimiter='\t'))
    mapping = {f['file']: f for f in figures}
    figure_rows, panel_rows, text_rows = [], [], []
    for f in figures:
        assert f['figure'] in FIG_NODES
        path = ROOT/'Fig'/f['file']
        assert sha(path) == f['sha256'], f['file']
        figure_rows.append(dict(figure=f['figure'], target='Fig/'+f['file'], action='rebuild' if FIG_NODES[f['figure']] else 'reuse_if_expression_hashes_unchanged', analyses=FIG_NODES[f['figure']].replace(' ', ';'), baseline_sha256=f['sha256'], title=f['title']))
    for p in source_files():
        if p.suffix != '.tex':
            continue
        t = p.read_text()
        for m in re.finditer(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}', t):
            f = mapping[Path(m[1]).name]
            # ContinuedFloat captions until the next image belong to the same asset.
            end = t.find('\\includegraphics', m.end())
            block = t[m.end():end if end >= 0 else len(t)]
            keys = sorted(set(re.findall(r'\\textbf\{([A-Z](?:,[A-Z]|[–-][A-Z])?)\}', block)))
            label = re.search(r'\\label\{([^}]+)\}', block)
            for key in keys or ['whole']:
                panel_rows.append(dict(figure=f['figure'], panel=key, file='Fig/'+f['file'], tex=str(p.relative_to(ROOT)), label=label[1] if label else '', analyses=FIG_NODES[f['figure']].replace(' ', ';'), replacement='new_run_panel_data;expand_combined_panel_letters_where_present'))
        # Whole paragraphs, with hashes rather than brittle automatic number replacement.
        for match in re.finditer(r'\S[\s\S]*?(?=\n\s*\n|\Z)', t):
            paragraph = match[0].strip()
            is_result = any(s in str(p.relative_to(ROOT)) for s in ['/results/', '/discussion/', '/conclusions/', '/supplementary-results/', '/supplementary-tables/', '/supplementary-figures/'])
            if not (re.search(r'PAS|APA|DWUI|MPRO|RWUI|PDUI|PPUI|polyadenyl|miRNA|RBP|Rad23b|Eif4ebp2|Hmgn2', paragraph, re.I) or (is_result and re.search(r'\d', paragraph))):
                continue
            text_rows.append(dict(id=f'T{len(text_rows)+1:04}', file=str(p.relative_to(ROOT)), line=t.count('\n', 0, match.start())+1, sha256=hashlib.sha256(paragraph.encode()).hexdigest(), action='review_or_replace_against_new_results', text=re.sub(r'\s+', ' ', paragraph)))
    write_tsv('figures.tsv', figure_rows)
    write_tsv('panels.tsv', panel_rows)
    write_tsv('text-targets.tsv', text_rows)
    old_panels = list(csv.DictReader((ROOT/'notes/2026-10-03-supplement-claude-dense/source/new-panel-index.tsv').open(), delimiter='\t'))
    for r in old_panels:
        n = int(re.search(r'\d+', r['panel'])[0]); r['current_panel'] = re.sub(r'^S\d+', 'S'+str(n+(n>=5)), r['panel'])
        r['warning'] = 'historical numeric/rendering trace only; new data required'
    write_tsv('supplement-panel-lineage.tsv', old_panels)
    baseline_paths = set(source_files()) | {ROOT/'Fig/figure-manifest.tsv', ROOT/'references.bib', ROOT/'Makefile', ROOT/'scripts/tex_sections.py'}
    baseline_paths.update(ROOT/r['target'] for r in figure_rows)
    baseline_paths.update(ROOT/p['file'] for p in programs if p['file'] not in SEPARATE_CALLER_LAUNCHERS)
    write_tsv('baseline-files.tsv', [dict(file=str(p.relative_to(ROOT)), sha256=sha(p)) for p in sorted(baseline_paths)])
    (HERE/'bindings.template.json').write_text(json.dumps({'schema_version': 1, 'nodes': {n['id']: {'argv': [], 'code_files': [], 'note': 'Bind an adapted command that reads PAS_INPUT_MANIFEST and PAS_DEPENDENCIES_JSON and writes PAS_NODE_DIR/result.json. Historical scripts cannot be run unchanged.'} for n in nodes}}, indent=2)+'\n')
    print(json.dumps(dict(analyses=len(all_ids), workflow_nodes=len(nodes), figures=len(figure_rows), panel_records=len(panel_rows), text_targets=len(text_rows), source_records=len(programs)), ensure_ascii=False))


if __name__ == '__main__':
    main()
