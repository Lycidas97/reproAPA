"""Verify the downloaded receipt, then expose the five-model results for review."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

root = Path(__file__).resolve().parent
snapshot = root / 'results'
manifest = json.loads((snapshot / 'HANDOFF-MANIFEST.json').read_text())
assert manifest['status'] == 'PASS'
assert manifest['payload'] == 'aggregate_statistics_and_verification_reports'
for source, identity in manifest['sources'].items():
    original = json.loads((root / (source + '.source.config.json')).read_text())
    assert identity == {key: original[key] for key in identity}, source
for item in manifest['files']:
    path = snapshot / item['path']
    assert path.stat().st_size == item['bytes'], item['path']
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    assert h.hexdigest() == item['sha256'], item['path']
labels = {'joint_measurement_PAS': 'Joint measurement/PAS',
          'independent_q5_reference_capacity': 'Quintile',
          'independent_range10_reference_capacity': 'Range',
          'Q5_expression_coverage': 'Quintile + expression/coverage',
          'Range_expression_coverage': 'Range + expression/coverage'}
parts = []
lines = ['# A1：五个跨物种共享近端候选 null 的修正重跑', '',
         '2026-10-07。直接修改原统计实现，先完成 polyApipe，随后完成 APADETECT 点过滤＋p50。', '',
         '在所有旧分层条件上增加实际可测效应研究数；保留 16,092 对背景、各来源的原效应、原候选标签及成员。'
         '人类和小鼠实际可测研究数小于二的基因，在五个 null 中的正事件概率均严格为零。', '',
         '原三个模型仍各自按两个对比计算 q₂；两扩展模型保留四检验 q₄，没有事后合并成一个十检验族。'
         'Joint 使用 5,000 次，其余模型各 20,000 次置换；扩展的 OR 1/2/4 校准与功效各用 3,000 次模拟。', '']
for source in ['polyapipe', 'apadetect-p50']:
    folder = snapshot / source
    for receipt in [folder / 'all-five-verification.json', folder / 'independent-verification.json',
                    snapshot / 'validation' / (source + '-grid') / 'sensitivity-verification.json',
                    snapshot / 'validation' / (source + '-exact') / 'exact-verification.json']:
        assert json.loads(receipt.read_text())['status'] == 'PASS', receipt
    frame = pd.read_csv(folder / 'all-five-null-comparison.tsv', sep='\t', float_precision='round_trip')
    frame.insert(0, 'source', source)
    exact = pd.read_csv(snapshot / 'validation' / (source + '-exact') / 'joint-exact-results.tsv',
                        sep='\t', float_precision='round_trip')
    exact_primary = exact[exact.layer.eq('supported50') & exact.threshold.eq(.02)].set_index('contrast')
    joint = frame.model.eq('joint_measurement_PAS')
    frame['joint_exact_p'] = np.nan
    frame['joint_exact_q2'] = np.nan
    frame.loc[joint, 'joint_exact_p'] = frame.loc[joint, 'contrast'].map(exact_primary.exact_p)
    frame.loc[joint, 'joint_exact_q2'] = frame.loc[joint, 'contrast'].map(exact_primary.exact_q2_primary)
    parts.append(frame)
    lines += ['## ' + source, '',
              '| 对比 | 模型 | 观察 | 原条件期望 | 修正期望 | 倍数 | p | q₂ | q₄ |',
              '|---|---|---:|---:|---:|---:|---:|---:|---:|']
    for row in frame.itertuples(index=False):
        q4 = '—' if pd.isna(row.q_four_tests) else f'{row.q_four_tests:.5g}'
        lines.append(f'| {row.contrast} | {labels[row.model]} | {row.observed} | {row.legacy_analytic_expected:.3f} | '
                     f'{row.analytic_expected:.3f} | {row.enrichment:.3f} | {row.p:.5g} | {row.q_two_contrasts:.5g} | {q4} |')
    pm = frame[frame.contrast.eq('PM')]
    significant = pm.apply(lambda row: (row.q_four_tests if pd.notna(row.q_four_tests)
                                      else row.q_two_contrasts) < .05, axis=1)
    lines += ['', f'P–M：五个模型中 {int(significant.sum())}/5 在相应既定检验族下显著。', '',
              f'Joint 主阈值的精确超几何卷积复核：P–M 精确 p={exact_primary.loc["PM", "exact_p"]:.6g}、'
              f'q₂={exact_primary.loc["PM", "exact_q2_primary"]:.6g}。'
              '精确分布的归一化、均值及 SD 均与原解析公式相符，排除了显著性边界附近的置换抽样波动。', '',
              '旧条件向不合格基因分配的期望事件（P–M；新条件两侧全部为 0）：', '',
              '| 模型 | 人类 | 小鼠 |', '|---|---:|---:|']
    for row in pm.itertuples(index=False):
        lines.append(f'| {labels[row.model]} | {row.h_legacy_expected_ineligible_events:.4f} | '
                     f'{row.m_legacy_expected_ineligible_events:.4f} |')
    mobility = pd.read_csv(folder / 'mobility-diagnostics.tsv', sep='\t', float_precision='round_trip')
    balances = pd.read_csv(folder / 'inference/balance.tsv', sep='\t', float_precision='round_trip')
    power = pd.read_csv(folder / 'inference/calibration-power.tsv', sep='\t', float_precision='round_trip')
    lines += ['', 'P–M 扩展模型的修正后诊断：', '',
              '| 模型 | 人候选可移动 | 鼠候选可移动 | 固定交叠 | 最大绝对 SMD | OR1 拒绝率 | OR2 功效 | OR4 功效 |',
              '|---|---:|---:|---:|---:|---:|---:|---:|']
    for row in mobility[mobility.contrast.eq('PM')].itertuples(index=False):
        select = balances[balances.contrast.eq('PM') & balances.model.eq(row.model) &
                          balances.matching.eq('expression_coverage')]
        maximum = select.standardized_difference.abs().max()
        rates = power[power.contrast.eq('PM') & power.model.eq(row.model) & power.alpha.eq(.025)].set_index('odds_ratio').rejection_rate
        lines.append(f'| {labels[row.model]} | {row.h_candidate_movable_fraction:.1%} | '
                     f'{row.m_candidate_movable_fraction:.1%} | {row.fixed_overlap} | {maximum:.3f} | '
                     f'{rates.loc[1]:.1%} | {rates.loc[2]:.1%} | {rates.loc[4]:.1%} |')
    lines += ['', '可移动比例以候选为分母；SMD 包括表达、覆盖和最低端点/库覆盖的原尺度及 log₂(1+x) 尺度。'
              '校准/功效在 α=0.025 下计算；OR 是条件注入的优势比，不等于观察/期望倍数，功效也不是 BH 联合拒绝功效。', '']
    grid_check = json.loads((snapshot / 'validation' / (source + '-grid') / 'sensitivity-verification.json').read_text())
    grid = pd.read_csv(folder / 'original-nulls/overlap_threshold_layer_sensitivity.tsv', sep='\t', float_precision='round_trip')
    grid = grid.merge(exact[['contrast', 'layer', 'threshold', 'exact_p', 'exact_q16']],
                      on=['contrast', 'layer', 'threshold'], validate='one_to_one')
    exact_sig = int((exact.contrast.eq('PM') & exact.exact_q16.lt(.05)).sum())
    lines += ['### Joint 阈值/层敏感性', '',
              f'P–M 显著设置（q₁₆<0.05）从旧条件 {grid_check["old_PM_significant_settings"]}/8 '
              f'变为修正后的 {grid_check["corrected_PM_significant_settings"]}/8。'
              f'消除 Monte Carlo 波动后，精确分布为 {exact_sig}/8。'
              '该检验族与主阈值各模型的 q₂、扩展模型的 q₄ 分开解释。', '',
              '| 对比 | 层 | 阈值 | 观察 | 修正期望 | MC p | MC q₁₆ | 精确 p | 精确 q₁₆ |',
              '|---|---|---:|---:|---:|---:|---:|---:|---:|']
    for row in grid.itertuples(index=False):
        lines.append(f'| {row.contrast} | {row.layer} | {row.threshold:g} | {row.observed} | '
                     f'{row.analytic_null_mean:.3f} | {row.empirical_p:.5g} | {row.q_16:.5g} | '
                     f'{row.exact_p:.5g} | {row.exact_q16:.5g} |')
    lines += ['', '网格的全部候选成员、分层事件边际、16 组置换数组及 q₁₆ 已独立重放核验。', '']
combined = pd.concat(parts, ignore_index=True)
combined.to_csv(root / 'all-sources-five-null-comparison.tsv', sep='\t', index=False, float_format='%.17g')
assert not combined[['h_expected_ineligible_events', 'm_expected_ineligible_events']].to_numpy().any()
lines[4:4] = ['主阈值 P–M：polyApipe 五个模型均支持超额交叠；p50 四个支持，joint 的精确 q₂=0.0597，未达 0.05。'
              'Joint 网格的精确 q₁₆ 仅 polyApipe 3/8、p50 1/8 显著，旧的跨全部阈值稳健性表述不能保留。', '']
lines += ['## 解释及交付边界', '',
          '“原条件期望”通过去掉新增 actual_effect_n 条件重建。三个原模型均复现原已执行结果；'
          'polyApipe 扩展复现前次扩展条件，p50 扩展此前未执行，因此其旧条件期望仅是未修正设计的解析对照。', '',
          '此次修正消除了违反候选最低可测资格的置换。分箱内的连续值残余差异仍需按诊断解释；'
          '实际研究数分层也不等同于逐研究组合模式匹配。不能把该检验解释为完全控制全部研究特异可测性。', '',
          '原 joint 的阈值/层敏感性随原函数一并重算，保存在各来源 original-nulls/overlap_threshold_layer_sensitivity.tsv。'
          '另外的 L0–L3 梯度控制未重跑，其旧校准/功效不能作为此次 A1 修正后的结果。', '',
          '本轮未改正文、图或 Table S5b；其中旧数值及将该问题称为“已声明限制”的文字仍需同步。'
          '旧统计包及硬链接绘图数据没有覆盖。', '',
          '## 验证与复现', '',
          '七项意图测试通过，跳过数为零。两来源分别通过原贡献者/表达覆盖复核，以及五模型的独立效应计数、'
          '背景/成员、零事件支持、解析期望、数组尾概率、BH 检验族与输入/输出 SHA 复核。', '',
          '另独立核验两来源共 32 个 joint 阈值/层设置：全部候选成员保持不变，置换数组和 q₁₆ 从分层边际准确重放。', '',
          '对同一 joint 零模型另做精确超几何分布卷积，核验归一化、解析均值/SD 及精确尾概率，'
          '没有改动原 5,000 次 Monte Carlo 结果。', '',
          f'本地汇总包 {len(manifest["files"])} 个清单文件的大小及 SHA-256 全部通过；'
          '两来源的 run_id、input_fingerprint、reference_variant 与原配置完全相符。', '',
          '- 全来源统一结果：all-sources-five-null-comparison.tsv。',
          '- 原接口结果：results/{polyapipe,apadetect-p50}/original-nulls/overlap_null_results.tsv。',
          '- 扩展结果：results/{polyapipe,apadetect-p50}/inference/expression-coverage-null-results.tsv。',
          '- 本地校准/功效、balance、分箱和独立验证报告在 results/ 对应来源目录；精确结果在 results/validation/。',
          '- 本地原代码备份：before/；本轮原代码修改：correction.patch；本轮执行代码冻结副本：code/。',
          '- 当前 analysis_code/pas_rerun/downstream_cross.py 的五模型相关函数与冻结副本一致（code-consistency.json）；'
          '当前文件另含 L0–L3 控制的修改，不属于本轮执行记录或七项测试的验证范围。',
          '- 完整配置、逐阶段命令/status.json、environment.json、提取记录、置换数组和失败启动记录只留在服务器。',
          '- 本地下载清单：results/HANDOFF-MANIFEST.json；完整服务器包的清单 SHA-256：'
          + manifest['complete_server_package_manifest_sha256'] + '。', '',
          '服务器：213docker:/root/nfsdata/apadetect/analysis_projects/reproapa_cross_overlap_null_a1_20261007_v1/。'
          '复现时使用新输出目录；原生 APA/GEX 矩阵仍在原服务器，路径和 SHA 在 extraction/audit.json 中。', '',
          '整包下载被自动审批拒绝，原因是用户授权重跑未包含贡献者数据等完整导出。'
          '本地仅交付经白名单筛选的汇总统计及验证报告，不含细胞/基因标识、贡献者数据、计数矩阵或置换数组。', '']
(root / 'REPORT.zh.md').write_text('\n'.join(lines))
(root / 'download-verification.json').write_text(json.dumps(dict(status='PASS', files=len(manifest['files']),
    bytes=sum(item['bytes'] for item in manifest['files']), missing=0, hash_errors=0,
    source_identities_match=True, payload=manifest['payload'],
    statistical_tests=20, sensitivity_settings=32, intent_tests=7, tests_skipped=0), indent=2) + '\n')
print(combined[combined.contrast.eq('PM')][['source', 'model', 'observed', 'analytic_expected', 'q_two_contrasts', 'q_four_tests']].to_string(index=False))
