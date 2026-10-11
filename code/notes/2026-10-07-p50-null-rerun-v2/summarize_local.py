"""Check the downloaded p50 aggregate receipt and expose the rerun for review."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

root = Path(__file__).resolve().parent
data = root / 'results'
read = lambda path: pd.read_csv(path, sep='\t', float_precision='round_trip')
manifest = json.loads((data / 'HANDOFF-MANIFEST.json').read_text())
assert manifest['status'] == 'PASS' and manifest['tests_skipped'] == 0
assert manifest['identity']['input_fingerprint'] == 'f297dfffd8459d3223fcad4654ce8eb0ff566668877895a41a4e01056ad855df'
for relative, fingerprint in manifest['code_sha256'].items():
    assert hashlib.sha256((root / relative).read_bytes()).hexdigest() == fingerprint, relative
for item in manifest['files']:
    path = data / item['path']
    assert path.stat().st_size == item['bytes'], item['path']
    assert hashlib.sha256(path.read_bytes()).hexdigest() == item['sha256'], item['path']
folder = data / 'apadetect-p50'
models = read(folder / 'all-five-null-comparison.tsv')
earlier = read(root.parent / '2026-10-07-cross-overlap-null-a1/results/apadetect-p50/all-five-null-comparison.tsv')
pd.testing.assert_frame_equal(models, earlier, check_dtype=False, rtol=0, atol=1e-14)
exact = read(data / 'validation/exact/joint-exact-results.tsv')
primary = exact[exact.layer.eq('supported50') & exact.threshold.eq(.02)].set_index('contrast')
models['joint_exact_p'] = np.nan
models['joint_exact_q2'] = np.nan
joint = models.model.eq('joint_measurement_PAS')
models.loc[joint, 'joint_exact_p'] = models.loc[joint, 'contrast'].map(primary.exact_p)
models.loc[joint, 'joint_exact_q2'] = models.loc[joint, 'contrast'].map(primary.exact_q2_primary)
models.to_csv(root / 'five-null-results.tsv', sep='\t', index=False, float_format='%.17g')
assert not models[['h_expected_ineligible_events', 'm_expected_ineligible_events']].to_numpy().any()

mobility = read(folder / 'mobility-diagnostics.tsv')
balance = read(folder / 'inference/balance.tsv')
power = read(folder / 'inference/calibration-power.tsv')
diagnostic_rows = []
for row in mobility.itertuples(index=False):
    subset = balance[balance.contrast.eq(row.contrast) & balance.model.eq(row.model) &
                     balance.matching.eq('expression_coverage')]
    rates = power[power.contrast.eq(row.contrast) & power.model.eq(row.model) &
                  power.alpha.eq(.025)].set_index('odds_ratio')
    diagnostic_rows.append(dict(contrast=row.contrast, model=row.model,
        fixed_overlap=row.fixed_overlap, h_candidate_movable_fraction=row.h_candidate_movable_fraction,
        m_candidate_movable_fraction=row.m_candidate_movable_fraction,
        maximum_absolute_SMD=subset.standardized_difference.abs().max(), alpha=.025,
        OR1_rejection_rate=rates.loc[1, 'rejection_rate'], OR2_power=rates.loc[2, 'rejection_rate'],
        OR4_power=rates.loc[4, 'rejection_rate']))
diagnostics = pd.DataFrame(diagnostic_rows)
diagnostics.to_csv(root / 'extension-diagnostics.tsv', sep='\t', index=False, float_format='%.17g')

ladder = read(folder / 'ladder/cross_species_power.tsv')
earlier_ladder = read(root.parent / '2026-10-07-cross-null-followup/ladder-results/apadetect-p50/ladder/cross_species_power.tsv')
pd.testing.assert_frame_equal(ladder, earlier_ladder, check_dtype=False, rtol=0, atol=1e-14)
chosen = ladder[ladder.generator.eq('L3_full') & ladder.tested_model.eq('L3_full') &
                ladder.alpha.eq(.05) & ladder.experiment.eq('conditional_injection')]
eligibility = read(data / 'eligibility-verification.tsv')
assert not eligibility[['h_expected_ineligible_events', 'm_expected_ineligible_events']].to_numpy().any()
labels = {'joint_measurement_PAS': 'Joint', 'independent_q5_reference_capacity': 'Quintile',
          'independent_range10_reference_capacity': 'Range', 'Q5_expression_coverage': 'Quintile + 表达/覆盖',
          'Range_expression_coverage': 'Range + 表达/覆盖'}
lines = ['# p50 五模型及梯度校准/功效重跑', '',
         '已在全新服务器目录完整重跑 APADETECT 点过滤＋p50。未改稿件表述。'
         '原背景16,092对、候选成员及原效应保持不变，P–M 仍为12个共享候选（人类84、小鼠154）。', '',
         '主阈值 P–M：两个 reference-capacity null 及其表达/覆盖扩展均显著，joint 未通过校正。'
         f'Joint MC q₂={primary.loc["PM", "monte_carlo_p"] * 2:.6g}，精确 q₂={primary.loc["PM", "exact_q2_primary"]:.6g}。', '',
         '| 对比 | 模型 | 观察 | 解析期望 | 倍数 | MC p | q₂ | q₄ |',
         '|---|---|---:|---:|---:|---:|---:|---:|']
for row in models.itertuples(index=False):
    q4 = '—' if pd.isna(row.q_four_tests) else f'{row.q_four_tests:.5g}'
    lines.append(f'| {row.contrast} | {labels[row.model]} | {row.observed} | {row.analytic_expected:.6f} | '
                 f'{row.enrichment:.4f} | {row.p:.5g} | {row.q_two_contrasts:.5g} | {q4} |')
lines += ['', '原三个模型保留各自两对比 BH 族；扩展另报告两对比 q₂ 和四检验 q₄。'
          'Joint 5,000次，其他四模型各20,000次。Joint 的精确超几何卷积另列，未覆盖 MC 结果。', '',
          '## 从服务器下载的扩展诊断', '',
          '| 对比 | 模型 | 固定交叠 | 人可移动 | 鼠可移动 | 最大绝对 SMD | OR1 拒绝率 | OR2 功效 | OR4 功效 |',
          '|---|---|---:|---:|---:|---:|---:|---:|---:|']
for row in diagnostics.itertuples(index=False):
    lines.append(f'| {row.contrast} | {labels[row.model]} | {row.fixed_overlap} | '
                 f'{row.h_candidate_movable_fraction:.1%} | {row.m_candidate_movable_fraction:.1%} | '
                 f'{row.maximum_absolute_SMD:.6f} | {row.OR1_rejection_rate:.3%} | '
                 f'{row.OR2_power:.4f} | {row.OR4_power:.4f} |')
lines += ['', '可移动比例以候选为分母；SMD 包括表达、APA覆盖及最低端点/库覆盖的原尺度和log₂(1+x)尺度。'
          '此表校准/功效为α=0.025；完整文件含各α及95%区间。OR是注入优势比，不是观察/期望倍数。', '',
          '## 带资格约束的梯度重跑', '',
          '每一级 L0–L3 和 L3 注入均加入各物种真实可测效应研究数。原背景保留，但不足两项研究者的事件概率为零。'
          '梯度各20,000次参考置换；L3注入 OR1/2/4/8 各3,000份数据，另执行每模型自校准。', '',
          '| 对比 | L3生成、L3检验：OR1拒绝率 | OR2功效 | OR4功效 | OR8功效 |',
          '|---|---:|---:|---:|---:|']
for contrast, part in chosen.groupby('contrast', sort=False):
    rates = part.set_index('odds_ratio').rejection_rate
    lines.append(f'| {contrast} | {rates.loc[1]:.3%} | {rates.loc[2]:.4f} | {rates.loc[4]:.4f} | {rates.loc[8]:.4f} |')
lines += ['', '以上梯度表为α=0.05。L3生成但用L1检验的OR1拒绝率属于交叉校准；不能代替L1自身生成的校准率。'
          '梯度与主图 Quintile 的分层条件一致，但置换种子不同，不可互换其MC p。', '',
          '## 验证与交付', '',
          '八项意图测试、五模型十项检验、16个joint阈值/层设置、十项梯度检验和90行梯度模拟均复核通过，跳过零项。'
          '由冻结效应独立重建真实研究数，核验完整旧特征、候选成员、零事件支持、解析期望、置换尾概率、BH族及输入输出哈希。', '',
          '本次五模型及梯度模拟表与前次已核验结果在1e-14容差内逐项一致。所有统计已在新目录重新执行；'
          '复用的仅是来自p50自身原生APA/GEX矩阵的已核验贡献者/协变量提取。', '',
          f'本地下载的{len(manifest["files"])}个汇总/诊断文件大小和SHA-256全部通过。'
          '不含基因/细胞标识、贡献者数据、矩阵或置换数组；完整配置、逐阶段命令、日志及数组留在服务器。', '',
          '- 主结果：five-null-results.tsv。', '- 扩展诊断摘要：extension-diagnostics.tsv。',
          '- 服务器原诊断：results/apadetect-p50/inference/{balance.tsv,calibration-power.tsv}。',
          '- 梯度原表：results/apadetect-p50/ladder/{cross_species_ladder.tsv,cross_species_power.tsv}。',
          '- 下载清单：results/HANDOFF-MANIFEST.json；本地核验：download-verification.json。', '',
          '服务器：213docker:' + manifest['server_directory'] + '/', '']
(root / 'REPORT.zh.md').write_text('\n'.join(lines))
(root / 'download-verification.json').write_text(json.dumps(dict(status='PASS',
    files=len(manifest['files']), bytes=sum(item['bytes'] for item in manifest['files']),
    missing=0, hash_errors=0, tests_skipped=0, five_null_results_match_prior=True,
    masked_ladder_results_match_prior=True, zero_events_outside_eligibility=True), indent=2) + '\n')
print(models[models.contrast.eq('PM')][['model', 'observed', 'analytic_expected', 'q_two_contrasts', 'q_four_tests']].to_string(index=False))
