"""Render the descriptive classification report from the verified tables."""
from pathlib import Path
import json
import pandas as pd

HERE = Path(__file__).resolve().parent


def read(name):
    return pd.read_csv(HERE / name, sep='\t', float_precision='round_trip')


def table(frame):
    columns = list(frame.columns)
    lines = ['| ' + ' | '.join(columns) + ' |', '| ' + ' | '.join(['---'] * len(columns)) + ' |']
    lines += ['| ' + ' | '.join(str(row[col]) for col in columns) + ' |' for row in frame.to_dict('records')]
    return '\n'.join(lines)


def main():
    pairs = read('event_class_pairs.tsv')
    species = read('species_class_counts.tsv')
    margins = read('paired_stratum_direction_marginals.tsv')
    cores = read('shared_proximal_core_all_thresholds.tsv')
    coverage = read('both_ALE_terminal_UTR_coverage_pairs.tsv')
    validation = json.loads((HERE / 'validation.json').read_text())
    rows = []
    for threshold in [.02, .05]:
        part = margins[margins.scope.eq('all') & margins.threshold.eq(threshold)]
        for class_ in ['both_tandem', 'both_ALE', 'other/mixed']:
            class_margins = part[part.pair_class.eq(class_)]
            n = int(class_margins.class_background.iloc[0])
            def count(sp, direction):
                return int(class_margins[class_margins.species.eq(sp) & class_margins.direction.eq(direction)].n_genes.iloc[0])
            shared = len(cores[cores.threshold.eq(threshold) & cores.pair_class.eq(class_)])
            rows.append({'阈值': threshold, '双方结构': class_, '背景对数': n,
                         '人近端/远端/未复现': '/'.join(str(count('human', label)) for label in ['proximal', 'distal', 'nonrecurrent']),
                         '鼠近端/远端/未复现': '/'.join(str(count('mouse', label)) for label in ['proximal', 'distal', 'nonrecurrent']),
                         '共同近端': shared})
    scan_rows = []
    for threshold, group in cores.groupby('threshold'):
        counts = group.pair_class.value_counts()
        scan_rows.append({'阈值': format(threshold, '.3f'), '全部共同近端': len(group),
                          '双方串联UTR': int(counts.get('both_tandem', 0)),
                          '双方ALE兼容': int(counts.get('both_ALE', 0)),
                          '其余/混合': int(counts.get('other/mixed', 0)),
                          '双方全轴串联UTR严格资格': int(group.strict_full_axis_sensitivity.sum())})
    effect_rows = []
    for (threshold, class_), group in cores[cores.threshold.isin([.02, .05])].groupby(['threshold', 'pair_class']):
        effect_rows.append({'阈值': threshold, '双方结构': class_, '共同近端对数': len(group),
                           '人gene等study均值的中位数': f'{group.human_delta_DWUI_equal_study_mean.median():.6f}',
                           '鼠gene等study均值的中位数': f'{group.mouse_delta_DWUI_equal_study_mean.median():.6f}'})
    report = '''# P–M共同近端响应的末端结构分层：描述统计

2026-10-10。主分类为双方所有可评估route的最大gain/minimum-loss位点对都兼容同一串联UTR、双方都兼容不同末端外显子（ALE），或其余/混合结构。本报告只描述基因方向和注释兼容性；类内条件匹配检验由另一个分支完成，不在这里作生物重复推断。

**751对背景中，315对双方串联UTR、52对双方ALE兼容、384对其他/混合。0.02下的27对共同近端核心由17/4/6对构成；0.05下的7对为6/1/0对。** 不能把27对全部称为3′UTR缩短事件。

## 范围与对账

固定背景为751对P_M一对一直系同源基因，读取原九阈值方向标签。3467条route事件与冻结study_events的所有实际study/series_pattern完全一致（人1965，鼠1502）；没有因方向不同删除route。各study内均值/最小/最大DWUI与冻结表最大绝对误差为2.0e−16。分类固定一次，不随阈值改变。浮点阈值使用round-trip读取，不把0.015/0.03/0.045拆成不同键。

结构标签来源于响应所确定的driver pair，因此是效应依赖的结构兼容性分层，不能称为与响应独立的结构变量。兼容注释不证明实际剪接或完整转录本；gene-wide轴近端/远端不等于核苷酸长度或相同PAS的跨物种保守。

## 每个物种的分类分母

'''
    report += table(species.rename(columns={'species': '物种', 'event_class': '结构分类', 'n_genes': '基因数'}))
    report += '''

单物种分母包含另一物种不同结构的同源对，因此人421个串联UTR、鼠432个串联UTR不等于双方共同315对。完整物种结构交叉矩阵见[species_class_matrix.tsv](species_class_matrix.tsv)。

## 主阈值及敏感性阈值

下面的边际限定在“双方相同结构”的共同背景内；未复现包括阈值以下、支持研究不足以及方向不能复现的基因。不是“无APA”。

'''
    report += table(pd.DataFrame(rows))
    report += '''

全部单物种结构边际见[recurrence_direction_marginals.tsv](recurrence_direction_marginals.tsv)，共同配对分层边际见[paired_stratum_direction_marginals.tsv](paired_stratum_direction_marginals.tsv)，每层九阈值3×3方向表见[recurrence_cross_direction.tsv](recurrence_cross_direction.tsv)。原复现规则仍允许第三项人研究相反，只要至少两项研究支持；每支持研究内部的所有可评估系列均须达到同向阈值。

## 九阈值共同近端核心构成

'''
    report += table(pd.DataFrame(scan_rows))
    report += '''

双方全轴串联UTR严格敏感性背景为67对（both_tandem与双方all_full_axis_tandem的交集），主阈值核心8对、0.05核心3对。完整轴兼容资格不要求跨study共用一个transcript ID，也不能当已确认的同一异构体。

## 核心成员与效应描述

'''
    for threshold in [.02, .05]:
        report += f'**阈值{threshold:.2f}**\n\n'
        for class_, group in cores[cores.threshold.eq(threshold)].groupby('pair_class'):
            report += f'- {class_}：' + '、'.join(group.human_gene_symbol) + '。\n'
        report += '\n'
    report += '各基因先平均study内route，再在study间等权平均；下表是共同近端集合内该gene描述值的中位数。不是物种间的效应显著性比较，也没有把不同PAS轴转成相同物理长度。\n\n'
    report += table(pd.DataFrame(effect_rows))
    report += '''

完整基因名单和逐gene效应分别见[genes_both_tandem.tsv](genes_both_tandem.tsv)、[genes_both_ALE.tsv](genes_both_ALE.tsv)、[genes_other_mixed.tsv](genes_other_mixed.tsv)与[core_002_005_members.tsv](core_002_005_members.tsv)。

## “UTR切换”的额外语义资格核对

ALE原分类只要求driver两个端点兼容不同末端外显子，不要求每端均在UTR。作为事后语义核对，我们另查每条route的两个retained bases是否分别位于某注释terminal UTR：52对both_ALE中15对双方所有route均满足，37对至少一端/route不满足。这个额外资格不改变主分类和检验家族。

0.02核心PCMT1、KLC1、POLR1D、PDLIM7四对全部通过双端UTR覆盖；0.05的ALE核心是PDLIM7。**每端被某个UTR注释覆盖，不保证它们属于同一CDS版本。CDS未核定，ALE仍可能改变编码序列，不能直接称为UTR-only切换。** 完整额外资格在[both_ALE_terminal_UTR_coverage_pairs.tsv](both_ALE_terminal_UTR_coverage_pairs.tsv)，坐标规则和输入SHA在[ALE_UTR_coverage_contract.json](ALE_UTR_coverage_contract.json)。

## Driver方向与全轴DWUI的关系

243/3467条route的最大gain/minimum-loss位点对方向与全轴DWUI符号相反：人174/1965（147个基因至少一条），鼠69/1502（60个基因至少一条）。所有route均有有效gain/loss且非零DWUI。

0.02核心中，人2/68条route反号，鼠0/54：ZC3H15和RAB10均在Panariello/GSE221739，其最大gain位点偏近端，但净DWUI分别为+0.000502和+0.014505，均未达到该研究0.02远端阈值。这两条保留在分类中，说明最大driver对不是整个多PAS轴净响应的替代指标，不重新定义复现方向。0.05核心人0/18、鼠0/14反号。详情见[route_sign_audit.tsv](route_sign_audit.tsv)和[driver_sign_summary.tsv](driver_sign_summary.tsv)。

## 验证与交付边界

'''
    report += f'独立CSV/标准库复核{validation["checks"]}项通过，0跳过：完整贡献范围、固定分类、九阈值study支持与方向、全部边际/交叉计数、核心成员、原生链方向、ALE端点覆盖和输入SHA均对账。[validation.json](validation.json)保存结果。主分类表SHA固定为`5ba39cd3a5439557c7c3edf604ec6c9ca03c88acb0249996f4748572c17842c3`，新统计不会重写该冻结表。\n\n'
    report += '本报告没有计算新的p值，没有把未复现记成无变化，没有强行二分其他/混合结构。旧封存分析和稿件均未修改。\n'
    (HERE / 'REPORT.zh.md').write_text(report)


if __name__ == '__main__':
    main()
