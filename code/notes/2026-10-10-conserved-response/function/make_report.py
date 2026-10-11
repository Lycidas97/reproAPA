"""Render the sealed tables without revising analysis definitions."""
from pathlib import Path
import json
import pandas as pd

HERE = Path(__file__).resolve().parent
results = pd.read_csv(HERE/'category_results.tsv', sep='\t')
primary = results[results.primary].sort_values(['set_name', 'category'])
source = pd.read_csv(HERE/'comparison_membership_features.tsv.gz', sep='\t', low_memory=False)
membership = pd.read_csv(HERE/'category_membership.tsv.gz', sep='\t')
annotations = pd.read_csv(HERE/'annotation_gene_membership.tsv.gz', sep='\t')
case_rows = []
for set_name in ['main48', 'hN79']:
    genes = source[source.set_name.eq(set_name)&source.threshold.eq(.02)&source.case]
    for row in genes.itertuples():
        record = dict(set_name=set_name, human_gene_id=row.human_gene_id, mouse_gene_id=row.mouse_gene_id,
            human_gene_symbol=row.human_gene_symbol, mouse_gene_symbol=row.mouse_gene_symbol)
        for category in ['GO:0003723', 'GO:0006397', 'GO:0006417']:
            z = membership[membership.set_name.eq(set_name)&membership.threshold.eq(.02)&membership.annotation.eq('human')&
                membership.model.eq('coarse')&membership.category.eq(category)&membership.human_gene_id.eq(row.human_gene_id)]
            assert len(z) == 1
            record[category+'_known'] = bool(z.annotation_known.iloc[0]); record[category+'_member'] = bool(z.category_member.iloc[0])
        case_rows.append(record)
case_table = pd.DataFrame(case_rows)
case_table.to_csv(HERE/'primary_case_categories.tsv', sep='\t', index=False)
table_lines = ['| 集合 | 类别 | 注释已知目标 / 背景 | 命中 / 条件期望 | 倍数 | p | 主 BH6 q |',
               '|---|---|---:|---:|---:|---:|---:|']
category_names = {'GO:0003723': 'RNA 结合', 'GO:0006397': 'mRNA 加工', 'GO:0006417': '翻译调控'}
for row in primary.itertuples():
    table_lines.append(f'| {row.set_name} | {category_names[row.category]} | {row.cases} / {row.N} | {row.hits} / {row.expected:.2f} | {row.matched_fold:.2f} | {row.p:.3f} | {row.q_primary6:.3f} |')
examples = []
for set_name in ['main48', 'hN79']:
    z = case_table[case_table.set_name.eq(set_name)]
    for cat, label in category_names.items():
        names = z[z[cat+'_known']&z[cat+'_member']].human_gene_symbol.tolist()
        examples.append(f'- {set_name}，{label}：'+', '.join(names)+'。')
review = HERE.parent/'review/final_validation.json'
qa_text = '独立数值审计进行中；下面为已封口计算结果，不能据此声称独立验收完成。'
if review.exists():
    payload = json.loads(review.read_text())
    qa_text = '独立审计结果文件已返回，状态：`'+str(payload.get('status', '见审计文件'))+'`；详见 `../review/final_validation.json`。'
text = '''# 跨物种同向基因：功能类别条件比较

结论：未支持“跨物种共享部分集中于转录后调控因子”的集合层面富集。六项主检验全部未通过 BH6；预设全部敏感性及 BP/MF 探索扫描亦无校正后显著项。名单中的个别 RNA 因子可以作为候选生物学实例，不能据此把共享集合改写成富集的 RNA 调控模块或 RNA 加工机器的反馈。

## 集合和实际比较对象

本轮协议在功能命中数前固定，但因已经看过基因名字后提出假说，属于事后探索，不称预注册。主 case 是四个 same-contrast 同向格并集，P–S 双近端7、双远端14，P–M 双近端27、双远端3，去重48。原用户名单的前三格并集45另作敏感性。主 case/control 库为48+539=587；control 是只一个物种在其双方实际可测 contrasts 内复现，另一物种所有这些可测 contrasts 均未复现。反向或不同 contrast 双边复现者不是单物种对照，未测也不是无响应。

hN 次集合是79个唯一同源对，control为单侧近端而另一侧不满足任一方向响应，共79+211=290。旧125行等于79加内部双参照46，没有125个不同基因；dual46为敏感性。hN侧是两个培养条件平均响应，鼠侧是研究级复现，不能把hN称人跨研究复现或同contrast事件保守。所有集合每ortholog一票；0.05完整重复。

## 六项主检验

统一用人GO注释，coarse匹配双物种表达、APA覆盖及PAS数的最小百分位箱，并固定可测contrast mask、响应方向类别。每类别只把相应MF/BP namespace已知注释者纳入；未知单列流失，不当非成员。GO资格过滤后保留原箱，固定每层case数，以精确分层超几何卷积求上尾P。

'''+ '\n'.join(table_lines)+'''

主48在RNA类别没有未知目标，BP类别缺2个目标；hN79在RNA类别缺1个目标，BP类别缺8个目标。因此78/71是类别已知资格分母，不是改变79的父集合。三个类别可以重叠，命中总数不能相加当不同基因数。

条件零分布的95%范围见 `category_results.tsv`，它是给定基因/层的参考范围，不是物种总体的生物学置信区间。上尾结果不构成“功能耗竭”检验；主48命中低于期望不能直接写成显著耗竭。

## 敏感性与 GO 探索

144项完整预设配置＝四集合×两阈值×coarse/fine×人/鼠/AND注释×三个类别；六项主检验BH6，另外138项作为独立敏感性族BH138，另报joint144。全部144项实际可检验，未有空目标/无方差配置被跳过。敏感性最低名义p为0.0804（hN、mRNA加工、AND、fine、0.02）；BH138全部q=1，joint144也全部q=1。

主48在0.02的coarse/fine和人/鼠/AND注释下均没有富集：RNA类别匹配倍数0.48–0.80，mRNA加工0.53–0.77，翻译调控0.63–0.70。hN有些mRNA加工点估计约1.3–1.4倍，最强名义检验仍p=0.0804，不构成显著模块身份。list45、dual46及0.05没有改变这一判断，不用其中配置替代主口径。

BP+MF全扫描只针对main48/hN79、人注释、coarse、0.02，按实际namespace-qualified背景term size10–500确定family，不用目标命中数筛term。main48共736 terms，hN79共442，合计1178；各集合及联合扫描q全部为1。名义最低p项包括hN的“cellular component organization or biogenesis”（p=0.00376），以及main48的“regulation of intracellular transport”（p=0.0271），均不能写成校正后通路富集。旧C4阴性仍保留；本轮比较的是共享部分相对单物种复现部分，是不同的问题。

## 匹配的实际能力和残差

未按GO过滤时，main48 coarse为87层，46/48 case可交换（95.8%）；main fine加两个物种各contrast精确可测study/series签名，43/48可交换（89.6%）。hN coarse为24层、79/79可交换；fine再加三个trait的双物种最大百分位，56/79可交换（70.9%）。注释过滤后的实际移动性均按类别输出，不能只引未过滤的最大值。

`covariate_balance.tsv.gz`使用冻结trait的原值，没有追加变换。main两物种的mean_GEX来源为端点GEX_CPM平均再跨贡献取中位数，是CPM proxy；hN human来源为各实际条件/端点的mean_log2CPM1等权合成，hN mouse仍来自复现表的CPM proxy。每species/context单独做百分位排序，未将两物种这些数值直接相加，也未把混合proxy当成相同表达单位；这里只匹配排序代理，不能声称表达量逐数值可比。APA覆盖是每qualified MC平均PAS分子数。初版把这些GEX代理统称为log2(CPM+1)是报告说明错误，本轮据冻结downstream_bundle.py:56及proximal-module-specificity/analyze.py:229纠正；协变量值、rank、层、balance数值和检验全部保持。standardized差的分母是整个namespace-qualified背景有限值的样本SD（ddof=1）。条件期望case均值用每gene权重n_case层/N层计算，表中残差是实际case均值减去这个条件期望再除背景SD；它不是case/control传统合并组内SD。

human注释、coarse、0.02下，三个匹配trait在hN的残差绝对值≤0.087；main48仍存在最高约0.166 SD的残差（小鼠覆盖偏低）。未直接匹配的方向容量残差在main达到约0.228 SD、hN约0.210 SD。fine未使main的三个trait全部平衡（最高约0.163 SD），也不能称完整六维逐变量匹配。功效、未知注释、基因功能注释完整度差异和未控制特征仍限制阴性解释；本结果是不支持当前富集假说，不是证明共享基因完全没有功能组织。

## 目标中的 GO 类别实例

下列是冻结类别自动分类的全体命中，不能依据文献印象手工给预设类别补成员；类别名称和实验功能不是同义词：

'''+ '\n'.join(examples)+'''

## 可写与不可写

可以写：两物种净方向分歧下仍有部分同源基因同向响应；其中包含转录后调控因子。事后提出的三个功能类别在测量匹配、单物种复现对照下未显示富集。具体候选因子按其真实PAS资格和外部实验讨论，基因身份、PAS事件对应和功能机制分开。

不能写：共享集合集中于RNA加工机器；这些RNA因子的自身APA构成跨物种功能反馈；hN共享模块具有naive特异性；GO未富集证明没有生物学共性；阴性结果排除了任何周期/CPA/UGUA机制。也不能把gene-wide DWUI同向自动称作同一个同源PAS事件受到选择。

## 验证及交付

'''+qa_text+'''

12项意图测试通过，0跳过。它们包含并集一票、真实可测性、异向/跨contrast双复现不是单物种控制、ties、未知namespace、is_a/part_of不走regulates、精确卷积全枚举、NE占family p=1、BH与全尾概率边界。初跑在BH范围断言捕获全尾求和1+机器误差；一行min(1,tail)修复，不改PMF或集合/层。初版脚本和失败日志保留在before-roundoff与initial_roundoff_failure.log，记录见numerical_repair.json。

输入/输出版本和协议哈希见membership_manifest.json、functional_manifest.json。统计结果为category_results.tsv；每gene未知/命中/层为category_membership.tsv.gz；可复算层N/K/n/hits为category_stratum_margins.tsv.gz；完整精确PMF为category_exact_pmf.tsv.gz；注释direct/ancestor证据与根后代为annotation_gene_membership.tsv.gz及root_descendants.tsv；探索全family为exploratory_GO.tsv.gz；原始比较库为comparison_membership_features.tsv.gz。不改正文或冻结旧分析。
'''
(HERE/'REPORT.zh.md').write_text(text)
