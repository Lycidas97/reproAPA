# A2/A8 连续周期条件对照：执行协议

2026-10-08。看到初步结果后的方法修订，非预注册验证。`../CYCLE-FREEZE-ADDENDUM.md`取代总协议的A2/A8范围、比较版本与正式检验行；本文件落实该修订。旧输出、候选、PAS、marker和正文保持只读。两PAS使用同一算法、验收及固定主候选。

## 权重与验收

冻结分数为既有18个marker的log1p(CP10K)均值（10个S、8个G2M）。先在purity≥0.8 metadata的每个真实同库PM context，或原P/M目标库加原来源S库的PMS context上拟合，与APA效应无关。PM和PMS条件目标分别命名，不能拼成三角。

连续主方法：每状态非负metacell权重和为1，Cycle_S和Cycle_G2M两加权均值在状态间相等，最小化`Σ_s n_s||w_s−1/n_s||²`。每状态观察二维cycle凸包包含共同目标，不外推。两个marker先按context全部合格metadata的均值和SD标准化以改善数值条件（常数marker不加均值约束）；不改变共同均值目标。`linprog(method='highs')`核对可行性并提供初值；`SLSQP`用解析梯度、ftol=1e−10、maxiter=300。仅移除完全数值冗余的线性约束行（QR pivot，rank tolerance=1e−10），全原始约束复查，均值残差≤1e−8。LP不可行、数值失败、拟合可行但验收失败分别报告；不换模型救可测率。

权重≤1e−8视为数值零，归一化并再次验收。每状态至少3个正权重MC、至少30个观察细胞、ESS=1/Σw²≥3。所有状态对、两marker均要求|SMD|≤0.1、加权population variance比[0.5,2]、weighted KS≤0.2；验收比较容差1e−10。SMD分母为对应**原始合格两组**未加权population variance的等权平均之平方根（库级用purity合格组，基因级用depth合格组）。若双方零方差且均值相同（绝对差≤1e−12），SMD=0、variance ratio=1；原始池SD为零但均值不同时SMD无限。一方加权variance≤1e−14、另一方更大时方差比无限并失败；双方常数且均值不同同样失败。另输出固定q10/q50/q90、ESS/N、最大权重、marker均值、方差与库级尺度SMD诊断，不事后加门槛。

每基因depth≥20限制冻结库权重，归一化后重新验收全部连续门槛、ESS、≥3MC/30cells，不逐基因重拟合。PS和PM各要求≥2 joint PAS（至少10分子、至少2个正权重MC检出）；不新加MS门槛。PMS输入仍来自旧weighted_triplets，上游已筛MS joint资格，不能恢复已被旧门槛排除的全背景。通过验收仅称两个marker边际达到预定条件，不能称完全控制增殖。

固定2×2敏感性：每路线purity合格metadata两个marker各自中位数为边界（等于中位数归高箱），共有箱目标比例为各状态箱MC数的最小值，箱内均匀分摊。metadata级拟合后逐基因限制并使用与连续主方法相同的连续验收。它不自动替代主方法，也不因为分箱相同而视为控制成功。

## 范围、候选和四种版本

PM：CIP/GSE103221、WOT/GSE115943、Xing/GSE118258、Liu/GSE150311、Panariello/GSE221739五个原真实同库研究。PMS：只用CIP、WOT、Xing原记录来源S三路线；来源S可跨库，不伪称四状态同库。来源重用不计新库，所有source在合格贡献者内保留旧相对triplet权重。

PM正式候选两PAS均固定主流程小鼠420/人类268原PM复现基因及方向；替代自身215/194只描述。PMS正式两PAS均固定主流程395/250 matched-measurable候选及原222/128双参照标签，替代自身203/174与112/85只描述。原0.02双参照对齐阈值固定。

四种版本明确区分：

1. `raw_all`：全部原合格context和其全部gene-depth合格MC的原效果。
2. `raw_context_common`：仅新方法可测gene/context，保留其中全部gene-depth合格MC均匀均值。
3. `raw_support_common`：同样gene/context，进一步限定调整正权重MC并均匀平均。
4. `adjusted`：同样gene/context及正权重MC，使用冻结连续或2×2权重。

主要配对变化为adjusted减raw_context_common，包含支持修剪与再加权；raw_support_common用于分解。missing不是效应零或失去对齐。PMS在每三联体保持同一P权重，PS=PM+MS；按原相对source权重先合并到基因/目标库。每库基因等权，路线物理库等权。

## 正式推断和描述

`A2-cycle-alignment-redesign`固定16位置，两PAS×(五PM＋三PMS)，仅连续主方法。PM为固定候选`sign×(adjusted_PM−raw_context_common_PM)`每库基因均值；PMS为`min(sign×PS,sign×PM)`两版本的变化，先合并source到基因/目标库后再每库基因均值。名称沿用任务，不称二元对齐率检验。

`A8-cycle-PM-redesign`固定10位置，两PAS×五PM路线，仅连续主方法，检验全部原可测PM基因adjusted库等权路线均值相对零。各raw版本和配对变化完整描述，不重复正式p。固定候选及替代own的retained/flipped/lost/unmeasurable、2×2完整描述，不另选p。

推断单位物理库，双侧sign-flip：n≤20精确枚举，更多库10000次固定SHA派生种子、plus-one；n<2为NA。库bootstrap10000次2.5/97.5百分位区间，n<2为NA。报告n、库值、可达最小p=2/2^n（非零一般位置）及研究内重复/日期/固定来源S边界，不用基因或细胞扩充n，不用t检验救小n。root各家族统一BH q<0.05，NA不当p=1。

## 输入、执行与验证

已有只读PAS缓存覆盖三PMS研究全部primary PAS轴；只读追加两PAS×Liu/Panariello的相同轴列和selected metadata，不重调PAS。沿用原两caller独立feature-index与原矩阵SHA，所有原效果逐条重放（误差<1e−10）。全7研究原周期分数表引用旧结果，本轮无需重复该family。

执行前UTC回执保存本协议、总协议、冻结修订及执行代码SHA；root授权后才提取追加PAS或计算APA。输入前后hash、非负/和1/均值约束、优化目标、weights、ESS、平衡/失败、PS=PM+MS、全部分母和16＋10位置必须验证。独立小例子验证weighted KS及sign-flip，独立重算拟合凸约束/优化一阶条件和主要库统计量。旧文件不改，成功不以显著或保留数量定义。
