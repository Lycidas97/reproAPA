# A2/A8执行前冻结修订

2026-10-08。本文件在任何新周期权重、APA效应或周期p计算前固定，替换总协议中A2/A8的候选角色、范围和正式检验位置。A1/A3的科学方法和家族不改变；它们的结果不进入本模块的方案选择。修订根据估计目标和输入覆盖决定，避免把多个二元保留率、候选口径及敏感性重复当作主要推断。

## 范围与权重

PM两态：CIP、WOT、Xing/GSE118258、Liu/GSE150311、Panariello/GSE221739五个真实同库PM研究，两PAS。后二者仅需只读提取冻结PAS矩阵列，核验原计数、metadata及feature-index SHA，不重新调用PAS。没有同库PM观测的研究不外推。

PMS：仍限定原CIP、WOT、Xing记录来源S三路线。PM主要对照与PMS双参照使用不同、独立命名的条件目标，不能拼成三角。

单一主要调整为在purity≥0.8 metadata上拟合的非负连续权重；按每物理PM库或PMS目标/来源context拟合一次，与基因APA效应无关。每状态权重和1，两个marker加权均值在状态间相等，最小化Σ n_s||w_s−1/n_s||²。linprog先核验共同均值支持，SLSQP用解析梯度、ftol=1e−10、maxiter=300；均值等式残差≤1e−8。只可去掉数值冗余的约束行，不能按结果换估计量或求解器模型。

固定2×2中位数分箱为敏感性，metadata级共有箱按最小状态MC数赋共享权重。两种方法都在库级及每基因depth≥20限制/归一化之后验收：正权重>1e−8、每状态≥3 MC、≥30观察细胞、ESS≥3；两marker所有状态对间SMD≤0.1、加权方差比[0.5,2]、加权KS≤0.2。SMD用相应原始合格组的两状态等权池内SD；零方差规则、分位诊断及数值公差在模块协议中明确。PS/PM各≥2 joint PAS（至少10分子与2正权重MC检出），不新加MS门槛。旧PMS输入已在上游筛过MS，本轮不恢复该全背景。

mean平衡可行但ESS或分布验收失败须区别报告；通过只能称两marker达到预定条件，不称完全控制增殖。支持外不估计。2×2绝不自动替代连续主方法，无论哪个保留数量更多。

## 固定候选、贡献者及主要配对比较

PM候选：两PAS都固定主流程小鼠420/人类268原PM复现基因及方向；替代自身215/194仅描述。PMS正式比较固定主流程395/250匹配可测候选及222/128原双参照标签；替代自身203/174及112/85另作描述，不新增正式筛选检验。

提供四种明确版本：

1. raw_all：全部原合格贡献者/MC的原始效果。
2. raw_context_common：只保留新方法可测的gene/context，但在该context内全部gene-depth合格MC仍均匀平均。
3. raw_support_common：同样gene/context，进一步限定调整中正权重MC，仍等权平均。
4. adjusted：同样gene/context及正权重MC，使用冻结调整权重。

主要配对变化为adjusted减raw_context_common，包含全部within-context周期调整。raw_support_common用于分解支持修剪与再加权；不把它与raw_context_common混称，不只估再加权而遗漏周期支持修剪。原triplet相对权重保留，合格贡献者内归一化；重复source不当新库。missing不作失去对齐或效应零。

## 正式家族

family=`A2-cycle-alignment-redesign`仅连续主方法16个计划位置：两PAS×(5 PM+3 PMS)。PM统计量为固定候选sign×(adjusted_PM−raw_context_common_PM)的每库基因均值；PMS统计量为固定候选min(sign×PS,sign×PM)的adjusted与raw_context_common差，先按原source相对权重合并到基因/目标库，再每库基因均值。统计量是连续效应/双参照弱项边距变化，家族名沿用任务名称，不把它叫二元对齐率检验。

family=`A8-cycle-PM-redesign`仅连续主方法10个计划位置：两PAS×五路线，检验全部原可测PM基因的adjusted库等权路线均值相对零。raw各版本及配对变化完整描述，不重复正式p。

二元retained/flipped/lost/unmeasurable、2×2所有效应及替代自身候选全部描述性，输出完整，不能择有利设置宣称正式验证。两个家族分别对全部有效p统一BH q<0.05；n<2为NA，不按新结果拆family。

推断单位为物理库；n≤20双侧sign-flip精确枚举，更多库10000次plus-one；不使用细胞/基因扩充n或换t检验救小n。库bootstrap 10000次区间作有经验库分布限制的描述，报告库值及可达p。CIP/Xing未确认独立生物重复与WOT日期异质性、来源S重用单列，研究内证据不冒充跨研究复制。

## 执行顺序与验收

最终`cycle/PROTOCOL.md`必须与本修订一致，保存其及本文件/原总协议/代码SHA与UTC执行前回执，root明确授权后再收取新增PAS和运行新效应。保留draft与修订记录。输入前后SHA、原效应重放、weights/ESS/平衡/失败原因、PMS加法关系及所有分母完整核验；独立验证主要权重的凸约束与优化目标，小样本检查weighted KS及sign-flip。不自动改正文。
