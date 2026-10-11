# 增殖对照与路线 P–M 均值：具体执行约定

本文件在新统计计算前写定。polyApipe contract-r2 是现稿主流程；APA_detect reject_original_point_terminal_p50 是替代流程。不混合两套 PAS 坐标、计数或候选。cycle archive 使用既有 `expression_inputs.npz`，必须与两套 measurement audit 的 SHA256 一致。

1. 分数使用既有 Cycle_S / Cycle_G2M marker 的 log1p(CP10K) 均值。先逐库、逐状态给出 metacell 等权均值、细胞数和 metacell 数；purity≥0.8、至少3个 metacell 且至少30细胞。P–M、P–S、M–S 只在同物理库同时具备两个状态时计算。每库只有一个观察均值，不进行细胞伪重复检验。路线汇总以物理库等权；至少2库才给 sign-flip p 与 bootstrap CI。缺状态须显式记录。两套流程共用冻结 GEX 因而不重复相同 cycle-score 检验。
2. 每套流程固定自己的 P–M WPP 0.02 exact-series recurrence 候选。polyApipe 原始 matched-measurable 需核对395小鼠、250人类，双参照222/128。另对 polyApipe 395/250候选在替代 PAS 流程作固定集合敏感性检查。
3. 沿用既有 measurement 的路线范围 S、G2M 三分位固定箱边界，生成3×3箱；相同分数导致边界折叠时沿用 numpy unique+searchsorted。每个基因原始真来源三联体内，三个状态均先应用原始基因深度≥20、purity≥0.8。只保留三个状态共有的箱，箱权重 ∝ min(nP,nM,nS)。每状态共有箱保留≥3 metacell、≥30细胞，并对PS/PM/MS各要求≥2个联合PAS（每PAS合计≥10分子且≥2 metacell检出）。三个状态使用相同箱权重，验证 PS=PM+MS。不可测单列，不能当作失去对齐。
4. 原始三联体的路线权重固定为既有 `triplet_weight`；匹配后仅对可测三联体重新归一化该权重，不按新效应调权。输出原始全部三联体、原始仅匹配可测三联体、匹配三个版本，以区分组成变化和周期匹配变化。路线对齐需PS、PM均沿固定复现方向≥0.02；反向对齐需二者均沿反向≥0.02，其他可测为失去对齐。全物种集合依原始必需路线全部成立；任何必需路线不可测则为不可测。
5. 对齐计数以基因为描述性分母。匹配前后推断以库为单位：在每库可测的原始候选中计算原始/匹配对齐比例的差，然后库等权 sign-flip。至少2库；共享来源S不构成新独立重复，因此该区间/检验仅为条件于冻结来源库的路线敏感性诊断。输出 `A2-cycle-alignment` 原始p，不因阴性或大量不可测调箱。
6. 第8项：主分析为所有原始同库P–M可测基因的库内平均DWUI，库等权路线均值。周期匹配采用同样路线固定3×3箱，两状态共同箱 min(nP,nM) 权重，遵守相同最小覆盖/联合PAS条件；同样提供raw-all、raw-common、cycle-matched，后两者固定为可测共同基因以避免组成混淆。基因不可作为路线重复。
7. sign-flip采用全部2^n组合（n≤20）；更大n用10000次Monte Carlo并加1修正。bootstrap对物理库有放回重采样10000次，种子20261008加SHA256键偏移；95%百分位CI。n=1仅给点值，CI/p为NA并标注。每项保留family/test_id/p字段供root统一BH，不自行挑显著结果。

新对照在初步结果之后设定。本分析衡量现有marker分数可比较范围内的效应保留，不能把匹配解释成排除全部增殖机制或因果验证。
