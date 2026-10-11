# 共同近端集合的 MET / QKI 外部扰动验证

本协议在只检查外部表结构、标识符、有限数值资格和覆盖数后固定；尚未计算 27 个基因的外部效应符号、靶标命中或富集统计。它是发现共同集合之后的探索性机制检验，不是预注册实验。

## 固定输入和假说

- 前景固定为 `target_members.tsv` 的 P_M 共同近端 27 对；背景固定为 `ortholog_roster.tsv.gz` 的 P_M 751 对。阈值 0.05 的 7 对仅作描述性敏感性，不替换主前景。
- 主外部资料为 Neumann 等（DOI `10.1080/15476286.2023.2294222`，PMCID `PMC10732628`）的 Supplementary Table 8：mesHMLE 中 QKI knockdown 相对 mesHMLE 的全测试基因 APA 表。辅助资料为同文 Table 7：mesHMLE 相对 HMLE 的 EMT APA 表。
- 检验预测：共享近端基因对 QKI knockdown 的近端响应富集。辅助预测：同一集合在 EMT 中远端响应富集。QKI knockdown 远端、EMT 近端是方向对照。基因重叠不能证明重编程中 QKI 的因果作用，也不能证明 M→P 就是这些外部细胞的 MET。

## 标识符与覆盖

外部 `name` 列为 Ensembl gene ID，`gene` 列为符号。仅精确 Ensembl ID 对接，不用符号或挑选旧名补命中。外部 `name` 若重复则整个 ID 排除；`effect` 和 `fdr_zero` 必须为有限值。当前预检：S7 / S8 各可覆盖 605 / 751 背景和 25 / 27 主前景，0.05 集合为 7 / 7。缺失保留 NA，不记为不响应。保存全部 751 对的资格和缺失原因。

主集合保留所有 APA 事件。另作两个预定结构敏感性集合：两个物种每个现有 route 的 driver pair 均为 `strict_tandem_UTR`；两个物种每个 route 的整个 primary PAS axis 均为 `full_primary_axis_tandem=True`。二者在外部交集前分别有 315 对 / 17 核心及 67 对 / 8 核心。标签是静态转录本结构兼容性，不等同于实际表达异构体。

## 外部响应和统计

外部响应由作者 `fdr_zero < 0.05` 且 `effect < 0`（近端）或 `effect > 0`（远端）定义；不把作者 effect 当作 DWUI，也不套用本项目 0.02 / 0.05 效应阈值。`effect=0` 不计方向。作者样本模型产生的 FDR 不重算。核对全表可否重现论文所报告的 S8 94 近端 / 25 远端、S7 27 近端 / 52 远端；若不一致，先查清阈值 / 表内容再报告，不根据 27 集合调整规则。

四个预定方向终点 × 三个结构面板（全部、driver strict、axis strict）× 三个零模型（未匹配、COARSE、FULL）共 36 格为一个 BH 家族；不可估格保守以 p=1 进入家族，但机器表原值 NA，注明原因。未匹配用单侧条件超几何尾概率，等价 Fisher 富集检验。COARSE / FULL 沿用冻结的双物种 assay covariates：COARSE 为两物种 `measured_pattern, proximal_capacity, distal_capacity, PAS_number_bin`；FULL 再加入两物种 `mean_GEX_bin, APA_coverage_bin, baseline_bin`。不对结果重新分箱或换模型。

每个匹配 stratum 固定 N、K（该外部响应数）和 n（核心数），以超几何分布卷积得到核心命中数零分布。显式报告随机化方差、point-mass、可变前景数，若总零分布方差为零则不给伪 p 值。匹配控制的是既有 assay covariates，不等于控制外部细胞背景。统计单位为基因；这是对外部扰动靶标的条件集合检验，不把 25 个基因视为 25 个独立生物学实验。

对全测量交集、核心及 0.05 核心另报告效应中位数和正/负比例、作者 FDR 阈值 0.01 / 0.05 / 0.10 下的响应计数；它们为描述性敏感性，不提供额外未经校正的显著性。列出所有核心，包含不响应和不可评估项。

## 表达和验证

QKI / Qki 的原 GEX 表达以原 NPZ 中每个 metacell 的 log2(CPM+1) 独立重算，同库 P/M 各状态至少 3 个 metacell 和 30 个细胞，匹配 library 后等权汇总；附时间点等权结果。四条现有路线并非 QKI RNA 普遍下降（WOT 上升），mRNA 也不能替代 QKI-5 蛋白或活性。其余路线原 RBP archive 未收录，不能记为零。

保留输入 SHA256、预检、冻结时间、来源说明和输出表。另实现确定性核验：独立读取原 Excel 的标识符/作者值、重新计数核心/背景/表内公开方向；用 SciPy 超几何 PMF / 手工卷积核验匹配分布和 SciPy Fisher 核验未匹配尾概率；从原 NPZ 核验表达聚合。所有未完成或不可估分析明示，不用阴性结果排除机制。
