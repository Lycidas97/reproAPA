# H1：周期标记与 APA 的条件关联

固定于 2026-10-10，写于新 GEX/PAS 配对结果提取之前。这是既有结果后的预设补充分析，不是预注册。主数据 v2_primary、supported50_p08、基因 DWUI；不改上游数据、正文或既有分析。

## 已有证据和缺口

正文 Supplementary Methods 已固定 10 个 S 和 8 个 G2/M 标记，并已有库配对周期差异及周期共同支持诊断：P 常高于 M，严格 P/M 平衡不可估。故 H1 并非完全没测。原始 GEX 新恢复并不创造不存在的周期共同支持。本次补各路线 S/M/P 固定口径量级、及同一物理库/状态内部的周期—APA 关联。周期转录评分代表周期程序，不等于细胞分裂速率。

## 评分和数据

S：MCM5, PCNA, TYMS, FEN1, MCM2, MCM4, RRM1, UNG, GINS2, MCM6。
G2M：CDK1, TOP2A, NUSAP1, CDC20, CDC45, UBE2C, BIRC5, TPX2。
来源：当前正文 supplement/supplementary-methods/within-species-recurrence.tex 固定标记表。通过 GEX gene symbol（大小写归一）、必要时冻结 10x GTF 解析占位符，每标记必须唯一映射，缺失显式记录。不新增/删减结果驱动标记。两物种用同一 18 个符号及相同公式。

每 MC 以全部基因 GEX UMIs 为分母，均值 ln(1+CP10K) 为 S/G2M 分数；cycle 为两个分数的等权均值。背景校正敏感性：每标记用全路线 S/M/P 中均值 ln(1+CP10K) 最接近的 20 个非周期基因（只依赖全体表达，排除 18 个标记；距离并列按 gene ID 排序），减去对应控制均值。保留原始与背景评分，避免把 CPM 分母现象当周期效应。背景选择不用 APA、状态对比或关联方向。

路线：CIP、WOT、Xing、Liu、Guan 两 series、Panariello。状态 S/M/P，purity≥0.8。状态/库至少 3 MC 和 30 cells；P−M 仅同一物理库双方合格，库等权；另报时间等权。S 是实际观测 S，并明确文档化 source 标志，Liu/Pan 的 observed S 不冒充起始态。跨物种只描述各固定路线量级，无物种推断 p 值。

## 库内 APA 汇总

每个合格的 route/library/state 中，预先要求 nMC≥12。PAS 位点、primary full-axis position 与当前 polyapipe_effects.py 完全同口径；保留 supported50；gene full primary 轴必须同染色体同链且至少 2 个不同 cut，supported 位点至少 2 个。固定 panel 由该 context 每个 MC 都有≥20 gene PAS UMI 的基因构成，且 context 合计至少 2 PAS 各≥10 UMIs、各至少 2 MC 非零。要求 panel≥20 genes；所有 MC 均用这同一 panel 的等基因平均 DWUI，防止 MC 基因组成差异。panel 不按周期或 APA 方向选择。基因 GEX 唯一匹配率要求 100%；缺失不隐式丢弃。另报 panel 中位数，及同一 panel/MC 每基因下采样至恰好20 PAS UMIs 的 20 个固定种子抽样均值（不改变 panel）。记录所有不合格 context 的原因。

## 关联与控制

每库/状态单独：周期分数与 mean DWUI 的 Pearson 相关、及控制后的相关和原始单位斜率；负值表示高周期更近端。covariates：panel 平均 ln(1+gene CP10K)、panel 平均 ln(1+gene PAS UMI)、ln(1+全 GEX UMI/n_cells)、ln(1+全 PAS UMI/n_cells)、ln(n_cells)。OLS 含截距、按数值秩计算可用自由度，残差 df≥5 才估计；常量周期/APA、秩不足均不转为零。score background adjusted、APA median、20 UMI 等深为预设敏感性。

不把 MC 当生物重复。不为各 MC 拟合 p 值。库内相关仅测量条件关联；按物理库独立展示、同一路线等权汇总同库可用状态后再库等权。不重复计同库 S/M/P 为三个独立样本。至少 2 库报告 10,000 library-block bootstrap 95%范围，以及按可用 collection date 等权、leave-one-library-out 稳定性；单库不提供重复推断。多个库/日期仍可能共享一次诱导/细胞系，bootstrap 不是生物总体置信区间。本分析不提供跨物种或因果 p 值。依据每条路线方向、数量、敏感性一致性回答 H1；不以“未显著”排除增殖机制。

## 验证

冻结协议 hash；远程只读提取，原始整数 GEX 与 PAS 分母守恒、MC ID/库/细胞数一致；复核已有18标记分数（可取得相应表时）与 current full-axis position。意图测试固定 panel 防组成混杂、covariate 调整能去除深度伪相关、library 等权防 MC 数驱动、20 UMI 无放回守恒。输入、协议、代码、表 SHA256 清单。
