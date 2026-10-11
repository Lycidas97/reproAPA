# 共享近端核心：事件分型与长异构体清除候选（结果前冻结）

本轮是新候选的探索性验证，不把已有因子阴性结果视为其正证据。固定前景为 P–M、|ΔDWUI|≥0.02 的 27 对同源共同近端基因，背景固定为原同一对比的 751 对可测基因。不修改原前景、阈值或测量设置。

## 事件定义

使用 supported50_p08 的原生 PAS 轴及 route_effects.primary P_M，每一个可评估 route 单独保留。driver 是 delta_usage 最大的 gain 和最小的 loss（并列按原生 PAS ID 字典序），记录其原生切点、Δusage、转录方向和 ΔDWUI；gain/loss 必须为正/负，否则 unknown。driver 不根据 motif 选择。

以冻结 Ensembl 116 的 transcript/exon/CDS 和 terminal-UTR 注释判断：同一转录本末端 UTR 覆盖两个 driver 为 strict_tandem_UTR；同一末端外显子但不满足 UTR 为 terminal_exon_nonUTR；两个 driver 属于不同末端外显子、没有任何共同外显子覆盖两个点且无 UTR 对应为 alternative_terminal_exon_compatible；有未被任何该基因外显子覆盖、但在转录本范围内的 driver 为 intronic_compatible；同时存在冲突结构解释或不同 route 分类为 mixed；其余为 unknown。分类表示静态结构兼容性，不证明实际转录本。每 study/route 的异质性均保留，不能将近端 gene-wide DWUI 改写为 UTR 缩短。

## 全区间序列资格与区间选择

主 motif 分析要求双方所有可评估 P_M study 的完整 primary PAS 轴均为 strict terminal UTR，且同一物种跨 study 的 compatible_transcripts 有非空交集。按坐标固定 union 轴最近端与最远端两个原生 PAS（不看 motif），两点之间完整 distal-only aUTR（正链 [prox,dist)，负链 [dist,prox)，反向互补）。对每 species/gene 选字典序第一 compatible transcript，仅用于核验整个区间确属 terminal UTR。另记录每 route 的 driver class，要求双方所有 driver 均 strict_tandem_UTR 才进入双物种检验。参考序列从冻结 GRCh38-2024-A / GRCm39-2024-A FASTA 读取并核对 SHA256；不做版本转换，不用 100 nt 窗口替代。完整区间含任何非 ACGT 或无法取得序列都不可估，unknown 不计成阴性。

## 候选与统计

ESCC 家族来自 Subramanyam 2011（DOI 10.1038/nbt.1862）及本地冻结 miRNA 注释：hsa-miR-302b-3p、hsa-miR-372-3p、mmu-miR-294-3p 的 nt2–8 共同 seed AAGUGCU。RNA 靶标侧 DNA 的 7mer-m8=AGCACTT、7mer-A1=GCACTTA、8mer=AGCACTTA。主结局为双方 aUTR 各含 ≥1 个 7mer-m8/7mer-A1；8mer 是严格敏感性。序列阳性仅表示潜在靶位点，不证明结合、RNA decay 或实际 miRNA 活性。

TIA1/U-rich 为预先列出的次要候选，仅使用 TTTT 重叠计数/有效 4-mer 起点数；不能将 U-rich 视为 TIA1 特异性结合，更不能推断 stress。其依据为 Zheng et al. 2018（DOI 10.1038/s41467-018-04730-7）。

双物种均合格背景上，9 个预设上尾统计为 human/mouse/both 的 ESCC ≥7 presence（3）、ESCC 8mer presence（3）、human/mouse 的 TTTT density 和两者最小 density（3）。gene 对为抽样单位。两个零模型均在资格背景中条件抽取与真实前景相同数量的 gene 对：coarse 按人鼠可评估 study pattern 精确分层，并加入两物种平均 log aUTR length、平均 GC 的各二分位；full 再加入两物种平均 M baseline DWUI、平均 PAS number 的各二分位。分位切点按完整合格背景，不按前景；同一 stratum 固定 foreground 数，包括不能交换的小层。100,000 次抽样，固定 seed 2026101007。所有 18 个 motif/statistic/model 检验统一 BH 校正。coarse/full 是同一数据的敏感性分析，不是独立复现。

报告零模型可交换前景数、singleton/全前景 strata、length/GC/测量/基线 SMD；模型对未直接分层的因子（GEX、APA coverage、capacity 等）的 max |SMD|>0.20 时明确残余混杂，不把其显著性升级为机制支持。合格前景不足 5、背景少于 20、或可交换前景不足 5 时取消正式富集声明，仅输出 descriptive 与离散置换值。不追加阈值、motif、事后扩大背景或重新选择 transcript 来改善结果。

## 验证与边界

分别保存逐 route 事件表、逐 gene 跨 study 结构/区间表、全序列/参考 hash、motif hit 坐标、9×2 统计及覆盖、不合格原因、输入/输出 hash。独立实现重算 motif（枚举 vs regex）、结构判定及序列对齐；真实场景测试须能发现负链/坐标 off-by-one、seed 正反互补错误、8mer 嵌套、unknown 被当作 0。机制最多 compatible：成熟 RNA 的比例变化仍不能区分 PAS 生成与长 RNA 清除，且 full-axis 资格子集可能与原 27 个核心不同。
