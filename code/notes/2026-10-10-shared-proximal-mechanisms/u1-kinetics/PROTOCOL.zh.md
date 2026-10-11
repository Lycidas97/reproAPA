# U1 可用量与转录动力学候选：探索性验证协议

2026-10-10。这是发现 27 个共同 P−M 近端基因之后固定的探索性协议；不是独立预注册。冻结前只查看文献、公开文件类型、表头、静态资格和样本量，未查看新的 U1 覆盖差异或共同集合的结构统计。

## 输入与问题

固定 P_M one-to-one ortholog 背景 751 对，foreground_002=27；foreground_005=7 仅作描述性敏感性。事件分型使用同目录父层 event_inventory.tsv.gz，逐研究 gain/loss PAS 的静态兼容性不是实际长读长异构体证据。背景之外的基因、未检出基因均不作为无响应。

假说是某些基因对有效 U1/新生转录比例或近端位点附近停顿具有保守易感性。poly(A) GEX 无法可靠量化 U1 snRNA，SNRPA/SNRNP70 表达也不等于 U1 活性；PAS 结构相关不能证明转录速率改变。

## A. 公开 U1 扰动：固定探针的 LECDS 重建

Oh et al. 2020（doi:10.1038/s41467-019-13993-7，GSE140543）公开六个 hg38 normalized bigWig：AMO Ctrl/12.5/62.5 pmol、OE Ctrl/1/1.5 µg。文章完整 LECDS/DaPars tested gene table 未公开；PDF 仅有精选癌症名单，Source XLSX 仅表型及精选表达。因此不能把附件未列出的基因当作未响应，不能复刻作者 adjusted P≤0.01 分类。

只读 HTTPS genomic range 查询六个轨道，使用冻结 Ensembl116 hg38 exon/CDS 结构，记录精确区间和查询结果。先按静态结构选择 transcript，后读取信号：人类 gene 的所有贡献研究 driver 均 strict_tandem_UTR，所有研究共有兼容 transcript；候选必须同时具有最后外显子 CDS 和完整末端 3′UTR。选择总 CDS 长度最长，再选择末端 UTR 最长，再 transcript_id 字典序；不按覆盖、方向或显著性挑选 transcript。plus 最后外显子是最高坐标 exon；minus 是最低坐标 exon。CDS 为最后 exon 内 CDS 交集，UTR 从最后 CDS 边界到 exon 转录末端。排除无末端 CDS、UTR/CDS 长度不足 50 nt、同一外显子内 CDS 多段导致解释歧义、不同 gene exon/CDS 有重叠（unstranded track 防污染）。

覆盖计算使用 exact bigWig intervals × overlap length，再除以完整查询区间长度，缺口视为 zero coverage；无染色体、查询失败单独记 missing，不能补零。主 gate：两个 control 的 CDS 和 UTR mean normalized coverage 均≥0.05；预设敏感 gate 0.01、0.1。gate 仅根据对照。零信号 treatment 的 log 比不用任意 pseudocount，记 nonfinite/不可估。

LECDS = UTR coverage/CDS coverage；各剂量 log2 fold = log2(LECDS_treatment/LECDS_control)。按两个 AMO 剂量均负、两个 OE 剂量均正记录同假说方向；阈值 0、0.1、0.2、0.5 log2（不是 DWUI 阈值，不能混用）。AMO/OE 平均剂量效应和 min(−AMO_mean, OE_mean) 同时报告，正反方向均保留。

本轮只做固定探针及合格 foreground/background 的描述性计数、连续均值/中位数和 coverage attrition；不把六个聚合轨道当独立生物重复，不做 treatment-level p 值。不从基因中挑出 RAB10 后升级整个集合。RAB10 仅作为原文已明确展示的机制先例，且重建 probe 与本研究 driver 不是同一 PAS 的证明。

## B. PAS 强度与间距结构

使用已冻结的 CFIm sequence-signature gene_features_and_labels.tsv.gz；该表以所有贡献研究 full-primary-axis 共同 terminal-UTR transcript 和完整 upstream100 资格选无响应的 native 最近/最远 anchors。仅连接751固定背景，不换背景。此处比 driverpair 更严格，不按响应挑位点；完全不适用于 ALE/intronic。

每物种预设三项：proximal_core_rank（高为有 canonical motif）、distal_core_rank−proximal_core_rank（正为近端更弱）、log(anchor_span)（双向，间距无预设单调机制）。core rank 来自 upstream [−35,−9)，AAUAAA=3、预定变体=2/1、无=0。共同 paired eligible 子集另外只作覆盖及描述性列表，样本过小不升级跨物种推断。

主零模型固定 measured_pattern、proximal_capacity、distal_capacity、PAS_number background median bin、baseline background median bin，块内打乱 foreground 二元标签，保留每块 foreground 个数。敏感零模型再加 mean_GEX、APA_coverage background median bin。不匹配被检验的 rank/span 特征本身。每 context 20000次，SHA256(context) 固定随机种子。检验统计为 foreground mean−nonforeground mean，双尾以相对块内 analytic expectation 的绝对超额计算。所有 2 species×3 features×2 models=12 项一起 BH；不从模型或物种中挑显著一项。另报告 null SD、可移动 foreground/背景；SD=0 记不可估 p=1，不称阴性。

## 验证与结论边界

验证 strand 和 CDS/UTR 半开区间规则；coverage exact interval 积分独立复核；gene/transcript 唯一选择与 hash；shuffle 每块边际保存；analytic expectation/Monte Carlo均值一致；p/BH 独立复核。保存来源sha、协议sha、机器表、覆盖/失败原因。最强结果也只是 U1 perturbation-compatible 或静态 PAS architecture 相关；没有本项目的 U1、新生 RNA 或 PolII 测量不能声称机制成立。负结果也不能排除 U1 蛋白装配、snRNA variants、局部暂停或其他机制。
