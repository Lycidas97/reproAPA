# CFIm sequence susceptibility：分析前协议

冻结时间：2026-10-10；本文件写入前仅检查表结构、代码、参考坐标与不带响应标签的可测背景计数（human 389、mouse 390 strict tandem-UTR；两边同时满足为 67 对）。未查看 UGUA 的响应组均值或富集结果。属于已有方向结果之后提出的探索性机制检验，不是独立验证。

## 问题与边界

检验 P−M 近端复现基因是否具有较高的远端减近端 upstream UGUA 数量，且是否强于远端复现/非复现基因。近端 shift 不等于 CFIm 活性下降；UGUA 是 PAS enhancer 的结构标志，不能把本检验称为活性测量。阴性结果不能排除 CFIm 蛋白量、磷酸化、复合物装配或其他靶标机制。

主文献：Zhu et al., 2018, Molecular Cell, DOI 10.1016/j.molcel.2017.11.031（https://pmc.ncbi.nlm.nih.gov/articles/PMC5756121/），指出 UGUA 的峰在 cleavage site 上游约 50 nt，且 UGUA 不是所有 PAS 的必要元素；Brumbaugh et al., 2018, Cell, DOI 10.1016/j.cell.2017.11.023（https://pmc.ncbi.nlm.nih.gov/articles/PMC5766360/），指出 Nudt21 knockdown 缩短靶基因的 distal PAS 周围 UGUA 相对 proximal PAS 更富集。窗口以实验/序列位置而不是本数据的富集峰选择。本轮不读取外部 knockdown 结果。

## 锁定输入与序列规则

- 当前 v2_primary / supported50_p08 / DWUI；P−M 同一合格 physical library；跨研究复现用冻结的 ≥2-study 同方向且研究内全部可测 series 一致规则，阈值主 0.02。
- ortholog：16092 高置信 one-to-one 原始全表，不能换成筛后的共享基因表。
- 每个物种保留 P−M 至少两研究可测，所有贡献研究 full-primary-axis tandem_utr_compatible 的基因；必须有共同的兼容 terminal-UTR transcript，能包含所有贡献研究 full primary PAS，并包含两个 anchor 的完整 100-nt upstream window。严格子集避免 alternative terminal exon / coding / intronic PAS 混入 3′UTR 缩短。
- Anchors 不按响应选：在该基因所有可测贡献研究的 native supported50 PAS 联合集中，按 transcript strand 排序取最 proximal 和最 distal；重复 genomic cuts 合并。至少两个不同 cut，单染色体/单链。保留原位 native PAS，annotation catalog 仅用于核对支持与候选，不用 response-dependent gain/loss driver，不作 liftover。
- 固定原始 GRCh38_2024_A / GRCm39_2024_A FASTA，config 与资源 hash 均核对。Plus upstream100 = genome[cut−100:cut]；minus = reverse-complement(genome[cut:cut+100])。cut 是 interbase 切点，不能再人为加 1。terminal UTR 在表中用 1-based inclusive，转换为 [utr_start−1, utr_end)。只接受完整 A/C/G/T 窗口。
- 主窗口为转录方向 [-100,−40)（60 nt，覆盖报道的 upstream40–100 region）；敏感窗口 [-100,0) 和 [-80,−40)。TGTA 对应 RNA UGUA；允许重叠起点，只有 motif 四碱基全部在窗口内才计数。
- 主特征 = distal TGTA count − proximal TGTA count；同时报告两个 site 的绝对计数、≥1 与 distal>proximal 的比例，不能仅凭 distal 高就推断 CFIm-specific。
- 预定 general-PAS 结构对照：相同窗口 TGAT（同单碱基组成）与 CGTA（相同 GC 比例）；core polyA hexamer rank 在 [-35,−9) 依当前 cis RANK 计算，作 covariate，不当 CFIm 读数。

## 背景、标签与检验

先完成无标签的 sequence/UTR/anchor 资格，后连接 response label。共同严格、完整序列 ortholog paired background 是直接可比核心；各自严格完整的所有 one-to-one genes 是预定物种内扩展，不能将各物种单边背景的 fold 直接解释为 species difference。

标签为 proximal recurrent / distal recurrent / nonrecurrent（最后一类不等于生物学完全不响应）。另列 all-contributor event=0 的严格 no-event 描述，不替换主比较。预定比较为 proximal vs nonrecurrent、proximal vs distal；方向假说为近端组主特征更高。每基因一票。

匹配协变量：实际 unit:series_pattern availability signature、negative/positive capacity、PAS number、baseline M DWUI、native-anchor distance、upstream100 GC mean、mean GEX、APA coverage、distal-minus-proximal core hexamer rank。前两者固定类别，其他连续项按无标签背景确定 quantile bins。

- 主 null `joint_covariate_median`：固定 actual availability + capacities + PAS-number tertile + baseline median + span median + GC median + GEX median + coverage median + core-rank-difference sign，块内同时打乱完整三方向标签，保留每块三方向边际。主结果 BH 调整 4 contexts ×2 comparison =8 tests。小背景可能产生很多 singleton，必须报告 mobile genes/labels 与 null SD；SD=0 时 p=1 且 informative=false，不能将不可检验写为阴性。
- 敏感 null `assay_architecture`：availability/capacity/PAS tertile + baseline median + span median；用于判断 joint null 是否因过度分层失去信息，不可独自升级主结论。
- 预定探索性 `residualized_assay_architecture`：先用不读取标签的 least-squares 模型从 motif 特征回归掉 continuous baseline、log1p PAS、log span、GC、log1p GEX、log1p coverage、core rank difference 与 availability/capacity indicators，再在 assay_architecture 相同块内打乱标签。检验 residual motif-score 的均值差。这只控制线性 covariate association；不是在连续协变量上完全条件化。
- 每 null/背景 20000 次，固定 seed 由 context string SHA256 确定；固定 feature 数值，重排标签。精确 expectation 用块内标签概率计算；one-sided tail (1+count(null>=observed))/(20001)，额外报告双尾用于发现相反关联。始终保存 draws / stratum margins / group feature summary / 不可检验标记。
- 主窗口 TGTA main-null 为主 8 项；其他窗口/对照/敏感模型的全部 tests 在一个独立 exploratory family BH，不能从若干窗口挑一个显著结果做主结论。
- 共同背景 direct interaction = (mouse proximal−control motif mean difference)−(human proximal−control difference)，按两个物种各自同一预定 null 独立重排标签，保留 ortholog序列对；只在两物种比较各有 mobile labels 且 null SD>0 时给条件检验，two-sided，以 2 comparisons BH。不是物种总体推断；不比较两个 p-values。
- 样本量 attrition、标签边际、null 变异、within-block exchangeable proportions、covariate distribution 均输出，不因显著性放宽 gate、窗口或模型。

## 成功标准与解释

正负链与 cut/window 边界通过合成 intent tests；UTR 完整窗口不跨 CDS/外显子；序列资格不依赖 response；标签与当前复现 roster 全同；shuffle 保留每块三方向边际；analytic expectation 与 Monte Carlo一致；BH/tail 独立复核；资源和产物 hash 入 manifest。报告完整主/敏感结果和信息量限制。即使结构相关通过，最多写“与 CFIm-related susceptibility 相容”；没有 binding/perturbation 和特异结构对照，就不能区分 CFIm activity 与一般 PAS architecture。
