# 连续 marker 参照敏感性：独立协议核对

2026-10-10。仅阅读已有注释、审计、固定基因映射和可测性代码；没有读取本轮重分类后的 APA 结果。本文件给新分析提供口径核对，不修改旧报告、标签或正文。

本轮能回答的是：同一批实际 metacell 上，换用预先固定的表达程序参照，原 P–M 的近端响应构成是否改变；改变来自重新选择 P、重新选择 M，还是覆盖范围。它不能认证功能多能性、统一 M 身份、跨研究状态等价或状态转换路径。

## 1. 建议锁定的最小可执行规则

**主多能性相关 panel：DPPA2、DPPA4、SALL4、TDGF1；ECM panel：COL1A1、COL3A1、DCN、LUM、FN1。** 按本物种原生 Ensembl ID 提取；不能仅靠大小写把符号当作一对一同源身份。每个 feature 必须有明确 ID 和覆盖记录；不存在的 feature 不能补成表达零值。

这里应称“共同内源 marker panel”，不能称“严格高置信一对一直系同源 panel”：本轮 16,092 对严格同源宇宙只包含 SALL4 和 ECM 五个；旧广义 Ensembl 116 表里的 DPPA2/DPPA4 为 documented one-to-one，但 confidence=0；Tdgf1 在该表没有人端匹配。保留四个原生程序 marker 用于参照敏感性是可行的，须把上述映射资格如实导出，不能继承主文同源重叠检验的身份资格。

预设 species panel 敏感性严格复用 09-26：人 DPPA4、DNMT3B、ESRG、L1TD1、TDGF1；鼠 Dppa4、Esrrb、Zfp42、Dppa5a、Tdgf1。它有物种及 naive/primed 程序差异，因此只作敏感性，不根据新结果替换共同主 panel。

主评分排除 POU5F1、SOX2、KLF4、MYC、NANOG、LIN28A 等可能引入的因子。尤其 Panariello 原研究转染 LIN28/NANOG，现有计数不能区分内外源；其高检出不能单独证明内源多能性。人 DNMT3B、ESRG 和鼠 naive marker 不进入共同主 panel。

**候选单元：** 本轮冻结目标路线内全部原有 metacell；不先限制 cluster 为 P/M，不用旧 label purity 筛候选。保留物理文库、路线、时间、供者和原标签，后者仅用于重分类覆盖与桥接对照。旧基线保持原 selected、purity≥0.8 资格；不能悄悄改变旧分母。

**计分：** 从整数 GEX 以全基因总 UMI 为分母计算每 metacell 每 marker 的 CP10K，再取 log1p；在每个 route×physical-library 内逐 marker 标准化，panel score 为固定 marker 数的等权均值。零方差 marker 保留在分母并记 z=0，另标不可变；不按库静默删 marker。每个 panel 至少两个非零方差 marker，否则该文库不可分类。两个 panel 各自要求 Q2/3>Q1/3，否则不可分类。z-score 不是身份概率。

**主 P\*：** P score≥Q2/3，ECM score≤Q1/3；P panel 至少两个 marker 的 CP10K≥1，ECM panel 至多一个 marker CP10K≥1。**主 M\*：** 对调这两项程序。CP10K≥1 的门槛是每个 marker 的名义表达资格；“至多一个”是 marker 数，不是 panel 总表达量上限。分位计算保存精确方法与阈值；ties 不拆分，不随机指定成员。混合、高高、低低、阈值/覆盖不足均保留 unassigned。

分位规则不能强制每个文库都生成两组。单态文库可以只有 P\*、只有 M\* 或两者均不可测；不拿相对低尾补另一组。绝对 marker 检测门槛降低这种风险，但没有经过外部生物标定，不能命名为分子“纯度”界限。

固定三项敏感性即可：共同 panel 四分位；共同 panel 三分位但取消对侧检测个数上限；species panel 三分位。相同公式、缺失处理及 PAS 口径，不根据 APA 方向另增阈值或切换主方案。

## 2. 效应重算与共同分母

沿用 v2_primary 的完整 primary 原生 PAS 位置轴，再选择 supported50 位点；至少两个不同有向 cut，同染色体同链。每个 metacell 每基因≥20 PAS UMI，文库内每组每基因至少三 metacell 和三十细胞；无资格记 NE，不记零效应。每库先按旧代码重现 baseline，route 的日期→sampling-unit→library 权重沿用当前冻结实现。

每个方案输出四种对照：original P−original M、P\*−original M、original P−M\*、P\*−M\*。它们分别识别重新选择 P 与 M 的影响。所有四种对照的直接变化必须在同一基因、同一实际贡献文库交集上计算；否则覆盖与表型选择混在一起。每个方案的扩大可测背景另列，不冒充 fixed-background 稳健性。

近端/远端主要以当前 |DWUI|≥0.02 规则计数，0.05 为敏感性。分母包括该共同背景全部可测基因，不能仅在改变基因中取比例。同时报告均值、中位数、近端和远端幅度及 subthreshold 比例。人、鼠每条实际路线单独列，Guan 两系列不算两个研究；不提供以 gene/metacell/library 为物种重复的 p 值。

重分类附表必须包括：全部候选数、各新组数、原标签→新组列联、原 P/M 保留比例、metacell/cell/library 覆盖、按时间与供者的覆盖、score 与原始 CP10K 分布、未分配/不可测原因、各基因实际贡献库和 PAS 深度。先保存这些表达结果及输入 hashes，再计算 APA。

## 3. 原 M 内梯度

仅在原 M 中，按主 P score 的同库三分位分层，ties 不拆；原 P 固定不变。每层都按原 PAS 门槛测量，所有层共用同一基因与物理库分母。原 M 的三层应输出其连续 P/ECM 表达、深度、细胞数和日期组成。若三层无法同时测量，该库/基因 NE，不能改为两层、移动门槛或补别库。

该梯度测试“参照群体的表达组成敏感性”；不是按 M 内分数证明一个 M→P 过程。marker 分层本身会选择表达、周期与深度背景，因此梯度与 APA 的对应仍不能证明状态纯度造成效应。

## 4. 判定与正文用语

不把“人约15%、鼠≥30%”定为通过标准。应预先比较：每条路线近端响应差异是否保留、固定背景变化量、变化主要由哪个桥接对照承担、方向是否跨固定敏感性一致、以及关键路线是否仍有覆盖。仅剩单条人路线或单条鼠路线，不称跨研究验证成功。

- 若 P\*−M\* 差异保留，且 original P−M\* 亦保留：可写“方向构成差异对这套连续表达程序参照定义不敏感”。不能写“标注错误已被排除”“M 已跨研究等价”或恢复 71.6% 为物种机制分解。
- 若 original P−M\* 明显改变而 P\*−original M 基本不变：支持原 M 的选择对 P–M 大小有贡献，保留 P–S 为主线并降低拆分解释。
- 若变化主要出现在 P\*−original M：提示 P 的重新选择/表型限制更重要；不能把 P\*−M\* 的变化全归因于 M。
- 若覆盖严重不足：该敏感性不可评估；不能把 NE 解释为方向消失、假说被反驳或差异稳健。

## 5. 已有证据中必须收紧的说法

1. **P marker 证明到程序层级。** 09-13 明确这些 marker 已参与原注释；新的组合评分仍来自同一 GEX。取消 cluster label 不等于获得独立身份真值，更不能靠 marker 高检出排除功能多能性、成熟度或 naive/primed 差异。人的 P–M 比例和 marker 完整度“没有关系”也没有正式可识别检验：路线数小，marker 不同、参照 M 不同，并有转染因素。
2. **label purity≠分子纯度。** metacell 纯度高表示其组成细胞原标签一致，不保证 ECM/多能性程序单一。提高 purity 到0.9检验 label 混合敏感性，不能据此排除真实混合表型。
3. **近邻不是等价性。** 09-26 的 15-NN+距离界限是描述性共同支持筛查；研究/技术效应与状态差异未被区分。Guan P 的邻居标签均 P、但距离不通过，并不能解释为同状态被验证或重编程失败。
4. **人 M 投票0是中位数。** 原表 Xing M 130个中1个、Guan M 433个中61个有非零M票；联合q95支持均0。准确写法为“同标签投票中位数为0，均未通过共同支持筛查”，不能写成所有单元完全没有同标签对应。
5. **P–S是不同资格的来源比较。** 正文主方向由有文档化起始态的路线支持；Liu observed S、Panariello D3不是严格诱导前来源。P–S不依赖M不等于所有P–S端点已达到跨研究等价、功能成熟或与物种标签独立。
6. **71.6%依赖中点选择。** P–S=P–M+M–S恒等式本身不能赋予M时间顺序。即使本轮重定义后比例保留，也只能作为采样路线及本参照定义下的描述性分解，不能升级为“物种差异发生在哪个重编程步骤”。

## 阅读来源

- `notes/2026-09-13-state-marker-supplement/README.zh.md`、`inputs/config.json`、原S5与原marker表。
- `notes/2026-09-14-state-grouping-robustness/README.zh.md`。
- `notes/2026-09-14-integration-reliability/A07-all110-program-marker-results.zh.md`、`A07-GSE221739-all45-review-results.zh.md`。
- `notes/2026-09-26-state-comparability-audit/PROTOCOL.md`、`README.zh.md`、`extract_remote.py`、`expression_support_summary.tsv`、`metacell_expression_support.tsv.gz`。
- `notes/2026-09-26-global-direction-audit/PROTOCOL.md`、`README.zh.md`。
- `table/production-v15/data/9cae42683c5e-ortholog_universe.tsv` 与旧广义Ensembl116原生gene映射。

未读新 APA 结果；未改正文、原标签或旧数据。
