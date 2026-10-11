# 连续 marker 参照敏感性：冻结协议

日期：2026-10-10。本协议在新分类的 APA 结果计算前固定，是既有观察后的探索性补充分析，不称预注册。分析目录为 `notes/2026-10-10-marker-reference-sensitivity/`；只读当前冻结 GEX/PAS，不修改旧分析或稿件。独立方法反查位于 `notes/2026-10-10-marker-reference-review/`，注释证据核对位于 `notes/2026-10-10-marker-annotation-audit/`。

## 问题与可识别范围

检验：不依赖原粗标签、以同库连续多能程序/ECM 表达重新定义参照时，人路线较弱、小鼠较强的 P–M 近端响应构成是否保留；以及在原 M 内改变多能程序评分，原 P 相对 M 的 APA 变化是否被系统性压小。

这是一种表达参照定义敏感性分析，不是重新发现功能性多能细胞、纠正标签的真值检验或跨研究状态等价检验。marker 与原标签来自同一 GEX，删除标签不构成独立实验验证。M 是替代状态；P–M 不代表实际 M→P 轨迹。即使结果保持，也不能恢复“71.6% 是可推广物种步骤效应”的解释。基因、metacell 和库不作独立物种重复，不新增物种总体 p 值。

## 候选池与冻结来源

- 6研究、7计划路线：Guo2019 CIP、Schiebinger2019 WOT/2i、Xing GSE118258、Liu GSE150311、Guan GSE178325 的0618/0809、Panariello GSE221739。OSK 不加入本次计划路线。
- 重分类池为当前 native `metacells.tsv.gz` 中计划路线的**全部**既有 metacell，精确按 MC ID 对 frozen GEX，主池不限制原 P/M、粗标签、是否进入旧 selected 或旧 purity。
- 原路线路由、物理文库、时间、采样单位来自冻结元数据，不能按新评分或 APA 重路由。Pan D9A_rep1 未入旧 selected，但研究/文库记录唯一归属该研究路线，可进入全 MC 候选；不得冒称它存在旧 P–M 基线。
- 原 selected/node/purity 单列保存。原基线与桥接的 original P/M 沿用 `selected_metacells.tsv.gz`、purity≥0.8及计划同库 P–M 设计，不用全 MC 的其他标签补基线。
- GEX 必须是原整数 counts，CP10K 分母为该 MC 全部基因计数。核查 MC ID、行号、细胞数、总 GEX 一致；基因 symbol 用冻结 GTF 补全，严格唯一，不把多个基因 ID 拼成一个 marker。
- 所有 primary PAS 原坐标及稀疏 counts 一次只读取回，先在完整 primary 轴上固定 position，再以 matched_source_count_50nt≥1 保留 supported50 PAS。须有单染色体、单链和至少两个不同原生 cut。gene-wide DWUI 不自动等于 UTR 长度。
- 采集、协议、代码、参考、元数据与整数计数保存哈希；不重计读段，不更换参考。

## 固定 marker 与评分

主多能程序 marker 为 **DPPA2、DPPA4、SALL4、TDGF1**（鼠同名大小写对应），主 ECM 为 **COL1A1、COL3A1、DCN、LUM、FN1**。前三个来自09-13冻结 acquisition discovery，TDGF1和ECM来自09-26冻结程序。避开 POU5F1/SOX2/KLF4/MYC，及 Panariello 引入的 NANOG/LIN28A。程序评分不是功能多能性测量。

主 panel 称“共同内源程序 marker panel”，**不能称严格高置信一对一同源 panel**：DPPA2/DPPA4在旧广义 documented one-to-one 表中的 confidence=0，TDGF1旧映射不是 documented one-to-one；三者不在当前16092对严格同源宇宙。映射表逐项保留这些信息。它们用于库内表型参照，不用于同源基因重叠检验。

唯一预设替代 panel 严格复用09-26五基因程序：人 DPPA4/DNMT3B/ESRG/L1TD1/TDGF1；鼠 Dppa4/Esrrb/Zfp42/Dppa5a/Tdgf1；ECM不变。它是物种内原程序敏感性，绝对评分不直接跨物种比较。

对每个 route×物理库，用全部候选 MC：

1. 每个 marker 计算 `CP10K = 10000 * integer_count / total_GEX`，再取 `log1p(CP10K)`。
2. 每个 marker 在该库内按总体均值和总体 SD（ddof=0）作 z-score。常量 marker 记 z=0，保留在原 panel 等权分母；不得逐库静默删 marker。
3. 分别取多能、ECM marker 的等权 z 均值，得到两个独立连续分数，不用二者差值替代双程序筛选。
4. 每个程序必须有≥2个 SD>1e−12 的 marker，且其高分位严格大于低分位；缺失/ambiguous marker、总 GEX≤0或不满足变异要求时该库该 panel 不可分类，记录原因。

## 分类阈值：固定后不按 APA 调参

主版本 `core_third_strict`：每程序分别用 `np.quantile(..., method='linear')` 的 Q1/3、Q2/3。

- P*：多能分数≥Q2/3，ECM分数≤Q1/3；多能 panel 中≥2个 marker 的 CP10K≥1，ECM panel 中 CP10K≥1 的 marker≤1。
- M*：ECM分数≥Q2/3，多能分数≤Q1/3；ECM panel 中≥2个 marker 的 CP10K≥1，多能 panel 中 CP10K≥1 的 marker≤1。
- 其余 MC 为 unassigned/mixed，记录是相对分位或绝对支持未通过。P*/M*互斥；不强制每库都有两组，不从未分配 MC 补足数量。并列分数不按 ID/随机排名拆分。

CP10K≥1及对侧≤1 marker 是固定的表达支持门槛，不是生物学校准的分子纯度。特别是 ECM 共表达真实存在时，strict P*可能筛掉原 P；结果不能隐去这种 attrition。

预定敏感版本仅三个，完整报告，不能挑“最像15%/30%”的一版：

1. `core_quartile_strict`：同主panel，分位改Q1/4、Q3/4，其余相同。
2. `core_third_relaxed`：同主panel/三分位/高程序≥2表达 marker，取消低侧≤1 marker门槛；仍须低侧连续分数≤Q1/3，仍不造组。
3. `native_third_strict`：09-26原物种panel，三分位，所有支持/变异要求同主。

expanded-primary PAS 仅对主分类做测量/覆盖敏感性；不另选 marker 或阈值。所有 setting 的 position 必须源自相同 full-primary 轴。

## 同库 APA：原流程重放与桥接对照

先从原 selected、purity≥0.8 重放 supported50_p08 的 original P−original M，逐 gene/library 及 route 对齐当前冻结 pair/route effects；keys、贡献库、位置、向量与数值相符才能继续。

每个新版本计算四个对照：original P−original M、P*−original M、original P−M*、P*−M*。原状态只使用原计划同库参照；新 stars 可在所有计划路线的同库中形成，不拿不同库补缺失。

沿用原资格与层级估计：每基因每 MC supported PAS≥20 UMI；每被比较组至少3个基因合格MC和30个合格cells；两组累计至少2个PAS各有≥10 pooled UMI且≥2 nonzero MC。usage为合格MC的完整 retained-PAS比例等权均值；所有对照使用同一原生完整位置轴。保留 gene/group metacell IDs、counts、PAS IDs、usage、cells及失败原因。

路线效应复用 `polyapipe_effects.weights` 的 date→sampling unit→captures层级（同库对比两侧 date/unit 一致），不改成按MC数或深度加权；每gene以其实际合格库计算。方向阈值统一 ΔDWUI≤−0.02/≥+0.02，|Δ|<0.02保留；0.05另列分类敏感性。

**公平比较必须有共同贡献库与共同基因背景。** 每个版本分别为四个对照取逐 gene×同库完整资格交集，再用该gene完全相同的库和权重重算四对照；路线的四个方向比例使用同一N，效应均值也同背景。另列各自全部可测范围，不能用不同分母的百分比认定参照效应。不同marker版本有不同共同背景，不能交换分母。

输出每路线及每库的 N、近/远/阈值内计数和比例、mean/median Δ、近端条件幅度、P*/M*覆盖与原标签混淆表。共同背景 additionally报告“仅换M”和“仅换P”的变化及二者同时更换，不能将P*−M*变化直接归因于M。

在主版本四对照共测背景中，原P−M为 D00，原P−M*为D01，P*−原M为D10，P*−M*为D11；平均DWUI严格满足D11−D00=(D01−D00)+(D10−D00)，逐gene/库核验。方向比例因为阈值非线性不应强迫相加。

## 原 M 内的多能程序梯度

保持原 P 不变，且只在旧 selected purity≥0.8 的原 M 中分层。主多能连续分数仍用全库候选池的标准化尺度，但三分位切点在该库原 M 内计算，ties不拆。

- low≤Q1/3、high≥Q2/3、middle为两者之间；要求Q2/3>Q1/3，三层均有足够MC/cells，否则记录不可测。使用分数层级，不称“纯度梯度”。
- 对每gene，仅当原P及三层M均符合基因资格、三组P−M均通过联合PAS支持时纳入；每gene/库的原P usage在三对照完全相同。
- 三层路线效应使用完全相同 gene×library贡献集合及层级权重；路线summary同一基因N。另列全原M对照，不能用它不同的库分母替代三层配对背景。
- 主要量为 Δ(P−M_high)−Δ(P−M_low) 与两层近端比例差。若前者为正，表示相对高多能评分M的效应在位置均值上较不近端；若为负，则相反。三个点是否单调单列，不把任意变化命名为系统性稀释。
- 在同PAS轴/同P/同库情况下，上述差严格等于 M_low−M_high，独立核验。它反映APA与原M内评分的关联，不识别原M是否混入实际P，也不推导成熟度因果。
- 09-26原物种程序评分的同样梯度作预设评分敏感性，分母另外记录。

## 判定与报告

不以“人≈15%、鼠≥30%”作为通过门槛，不因覆盖少而修改marker、分位或PAS支持。

- 若共同背景的方向构成差异在两条鼠路线、可测人路线及预定敏感性中保留，报告其对这些表达参照规则稳健。仍不证明M等价、普遍物种机制或功能多能性。
- 若仅换M显著改变人路线的近端比例/均值，而仅换P改变较小，且M内梯度方向相符，可报告原M参照选择影响结果；“显著改变”这里只报告具体幅度和路线一致性，不新增生物学p或任意成功阈值。
- 若变化主要来自换P、版本间不一致或近远端同时变化，按实际桥接结果描述，不能套预设M稀释故事。
- 若 P*/M*或梯度覆盖不足，明确not estimable及缺失位置；不得把剩余少数人路线称“所有人研究通过”，也不能由未检验恢复R3强度。

输出完整分类、每库资格、PAS支持、背景流失、所有版本结果、原基线重放、意图测试、独立数值重建和输入/产物哈希。不生成MC伪重复p，不解释“鼠显著而人不显著”为物种交互。主文和旧报告暂不改写；本轮报告明确旧71.6%拆分仍依赖原M定义。
