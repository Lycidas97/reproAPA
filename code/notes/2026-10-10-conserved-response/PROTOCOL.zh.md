# 跨物种同向响应的功能与位点证据：计算前协议

2026-10-10。用户已看过候选基因名字后提出功能假说，本分析属于事后、探索性补充；此协议在读取集合的功能命中、富集与位点组间结果前固定，不称预注册。目标是检验共享部分是否具有转录后调控相关的生物学身份，并区分同源基因响应、同源 PAS 变化和预测区段后果。

仅新增 `notes/2026-10-10-conserved-response/` 下的协议、代码、数据表、报告及验证；不覆盖冻结分析或稿件。功能、位点和文献由三个 agent 分工，根独立复核计数和统计。既有本地资源优先；若需要追加固定窗口/seed 映射，可对 `213docker` 上已存在的参考只读查询，不重计 reads、不更换版本。

## 1. 集合与术语

主口径为 polyApipe、v2 primary、supported50_p08、原研究级 recurrence 定义和 ΔDWUI 阈值 0.02。gene-wide PAS 位置变化不能直接称 3′UTR 长度变化；UTR 解释另设严格资格。

- 主集合为四个 same-contrast 同向单元并集：人鼠 P−S 均近端、均远端，及 P−M 均近端、均远端。四格计数为 7、14、27、3；去重后 48 对。一对一同源关系为当前 16092 对高置信、双向唯一映射。原名单遗漏的 P−M 双远端为 SNX5、FGFR1OP2、SAT1；不能因名单印象排除它们。用户原列前三组并集 45 作敏感性。
- 次级集合为人两个 hN 条件相对 primed 的等权平均近端响应与鼠 P−M 近端复现共享的 79 个唯一同源对。原 125 行是 79 加其内部双参照子集 46，不能作为 125 个不同基因。双参照 46 作敏感性。hN 不称跨研究复现或 same-contrast 保守事件，两个培养条件不当同质生物学重复。
- 0.05 阈值为预定敏感性；沿用原研究票与真实可测性，不将未测或缺失当无响应。主集合对应 12 对；hN 主集对应 16 对。
- 每个 ortholog 一票，跨多个方向/参照重复出现不增加样本量。每一格的名单另表保存，包括反向响应的描述性计数，不能将“不富集反向”写成“没有反向基因”。

## 2. 对照与测量匹配

主功能比较库：对每个同源对，保留 P−S/P−M 中实际双方各≥2研究可测的 contrast。case 为在其中任一 contrast 双方同向复现；control 为一物种在至少一个可测 contrast 复现、另一物种在其全部可测 contrast 都不复现。双物种异向或分别在不同 contrast 复现而未共享者不纳入单物种 control；无复现者也不纳入。只按实际可测 contrast 判定，不对缺失 contrast 填零。0.02 口径为 48 cases+539 controls，比较库 587；完整可测背景另外保留。

hN 比较库要求双方实际可测。case 为双方近端响应；control 为只一侧近端，另一侧既非近端也未达到远端响应定义。它是单侧响应 control，不是两个物种的单侧跨研究复现。主为 79+211=290，dual 为 46+189=235。明确记录 hN pooled 和 mouse recurrence 的不同含义。

协变量来源为冻结研究/库贡献表，不从基因名字补值：每物种 GEX、APA coverage、supported PAS 数在双方可测 contrast 间取中位数 proxy；hN 使用已补人 GEX 的 joint naive-direction 表。每个 trait 分别在当前 case/control 库内作保留 average ties 的物种内 percentile rank，两物种取 min，再按线性三分位切点分箱；不以随机排名拆 ties。分箱在任何 GO 资格过滤前固定。

主 coarse 层为上述三个 bin 加实际可测 contrast mask、已复现方向类别（proximal only/distal only/both）。main fine 保留 coarse 层并加两物种、两 contrast 的 exact measured-pattern；hN fine 保留 coarse，并加三个 trait 的双物种 max-rank 三分位。不得根据功能结果自适应松箱。协变量缺失独立 NA 层；本次已盘点主要协变量完整。

报告所有物种原始 GEX/APA/PAS、参考位置/方向容量的匹配前后 balance。核查匹配是否仍有明显残差，不把 min-rank 匹配称完整六维逐变量精确匹配。报告 case 可交换比例、预期换标签比例、singleton/frozen 层和方差；功效受限或无法交换时标 NE，不事后选择更显著的模型。

## 3. 固定功能类别和注释资格

三个固定根类别：RNA binding `GO:0003723`（MF），mRNA processing `GO:0006397`（BP），regulation of translation `GO:0006417`（BP）。使用冻结 GO-basic 2026-07-26 与 Ensembl 116 gene→GO 资源；成员为根及沿 `is_a` / `part_of` 的后代，解析 alt_id、排除 obsolete，不沿 `regulates` 边传播。类别可重叠，不将三个命中数相加当独立基因总数。

主注释统一为人注释；鼠注释及双物种 AND 是预设敏感性。每类别要求相应物种有至少一个有效 namespace 注释：无 MF 的基因对 RNA 类别是未知，无 BP 的基因对另两类别是未知；不把未知计为非成员。AND 要求双方 namespace 均已注释。过滤之后保留原分箱，重新计算该类别各层的 case 数，不重新定边界。报告 annotation attrition 和命中来源。

冻结注释只有 gene/GO 两列，没有 evidence code、qualifier/reference，不能声称实验性注释或 GAF 级 NOT 审核。旧 C4 loader 只读 BP，不能直接用于 RNA binding。三个固定根类别不套 10–500 的 term-size 筛选；即使大、空或不变也完整报告。三个类别适合“转录后调控相关因子”，不自动支持“RNA 加工机器的功能性反馈”。

## 4. 功能统计与多重检验

每层条件置换 case/control 标签，固定层总数及 case 数。类别命中总数 X 的零分布为分层超几何分布的精确卷积：`Hypergeom(N_s, K_s, n_case_s)`，主检验为 X 的富集上尾。输出观察命中、条件期望、观察/期望、95%条件零分布范围、p、方差和移动性；该范围不是生物学总体置信区间。空目标、无类别变异、全冻结导致不能检验时写 NE，BH 中占 p=1，不假称阴性。

主正式族仅 6 项：48 主集合及 hN79 次集合，各3类别，人注释/coarse/0.02，BH6。完整敏感性网格为四个集合（main48/list45/hN79/dual46）×两阈值×两模型×三注释口径×三类别，共144项；主6之外138项独立 BH，同时报告全144项联合 q。不得按敏感性结果替换主定义。

GO 全扫描作为探索：仅 BP+MF，main48/hN79、人注释/coarse/0.02，按每个实际 namespace-qualified 库的背景 term size 10–500 过滤，不按 case 命中数过滤。报告各集合及跨两集合的探索 BH，不将探索术语升级为固定假说。旧整个复现集合 C4 阴性与本条件比较是不同问题，不抹去旧阴性。

## 5. 位点层证据：严格资格与完整流失

位点分析另用未按 GO namespace 过滤的主/次及对照库。先保留各 contrast 的原始贡献；P−S 用标签对应的 source 设计，P−M 必须使用实际同库库对向量。旧 cis node-difference `driver_sites` 不能代替 paired P−M。先将原生 PAS 按冻结 annotation identity 合并 usage，保留未映射 mass；不将合并前多个原生位点当独立事件。

严格 3′UTR 解释要求每物种该 contrast 的全部实际贡献 study 的完整 primary PAS 轴都 tandem-UTR-compatible，且有跨 study 的共同 terminal-UTR transcript。不满足者标不适用，不称无保守响应。hN 按完整实际条件/库对检查，缺冻结 UTR 资格时不从别的状态补。

对应 PAS 使用冻结的双向唯一 1bp liftOver、同源 gene/strand 及全候选上已固定的 greedy 一对一分配；主双侧容差25nt，50nt敏感性，不能在名单内重分配。至少2对应PAS、转录方向次序一致；每个 study/state 的匹配 usage mass≥0.8。保留原轴 Δ 和匹配子轴描述，不用子轴重归一结果筛选父集合。

在每物种内，将实际贡献 route/series 向量按原规则合成 study，再等权平均 study；最大 gain 与最大 loss 从完整 usage 差独立确定，annotation ID 字典序处理平局。阈值 gain≥0.02、loss≤−0.02；最大 driver 不可映射时标不可判定，不改选第二大 driver。只有 gain↔gain、loss↔loss 对应，才记“双 driver 对应”。主 gene-response 并集按不同 contrast 分行，union 描述每 gene 一票。

“跨研究支持同一个对应切换”还须每物种≥2独立study支持同一对 driver 的同符号、达位点阈值变化；Guan系列不能作两研究。hN 次级只有固定条件，可报告条件内 driver 一致性，不称人跨研究事件复现。事件覆盖和对应性先作描述，不能把与父集合方向定义相关的事件比例 p 当独立保守性证明。

输出 gene共享→严格UTR→对应PAS→mass→driver可判定→双driver→逐研究支持的全部流失、失败原因和原/匹配效应。另列 full gene-wide/ALE 等候选描述，但不混入严格UTR分母，不因小样本放宽门槛。

## 6. 序列保守度与 miRNA 后果

现有 phastCons 来自同一个 mm10 60way track；两物种位点均映到 mm10，在 genomic `[cut−100,cut+50)` 的可用碱基上取均值。它们不是两个独立的物种保守度观测；负链也不改成转录方向窗口。高分不能证明 APA 响应受到选择。

若只读服务器可取得固定窗可用碱基比例及相同 mm10窗 GC，并能回放原分数，再在严格、可映射 driver 背景上开展探索性匹配比较。case 限于上节双 driver 对应，control 使用复现一侧的相同最大driver规则，要求其driver也有合格对应但不强迫另一侧有响应；两组均至少2可判定driver、每窗≥80%有限phastCons覆盖。每 gene 一票，先对一个对应PAS的人鼠分数合成一项，再平均独特driver位点，最后平均可测contrast；使用同一track，不能把两侧当独立观测。

sequence coarse 在固定功能 coarse 层上，加 gene-level GC 三分位、有限覆盖二类（0.8至小于1、完全覆盖）和对应PAS数三类（2、3–4、≥5）。sequence fine 再加完整native span的两物种min-percentile三分位、matched-site fraction三分位。新序列feature箱只在此可测case/control库中按上述固定规则形成、保留ties，不按phastCons数值或关联结果调整。主0.02、两集合（main/hN）、两容差（25/50）、两模型，共8项探索族BH8；未可测占p=1。统计量为每gene score的case均值，层内固定case数标签置换20000次，种子由固定前缀20261010及集合/容差/模型的SHA256派生，p=(1+上尾次数)/20001。输出条件期望、95% null范围和可交换比例；case少于3或零方差标NE。此比较是条件序列保守度关联，非事件选择检验。

若无法补齐这些控制，只报告冻结分数的描述性分布和缺失，不给“充分测量匹配的序列富集”p。追加窗口资源须记录原参考路径、sha、窗口、分子/非零支持无关、finite fraction和GC，不换新track。序列比较即使通过，只能支持该映射邻域的保守度关联。

miRNA 后果为候选个案及结构描述：只在可判定严格UTR事件区段内使用冻结 physical seed 坐标，报告各真实study区段 inclusion 改变、完整seed是否位于区段、家族。两物种同seed家族不等于同源位置；需要 entire seed 的唯一、gap-free、reciprocal 映射才称位置对应。无法做映射时保留预测、未判定。同一基因多seed去重，site不存在不可当证明无调控；结合/扰动缺失不称功能验证。EIF4EBP2 例子单独核对，不能因已有展示而降低本轮资格。

## 7. 文献、验证与叙事

独立文献核查 DDX6 多能可塑性、CSTF3/NUDT21 自身 APA、EIF4EBP2 miR-34/449：分别记录原始论文、物种/细胞、位点/转录本、直接实验及是否与本项目事件一致。RNA因子表达/扰动、其靶标机制、其自身APA功能不能相互替代。

意图测试验证祖先边/namespace与未知注释、并集去重、真实可测control、分层标签数、精确超几何及BH、位点未映射/复杂轴/series非研究、最大driver规则与usage mass。根另独立重建集合、类别成员/期望/p/q，并核查位点覆盖和实例坐标；完整哈希、输入版本、环境和所有 NE 归档。不将基因、MC或库作为物种总体独立重复。

v11 的“分歧之下的共性”结构按本轮作者意见保留；本轮不把负结果自动写成已排除周期/CPA/UGUA机制，也不把“不富集反向”写成所有共享基因不反向。若新增功能/事件证据成立，拟仅更新 R4、TH_conserved、K0 的证据与措辞；若不成立，保留功能候选在 TH_function。正文升级还须区分基因身份、事件结构和功能机制三个证据层级，不承诺RNA机器反馈或事件受到选择。
