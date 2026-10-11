# 周期与来源状态归因：冻结补充协议

2026-10-10。本协议写于覆盖/设计预检之后、新 S−M 效应和本轮周期—APA关联读取之前。用户给出的粗略关联、旧 H1 的评分差异和既有 WOT 检验的设计/覆盖已经已知，因此本分析是探索性补充，不是预注册。仅分析冻结原始 selected、purity≥0.8 的 S/M/P；不使用 marker 重定义状态，不改正文，不提交 v11，不进入 W2。

## 问题和共同口径

检验原先提出的周期解释能否同时解释不涉及 P 的 S−M、文档化来源的 P−S 和 WOT 的同日期 P−M 剂量关系。ΔDWUI<0 表示近端；周期转录评分不能直接解释为分裂速率。M 是研究内替代状态，S−M/P−M 不是已证明的沿途转换。来源、物种、培养和时间仍混杂，任何回归系数都不能识别物种固有或因果效应。

沿用 proliferation-apa 冻结的 10 S/8 G2M 标记及公式：MC 全 GEX 为 CP10K 分母，模块 mean ln(1+CP10K)，Cycle 为两模块等权平均；Adjusted_Cycle 沿用已冻结的表达匹配背景。MC 等权。APA 采用 supported50_p08、冻结完整基因位置轴，不按当前支持 PAS 重缩放；每 gene/MC≥20 PAS UMI、每状态≥3合格 MC且≥30 cells，联合≥2 PAS 各合计≥10 UMI且≥2 MC非零。状态本身先通过≥3 MC/30cells。缺失不补零，不放宽门槛。

新跨状态分析主评分与每个 APA pair 的 gene-qualified MC 完全一致，辅助评分使用同库原始 pre-gene 状态 cohort。每基因先按原流程 time_a→sample_unit_a→time_b→sample_unit_b→capture 的层级权重汇总，再等基因汇总。保留原始与背景校正周期，输出具体贡献库、时间、合格 MC 和权重。既有 WOT core 正式检验保持其原 pre-gene Cycle 定义，不能悄悄改成新的 gene-qualified 评分。

## A：S−M 核心响应

固定 core27：recurrence-utr-ortholog/measured_ortholog_members.tsv.gz 中 P_M、threshold=0.02、人鼠均 proximal 的 27 对一对一直系同源基因。固定完整集合：directional-recurrence/observed_members.tsv.gz 的 P_M、threshold=0.02、proximal；human144/mouse372，按原主文研究支持定义，不另改严格程度。它们由 P−M 结果选择，WOT等参与发现，因此 S−M 结果不是独立验证，也不证明周期模块身份。

主设计按 route 使用所有文档化来源 S × 全部合格 M（Xing、Guan0618、CIP、WOT）；没有文档化来源的 Liu/Pan 使用所有 observed S × M，单独标明。Guan0809无 M 记不可估。计算时 M 为层级加权 target、S 为 reference，最后取负得到 S−M，因此 M 日期按原层级等权。每个基因只汇总其可测 pair，记录覆盖和时间组成，核心名单不按结果修改。

报告每路线/集合的 N、S−M mean/median DWUI、±0.02 近端/远端/阈下数、同贡献 ΔCycle/Adjusted_Cycle、逐库/日期稳定性。方向预测用正 ΔCycle 应对应负 ΔDWUI；WOT/Liu 在用户非配对均值中接近零只是已有假说，不能先验称为等效零或纯负对照。没有生物重复支持的 gene-sign p、MC p 或物种 p 不计算。

固定敏感性：仅既有同库 P−M 设计所包含的 M 库（late_M）桥接；同物理库 S−M；同时间标签不同库 S−M（Liu endpoint是条件标签、D0_D7_pool为混合时间，不称同步）；pre-gene Cycle；背景校正 Cycle；在同一设计所有候选 pair 均可测的固定基因交集。交集为空明确记不可估，不改窗口。每设计单独报告，不能把一物理对的多种标签当重复。

另为 B 的所有路线补充 observed P−S：只在 Liu/Pan 用固定五条文档化来源 P−S 共同的1499对高置信一对一同源背景，按原始 S×P 合格群组和同一全轴/资格重算，P 为 target。记录每基因贡献，原始周期及背景周期；返回给 B。这是 observed S 对比，不能并入文档化来源主模型。覆盖可在读取新效应前决定可测集合；不按效应选择。

## B：路线/对比剂量和来源差异

主分析是五条文档化来源 P−S 路线（Xing、Guan0618、Guan0809、CIP、WOT）共同的1499对高置信一对一同源背景，不用方向或复现筛选。按冻结 pair membership/weights 重建 ΔDWUI 和 exact-MC ΔCycle，核对 route_effects。每个路线是一个描述性点；Guan 两 series 各权0.5，其他三个研究各权1。模型 Y=intercept+ΔCycle，以及 Y=intercept+ΔCycle+mouse indicator。输出系数、拟合/残差、两物种固定研究加权差、秩、条件数、leave-one-study-out 稳定性；原始/背景校正及 pre-gene Cycle 同协议敏感性。没有 p/生物总体 CI，不把系数变化转成“解释了百分之多少”的因果归因。

补充：native 全可测基因背景；加入 Liu/Pan observed S 对比时采用七路线可测共同同源交集并保留 observed 标志，另表展示且不当主来源模型。P−M 五条同库路线（Xing、Liu、Pan、CIP、WOT）用其共同可测一对一同源背景以及 native 背景描述，Guan无同库 P−M，明确记录。重点核查 Guan0618/CIP 在完全匹配评分后是否仍构成反例，不事后丢弃。

既有三路678对共同基因的 P/S/M matched triplets 单独重建 exact-MC 三状态周期，P/S来自triplet MC IDs，M来自同 gene/同 P 库的原 P−M reference IDs，权重与triplets完全一致。验证 ΔCycle_PS=ΔCycle_MS+ΔCycle_PM 及 APA 恒等式。九点相关只能作为描述性重算，不能当9个独立观察拟合 p，也不以这个 r 判定归因。另输出两个代数独立对比 MS/PM，不增加研究数。

判读不以单个回归系数的符号或 R² 阈值决定叙事。看 matched 数据中预测的方向、反例、背景/评分敏感性和 WOT 同日期结果能否共同支持周期解释。五条路线/四研究不支持因果物种调整；即使残差排序仍在，也只能说周期评分未充分解释所采样路线差异。

## C：WOT24同库剂量

覆盖预检：24库、13日期（11双库日期，D12/D12.5单库）；固定27在所有24库共同可测为0，因此“全24完整27固定模块检验”不可估。不使用有利窗口替代。全24共有固定 native assay85基因，751同源背景共有45基因，完整mouse372近端集合共有8基因，均由资格而非效应筛选。

新描述性全24分析主 outcome 为固定85基因 mean P−M DWUI 的负值（近端强度）。在全部24列表保留两个 singleton，但同日期消除 date 后系数仅有11个双库日期贡献。每日期等权，C1−C2 oriented ΔCycle 和近端强度之差，r=sum(XY)/sqrt(sum(X²)sum(Y²))、beta=sum(XY)/sum(X²)，不再次居中 oriented differences。主 Cycle 为85基因相同 APA-qualified MC/贡献的等基因周期差；原 pre-gene Cycle 与 Adjusted_Cycle为固定敏感性。另在固定45同源和8近端panel上同式重算，报告 mean/median outcome、批次常量差消除、leave-one-date-out范围。不增新 p，不用最有利panel作结论。

既有 core 次级正式检验严格独立重放 cpa-cycle-core-dose/dose_protocol.json（及其书面 clarification）：7日期、14库、每日期≥5共同core；每日期两库同gene panel，跨日期panel不同，不估计完整27均值。由原始 GEX计数独立重建18marker Cycle 和11CPA因子 score，重建 pair outcomes，128 within-date predictor joint sign-flips，双侧 p=tail/128（含observed），保留原 Cycle/CPA family2 的BH q。原≥10基因敏感性实际3日期，既有protocol的4日期是已澄清计数错误，不改阈值。

任何条件检验依赖 C1/C2 交换性及日期独立假设；同实验/细胞系和可能持续batch仍是限制，constant-batch-centering仅诊断常量偏移。WOT参与发现core，结果探索性，不能单独证明增殖机制或解释人鼠分歧。

## 验证与交付

三个子目录独立写代码/数据/报告；总协议先算 SHA256 再授权读取效应。必须核对身份、coverage、score分母、原权重、冻结 effect replay、三状态加和、同日期固定gene panel、缺失无零填、singleton零贡献、WOT exact128/BH、符号和层级权重意图检查。root独立复算关键聚合/模型和统计，保存输入/代码/协议/输出 hashes、验证结果和总 REPORT.zh.md。不得因未显著而排除周期机制，不得把不可估写成阴性。v11保持未提交，W2保持未开始。
