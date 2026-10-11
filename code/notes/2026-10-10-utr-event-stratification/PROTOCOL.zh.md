# 串联 UTR 与末端外显子/UTR 切换：分层方向统计

2026-10-10。用户授权在现有数据上分别统计两类事件。本协议在本轮分层推断前固定；已有 27 对核心的 17/4/6 分型计数已知，因此属于补充探索性分析，不称预注册或独立发现。仅新增分析记录，不修改冻结来源或稿件。

## 问题与固定输入

回答 P 相对研究内 M 的共同近端响应，分别在串联 UTR 与可变末端外显子（ALE，作为本轮“UTR 切换”的操作定义）中是否超过条件匹配预期。M 不是已证明的 P 前体。ALE 可能改变编码序列，不能全部称为保持编码序列的 UTR 长度变化。

输入固定为 751 对双物种至少两项研究可评估的一对一直系同源基因、已封存的 3467 条 route driver inventory、九阈值复现成员与 assay features。使用 supported50_p08 原生 PAS 轴；不重新调用 PAS、改变状态或缺失资格。

## 分型

以 notes/2026-10-10-shared-proximal-mechanisms/event_inventory.tsv.gz 为源。核对 inventory 与冻结 recurrence 的实际贡献 study/route；如不一致，先报告并按冻结可评估贡献范围处理，不因响应方向删 route。

每个 species/gene 的所有可评估 route 均为 strict_tandem_UTR 才归 tandem；均为 alternative_terminal_exon_compatible 才归 ALE。其他单一类别保留其原名；跨 route 不同类别归 mixed。两个物种均 tandem 为 both_tandem，均 ALE 为 both_ALE，其余单列。全部背景和前景采用相同分型规则，分类一次固定，不随阈值改变。

driver 为既有最大 usage gain/minimum loss 位点对，与响应效应有关；这里是效应依赖的结构兼容性分层，不能称独立于响应的结构特征，也不能证明表达中的完整异构体。driver 符号不一定等于多 PAS 轴的 DWUI 符号，分别报告，不用 driver 符号重新定义方向。

另定义双方所有贡献 study 的完整原生 PAS 轴均兼容同一 terminal UTR 的 both_full_axis_tandem，并与 both_tandem 交集作为更严格敏感性。该资格不要求跨研究同一个 transcript ID，不能写为同一完整异构体。

## 复现与描述统计

沿用冻结成员：主阈值 |delta_DWUI|>=0.02，敏感性 0.05；扫描 0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05。每项研究所有可评估系列同向并达到阈值，至少两项研究支持才记复现；中间或未复现不计方向。直接读取冻结标签并以 study events 的 extrema/support 独立复核。

输出每个物种的分型分母、近端/远端复现数及比例，双物种分型交叉表、各类 3x3 方向表、核心 27/7 分型与成员名单。DWUI 描述采用 route 先平均到 study，再等权平均到 gene，不把 gene/route 当独立研究重复。近端/远端表示 PAS 轴位置；ALE 不直接解释为 UTR 碱基长度的缩短/延长。

## 类内匹配零模型

分别过滤 both_tandem 与 both_ALE 的前景及背景。保留先前 751 对总体产生的 quantile bins，不在子集中重新分箱。各类内部，两个物种各自在以下 frozen strata 内独立置换 proximal/distal/nonrecurrent 标签，保留所有 stratum 的三类边际。100000 次固定 seed 置换；记录解析期望、null 分布、上尾 p=(tail+1)/(100000+1) 及 Monte Carlo 精度。

- measurement_PAS_capacity：measured_pattern、proximal_capacity、distal_capacity、PAS_number_bin。
- full_assay：以上加 mean_GEX_bin、APA_coverage_bin、baseline_bin。

保存每类、每阈值、模型、物种的 singleton/homogeneous fractions、可改变标签数、null variance 与有效支持点。零方差保留计划家族位置 p=1、标为 nonmobile，不能作阴性证据。稀疏/低流动性与小 ALE 背景必须明确报告；不因结果再合并 bins 或类别。

主 0.02 的 proximal/proximal，2 类 x 2 模型 = BH4（q_primary4）。0.05 同样 BH4（q_sensitivity4）。九阈值上述全部 36 检验 BH36（q_scan36）。

另对四种有方向配对（包括 proximal/proximal）报告 secondary BH16 每阈值与 BH144 全扫描；其余 nonrecurrent 格仅描述。不能用 secondary 家族替换主家族或挑选有利阈值。两个模型是敏感性分析，九阈值不是独立复现。

完整 PAS 轴严格串联 UTR 交集只报告计数、类内匹配解析期望、mobility 及效应描述，不新增主推断。无正式类别 interaction：一类显著另一类不显著不代表两类之间差异显著。分析只检验分类条件下的基因方向对应，不推断同一 PAS、统一机制、APA 因果作用或物种层面的统计推广。

## 完成标准

751 对与 3467 route 输入完整对账；所有九阈值分类固定；0.02 核心 27 与 0.05 核心 7 对账；各类计数加总与旧背景一致。独立复核核心成员、方向/分母、类内解析期望及主 p、校正家族、输入 SHA。输出代码、结果表、中文报告和独立 QA，明确任何未通过或无法检验项。
