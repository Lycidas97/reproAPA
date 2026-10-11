# 协议澄清（不改变集合选择规则）

记录时间 UTC：2026-10-10T16:13:41.775793+00:00

1. coverage_protocol.json 是只读覆盖可行性预检查，提出原全24库固定面板的 date-blocked rank correlation。覆盖显示两预指定窗口共同可测基因均为0后，dose_protocol.json 在加载新的 APA outcome/关联之前冻结了另一个明确的探索性 estimand：同日期、同基因 panel 的 C1−C2 强度差与 predictor 差，以 date-fixed Pearson coefficient（不再次中心化 oriented differences）及单向量日内交换作检验。因此最终分析按更晚、更具体的 dose_protocol.json 执行，不报告初始 rank 方案已经完成。原全24固定核心检验仍 unavailable。

2. dose_protocol.json sensitivities 中“≥10 common-gene blocks (4 dates)”为纯计数笔误。明确选择规则始终是 ≥10，覆盖输入实际只有 D13=12、D14=22、D14.5=14，共3个日期；D18=9不纳入。执行 ≥10 原规则，不通过降低阈值补成4日，不改动带哈希的原协议。该敏感性仍只报告描述性系数，无新增 p/q。

澄清由覆盖数据和独立源表重建触发，没有按关联的方向或显著性修改选择规则。
