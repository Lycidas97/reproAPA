# 第4–6项具体执行约定（运行新检验前）

输入仅用当前 polyApipe contract-r2 冻结表。第4–6项不是第1–3项强制替代调用验证家族。

A4：先按 supported50_p08 合格状态库与 frozen state_libraries 的记录日期检查 Guan2022/Liu2020。合格阈值维持每状态>=3 metacells、>=30 cells、基因PAS深度20及>=2 joint PAS；不能放宽状态阈值来获取第二研究。Guan同route同time的P/M无交集即不可行。Liu无 documented source，若存在同库P/M/observed-S，则构造observed-S敏感性三角，限定为改变S含义的敏感性而非documented-source重复。固定人类268个原P–M复现候选，使用同库全部三条边合格的triangle，每sample_unit再日期等权，|P-S|与|P-M|>=0.02且方向符合候选。提供人数/库/可测/对齐与Xing集合保留表；不作基因级显著性筛选。

A5：双参照组与同物种P–M复现候选互补组（不把整个420/268与子集当独立组）比较。先核对全候选预查71/222、5/128，再区分未建模、未合格、可测。区间率分母限定两指定路线都有 common-library MP interval eligibility，且 exact_site_cross_study 中至少一条 both_observed 可变注释区间的基因；分子为现有 recurrent_region_pairs 的MP严格复现基因。Fisher双侧，单位基因。位点净保留定义为共同可评估且物理坐标去重的预测位点，每基因先平均(P-M inclusion effect20、effect50、两路线)，再比较基因级均值（Mann–Whitney双侧）；常量保留位点记效应0是已测构造结果，不把不可测置0。类别/因子链接不伪作独立样本。只对至少3基因/组作推断；区间含位点且仅存在于既有union pool，结论限于该背景，并报告筛选所致分母损失。A5 family包含两个物种×区间率/净保留共4个预定检验，无法估计行p=1标明descriptive。

A6：固定k3 profile候选，WPP=.02且series_matching=True的P–S recurrence。可测分母要求该基因在GSE118258与GSE178325（人类）或GSE103221与GSE115943（小鼠）的 primary study_event_covariates 中均有该endpoint的有效记录，不能要求其跨阈值。C2与C1、C2与C3作Fisher双侧，4检验为A6家族。逐基因保存class、可测研究数、正负票与复现标记。报告未校正计数、可测分母以及复现率的方向，基因级率差不被解释成时间因果。没有改变profile筛选或聚类。

置换不用于这三个模块；所有输入和代码保存SHA256；raw_p/family输出供root统一BH。若可测不足，报告未完成推断，保持p=1行计入预定家族。

审计澄清（2026-10-08，首轮后）：experimental_time中的endpoint是采样阶段标签，不代表已知日历日期。A4“同日期”计数改标同记录时间标签；Liu敏感性实际仅用同物理库三角，未把不同库endpoint标签冒充同日配对。Guan可确认没有合格P/M的同记录采样时间，未记录精确日期不能自行推断。统计家族与A5/A6检验不变。
