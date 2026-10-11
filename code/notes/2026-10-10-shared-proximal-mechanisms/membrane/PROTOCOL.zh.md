# 膜系统功能候选：探索性检验协议

冻结时间：2026-10-10T16:39:24.760602+00:00。27基因身份此前已被看过，本协议是在提出功能解释后、尚未计算新类别命中/显著性前固定，不能称预注册或独立发现。

## 问题与边界
检验共享P–M近端集合是否比同样可测同源背景更常注释到膜运输、ER和细胞骨架类别。这是响应基因功能的检验，不是APA上游机制检验，不能解释近端使用如何产生。SURF4/YIPF5的已知相互作用只可作为事后个案，不能计算挑选这一对的PPI显著性。不得追加类别或挑较有利来源。

## 固定集合与资源
共同背景751个P_M一对一直系同源基因，目标0.02阈值27个，0.05阈值7个仅描述。来自2026-10-10-rna-module-enrichment/ortholog_roster.tsv.gz与target_members.tsv；独立从原measured_ortholog_members.tsv.gz核验。人GO主注释，鼠GO敏感性，使用项目冻结Ensembl116与GO2026-07-26的GO ID表、本体、manifest，按is_a/part_of闭包；obsolete/别名处理沿用已审计工具。无注释不是无功能。

固定四类（不合并）：ER membrane GO:0005789；ER to Golgi vesicle-mediated transport GO:0006888；vesicle-mediated transport GO:0016192；cytoskeleton organization GO:0007010。宽类别彼此可重叠，祖先传播记录证据路径。

## 零模型与校正
复用冻结joint human+mouse匹配分层：COARSE measured_pattern、proximal_capacity、distal_capacity、PAS_number_bin；FULL另外双方mean_GEX_bin、APA_coverage_bin、baseline_bin。固定每层N和目标n，各层命中超几何精确卷积单侧上尾。COARSE为主、FULL为低功效敏感性；零方差明确标不可区分。每注释×模型四类q4，两个注释×两个匹配模型16项q16共同给出，主要判断要求COARSE人注释q16≤.05且不是零方差；鼠注释支持程度另报。未匹配Fisher仅描述，另q8，不用于替代匹配失败。所有零命中/退化项保留家族。不给7基因新的功能显著性。

## 外部ER关联预查和可行性
检索Cheng2021（PMID34289366，PMC8501909，GSE162971）的全量PAS/分级读数或作者gene-level表。先下载资源、核实物种/组别/数值定义及tested universe，不提前打印27命中结果。若有完整可评估表，可另冻结分析增补协议再做外部ER相关异构体检验；若只有原始raw counts，且需要新原始测序流程才能获得可比较APA，报告当前不能直接检验，不把查不到作为阴性。此可行性预查不进入四类p值。

## 验证与报告
重算751/27/7与固定输入hash；独立解析GO父图与注释，独立统计卷积并逐项比较观察数、期望、p/q及target名单；toy exact与退化边界测试。有缺失/覆盖不足全部记录，发布机器表与REPORT.zh.md。来源与协议hash在freeze.json。无稿件修改。
