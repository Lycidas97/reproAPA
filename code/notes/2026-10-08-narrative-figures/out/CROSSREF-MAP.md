# 面板去向对照（供第 ③ 步改文字与替换 Fig/ 时使用）2026-10-09

## 新旧面板
| 原面板 | 新位置 | 说明 |
|---|---|---|
| 图 3c（示例表、散点、计数） | 图 3c | 示例表删去（逐基因数据在 Table S4e）；新增两个零模型的对比；散点保留 |
| 图 4a 研究对共享 | **S13a** | 替换原 S13b |
| 图 4b 全基因组相关 | 删去 | 研究对层面的 ρ 见新 S13d；均值写进正文和 Table S5a |
| 图 4c 近端共享 | **图 5a** | |
| 图 4d 同源 PAS | **图 5b** | |
| 图 4e *PDLIM7* | **图 5c** | |
| 图 5a 区间一致性 | **图 4a** | |
| 图 5b 位点净保留 | **图 4f** | |
| 图 5c 低净位移筛选 | **图 4c** | |
| 图 5d 部分保留 | 删去 | 同一组数据见 S17c、d |
| 图 5e *Eif4ebp2* | 删去 | 同类例子见 S17e、f；*Eif4ebp2* 留在 Table S6f |
| 图 5f *Riok3* | **图 4d** | |
| —— | **图 4b**（新） | 严格复现区间基因的全基因 ΔDWUI，只用 polyApipe |
| —— | **图 4e**（新） | *Fgfr1op2* 三项分解与区间 ΔI |
| S13a 阈值/层敏感性 | **只放进 Table S5b**（块 S5b.5，16 个设置，含 q_16，已在表中，无需改表） | 用户 2026-10-09 决定 |
| S13b 研究对同向组成 | 删去 | 新 S13a（原主图 4a）已包含 |
| S13c / d / e / f | S13b / c / d / e | 只平移，字母改动 |

## 需要同步改的引用（替换 Fig/ 时一起改；现在的源文件仍对应生产版图）
- **主文 `\ref`：** 正文的 `fig:cross-species-apa` 和 `fig:apa-regulatory-sites` 两个 label 不用改。调换 R4/R5 小节后，浮动体顺序会自动给出新编号。但面板字母要逐处改：
  - 跨物种一节：a、b → S13a/S13d；c → a；d → b；e → c；
  - 区间一节：b → f；c 不变；d、e → S17；f → d；新增 b、e；
  - R3 中的 3c；
  - R6 中的 `apa-regulatory-sites}c`（不变）。
- **主文中 S13 的引用：** "S13b--e" → "S13a--d"；"S13a,f" → "Table S5b、S13e"。
- **补充图图注：**
  - `supplementary-figures/main.tex`，第 203、223、232 行（Figure 4c/d/e → 5a/b/c）；
  - 第 275–378 行（Figure 5a–f → 4a/f/c/d 以及 S17；Figure 5e 的例子已删去）；
  - S13 图注整体替换为 `captions.tex` 中的新版本。
- **补充方法与结果：**
  - `candidate-intervals-and-predicted-sites.tex`：Figure 5、5a、5b、5c、5f；
  - `cross-species-comparisons.tex`（方法与结果两处）：Figure 4c、4e；
  - `display-conventions-and-reproducibility.tex`：Figure 5e、5f；
  - `alternative-pas-calling.tex`、`interval-inclusion-and-predicted-sites.tex`：Figure 5e、5f。
- **补充表索引 `supplementary-tables/main.tex`：**
  - 第 49 行 S5a "(Figure S13b--f)" → "(Figure S13)"；
  - 第 50 行 S5b "(Figure S13a)" → "(Figure 5a)"；
  - 第 71 行 Figure 5f → 4d。
