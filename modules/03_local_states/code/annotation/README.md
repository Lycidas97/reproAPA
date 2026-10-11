# 统一注释与状态对齐 v1

本目录执行[六研究统一注释与状态对齐计划](../2026-08-31-unified-annotation-state-alignment-plan-v1.md)的输入冻结和 anchor 裁决阶段。

当前阶段只做以下工作：

1. 合并五个 R9 与 WOT R10B 的权威 cluster 注释；
2. 保留所有冻结 cluster、细胞计数和限制字段；
3. 生成 A/B/C、reference-only 和 validation-only 的 anchor 裁决底表；
4. 冻结第一版同物种状态族 registry。

当前阶段不运行 scVI/scANVI、moscot、RNA velocity 或 APA，也不修改任何 cluster membership。

运行：

```bash
python3 scripts/materialize_annotation_inputs.py
python3 -m unittest discover -s tests -v
```

权威输出位于 `release/`。

