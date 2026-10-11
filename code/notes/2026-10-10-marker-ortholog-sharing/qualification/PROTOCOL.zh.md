# Marker共享复现：资格与测量机会审计

2026-10-10。只读取旧复现定义、冻结同源背景与marker桥接的资格/贡献字段；不读取新方向交集或显著性，不运行机制分析。

主输入为 `unique_support_bridge_route_effects.tsv.gz` 的 `common_M_bridge`（原P固定、D00/D01共同gene×library贡献），其余原桥接用于一致性检查。主marker版本 `core_third_strict` / `supported50`；其余预设版本及expanded仅作资格敏感性。按旧751个P–M、1656个P–S同源背景固定身份，列出主背景内至少2项研究双方可测、全旧study/series机会完整及旧P–S/P–M交集资格的分母，不依据效应选背景。

成功标准：D00/D01各route贡献库与层级权重完全相同；缺测不填neutral；旧study内所有可测series的阈值规则保留；至少2studies同向且允许第3项反向的旧复现定义不变。Pan库/采样单位不增造study或series。明确旧测量proxy bins与新实际参照容量/机会signature的区别，提出输入字段与NE判定方案。只形成审计建议，不产出方向标签或p值。
