# QC v3.2 到 study-wide 分群生产说明

## 前置条件

六研究 QC 生产必须先完成，且每个研究都存在：

```text
independent_validation.tsv
primary_membership.tsv.gz
sensitivity_membership.tsv.gz
technical_membership.tsv.gz
tier_counts_and_hashes.tsv
```

六个 study 的 contact sheet 必须经人工视觉复核，并在 QC 项目中先执行：

```bash
bash production/run_96core.sh acceptance
```

只有 `runs/_production_acceptance/acceptance.json` 为六研究全 `PASS`，
且该报告不早于 validator、membership 和视觉复核台账时，才能生成
分群配置。

构建器会拒绝缺失文件、非全 PASS 独立验证、tier 数量/hash 不一致、重叠
membership、文库元数据不闭合或表达矩阵缺失。GSE225142 不在范围内。

## 服务器目录

```text
QC=/root/nfsdata/apadetect/analysis_projects/reprogramming_unified_rna_qc_v3.2_release
ROOT=/root/nfsdata/apadetect/analysis_projects/reprogramming_qc_v3_2_clustering_release_v1
PYTHON=/root/miniconda3/envs/single-cell-qc/bin/python
SNAKEMAKE=/root/miniconda3/envs/single-cell-qc/bin/snakemake
```

## 生成冻结输入合同

```bash
"$PYTHON" "$ROOT/scripts/build_qc_clustering_inputs.py" \
  --qc-root "$QC" \
  --clustering-root "$ROOT" \
  --baseline-config-dir "$ROOT/configs" \
  --output-dir "$ROOT/production_generated"
```

检查：

```text
production_generated/build.ok
production_generated/qc_to_clustering_registry.tsv
production_generated/config.production.yaml
production_generated/configs/*.yaml
```

registry 会根据新旧 primary cell-set hash 给出 `rerun` 或
`reuse_only_after_new_release_gates`。本 release 使用新输出根，不覆盖历史分群。

## 生产 dry-run

```bash
"$SNAKEMAKE" \
  --directory "$ROOT" \
  --snakefile "$ROOT/Snakefile" \
  --configfile "$ROOT/production_generated/config.production.yaml" \
  --cores 96 --jobs 96 \
  --rerun-incomplete --latency-wait 120 --dry-run \
  --resources mem_mb=512000
```

六研究完整 DAG 应包含 242 个 job，其中 120 个是重采样 shard。每个 study
的每个步骤独立调度；没有 agent 负责后台运行。

## 正式运行

```bash
cd "$ROOT"
nohup ./run_96core.sh > clustering_qc_v3_2.log 2>&1 &
```

资源上限为 96 核、512,000 MB。`--rerun-incomplete` 支持断点续跑，
`--keep-going` 允许其他研究在单个研究失败时继续。

## 分群合同

- primary 细胞用于 HVG、PCA、kNN、chooseR 和全量 Leiden。
- sensitivity 细胞不参与拟合，以 30 个 primary 近邻输出簇概率、熵和
  out-of-support 标志；不删除细胞。
- technical 细胞保留在 QC ledger，不进入生物学聚合。
- 时间、条件、作者标签、表达程序和 APA 均不参与分群参数选择。
- 真实时间、条件、cell-cycle、stress 和 interferon 仅在冻结分群后叠加。
- 评估子样与全量分群通过 ARI、簇数漂移和逐簇匹配核验；未达门禁的簇标为
  `rare_or_unresolved`，细胞仍保留。

## 完成与审查

计算完成标记：

```text
results_qc_v3_2_studywide_v1/release/.ok
results_qc_v3_2_studywide_v1/release/project_release.tsv
```

`.ok` 只表示计算合同闭合。任何 scope 的状态为 `REVIEW` 时，必须视觉审查：

```text
visual/partition_umap.png
tree/resolution_tree_selected_window.png
full_partition_audit/evaluation_to_full_transfer.png
program_audit/program_overlay.png
sensitivity/projection/sensitivity_projection.png
```

最终生物学验收还需要核对 marker、时间/条件组成、单文库主导簇、细胞周期/
应激程序和 sensitivity 投影是否共同支持该状态划分。
