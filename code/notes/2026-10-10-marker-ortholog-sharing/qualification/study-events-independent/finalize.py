"""Seal the independent opportunity/usage/events/member audit; no null inference."""
import hashlib
import json
from pathlib import Path

import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
TASK=HERE.parents[1]
read=lambda p:pd.read_csv(p,sep="\t",float_precision="round_trip")


def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):h.update(block)
    return h.hexdigest()


def main():
    replay=json.loads((HERE/"QA.json").read_text())
    comparisons=json.loads((HERE/"comparison_QA.json").read_text())
    groups=json.loads((HERE/"group_request_QA.json").read_text())
    scalars=json.loads((HERE/"scalar_GEX_QA.json").read_text())
    assert comparisons["status"]=="PRODUCTION_EVENTS_AND_ROOT_LABELS_PASS"
    assert set(comparisons["comparisons"])=={"root_labels","root_3x3","production_labels","production_core27_flow","production_primary_shared_members","production_events"}
    frozen={}
    for filename in ["input_sha256.json","group_input_sha256.json"]:
        for path,digest in json.loads((HERE/filename).read_text()).items():
            assert sha(ROOT/path)==digest,path
            frozen[path]=digest
    for path,digest in scalars["input_sha256"].items():
        assert sha(ROOT/path)==digest,path
        frozen[path]=digest
    for path,digest in comparisons["input_sha256"].items():
        assert sha(TASK/path)==digest,path
        frozen[str((TASK/path).relative_to(ROOT))]=digest

    labels=read(HERE/"paired_labels_independent.tsv")
    core=read(HERE/"original_core27_flow_independent.tsv")
    old=read(ROOT/"notes/2026-10-10-ortholog-threshold-scan/measured_ortholog_members.tsv.gz")
    old=old[old.contrast.eq("P_M")&old.threshold.eq(.02)]
    symbols=old.set_index("human_gene_id").human_gene_symbol
    shared=labels[labels.human_direction.eq("proximal")&labels.mouse_direction.eq("proximal")]
    count_table=shared.groupby(["threshold","definition"]).size().unstack().reset_index()
    count_table.to_csv(HERE/"shared_count_summary_independent.tsv",sep="\t",index=False)
    primary=labels[labels.threshold.eq(.02)]
    original_ids=set(core[core.joint_measurable].human_gene_id)
    flows=[]
    for definition,q in primary.groupby("definition"):
        selected=q[q.human_direction.eq("proximal")&q.mouse_direction.eq("proximal")]
        selected_ids=set(selected.human_gene_id)
        flows.append(dict(definition=definition,shared_proximal=len(selected),old_core27_measurable=len(original_ids),old_core_measurable_retained=len(original_ids&selected_ids),old_core_measurable_changed=len(original_ids-selected_ids),new_members=len(selected_ids-original_ids),human_proximal=int(q.human_direction.eq("proximal").sum()),mouse_proximal=int(q.mouse_direction.eq("proximal").sum()),lost_direction_symbols=sorted(symbols.loc[list(original_ids-selected_ids)]),new_member_symbols=sorted(symbols.loc[list(selected_ids-original_ids)])))
    main=next(x for x in flows if x["definition"]=="human_Mstar")
    assert main["shared_proximal"]==31 and main["old_core_measurable_retained"]==23
    assert main["lost_direction_symbols"]==["YWHAB","ZC3H15"]
    assert main["new_member_symbols"]==["ARPP19","ASPH","GLUL","MTDH","OAZ2","RTN3","TPM1","ZMAT2"]
    qa=dict(status="INDEPENDENT_QUALIFICATION_USAGE_EVENTS_MEMBER_AUDIT_PASS",failed=0,scope="common background, original usage/library weights, study events, recurrence/support/member tables, local group identities/denominators/scalar CPM",source_replay=replay,production_comparison=comparisons,group_request_review=groups,scalar_review=scalars,frozen_input_files_verified=len(frozen),frozen_input_sha256=frozen,primary_flow=flows,new_null_computed=False,raw_GEX_numerators_reexecuted_here=False,remote_GEX_spotcheck_delegated="cis reviewer 15 groups",covariate_null_and_multiplicity_audit_delegated="cis reviewer",interpretation="The 31-member identity overlap is independently verified; enrichment significance is outside this audit.")
    (HERE/"FINAL_QA.json").write_text(json.dumps(qa,indent=2,ensure_ascii=False)+"\n")
    table_lines=["| 阈值 | matched原M | 人M* | 双物种M* |","|---|---:|---:|---:|"]
    for r in count_table.itertuples(index=False):table_lines.append(f"| {r.threshold:g} | {r.original_matched} | {r.human_Mstar} | {r.both_Mstar} |")
    report="""# 共同背景、usage、study事件和成员独立验收

2026-10-10。结论：本分支职责范围全部PASS，未发现资格、方向、贡献权重或成员分流错误。只写本子目录，不改封存的资格报告、marker输出或正文；不运行或验收null显著性。

## 数据链与生产逐值比较

- 在封存728对共同背景上，独立读取30,470条library对比的原P/参照PAS usage和固定position，重算ΔDWUI；library层最大误差0。
- 按target日期、target样本、reference日期、reference样本、capture逐级等权独立重建权重；D00/D01同gene共同库集合、原P usage、权重保持相同，权重最大误差0。
- 重建6,638条route效应、59,742条source D00/D01 study阈值事件，源数据链139,953项检查通过。route与已封marker效应最大差2.47×10⁻¹⁶。
- 三定义展开后89,613条生产study事件的event、全部series pattern、effect mean/min/max、baseline和正负理论容量逐值一致，最大浮点差2.47×10⁻¹⁶。投票由所有可测series跨阈值定义，未用series均值替代，未把Pan capture增加成独立研究，缺失study未填neutral。
- 19,656条配对标签与根独立表及生产表完全一致；生产两物种方向、近端/远端支持研究数、可测研究数、opposing标记共10列一致。243格3×3计数与根表一致。
- 729条原27分流、87条主阈值共享成员记录与生产逐值一致。GNPNAT1和KLC1在所有27个定义×阈值context中方向均为NA、状态为not_estimable_background_loss，不是nonrecurrent。

## 主阈值0.02的成员分流

旧27中25在共同背景可测，GNPNAT1/KLC1失去跨研究可测资格。matched原M保留这25个全部原共享成员；换人M*或双物种M*后，其中23个保持共享、2个改变方向类别。ZC3H15与YWHAB均是人的复现类别变成nonrecurrent，小鼠仍近端；这不是测量缺失。

两种M*定义的主点共享集均31个且身份完全相同，新增8个ARPP19、ASPH、GLUL、MTDH、OAZ2、RTN3、TPM1、ZMAT2。PDLIM7与EIF4EBP2在三个定义中均保留。主点人近端数为53/62/62，鼠近端为192/192/191（matched原M/人M*/双M*）。

下表只列独立核实的交集计数，不据此作富集推断。

"""+"\n".join(table_lines)+"""

## GEX映射与分工边界

独立从封存gene_group_profiles及library group ID恢复45,705个gene×group，与生产qualified_groups和GEX请求逐行一致。qualifiedMC身份、同库归属、nMC/nCells、每MC APA深度≥20及pooled PAS总数均通过；全gene GEX分母由同一qualifiedMC的旧metacell_metadata独立求和恢复。45,705条小scalar的分母与身份匹配，无缺失；integer gene sum / allgene sum ×10⁶的CPM公式最大误差0，无零gene计数组。

根最新冻结规则采用真实qualifiedMC组CPM，优先于早期资格报告的legacy proxy建议。本分支没有重新载入服务器GEX矩阵、没有独立重算gene numerator，也不声称完成fresh numerator、完整covariates或null验收。原计划30组远程抽查被根明确改派cis的15组独立抽查后取消，未执行重复远程提取；covariates/null/校正亦由cis审查。

## 解释边界与记录

本验收确认改变M定义后近端同源共享成员仍在；31个交集的富集显著性须以另行验收的新null结果判断。旧P–S与P–M人的复现研究构成不同（Xing/Guan vs Xing/Liu/Pan），共同人研究仅Xing，因此不把旧两个对比的保守性差异称为纯粹由参照改变造成。本次固定P与共同gene×library桥接才是原M→M*的公平比较；M及M*仍是研究内参照，不据此称共同培养的未重编程细胞。

比较脚本扩展布尔列时首次运行遇到numpy不支持bool减法的QA实现错误；改为布尔精确相等后全部重跑通过，未更改分析数据或类别规则。早期源metadata读取的mixed dtype警告涉及未用辅助列；后续reader限定所需五列，实际分母/身份检查已通过。阶段QA保留为当时状态，最终结论以FINAL_QA.json为准。所有固定输入结束时复核SHA256未变；本目录manifest列出审查代码、表、日志、报告与QA的SHA256。
"""
    (HERE/"REPORT.zh.md").write_text(report)
    # Do not replace phase snapshots; final comparison stdout is archived separately.
    manifest={"status":"SEALED_SCOPE_PASS","files":{str(p.relative_to(HERE)):sha(p) for p in sorted(HERE.iterdir()) if p.is_file() and p.name!="manifest.json"},"new_null_computed":False,"raw_GEX_reexecution":False}
    (HERE/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print(json.dumps({"status":qa["status"],"frozen_input_files_verified":len(frozen),"primary_flow":flows},indent=2,ensure_ascii=False))


if __name__=="__main__":main()
