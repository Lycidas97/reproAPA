"""Validate the four independently audited branches before sealing the synthesis."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def main():
    checked = []

    def check(condition, label):
        assert condition, label
        checked.append(label)

    original = pd.read_csv(PROJECT / 'notes/2026-10-10-recurrence-utr-ortholog/measured_ortholog_members.tsv.gz', sep='\t')
    base = original[original.contrast.eq('P_M') & original.threshold.eq(.02)].reset_index(drop=True)
    foreground = base.human_direction.eq('proximal') & base.mouse_direction.eq('proximal')
    core = base[foreground]
    check(len(base) == 751 and len(core) == 27, 'common discovery universe')
    check(base.human_gene_id.is_unique and base.mouse_gene_id.is_unique, 'one-to-one unique genes')
    ids = set(core.human_gene_id)
    decay = pd.read_csv(HERE / 'decay/paired_coverage.tsv.gz', sep='\t')
    check(set(decay.human_gene_id) == set(base.human_gene_id), 'decay retains all eligible and ineligible genes')
    check(set(decay[decay.foreground_002].human_gene_id) == ids, 'decay unchanged core')
    qki = pd.read_csv(HERE / 'met-qki/ortholog_external_roster.tsv.gz', sep='\t')
    for comp, part in qki.groupby('comparison'):
        check(set(part.human_gene_id) == set(base.human_gene_id), comp + ' retains all genes')
        check(set(part[part.foreground_002].human_gene_id) == ids, comp + ' unchanged core')
        check(part.loc[~part.external_evaluable, 'effect'].isna().all(), comp + ' missing not zero')
    u1 = pd.read_csv(HERE / 'u1-kinetics/external_probe_inventory.tsv.gz', sep='\t')
    check(set(u1.gene_id) == set(base.human_gene_id), 'U1 all static eligibility records retained')
    check(set(u1[u1.foreground_002].gene_id) == ids, 'U1 unchanged core')
    er = pd.read_csv(HERE / 'membrane/ER_proxy_roster.tsv', sep='\t')
    check(set(er.human_gene_id) == set(base.human_gene_id), 'ER all gene records retained')
    check(set(er[er.foreground_002].human_gene_id) == ids, 'ER unchanged core')
    check(er.loc[~er.eligible, 'ER_region_proxy'].isna().all(), 'ER missing scores not zero')
    check(decay.loc[~decay.paired_eligible].shape[0] > 0, 'sequence coverage explicitly restricted')
    class_pairs = pd.read_csv(HERE / 'core_event_class_pairs.tsv', sep='\t')
    check(set(class_pairs.human_gene_id) == ids, 'class synthesis unchanged core')
    events = pd.read_csv(HERE / 'event_inventory.tsv.gz', sep='\t')
    for sp in ['human', 'mouse']:
        check(set(events[events.species.eq(sp)].gene_id) == set(base[sp + '_gene_id']), sp + ' full event roster')
    contract = json.loads((HERE / 'event_contract.json').read_text())
    check(sha(HERE / 'event_inventory.tsv.gz') == contract['inventory_sha256'], 'event contract hash')
    check(contract['gene_classes_are_structural_compatibility'], 'classification not expressed isoform proof')
    checks = {
        'decay/validation.json': ('status', 'passed'),
        'decay/independent_reference_validation.json': ('status', 'passed'),
        'met-qki/validation.json': ('status', 'PASS'),
        'membrane/independent_validation.json': ('status', 'PASS'),
        'membrane/ER_proxy_independent_validation.json': ('status', 'PASS'),
        'u1-kinetics/verification.json': ('status', 'complete'),
    }
    for name, (key, value) in checks.items():
        info = json.loads((HERE / name).read_text())
        check(info[key] == value, name + ' passed/completed')
        check(info.get('skipped', 0) == 0 and info.get('failed', 0) == 0, name + ' no failed/skipped checks')
    for branch in ['decay', 'met-qki', 'u1-kinetics', 'membrane']:
        check((HERE / branch / 'REPORT.zh.md').is_file(), branch + ' final report exists')
    # The six aggregate tracks are not six biological replicate pairs.
    cov = pd.read_csv(HERE / 'u1-kinetics/external_gene_effects.tsv.gz', sep='\t')
    static = u1[u1.status.eq('static_probe_eligible')]
    check(set(cov.gene_id) == set(static.gene_id), 'U1 covers complete static universe')
    counts = pd.read_csv(HERE / 'u1-kinetics/external_direction_counts.tsv', sep='\t')
    check(set(counts.group) == {'core002', 'core005', 'other', 'all_background'}, 'all U1 foreground/background views')
    check(len(counts) == 3 * 4 * 4, 'all U1 gates/groups/thresholds retained')
    for _, part in counts.groupby(['control_coverage_gate', 'log2_threshold']):
        part = part.set_index('group')
        for column in ['n', 'AMO_both_proximal', 'OE_both_distal', 'reciprocal_expected_direction', 'reciprocal_reverse_direction']:
            check(part.loc['all_background', column] == part.loc['core002', column] + part.loc['other', column],
                  'U1 full background partitions ' + column)
    check((counts.reciprocal_expected_direction <= counts.n).all(), 'U1 reciprocal counts bounded')
    check((HERE / 'membrane/ER_PROXY_INTERPRETATION_CORRECTION.zh.md').is_file(), 'ER library/pseudocount correction retained')
    # Verify sealed branch outputs without relying on prose claims.
    for branch in ['decay', 'met-qki', 'membrane']:
        branch_manifest = json.loads((HERE / branch / 'manifest.json').read_text())
        for relative, digest in branch_manifest['files'].items():
            check(sha(HERE / branch / relative) == digest, branch + ' sealed artifact ' + relative)
    u1_manifest = json.loads((HERE / 'u1-kinetics/manifest.json').read_text())
    for relative, record in u1_manifest['products'].items():
        check(sha(HERE / 'u1-kinetics' / relative) == record['sha256'], 'U1 sealed artifact ' + relative)
    report = {'status': 'PASS', 'checks': len(checked), 'failed': 0, 'skipped': 0,
              'background': 751, 'core': 27, 'validation_checks': checked,
              'validator_sha256': sha(__file__)}
    (HERE / 'bundle_validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'validation_checks'}, indent=2))


if __name__ == '__main__':
    main()
