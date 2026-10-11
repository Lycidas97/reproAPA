"""Spot-check of plotted values against the numbers printed in the analysis REPORTs (independent of the drawing code),
plus the layout-check summary of every new figure. Usage: python3 verify.py"""
import json
from pathlib import Path
O = Path(__file__).resolve().parent/'out'
V = {n: json.loads((O/f'SuppFig.S{n}.values.json').read_text()) for n in (10, 12, 18, 19)}
close = lambda a, b, t=5e-4: abs(a-b) <= t
checks = [
 ('S10a tandem|tandem', V[10]['a tandem|tandem'] == dict(pairs=315, shared_proximal_002=17)),
 ('S10a ALE|ALE', V[10]['a ALE|ALE'] == dict(pairs=52, shared_proximal_002=4)),
 ('S10b t=0.020', V[10]['b t=0.020']['total'] == 27 and V[10]['b t=0.020']['both_tandem'] == 17 and V[10]['b t=0.020']['both_ALE'] == 4),
 ('S10b t=0.050', V[10]['b t=0.050']['total'] == 7 and V[10]['b t=0.050']['strict_full_axis'] == 3),
 ('S10c tandem full assay', V[10]['c both_tandem|full_assay|proximal-proximal']['observed'] == 17 and close(V[10]['c both_tandem|full_assay|proximal-proximal']['expected'], 8.9811)),
 ('S10c tandem q4', close(V[10]['c both_tandem|measurement_PAS_capacity|proximal-proximal']['q'], 0.000160, 1e-6)),
 ('S10c ALE full assay q4', close(V[10]['c both_ALE|full_assay|proximal-proximal']['q'], 0.1782, 1e-3)),
 ('S10d inclusion pooled', close(V[10]['d segment']['pooled_inclusion']['human'][1]-V[10]['d segment']['pooled_inclusion']['human'][0], 0.307851, 1e-6)),
 ('S10e 25 nt', [V[10][f'e 25 {s}'] for s in ('all_contexts', 'strict_UTR', 'order_eligible', 'mass_eligible', 'driver_judgeable', 'double_driver_correspondence', 'two_study_exact_event_both')] == [48, 9, 7, 2, 1, 1, 1]),
 ('S10e 50 nt', [V[10][f'e 50 {s}'] for s in ('all_contexts', 'strict_UTR', 'order_eligible', 'mass_eligible', 'driver_judgeable', 'double_driver_correspondence', 'two_study_exact_event_both')] == [48, 9, 8, 4, 3, 2, 2]),
 ('S12a Liu P7-S5 naive', V[12]['a P_M_recurrent|primed_P7_S5_pool']['naive'][0] == 78 and close(V[12]['a P_M_recurrent|primed_P7_S5_pool']['naive'][1], 55.864, 1e-3)),
 ('S12b Liu P7-S5 q', close(V[12]['b P_M_recurrent|primed_P7_S5_pool']['q14'], 0.387781, 1e-6)),
 ('S12b P-M recurrent dual q', close(V[12]['b P_M_PS_same_contributor_dual|human_P_M_recurrent']['q14'], 0.945303, 1e-6)),
 ('S18a CIP P-M', close(V[18]['a CIP']['P_M'], 0.4389)),
 ('S18a WOT P-M', close(V[18]['a WOT']['P_M'], 0.4668)),
 ('S18c CIP P-S NUDT21', close(V[18]['c P_S|CIP|NUDT21'], 0.36720)),
 ('S18c Guan0618 P-S CPSF6', close(V[18]['c P_S|178325_0618|CPSF6'], 1.43096)),
 ('S18d mouse common q8', close(V[18]['d paired_strict|mouse|proximal_vs_nonrecurrent']['q8'], 0.08340, 1e-4)),
 ('S18e Xing P', close(V[18]['e GSE118258|P']['r'], 0.427, 1e-3)),
 ('S19a WOT', close(V[19]['a WOT']['cycle_S_minus_M'], 0.0343) and close(V[19]['a WOT']['core_mean_dDWUI_S_minus_M'], -0.0437)),
 ('S19a Liu', V[19]['a GSE150311']['proximal'] == 1 and V[19]['a GSE150311']['distal'] == 7),
 ('S19b core subset', close(V[19]['b Core subset']['r'], 0.664, 1e-3) and close(V[19]['b Core subset']['q2'], 0.71875, 1e-5)),
 ('S19b native85', close(V[19]['b 85 fixed genes']['r'], -0.120, 1e-3)),
 ('S19c EMT', close(V[19]['c External targets, 605 genes|EMT, distal']['q'], 0.1167, 1e-3)),
 ('S19d U1 t0.2', V[19]['d U1_reciprocal_log2_t0.2 core']['hits'] == 5 and V[19]['d U1_reciprocal_log2_t0.2 other']['hits'] == 81),
]
bad = [c for c, ok in checks if not ok]
print(f'{len(checks)-len(bad)}/{len(checks)} value checks pass', bad)
for n in (10, 12, 18, 19):
    c = json.loads((O/f'SuppFig.S{n}.receipt.json').read_text())['checks']
    print(f'S{n}', {k: len(c[k]) for k in ('outside', 'text_overlaps', 'letter_zone', 'line_aspect_bad', 'colour_unnamed')}, 'sizes', c['sizes'])
