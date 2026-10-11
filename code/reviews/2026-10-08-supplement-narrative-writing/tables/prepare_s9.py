"""Copy frozen result rows into S9 blocks; no inference or statistical recomputation."""
from pathlib import Path
import csv, gzip, hashlib, json, math, re

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
DATA = HERE / 'data'
DATA.mkdir(parents=True, exist_ok=True)
blocks = []

def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def read(rel):
    path = ROOT / rel
    with (gzip.open(path, 'rt') if str(path).endswith('.gz') else path.open()) as f:
        reader = csv.DictReader(f, delimiter='\t')
        rows = list(reader)
        return reader.fieldnames, rows

def save(name, columns, rows):
    path = DATA / name
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter='\t')
        writer.writeheader(); writer.writerows(rows)
    return str(path.relative_to(ROOT))

def add(sheet, rel, title, notes='', predicate=None, filter_text='All rows; no statistical recomputation.', columns=None):
    fields, original = read(rel)
    rows = [r for r in original if predicate is None or predicate(r)]
    cols = columns or fields
    assert len(set(c.lower() for c in cols)) == len(cols), rel
    assert rows, rel
    blocks.append(dict(sheet=sheet, source=rel, title=title, notes=notes,
                       filter=filter_text + ((' Columns: ' + ', '.join(cols) + '.') if columns else ''),
                       source_rows=len(original), rows=len(rows), columns=cols,
                       values=[[r[c] for c in cols] for r in rows],
                       inputs=[dict(path=rel, sha256=sha(ROOT/rel))]))

def derived(sheet, name, cols, rows, inputs, title, notes, rule):
    rel = save(name, cols, rows)
    add(sheet, rel, title, notes, filter_text=rule)
    blocks[-1]['inputs'] += [dict(path=p, sha256=sha(ROOT/p)) for p in inputs]

RED = 'reviews/2026-10-08-plan-a-redesign/'
OLD = 'reviews/2026-10-08-plan-a-analysis/'
NET = 'reviews/2026-10-08-narrative-logic-chain/netshift-local/'

a1cols = ['test_id','family','pas_source','species','test','observed','candidates','pm_eligible_candidates',
          'background','raw_p','q_BH_family','rejected_BH_05','status_for_correction','permutations','seed',
          'analytic_null_mean','analytic_null_variance','simulated_mean','null_p025','null_p975',
          'attainable_min','attainable_max','tail_probability_ci_low','tail_probability_ci_high',
          'odds_ratio','dual_proximal','dual_distal','complement_proximal','complement_distal','interpretation']
add('S9a', RED+'corrected-tests.tsv', 'Dual-reference tests and correction families',
    'Both conditional null models are reported. The eight conditional tests and four direction-composition tests use separate BH families.',
    lambda r:r['family'].startswith('A1-'), "family begins A1-; retain both PAS workflows and both conditional models.", a1cols)
for stem,title in [('structure-summary.tsv','Conditional background and attainable counts'),
                   ('strata-structure.tsv.gz','Stratum structure'),
                   ('layer-null-diagnostics.tsv.gz','Stratum null moments'),
                   ('candidate-null-probabilities.tsv.gz','Every fixed candidate and its conditional alignment probability'),
                   ('background-input.tsv.gz','Complete measurable background, fixed candidates and directions')]:
    add('S9a', RED+'background/results/'+stem, title)

scorecols = ['test_id','family','species','study_id','route','contrast','score','n_libraries','mean',
             'ci_low','ci_high','raw_p','q_BH_family','status','status_for_correction','rejected_BH_05','null_method']
add('S9b', OLD+'corrected-tests.tsv', 'Raw cell-cycle score tests',
    'Scores use the same frozen expression data for both PAS workflows; these are not duplicated PAS-specific tests.',
    lambda r:r['family']=='A2-cycle-scores', "family == A2-cycle-scores; retain all 102 positions, including NA.", scorecols)
for stem,title in [('cycle_score_library_states.tsv','Raw library-state cell-cycle scores'),
                   ('cycle_score_library_differences.tsv','Raw within-library cell-cycle score differences')]:
    add('S9b', OLD+'cycle/'+stem, title)
cycols = ['test_id','family','pas_source','scope','study_id','route','cohort','statistic','n_libraries',
          'mean','ci_low','ci_high','raw_p','q_BH_family','status','status_for_correction','null_method']
add('S9b', RED+'corrected-tests.tsv', 'Cycle-conditioned dual-reference positions',
    'All 16 formal A2 positions are unavailable (n_libraries = 0, p/q = NA). Zero category counts do not indicate disappearance of APA effects.',
    lambda r:r['family']=='A2-cycle-alignment-redesign',
    "family == A2-cycle-alignment-redesign; exclude A8 route-mean tests.", cycols)
for stem,title in [('context_failure_summary.tsv','Common support and balance failure summary'),
                   ('context_failure_reasons.tsv','Context failure reasons'),
                   ('context_raw_marker_ranges.tsv','Raw marker ranges in each comparison context')]:
    add('S9b', RED+'cycle/'+stem, title,
        'PMS contexts pair a target P/M library with a documented-source S library; source reuse does not create independent biological replicates.')
for caller in ['polyapipe','apadetect_p50']:
    for stem,title in [('context_solvers.tsv','Common-mean feasibility and balance acceptance'),
                       ('context_state_balance.tsv','State-level balance diagnostics'),
                       ('context_pair_balance.tsv','Pairwise marker-balance diagnostics'),
                       ('pm_library_identities.tsv','Physical P–M library identities'),
                       ('dual_species_summary.tsv','Dual-reference candidate measurability by species'),
                       ('dual_route_summary.tsv','Dual-reference candidate measurability by route'),
                       ('dual_species_calls.tsv','All dual-reference candidates including unmeasurable genes'),
                       ('dual_candidate_route_calls.tsv','All route-level dual-reference candidate calls'),
                       ('pm_candidate_summary.tsv','P–M candidate measurability summary'),
                       ('pm_candidate_route_calls.tsv','All P–M candidate calls including unmeasurable genes')]:
        add('S9b', RED+'cycle/'+caller+'/'+stem, title+' ('+caller+')',
            'Continuous weighting and fixed 2×2 descriptive sensitivity are distinct methods. No real-data context passed the prespecified balance gate; no adjusted gene-level effect is available.')

for stem,title in [('corrected-tests.tsv','Dual-reference cross-species overlap and direction tests'),
                   ('polyapipe-shared-members.tsv','Shared dual-reference orthologs (polyApipe)'),
                   ('apadetect-p50-shared-members.tsv','Shared dual-reference orthologs (alternative PAS calls)'),
                   ('direction-exact-pmf.tsv','Exact direction null support'),
                   ('stratum-diagnostics.tsv','Cross-species null strata'),
                   ('power-curves.tsv','Specified alternatives and power bounds'),
                   ('detectable-effects.tsv','Detectable effect summaries')]:
    note = ('Nominal alpha 0.05 is a marginal BH-power upper bound; Bonferroni eight-test power is a conservative lower bound. '
            'conditional_BH_fixed_other_7_observed_p is a conditional diagnostic, not joint BH power.') if 'power' in stem or 'detectable' in stem else ''
    add('S9c', RED+'cross-species/'+stem, title, note)
for caller in ['polyapipe','apadetect-p50']:
    add('S9c', OLD+'cross-species/results/'+caller+'-universe.tsv.gz',
        'Complete jointly measurable one-to-one ortholog background ('+caller+')',
        'The redesign reuses this frozen background and labels; its older test results are not included.')

for stem,title in [('A6-measurable-endpoint-recurrence.tsv','All three profile classes, measurable denominators and missing genes'),
                   ('A6-profile-member-audit.tsv','Every profile candidate and endpoint measurability')]:
    add('S9d', OLD+'bridges/'+stem, title,
        'polyApipe only. Class-entry rates extend the profile-geometry result; they are not independent evidence for unobserved trajectories.')
a6cols=['family','species','test','odds_ratio','rate_difference','raw_p','q_BH_family',
        'family_valid_tests','status','status_for_correction','rejected_BH_05']
add('S9d', OLD+'corrected-tests.tsv', 'All four class-entry comparisons',
    'Four two-sided Fisher tests form one BH family; nonsignificant class comparisons are retained.',
    lambda r:r['family']=='A6-shape-endpoint', "family == A6-shape-endpoint; all four tests.", a6cols)

for caller in ['polyapipe','p50']:
    for stem,title in [('summary.tsv','Gene-wide classification counts'),('gene_classes.tsv','Every strict-interval gene, admission path and gene-wide measurability')]:
        add('S9e', NET+caller+'/'+stem, title+' ('+caller+')',
            'Descriptive selected-pool analysis; no null test. P–S local intervals use observed S whereas gene-wide effects use documented-source S. P–S rows are reference-only. P–M uses the same reference state. Gene-wide unavailable rows remain separate.')

orig = 'notes/2026-10-08-narrative-figures/opposite-sign/opposite-sign.json'
ep = 'notes/2026-10-06-supplement-apadetect-p50/contract-r2/bundle/tables/endpoints/route_effects.tsv.gz'
obj = json.loads((ROOT/orig).read_text())
fields, endpointrows = read(ep)
selected = [r for r in endpointrows if r['setting']=='supported50_p08' and r['contrast']=='P_M'
            and r['design']=='same_library' and r['route'] in ['CIP','WOT']
            and r['gene_id'] in {v['gene_id'] for v in obj.values()}]
assert len({(r['gene_id'],r['route']) for r in selected}) == len(selected)
for symbol,g in obj.items():
    rr=[r for r in selected if r['gene_id']==g['gene_id']]
    g['p50']['dDWUI']={r['route']:float(r['delta_DWUI']) for r in rr}
    g['p50']['setting']='supported50_p08'
    g['p50']['design']='same_library'
    g['p50']['n_pairs']={r['route']:int(r['n_pairs']) for r in rr}
corrected = DATA/'opposite-sign.primary-corrected.json'
corrected.write_text(json.dumps(obj,indent=2)+'\n')
decomp, intervals, avail = [],[],[]
for symbol,g in obj.items():
    for r in g['decomposition']:
        row=dict(gene_symbol=symbol,gene_id=g['gene_id'],pas_source='polyapipe',setting='supported50_p08',
                 design='same_library',dwui_pm_recurrent=g['dwui_pm_recurrent'])
        row.update(r)
        row['family_start0']=row['family'][0];row['family_end0']=row['family'][1];del row['family']
        for k in ['terminal_shift_nt','within','mass','outside']: row.setdefault(k,'NA')
        decomp.append(row)
    for iv in g['strict_intervals']:
        for route in ['CIP','WOT']:
            for cutoff in ['20','50']:
                x=iv[route][cutoff]
                intervals.append(dict(gene_symbol=symbol,gene_id=g['gene_id'],pas_source='polyapipe',
                                      route=route,interval_start0=iv['interval'][0],interval_end0=iv['interval'][1],
                                      interval_length_nt=iv['length'],family_UMI_cutoff=int(cutoff),
                                      n_paired_libraries=x['n'],dI_pp=x['dI_pp']))
    for route in ['CIP','WOT']:
        avail.append(dict(gene_symbol=symbol,gene_id=g['gene_id'],pas_source='apadetect-p50',route=route,
                          setting='supported50_p08',design='same_library',
                          delta_DWUI=g['p50']['dDWUI'].get(route,'NA'),
                          n_pairs=g['p50']['n_pairs'].get(route,'NA'),
                          dwui_pm_recurrent=g['p50']['dwui_pm_recurrent'],interval_pairs=g['p50']['interval_pairs'],
                          strict_pairs=g['p50']['strict_pairs'],
                          local_evaluable='yes' if g['p50']['interval_pairs'] else 'no_interval_pair'))
derived('S9e','opposite-sign-decomposition.tsv',list(decomp[0]),decomp,[orig],
        'Native-axis decomposition of the opposite-sign examples',
        'The three terms sum to the native effect and depend on the family origin; they are not independent causal mechanisms. Unavailable family terms remain NA.',
        'Flatten polyApipe decomposition records from the frozen JSON, retaining all five genes and both routes; no recomputation.')
derived('S9e','opposite-sign-intervals.tsv',list(intervals[0]),intervals,[orig],
        'All strict intervals and library contributors of the opposite-sign examples',
        'Both Fgfr1op2 intervals (970 and 492 nt) and both 20/50-UMI effects are retained.',
        'Flatten strict_intervals across gene, interval, route and UMI cutoff; no recomputation.')
derived('S9e','opposite-sign-p50-primary.tsv',list(avail[0]),avail,[orig,ep],
        'Alternative-caller primary effects and local availability of the examples',
        'Primary supported50_p08 only. The old JSON combined primary and supported50_p08_mc5; these values use an exact primary-setting filter.',
        "Gene IDs from frozen opposite-sign JSON; exact setting == supported50_p08, contrast == P_M, design == same_library, route in CIP/WOT. Unique gene/route; missing routes remain NA.")

titles={'S9a':'Dual-reference conditional null models','S9b':'Cell-cycle scores and common-support diagnostics',
        'S9c':'Dual-reference cross-species overlap, direction and power','S9d':'Profile-class entry into endpoint recurrence',
        'S9e':'Strict intervals and gene-wide shifts'}
analyses={'S9a':'Dual-reference comparison','S9b':'P–M and P/M/source-S cycle diagnostics',
          'S9c':'Cross-species dual-reference comparison','S9d':'Profile-class endpoint comparison',
          'S9e':'Interval and gene-wide comparison'}
order={s:i for i,s in enumerate(titles)}
blocks.sort(key=lambda b:(order[b['sheet']],b['rows']))
for sh in titles:
    for i,b in enumerate([b for b in blocks if b['sheet']==sh],1):b['block']=f'{sh}.{i}'
payload={'baseline':'notes/2026-10-06-supplement-tables/Supplementary_Tables_S1-S8.xlsx',
         'baseline_sha256':sha(ROOT/'notes/2026-10-06-supplement-tables/Supplementary_Tables_S1-S8.xlsx'),
         'titles':titles,'analyses':analyses,'blocks':blocks}
(HERE/'s9-payload.json').write_text(json.dumps(payload,ensure_ascii=False,separators=(',',':'))+'\n')
cols=['sheet','block','title','source','source_rows','rows','filter','sha256','inputs']
with (HERE/'table-sources-S9.tsv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader()
    for b in blocks:
        w.writerow({k:(b[k] if k in b else '') for k in cols if k not in ['sha256','inputs']}|
                   {'sha256':sha(ROOT/b['source']),'inputs':json.dumps(b['inputs'],ensure_ascii=False)})
print(json.dumps({s:{'blocks':sum(b['sheet']==s for b in blocks),'rows':sum(b['rows'] for b in blocks if b['sheet']==s)} for s in titles},indent=2))
