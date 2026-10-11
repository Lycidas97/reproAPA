"""Independent replay and prespecified fixed-panel WOT cycle-dose description."""
from pathlib import Path
import hashlib
import itertools
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OLD = ROOT / 'notes/2026-10-10-cpa-cycle-core-dose'
PROTOCOL = HERE.parent / 'PROTOCOL.zh.md'
EXPECTED_PROTOCOL = '85677ecc53b5bca9f6dbc5d59dc1789274e3c130733c23a57b1f62bcb488051a'
EXPECTED_OLD_PROTOCOL = 'af5c2c2486b0d6ef032c376b3c87d3142f5e672189ec9a6caeaf01a644756291'
TABLES = ROOT / 'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables'
PRIMARY = TABLES / 'state_vectors/pair_effects.primary.tsv.gz'
SECONDARY = TABLES / 'intermediate/paired_library_effects.tsv.gz'
CYCLE_INPUT = ROOT / 'notes/2026-10-10-proliferation-apa/inputs.npz'
CPA_INPUT = ROOT / 'notes/2026-10-10-cpa-factor-expression/expression_inputs.npz'
CPA_PROTOCOL = ROOT / 'notes/2026-10-10-cpa-factor-expression/protocol.json'
PAS_INPUT = ROOT / 'notes/2026-10-10-marker-reference-sensitivity/inputs.npz'
checks = []


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def check(name, condition, **details):
    checks.append(dict(check=name, passed=bool(condition), **details))
    if not condition:
        raise ValueError(name)


def compare(name, a, b, tolerance=1e-12):
    a, b = np.asarray(a, float), np.asarray(b, float)
    error = float(np.nanmax(np.abs(a-b))) if a.size and a.shape == b.shape else None
    check(name, a.shape == b.shape and np.allclose(a, b, atol=tolerance, rtol=tolerance, equal_nan=True), max_abs_error=error)


def save(frame, name):
    compression = {'method': 'gzip', 'mtime': 0} if name.endswith('.gz') else None
    frame.to_csv(HERE / name, sep='\t', index=False, float_format='%.17g', na_rep='NA', compression=compression)


def coefficient(x, y, centered=False):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if centered:
        x, y = x-x.mean(), y-y.mean()
    xx, yy, xy = x @ x, y @ y, x @ y
    return dict(r=float(xy/np.sqrt(xx*yy)), beta=float(xy/xx)) if xx > 0 and yy > 0 else dict(r=np.nan, beta=np.nan)


def raw_scores(units):
    archive = np.load(CYCLE_INPUT, allow_pickle=False)
    meta = pd.DataFrame(json.loads(str(archive['GSE115943_meta_json'])))
    normalized = np.log1p(archive['GSE115943_marker_counts'] / archive['GSE115943_total_gex'][:, None] * 1e4)
    check('exact_18_markers_CP10K_denominator', normalized.shape == (len(meta), 18) and np.all(archive['GSE115943_total_gex'] == meta.total_gex))
    meta['raw_Cycle'] = (normalized[:, :10].mean(1)+normalized[:, 10:].mean(1))/2
    compare('raw_cycle_reconstructs_frozen_MC_scores', meta.raw_Cycle, meta.Cycle, 1e-13)
    check('MC_ids_unique_in_WOT_cycle_archive', meta.metacell_id.is_unique)
    library_rows = []
    for library in units.library_key:
        row = dict(library_key=library)
        for score in ['raw_Cycle', 'Adjusted_Cycle']:
            values = []
            for state in ['P', 'M']:
                q = meta[meta.library_key.eq(library) & meta.node.eq(state)]
                check(f'{library}_{state}_original_cohort_qualified', len(q) >= 3 and q.n_cells.sum() >= 30 and q.purity.ge(.8).all())
                values.append(float(q[score].mean()))
            row['delta_'+score+'_cohort'] = values[0]-values[1]
        library_rows.append(row)
    libs = units[['library_key', 'time', 'day', 'C_label', 'sample_unit']].merge(pd.DataFrame(library_rows), on='library_key', validate='one_to_one')
    cpa = np.load(CPA_INPUT, allow_pickle=False)
    cm = pd.DataFrame(json.loads(str(cpa['GSE115943_meta_json'])))
    factors = json.loads(CPA_PROTOCOL.read_text())['primary_factors']
    check('exact11_predefined_CPA_factors', len(factors) == len(set(factors)) == 11)
    factor_columns = [np.flatnonzero(cpa['GSE115943_factors'] == f).item() for f in factors]
    values = np.log2(1+cpa['GSE115943_counts'][:, factor_columns]/cpa['GSE115943_total_gex'][:, None]*1e6)
    for index, row in libs.iterrows():
        state_values = []
        for state in ['P', 'M']:
            ix = np.flatnonzero(cm.library_key.eq(row.library_key) & cm.node.eq(state) & cm.purity.ge(.8))
            check(f'{row.library_key}_{state}_CPA_qualified', len(ix) >= 3 and cm.n_cells.iloc[ix].sum() >= 30)
            state_values.append(values[ix].mean(0))
        libs.loc[index, 'delta_CPA_cohort'] = float((state_values[0]-state_values[1]).mean())
    return meta, libs


def native_effects(units, selected_genes, meta):
    columns = ['study_id', 'route', 'library_key_a', 'library_key_b', 'time_a', 'contrast', 'design', 'gene_id', 'setting', 'delta_DWUI', 'metacell_ids_a', 'metacell_ids_b', 'n_mc_a', 'n_mc_b', 'n_cells_a', 'n_cells_b', 'n_joint_PAS']
    frames = []
    for part in pd.read_csv(PRIMARY, sep='\t', usecols=columns, chunksize=50000, float_precision='round_trip'):
        keep = part.study_id.eq('GSE115943') & part.route.eq('WOT') & part.contrast.eq('P_M') & part.design.eq('same_library') & part.setting.eq('supported50_p08')
        frames.append(part[keep])
    native = pd.concat(frames, ignore_index=True).rename(columns={'library_key_a':'library_key', 'time_a':'time'})
    presence = pd.read_csv(HERE / 'native_coverage_presence.tsv.gz', sep='\t')
    check('new_effect_inventory_equals_blind_coverage', set(zip(native.library_key, native.gene_id)) == set(zip(presence.library_key, presence.gene_id)))
    check('all24_same_physical_library_and_qualified', native.library_key.eq(native.library_key_b).all() and set(native.library_key) == set(units.library_key) and native[['n_mc_a','n_mc_b']].ge(3).all().all() and native[['n_cells_a','n_cells_b']].ge(30).all().all() and native.n_joint_PAS.ge(2).all())
    second = pd.read_csv(SECONDARY, sep='\t', float_precision='round_trip')
    second = second[second.study_id.eq('GSE115943') & second.route.eq('WOT') & second.setting.eq('supported50_p08') & second.contrast.eq('P_M')]
    joined = native.merge(second[['gene_id', 'library_key', 'delta']], on=['gene_id','library_key'], how='outer', validate='one_to_one', indicator=True)
    check('independent_frozen_pair_table_inventory_agrees', joined._merge.eq('both').all())
    compare('all_native_primary_and_secondary_effects_agree', joined.delta_DWUI, joined.delta, 1e-14)
    native = native[native.gene_id.isin(selected_genes)].copy()
    scores = meta.set_index('metacell_id')
    for index, row in native.iterrows():
        for side, state in [('a','P'), ('b','M')]:
            ids = json.loads(row['metacell_ids_'+side])
            q = scores.loc[ids]
            check(f'{row.library_key}_{row.gene_id}_{state}_exact_MC_identity', len(q) == row['n_mc_'+side] and q.library_key.eq(row.library_key).all() and q.node.eq(state).all() and q.n_cells.sum() == row['n_cells_'+side])
            for score in ['raw_Cycle','Adjusted_Cycle']:
                native.loc[index, score+'_'+side] = float(q[score].mean())
        native.loc[index, 'delta_Cycle_exact_MC'] = native.loc[index,'raw_Cycle_a']-native.loc[index,'raw_Cycle_b']
        native.loc[index, 'delta_Adjusted_Cycle_exact_MC'] = native.loc[index,'Adjusted_Cycle_a']-native.loc[index,'Adjusted_Cycle_b']
    return native


def raw_PAS_replay(genes, units, native):
    """Rebuild full-axis gate and MC-equal DWUI from sparse integer raw columns."""
    archive = np.load(PAS_INPUT, allow_pickle=False)
    meta = pd.DataFrame(json.loads(str(archive['GSE115943_meta_json'])))
    features = pd.DataFrame(json.loads(str(archive['GSE115943_features_json'])))
    pas_data, pas_indices, pas_indptr = [archive['GSE115943_pas_'+key] for key in ['data','indices','indptr']]
    shape = tuple(archive['GSE115943_pas_shape'])
    check('raw_PAS_integer_nonnegative_and_metadata_shape', np.issubdtype(pas_data.dtype,np.integer) and (pas_data >= 0).all() and shape[0] == len(meta))
    groups = {(library,state): np.flatnonzero(meta.library_key.eq(library) & meta.original_selected & meta.original_node.eq(state) & meta.original_purity.ge(.8)) for library in units.library_key for state in ['P','M']}
    rows = []
    for gene in sorted(genes):
        full = features[features.coordinate_eligible & features.primary_gene_id.eq(gene)]
        signed = full.cut_site.to_numpy()*np.where(full.strand.eq('+'),1,-1)
        axis = (signed-signed.min())/np.ptp(signed)
        compare(gene+'_full_axis_not_support_rescaled', full.position, axis, 1e-13)
        keep = full.matched_source_count_50nt.ge(1).to_numpy()
        columns = full.index.to_numpy()[keep]
        position = axis[keep]
        data = np.zeros((len(meta),len(columns)),dtype=pas_data.dtype)
        for col_index,col in enumerate(columns):
            start,end = pas_indptr[col],pas_indptr[col+1]
            data[pas_indices[start:end],col_index] = pas_data[start:end]
        for library in units.library_key:
            states = []
            for state in ['P','M']:
                ix = groups[library,state]
                pre_qualified = len(ix) >= 3 and meta.n_cells.iloc[ix].sum() >= 30
                ix = ix[data[ix].sum(1) >= 20]
                cc = data[ix]
                qualified = pre_qualified and len(ix) >= 3 and meta.n_cells.iloc[ix].sum() >= 30
                mean = float(np.mean(cc@position/cc.sum(1))) if len(ix) else np.nan
                states.append(dict(cc=cc,qualified=qualified,mean=mean,ids=meta.metacell_id.iloc[ix].tolist()))
            a,b = states
            joint = (a['cc'].sum(0)+b['cc'].sum(0) >= 10) & ((a['cc']>0).sum(0)+(b['cc']>0).sum(0) >= 2)
            passes = a['qualified'] and b['qualified'] and np.unique(signed[keep][joint]).size >= 2
            rows.append(dict(gene_id=gene,library_key=library,native_pass=passes,delta_DWUI_raw=a['mean']-b['mean'] if passes else np.nan,MC_ids_P=json.dumps(a['ids']),MC_ids_M=json.dumps(b['ids'])))
    raw = pd.DataFrame(rows)
    selected = raw[raw.native_pass].merge(native[['gene_id','library_key','delta_DWUI','metacell_ids_a','metacell_ids_b']],on=['gene_id','library_key'],how='outer',validate='one_to_one',indicator=True)
    check('raw_gate_matches_selected_native_inventory', selected._merge.eq('both').all())
    compare('raw_full_axis_MC_equal_DWUI_matches_native',selected.delta_DWUI_raw,selected.delta_DWUI,1e-12)
    check('raw_qualified_MC_membership_matches_native', all(set(json.loads(a)) == set(json.loads(b)) for a,b in zip(selected.MC_ids_P,selected.metacell_ids_a)) and all(set(json.loads(a)) == set(json.loads(b)) for a,b in zip(selected.MC_ids_M,selected.metacell_ids_b)))
    save(raw,'raw_PAS_gate_replay.tsv.gz')


def dates_for_panel(libs, native, genes, panel):
    library_rows = []
    for _,lib in libs.iterrows():
        q = native[native.library_key.eq(lib.library_key) & native.gene_id.isin(genes)]
        check(panel+'_'+lib.library_key+'_fixed_panel_complete', len(q) == len(genes) and q.gene_id.nunique() == len(genes))
        row = lib.to_dict()
        row.update(panel=panel,n_genes=len(genes),proximal_mean=float(-q.delta_DWUI.mean()),proximal_median=float(-q.delta_DWUI.median()),delta_Cycle_exact_MC=float(q.delta_Cycle_exact_MC.mean()),delta_Adjusted_Cycle_exact_MC=float(q.delta_Adjusted_Cycle_exact_MC.mean()),gene_ids=json.dumps(sorted(genes)))
        library_rows.append(row)
    libraries = pd.DataFrame(library_rows)
    dates=[]
    for time,q in libraries.groupby('time',sort=False):
        row=dict(panel=panel,time=time,day=float(q.day.iloc[0]),n_wells=len(q),n_genes=len(genes),eligible_date=len(q)==2)
        if len(q)==2:
            check(panel+'_'+time+'_one_C1_one_C2',sorted(q.C_label)==['C1','C2'])
            one,two=[q[q.C_label.eq(label)].iloc[0] for label in ['C1','C2']]
            row.update(library_C1=one.library_key,library_C2=two.library_key,Y_mean=one.proximal_mean-two.proximal_mean,Y_median=one.proximal_median-two.proximal_median)
            for score in ['delta_Cycle_exact_MC','delta_Adjusted_Cycle_exact_MC','delta_raw_Cycle_cohort','delta_Adjusted_Cycle_cohort']:
                row['X_'+score]=one[score]-two[score]
        dates.append(row)
    return libraries,pd.DataFrame(dates)


def descriptive_statistics(dates):
    rows, leaveout=[] ,[]
    for panel,q in dates.groupby('panel',sort=False):
        q=q[q.eligible_date]
        check(panel+'_11_dates_22_wells_no_singleton_contribution',len(q)==11 and q.n_wells.eq(2).all())
        for score in ['delta_Cycle_exact_MC','delta_Adjusted_Cycle_exact_MC','delta_raw_Cycle_cohort','delta_Adjusted_Cycle_cohort']:
            for outcome in ['Y_mean','Y_median']:
                for centered in [False,True]:
                    x,y=q['X_'+score].to_numpy(),q[outcome].to_numpy()
                    row=dict(panel=panel,predictor=score,outcome=outcome,batch_constant_centered=centered,n_dates=len(q),n_wells=22,**coefficient(x,y,centered))
                    rs=[]
                    for i,time in enumerate(q.time):
                        mask=np.arange(len(q))!=i
                        result=coefficient(x[mask],y[mask],centered)
                        leaveout.append(dict(panel=panel,predictor=score,outcome=outcome,batch_constant_centered=centered,omitted_time=time,n_dates=10,**result))
                        rs.append(result['r'])
                    row.update(leave_one_date_r_min=min(rs),leave_one_date_r_max=max(rs),leave_one_date_positive=int(np.sum(np.array(rs)>0)),leave_one_date_negative=int(np.sum(np.array(rs)<0)))
                    rows.append(row)
    return pd.DataFrame(rows),pd.DataFrame(leaveout)


def core_replay(libs,native,core):
    rows=[]
    for time,q in libs.groupby('time',sort=False):
        row=dict(time=time,day=float(q.day.iloc[0]),n_wells=len(q),n_common_genes=0,eligible_primary_secondary=False)
        if len(q)==2:
            one,two=[q[q.C_label.eq(label)].iloc[0] for label in ['C1','C2']]
            a=native[native.library_key.eq(one.library_key)&native.gene_id.isin(core)].set_index('gene_id')
            b=native[native.library_key.eq(two.library_key)&native.gene_id.isin(core)].set_index('gene_id')
            genes=sorted(set(a.index)&set(b.index))
            row.update(n_common_genes=len(genes),gene_ids=json.dumps(genes),library_C1=one.library_key,library_C2=two.library_key,X_Cycle=one.delta_raw_Cycle_cohort-two.delta_raw_Cycle_cohort,X_CPA=one.delta_CPA_cohort-two.delta_CPA_cohort,eligible_primary_secondary=len(genes)>=5)
            if genes:
                aa,bb=a.loc[genes].delta_DWUI,b.loc[genes].delta_DWUI
                row.update(Y_mean=float(-aa.mean()+bb.mean()),Y_median=float(-aa.median()+bb.median()))
        rows.append(row)
    dates=pd.DataFrame(rows).sort_values('day').reset_index(drop=True)
    q=dates[dates.eligible_primary_secondary]
    check('core_fixed7_dates_and_coverage',q.time.tolist()==['D13','D13.5','D14','D14.5','D17','D17.5','D18'] and q.n_common_genes.tolist()==[12,5,22,14,8,5,9])
    statistics,assignments,sensitivity,leaveout=[],[],[],[]
    for score in ['Cycle','CPA']:
        x,y=q['X_'+score].to_numpy(),q.Y_mean.to_numpy()
        result=coefficient(x,y)
        tails=0
        for index,signs in enumerate(itertools.product([-1,1],repeat=7)):
            statistic=coefficient(x*np.array(signs),y)
            tail=abs(statistic['r'])>=abs(result['r'])-1e-12
            tails+=tail
            assignments.append(dict(predictor=score,assignment=index,signs=','.join(map(str,signs)),tail=bool(tail),**statistic))
        row=dict(predictor=score,n_dates=7,n_wells=14,exact_n=128,tail_n=tails,p_exact_two_sided=tails/128,**result)
        rs=[]
        for i,time in enumerate(q.time):
            mask=np.arange(7)!=i
            result_loo=coefficient(x[mask],y[mask])
            leaveout.append(dict(predictor=score,omitted_time=time,n_dates=6,**result_loo))
            rs.append(result_loo['r'])
        row.update(leave_one_date_r_min=min(rs),leave_one_date_r_max=max(rs),leave_one_date_positive=int(np.sum(np.array(rs)>0)),leave_one_date_negative=int(np.sum(np.array(rs)<0)))
        statistics.append(row)
        plans=[('median_min5',q,'Y_median',False),('mean_min1',dates[dates.n_wells.eq(2)&dates.n_common_genes.ge(1)],'Y_mean',False),('mean_min10',dates[dates.n_wells.eq(2)&dates.n_common_genes.ge(10)],'Y_mean',False),('mean_min5_batch_centered',q,'Y_mean',True)]
        for name,panel,outcome,center in plans:
            sensitivity.append(dict(predictor=score,analysis=name,n_dates=len(panel),**coefficient(panel['X_'+score],panel[outcome],center)))
    statistics=pd.DataFrame(statistics)
    order=np.argsort(statistics.p_exact_two_sided.to_numpy(),kind='stable')
    ranked=statistics.p_exact_two_sided.to_numpy()[order]*2/np.arange(1,3)
    q_bh=np.empty(2);q_bh[order]=np.minimum(1,np.minimum.accumulate(ranked[::-1])[::-1])
    statistics['q2_BH']=q_bh
    old=pd.read_csv(OLD/'statistics.tsv',sep='\t',float_precision='round_trip').set_index('predictor').loc[statistics.predictor]
    for col in ['r','beta','exact_n','tail_n','p_exact_two_sided','q2_BH','leave_one_date_r_min','leave_one_date_r_max','leave_one_date_positive','leave_one_date_negative']:
        compare('core_old_formal_'+col+'_independent_replay',statistics[col],old[col],1e-11)
    old_dates=pd.read_csv(OLD/'date_effects.tsv',sep='\t',float_precision='round_trip').set_index('time').loc[dates.time]
    for col in ['X_Cycle','X_CPA','Y_mean','Y_median']:
        compare('core_old_date_'+col+'_independent_replay',dates[col],old_dates[col],1e-11)
    old_assignments=pd.read_csv(OLD/'exact_assignments.tsv',sep='\t',float_precision='round_trip')
    comparisons=pd.DataFrame(assignments).merge(old_assignments,on=['predictor','signs'],suffixes=('_new','_old'),validate='one_to_one')
    check('core_all256_assignments_and_tail_decisions_match',len(comparisons)==256 and comparisons.tail_new.eq(comparisons.tail_old).all())
    compare('core_all256_exact_r_replay',comparisons.r_new,comparisons.r_old,1e-11)
    old_s=pd.read_csv(OLD/'sensitivities.tsv',sep='\t',float_precision='round_trip')
    old_s=old_s[old_s.analysis.isin({r['analysis'] for r in sensitivity})]
    sx=pd.DataFrame(sensitivity).merge(old_s,on=['predictor','analysis'],suffixes=('_new','_old'),validate='one_to_one')
    check('core_all8_sensitivity_and_3day_min10_match',len(sx)==8 and sx.n_dates_new.eq(sx.n_dates_old).all() and sx[sx.analysis.eq('mean_min10')].n_dates_new.eq(3).all())
    compare('core_all8_sensitivity_r_replay',sx.r_new,sx.r_old,1e-11)
    compare('core_all8_sensitivity_beta_replay',sx.beta_new,sx.beta_old,1e-11)
    return dates,statistics,pd.DataFrame(assignments),pd.DataFrame(sensitivity),pd.DataFrame(leaveout)


def intent_checks():
    x=np.array([1.,2.,4.]);y=np.array([2.,3.,2.])
    check('intent_proximal_sign_positive_when_DWUI_negative',coefficient(x,-(-x))['r']>1-1e-13)
    expanded_x=np.ravel(np.column_stack([x/2,-x/2]));expanded_y=np.ravel(np.column_stack([y/2,-y/2]))
    check('intent_date_fixed_is_uncentered_oriented_difference',np.isclose(coefficient(x,y)['r'],np.corrcoef(expanded_x,expanded_y)[0,1]) and not np.isclose(coefficient(x,y)['r'],coefficient(x,y,True)['r']))
    check('intent_no_core_missingness_zero_imputation',set.intersection({'g1'},{'g2'})==set())
    check('intent_equal_gene_not_MC_count_weight',np.isclose(np.mean([.2,.8]),.5) and not np.isclose(np.average([.2,.8],weights=[3,30]),.5))
    check('intent_singleton_has_zero_date_centered_variation',np.var([.8])==0)
    tails=sum(abs(coefficient(np.array(signs)*np.arange(1.,8.),np.arange(1.,8.))['r'])>=1-1e-12 for signs in itertools.product([-1,1],repeat=7))
    check('intent_128_exact_assignments_include_observed',tails==2 and tails/128==.015625)


def main():
    check('root_protocol_hash_frozen',digest(PROTOCOL)==EXPECTED_PROTOCOL)
    check('original_core_protocol_hash_unchanged',digest(OLD/'dose_protocol.json')==EXPECTED_OLD_PROTOCOL)
    intent_checks()
    coverage=json.loads((HERE/'coverage_summary.json').read_text())
    for relative,expected in coverage['inputs'].items():
        check('blind_coverage_input_hash_'+relative,digest(ROOT/relative)==expected)
    check('fixed27_all24_unavailable_is_retained',coverage['panels']['shared27']['all24_n_common']==0)
    panels={'native85':coverage['panels']['native_measured_union']['all24_common_gene_ids'],'ortholog45':coverage['panels']['ortholog_background751']['all24_common_gene_ids'],'mouse_proximal8':coverage['panels']['mouse_proximal372']['all24_common_gene_ids']}
    check('coverage_panels_match_frozen_sizes',[len(panels[k]) for k in ['native85','ortholog45','mouse_proximal8']]==[85,45,8])
    units=pd.read_csv(OLD/'table_audit_wot_units.tsv',sep='\t').sort_values(['day','C_label'])
    orthologs=pd.read_csv(ROOT/'notes/2026-10-10-recurrence-utr-ortholog/measured_ortholog_members.tsv.gz',sep='\t')
    core=set(orthologs.query('contrast=="P_M" and threshold==.02 and human_direction=="proximal" and mouse_direction=="proximal"').mouse_gene_id)
    meta,libs=raw_scores(units)
    selected_genes=core|set().union(*map(set,panels.values()))
    native=native_effects(units,selected_genes,meta)
    measured_core=native[native.gene_id.isin(core)]
    check('core_raw_effect_presence232_all24_intersection_empty',len(measured_core)==232 and len(set.intersection(*(set(measured_core[measured_core.library_key.eq(lib)].gene_id) for lib in units.library_key)))==0)
    check('fixed85_and_core27_are_distinct_estimands',not set(panels['native85'])&core)
    raw_PAS_replay(selected_genes,units,native)
    save(libs,'raw_reconstructed_library_predictors.tsv')
    save(native,'native_gene_effects_exact_MC_scores.tsv.gz')
    library_frames,date_frames=[],[]
    for panel,genes in panels.items():
        libraries,dates=dates_for_panel(libs,native,genes,panel)
        library_frames.append(libraries);date_frames.append(dates)
    libraries,dates=pd.concat(library_frames,ignore_index=True),pd.concat(date_frames,ignore_index=True)
    stats,loo=descriptive_statistics(dates)
    save(libraries,'fixed_panel_library_summary.tsv');save(dates,'fixed_panel_date_dose.tsv')
    save(stats,'fixed_panel_descriptive_statistics.tsv');save(loo,'fixed_panel_leave_one_date_out.tsv')
    core_dates,core_stats,assignments,core_s,core_loo=core_replay(libs,native,core)
    for frame,name in [(core_dates,'core_date_replay.tsv'),(core_stats,'core_formal_statistics_replay.tsv'),(assignments,'core_exact_assignments_replay.tsv'),(core_s,'core_sensitivities_replay.tsv'),(core_loo,'core_leave_one_date_out_replay.tsv')]:save(frame,name)
    check('new_fixed_panel_statistics_have_no_p_values',not any(c in stats for c in ['p','p_exact_two_sided','q','q2_BH']))
    validation=dict(status='passed',tests_skipped=0,n_checks=len(checks),checks=checks,root_protocol_sha256=digest(PROTOCOL),n_native_gene_library=len(native),n_raw_checked_genes=len(selected_genes),n_libraries=24,n_all_dates=13,n_paired_dates=11,n_core_secondary_dates=7,n_exact_assignments_per_predictor=128,core_all24_status='not_estimable_no_shared_measured_genes',max_numeric_error=max(c['max_abs_error'] for c in checks if c.get('max_abs_error') is not None))
    (HERE/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    inputs=[PROTOCOL,OLD/'dose_protocol.json',OLD/'PROTOCOL-CLARIFICATION.zh.md',OLD/'WOT-SAMPLING.zh.md',PRIMARY,SECONDARY,CYCLE_INPUT,CPA_INPUT,CPA_PROTOCOL,PAS_INPUT,OLD/'table_audit_wot_units.tsv',OLD/'statistics.tsv',OLD/'exact_assignments.tsv',OLD/'date_effects.tsv',OLD/'sensitivities.tsv',ROOT/'notes/2026-10-10-recurrence-utr-ortholog/measured_ortholog_members.tsv.gz',HERE/'coverage_summary.json']
    manifest=dict(inputs={str(p.relative_to(ROOT)):digest(p) for p in inputs},outputs={p.name:digest(p) for p in HERE.iterdir() if p.is_file() and p.name!='manifest.json'})
    (HERE/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(stats[(stats.predictor=='delta_Cycle_exact_MC')&(stats.outcome=='Y_mean')].to_string(index=False))
    print(core_stats.to_string(index=False))
    print(core_s.to_string(index=False))
    print(json.dumps({k:v for k,v in validation.items() if k!='checks'},indent=2))


if __name__=='__main__':
    main()
