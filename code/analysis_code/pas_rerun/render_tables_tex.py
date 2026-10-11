"""Build candidate supplementary summaries from verified current-run tables."""
import argparse
import json
from pathlib import Path
import re
import zipfile

import numpy as np
import pandas as pd

from downstream_runtime import Context
from polyapipe_prepare import digest


def escape(value):
    return re.sub(r'[\\&%$#_{}~^]',lambda m:{'\\':r'\textbackslash{}','~':r'\textasciitilde{}','^':r'\textasciicircum{}'}.get(m[0],'\\'+m[0]),str(value))


def display(value):
    if pd.isna(value):return 'NA'
    if isinstance(value,(bool,np.bool_)):return 'Yes' if value else 'No'
    if isinstance(value,(float,int,np.number)):
        if value==0:return '0'
        if isinstance(value,(int,np.integer)):return f'{value:,}'
        if abs(value)<.001:
            mantissa,power=f'{value:.3e}'.split('e')
            return '$'+mantissa+r'\times10^{'+str(int(power))+'}$'
        return f'{value:.4g}'
    return escape(value)


def endpoint_site_union(sites):
    """Count the P--S/P--M union without importing intermediate-state tests."""
    keys=['species','category','gene_id','chrom','strand','start0','end0']
    sites=sites[sites.contrast.isin(['SP','MP'])].drop_duplicates(keys)
    return sites.groupby(['species','category']).agg(
        unique_physical_sites=('gene_id','size'),genes=('gene_id','nunique')).reset_index()


class Tables:
    def __init__(self,ctx):self.ctx=ctx;self.parts=[];self.sources=[];self.blocks=[]

    def read(self,relative):
        path=self.ctx.analysis/relative;node=Path(relative).parts[0]
        receipt=json.loads(self.ctx.register(self.ctx.analysis/node/'result.json').read_text())
        assert receipt['node']==node and receipt['status']=='complete' and all(receipt['checks'].values())
        assert receipt['run_id']==self.ctx.run_id and receipt['input_fingerprint']==self.ctx.config['input_fingerprint']
        known={a['path']:a['sha256'] for a in receipt['artifacts']}
        assert digest(path)==known[str(Path(relative).relative_to(node))]
        self.sources.append(relative);return self.ctx.read(path)

    def table(self,number,caption,label,frame,columns,note,continued=False,rows_per_page=18):
        keys=[k for k,h in columns];missing=set(keys)-set(frame)
        assert not len(frame) or not missing,(number,missing)
        frame=frame.reindex(columns=keys)
        self.blocks.append(dict(table=number,label=label,rows=len(frame),columns=keys,sources=sorted(set(self.sources))))
        self.ctx.table(frame,f'PDF-table-data/S{number}-{label}.tsv')
        for page,start in enumerate(range(0,max(1,len(frame)),rows_per_page)):
            more=continued or page>0;chunk=frame.iloc[start:start+rows_per_page]
            spec='@{}'+''.join('X' if frame[k].dtype.kind not in 'bifu' else 'r' for k in keys)+'@{}'
            text=[r'\clearpage',r'\setcounter{table}{'+str(number-1)+'}' if not more else '',r'\begin{table}[p]']
            if more:text.append(r'\ContinuedFloat')
            text.extend([r'\centering\footnotesize',r'\setlength{\tabcolsep}{3pt}',
                r'\caption'+('[]' if more else '')+'{'+(r'\textit{Cont.} ' if more else '')+caption+'}',
                r'\label{'+label+'}' if not page else '',r'\renewcommand{\tabularxcolumn}[1]{>{\raggedright\arraybackslash}p{#1}}',
                r'\begin{tabularx}{\textwidth}{'+spec+'}',r'\toprule',' & '.join(h for k,h in columns)+r' \\',r'\midrule'])
            text.extend(' & '.join(display(v) for v in row)+r' \\' for row in chunk.itertuples(index=False,name=None))
            if chunk.empty:text.append(r'\multicolumn{'+str(len(keys))+r'}{l}{No eligible records under the stated criteria.} \\')
            text.extend([r'\bottomrule',r'\end{tabularx}',r'\par\smallskip',r'\begin{minipage}{0.98\textwidth}\footnotesize\raggedright',note,
                r'\end{minipage}',r'\end{table}'])
            self.parts.append('\n'.join(text))


def build(ctx,manuscript):
    t=Tables(ctx);original=ctx.register(manuscript/'supplement/supplementary-tables/main.tex').read_text()
    frozen=original.split('The accompanying workbook contains')[0]
    assert frozen.count(r'\caption{')==2
    start=frozen.index('The numbered Tables')
    end=frozen.index(r'\begin{table}',start)
    frozen=frozen[:start]+r'''Tables S1--S7 retain the frozen expression framework. Tables S8--S30 provide the recomputed native polyApipe analyses, their complete denominators and correction families. Historical original-pipeline records are labeled separately and do not enter current inference. The following two tables summarize the unchanged dataset inventory and local expression partitions.

'''+frozen[end:]
    t.parts=[frozen,r'''The accompanying workbook and full-precision TSV collection cover Tables S1--S30. Tables S1--S7 retain the frozen expression-atlas metadata, library inventory, state annotations, integration orthologs and annotation-program features. Tables S8--S30 use the current native polyApipe rerun with the integrated v2 primary reference, except files explicitly marked as historical original-pipeline audits. The v3 stringency-at-least-62 reference is an independently recomputed control.

Table S8 supplies continuous cross-species correspondence, matched-state comparisons and selected-set nulls. Table S9 retains all newly selected shared orthologs, state effects, supporting records and the historical-to-current membership ledger. Tables S10--S13 recompute direction agreement, gene-equal summaries, partial retention and exact-site recurrence on the freshly selected DWUI pool; this is a separate estimand from the union pool in Table S19. Tables S14 and S15 supply current state-family coverage, inclusion profiles, first-crossing labels and family migration. Complete denominators and unavailable measurements remain in the source tables.

The workbook index identifies each part, evidence scope and its lossless TSV. Large gene-level files are supplied as TSV rather than truncated worksheets. NA means unavailable or non-estimable, never zero. Printed summaries round for readability; full-precision values, all predefined test positions and correction-family identifiers accompany the tables.''']
    q=t.read('cross_sharing/sharing_summary.tsv');q=q[q.metric.eq('WPP')&q.test.eq('same_direction')].copy();q['model']='Q5'
    r=t.read('cross_sharing_range10/sharing_summary.tsv');r=r[r.test.eq('same_direction')].copy()
    f=pd.concat([q,r]);names={'GSE118258':'Xing','GSE178325':'Guan','GSE103221':'Guo--CIP','GSE115943':'Schiebinger'}
    f['pair']=f.human_study.map(names)+' / '+f.mouse_study.map(names)
    t.table(16,'Documented-source cross-species DWUI sharing.','tab:opposite-background-sharing',f,
        [('pair','Study pair'),('threshold',r'$\tau$'),('model','Model'),('common_orthologs','$N$'),('observed','$O$'),('null_mean','$E_0$'),('empirical_p','$P$'),('q_primary_8','$q_8$')],
        r'Each model has its own eight-position same-direction family. Q5 matches source-position quintiles; Range10 matches fixed-width 0.1 bins for the source minimum and maximum. Both retain 10,000 permutations and all jointly measurable orthologs. The full TSVs additionally distinguish the primary cross-metric 24-test, broader 56-test and Range10 DWUI-only 24-test families. Counts are responses, not PAS-calling accuracy.')
    definitions=[
      ('Observed-state profiles','Observed S; equal library means','Current profile-selected genes; unchanged full primary-axis position weights'),
      ('Documented-source P--S','Hierarchical source/target sampling-unit and date weights','Current qualified library pairs; one directional vote per study'),
      ('Paired P--M','Same physical library; hierarchical aggregation','Current independently qualified P--M background'),
      ('Source-connected P/M/S','Identical complete contributor records and weights','All three contrasts independently meet joint PAS support; human Xing, mouse CIP and WOT'),
      ('Cross-species study pairs','Documented source; current ortholog background','Q5 primary and Range10 control; distinct correction families'),
      ('Selected proximal sets','Observed S or paired P--M; M-node context proxy','Fresh species-specific memberships; fixed only within each null sensitivity'),
      ('Interval inclusion','Observed S or paired P--M; annotated terminal-family masks','Fresh DWUI pool plus MPRO recurrence; inclusion is distinct from DWUI'),
      ('Historical six counting rules','Original archived workflow','Labeled audit only; excluded from current polyApipe inference')]
    t.table(17,'References, estimators and analysis sets.','tab:analysis-definitions',pd.DataFrame(definitions,columns=['analysis','estimator','scope']),
        [('analysis','Analysis'),('estimator','Reference and weighting'),('scope','Population and interpretation')],
        r'Support is recomputed against the integrated v2 reference within 50 nt on the same chromosome and strand. The four settings vary PAS support, purity and metacell requirements. GEX cells, states, metacell membership and embeddings remain frozen. Definitions in the Supplementary Methods distinguish cells, metacells, physical libraries, capture units, dates and independent biological units.',rows_per_page=4)
    f=t.read('cross_overlap/overlap_null_results.tsv');f=f[f.model.ne('joint_measurement_PAS')]
    t.table(18,'Reference-conditioned selected-set overlap and reference sensitivity.','tab:core-sensitivity-summary',f,
        [('contrast','Contrast'),('model','Model'),('observed','Observed'),('null_median','Null median'),('null_q025',r'2.5\%'),('null_q975',r'97.5\%'),('q_within_model_two_contrasts','$q_2$')],
        r'Each model uses the complete current ortholog universe and its newly selected memberships, with 20,000 independent species-label permutations. M-node reference values are contextual capacity proxies. Historical six-rule consensus is retained separately in explicitly labeled S18 source files and Supplementary Fig.~S9B; it does not enter these current tests.')
    f=t.read('reference_comparison/PAS_support_comparison.tsv');f=f[f.window_nt.eq(50)&f.layer.eq('primary')]
    t.table(18,'Primary native PAS support under the two integrated references.','tab:integrated-reference-support',f,
        [('study_id','Study'),('species','Species'),('denominator','Native PASs'),('supported_v2','v2'),('supported_v3','v3'),('supported_both','Both')],
        r'Denominators are native primary PASs and differ from the reference catalog sizes. Counts at 25 and 100 nt, full-catalog support, independent run identities, effect agreement and membership migration are supplied in S18. Both references use the same native counts and primary axes. These are reference-support comparisons, not sensitivity or specificity against biological truth.',continued=True)
    f=t.read('site_statistics/region_consistency_statistics.tsv');f=f[f.unit.eq('region_pair')&f.category.eq('all')].copy()
    for key in ['cobs','c0','delta_c']:f[key]=100*f[key]
    t.table(19,'Union-pool changed-region direction agreement.','tab:union-pool-regions',f,
        [('species','Species'),('contrast','Contrast'),('genes','Genes'),('cobs',r'$C_{obs}$ (\%)'),('c0',r'$C_0$ (\%)'),('delta_c',r'Excess (pp)'),('gene_label_p',r'$P_g$')],
        r'All matched changed-region pairs enter, including discordant effects. P--M uses common qualified physical libraries; P--S uses observed S. Gene-label magnitude tests require at least three genes and variable gene medians. Smaller groups remain descriptive. Full gene-block intervals, partial retention, exact-site recurrence and counting-unit definitions accompany S19. Physical sites, region pairs and gene--factor relationships are distinct units.')
    f=t.read('site_statistics/recurrence_summary.tsv');f=f[f.scope.eq('observed')]
    t.table(19,'Exact physical-site recurrence in the union pool.','tab:union-pool-exact-sites',f,
        [('species','Species'),('contrast','Contrast'),('category','Class'),('shared_physical_sites','Shared sites'),('both_changed_sites','Changed in both'),('recurrent_sites','Concordant'),('recurrent_genes','Genes')],
        r'Sites require the same genomic interval and strand, sufficient observation in both studies and consistent changes at both coverage thresholds. Physical sites are deduplicated within gene and class. Multiple compatible factors at one physical site do not create additional physical sites. Full site--factor relationships and uncertainty are supplied separately.',continued=True)
    f=endpoint_site_union(t.read('site_statistics/recurrent_sites.tsv.gz'))
    t.table(19,'Physical-site union across the two contrasts.','tab:union-pool-site-union',f,
        [('species','Species'),('category','Class'),('unique_physical_sites','Unique sites'),('genes','Genes')],
        r'Physical sites shared by P--S and P--M are counted once within species, gene and site class. The two classes can overlap at genomic coordinates and must not be summed as a unique all-class inventory.',continued=True)
    members=t.read('matched_pms/metric_membership.tsv.gz');members=members[members.contrast.eq('P_M')&members.DWUI_threshold.eq(.02)]
    vectors=t.read('matched_pms/matched_state_vectors.tsv.gz');v=vectors[vectors.setting.eq('supported50_p08')]
    rows=[]
    for sp in ['human','mouse']:
        routes=['GSE118258'] if sp=='human' else ['CIP','WOT']
        for cohort in ['DWUI_only','shared','MPRO_only']:
            genes=set(members[members.species.eq(sp)&members.cohort.eq(cohort)].gene_id);complete=ps=pm=both=0
            for gene in genes:
                z=v[v.species.eq(sp)&v.gene_id.eq(gene)&v.route.isin(routes)]
                if set(z.route)!=set(routes):continue
                complete+=1;a=z.MPRO_PS.ge(.1).all();b=z.MPRO_PM.ge(.1).all();ps+=int(a);pm+=int(b);both+=int(a and b)
            rows.append(dict(species=sp,cohort=cohort,total=len(genes),complete=complete,unavailable=len(genes)-complete,PS=ps,PM=pm,both=both))
    t.table(20,'Source-connected state summaries in current P--M members.','tab:matched-state-metrics',pd.DataFrame(rows),
        [('species','Species'),('cohort','Membership'),('total','Current'),('complete','Complete'),('unavailable','NA'),('PS','MPRO PS'),('PM','MPRO PM'),('both','Both')],
        r'Membership compares same-direction DWUI 0.02 with unsigned MPRO 0.10 recurrence. Complete matched P/M/S measurements require Xing2020 in human or both mouse studies. MPRO columns count complete genes meeting 0.10 in the indicated contrasts in every required study. Unavailable members remain in the original current-cohort denominator. Shared P contributors make these reference relationships dependent.')
    f=vectors[vectors.gene_symbol.str.upper().eq('NUDT21')]
    t.table(20,r'Matched \textit{Nudt21} effects across all four settings.','tab:matched-nudt21-settings',f,
        [('route','Route'),('setting','Setting'),('DWUI_PS','DWUI PS'),('MPRO_PS','MPRO PS'),('DWUI_PM','DWUI PM'),('MPRO_PM','MPRO PM'),('n_target_libraries','Target libs')],
        r'Each setting requalifies all three contrasts and uses common complete P/M/S contributors. These are descriptive effects; library counts do not establish independent biological replication. Full vectors, weights, source libraries and dates are in S20.',continued=True)
    f=t.read('well_inference/detection_summary.tsv')
    t.table(21,'Joint five-metric independent-well inference.','tab:metric-candidate-differences',f,
        [('metric','Metric'),('genes','Tested genes'),('significant_joint',r'$q_{joint}<0.05$'),('minimum_p','Minimum $P$'),('minimum_joint_q','Minimum joint $q$')],
        r'The family includes every tested gene across DWUI, RWUI, PDUI, PPUI and MPRO. Independent culture wells within dates are the inference units. Detection counts do not compare metric sensitivity; effect thresholds and current recurrence membership are supplied separately.')
    tests=t.read('well_inference/gene_metric_tests.tsv.gz')
    a=tests[tests.metric.eq('DWUI')];b=tests[tests.metric.eq('MPRO')]
    fields=['gene_id','gene_symbol','effect','q_BH_joint','n_units','n_dates']
    f=a[fields].merge(b[['gene_id','effect','q_BH_joint']],on='gene_id',suffixes=('_DWUI','_MPRO'),validate='one_to_one')
    f=f[f.q_BH_joint_DWUI.lt(.05)!=f.q_BH_joint_MPRO.lt(.05)].sort_values('gene_symbol')
    t.table(21,'Complete DWUI/MPRO joint-q discordant gene set.','tab:metric-candidate-significance',f,
        [('gene_symbol','Gene'),('effect_DWUI','DWUI'),('q_BH_joint_DWUI','$q_D$'),('effect_MPRO','MPRO'),('q_BH_joint_MPRO','$q_M$'),('n_units','Wells'),('n_dates','Dates')],
        r'All genes significant by exactly one of the two joint-family criteria are listed, without selecting for effect size. MPRO is unsigned. The source tables also provide simultaneous Monte Carlo error bounds and identify joint-q decisions that remain uncertain at the boundary.',continued=True)
    mc=t.read('well_inference/mc_stability_summary.tsv')
    t.table(21,'Monte Carlo decision stability.','tab:mc-stability',mc,
        [('metric','Metric'),('tests','Tests'),('stable_reject','Stable reject'),('stable_nonreject','Stable nonreject'),('boundary_uncertain','Uncertain')],
        r'Simultaneous Clopper--Pearson bounds cover the simulated tail probabilities across the full Monte Carlo family. Exact sign-enumeration probabilities remain fixed. Applying BH to lower and upper bounds classifies joint-q decisions without changing the primary probabilities.',continued=True)
    f=t.read('measurement/nulls/recurrence_summary.tsv');f=f[f.metric.eq('WPP')&f.threshold.eq(.02)&f.series_matching]
    t.table(22,'DWUI 0.02 recurrence under measurement controls.','tab:measurement-controls',f,
        [('variant','Setting'),('species','Species'),('contrast','Contrast'),('observed','Observed'),('null_mean','Null mean'),('null_p025','2.5 percentile'),('null_p975','97.5 percentile'),('status','Status')],
        r'Each setting rebuilds eligibility and its own conditional null. Unavailable configurations remain explicit and enter the stated full family conservatively. Null ranges are permutation ranges, not biological confidence intervals. The full files include all metrics, thresholds, thinning draws, member retention and same-state controls.')
    f=t.read('motif_background/density_summary.tsv')
    t.table(23,'Sequence-matched direct motif-density differences.','tab:sequence-background-external',f,
        [('species','Species'),('contrast','Contrast'),('category','Class'),('n_genes','Genes'),('mean_difference','Sites/kb'),('lo','Lower'),('hi','Upper'),('q_8_tests','$q_8$')],
        r'Changed and stable windows are matched within genes for length, GC and relative position. Gene-equal differences use 5,000 gene-bootstrap draws and 10,000 sign permutations. Fewer than ten genes remains descriptive with P=1 in the eight-position family; matching attrition remains explicit. Compatibility does not imply occupancy.')
    f=t.read('external_pas/vector_summary.tsv');f=f[f.statistic.eq('tvd')|f.statistic.isna()]
    for field in ['retention','design_gate_pass']:
        if field not in f:
            assert f.status.eq('no_complete_matched_sets').all()
            f[field]=np.nan
    t.table(23,'External iPSC--MEF usage on compatible shared PASs.','tab:external-pas-current',f,
        [('configuration','Design'),('candidate_n','Current candidates'),('measurable_n','Measurable'),('matched_n','Matched'),('retention','Retention'),('estimate','TVD excess'),('design_gate_pass','Gate'),('status','Status')],
        r'The candidate roster is rebuilt from current mouse P--S recurrence; the former counting-rule consensus is not used as a current roster. Controls are globally nonreused. Two external libraries per endpoint permit six whole-library label assignments, so these diagnostics do not establish significant validation. Complete missingness and weight sensitivity are supplied in S23.',continued=True)
    exact=pd.concat([t.read('functional/go_exact.tsv.gz'),t.read('functional/kegg_exact.tsv.gz')])
    f=exact[exact.q_comparison.lt(.05)].copy();f['comparison_label']=f.species+' '+f.contrast+' '+f.comparison+' '+f.target
    f['term']=f.term_name+' ('+f.term_id+')'
    t.table(24,'Exact conditional functional associations.','tab:exact-functional',f,
        [('comparison_label','Comparison'),('term','Term'),('hits','Hits'),('expected','Expected'),('q_comparison','$q_{comparison}$'),('q_database_model','$q_{database}$')],
        r'All term positions passing within-comparison q<0.05 under complete conditioning are shown. G0 is expression-measurable, G1 APA-testable and G2 recurrent; their nested memberships and attrition are rebuilt. Exact tails convolve stratum-specific hypergeometric distributions. Each database retains 28 analyses; terms recurring in different analyses remain distinct tests. Annotation names do not establish pathway activity.',rows_per_page=10)
    f=t.read('cross_controls/cross_species_ladder.tsv')
    t.table(25,'Cross-species sharing under nested conditioning.','tab:null-controls',f,
        [('contrast','Contrast'),('model','Conditions'),('observed','Observed'),('expected','Expected'),('p','$P$'),('q_two_contrasts','$q_2$')],
        r'All models retain the newly selected memberships and complete current ortholog universe. Independent species-label shuffles use 20,000 draws; q2 adjusts both contrasts per model. Measurement-only conditioning is diagnostic. Conditioning, calibration and signal-injection power are supplied separately; these sensitivity analyses do not replace the primary tests.')
    candidates=t.read('rbp_expression/candidates.tsv');selected=candidates[candidates.display_selected]
    f=t.read('rbp_expression/cross_study_summary.tsv').merge(selected[['species','factor','target_genes']],on=['species','factor'],how='inner',validate='many_to_one')
    t.table(26,'Current candidate-RBP expression context.','tab:rbp-expression',f,
        [('species','Species'),('factor','Factor'),('target_genes','Targets'),('contrast','Contrast'),('effect_a',r'$\Delta_A$'),('effect_b',r'$\Delta_B$'),('direction','Direction')],
        r'All current mouse direct-motif candidates and the top six human candidates by recurrent target-gene count are shown. Factors are selected before expression is examined. The full candidate set, including transferred motifs and low-expression factors, remains in S26. Effects are differences in mean log2(CPM+1); P--M uses common physical libraries. Study identities and weights are supplied in each source record; these are descriptive effects.')
    f=t.read('longread/structural_summary.tsv').merge(t.read('longread/coordinate_liftover.tsv')[
        ['feature_id','chrom','strand','start0','end0','factor','contrast']],on='feature_id',validate='one_to_one');f['feature']=np.arange(1,len(f)+1)
    ctx.table(f,'table/S27-feature-key.tsv')
    unique=f.drop_duplicates(['gene','kind','chrom','strand','start0','end0'])
    summary=unique.assign(mapped=unique.mapping_status.eq('unique_gapfree_reciprocal'),
        strict=unique.strict_pairs.gt(0),boundary=unique.boundary1nt_pairs.gt(0)).groupby(['gene','kind']).agg(
            queried=('feature_id','size'),mapped=('mapped','sum'),strict=('strict','sum'),boundary=('boundary','sum')).reset_index()
    t.table(27,'Long-read compatibility with current candidate features.','tab:external-long-read',summary,
        [('gene','Gene'),('kind','Kind'),('queried','Queries'),('mapped','Mapped'),('strict','Strict pair'),('boundary','Boundary pair')],
        r'Counts are distinct physical query features, deduplicated across factors and contrasts within gene and feature kind. Strict/boundary columns count features with at least one compatible pair, not transcript-pair totals. Strict pairs share the complete intron chain and terminal boundary, with inclusion in one isoform and full exclusion in the other; one-base boundary overlaps are separate. All query identities, unsuccessful mappings, transcript pairs and sample columns remain in S27 and its feature-key TSV. Archived models were queried without reprocessing raw reads.')
    tracking=t.read('longread/historical_three_gene_tracking.tsv')
    t.table(27,'Coverage of all three prespecified long-read loci.','tab:external-three-locus-coverage',tracking,
        [('gene','Gene'),('status','Current eligibility'),('new_features','Query records'),('filtered_isoforms','Archived models')],
        r'The three historical gene identities are retained even when no current recurrent feature is available. Query records distinguish contrasts and factors and therefore differ from the physical-feature counts above. TPM columns are sequencing samples, not inferred independent biological replicates.',continued=True)
    f=t.read('ago2/overlap_summary.tsv').copy();ids={v:i+1 for i,v in enumerate(f.feature_id.drop_duplicates())};f['feature']=f.feature_id.map(ids)
    ctx.table(pd.DataFrame(list(ids.items()),columns=['feature_id','feature']),'table/S28-feature-key.tsv')
    t.table(28,'Current Rad23b query overlap in the AGO2 archive.','tab:rad23b-ago2',f,
        [('feature','Feature'),('label','Archived track'),('biological_unit','Biological unit'),('overlap_clusters','Clusters'),('overlap_bp_sum','Shared bp')],
        r'Queries use current intervals and seed sites, with unique reciprocal conversion through mm39--mm10--mm9 and same-chromosome/strand overlap. The S28 feature-key TSV, full sequences and absent-query eligibility accompany this table. WT1A and WT1B are technical libraries of one preparation; merged tracks do not retain replicate membership. Cluster boundaries do not measure nucleotide-resolved binding.')
    f=t.read('tables/table/S29-source-data-index.tsv');index=f.groupby(['node','evidence_scope']).agg(files=('source','size'),bytes=('bytes','sum')).reset_index()
    t.table(29,'Complete numerical source-data inventory.','tab:supp-revision-data',index,
        [('node','Analysis'),('evidence_scope','Evidence scope'),('files','Files'),('bytes','Bytes')],
        r'The archive contains full-precision tables, null distributions, plotted numerical arrays, coordinates and audit records from this run. Each record has its path, SHA-256 and evidence scope. Historical original-pipeline records are identified explicitly and cannot be interpreted as current polyApipe results. Frozen GEX resources retain their original memberships. Exact table and panel mappings are supplied with the workbook index.')
    f=t.read('caller_comparison/caller_availability.tsv')
    t.table(30,'PAS-caller comparison coverage.','tab:pas-caller-comparison',f,
        [('study_id','Study'),('caller','Caller'),('status','Status'),('native_result','Current result')],
        r'This batch compares native polyApipe with the original PAS pipeline on explicit native and common-measurement axes. Complete bidirectional PAS denominators, effect agreement, original recurrent-gene migration, named examples and missingness accompany S30. SCAPE is reserved for its separate comparison; availability of an intermediate SCAPE file does not imply a completed seven-study benchmark.')
    return t


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);p.add_argument('--manuscript',type=Path,required=True)
    a=p.parse_args();ctx=Context(a.config,'text',a.out);t=build(ctx,a.manuscript)
    target=ctx.out/'split_sources/supplement/supplementary-tables/main.tex'
    target.parent.mkdir(parents=True,exist_ok=True);target.write_text('\n\n'.join(t.parts)+'\n\\clearpage\n')
    assert {b['table'] for b in t.blocks}==set(range(16,31))
    expected=set(re.findall(r'\\label\{([^}]+)\}',(a.manuscript/'supplement/supplementary-tables/main.tex').read_text()))
    assert expected<=set(re.findall(r'\\label\{([^}]+)\}',target.read_text()))
    (ctx.out/'table-text-review.json').write_text(json.dumps(dict(status='generated; compile_and_visual_review_pending',
        run_id=ctx.run_id,input_fingerprint=ctx.config['input_fingerprint'],source=str(target),sha256=digest(target),
        blocks=t.blocks,input_sha256=ctx.inputs),indent=2)+'\n')
    archive=ctx.out/'table/S29-PDF-summary-source-data.zip'
    archive.parent.mkdir(exist_ok=True)
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=3) as bundle:
        for source in sorted((ctx.out/'PDF-table-data').glob('*.tsv')):bundle.write(source,str(source.relative_to(ctx.out)))
        for source in sorted((ctx.out/'table').glob('*feature-key.tsv')):bundle.write(source,source.name)
        bundle.write(Path(__file__),'render_tables_tex.py')
        bundle.write(ctx.out/'table-text-review.json','table-text-review.json')
    print(target)
