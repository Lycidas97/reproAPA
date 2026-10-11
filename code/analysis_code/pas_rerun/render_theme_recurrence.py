from render_theme import *
DRAW={}
def D09(r):
    r.panel('A');x=r.read('recurrence_controls/recurrence_ladder.tsv');mods=['L0_feasible_series','L1_measurement','L2_baseline','L2_PAS_first','L3_full'];labs=[context(sp,co)+' '+['L0','L1','L2','L2-PAS','L3'][j] for sp,co in CTX for j in range(5)];q=x.copy();q['_row']=[CTX.index((a.species,a.contrast))*5+mods.index(a.model) for a in q.itertuples()];q['low']=(q.expected-q.null_sd)/q.expected;q['high']=(q.expected+q.null_sd)/q.expected;q['ratio']=q.observed/q.expected;ax=r.dots(0,q,'ratio',labs,lo='low',hi='high',series='threshold',xlabel='Observed / expected',reference=1,q='q_eight');r.note('Markers observed/expected. Bands show null mean ± 1 null SD divided by null mean; these are not confidence intervals. Forty settings, original eight-test q per model; all joint-40 q retained in data.')
    m=r.read('recurrence/recurrence_members.tsv.gz');m=m[m.metric.eq('WPP')&m.threshold.eq(.02)&m.series_matching&m.variant.eq('primary')];f=m.groupby(['species','contrast','direction']).size().unstack(fill_value=0).reindex(CTX);f.index=CL;f=f.rename(columns={'negative':'Proximal','positive':'Distal'});r.bars(1,f,'Recurrent genes',colorset=[BLUE,ORANGE]);r.note('Distal/proximal counts are shown separately for each cohort; their split is not a study vote.')
    x=r.read('measurement/same-state/recurrence_draws.tsv.gz');x=x[x.metric.eq('WPP')&x.threshold.eq(.02)];groups=[]
    for (sp,co),g in x.groupby(['species','contrast']):groups.append((context(sp,co),g.pseudo_recurrent))
    r.ecdf(2,groups,'Pseudo-recurrent genes');r.record('all-same-state-eligible-backgrounds',x);r.note('Pseudo recurrence distribution pools the two labeled reference/target diagnostic splits only for display; all replicates and eligible/observed backgrounds remain separate in data.')
    r.panel('B');q=r.read('measurement/nulls/member_retention.tsv');q=q[q.metric.eq('WPP')&q.threshold.eq(.02)&q.contrast.eq('P_S')];q=q[q.variant.eq('interval_supported')];r.record('bootstrap-source',q)
    if len(q)>2:q=q.drop_duplicates(['species','same_direction_intersection','original'])
    f=pd.DataFrame({'Retained':q.same_direction_intersection,'Other':q.original-q.same_direction_intersection});f.index=q.species.str.title();r.bars(0,f,'Original recurrent genes',colorset=[BLUE,GREY])
    x=r.read('recurrence/leave_one_study_out_nulls.tsv');x=x[x.metric.eq('WPP')&x.species.eq('human')&x.contrast.eq('P_M')&x.threshold.eq(.02)&x.series_matching].copy();x['_row']=range(len(x));x['low']=x.null_p025/x.null_mean;x['high']=x.null_p975/x.null_mean;r.dots(1,x,'enrichment',['Omit '+short(a) for a in x.omitted_study],lo='low',hi='high',xlabel='Human P–M observed / expected',reference=1)
    g=r.read('recurrence_controls/guan_series_unit_sensitivity.tsv');g=g[g.metric.eq('WPP')];ax=r.ax(2,kind='line',pad=(12,3,1,10))
    for j,(co,z) in enumerate(g.groupby('contrast')):
        if z.enrichment.notna().any():ax.plot(z.threshold,z.enrichment,'o-',color=PAL[j],ms=3,label='H '+co.replace('_','–'))
    ax.text(.05,.13,'P–M NA',transform=ax.transAxes,fontsize=7,color=MUTED)
    ax.axhline(1,color=GREY,ls=':');ax.set(xlabel='DWUI threshold',ylabel='Pseudounit O/E');legend(ax,loc='best');r.record('Guan-series-sensitivity',g);r.note('Guan series are deliberately treated as pseudounits in this sensitivity; they are not independent studies. The plotted ratio does not replace the study-level result.');r.note('Study, series and library levels are not interchangeable replication units.')
def D10(r):
    r.panel('A');q=r.read('well_inference/gene_metric_tests.tsv.gz');z=q.pivot(index='gene_id',columns='metric',values='q_BH_joint');ax=r.ax(0,kind='scatter');sig=(z.DWUI<.05)&(z.MPRO<.05);ax.scatter(-np.log10(z.DWUI),-np.log10(z.MPRO),s=8,c=np.where(sig,BLUE,GREY),lw=0,alpha=.7);ax.axvline(-np.log10(.05),color=GREY,ls=':');ax.axhline(-np.log10(.05),color=GREY,ls=':');ax.set(xlabel='−log10 joint q · DWUI',ylabel='−log10 joint q · MPRO');legend(ax,handles=[Line2D([],[],ls='',marker='o',ms=3,c=BLUE,label='Both q < .05'),Line2D([],[],ls='',marker='o',ms=3,c=GREY,label='Other')],ncol=2);r.record('all-joint-q',z)
    m=q.groupby('metric').joint_q05.agg(['sum','count']).reindex(['DWUI','RWUI','MPRO','PDUI','PPUI']);f=pd.DataFrame({'q < .05':m['sum'],'Other':m['count']-m['sum']});r.bars(1,f,'WOT genes',colorset=[BLUE,GREY]);r.note(f'{q.gene_id.nunique()} genes, {len(q)} metric tests; q_BH_joint retains the producer correction family. Blue scatter points pass both plotted metrics; non-estimable tests are not zeros.')
    r.panel('B');x=r.read('well_inference/source_endpoint_WPP_tests.tsv.gz');d=r.read('axis:well_six_unit_downsample_summary.tsv');groups=[('Source fraction · count',x[x.metric.eq('count_WPP')].source_variance_fraction),('Source fraction · equal MC',x[x.metric.eq('equal_MC_WPP')].source_variance_fraction),('Six-well minimum q',d.min_joint_q)];
    for j,gg in enumerate([groups[:2],groups[2:]]):
        ax=r.ax(0,sub=j,n=2,kind='distribution',pad=(12 if j==0 else 2,3,2,10))
        for k,(lab,values) in enumerate(gg):
            vv=np.sort(np.asarray(values,float));vv=vv[np.isfinite(vv)];ax.step(vv,np.arange(1,len(vv)+1)/len(vv),where='post',c=PAL[k],lw=.9,label=('Count' if k==0 else 'Equal MC') if j==0 else '100 draws')
        ax.set(xlabel='Source variance fraction' if j==0 else 'Minimum joint q',ylabel='Cumulative fraction',ylim=(0,1),yticks=[0,.5,1]);legend(ax,loc='best')
    r.record('100-six-unit-draws',d);r.record('source-statuses',x)
    p=r.read('axis:pool_resolution.tsv');ax=r.ax(1,kind='line');ax.plot(p.n_units,p.exact_two_sided_sign_min_p,'o-',color=BLUE,ms=3);ax.axhline(.05,color=GREY,ls=':');ax.set(xlabel='Independent sign units',ylabel='Minimum two-sided exact P',xticks=[4,5,6]);r.record('finite-assignment-resolution',p);r.note('Six-well draws are a resolution diagnostic; the source table retains significant-gene counts for every draw. Source variance fractions and minimum q use separate x axes and a shared cumulative-probability y axis.')
    r.panel('C');r.note('External design gate:4 physical libraries,6 whole-library assignments; primary design fails. Definition moved from a standalone schematic into this caption.');x=r.read('external_pas/vector_summary.tsv');q=x[x.statistic.eq('tvd')].copy();q['_row']=range(len(q));ax=r.dots(1,q,'estimate',q.configuration.tolist(),lo='sensitivity_min',hi='sensitivity_max',xlabel='External PAS TVD contrast');ax.set_yticklabels([a.configuration+(' +' if a.design_gate_pass else ' ×') for a in q.itertuples()]);ax.text(0,1.04,'Matching gate\n6 assignments',transform=ax.transAxes,fontsize=7);r.record('all-external-diagnostics',x);r.note('Row suffix + passes the matching gate; × fails. Every configuration has only six whole-library assignments. Assignment ranges are descriptive. controls2 passes matching qualification but still has only six assignments; no external validation significance is claimed.');r.record('date-conditioned-summary',r.read('axis:date_conditioned_summary.tsv'))
def D11(r):
    from render_theme_panels import pas_curves
    r.panel('A');q=r.read('base:validation/narrative-v4-review/Figure3C-counts.tsv');ax=r.ax(0,kind='line',pad=(13,4,3,10))
    for i,a in enumerate(q.itertuples()):
        vals=[a.candidate_genes,a.measurable_genes,a.dual_genes];ax.plot(range(3),vals,'o-',color=PAL[i],ms=3,label='Xing' if a.species=='human' else 'CIP/WOT')
        for j,v in enumerate(vals):ax.annotate(str(v),(j,v),xytext=(0,4 if i else -10),textcoords='offset points',ha='center',fontsize=7)
    ax.set(xticks=[0,1,2],xticklabels=['Candidate','Measured','Dual'],ylabel='Genes',ylim=(0,q.candidate_genes.max()*1.15),xlim=(-.15,2.2));legend(ax,loc='best');r.note('Counts are genes, not gene×route records. Human uses source-qualified Xing; mouse requires both routes.')
    r.record('dual-qualification-flow',q);m=r.read('recurrence/recurrence_members.tsv.gz');m=m[m.series_matching&m.variant.eq('primary')];rows=[]
    for sp,co in CTX:
        z=m[m.species.eq(sp)&m.contrast.eq(co)];a=set(z.loc[z.metric.eq('WPP')&z.threshold.eq(.02),'gene_id']);b=set(z.loc[z.metric.eq('MPRO')&z.threshold.eq(.1),'gene_id']);rows.append([len(a-b),len(a&b),len(b-a)])
    r.bars(1,pd.DataFrame(rows,index=CL,columns=['DWUI','Both','MPRO']),'Unique genes');r.note('Set memberships are DWUI-only, both and MPRO-only. Dual-reference qualification is a separate within-gene criterion.')
    r.panel('B')
    for j,sp in enumerate(['human','mouse']):
        x=r.read(f'base:validation/narrative-v4-review/Figure3C-{sp}-members.tsv');ax=r.ax(0,sub=j,n=2,kind='scatter')
        for k,(route,g) in enumerate(x.groupby('route')):
            col=dataset_color(short(route),PAL[k])
            for qualified in [False,True]:
                z=g[g.dual_qualified_gene.eq(qualified)];ax.scatter(z.DWUI_PS,z.DWUI_PM,facecolor=col if qualified else 'none',edgecolor=col,marker=['o','^'][k],s=6,lw=.35,alpha=.65,rasterized=True)

        zeros(ax,True);ax.set(xlabel=('Xing' if sp=='human' else 'CIP/WOT')+' P–S ΔDWUI',ylabel='Matched P–M ΔDWUI');legend(ax,handles=[Line2D([],[],ls='',marker='o',ms=3,mfc='none',mec=INK,label='Other measured'),Line2D([],[],ls='',marker='o',ms=3,c=INK,label='Dual qualified')],ncol=2);ax.text(.02,1.04,f'{sp.title()} · {x.gene_id.nunique()} genes',transform=ax.transAxes,fontsize=8);r.record('complete-dual-'+sp,x)
    r.panel('C');x=r.read('matched_pms/matched_state_vectors.tsv.gz');x=x[x.setting.eq('supported50_p08')&x.gene_symbol.eq('Pdlim7')&x.route.isin(['CIP','WOT'])]
    for j,ro in enumerate(['CIP','WOT']):pas_curves(r,0,x[x.route.eq(ro)],'Pdlim7',j,2,matched=True)
    rows=[]
    for i,ro in enumerate(['CIP','WOT']):
        a=x[x.route.eq(ro)].iloc[0]
        for j,co in enumerate(['PS','PM','MS']):rows.append(dict(_row=i*3+j,value=a['DWUI_'+co]))
    ax=r.dots(1,pd.DataFrame(rows),'value',[ro+' '+co for ro in ['CIP','WOT'] for co in ['PS','PM','MS']],xlabel='ΔDWUI',pad=(15,5,2,11));ax.text(0,1.04,'Pdlim7',transform=ax.transAxes,fontstyle='italic',fontsize=8);r.note('All panels use matched common-contributor P/M/S vectors. Step curves retain all PAS usage fractions.')
def D13(r):
    r.panel('A');allq=r.read('intermediate/process/coverage_flow.tsv');q=allq[allq.variant.eq('shared_p08')&allq.scope.eq('requalified')&allq.endpoint_threshold.eq(.02)&allq.intermediate_threshold.eq(.02)];ax=r.ax(0,kind='structure');ax.set(xlim=(0,1),ylim=(0,1))
    from matplotlib.transforms import ScaledTranslation
    labels=['Candidate','Measured','Retained','R','E'];xx=np.array([0,.25,.5,.75,1.0]);ax.set(xlim=(-.09,1.10),ylim=(-.5,2.15));ax._row_alignable=True
    for j,label in enumerate(labels):ax.text(xx[j],2.05,label,ha='center',va='bottom',fontsize=7)
    for i,a in enumerate(q.itertuples()):
        yy=1-i;vals=[a.original_candidates,a.endpoint_measurable,a.endpoint_retained,a.R,a.E]
        ax.text(-.09,yy,a.species[0].upper(),ha='right',va='center',fontsize=8)
        ax.plot(xx[:3],[yy]*3,c=PAL[i],lw=.8)
        for k in [3,4]:ax.plot([xx[2],xx[k]],[yy,yy+.17 if k==3 else yy-.17],c=PAL[i],lw=.6)
        for j,v in enumerate(vals):
            y=yy if j<3 else yy+(.17 if j==3 else -.17)
            ax.scatter(xx[j],y,s=18,marker='o' if np.isfinite(v) else 'x',c=PAL[i] if np.isfinite(v) else GREY,zorder=3)
            ax.text(xx[j],y,str(int(v)) if np.isfinite(v) else 'NA',ha='center',va='bottom' if j==3 else 'top',fontsize=7,transform=ax.transData+ScaledTranslation(0,(4 if j==3 else -4)/72,ax.figure.dpi_scale_trans))
    r.record('all-process-coverage',allq);r.note('Flow: original recurrent → endpoint measured → endpoint retained. R and E are parallel branches, not successive attrition.')
    x=r.read('intermediate/process/frozen_expression_support_summary.tsv');f=x.groupby('node')[['support_q90','support_q95']].sum().reindex(list('SREP'));ax=r.ax(1,labels=f.index.tolist(),kind='bar')
    for j,c in enumerate(f):ax.barh(np.arange(4)+(j-.5)*.32,f[c],height=.28,color=PAL[j],label=c[-3:])
    ax.set(yticks=range(4),yticklabels=f.index,xlabel='Supported metacells');legend(ax,loc='best');r.record('frozen-expression-support',x)
    r.panel('B');allx=r.read('intermediate/process/same_series_process_classes.tsv.gz');x=allx[allx.variant.eq('shared_p08')&allx.scope.eq('requalified')&allx.intermediate_threshold.eq(.02)&allx.endpoint_retained&allx.pair_measurable];cats=['endpoint_aligned','reverse','within_threshold','direction_difference','mixed_or_threshold_undetermined'];names=['Aligned','Reverse','Small','Discordant','Mixed']
    for j,sp in enumerate(['human','mouse']):
        ax=r.ax(0,sub=j,n=2,kind='scatter',pad=(12,6,1,11));z=x[x.species.eq(sp)]
        for state,mk in [('R','o'),('E','^')]:
            for k,cat in enumerate(cats):
                g=z[z.state.eq(state)&z.pair_class.eq(cat)];ax.scatter(g.aligned_delta_1,g.aligned_delta_2,s=4,marker=mk,c=PAL[k],alpha=.5,lw=0,rasterized=True)
        zeros(ax,True);assert z.route_1.nunique()==z.route_2.nunique()==1;ax.set(xlabel=short(z.route_1.iloc[0])+' aligned Δ',ylabel=short(z.route_2.iloc[0])+' aligned Δ');ax.text(.02,1.05,sp.title(),transform=ax.transAxes,fontsize=7);legend(ax,handles=[Line2D([],[],ls='',marker='o',ms=3,c=PAL[k],label=name) for k,name in enumerate(names)]+[Line2D([],[],ls='',marker=mk,ms=3,c=INK,label=st) for st,mk in [('R','o'),('E','^')]],ncol=3);r.record('complete-process-'+sp,z)
    r.share_y(r.panel_axes[r.active][:2],'Aligned intermediate ΔDWUI in each named target study; species retain separate cohorts',hide=False)
    f=x.groupby(['species','state','pair_class']).size().unstack(fill_value=0).reindex(columns=cats,fill_value=0);f.index=[sp[0].upper()+' '+st for sp,st in f.index];f.columns=names;r.bars(1,f,'Genes',pad=(10,11,1,11))
    r.panel('C');g=allx[allx.scope.eq('requalified')&allx.intermediate_threshold.eq(.02)&allx.endpoint_retained&allx.pair_measurable];settings=sorted(g.variant.unique());ax=r.ax(0,kind='line',labels=['Omit C012','Omit L00','Omit both','Native PAS','Depth40','p08','p09'],pad=(22,8,2,10))
    for j,(sp,st) in enumerate([('human','R'),('mouse','R'),('mouse','E')]):
        a=g[g.species.eq(sp)&g.state.eq(st)].groupby('variant').pair_class.agg(lambda a:np.mean(a.eq('endpoint_aligned'))).reindex(settings);ax.plot(a,range(len(settings)),'o-',color=PAL[j],ms=3,label=sp[0].upper()+' '+st)
    ax.set(yticks=range(len(settings)),yticklabels=[{'exclude_C012':'Omit C012','exclude_L00':'Omit L00','exclude_both':'Omit both','native_PAS_coordinates':'Native PAS','shared_depth40':'Depth40','shared_p08':'p08','shared_p09':'p09'}[v] for v in settings],ylabel='',xlabel='Endpoint-aligned / requalified',xlim=(0,1),ylim=(len(settings)-.5,-.5));legend(ax,loc='lower right',bbox_to_anchor=(1,1.01),ncol=3);r.record('all-qualified-setting-records',g);r.record('fixed-cohort-outcomes',r.read('intermediate/process/fixed_cohort_outcomes.tsv'))
    fc=r.read('intermediate/process/fixed_cohort_outcomes.tsv');fc=fc[fc.scope.eq('fixed_primary_admission')&fc.variant.isin(settings)].copy();ax=r.ax(1,kind='point_interval',labels=[str(i+1) for i in range(len(settings))],pad=(9,8,2,10))
    for j,(sp,st) in enumerate([('human','R'),('mouse','R'),('mouse','E')]):
        z=fc[fc.species.eq(sp)&fc.state.eq(st)].set_index('variant').reindex(settings);ax.scatter(z.stable/z.fixed_cohort,np.arange(len(settings))+(j-1)*.19,s=13,c=PAL[j],marker=['o','s','^'][j],label=sp[0].upper()+' '+st)
    ax.set(yticks=range(len(settings)),yticklabels=[str(i+1) for i in range(len(settings))],ylim=(len(settings)-.5,-.5),xlim=(0,1),xlabel='Stable class / fixed cohort');legend(ax,loc='lower right',bbox_to_anchor=(1,1.01),ncol=3);r.share_y(r.panel_axes[r.active][:2],'Seven qualification settings in identical order');r.record('fixed-cohort-outcomes-plotted',fc);r.record('setting-key',pd.DataFrame({'index':range(1,len(settings)+1),'setting':settings}));r.note('Settings1–7 follow the linked key: '+', '.join(f'{i+1} {v}' for i,v in enumerate(settings))+'. Left requalifies its denominator; right holds primary admission fixed. Threshold and GEX-distance qualification do not prove a lineage trajectory.')

DRAW.update(D09=D09,D10=D10,D11=D11,D13=D13)
