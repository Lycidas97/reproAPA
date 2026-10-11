from render_theme import *
DRAW={}

def flow(ax,frame,columns,orders=None,titles=None):
    """Exact-count ribbons, with each intermediate node mass conserved."""
    total=len(frame);orders=orders or [sorted(frame[c].unique()) for c in columns];positions=[]
    for j,(c,order) in enumerate(zip(columns,orders)):
        counts=frame[c].value_counts();y=0;pos={}
        for i,v in enumerate(order):
            h=counts.get(v,0)/total*.86;pos[v]=[y,h];ax.add_patch(Rectangle((j-.028,y),.056,h,color=dataset_color(v,PAL[i%6])));y+=h+.025
        positions.append(pos)
        from matplotlib.transforms import ScaledTranslation
        ax.text(j,.94,(titles or columns)[j],ha='center',va='bottom',fontsize=8,transform=ax.transData+ScaledTranslation(0,29/72,ax.figure.dpi_scale_trans))
        for i,v in enumerate(order):
            ax.text(j,.94,f'{v}  {counts.get(v,0)}',ha='center',va='bottom',fontsize=7,color=dataset_color(v,PAL[i%6]),transform=ax.transData+ScaledTranslation(0,(20-i*8)/72,ax.figure.dpi_scale_trans))
    for j in range(len(columns)-1):
        usedleft={v:0 for v in orders[j]};usedright={v:0 for v in orders[j+1]}
        for (v,w),g in frame.groupby([columns[j],columns[j+1]],sort=True):
            h=len(g)/total*.86;a=positions[j][v][0]+usedleft[v];b=positions[j+1][w][0]+usedright[w];usedleft[v]+=h;usedright[w]+=h
            verts=[(j+.029,a),(j+.48,a),(j+.52,b),(j+.971,b),(j+.971,b+h),(j+.52,b+h),(j+.48,a+h),(j+.029,a+h),(j+.029,a)]
            ax.add_patch(PathPatch(MPath(verts,[1,4,4,4,2,4,4,4,79]),facecolor=dataset_color(v,PAL[orders[j].index(v)%6]),alpha=.35,edgecolor='none'))
    upper=max(a+h for pos in positions for a,h in pos.values())
    ax.set(xlim=(-.36,len(columns)-.6),ylim=(0,upper));ax.set_axis_off();ax._row_alignable=True

def D04(r):
    a=r.pdefs['D04A'];a['question']='最终参考下原生PAS及UMI获得多少支持？';a['selection']='仅使用PolyASite v2.0整合主版；原生primary PAS/UMI分别为完整分母，支持定义为同链±50 nt，不将支持率称为准确率。'
    a['components'][1]['label']='Primary PAS count and UMI denominator'
    c=r.pdefs['D04C'];c['question']='匹配窗口与UMI深度门槛如何改变覆盖？';c['selection']='固定v2整合主版和同一primary PAS分母，比较±25/50/100 nt支持窗口；supported50_p08下比较20/40/80 UMI深度。'
    c['sources']=[v for v in c['sources'] if 'v3_ge62_comparison' not in v and v!='reference_comparison/reference_runs.tsv']
    c['components'][0]['label']='v2 primary support-window sensitivity'
    r.panel('A');q=r.read('pas_qc/catalog_summary.tsv');q=q[q.layer.eq('primary')].sort_values('study_id');f=pd.DataFrame({'PAS':q.supported50_sites/q.sites,'UMI':q.supported50_UMI/q.UMI});f.index=q.study_id.map(short)
    ax=r.ax(0,labels=f.index.tolist(),kind='bar');
    for j,c in enumerate(f):ax.barh(np.arange(len(f))+(j-.5)*.32,f[c],height=.29,color=PAL[j],label=c)
    ax.set(yticks=range(len(f)),yticklabels=f.index,ylim=(len(f)-.5,-.5),xlim=(0,1),xlabel='Reference-supported fraction');legend(ax,loc='lower right',bbox_to_anchor=(1,1.015),ncol=2);r.record('supported-PAS-and-UMI',f)
    left=ax;ax=r.ax(1,kind='point_interval',labels=f.index.tolist(),pad=(2,3,22,10))
    ax.scatter(q.sites,np.arange(len(q)),s=18,c=BLUE)
    for i,a in enumerate(q.itertuples()):ax.annotate(f'{a.UMI/1e6:.0f}M',(a.sites,i),xytext=(5,0),textcoords='offset points',va='center',fontsize=7)
    ax.set(xscale='log',xlabel='Primary PAS · label = UMIs',yticks=range(len(q)),yticklabels=f.index,ylim=(len(q)-.5,-.5));ax.xaxis.set_major_locator(ticker.FixedLocator([50000,100000,200000]));ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda v,p:f'{v/1000:.0f}k'));ax.xaxis.set_minor_locator(ticker.NullLocator());ax.xaxis.set_minor_formatter(ticker.NullFormatter());r.share_y([left,ax],'Seven studies in identical order');r.record('catalog-size-and-count-support',q)
    r.note('One study per shared row: supported PAS/UMI fractions at left, catalog PAS count and UMI denominator at right. DWUI is distance-weighted usage; RWUI is rank-weighted usage within a gene, not absolute UTR length.')
    r.panel('B');curves=[];reads=[]
    for row in q.itertuples():
        folder='WOT_GSE115943' if row.study_id=='GSE115943' else row.study_id
        curves.append(r.read(f'pas_qc/{folder}/sequence_curves.tsv.gz').assign(species=row.species));reads.append(r.read(f'pas_qc/{folder}/read_profile_curves.tsv.gz').assign(species=row.species))
    curves=pd.concat(curves);reads=pd.concat(reads)
    for j,kind in enumerate(['Canonical signal','Read ends']):
        ax=r.ax(0,sub=j,n=2,kind='line')
        for sp,ls in [('human','-'),('mouse','--')]:
            for k,cat in enumerate(['supported50','unsupported50']):
                if j==0:g=curves[curves.species.eq(sp)&curves.category.eq(cat)&curves.alignment.eq('representative')&curves.feature.eq('canonical_signal_start')].groupby('position_nt')[['count','denominator']].sum();y=g['count']/g.denominator
                else:g=reads[reads.species.eq(sp)&reads.category.eq(cat)].groupby('position_nt')[['normalized_site_sum','selected_sites']].sum();y=g.normalized_site_sum/g.selected_sites
                ax.plot(g.index,y,color=PAL[k],ls=ls,lw=.85,label=sp[0].upper()+' '+('supported' if k==0 else 'unsupported'));r.record(f'position-{j}-{sp}-{cat}',g.assign(plotted=y))
        ax.set(xlabel='Position relative to PAS (nt)',ylabel=kind+' fraction',xlim=(-80,50) if j==0 else (-400,400),ylim=(0,None));legend(ax,loc='best',ncol=1)
    r.panel('C');x=r.read('reference_comparison/PAS_support_comparison.tsv');x=x.loc[x.layer.eq('primary'),['study_id','species','window_nt','layer','denominator','supported_v2']].copy();ids=q.study_id.tolist();x.study_id=x.study_id.replace({'WOT_GSE115943':'GSE115943'});ax=r.ax(0,labels=[short(v) for v in ids],kind='point_interval',pad=(20,10,2,11))
    for window,mk,offset in [(25,'o',-.2),(50,'s',0),(100,'^',.2)]:
        g=x[x.window_nt.eq(window)].set_index('study_id').reindex(ids);ax.scatter(g.supported_v2/g.denominator,np.arange(7)+offset,s=11,marker=mk,c=BLUE,label=f'±{window} nt')
    # The main ±50 nt support and the sensitivity view must use the same PAS denominator.
    g=x[x.window_nt.eq(50)].set_index('study_id').reindex(ids)
    assert np.array_equal(g.denominator.to_numpy(),q.sites.to_numpy())
    assert np.array_equal(g.supported_v2.to_numpy(),q.supported50_sites.to_numpy())
    ax.set(yticks=range(7),yticklabels=[short(v) for v in ids],ylim=(6.6,-.6),xlabel='v2-supported / primary PAS',xlim=(0,1));legend(ax,loc='lower right',bbox_to_anchor=(1,1.02),ncol=3);r.record('primary-support-windows',x,'Primary-layer v2 columns only; unchanged counts and denominators from the archived reference comparison.')
    n=r.read('profiles/all_route_node_positions.tsv.gz');n=n[n.setting.eq('supported50_p08')];w=n.pivot(index=['species','route','depth','gene_id'],columns='node',values='position_index');z=w[['S','P']].dropna().reset_index();counts=[len(z[z.depth.eq(d)]) for d in [20,40,80]]
    left=ax;z['study']=z.route.replace({'CIP':'GSE103221','WOT':'GSE115943','178325_0618':'GSE178325','178325_0809':'GSE178325'});f=z.groupby(['study','depth']).size().unstack().reindex(index=ids,columns=[20,40,80]);f.index=[short(a) for a in f.index];ax=r.ax(1,kind='bar',labels=f.index.tolist(),pad=(2,8,1,10))
    for j,d in enumerate(f):ax.barh(np.arange(len(f))+(j-1)*.23,f[d],height=.2,color=PAL[j],label=f'{d} UMI')
    for i in np.flatnonzero(f.isna().all(axis=1)):ax.text(.03,i,'Not in route analysis',transform=ax.get_yaxis_transform(),va='center',fontsize=7)
    ax.set(yticks=range(len(f)),yticklabels=f.index,ylim=(len(f)-.5,-.5),xlabel='Measurable gene–route profiles');legend(ax,loc='lower right',bbox_to_anchor=(1,1.01),ncol=3);r.share_y([left,ax],'Seven studies; Guan series profiles summed, Jaber unavailable');r.record('depth-by-study-qualification',f)
    r.note('Only the final PolyASite v2.0 + PolyA_DB v4.1 Main + GENCODE v40/M25 manual polyA_site reference is used; no between-version comparison. ±25/50/100 nt are matching tolerances, not confidence intervals. Depth20/40/80UMI counts use the current measurable gene×route profiles. The right display shares study rows with v2 support-window sensitivity; Guan series are summed as gene×route records, not unique genes. Jaber has no trajectory route profile and is not zero.')


def D06(r):
    r.panel('A');q=r.read('profiles/resolution_comparison.tsv');ax=r.ax(0,kind='line');ax.plot(q.k,q.species_balanced_silhouette,'o-',color=BLUE,ms=3);ax.axvline(3,color=GREY,ls=':');ax.set(xlabel='Number of classes, k',ylabel='Balanced silhouette',xticks=[2,3,4,5,6]);r.record('all-k-diagnostics',q)
    f=r.read('profile_controls/trajectory_fixed_class_T.tsv');n=r.read('profile_controls/trajectory_reselected_draws.tsv.gz');r.ecdf(1,[('Observed class medians',f.observed_median_T),('Fixed-class null medians',f.null_median_T),('Reselected draws',n.median_T)],xlabel='Transient fraction T');r.note('T ECDFs distinguish nine observed/fixed-class summaries from the full candidate-reselection draw distribution; they are not equally replicated experiments.')
    r.panel('B');x=r.read('prepared:E04-LOO-gene-assignments.tsv.gz');ax=r.ax(0,kind='structure',pad=(10,1,4,1));f=x[['heldout_route','full_cluster','projected_cluster']].copy();f.heldout_route=f.heldout_route.map(short);flow(ax,f,list(f),orders=[['Xing','CIP','WOT'],['C1','C2','C3'],['C1','C2','C3']],titles=['Held out','Full class','Projected']);r.record('all-assignments',x)
    q=r.read('profiles/leave_route_projection.tsv');rows=[]
    for i,row in enumerate(q.itertuples()):
        for label,col in [('Agreement','assignment_agreement_vs_full'),('ARI','ARI_vs_full')]:rows.append(dict(_row=i,metric=label,value=getattr(row,col),n=row.n_heldout))
    r.dots(1,pd.DataFrame(rows),'value',q.heldout_route.map(short).tolist(),series='metric',xlabel='Agreement / ARI',reference=1);legend(r.panel_axes[r.active][-1],loc='lower right',bbox_to_anchor=(1,1.015),ncol=2);r.note('Candidates were selected using full data. The flow has current gene×route assignments; holding out a route does not reselect the candidate set.')
    r.panel('C');x=r.read('prepared:E04-center-display.tsv');q=r.read('profiles/leave_route_centroids.tsv').sort_values(['heldout_route','cluster']);labels=[short(a.heldout_route)+' '+a.cluster for a in q.itertuples()];ax=r.ax(0,kind='line',labels=labels,pad=(20,5,2,10))
    for i,a in enumerate(q.itertuples()):
        vals=x[x.heldout_route.eq(a.heldout_route)&x.cluster.eq(a.cluster)].set_index('state').reindex(list('SREP')).difference;ax.hlines(i,0,3,color='#DDE3E7',lw=.6);ax.plot(range(4),i-vals*.7,color=PAL[['C1','C2','C3'].index(a.cluster)],lw=.9,marker='o',ms=2)
    ax.set(yticks=range(9),yticklabels=labels,ylim=(8.5,-.5),xticks=range(4),xticklabels=list('SREP'),xlabel='Centroid difference · local ±0.5');left=ax
    aa=r.ax(1,sub=0,n=2,kind='point_interval',pad=(2,5,1,10));aa.scatter(q.centroid_correlation,range(9),s=15,c=INK);aa.set(ylim=(8.5,-.5),yticks=range(9),xlabel='Shape r',xlim=(.7,1.02),xticks=[.7,1]);bb=r.ax(1,sub=1,n=2,kind='bar',pad=(2,5,1,10));bb.barh(range(9),q.centroid_RMSE,color=BLUE,height=.5);bb.set(ylim=(8.5,-.5),yticks=range(9),xlabel='RMSE',xlim=(0,.3),xticks=[0,.3]);r.share_y([left,aa,bb],'Nine held-out route × class rows');r.record('nine-centroid-recovery-metrics',q);r.record('all-state-centroid-differences',x);r.note('One shared row per route/class: state-resolved centroid differences, centroid correlation, then RMSE. Sparklines use the same local ±0.5 difference range around each row baseline; positive is above baseline. No assignment agreement is substituted for shape recovery.')


def D07(r):
    r.panel('A');x=r.read('hyp:B_shape/ortholog_shape_members.tsv.gz');ax=r.ax(0,kind='structure');flow(ax,x,['human_cluster','mouse_cluster'],titles=['Human','Mouse']);r.record('ortholog-pairs',x)
    null=r.read('supplement_statistics/shape-null-PMF.tsv');ax=r.ax(1,kind='distribution',pad=(10,5,1,11))
    for j,(model,z) in enumerate(null.groupby('model')):
        ax.plot(z.same_class_genes,z.probability,color=PAL[j],label='Margins' if model=='selected_overlap_margins' else 'Measurement');r.record('null-'+model,z)
    obs=int((x.human_cluster==x.mouse_cluster).sum());ax.axvline(obs,color=INK,ls=':');ax.text(obs+.3,.96,str(obs),transform=ax.get_xaxis_transform(),va='top',fontsize=7);ax.set(xlabel='Same-class pairs',ylabel='Null probability');legend(ax,loc='best')
    m=r.read('profiles/gene_membership_k3.tsv.gz');f=m.groupby(['species','cluster']).size().unstack().reindex(['human','mouse']);ax=r.bars(2,f.div(f.sum(axis=1),axis=0),'Class fraction',fraction=True,pad=(9,9,1,11));ax.set_xticks([0,.5,1]);ax.set_yticklabels([f'{sp[0].upper()}\nn={int(f.loc[sp].sum())}' for sp in ['human','mouse']]);r.note('The measured ortholog pairs are compared to both producer conditional nulls. Class proportions use species-specific candidate backgrounds and do not explain global route means.')
    r.panel('B');g=r.read('profiles/architecture_gene_metrics.tsv.gz');grouped=list(g.groupby(['route','cluster']));labels=[short(ro)+' '+cl for (ro,cl),z in grouped]
    p=r.pdefs[r.active];x0,y0,_,hh=p['components'][0]['rectangle_mm'];template=copy.deepcopy(p['components'][0]);p['components']=[copy.deepcopy(template) for _ in range(3)]
    for comp,x,w in zip(p['components'],[1.5,58,110],[53,48.5,48.5]):comp['rectangle_mm']=[x,y0,w,hh];comp['axes']=1
    effects=np.array([[np.mean(z.delta_R-z.delta_S),np.mean(z.delta_E-z.delta_R),np.mean(z.delta_P-z.delta_E)] for _,z in grouped]);lim=float(np.max(np.abs(effects)))
    ax=r.ax(0,labels=labels,kind='matrix',pad=(17,6,2,11))
    for i in range(9):
        if i%3==0:ax.axhspan(i-.5,i+2.5,color='#F2F4F5',zorder=-3)
        for j in range(3):
            v=effects[i,j];ax.plot([j,j+.39*v/lim],[i,i],color=INK,lw=2.2,solid_capstyle='butt');ax.plot(j+.39*v/lim,i,'o',ms=2.6,color=INK)
        if i in [2,5]:ax.axhline(i+.5,color=GREY,lw=.5,ls=':')
    for j in range(3):ax.vlines(j,-.5,8.5,color='#CBD3D8',lw=.5,zorder=-1)
    ax.set(yticks=range(9),yticklabels=labels,xticks=range(3),xticklabels=['S→R','R→E','E→P'],xlim=(-.5,2.5),ylim=(8.5,-.5),xlabel='Signed mean ΔDWUI');ax.text(.5,1.04,f'Each cell ±{lim:.3f}',transform=ax.transAxes,ha='center',fontsize=7)
    m=g.groupby(['route','cluster'])[['contribution_SR','contribution_RE','contribution_EP']].mean();m.index=labels;m.columns=['S→R','R→E','E→P'];middle=r.bars(1,m,'Mean stage contribution',fraction=True,pad=(2,6,1,11))
    groups=[(label,z.endpoint_displacement_D/z.path_length_L.replace(0,np.nan)) for label,(_,z) in zip(labels,grouped)];right=r.row_distribution(2,groups,'D/L',pad=(2,6,1,11));right.set_xlim(0,1)
    r.share_y([ax,middle,right],'Nine identical route × class rows; direction, stage contribution and per-gene D/L')
    r.record('mean-signed-stage-effects',pd.DataFrame(effects,index=labels,columns=['S_R','R_E','E_P']));r.record('all-geometry-records',g)
    r.note('Signed mean transition effects, mean absolute stage contributions and all finite gene D/L ratios share the same nine rows. Every glyph cell has the same symmetric ΔDWUI scale; direction is encoded by displacement about each zero. Stage fractions use their own categorical legend. Distribution width is normalized within each row, not a count. Mean contributions are not medians; these stages are not used to explain global route means.')
    r.panel('C');x=r.read('profiles/trajectory_amplitude_sensitivity.tsv');ax=r.ax(0,kind='line')
    for (cl,measure),z in x.groupby(['cluster','measure']):
        a=z.groupby('cutoff').median_T.median();ax.plot(a.index,a,color=PAL[['C1','C2','C3'].index(cl)],ls='-' if measure=='path_length_L' else '--',label=cl+(' L' if measure=='path_length_L' else ' D'))
    ax.set(xlabel='Amplitude cutoff',ylabel='Median T\nacross route medians');legend(ax,loc='best',ncol=3);r.record('all-amplitude-settings',x)
    q=r.read('axis:trajectory_mouse_paired.tsv');q=q[q.measure.eq('dominant_stage_agreement')].copy();q['_row']=range(len(q));ax=r.dots(1,q,'estimate',q.cluster.tolist(),lo='low',hi='high',xlabel='Mouse dominant-stage agreement',reference=None)
    for row in q.itertuples():
        paired=g[g.route.isin(['CIP','WOT'])&g.cluster.eq(row.cluster)].pivot(index='gene_id',columns='route',values='dominant_stage').dropna()
        assert len(paired)==row.n_genes
        a=paired.CIP.value_counts(normalize=True);b=paired.WOT.value_counts(normalize=True);expected=float((a*b).fillna(0).sum())
        ax.scatter(expected,row._row if hasattr(row,'_row') else list(q.cluster).index(row.cluster),marker='D',s=16,facecolor='white',edgecolor=ORANGE,lw=.8,zorder=5)
    legend(ax,handles=[Line2D([],[],marker='o',ls='',c=INK,label='Observed'),Line2D([],[],marker='D',ls='',mfc='white',mec=ORANGE,label='Marginal expected')],ncol=1)
    r.note('Hollow diamonds give the descriptive independent-marginal expectation: sum of CIP and WOT dominant-stage proportions multiplied within the same paired class. No new p value is calculated.');r.note('Mouse CIP/WOT gene pairing is descriptive. Intervals retain the original gene-resampling definition; agreement does not show conservation across species.')

def D08(r):
    r.panel('A');x=r.read('moran/moran_results.tsv');q=x[x.k.eq(20)&(x.threshold.isna()|np.isclose(x.threshold,.8))].copy();q['test_label']=q.test+q.threshold.fillna(0).map(lambda v: '' if v==0 else f' {v:.1f}');columns=list(q.test_label.drop_duplicates());q=q.drop_duplicates(['study_id','cluster','test_label']);m=q.pivot(index=['study_id','cluster'],columns='test_label',values='excess_i');qs=q.pivot(index=['study_id','cluster'],columns='test_label',values='q_BH');m=m.reindex(columns=columns);qs=qs.reindex_like(m);labels=[short(a)+' '+b for a,b in m.index];order=m.index.copy();m.index=labels;qs.index=labels;m.columns=['Coverage','Total','State res.','Lib. res.','Same nodes','Fixed .8'];qs.columns=m.columns
    ax=r.ax(0,kind='point_interval',labels=labels,pad=(24,9,2,10))
    for j,c in enumerate(m):
        yy=np.arange(len(m))+(j-2.5)*.10;v=m[c].to_numpy();ok=np.isfinite(v);sig=qs[c].to_numpy()<.05
        ax.scatter(v[ok],yy[ok],s=11,marker=['o','s','^','D','v','>'][j],facecolor=[PAL[j] if a else 'white' for a in sig[ok]],edgecolor=PAL[j],lw=.6,label=c)
        for i in np.flatnonzero(~ok):ax.text(.99,yy[i],'×',transform=ax.get_yaxis_transform(),ha='right',va='center',color=PAL[j],fontsize=7)
    for i in range(len(m)):
        if i%2==0:ax.axhspan(i-.5,i+.5,color='#F2F4F5',zorder=-3)
    ax.axvline(0,color=GREY,ls=':',lw=.6);ax.set(yticks=range(len(m)),yticklabels=labels,ylim=(len(m)-.5,-.5),xlabel='Excess Moran I');legend(ax,loc='lower right',bbox_to_anchor=(1,1.01),ncol=3);r.record('all-six-controls',m);r.record('original-q-families',qs);r.record('null-column-key',pd.DataFrame({'code':m.columns,'test':columns}));r.note('Coverage, unrestricted total, within-state residual, within-library/state residual, same-node and fixed-panel tests are named in the linked null-column key; filled points use each original BH family. Colored × means not estimable, never zero.')
    f=r.read('moran/fixed_panel_feasibility.tsv');f=f[np.isclose(f.threshold,.8)].set_index(['study_id','cluster']).reindex(order);a=f[['selected_genes']].copy();a.index=[short(s)+' '+c for s,c in a.index];a.columns=['Genes'];r.bars(1,a,'Eligible genes',pad=(19,9,1,10));r.record('complete-metacells',f);r.share_y(r.panel_axes[r.active][:2],'Twenty-one study × class rows')
    r.panel('B');q=x[x.test.eq('state_coverage_library_state')];ax=r.ax(0,kind='line')
    for j,((sid,cl),z) in enumerate(q.groupby(['study_id','cluster'])):ax.plot(z.k,z.excess_i,color=PAL[['C1','C2','C3'].index(cl)],alpha=.55,lw=.65)
    ax.set(xlabel='Neighborhood size k',ylabel='Excess Moran I',xticks=[10,15,20,30]);ax.axhline(0,color=GREY,ls=':');legend(ax,handles=[Line2D([],[],c=PAL[j],lw=.8,label=cl) for j,cl in enumerate(['C1','C2','C3'])],ncol=3);r.record('all-k-sensitivity',q);ff=r.read('moran/fixed_panel_feasibility.tsv');ax=r.ax(1,kind='scatter',pad=(12,3,2,10))
    for j,(cl,g) in enumerate(ff.groupby('cluster')):
        for ro,z in g.groupby('study_id'):ax.plot(z.selected_genes,z.complete_metacells,color=PAL[j],lw=.6,alpha=.45)
        ax.scatter(g.selected_genes,g.complete_metacells,s=10,color=PAL[j],label=cl)
    ax.set(xlabel='Fixed-panel genes',ylabel='Complete metacells',yscale='symlog',ylim=(0,None),yticks=[0,10,100,1000]);ax.set_yticklabels(['0','10','100','1000']);legend(ax,loc='best',ncol=3);r.record('coverage-threshold-feasibility',ff);r.note('C1 blue, C2 orange, C3 purple; every study×class curve remains visible. Coverage and fixed-panel tests answer distinct questions.')
    r.panel('C');g=r.read('profile_gex/APA_GEX_gene_correlations.tsv.gz');a=r.read('profile_gex/adjacent_direction_summary.tsv');a=a[a.route.isin(['GSE118258','CIP','WOT'])&a.transition.isin(['S_R','R_E','E_P'])].copy();keys=list(a[['route','cluster']].drop_duplicates().itertuples(index=False,name=None));labs=[short(ro)+' '+cl for ro,cl in keys];groups=[(label,g[g.route.eq(ro)&g.cluster.eq(cl)].APA_GEX_r) for label,(ro,cl) in zip(labs,keys)];left=r.row_distribution(0,groups,'Per-gene APA–GEX profile r');a['_row']=[keys.index((ro,cl)) for ro,cl in zip(a.route,a.cluster)];a['fraction']=a.same_direction/a.genes;right=r.dots(1,a,'fraction',labs,series='transition',xlabel='Same-sign fraction',reference=.5,pad=(2,9,2,10));legend(right,loc='lower right',bbox_to_anchor=(1,1.01),ncol=3);r.share_y([left,right],'Nine route × class rows');r.record('all-profile-correlations',g);r.note('Gene-level profile correlations and adjacent-stage direction fractions share route/class rows. Correlation distributions use finite measured genes; transition fractions retain their own gene denominators. Full additional-study correlations are retained in source data.')

DRAW.update(D04=D04,D06=D06,D07=D07,D08=D08)

def D05(r):
    r.panel('A');x=endpoint(r.read('global/global_direction_summary.tsv'));q=x[x.contrast.eq('P_S')].copy();routes=list(q.route.drop_duplicates());settings=list(q[['setting','background']].drop_duplicates().itertuples(index=False,name=None));ax=r.ax(0,labels=[short(a) for a in routes],kind='point_interval',pad=(18,9,2,10))
    for j,(setting,bg) in enumerate(settings):
        g=q[q.setting.eq(setting)&q.background.eq(bg)];yy=np.array([routes.index(a) for a in g.route])+(j-(len(settings)-1)/2)*.04;ax.scatter(g['mean'],yy,c=PAL[j//4%6],marker=['o','s','^','D'][j%4],s=9,lw=.4)
    ax.axvline(0,color=GREY,ls=':');ax.set(yticks=range(len(routes)),yticklabels=[short(a) for a in routes],ylim=(len(routes)-.6,-.6),xlabel='P–S route mean ΔDWUI');legend(ax,handles=[Line2D([],[],c=PAL[j],ls='',marker='o',ms=3,label=label) for j,label in enumerate(['50/p08','50/p09','50/MC5','Expanded'])],loc='lower right',bbox_to_anchor=(1,1.01),ncol=4);r.record('all-route-setting-means',q);r.record('sixteen-setting-key',pd.DataFrame(settings,columns=['setting','background']));r.note('Colors in setting order: blue50nt/p08, orange50nt/p09, purple50nt/p08/MC5, green expanded/p08. Shapes circle available, square fixed within route, triangle common within species, diamond common one-to-one orthologs.');r.note('Every source-qualified P–S route mean is shown across all16 setting/background combinations. It is a route ordering, not an induction-method effect.')
    left=ax;d=r.read('axis:layer_paired_decomposition.tsv');d=d[d.contrast.eq('P_S')].set_index('route').reindex(routes).reset_index();ax=r.ax(1,labels=d.route.map(short).tolist(),kind='bar')
    for j,col in enumerate(['common_reestimate_contribution','added_background_contribution']):ax.barh(np.arange(len(d))+(j-.5)*.3,d[col],height=.27,color=PAL[j],label=['Common genes','Added genes'][j])
    ax.set(yticks=range(len(d)),yticklabels=d.route.map(short),ylim=(len(d)-.5,-.5),xlabel='Contribution to mean change');ax.axvline(0,color=GREY,lw=.6);legend(ax,loc='lower right',bbox_to_anchor=(1,1.01));r.share_y([left,ax],'Five source-qualified routes');r.record('same-gene-versus-added-background',d)
    r.panel('B');z=endpoint(r.read('state_vectors/route_effects.primary.tsv.gz'));z=z[z.setting.eq('supported50_p08')&z.contrast.eq('P_S')];left=r.row_distribution(0,[(short(ro),z[z.route.eq(ro)].delta_DWUI) for ro in routes],'Gene P–S ΔDWUI');s=r.read('global/single_source_library_summary.tsv');s=s[s.contrast.eq('P_S')];q=s.groupby('route')['mean'].agg(value='median',lo='min',hi='max').reindex(routes).reset_index();q['_row']=range(len(q));right=r.dots(1,q,'value',q.route.map(short).tolist(),lo='lo',hi='hi',xlabel='Source-setting route mean');r.share_y([left,right],'The same five routes; distributions and source-choice ranges');r.record('all-source-settings',s);r.record('target-date-omission',r.read('global/leave_target_date_summary.tsv'));r.note('Source-setting bands are min–max ranges, not confidence intervals. Full gene distributions are shown beside route means.')
    r.panel('C');x=r.read('hyp:A_naive/RSeT_same_library_triangle_gene_effects.tsv.gz');x=x[x.metric.eq('DWUI')&x.setting.eq('supported50_p08')];r.ecdf(0,[(lab,x[col]) for lab,col in [('hN − P','hN_P'),('hN − S','hN_S'),('P − S','P_S')]],xlabel='One RSeT library · ΔDWUI');r.note(f'Same-library comparisons use {x.gene_id.nunique()} measurable genes. This is a state refinement within one library, not a general naive-state effect.')
    x=r.read('hyp:A_naive/global_effect_summary.tsv');cons=['hN_P_primed_condition','hN_S_D0_D7_pool','P_primed_S_D0_D7_pool'];q=x[x.setting.eq('supported50_p08')&x.background.eq('common_hN_primed_P_pool_S')&x.subset.eq('all_primary_APA')&x.metric.eq('DWUI')&x.contrast.isin(cons)];ax=r.ax(1,kind='structure');ax.set(xlim=(0,1),ylim=(0,1));pos={'hN':(.5,.88),'P':(.85,.15),'S':(.15,.15)}
    for lab,(xx,yy) in pos.items():ax.text(xx,yy,lab,ha='center',fontsize=8,fontweight='bold')
    for con,(aa,bb),(xx,yy) in zip(cons,[('hN','P'),('hN','S'),('P','S')],[(.78,.57),(.21,.57),(.5,.09)]):
        a,b=pos[aa],pos[bb];ax.plot([a[0],b[0]],[a[1]-.02,b[1]+.08],color=GREY,lw=.6);g=q[q.contrast.eq(con)].sort_values('version');ax.text(xx,yy,' / '.join(f'{v:.3f}' for v in g['mean']),ha='center',fontsize=7)
    ax.text(.5,.99,'Separate cross-condition cohort',ha='center',fontsize=7);r.record('cross-condition-triangle',q);r.note('The separate cross-condition cohort compares hN and primed-P conditions with pooled S; all displayed values use v2 only. Its gene denominator is retained in the linked source and differs from the same-library comparison.')
DRAW.update(D05=D05)
