"""Main Figure 2 and profile supplements, with newly selected class membership."""
import numpy as np
import pandas as pd

from render_common import CLASSES,COLORS,INK,GREY,BLUE,ORANGE,LABELS,DIVERGE,empty,zero,matrix,ecdf,medians,embedding,join_frozen_scores


FIT_ROUTES=['WOT','CIP','GSE118258']
ROUTE_SHORT={'GSE118258':'Xing','GSE150311':'Liu','178325_0618':'Guan0618','CIP':'CIP','OSK':'OSK','WOT':'WOT'}
ROUTE_MARKERS={'GSE118258':'o','GSE150311':'D','178325_0618':'^','CIP':'o','WOT':'s'}


def profile_data(r):
    nodes=r.read('profiles/all_route_node_positions.tsv.gz')
    members=r.read('profiles/gene_membership_k3.tsv.gz')
    q=nodes[nodes.setting.eq('supported50_p08')&nodes.depth.eq(20)].merge(members,on=['species','gene_id'],validate='many_to_one')
    return q[q.route.isin(FIT_ROUTES+['GSE150311','178325_0618'])],members


def centered_nodes(q):
    keys=['species','route','cluster','gene_id']
    wide=q.pivot(index=keys,columns='node',values='position_index').reindex(columns=list('SREPM'))
    return wide.sub(wide.S,axis=0).reset_index()


def figure2(r):
    from render_common import HOUSE_MAP,HOUSE_INK,STATE_COLORS,STUDY_STYLE
    from matplotlib.lines import Line2D
    from matplotlib.colors import ListedColormap
    r.main(2)
    nodes=r.read('profiles/all_route_node_positions.tsv.gz');members=r.read('profiles/gene_membership_k3.tsv.gz')
    q=nodes[nodes.setting.eq('supported50_p08')&nodes.depth.eq(20)].merge(members,on=['species','gene_id'],validate='many_to_one')
    d=centered_nodes(q);g=r.read('profiles/architecture_gene_metrics.tsv.gz')
    styles=dict(STUDY_STYLE)
    extras={'171127_GETM':('#88719A','D'),'171127_OSKM':('#AAA0B7','v'),'171127_late':('#A5A8AA','P'),
        '178325_0330_perturbation':('#AB9270','v'),'178325_1117':('#818670','<'),
        '178325_1230':('#98917D','h'),'178325_reference':('#787B7C','X')}
    styles.update(extras)
    labels=dict(LABELS);labels.update({'171127_GETM':'Jaber2022–GETM','171127_OSKM':'Jaber2022–OSKM',
        '171127_late':'Jaber2022–late','178325_0330_perturbation':'Guan2022–0330 pert.',
        '178325_1117':'Guan2022–1117','178325_1230':'Guan2022–1230','178325_reference':'Guan2022–reference'})
    route_order={'mouse':['WOT','CIP','OSK','171127_GETM','171127_OSKM','171127_late'],
        'human':['GSE118258','GSE150311','GSE221739','178325_0330','178325_0330_perturbation','178325_0605','178325_0618','178325_0809','178325_1117','178325_1230','178325_reference']}
    r.label('A',0,0);centers=[]
    for col,sp in enumerate(['mouse','human']):
        x=col*80;r.text(x+24,5,sp.title(),ha='center')
        r.text(x+67,5,'Analysis\ngroup',ha='left',fontsize=7.2)
        for j,c in enumerate(CLASSES):
            y=10+j*11.7;ax=r.axes((x+8,y,21,9));f=d[d.species.eq(sp)&d.cluster.eq(c)]
            allvals=[]
            for route in route_order[sp]:
                ff=f[f.route.eq(route)];curve=ff[list('SREPM')].mean();color,marker=styles[route]
                vals=curve.to_numpy();allvals.extend(vals[np.isfinite(vals)])
                good=np.flatnonzero(np.isfinite(vals[:4]));
                if len(good):ax.plot(good,vals[good],color=color,marker=marker,ms=2,ls='-' if len(good)==4 else '--',lw=.75)
                if np.isfinite(vals[4]):ax.plot(4.35,vals[4],marker=marker,color=color,ms=2,ls='none')
                centers.append(dict(species=sp,cluster=c,route=route,**curve.to_dict()))
            ax.axhline(0,color=GREY,lw=.4,ls=':');ax.axvline(3.7,color=GREY,lw=.4,ls=':')
            ax.set_xlim(-.1,4.55);ax.set_xticks([0,1,2,3,4.35],list('SREPM') if j==2 else ['']*5)
            ax.tick_params(pad=1);ax.set_yticks([-.1,0,.1]);
            if allvals:ax.set_ylim(min(-.03,min(allvals)-.012),max(.03,max(allvals)+.012))
            for t in ax.get_xticklabels():t.set_color(STATE_COLORS.get(t.get_text(),HOUSE_INK))
            # State blocks contain every route; unavailable S-referenced contexts stay grey.
            wide=f.set_index(['gene_id','route'])[list('REPM')].rename_axis('state',axis=1).stack(future_stack=True).rename('delta').reset_index()
            mat=wide.pivot(index='gene_id',columns=['state','route'],values='delta')
            columns=pd.MultiIndex.from_product([list('REPM'),route_order[sp]])
            mat=mat.reindex(columns=columns);mat=mat.loc[mat.mean(axis=1).sort_values(na_position='last').index]
            ax=r.axes((x+32,y,31,9));im=ax.imshow(np.ma.masked_invalid(mat.to_numpy()),aspect='auto',cmap=HOUSE_MAP,vmin=-.10,vmax=.10,rasterized=True)
            ax.set_yticks([]);nr=len(route_order[sp]);ax.set_xticks(np.arange(4)*nr+(nr-1)/2,list('REPM') if j==2 else ['']*4);ax.tick_params(length=0,pad=2)
            for z in np.arange(nr,len(columns),nr)-.5:ax.axvline(z,color='white',lw=1)
            for spine in ax.spines.values():spine.set_visible(False)
            for t in ax.get_xticklabels():t.set_color(STATE_COLORS.get(t.get_text(),HOUSE_INK))
            if j==0:
                strip=r.axes((x+32,y-2.5,31,1.5));strip.imshow([list(range(len(columns)))],aspect='auto',cmap=ListedColormap([styles[v][0] for _,v in columns]));strip.axis('off')
            n=int(((members.species==sp)&(members.cluster==c)).sum());r.text(x+65,y+.8,f'{c}\n{n}',fontsize=7.2)
    r.display('A-centroids',pd.DataFrame(centers));r.display('A-members',members)
    # Same legend band as the original figure; study identities, not result titles.
    handles=[Line2D([],[],marker=styles[v][1],color=styles[v][0],ls='none',ms=3,label=labels[v]) for sp in ['mouse','human'] for v in route_order[sp]]
    leg=r.fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.03,1-53/r.height),ncol=3,frameon=False,
                     columnspacing=.7,handletextpad=.35,handlelength=.9,borderaxespad=0,labelspacing=.35)
    r.panel_texts['A'].append(leg)
    r.text(5,49,'Grey  Not estimable    Dashed line  Unobserved intermediate state',fontsize=7.2)
    cb=r.axes((149,54,1.5,19));r.fig.colorbar(im,cax=cb,ticks=[-.1,0,.1],extend='both');cb.set_ylabel('ΔDWUI vs S',labelpad=1)
    r.label('B',0,0)
    ax=r.axes((10,88,25,12));v=[0,.7,.6,.2];ax.plot(range(4),v,color=HOUSE_INK,lw=1)
    for k,st in enumerate('SREP'):ax.plot(k,v[k],'o',ms=3,color=STATE_COLORS[st])
    ax.axhline(0,color=GREY,lw=.4,ls=':');ax.set_xticks(range(4),list('SREP'));ax.set_yticks([]);ax.set_ylabel('DWUI');ax.set_ylim(-.05,.85)
    for t in ax.get_xticklabels():t.set_color(STATE_COLORS[t.get_text()])
    r.text(41,83,'L = |R−S| + |E−R| + |P−E|\nD = |P−S|\nT = 1 − D/L',fontsize=7.2)
    for j,c in enumerate(CLASSES):
        ax=r.axes((10+j*24,111,18,12));f=g[g.cluster.eq(c)]
        for route in FIT_ROUTES:
            v=np.sort(f.loc[f.route.eq(route),'transient_fraction_T'].dropna().to_numpy());color,_=styles[route]
            ax.step(np.r_[0,v,1],np.r_[100,100*(1-np.arange(1,len(v)+1)/len(v)),0],where='post',color=color,lw=.8)
        ax.set(xlim=(0,1),ylim=(0,100));ax.set_xticks([0,.5,1]);ax.set_yticks([0,50,100] if j==0 else []);ax.axvline(.5,color=GREY,lw=.4,ls=':');r.text(19+j*24,107,c,ha='center')
        if j==0:ax.set_ylabel('P(T ≥ t) (%)',labelpad=1)
    r.text(41,127,'Threshold t',ha='center')
    r.label('C',0,0);row_order=[(c,route) for c in CLASSES for route in FIT_ROUTES]
    route_stats=r.read('profiles/architecture_route_summary.tsv');med=[];fractions=[];ret=[]
    for c,route in row_order:
        f=g[g.cluster.eq(c)&g.route.eq(route)];med.append([f.delta_R.mean(),(f.delta_E-f.delta_R).mean(),(f.delta_P-f.delta_E).mean()])
        fractions.append(f[['contribution_SR','contribution_RE','contribution_EP']].mean().to_numpy())
        ret.append(f.endpoint_retention_D_over_L.quantile([.25,.5,.75]).to_numpy())
    ax=r.axes((99,89,22,24));step_image=ax.imshow(med,cmap=HOUSE_MAP,vmin=-.08,vmax=.08,aspect='auto');ax.set_xticks(range(3),['S→R','R→E','E→P'],fontsize=7.2);ax.xaxis.tick_top();ax.tick_params(length=0,pad=2);ax.set_yticks([])
    for i,(c,route) in enumerate(row_order):
        if i%3==1:r.text(86,89+(i+.5)*24/9,c,ha='center',va='center')
        ax.plot(-1.15,i,marker=styles[route][1],color=styles[route][0],ms=3,clip_on=False)
        for j in range(3):ax.text(j,i,'↑' if med[i][j]>0 else '↓',ha='center',va='center',fontsize=7.2,color=HOUSE_INK)
    for sp in ax.spines.values():sp.set_visible(False)
    ax=r.axes((125,89,18,24));left=np.zeros(9)
    for j,st in enumerate('REP'):v=np.asarray(fractions)[:,j];ax.barh(range(9),v,left=left,height=.75,color=STATE_COLORS[st]);left+=v
    ax.set(ylim=(8.6,-.6),xlim=(0,1));ax.set_yticks([]);ax.set_xticks([0,1]);r.text(134,83,'|ΔDWUI| / L',ha='center',fontsize=7.2)
    ax=r.axes((147,89,10,24));
    for i,(lo,mid,hi) in enumerate(ret):ax.plot([lo,hi],[i,i],color=styles[row_order[i][1]][0],lw=.9);ax.plot(mid,i,'o',ms=2.5,color=styles[row_order[i][1]][0])
    ax.set(ylim=(8.6,-.6),xlim=(0,1));ax.set_yticks([]);ax.set_xticks([0,1]);r.text(152,83,'D/L',ha='center',fontsize=7.2)
    r.legend(85,117,[Line2D([],[],color=styles[v][0],marker=styles[v][1],ls='',ms=3,label=labels[v]) for v in FIT_ROUTES],ncol=2)
    cb=r.axes((126,125,26,1.3));r.fig.colorbar(step_image,cax=cb,orientation='horizontal',ticks=[-.08,0,.08]);r.text(119,124,'ΔDWUI',ha='right',fontsize=7.2)
    r.display('C-route-architecture',route_stats)
    r.label('D',0,0);scores=r.read('moran/k3_group_wpp_metacells.tsv.gz');cells=pd.read_csv(r.resource('frozen_cell_embedding'),sep='\t');cells=join_frozen_scores(cells,scores)
    for col,sp in enumerate(['mouse','human']):
        r.text(40+col*80,133,sp.title(),ha='center')
        for j,c in enumerate(CLASSES):
            ax=r.axes((3+col*80+j*25,142,24,21));e=embedding(ax,cells[cells.species_key.eq(sp)],c+'_centered',.10);e.set_cmap(HOUSE_MAP)
            r.text(5+col*80+j*25,138,c,fontsize=7.2)
        cb=r.axes((33+col*80,166,37,1.4));r.fig.colorbar(e,cax=cb,orientation='horizontal',ticks=[-.1,0,.1],extend='both')
        r.text(3+col*80,165,'Gene-centered DWUI\nwithin study',fontsize=7.2)
    r.label('E',0,0);corr=r.read('profile_gex/APA_GEX_gene_correlations.tsv.gz');routes=['WOT','CIP','GSE118258','OSK','GSE150311','178325_0618']
    vals=corr.groupby(['route','cluster']).APA_GEX_r.median().unstack().reindex(index=routes,columns=CLASSES)
    ax=r.axes((33,182,33,24));e=ax.imshow(vals,cmap=HOUSE_MAP,vmin=-1,vmax=1,aspect='auto');ax.set_yticks(range(6),[LABELS[v] for v in routes],fontsize=7.2);ax.set_xticks(range(3),CLASSES);ax.xaxis.tick_top();ax.tick_params(length=0,pad=3)
    for i in range(6):
        for j in range(3):ax.text(j,i,f'{vals.iloc[i,j]:.2f}',ha='center',va='center',fontsize=7.2)
    for t,route in zip(ax.get_yticklabels(),routes):t.set_color(styles[route][0])
    ax.axhline(2.5,color=HOUSE_INK,ls=':',lw=.7)
    ax=r.axes((70,182,39,24));
    for i,route in enumerate(routes):
        f=corr[corr.route.eq(route)].APA_GEX_r.dropna();lo,mid,hi=f.quantile([.25,.5,.75]);color,marker=styles[route]
        ax.plot([f.min(),f.max()],[i,i],color=color,lw=.55);ax.barh(i,hi-lo,left=lo,height=.55,color=color,alpha=.17);ax.plot(mid,i,marker=marker,color=color,ms=3)
    ax.set(ylim=(5.5,-.5),xlim=(-1.05,1.05));ax.axvline(0,color=GREY,ls=':',lw=.5);ax.set_yticks([]);ax.set_xticks([-1,0,1]);ax.set_xlabel('r(ΔDWUI, ΔGEX)',labelpad=1);r.text(89.5,178,'All classes',ha='center',fontsize=7.2)
    r.display('E-correlations',corr);r.display('E-class-medians',vals.reset_index())
    r.label('F',0,0);tr=r.read('profile_gex/adjacent_direction_summary.tsv');tr=tr[tr.route.isin(FIT_ROUTES)]
    ax=r.axes((128,183,28,23));out=[]
    for i,c in enumerate(['All']+CLASSES):
        for sp,shift,marker,color in [('mouse',-.12,'o','#207450'),('human',.12,'D','#753399')]:
            f=tr[tr.species.eq(sp)];f=f if c=='All' else f[f.cluster.eq(c)]
            den=f.genes.sum();n=f.same_direction.sum();ax.plot(100*n/den,i+shift,marker=marker,color=color,ms=3)
            out.append(dict(species=sp,cluster=c,same_direction=n,measured_transitions=den,percent=100*n/den))
    ax.set_yticks(range(4),['All']+CLASSES);ax.set_ylim(3.5,-.5);ax.set_xlim(25,70);ax.set_xticks([30,50,70]);ax.axvline(50,color=GREY,lw=.5,ls=':');ax.set_xlabel('Direction agreement (%)',labelpad=1,fontsize=7.2)
    r.legend(120,177,[Line2D([],[],color=color,marker=m,ls='',ms=3,label=sp) for sp,color,m in [('Mouse','#207450','o'),('Human','#753399','D')]],ncol=2);r.display('F-direction-agreement',pd.DataFrame(out))
    r.save('Fig/Figure.2.pdf','Original A–F layout; current fixed-candidate k3 profiles. M is disconnected from the SREP process. Grey heatmap cells are unavailable S-referenced contexts. E shows route-level gene-correlation medians/IQR and full ranges. F counts transitions in the three complete routes, retaining zero changes in denominators. Local graph structure is documented in the caption and full Moran source.')

def supplement7(r):
    r.start('Figure S7',150,207)
    resolution=r.read('profiles/resolution_comparison.tsv');robust=r.read('profile_controls/clustering_robustness.tsv')
    r.label('A',2,2,'Resolution diagnostics')
    ax=r.axes((14,13,51,27));ax.plot(resolution.k,resolution.species_balanced_silhouette,'o-',color=INK,ms=3)
    ax.axvspan(2.8,3.2,color=GREY,alpha=.2);ax.set(xlabel='Number of classes k',ylabel='Balanced silhouette')
    r.label('B',78,2,'Refitted class stability')
    f=robust[robust.k.eq(3)&robust.method.eq('centered')&robust.status.eq('refitted')].copy()
    ax=r.axes((90,13,42,27));order=list(f.kind.drop_duplicates());medians(ax,f,'kind','ARI',order)
    short={'80pct_gene_subsample':'80% genes','library_reweight':'Library weights','seed':'Seed',
        'expanded_primary_p08_depth20':'Expanded','supported50_p08_depth40':'40 UMI','supported50_p08_depth80':'80 UMI',
        'leave_route_CIP':'−CIP','leave_route_WOT':'−WOT'}
    assert set(order)<=set(short)
    ax.set_xticklabels([short[s] for s in order],rotation=50,ha='right',fontsize=6.5);ax.set_ylabel('ARI')
    secondary=ax.twinx();r.panel_axes['B'].append(secondary)
    v=f.groupby('kind').centroid_RMSE.median().reindex(order);secondary.plot(range(len(order)),v,'s--',color=ORANGE,ms=2.5);secondary.set_ylabel('Centroid RMSE',fontsize=7,color=ORANGE)
    r.label('C',2,57,'Leave-route-out class projection')
    loo=r.read('profiles/leave_route_projection.tsv');ax=r.axes((14,68,51,27))
    ax.plot(range(len(loo)),loo.assignment_agreement_vs_full,'o-',color=BLUE,ms=3,label='Agreement');ax.plot(range(len(loo)),loo.ARI_vs_full,'s--',color=ORANGE,ms=3,label='ARI')
    ax.set_xticks(range(len(loo)),[ROUTE_SHORT[v] for v in loo.heldout_route]);ax.set_ylim(0,1)
    centroids=r.read('profiles/leave_route_centroids.tsv')
    for j,route in enumerate(loo.heldout_route):
        q=centroids[centroids.heldout_route.eq(route)];ax.scatter(j+np.linspace(-.12,.12,len(q)),q.centroid_correlation,s=8,marker='^',color=INK,label='Centroid r' if j==0 else None)
    ax.legend(frameon=False,loc='lower left')
    r.label('D',78,57,'Path and endpoint geometry')
    geometry=r.read('profiles/architecture_gene_metrics.tsv.gz');ax=r.axes((90,68,51,27))
    for i,key in enumerate(['path_length_L','endpoint_displacement_D','transient_fraction_T']):
        values=[geometry.loc[geometry.cluster.eq(c),key].dropna().to_numpy() for c in CLASSES]
        ax.boxplot(values,positions=np.arange(3)+(i-1)*.22,widths=.18,showfliers=False,manage_ticks=False,
            patch_artist=True,boxprops=dict(facecolor=[BLUE,ORANGE,GREY][i],lw=.5),medianprops=dict(color=INK,lw=.7),whiskerprops=dict(lw=.5),capprops=dict(lw=.5))
    ax.set_xticks(range(3),CLASSES);ax.set_ylabel('L / D / T distributions')
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=color,label=label) for color,label in [(BLUE,'L'),(ORANGE,'D'),(GREY,'T')]],frameon=False,ncol=3,loc='lower right',bbox_to_anchor=(1,1.01),handlelength=.8,columnspacing=.6)
    r.label('E',2,112,'Transient fraction across amplitude cutoffs')
    sens=r.read('profiles/trajectory_amplitude_sensitivity.tsv');ax=r.axes((14,123,51,27))
    for (c,measure),g in sens.groupby(['cluster','measure']):
        s=g.groupby('cutoff').median_T.median();ax.plot(s.index,s,color=COLORS[CLASSES.index(c)],ls='-' if measure=='path_length_L' else '--',label=c if measure=='path_length_L' else None)
    ax.set(xlabel='Minimum amplitude',ylabel='Median T');ax.legend(frameon=False,ncol=3)
    r.label('F',78,112,'Fixed-candidate PAS/depth sensitivity')
    ax=r.axes((90,123,51,27));s=f[f.kind.str.contains('depth|expanded')].copy()
    if len(s):
        ax.scatter(range(len(s)),s.assignment_agreement,color=BLUE,s=12,label='Agreement');ax.scatter(range(len(s)),s.ARI,color=ORANGE,s=12,label='ARI')
        ax.set_xticks(range(len(s)),[short[row.kind]+f'\nn={row.genes}' for row in s.itertuples()],fontsize=6.5);ax.set_ylabel('Agreement / ARI');ax.set_ylim(0,1)
        ax.legend(frameon=False,ncol=2,loc='lower center',bbox_to_anchor=(.5,1.01),columnspacing=.6,handletextpad=.3)
    else:empty(ax,'No jointly complete sensitivity cohort')
    fixed=r.read('profile_controls/trajectory_fixed_class_T.tsv');null=r.read('profile_controls/trajectory_reselected_draws.tsv.gz');null=null[null.cluster.eq('C2')]
    for letter,x,measure,label in [('G',2,'median_T','C2 transient fraction'),('H',78,'n','C2 class size')]:
        r.label(letter,x,163,label+' after null reselection');ax=r.axes((x+12,174,51,15))
        for i,route in enumerate(FIT_ROUTES):
            vals=null.loc[null.route.eq(route),measure].dropna()
            if len(vals):
                lo,mid,hi=vals.quantile([.025,.5,.975]);ax.errorbar(i,mid,yerr=[[mid-lo],[hi-mid]],fmt='s',ms=2.5,color=GREY,capsize=2)
            obs=fixed[fixed.route.eq(route)&fixed.cluster.eq('C2')]
            if len(obs):ax.scatter(i,obs.observed_median_T.iloc[0] if measure=='median_T' else obs.genes.iloc[0],s=12,color=INK,zorder=3)
        ax.set_xticks(range(3),['WOT','CIP','Xing']);ax.set_ylabel('T' if measure=='median_T' else 'Genes')
    from matplotlib.lines import Line2D
    r.fig.legend(handles=[Line2D([],[],marker='o',color=INK,ls='',ms=3,label='Observed'),
        Line2D([],[],marker='s',color=GREY,ls='',ms=3,label='Reselected null: median and 95% range')],
        loc='upper center',bbox_to_anchor=(.5,1-200/r.height),ncol=2,frameon=False)
    r.save('Fig/SuppFig.S6.pdf','All profile diagnostics recomputed with newly selected candidates; null bars are 95% simulation ranges, not effect intervals.')


def supplement8(r):
    r.start('Figure S8',150,224);m=r.read('moran/moran_results.tsv')
    r.label('A',2,2,'Expression-neighborhood organization')
    ax=r.axes((25,16,115,30));base=m[m.k.eq(20)&m.test.isin(['total_unrestricted','state_residual_within_state','state_coverage_library_state'])].copy()
    base['label']=base.study_id+' · '+base.cluster
    order=list(base.label.drop_duplicates());methods=list(base.test.drop_duplicates())
    for j,method in enumerate(methods):
        g=base[base.test.eq(method)].set_index('label').reindex(order)
        # Output column names are asserted at rendering; no fallback to historical values.
        color=[BLUE,ORANGE,INK][j]
        label={'total_unrestricted':'Unadjusted','state_residual_within_state':'State only','state_coverage_library_state':'State + coverage'}[method]
        ax.scatter(np.arange(len(order))+(j-1)*.2,g.excess_i,facecolors=[color if v<.05 else 'white' for v in g.q_BH],edgecolors=color,s=8,label=label,marker=['o','s','^'][j])
    zero(ax);ax.set_xticks(range(len(order)),[v.split(' · ')[1] for v in order],fontsize=6.5);ax.set_ylabel('Excess Moran I')
    ax.legend(frameon=False,ncol=3,fontsize=6.5,loc='lower right',bbox_to_anchor=(1,1.02))
    for study in base.study_id.drop_duplicates():
        positions=[i for i,v in enumerate(order) if v.startswith(study+' · ')]
        ax.text(np.mean(positions),-.32,study.replace('GSE',''),transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=6.5)
        if min(positions):ax.axvline(min(positions)-.5,color=GREY,lw=.35)
    r.label('B',2,64,'Fixed-panel and neighborhood sensitivity')
    feasibility=r.read('moran/fixed_panel_feasibility.tsv')
    for j,(subset,dimension,title) in enumerate([(m[m.test.eq('fixed_panel')],'threshold','Fixed panels: color = excess I; text = genes'),
            (m[m.test.eq('state_coverage_library_state')],'k','Neighborhood size')]):
        cov=subset.pivot(index='study_id',columns=[dimension,'cluster'],values='excess_i').sort_index(axis=1)
        qs=subset.pivot(index='study_id',columns=[dimension,'cluster'],values='q_BH').reindex(index=cov.index,columns=cov.columns)
        ax=r.axes((19 if j==0 else 102,77,76 if j==0 else 40,28))
        im=ax.imshow(np.ma.masked_invalid(cov.to_numpy()),aspect='auto',cmap=DIVERGE,vmin=-.3,vmax=.3)
        yy,xx=np.where(qs.to_numpy()<.05);ax.scatter(xx+.28,yy-.28,s=3,color=INK)
        yy,xx=np.where(~np.isfinite(cov.to_numpy()));ax.scatter(xx,yy,s=6,color=GREY,marker='x')
        ax.set_yticks(range(len(cov)),[sid.replace('GSE','') for sid in cov.index] if j==0 else ['']*len(cov),fontsize=6.5)
        values=list(cov.columns.get_level_values(0).unique())
        ax.set_xticks(np.arange(len(values))*3+1,[f'{v:g}' for v in values],fontsize=6.5)
        for line in np.arange(2.5,len(cov.columns),3):ax.axvline(line,color='white',lw=.5)
        ax.set_title(title,fontsize=6.5);ax.set_xlabel('Coverage fraction · C1/C2/C3' if j==0 else 'k · C1/C2/C3',fontsize=6.5)
        if j==0:
            counts=feasibility.pivot(index='study_id',columns=['threshold','cluster'],values='selected_genes').reindex(index=cov.index,columns=cov.columns)
            for row in range(len(cov)):
                for col in range(len(cov.columns)):
                    if np.isfinite(counts.iloc[row,col]):ax.text(col,row,int(counts.iloc[row,col]),ha='center',va='center',fontsize=6.5,color=INK)
    r.fig.colorbar(im,cax=r.axes((53,119,45,1.5)),orientation='horizontal',ticks=[-.3,0,.3])
    r.fig.text(102/150,1-120/r.height,'Excess I',fontsize=6.5,va='center')
    corr=r.read('profile_gex/APA_GEX_gene_correlations.tsv.gz')
    r.label('C',2,132,'APA–GEX correlations');ax=r.axes((15,144,26,22))
    for i,c in enumerate(CLASSES):
        values=corr.loc[corr.cluster.eq(c),'APA_GEX_r'].dropna()
        if len(values):
            box=ax.boxplot([values],positions=[i],widths=.55,showfliers=False,patch_artist=True,
                medianprops=dict(color=INK,lw=.8),whiskerprops=dict(color=COLORS[i],lw=.6),capprops=dict(color=COLORS[i],lw=.6))
            box['boxes'][0].set(facecolor=COLORS[i],alpha=.3,edgecolor=COLORS[i])
        ax.text(i,1.03,str(len(values)),ha='center',va='bottom',fontsize=6.5)
    ax.set_xticks(range(3),CLASSES);zero(ax);ax.set_ylim(-1,1);ax.set_ylabel('Pearson r')
    r.label('D',48,132,'State-count context');ax=r.axes((61,144,26,22))
    for i,c in enumerate(CLASSES):
        g=corr[corr.cluster.eq(c)].groupby('n_states').APA_GEX_r.quantile([.25,.5,.75]).unstack()
        if len(g):ax.errorbar(g.index+(i-1)*.08,g[.5],yerr=[g[.5]-g[.25],g[.75]-g[.5]],fmt='o-',ms=2.5,color=COLORS[i],label=c,lw=.6)
    zero(ax);ax.set_ylim(-1,1);ax.set_xticks(sorted(corr.n_states.unique()));ax.set_xlabel('Observed states');ax.set_ylabel('Median / IQR r')
    r.label('E',95,132,'Pearson / Spearman');ax=r.axes((108,144,30,22))
    for i,key in enumerate(['APA_GEX_r','APA_GEX_r_spearman']):
        g=corr.groupby('cluster')[key].quantile([.25,.5,.75]).unstack().reindex(CLASSES)
        ax.errorbar(np.arange(3)+(i-.5)*.12,g[.5],yerr=[g[.5]-g[.25],g[.75]-g[.5]],fmt='o-' if i==0 else 's--',ms=2.5,label=['Pearson','Spearman'][i],color=[BLUE,ORANGE][i],lw=.6)
    zero(ax);ax.set_ylim(-1,1);ax.set_xticks(range(3),CLASSES);ax.legend(frameon=False,fontsize=6.5,ncol=2,loc='lower center',bbox_to_anchor=(.5,1.01),columnspacing=.6,handlelength=.9)
    r.label('F',2,175,'Adjacent APA–GEX direction agreement')
    tr=r.read('profile_gex/adjacent_direction_summary.tsv');assert tr.same_direction.le(tr.genes).all()
    ax=r.axes((23,186,115,22));tr['label']=tr.route.map(ROUTE_SHORT)+' '+tr.transition.str.replace('_','→');assert tr.label.notna().all();order=list(tr.label.drop_duplicates())
    for i,c in enumerate(CLASSES):
        f=tr[tr.cluster.eq(c)].set_index('label').reindex(order)
        ax.scatter(np.arange(len(order))+(i-1)*.18,f.same_direction/f.genes,s=5+18*f.genes/max(1,tr.genes.max()),color=COLORS[i],label=c)
    ax.set_xticks(range(len(order)),order,rotation=40,ha='right',fontsize=6.5);ax.set_ylim(0,1);ax.set_ylabel('Same sign / measured');ax.legend(frameon=False,ncol=3)
    r.save('Fig/SuppFig.S7.pdf','New class Moran statistics and gene-level APA/GEX correlations. Filled A symbols and B dots indicate q<0.05; grey crosses are unavailable. Adjacent-direction denominators include zero changes.')
