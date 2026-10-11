"""Native interval/site inference, matched backgrounds and RBP expression context."""
import numpy as np
import pandas as pd

from render_common import BLUE,ORANGE,INK,GREY,DIVERGE,empty,zero

CONTEXTS=[('mouse','SP'),('mouse','MP'),('human','SP'),('human','MP')]
CONTEXT_LABELS=['Mouse P–S','Mouse P–M','Human P–S','Human P–M']


def gene_points(ax,data,value,category):
    f=data[data.category.eq(category)];rng=np.random.default_rng(20261003)
    for i,(sp,co) in enumerate(CONTEXTS):
        q=f[f.species.eq(sp)&f.contrast.eq(co)];x=q[value].to_numpy()
        ax.scatter(x,i+rng.uniform(-.15,.15,len(x)),s=7,color=[BLUE if a<0 else ORANGE if a>0 else GREY for a in x],marker='s' if sp=='mouse' else 'o',alpha=.65,lw=0)
        if len(x):ax.plot(np.median(x),i,'|',ms=6,color=INK)
        ax.text(1.01,i,f'{(x<0).sum()}/{(x==0).sum()}/{(x>0).sum()}',transform=ax.get_yaxis_transform(),va='center',fontsize=6.5)
    ax.set_yticks(range(4),CONTEXT_LABELS);ax.invert_yaxis();zero(ax,False)


def rbp_heat(r,ax,species,n=6,expression=False):
    candidates=r.read('rbp_expression/candidates.tsv');candidates=candidates[candidates.species.eq(species)&candidates.evidence.eq('direct_species')]
    selected=candidates.sort_values(['target_genes','factor'],ascending=[False,True]);selected=selected.head(n) if species=='human' else selected
    factors=selected.factor.tolist()
    if not factors:empty(ax,'No new direct-motif candidates');return None
    if expression:
        table=r.read('rbp_expression/state_expression.tsv');table=table[table.species.eq(species)&table.evidence.eq('direct_species')&table.factor.isin(factors)]
        table=table.drop_duplicates(['factor','route','state']);value='mean_log2CPM1';keys=['route','state'];flag='expressed_any_2CPM'
    else:
        table=r.read('rbp_expression/expression_contrasts.tsv');table=table[table.species.eq(species)&table.evidence.eq('direct_species')&table.factor.isin(factors)]
        table=table.drop_duplicates(['factor','route','contrast']);value='delta_log2CPM1';keys=['route','contrast'];flag='expression_in_either_state'
    frame=table.pivot(index='factor',columns=keys,values=value).reindex(factors)
    support=table.pivot(index='factor',columns=keys,values=flag).reindex(index=factors,columns=frame.columns)
    cmap='viridis' if expression else DIVERGE;im=ax.imshow(np.ma.masked_invalid(frame.to_numpy()),cmap=cmap,vmin=0 if expression else -4,vmax=12 if expression else 4,aspect='auto')
    for i,j in zip(*np.where((support==False).to_numpy()&np.isfinite(frame.to_numpy()))):ax.plot(j,i,'o',ms=2.5,mfc='none',mec='white',mew=.6)
    ax.set_yticks(range(len(factors)),factors,fontsize=6.5);ax.set_xticks(range(len(frame.columns)),[('1' if route in ['CIP','GSE118258'] else '2')+' '+state for route,state in frame.columns],rotation=90,fontsize=6.5);ax.set_title(species.title(),fontsize=7)
    return im


def inclusion_example(r,ax,gene,paired=False,legend=True):
    from matplotlib.lines import Line2D
    table=r.read('examples/inclusion_paired_libraries.tsv.gz' if paired else 'examples/inclusion_state_summary.tsv')
    table=table[table.species.eq('mouse')&table.gene.eq(gene)&table.cutoff.eq(20)]
    if gene=='EIF4EBP2':table=table[table.kind.eq('site')&table.label.str.upper().eq('PCBP1')]
    else:
        table=table[table.kind.eq('interval')]
        if gene=='RAD23B':
            queries=r.read('examples/external_query_features.tsv');queries=queries[queries.gene.str.upper().eq(gene)]
            seeds=queries[queries.kind.eq('site')&queries.factor.eq('AAGTGCT')]
            regions=queries[queries.kind.eq('shared_region')]
            regions=regions.loc[np.array([any(a.start0<=s.start0 and a.end0>=s.end0 for s in seeds.itertuples()) for a in regions.itertuples()],dtype=bool)]
            table=table.merge(regions[['start0','end0']].drop_duplicates(),on=['start0','end0'],validate='many_to_one')
    if table.empty:empty(ax,'No newly qualified '+gene+' interval');return
    features=table[['start0','end0']].drop_duplicates().sort_values(['start0','end0']);colors=[BLUE,ORANGE,'#8064A2',INK]
    if len(features)>4:raise ValueError(gene+' has >4 current intervals; reflow example without dropping features')
    for j,(a,b) in enumerate(features.itertuples(index=False,name=None)):
        for i,route in enumerate(['CIP','WOT']):
            q=table[table.route.eq(route)&table.start0.eq(a)&table.end0.eq(b)]
            if paired:
                x=q.delta.to_numpy()*100;ax.scatter(x,np.full(len(x),i)+(j-(len(features)-1)/2)*.16,s=8,color=colors[j],alpha=.55)
                if len(x):ax.plot(x.mean(),i+(j-(len(features)-1)/2)*.16,'D',ms=4,color=colors[j])
            else:
                q=q.set_index('state').reindex(list('SREPM'));x=np.arange(5)+i*6
                ax.plot(x,q.equal_library*100,'o-',color=colors[j],ms=2.5,label=f'{a}–{b}' if i==0 else None)
                ax.vlines(x,q.minimum*100,q.maximum*100,color=colors[j],lw=.6)
    if paired:zero(ax,False);ax.set_yticks([0,1],['Guo–CIP','Schiebinger']);ax.set_xlabel('Common-library Δ inclusion (pp)')
    else:ax.set_xticks(list(range(5))+list(range(6,11)),list('SREPM')*2,fontsize=6.5);ax.set_ylabel('Inclusion (%)');ax.set_xlabel('Guo–CIP                  Schiebinger')
    handles=[Line2D([0],[0],color=colors[j],marker='o',ms=2.5,label=f'{a}–{b}')
        for j,(a,b) in enumerate(features.itertuples(index=False,name=None))]
    if legend:
        ax.legend(handles=handles,frameon=False,fontsize=6.5,loc='upper left',bbox_to_anchor=(0,-.75),
            title='Genomic '+('site' if gene=='EIF4EBP2' else 'interval')+' [start, end)',title_fontsize=6.5,borderaxespad=0,handlelength=1)


def figure5(r):
    import json
    from matplotlib.lines import Line2D
    from render_common import HOUSE_INK,STATE_COLORS,STUDY_STYLE,LABELS
    from matplotlib.colors import LinearSegmentedColormap
    r.main(5);r.label('A',0,0)
    stats=r.read('site_statistics/region_consistency_statistics.tsv');stats=stats[stats.unit.eq('region_pair')&stats.category.eq('all')]
    equal=r.read('site_statistics/gene_equal_consistency.tsv');equal=equal[equal.unit.eq('region_pair')&equal.category.eq('all')]
    ax=r.axes((24,7,48,17));ex=r.axes((84,7,27,17))
    for i,(sp,co) in enumerate(CONTEXTS):
        a=stats[stats.species.eq(sp)&stats.contrast.eq(co)].iloc[0];g=equal[equal.species.eq(sp)&equal.contrast.eq(co)].iloc[0]
        ax.plot([a.cobs_ci_low*100,a.cobs_ci_high*100],[i,i],color=HOUSE_INK,lw=.8);ax.plot(a.cobs*100,i,'o',color=HOUSE_INK,ms=3.5);ax.plot(a.c0*100,i,'x',color='#899BA5',ms=4,mew=1.3)
        ex.plot([a.delta_c_ci_low*100,a.delta_c_ci_high*100],[i,i],color=HOUSE_INK,lw=.8);ex.plot(a.delta_c*100,i,'o',color=HOUSE_INK,ms=3.5)
        if co=='SP':ex.plot([g.delta_C_lo*100,g.delta_C_hi*100],[i+.18,i+.18],color=HOUSE_INK,lw=.6);ex.plot((g.Cobs-g.C0)*100,i+.18,'D',mfc='white',mec=HOUSE_INK,ms=3)
    ax.set_yticks(range(4),CONTEXT_LABELS,fontsize=7.2);ax.set_xticks([0,50,100]);ax.set_xlim(0,102);ax.set_ylim(3.5,-.5);ax.set_xlabel('Direction agreement (%)',labelpad=1,fontsize=7.2)
    ex.set_yticks([]);ex.set_ylim(3.5,-.5);ex.axvline(0,color=GREY,ls=':',lw=.5);ex.set_xlim(-5,55);ex.set_xticks([0,20,50]);ex.set_xlabel('Cobs−C0 (pp)',labelpad=1,fontsize=7.2)
    r.text(24,1.5,'● Observed    × Marginal C0',fontsize=7.2);r.legend(84,1,[Line2D([],[],color=HOUSE_INK,mfc='white',marker='D',ms=3,ls='',label='Gene-weighted')])
    r.display('A-interval-statistics',stats);r.display('A-equal-gene-statistics',equal)
    r.label('B',0,0);retention=r.read('weighted_retention/retention_summary.tsv')
    for j,cat in enumerate(['miRNA','RBP']):
        ax=r.axes((26+j*47,42,36,17));f=retention[retention.category.eq(cat)].set_index(['species','contrast']).reindex(CONTEXTS)
        ax.barh(range(4),-f.negative,height=.6,color=HOUSE_INK);ax.barh(range(4),f.positive,height=.6,color='#0072B2');ax.axvline(0,color=GREY,lw=.4)
        for i,a in enumerate(f.itertuples()):ax.text(-a.negative-2,i,str(a.negative),ha='right',va='center',fontsize=7.2);ax.text(a.positive+2,i,str(a.positive),ha='left',va='center',fontsize=7.2)
        ax.set_yticks(range(4),CONTEXT_LABELS if j==0 else ['']*4,fontsize=7.2);ax.set_ylim(3.5,-.5);ax.set_xlim(-95,65);ax.set_xticks([-80,-40,0,40],['80','40','0','40']);r.text(44+j*47,36,cat+'-compatible sites',ha='center',fontsize=7.2)
    r.legend(34,64,[Line2D([],[],color=color,marker='s',ls='',ms=4,label=label) for color,label in [(HOUSE_INK,'Net loss'),('#0072B2','Net gain')]],ncol=2)
    r.display('B-retention-counts',retention)
    r.label('C',0,0);cand=r.read('rbp_expression/candidates.tsv');expression=r.read('rbp_expression/expression_contrasts.tsv');expmap=LinearSegmentedColormap.from_list('GEX_house',['#2166AC','#F7F7F7','#B2182B'])
    selected=[]
    for si,sp in enumerate(['mouse','human']):
        f=cand[cand.species.eq(sp)&cand.evidence.eq('direct_species')].sort_values(['target_genes','factor'],ascending=[False,True]).head(6);selected.append(f)
        routes=['CIP','WOT'] if sp=='mouse' else ['GSE118258','GSE150311'];cols=[(route,co) for co in ['SP','MP'] for route in routes]
        q=expression[expression.species.eq(sp)&expression.evidence.eq('direct_species')&expression.factor.isin(f.factor)].drop_duplicates(['factor','route','contrast'])
        mat=q.pivot(index='factor',columns=['route','contrast'],values='delta_log2CPM1').reindex(index=f.factor,columns=pd.MultiIndex.from_tuples(cols));flags=q.pivot(index='factor',columns=['route','contrast'],values='expression_in_either_state').reindex(index=f.factor,columns=mat.columns)
        y=11+si*24;ax=r.axes((135,y,22,19));im=ax.imshow(np.ma.masked_invalid(mat.to_numpy()),cmap=expmap,vmin=-4,vmax=4,aspect='auto');ax.set_yticks(range(len(f)),f.factor,fontsize=7.2,fontstyle='italic');ax.set_xticks(range(4),['1','2','1','2']);ax.xaxis.tick_top();ax.tick_params(length=0,pad=2)
        for a,b in zip(*np.where((flags==False).to_numpy()&np.isfinite(mat.to_numpy()))):ax.plot(b,a,'o',ms=2.3,mfc='white',mec=HOUSE_INK,mew=.5)
        for spine in ax.spines.values():spine.set_visible(False)
        r.text(119,y-7,sp.title(),fontsize=7.2)
    r.text(140,2,'P–S',ha='center',fontsize=7.2);r.text(151,2,'P–M',ha='center',fontsize=7.2)
    cb=r.axes((133,59,24,1.4));r.fig.colorbar(im,cax=cb,orientation='horizontal',ticks=[-4,0,4]);r.text(119,58,'Δ log2\n(CPM+1)',fontsize=7.2);r.display('C-selected-factors',pd.concat(selected))
    r.label('D',0,0)
    ax=r.axes((7,75,29,7));ax.plot([0,1],[.6,.6],color='#C4CDD3',lw=2)
    ax.plot([.12,.65,.88],[.6]*3,'s',color='#0072B2',ms=4);ax.plot(.12,.6,'s',color='#929FA7',ms=4);ax.axis('off')
    r.text(8,84,'Retained       Changed',fontsize=7.2)
    exact=r.read('site_statistics/exact_partial_recurrence_summary.tsv');exact=exact[exact.contrast.eq('SP')&exact.category.isin(['miRNA','RBP'])]
    r.text(44,72,'P–S pairs (genes)',fontsize=7.2)
    r.text(52,76,'miRNA',ha='center',fontsize=7.2);r.text(72,76,'RBP',ha='center',fontsize=7.2)
    for i,sp in enumerate(['mouse','human']):
        r.text(40,80+i*4,sp.title(),fontsize=7.2)
        for j,cat in enumerate(['miRNA','RBP']):a=exact[exact.species.eq(sp)&exact.category.eq(cat)].iloc[0];r.text(56+j*19,80+i*4,f'{a.exact_partial_pairs} ({a.genes})',ha='center',fontsize=7.2)
    partial=r.read('site_statistics/partial_retention_summary.tsv');routes=['CIP','WOT','GSE118258','GSE150311']
    for j,cat in enumerate(['miRNA','RBP']):
        ax=r.axes((13+j*34,96,30,18));r.text(28+j*34,91,cat,ha='center')
        for i,route in enumerate(routes):
            for k,co in enumerate(['SP','MP']):
                a=partial[partial.route.eq(route)&partial.category.eq(cat)&partial.contrast.eq(co)].iloc[0];y=i+(k-.5)*.22;color,marker=STUDY_STYLE[route]
                ax.plot([a.ci_low*100,a.ci_high*100],[y,y],color=color,lw=.8);ax.plot(a.partial_fraction*100,y,'o' if k==0 else 'D',mfc=color if k==0 else 'white',mec=color,ms=3.2)
        ax.set(ylim=(3.6,-.6),xlim=(0,40));ax.set_xticks([0,20,40]);ax.set_yticks([])
        if j==0:
            for i,route in enumerate(routes):ax.plot(-6,i,marker=STUDY_STYLE[route][1],color=STUDY_STYLE[route][0],ms=4,clip_on=False)
    r.text(43,118,'Partial retention (%)',ha='center');r.legend(7,124,[Line2D([],[],color=HOUSE_INK,mfc=fc,marker=m,ms=3,ls='',label=co) for m,fc,co in [('o',HOUSE_INK,'P–S'),('D','white','P–M')]],ncol=2)
    r.legend(39,122,[Line2D([],[],color=STUDY_STYLE[v][0],marker=STUDY_STYLE[v][1],ms=3,ls='',label=label) for v,label in [('CIP','Guo–CIP'),('WOT','Sch.–2i'),('GSE118258','Xing2020'),('GSE150311','Liu2020')]],ncol=2)
    r.display('D-partial-retention',partial);r.display('D-exact-shared-pairs',exact)
    r.label('E',0,0);r.text(94,69,'Eif4ebp2',fontstyle='italic',fontweight='bold');r.text(146,69,'Pcbp1',ha='center')
    features=r.read('examples/inclusion_features.tsv');features=features[features.species.eq('mouse')&features.gene.eq('EIF4EBP2')&features.kind.eq('site')&features.label.str.upper().eq('PCBP1')]
    sites=features[['start0','end0']].drop_duplicates().sort_values('start0',ascending=False);assert len(sites)==3
    pas=r.read('examples/annotated_example_PAS.tsv');pas=pas[pas.species.eq('mouse')&pas.gene.eq('EIF4EBP2')];base=pas.cut_interbase.max();span=base-pas.cut_interbase.min()
    ax=r.axes((94,77,62,8));
    for i,route in enumerate(['CIP','WOT']):
        f=pas[pas.route.eq(route)];ax.plot([0,span],[i,i],color='#BBC6CD',lw=2);ax.vlines(base-f.cut_interbase,i-.12,i+.12,color=STUDY_STYLE[route][0],lw=.7)
        for j,a in enumerate(sites.itertuples()):ax.plot(base-(a.start0+a.end0)/2,i,'s',ms=4,color='#929FA7' if j==0 else '#0072B2')
    ax.set_ylim(1.4,-.7);ax.set_xlim(-span*.03,span*1.03);ax.axis('off')
    for j,a in enumerate(sites.itertuples()):ax.text(base-(a.start0+a.end0)/2,-.45,str(j+1),ha='center',fontsize=7.2)
    summary=r.read('examples/inclusion_state_summary.tsv');summary=summary[summary.species.eq('mouse')&summary.gene.eq('EIF4EBP2')&summary.cutoff.eq(20)&summary.kind.eq('site')&summary.label.str.upper().eq('PCBP1')]
    # Sites 2 and 3 share exactly the same inclusion vector in each study; show one curve after checking.
    variable=sites.iloc[1:];shown=[]
    for i,route in enumerate(['CIP','WOT']):
        f=summary[summary.route.eq(route)&summary.start0.isin(variable.start0)]
        wide=f.pivot(index='state',columns='start0',values='equal_library');assert np.allclose(wide.iloc[:,0],wide.iloc[:,1],equal_nan=True)
        f=f[f.start0.eq(variable.start0.iloc[0])].set_index('state').reindex(list('SREPM'));shown.append(f.reset_index())
        ax=r.axes((94+i*35,96,25,19));color,marker=STUDY_STYLE[route];x=[0,1,2,3,4.4];ax.plot(x[:4],f.equal_library.iloc[:4]*100,color=color,marker=marker,ms=3)
        ax.plot(4.4,f.equal_library.iloc[4]*100,color=color,marker=marker,ms=3);ax.vlines(x,f.minimum*100,f.maximum*100,color=color,lw=.7);ax.axvline(3.7,color=GREY,ls=':',lw=.5)
        ax.set(xlim=(-.3,4.7),ylim=(0,105));ax.set_xticks(x,list('SREPM'));ax.set_yticks([0,50,100]);ax.set_xlabel('State',fontsize=7.2,labelpad=1)
        if i==0:ax.set_ylabel('I (%)',labelpad=1)
        r.text(106.5+i*35,123,LABELS[route],ha='center',fontsize=7.2)
    r.text(125,128,'Pcbp1-compatible sites 2–3',ha='center',fontsize=7.2);r.display('E-site-features',features);r.display('E-inclusion-curves',pd.concat(shown))
    r.label('F',0,0);r.text(12,135,'Rad23b',fontstyle='italic',fontweight='bold')
    r.legend(42,135,[Line2D([],[],color='#8D79B0',marker='s',ls='',ms=4,label='Shared interval'),Line2D([],[],color='#0072B2',marker='v',ls='',ms=4,label='miR-290/302 (3p)'),Line2D([],[],color='#899BA5',ls='--',label='M'),Line2D([],[],color=HOUSE_INK,label='P')],ncol=4)
    native=r.read('examples/example_endpoint_effects.tsv.gz');native=native[native.setting.eq('supported50_p08')&native.gene_symbol.eq('Rad23b')&native.contrast.eq('P_M')]
    paired=r.read('examples/inclusion_paired_libraries.tsv.gz');paired=paired[paired.gene.eq('RAD23B')&paired.cutoff.eq(20)&paired.kind.eq('interval')&paired.start0.eq(55390760)&paired.end0.eq(55391205)]
    query=r.read('examples/external_query_features.tsv');seed=query[query.gene.str.upper().eq('RAD23B')&query.kind.eq('site')&query.factor.eq('AAGTGCT')&query.start0.ge(55390760)&query.end0.le(55391205)].drop_duplicates(['start0','end0']);assert len(seed)==1
    family_records=[]
    for i,route in enumerate(['CIP','WOT']):
        x=11+i*80;color,marker=STUDY_STYLE[route];r.text(x+31,141,LABELS[route],ha='center',color=color,fontsize=7.2)
        a=native[native.route.eq(route)].iloc[0];r.text(x+12,147,f'Native MPRO {a.MPRO:.3f}',ha='center',fontsize=7.2);r.text(x+48,147,f'ΔDWUI {a.delta_DWUI:+.4f}',ha='center',fontsize=7.2)
        models=json.loads(r.artifact(f'intervals/annotation/mouse/{route}/models.json').read_text());model=next(m for m in models.values() if m['gene_symbol']=='Rad23b')
        lib={state:{a['library']:a for a in model['profiles'][state]['libraries'] if a['total']>=20 and a['coverage']>=.5} for state in ['M','P']};common=sorted(set(lib['M'])&set(lib['P']));f=paired[paired.route.eq(route)].sort_values('library_key');assert common==sorted(f.library_key)
        cuts=np.array([a['cut_interbase'] for a in model['pas']]);order=np.argsort(cuts);pos=cuts[order]-cuts.min();ax=r.axes((x,157,30,20))
        for st,ls in [('M','--'),('P','-')]:
            usage=np.mean([np.array(lib[st][k]['counts'])/lib[st][k]['total'] for k in common],axis=0)[order]
            ax.step(np.r_[0,pos,pos[-1]+10],np.r_[100,100*(1-np.cumsum(usage)),0],where='post',color=color if st=='P' else '#899BA5',ls=ls,lw=1)
            family_records.extend(dict(route=route,state=st,cut=int(cut),usage=float(u),libraries=len(common)) for cut,u in zip(cuts[order],usage))
        ax.axvspan(55390760-cuts.min(),55391205-cuts.min(),color='#8D79B0',alpha=.18);sx=seed.start0.iloc[0]-cuts.min();ax.axvline(sx,color='#0072B2',ls=':',lw=.65);ax.plot(sx,103,'v',color='#0072B2',ms=4,clip_on=False)
        ax.set(xlim=(0,pos[-1]+10),ylim=(0,105));ax.set_yticks([0,50,100]);ax.set_xlabel('Distance (nt)',labelpad=1,fontsize=7.2);ax.set_ylabel('Family PAS usage (%)',labelpad=1,fontsize=7.2)
        ax=r.axes((x+45,157,21,20));
        for row in f.itertuples():ax.plot([0,1],[row.reference_inclusion*100,row.target_inclusion*100],color=color,lw=.5,alpha=.35,marker='o',ms=2)
        ax.plot([0,1],[f.reference_inclusion.mean()*100,f.target_inclusion.mean()*100],color=color,lw=1.2,marker=marker,ms=3.5)
        ax.set_xticks([0,1],['M','P']);ax.set_xlim(-.2,1.2);ax.set_ylim(max(0,min(f.reference_inclusion.min(),f.target_inclusion.min())*100-5),min(105,max(f.reference_inclusion.max(),f.target_inclusion.max())*100+5));ax.set_ylabel('Interval inclusion (%)',fontsize=7.2,labelpad=1)
        r.text(x+55.5,182,f'ΔI {100*f.delta.mean():+.1f} pp; n={len(f)}',ha='center',fontsize=7.2)
    r.display('F-family-usage',pd.DataFrame(family_records));r.display('F-paired-inclusion',paired);r.display('F-native-metrics',native);r.display('F-seed',seed)
    r.save('Fig/Figure.5.pdf','Original A–F layout. A reports four selected interval cohorts, marginal C0 and gene-block excess intervals; gene-label correlation p values are not reused. B shows strict-cohort 20-UMI weighted-site net gain/loss counts with incomplete-site denominators in source. C top six direct factors per species by target count, not by GEX. D retains changed-multisite gene-factor denominators and exact shared P–S pair counts. E all three current Pcbp1-compatible sites; sites 2/3 share one verified curve, first site is not universally constitutive. F native metrics are explicitly distinct from annotation-family curves and paired 445-nt interval inclusion; no AGO2 binding claim.')

def supplement19(r):
    r.start('Figure S19',150,214)
    gene=r.read('weighted_retention/gene_retention.tsv')
    r.label('A',2,2,'All strict-change genes, including zero-site genes')
    for i,cat in enumerate(['miRNA','RBP']):
        ax=r.axes((32+i*65,15,42,24));gene_points(ax,gene,'weighted_sites20',cat)
        if i:ax.set_yticklabels([])
        ax.set_xlabel(cat+' weighted site change')
    r.label('B',2,55,'The same genes at 20 and 50 family UMIs')
    for i,cat in enumerate(['miRNA','RBP']):
        ax=r.axes((22+i*69,68,48,28));f=gene[gene.category.eq(cat)]
        for sp,color in [('human',BLUE),('mouse',ORANGE)]:
            q=f[f.species.eq(sp)];ax.scatter(q.weighted_sites20,q.weighted_sites50,s=7,color=color,marker='o' if sp=='human' else 's',alpha=.65,label=sp)
        z=f[['weighted_sites20','weighted_sites50']].to_numpy();finite=z[np.isfinite(z)]
        if len(finite):lo,hi=finite.min(),finite.max();ax.plot([lo,hi],[lo,hi],':',color=GREY,lw=.5)
        ax.set(xlabel='20-UMI weighted change',ylabel='50-UMI weighted change',title=cat)
        if i==0:ax.legend(frameon=False,fontsize=6.5)
    r.label('C',2,111,'Complete MPRO-added pool and structural attrition')
    admission=r.read('union_pool/figure5_MPRO_admission.tsv');models=r.read('intervals/model_effects.tsv');regions=r.read('site_statistics/all_changed_region_pairs.tsv.gz')
    rows=[]
    for sp,co in CONTEXTS:
        genes=set(admission.loc[admission.species.eq(sp)&admission.contrast.eq(co)&~admission.already_DWUI_gene_contrast,'gene_id'])
        q=models[models.species.eq(sp)&models.contrast.eq(co)&models.eligible&models.annotated_pas.ge(2)&models.n_ref50.ge(1)&models.n_target50.ge(1)]
        n=q.groupby('gene_id').route.nunique();valid=set(n[n.eq(2)].index)&genes
        strict=regions[regions.species.eq(sp)&regions.contrast.eq(co)&regions.effect20_a.mul(regions.effect20_b).gt(0)].gene_id
        shared=set(strict)&valid
        rows.append(dict(species=sp,contrast=co,total=len(genes),evaluable=len(valid),strict=len(shared)))
    f=pd.DataFrame(rows);f.to_csv(r.out/'S19-MPRO-added-full-denominators.tsv',sep='\t',index=False);ax=r.axes((32,124,109,21));left=np.zeros(4)
    for value,color,label in [(f.strict,BLUE,'Shared strict'),(f.evaluable-f.strict,ORANGE,'Other evaluable'),(f.total-f.evaluable,GREY,'Unevaluable')]:
        v=value/f.total.replace(0,np.nan);ax.barh(range(4),v,left=left,color=color,label=label);left+=v.fillna(0)
    ax.set_yticks(range(4),CONTEXT_LABELS,fontsize=6.5);ax.invert_yaxis();ax.set_xlim(0,1);ax.legend(frameon=False,ncol=3,fontsize=6.5,loc='upper left',bbox_to_anchor=(-.1,-.2))
    r.label('D',2,167,'All union-cohort region-pair effect magnitudes')
    for i,(sp,co) in enumerate(CONTEXTS):
        ax=r.axes((14+i*35,180,24,24));q=regions[regions.species.eq(sp)&regions.contrast.eq(co)]
        color=np.where(q.effect20_a.mul(q.effect20_b).le(0),GREY,np.where(q.effect20_a.lt(0),BLUE,ORANGE))
        ax.scatter(q.effect_mean_a*100,q.effect_mean_b*100,s=7,color=color,alpha=.6)
        zero(ax);zero(ax,False);ax.set_title(CONTEXT_LABELS[i],fontsize=6.5);ax.set_xlabel('Study A ΔI (pp)',fontsize=6.5)
        if i==0:ax.set_ylabel('Study B ΔI (pp)',fontsize=6.5)
    r.save('Fig/SuppFig.S18.pdf','A/B preserve the full strict-change gene denominator including genes with zero compatible sites. C reports all MPRO additions and fresh terminal-family coverage. D shows every matched changed region pair including discordant effects; it is distinct from the strict shared motif inventory.')


def supplement20(r):
    r.start('Figure S20',150,201)
    gene=r.read('motif_background/gene_density.tsv');summary=r.read('motif_background/density_summary.tsv')
    for i,cat in enumerate(['miRNA','RBP']):
        r.label(['A','B'][i],2+i*77,2,cat+' changed − matched stable density')
        ax=r.axes((31+i*66,15,40,31));gene_points(ax,gene,'density_difference',cat)
        if i:ax.set_yticklabels([])
        for j,(sp,co) in enumerate(CONTEXTS):
            q=summary[summary.species.eq(sp)&summary.contrast.eq(co)&summary.category.eq(cat)];a=q.iloc[0]
            if a.n_genes>=10:ax.plot([a.lo,a.hi],[j,j],color=INK,lw=1)
            if a.n_genes:ax.plot(a.mean_difference,j,'D',ms=3,mfc='white',mec=INK)
        ax.set_xlabel('Sites / kb difference')
    r.label('C',2,63,'Full changed-interval attrition ledger')
    attr=r.read('motif_background/motif_attrition.tsv');ax=r.axes((28,77,111,25))
    for i,(sp,co) in enumerate(CONTEXTS):
        q=attr[attr.species.eq(sp)&attr.contrast.eq(co)]
        values=[len(q),int(q.stable_bp.gt(0).sum()),int(q.equal_length_windows.gt(0).sum()),int(q.gc_alone_windows.gt(0).sum()),int(q.joint_windows.gt(0).sum())]
        ax.plot(range(5),values,'o-',ms=3,color=BLUE if sp=='human' else ORANGE,ls='-' if co=='SP' else '--',label=CONTEXT_LABELS[i])
    ax.set_xticks(range(5),['Changed','Stable sequence','Equal length','GC','GC + position']);ax.set_ylabel('Intervals');ax.legend(frameon=False,ncol=2,fontsize=6.5)
    r.label('D',2,119,'P–S miRNA density under relaxed matching')
    grid=r.read('motif_controls/grid_summary.tsv');grid=grid[grid.cohort.eq('dynamic')&grid.category.eq('miRNA')&grid.contrast.eq('SP')]
    ax=r.axes((31,133,42,26))
    for i,sp in enumerate(['mouse','human']):
        q=grid[grid.species.eq(sp)&grid.setting.isin(['GC0.05_POS0.2','GCinf_POSinf'])]
        for j,a in enumerate(q.itertuples()):
            ax.plot([a.lo,a.hi],[i+(j-.5)*.2]*2,color=GREY,lw=1);ax.plot(a.mean_difference,i+(j-.5)*.2,['o','s'][j],ms=3,color=BLUE if sp=='human' else ORANGE)
    zero(ax,False);ax.set_yticks([0,1],['Mouse','Human']);ax.set_xlabel('Mean sites / kb difference')
    r.label('E',80,119,'Conditional detection power')
    power=r.read('motif_controls/motif_power.tsv');power=power[power.contrast.eq('SP')&power.category.eq('miRNA')&power.alpha.eq(.05)]
    ax=r.axes((99,133,40,26))
    for sp,color in [('human',BLUE),('mouse',ORANGE)]:
        for prevalence,ls in [('all','-'),('quarter','--')]:
            q=power[power.species.eq(sp)&power.prevalence.eq(prevalence)];ax.plot(q.injected_sites_per_kb,q.rejection_rate,ls=ls,color=color,marker='o',ms=2.5,label=sp+' '+prevalence)
    ax.set(ylim=(0,1),xlabel='Injected sites / kb',ylabel='Detection probability');ax.legend(frameon=False,fontsize=6.5,loc='upper left',bbox_to_anchor=(-.4,-.5))
    r.save('Fig/SuppFig.S19.pdf','Matching windows are rebuilt from current terminal families and selected without motif counts. All eight primary tests, direct96/transfer24 sensitivity families, unmatched intervals, fixed and dynamic cohorts, and conditional power retain their own denominators. Small groups are descriptive; the caption must be determined from new q values.')


def supplement21(r):
    r.start('Figure S21',150,212)
    candidates=r.read('rbp_expression/candidates.tsv');states=r.read('rbp_expression/state_expression.tsv')
    mouse=candidates[candidates.species.eq('mouse')&candidates.evidence.eq('direct_species')]
    rows_h=min(12,len(candidates[candidates.species.eq('human')&candidates.evidence.eq('direct_species')]));rows_m=len(mouse)
    cell_height=3.2;total_height=max(10,(rows_h+rows_m)*cell_height)
    if total_height>100:raise ValueError('New direct-motif factors require S21 reflow; none may be silently omitted')
    r.label('A',2,2,'State expression of the highest-ranked human and all mouse direct factors')
    y=15;im=None
    for sp,n in [('human',rows_h),('mouse',rows_m)]:
        h=max(6,n*cell_height);im=rbp_heat(r,r.axes((35,y,94,h)),sp,12,True);y+=h+19
    if im is not None:
        cb=r.fig.add_axes([35/r.width,1-112/r.height,94/r.width,3/r.height]);r.fig.colorbar(im,cax=cb,orientation='horizontal',label='Mean log2(CPM+1)')
    cross=r.read('rbp_expression/cross_study_summary.tsv');sites=r.read('rbp_expression/candidate_sites.tsv.gz')
    for i,co in enumerate(['SP','MP']):
        r.label(['B','C'][i],2+i*77,133,'P–'+co[0]+' candidate expression changes')
        eligible=set(sites.loc[sites.species.eq('human')&sites.contrast.eq(co),'factor']);q=cross[cross.species.eq('human')&cross.contrast.eq(co)&cross.factor.isin(eligible)]
        ax=r.axes((17+i*76,147,50,48))
        for supported in [True,False]:
            z=q[q.expressed_in_both_studies.eq(supported)];ax.scatter(z.effect_a,z.effect_b,s=10,facecolors=INK if supported else 'white',edgecolors=INK,lw=.5,label=f'{"Supported" if supported else "Low / missing"} · {len(z)}')
        zero(ax);zero(ax,False);ax.set(xlabel='Xing Δ log2(CPM+1)',ylabel='Liu Δ log2(CPM+1)');ax.legend(frameon=False,fontsize=6.5)
    r.save('Fig/SuppFig.S20.pdf','Factors are selected from newly recurrent sites before querying unchanged GEX matrices. The heatmap uses the top twelve human direct factors by current target-gene count and every mouse direct factor, with equal cell heights. Scatterplots retain all current human factors and distinguish expression support without applying a post hoc expression filter.')
