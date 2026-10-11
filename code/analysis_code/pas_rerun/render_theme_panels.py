"""Data-specific compound panels; all entry calls pass through render_candidate."""
from render_theme import *
DRAW={}

def pas_curves(r,k,frame,gene,sub=0,n=1,matched=False):
    q=frame[frame.gene_symbol.str.upper().eq(gene.upper())].copy();assert len(q)
    ax=r.ax(k,sub,n,kind='genomic',pad=(12,13,2,11));allcoords=[]
    axes=r.read('state_vectors/PAS_axes.primary.tsv.gz');axes=axes[axes.setting.eq('supported50_p08')]
    lookup=axes.assign(oriented=axes.cut_site.where(axes.strand.eq('+'),-axes.cut_site)).set_index(['study_id','primary_gene_id','atlas_pas_id']).oriented
    def coordinates(a):
        study={'CIP':'GSE103221','WOT':'WOT_GSE115943'}.get(a.route,a.study_id if hasattr(a,'study_id') else a.route)
        if study=='WOT_GSE115943' and study not in axes.study_id.unique():study='GSE115943'
        return np.array([lookup.loc[(study,a.gene_id,p)] for p in json.loads(a.PAS_ids)])
    for a in q.itertuples():
        ids=json.loads(a.PAS_ids);coords=coordinates(a);allcoords.extend(coords)
    origin=min(allcoords)
    for i,a in enumerate(q.itertuples()):
        ids=json.loads(a.PAS_ids);coords=coordinates(a);order=np.argsort(coords);x=(coords[order]-origin)/1000
        for st,col in ([('S','S_usage'),('P','P_usage'),('M','M_usage')] if matched else [('S','reference_usage'),('P','target_usage')]):
            v=np.asarray(json.loads(getattr(a,col)));assert len(v)==len(x) and np.isclose(v.sum(),1)
            xx=np.r_[x[0]-.03,x,x[-1]+.03];yy=np.r_[0,np.cumsum(v[order]),1]
            ax.step(xx,yy,where='post',color=STUDY.get(a.route,PAL[i%6]),ls={'S':'--','P':'-','M':'-.'}[st],lw=.85,label=short(a.route)+' '+st)
    ax.set(xlabel='5′ → 3′ position (kb)',ylabel='Cumulative PAS usage',ylim=(-.03,1.03));ax.text(0,1.08,gene,transform=ax.transAxes,fontstyle='italic',fontweight='bold',fontsize=8);
    legend(ax,handles=[Line2D([],[],color=INK,ls={'S':'--','P':'-','M':'-.'}[st],label=st) for st in ('SPM' if matched else 'SP')],loc='lower right',bbox_to_anchor=(1,1.01),ncol=3)
    r.record(gene+f'-PAS-curves-{sub}',q,'All qualified PAS fractions; signed genomic positions oriented 5′ to 3′, shared gene origin. Step curves are cumulative usage.');return ax

def D12(r):
    r.panel('A');family=r.read('apa_gex/prespecified_eight_test_family.tsv');nulls=r.read('apa_gex/four_class_nulls.tsv.gz');q=family[family.baseline_bins.eq(5)].sort_values(['species','threshold','null_model']).copy();assert len(q)==8
    comp=r.pdefs[r.active]['components'][0];comp['rectangle_mm'][0]=1.5;comp['rectangle_mm'][2]=157;r.pdefs[r.active]['components']=[comp]
    labs=[('H' if a.species=='human' else 'M')+f' {a.threshold:.2f} '+('ID' if a.null_model=='response_identity' else 'Dir') for a in q.itertuples()];ax=r.ax(0,labels=labs,kind='point_interval',pad=(24,9,36,11));exports=[]
    for i,a in enumerate(q.itertuples()):
        if i%2==0:ax.axhspan(i-.5,i+.5,color='#F2F4F5',zorder=-3)
        z=nulls[nulls.species.eq(a.species)&nulls.threshold.eq(a.threshold)&nulls.null_model.eq(a.null_model)&nulls.baseline_bins.eq(5)&nulls.cell.eq('B')].sort_values('shuffled_study');assert len(z)==2
        for j,v in enumerate(z.itertuples()):
            yy=i+(j-.5)*.25;mu=v.null_mean;assert mu>0
            ax.hlines(yy,v.null_p025/mu,v.null_p975/mu,color=PAL[j],lw=2.5,alpha=.35);ax.scatter(v.observed/mu,yy,s=16,marker=['o','s'][j],facecolor=PAL[j] if a.q_family<.05 else 'white',edgecolor=PAL[j],lw=.7,zorder=3);exports.append(dict(row=i,orientation=j+1,**v._asdict()))
        b3=family[family.species.eq(a.species)&family.threshold.eq(a.threshold)&family.null_model.eq(a.null_model)&family.baseline_bins.eq(3)].iloc[0]
        for xx,v in [(1.10,b3.q_family),(1.31,a.q_family)]:ax.text(xx,i,f'{v:.2g}',transform=ax.get_yaxis_transform(),va='center',ha='center',fontsize=7)
    for xx,label in [(1.10,'q B3'),(1.31,'q B5')]:ax.text(xx,1.04,label,transform=ax.transAxes,ha='center',fontsize=8)
    ax.axvline(1,color=GREY,ls=':',lw=.6);ax.set(yticks=range(8),yticklabels=labs,ylim=(7.5,-.5),xlabel='Observed / null mean · B5');legend(ax,handles=[Line2D([],[],marker='o',color=BLUE,lw=2,alpha=.7,label='Shuffle study 1'),Line2D([],[],marker='s',color=ORANGE,lw=2,alpha=.7,label='Shuffle study 2')],loc='lower left',bbox_to_anchor=(0,1.015),ncol=2)
    r.record('B5-both-orientations',pd.DataFrame(exports));r.record('eight-test-families-B3-B5',family)
    r.note('ID=response-identity shuffle; Dir=sign shuffle with responder identities fixed. Observed counts and each original null 95% range are divided by that orientation’s null mean. The point is observed, not the null mean. Original eight-test q is printed separately for baseline B3/B5; filled markers use B5 q<.05. These null ranges are not observed-effect confidence intervals. Relationship cell B is APA-concordant/GEX-opposite; the same response category is used throughout.')
    r.panel('B');x=r.read('apa_gex/four_class_nulls.tsv.gz');q=x[x.baseline_bins.eq(5)&x.null_model.eq('direction_specific')].drop_duplicates(['species','threshold','cell']);m=q.pivot(index=['species','threshold'],columns='cell',values='observed');m.index=[('H' if sp=='human' else 'M')+f' {t:.2f}' for sp,t in m.index];m.columns=['APA+/GEX+','APA+/GEX−','APA−/GEX+','APA−/GEX−']
    r.bars(0,m,xlabel='Genes',pad=(12,12,2,11));a=x[x.cell.eq('B')&x.baseline_bins.eq(5)&x.null_model.eq('direction_specific')].groupby(['species','threshold']).agg(all=('comparable_genes','max'),movable=('shuffleable_genes','min'))
    f=pd.DataFrame({'Movable':a.movable,'Fixed':a['all']-a.movable});f.index=[('H' if sp=='human' else 'M')+f' {t:.2f}' for sp,t in f.index];ax=r.bars(1,f,'Comparable genes',colorset=[INK,GREY],pad=(13,7,2,10))
    for i,(fixed,movable) in enumerate(zip(a['all'],a.movable)):ax.text(fixed,i,f' {int(movable)}/{int(fixed)}',va='center',fontsize=7)
    ax.set_xlim(0,a['all'].max()*1.4);r.share_y(r.panel_axes[r.active][:2],'Species × threshold, four rows; distinct printed denominators');r.record('all-null-orientations',x);r.note('Legend + means same direction across studies, − means opposite direction. A=APA+/GEX+, B=APA+/GEX−, C=APA−/GEX+, D=APA−/GEX−; no O/E substitution for relationship counts. Complete B3/B5 and all shuffle orientations are linked.')
    r.panel('C');x=endpoint(r.read('state_vectors/route_effects.primary.tsv.gz'));q=x[x.setting.eq('supported50_p08')&x.contrast.eq('P_S')&((x.species.eq('human')&x.gene_symbol.eq('AP1S2'))|(x.species.eq('mouse')&x.gene_symbol.eq('Rtf2')))].sort_values(['gene_symbol','route']);pas_curves(r,0,q,'AP1S2',0,2);pas_curves(r,0,q,'Rtf2',1,2)
    q=q.copy();q['_row']=range(len(q));ax=r.dots(1,q,'GEX_log2FC',[short(a.route) for a in q.itertuples()],xlabel='GEX log2 FC',pad=(19,5,2,11));ax.set_yticklabels([a.gene_symbol+' '+short(a.route) for a in q.itertuples()],fontstyle='italic');r.note('PAS and expression columns are the same source-qualified gene×route comparison. AP1S2 human then Rtf2 mouse rows, in the exported order. S dashed, P solid; color denotes study; no common-contributor values are borrowed from another cohort.')

def D14(r):
    r.panel('A');ora=r.read('profile_gex/enrichment_all_tests.tsv.gz');rank=r.read('profile_gex/rank_enrichment_all_tests.tsv.gz')
    r.note('Class identity is tested by ORA; amplitude by rank enrichment. This definition is in the caption, with the complete test families plotted.')
    counts=pd.DataFrame({'Tests':[len(ora),int(((rank.layer=='supported50')&(rank.measure=='APA')).sum()),int(((rank.layer=='expanded_primary')&(rank.measure=='APA')).sum())]},index=['Class ORA','APA supported','APA expanded']);ax=r.bars(1,counts,'Test positions',pad=(25,4,2,10));ax.xaxis.set_major_locator(ticker.MaxNLocator(3));ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda v,p:f'{v/1000:g}k' if v else '0'))
    r.note('Class and amplitude are separate questions and gene universes. Each rank layer/measure retains its original route-specific correction across both databases.')
    r.panel('B');ax=r.ax(0,0,2,kind='distribution',pad=(12,6,2,11))
    for j,(sp,g) in enumerate(ora.groupby('species')):
        v=np.sort(g.q_value);ax.step(v,np.arange(1,len(v)+1)/len(v),color=PAL[j],label=sp.title())
    ax.set(xlabel='Class ORA q',ylabel='Cumulative fraction',xlim=(0,1),ylim=(0,1));ax.axvline(.05,color=GREY,ls=':',lw=.6);legend(ax,loc='upper left');r.record('complete-ORA',ora)
    ax=r.ax(0,1,2,kind='distribution',pad=(12,6,2,11));summary=[]
    for j,((layer,measure),g) in enumerate(rank.groupby(['layer','measure'])):
        v=np.sort(g.q_both_databases);label=('S' if layer=='supported50' else 'E')+' '+measure
        ax.step(v,np.arange(1,len(v)+1)/len(v),color=BLUE if measure=='APA' else ORANGE,ls='-' if layer=='supported50' else '--',label=label,lw=.9);summary.append([label,len(g),int((v<.05).sum())])
    ax.set(xlabel='Rank q',ylabel='Cumulative fraction',xlim=(0,1),ylim=(0,1));ax.axvline(.05,color=GREY,ls=':',lw=.6);legend(ax,loc='upper left')
    rows=[['ORA',len(ora),int(ora.q_value.lt(.05).sum())]]+summary;ax=r.ax(1,kind='point_interval',labels=[a[0] for a in rows],pad=(14,3,2,10))
    for i,(label,total,hits) in enumerate(rows):
        ax.hlines(i,hits,total,color=GREY,lw=.7);ax.scatter([hits,total],[i,i],s=17,facecolors=[BLUE,'white'],edgecolors=BLUE,lw=.7);ax.text(hits,i-.23,str(hits),fontsize=7,ha='center')
    ax.set(xscale='symlog',xticks=[0,100,10000],yticks=range(len(rows)),yticklabels=[a[0] for a in rows],ylim=(len(rows)-.5,-.7),xlabel='Hits / tested positions');ax.xaxis.set_major_formatter(ticker.ScalarFormatter());r.record('hits-and-test-denominators',pd.DataFrame(rows,columns=['family','tests','q05']));r.record('all-rank-tests',rank)
    r.panel('C');a=rank[rank.measure.eq('APA')];key=['species','route','database','term_id'];w=a.pivot(index=key,columns='layer',values='rank_effect').dropna();ax=r.ax(0,kind='scatter')
    for j,(route,g) in enumerate(w.groupby(level='route')):ax.scatter(g.supported50,g.expanded_primary,s=3,color=STUDY.get(route,PAL[j]),alpha=.4,lw=0,rasterized=True,label=short(route))
    lim=max(abs(w.to_numpy()).max()*1.03,.01);ax.plot([-lim,lim],[-lim,lim],color=GREY,ls=':',lw=.6);ax.set(xlabel='Supported APA rank effect',ylabel='Expanded APA rank effect',xlim=(-lim,lim),ylim=(-lim,lim));legend(ax,loc='best');r.record('all-paired-term-effects',w)
    z=a[a.layer.eq('expanded_primary')].groupby('route').q_both_databases.apply(lambda s:int(s.lt(.05).sum())).to_frame('q<.05');z.index=[short(s) for s in z.index];r.bars(1,z,'Expanded significant positions');r.note('All common terms plotted, including negative and nonsignificant terms. S=supported50, E=expanded; rank q is not a matched-rank test q. The supported and expanded significant counts are shown separately; no absence of rank-enrichment is assumed.')

def D17(r):
    members=None
    for letter,co in [('A','PS'),('B','PM')]:
        r.panel(letter);members=r.read('cross_overlap/new_shared_members.tsv.gz');members=members[members.layer.eq('supported50')&members.threshold.eq(.02)&members.contrast.eq(co)].sort_values('human_gene_symbol')
        nodes=r.read('cross_overlap/all_mapped_nodes.tsv.gz');q=nodes[nodes.layer.eq('supported50')&nodes.measurement_pass&nodes.human_gene_id.isin(members.human_gene_id)]
        codes=['H1','H2','H3','H4','C','W'];states=list('REPM');wide=q.pivot(index='human_gene_id',columns=['code','node'],values='position_index');parts=[]
        for code in codes:
            block=wide[code] if code in wide.columns.get_level_values(0) else pd.DataFrame(index=wide.index)
            for state in states:parts.append((block.get(state,pd.Series(index=wide.index,dtype=float))-block.get('S',pd.Series(index=wide.index,dtype=float))).rename(code+' '+state))
        m=pd.concat(parts,axis=1).reindex(members.human_gene_id);m.index=members.human_gene_symbol
        if co=='PM':
            paired=r.read('cross_overlap/paired_pm_study_effects.tsv.gz');paired=paired[paired.layer.eq('supported50')&paired.quality_pass&paired.human_gene_id.isin(members.human_gene_id)];m=paired.pivot(index='human_gene_id',columns='code',values='mean_delta').reindex(index=members.human_gene_id,columns=codes);m.index=members.human_gene_symbol
        ax=r.matrix(0,m,'ΔDWUI',limit=max(.1,float(np.nanmax(abs(m.to_numpy())))),pad=(24,17 if co=='PS' else 11,2,10));ax.set_xticklabels(m.columns,rotation=90 if co=='PS' else 0,fontsize=7)
        for t in ax.get_yticklabels():t.set_fontstyle('italic')
        values=m.to_numpy(float);iy,ix=np.where(np.isfinite(values)&(abs(values)>=.02));ax.scatter(ix,iy,s=4,c=INK,lw=0,zorder=3)
        for iy,ix in zip(*np.where(np.isfinite(values)&(values<=-.02))):
            if co=='PM' or m.columns[ix].endswith(' P'):ax.add_patch(Rectangle((ix-.48,iy-.48),.96,.96,fill=False,edgecolor=INK,lw=.55,zorder=4))
        ax.text(0,1+(9 if co=='PS' else 5.5)/(ax.get_position().height*r.height),'• |Δ| ≥ .02     □ P contrast ≤ −.02     × unmeasured',transform=ax.transAxes,fontsize=7)
        r.record('cell-threshold-markers',pd.DataFrame(values,index=m.index,columns=m.columns),'Dots: absolute measured effect >=.02; box: P–S or paired P–M <=−.02. Boxes do not assert a separate selection-qualification test.')
        r.record('complete-selected-members',members)
        if co=='PS':
            support=members[['h_PS_effect_n','m_PS_effect_n']].copy();support.index=members.human_gene_symbol;support.columns=['Human','Mouse'];r.bars(1,support,'Studies',pad=(2,9,2,10));r.panel_axes[r.active][-1].set_yticklabels(['']*len(support));r.share_y(r.panel_axes[r.active][:1]+[r.panel_axes[r.active][-1]],f'All {len(members)} shared P–S ortholog genes')
        else:
            native=[]
            for code,g in q.groupby('code'):
                w=g.pivot(index='human_gene_id',columns='node',values='position_index')
                if all(s in w for s in ['P','M']):native.append((w.P-w.M).rename(code))
            w=pd.concat(native,axis=1);z=paired.copy();z['separate']=[w.loc[a.human_gene_id,a.code] if a.human_gene_id in w.index and a.code in w else np.nan for a in z.itertuples()];ax=r.ax(1,kind='point_interval',pad=(2,9,2,10));keys=members.human_gene_id.tolist()
            for j,(code,g) in enumerate(z.groupby('code')):ax.scatter(g.mean_delta-g.separate,[keys.index(a) for a in g.human_gene_id],s=13,marker=['o','s','^','D','v','>'][j],color=STUDY.get({'C':'CIP','W':'WOT','H1':'GSE118258','H2':'GSE150311','H3':'178325_0618','H4':'GSE221739'}[code],PAL[j]),label=short(code),alpha=.8)
            ax.set(yticks=range(len(keys)),ylim=(len(keys)-.5,-.5),xlabel='Paired − separate ΔDWUI');ax.axvline(0,color=GREY,ls=':');legend(ax,loc='best',ncol=1);r.share_y([r.panel_axes[r.active][0],ax],f'All {len(members)} shared P–M ortholog genes');r.record('all-paired-vs-separate',z)
    r.note('Continuous colors retain the full observed range. Threshold dots and proximal-contrast outlines restore visibility near .02 without clipping color limits. Complete current shared proximal member sets; the source table gives each denominator. Human symbols identify one-to-one orthologs. Columns H1=Xing,H2=Liu,H3=Guan,H4=Panariello,C=CIP,W=WOT. Grey × is unmeasured; white is measured near zero. P–M uses paired values.')

def D26(r):
    gene='Trip12'
    r.panel('A');regions=r.read('examples/recurrent_region_pairs.tsv.gz');regions=regions[regions.gene_symbol_a.eq(gene)&regions.contrast.eq('MP')].drop_duplicates(['shared_start0','shared_end0']).sort_values('shared_start0');assert len(regions)==2
    intervals=list(zip(regions.shared_start0.astype(int),regions.shared_end0.astype(int)));pairs=r.read('examples/inclusion_paired_libraries.tsv.gz');pairs=pairs[pairs.gene.eq(gene.upper())&pairs.kind.eq('interval')&pairs.cutoff.eq(20)];pairs=pairs[[ (a.start0,a.end0) in intervals for a in pairs.itertuples() ]];ax=r.ax(0,kind='distribution',pad=(12,7,1,11));groups=[(ro,k) for ro in ['CIP','WOT'] for k in range(len(intervals))];top=pairs[['reference_inclusion','target_inclusion']].max().max()*100*1.18
    for j,(route,k) in enumerate(groups):
        start,end=intervals[k];z=pairs[pairs.route.eq(route)&pairs.start0.eq(start)&pairs.end0.eq(end)];assert len(z)>0
        ax.plot([j-.16,j+.16],z[['reference_inclusion','target_inclusion']].to_numpy().T*100,c=STUDY[route],alpha=.35,lw=.6);ax.scatter(np.repeat([j-.16,j+.16],len(z)),z[['reference_inclusion','target_inclusion']].to_numpy().T.ravel()*100,s=5,c=STUDY[route],lw=0)
        ax.plot([j-.16,j+.16],z[['reference_inclusion','target_inclusion']].mean()*100,c=INK,lw=1.5,marker='o',ms=3);ax.text(j,top*.92,f'n={len(z)}',ha='center',fontsize=7)
    ax.set(xticks=range(len(groups)),xticklabels=[ro+f' I{k+1}\nM → P' for ro,k in groups],ylim=(-2,top),ylabel='Interval inclusion (%)');r.record('two-interval-all-paired-libraries',pairs);r.record('strict-intervals',regions)
    x=r.read('examples/example_endpoint_effects.tsv.gz');x=x[x.gene_symbol.eq(gene)&x.setting.eq('supported50_p08')&x.contrast.eq('P_M')&x.route.isin(['CIP','WOT'])];v=pd.concat([x.assign(metric='DWUI',value=x.delta_DWUI),x.assign(metric='MPRO',value=x.MPRO)]).sort_values(['metric','route']);v['_row']=range(len(v));r.dots(1,v,'value',[a.metric+' '+a.route for a in v.itertuples()],xlabel='Native metric',pad=(20,6,2,11));r.note('Trip12 has two strict recurrent intervals. Lines join M/P in the same physical library at 20 UMI; black lines show equal-library means. Native DWUI and MPRO have their own contributors, distinct from interval contributors.')
    r.panel('B');lr=r.read('longread/compatible_isoform_pairs.tsv');lr=lr[lr.gene.eq(gene)&lr.strict_full_exclusion&lr.same_complete_intron_chain];key=f'|{intervals[0][0]}|{intervals[0][1]}|';lr=lr[lr.feature_id.str.contains(key,regex=False)];assert len(lr)>0;pair=lr.sort_values(['included_transcript','proximal_transcript']).iloc[0];models=r.read('longread/transcript_audit.tsv');models=models[models.transcript_id.isin([pair.included_transcript,pair.proximal_transcript])];assert len(models)==2;origin=intervals[0][0];ax=r.ax(0,kind='genomic',pad=(25,7,1,11));coords=[]
    for j,tx in enumerate([pair.included_transcript,pair.proximal_transcript]):
        a=models[models.transcript_id.eq(tx)].iloc[0];exons=json.loads(a.exons_mm39);lo,hi=min(exons,key=lambda a:a[0]);coords.extend([lo,hi]);ax.broken_barh([((lo-origin)/1000,(hi-lo)/1000)],(j-.13,.26),facecolors='#CBD5DB',edgecolors=INK,linewidths=.5)
        for k,(start,end) in enumerate(intervals):
            left=max(lo,start);right=min(hi,end)
            if right>left:ax.broken_barh([((left-origin)/1000,(right-left)/1000)],(j-.13,.26),facecolors=PAL[k],alpha=.7)
    for k,(start,end) in enumerate(intervals):ax.axvspan((start-origin)/1000,(end-origin)/1000,color=PAL[k],alpha=.06);ax.text(((start+end)/2-origin)/1000,-.35,f'I{k+1}',ha='center',fontsize=7)
    ax.set(yticks=[0,1],yticklabels=[pair.included_transcript.split('.',1)[1],pair.proximal_transcript.split('.',1)[1]],ylim=(1.5,-.6),xlabel=f'Trip12 · chr1 (−)\nCoordinate − {origin:,} (kb)');r.record('strict-structural-pair',lr);r.record('displayed-transcript-models',models)
    eff=r.read('examples/inclusion_effect_sensitivity.tsv');eff=eff[eff.gene.eq(gene.upper())&eff.kind.eq('interval')&eff.contrast.eq('MP')&eff.cutoff.isin([20,50])];eff=eff[[ (a.start0,a.end0) in intervals for a in eff.itertuples() ]].copy();eff['interval']=[intervals.index((a.start0,a.end0))+1 for a in eff.itertuples()];eff=eff.sort_values(['interval','route','cutoff']);eff['_row']=range(len(eff));eff[['equal_library','leave_one_min','leave_one_max']]*=100;r.dots(1,eff,'equal_library',[f'I{a.interval} {a.route} {a.cutoff}' for a in eff.itertuples()],lo='leave_one_min',hi='leave_one_max',xlabel='Interval ΔI (pp)',pad=(24,5,1,10));r.note('Terminal-exon models illustrate the one strict same-intron-chain inclusion/exclusion pair for interval 1. Every model coordinate is lifted and audited by the producer. Influence ranges describe deleting one library, not confidence intervals; both intervals and both UMI cutoffs are retained.')
    r.panel('C');struct=r.read('longread/structural_summary.tsv');struct=struct[struct.gene.eq(gene)&struct.kind.eq('shared_region')];struct=struct[struct.feature_id.str.contains('|MP|',regex=False)].copy();struct['interval']=[next((k+1 for k,(start,end) in enumerate(intervals) if f'|{start}|{end}|' in f),None) for f in struct.feature_id];struct=struct.dropna(subset=['interval']).sort_values('interval');assert len(struct)==2
    counts=struct.set_index('interval')[['included','excluded']].copy();counts['Other']=struct.set_index('interval').isoforms-counts.sum(axis=1);counts.columns=['Included','Excluded','Other'];counts.index=[f'I{int(k)}' for k in counts.index];ax=r.bars(0,counts,'Audited transcript models',colorset=[BLUE,ORANGE,GREY]);
    for j,a in enumerate(struct.itertuples()):ax.text(1.02,j,f'{int(a.strict_pairs)} strict pair',transform=ax.get_yaxis_transform(),fontsize=7,va='center')
    r.record('both-interval-structural-admission',struct);r.note('Structural compatibility is interval-specific: inclusion/exclusion counts do not themselves imply a compatible isoform pair, expression replication or AGO2 binding. No Trip12 AGO2 query is present in the completed bundle; the former Rad23b AGO2 panel is not transferred to this gene.')

DRAW.update(D12=D12,D14=D14,D17=D17,D26=D26)
from render_theme_foundation import DRAW as FOUNDATION
DRAW.update(FOUNDATION)
from render_theme_cis import DRAW as CIS
DRAW.update(CIS)
from render_theme_recurrence import DRAW as RECURRENCE
DRAW.update(RECURRENCE)
from render_theme_functional import DRAW as FUNCTIONAL
DRAW.update(FUNCTIONAL)
from render_theme_regulatory import DRAW as REGULATORY
DRAW.update(REGULATORY)
from render_theme_gex import DRAW as GEX
DRAW.update(GEX)
from render_theme_reference import DRAW as REFERENCE
DRAW.update(REFERENCE)
