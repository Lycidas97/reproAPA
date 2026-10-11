"""Review-requested additions, invoked only by the existing theme entry point."""
from render_theme import *
from matplotlib.transforms import ScaledTranslation

def data(r,name):
    return r.read('base:'+r.out.name+'/review-data/'+name)

def extra05(r):
    r.panel('D');x=data(r,'E03b-tandem-comparison.tsv');routes=['GSE118258','178325_0618','178325_0809','CIP','WOT'];ax=r.ax(0,kind='bar',pad=(12,7,1,10))
    for j,subset in enumerate(['All genes','Tandem-compatible']):
        z=x[x.subset.eq(subset)].set_index('route').reindex(routes);y=np.arange(5)+(j-.5)*.55
        ax.barh(y,z['mean'],height=.32,color=[INK,BLUE][j],label=subset)
        for yy,a in zip(y,z.itertuples()):ax.text(a.mean+(0.002 if a.mean>=0 else -.002),yy,str(a.n),ha='left' if a.mean>=0 else 'right',va='center',fontsize=7)
    ax.axvline(0,c=GREY,ls=':',lw=.6);ax.set(yticks=range(5),yticklabels=routes,ylim=(4.6,-.6),xlim=(-.11,.045),xlabel='Mean P–S ΔDWUI · n beside bar');legend(ax,ncol=1)
    r.record('tandem-versus-all',x);r.note('The same annotated terminal 3′UTR must contain all primary PASs for tandem compatibility. This is a restricted-gene comparison, not a paired estimator on the all-gene denominator. All nt summaries and P–M extensions are in E03b-all-tandem-settings. Counts beside each bar are genes.')
    r.panel('E');x=data(r,'E03c-state-reference-summary.tsv');q=x;routes=q[q.n.gt(0)].route.drop_duplicates().tolist();names={'171127_GETM':'D3 GETM','171127_OSKM':'D3 OSKM','178325_0330':'Guan2022–0330','178325_0605':'Guan2022–0605'};aa=[]
    for j,(bg,label) in enumerate([('available','Available genes'),('fixed_common_observed_states_this_route','Fixed common genes')]):
        z=q[q.background.eq(bg)].pivot(index='route',columns='state',values='mean').reindex(index=routes,columns=list('REPM'));n=q[q.background.eq(bg)].pivot(index='route',columns='state',values='n').reindex(index=routes,columns=list('REPM'));ax=r.ax(j,kind='matrix',pad=(12 if j==0 else 2,7,2,10));im=ax.imshow(np.ma.masked_invalid(z),cmap=DIV,vmin=-.08,vmax=.08,aspect='auto')
        for i in range(len(routes)):
            for k in range(4):ax.text(k,i,str(int(n.iloc[i,k])) if n.iloc[i,k]>0 else '×',ha='center',va='center',fontsize=7,color='white' if abs(z.iloc[i,k])>.05 else INK)
        ax.set(yticks=range(len(routes)),yticklabels=[names.get(a,a) for a in routes],xticks=range(4),xticklabels=list('REPM'),xlabel='');ax.xaxis.tick_top();ax.text(.5,1,label,transform=ax.transAxes+ScaledTranslation(0,16/72,r.fig.dpi_scale_trans),ha='center',va='bottom',fontsize=8);aa.append(ax)
    r.share_y(aa,'Identical observed-S route rows; n printed per state and background');pos=aa[1].get_position();cb=r.fig.add_axes([pos.x0,pos.y0-4/r.height,min(pos.width*.4,.12),1.4/r.height]);r.fig.colorbar(im,cax=cb,orientation='horizontal',ticks=[-.08,.08]);cb.tick_params(labelsize=7,length=1,pad=1);cb.text(1.12,.5,'Mean ΔDWUI',transform=cb.transAxes,fontsize=7,va='center');r.panel_axes[r.active].append(cb);r.colorbars.append((aa[1],cb))
    r.record('all-state-reference-summaries',x);r.note('Cell colour is mean target minus observed-S DWUI; text is the complete gene denominator, × is unmeasured. Fixed genes intersect all available R/E/P/M contrasts in a route under the current annotation; absent states are not forced into the intersection. All 17 routes, including five routes with no S-relative effect, in the current v2 data are in the linked table. No v3 result is displayed.')

def extra18(r):
    r.panel('E');tx=r.read('geometry/transcript-models.tsv.gz');ex=r.read('supplement_statistics/same-PAS-example-exon-CDS.tsv.gz');pas=r.read('axis:homolog_site_pairs.tsv.gz');pas=pas[pas.human_gene_id.eq('ENSG00000196923')&pas.tolerance.eq(25)&pas.design.eq('qualified_P_M')&pas.pair.eq('H1_C')];aa=[]
    ids={'human':['ENST00000355572','ENST00000355841'],'mouse':['ENSMUST00000069968','ENSMUST00000046246']}
    for j,sp in enumerate(['human','mouse']):
        z=tx[tx.transcript_id.isin(ids[sp])].set_index('transcript_id').loc[ids[sp]].reset_index();assert len(z)==2;origin=z.end0.max();ax=r.ax(j,kind='structure');ax.set_axis_on();ax.spines[['top','right']].set_visible(False)
        for row,a in enumerate(z.itertuples()):
            seg=ex[ex.transcript_id.eq(a.transcript_id)];assert len(seg)>0;ax.hlines(row,(origin-a.end0)/1000,(origin-a.start0)/1000,color=GREY,lw=.7)
            for v in seg.itertuples():
                lo=(origin-v.end0)/1000;hi=(origin-v.start0)/1000;h=.20 if v.kind=='exon' else .4
                ax.add_patch(Rectangle((lo,row-h/2),hi-lo,h,facecolor='#CBD5DB' if v.kind=='exon' else INK,edgecolor='none'))
            ax.text(.02,row-.30,a.transcript_id,transform=ax.get_yaxis_transform(),fontsize=7)
        cuts=abs(pas['human_signed_cut' if sp=='human' else 'mouse_signed_cut'].to_numpy());xx=np.sort((origin-cuts)/1000)
        for k,v in enumerate(xx):ax.vlines(v,-.15,1.25,color=PAL[k%len(PAL)],ls=':',lw=.85);ax.plot(v,1.45,'v',ms=4,color=PAL[k%len(PAL)])
        ax.set(ylim=(1.7,-.65),yticks=[0,1],yticklabels=['Short CDS','Long CDS'],xlabel=('PDLIM7 · chr5' if sp=='human' else 'Pdlim7 · chr13')+' (−)\n5′ → 3′ distance (kb)');ax.xaxis.set_major_locator(ticker.MaxNLocator(4));aa.append(ax)
        legend(ax,handles=[Rectangle((0,0),1,1,fc='#CBD5DB',label='Exon / UTR'),Rectangle((0,0),1,1,fc=INK,label='CDS'),Line2D([],[],c=BLUE,ls=':',label='Matched PAS')],ncol=2)
        r.record(sp+'-displayed-transcripts',z)
    r.share_y(aa,'Two selected annotated transcript rows per species; not protein measurements',hide=False);r.record('all-exon-CDS-segments',ex);r.record('actual-homologous-cuts',pas)
    r.note('Ensembl release 116 exon/CDS coordinates are 0-based half-open; PASs are interbase cuts. The displayed transcript identifiers are selected to illustrate alternative terminal structures and checked against the frozen annotation. PAS usage does not establish protein isoform switching. Full transcript models are linked.')

def extra26(r):
    r.panel('D');x=data(r,'E14d-native-sites.tsv.gz');x=x[x.gene_symbol.eq('Trip12')];a=data(r,'E14d-axis-decomposition-128.tsv');a=a[a.gene_symbol.eq('Trip12')];origin=x.cut_interbase.min();aa=[];family_lo=a.family_start0.min();family_hi=a.family_end0.max();regions=r.read('examples/recurrent_region_pairs.tsv.gz');regions=regions[regions.gene_symbol_a.eq('Trip12')&regions.contrast.eq('MP')].drop_duplicates(['shared_start0','shared_end0'])
    for j,local in enumerate([False,True]):
        ax=r.ax(j,kind='scatter',pad=(10,4,1,11))
        for route,z in x.groupby('route'):
            for state,mk in [('P','o'),('M','s')]:ax.scatter((z.cut_interbase-origin)/1000,z[state+'_usage'],s=17,marker=mk,facecolor=STUDY[route] if state=='P' else 'white',edgecolor=STUDY[route],lw=.7,label='_nolegend_')
        for route in ['CIP','WOT']:ax.plot([],[],c=STUDY[route],label=route)
        if not local:ax.axvspan((family_lo-origin)/1000,(family_hi-origin)/1000,color='#DDE3E7',alpha=.5)
        else:
            for row in regions.itertuples():ax.axvspan((row.shared_start0-origin)/1000,(row.shared_end0-origin)/1000,color='#DDE3E7',alpha=.5)
        ax.set(ylim=(-.02,1.02),xlim=((family_lo-origin)/1000-.2,(family_hi-origin)/1000+.2) if local else (-.7,(x.cut_interbase.max()-origin)/1000+.7),ylabel='Native PAS usage',xlabel=('Terminal family' if local else 'Full PAS axis')+f'\nchr1 (−), coordinate − {origin:,.0f} (kb)');aa.append(ax);legend(ax,handles=[Line2D([],[],ls='',marker='o',c=INK,label='P'),Line2D([],[],ls='',marker='s',mfc='white',mec=INK,label='M')],ncol=4)
    r.share_y(aa,'The same native P/M PAS usage fractions, whole gene and terminal-family zoom');r.record('Trip12-all-native-sites',x);r.note('All current native PASs are retained, including zero-usage sites. The full-axis shade is the frozen terminal UTR family; local shades are the two strict recurrent intervals. Both views use the same genomic origin and usage scale. Native contributors differ from interval contributors.')
    r.panel('E');ax=r.ax(0,kind='bar',pad=(12,7,1,10));cols=['within_family_DWUI','family_mass_DWUI','nonterminal_DWUI'];labs=['Local shift','Family mass','Outside family'];colorset=[BLUE,'#9A7FB4',ORANGE]
    for i,route in enumerate(['CIP','WOT']):
        z=a[a.route.eq(route)].iloc[0];pos=neg=0
        for k,col in enumerate(cols):
            val=z[col]*1000;left=pos if val>=0 else neg;ax.barh(i,val,left=left,height=.55,color=colorset[k],label=labs[k] if i==0 else None)
            if val>=0:pos+=val
            else:neg+=val
        ax.scatter(z.native_delta_DWUI*1000,i,marker='D',s=20,facecolor='white',edgecolor=INK,zorder=5)
    ax.axvline(0,color=GREY,ls=':',lw=.6);ax.set(yticks=[0,1],yticklabels=['CIP','WOT'],ylim=(1.6,-.6),xlabel='Contribution to ΔDWUI × 1,000');legend(ax,handles=[Rectangle((0,0),1,1,fc=c,label=l) for l,c in zip(labs,colorset)]+[Line2D([],[],ls='',marker='D',mfc='white',mec=INK,label='Native net')],ncol=1)
    r.record('Trip12-exact-decomposition',a);r.note('With x measured from the fixed family origin, within-family position, family mass and outside-family contributions sum to native ΔDWUI. The producer residual is verified separately; this algebraic decomposition does not imply a molecular mechanism. Local shift estimates and contributor counts are exported in full.')

def D27(r):
    r.panel('A');x=data(r,'E14d-axis-decomposition-128.tsv');wide=x.pivot(index='gene_id',columns='route',values='axis_span_nt');valid=x.assign(ok=x.family_status.eq('measurable_native_family')).groupby('gene_id').ok.sum().reindex(wide.index);strict=x.groupby('gene_id').strict_region_pairs.first().reindex(wide.index).gt(0);ax=r.ax(0,kind='scatter')
    for k,label,col in [(2,'Both routes',BLUE),(1,'One route',ORANGE),(0,'Neither route',GREY)]:
        z=wide[valid.eq(k)];ax.scatter(z.CIP/1000,z.WOT/1000,c=col,s=17,lw=0,label=label+f' ({len(z)})')
    z=wide[strict];ax.scatter(z.CIP/1000,z.WOT/1000,s=40,facecolor='none',edgecolor=INK,lw=.7,label=f'Strict interval ({int(strict.sum())})');ax.set(xscale='log',yscale='log',xlabel='CIP native span (kb)',ylabel='WOT native span (kb)');bound=wide.to_numpy().max()/1000*1.2;ax.plot([.8,bound],[.8,bound],':',c=GREY,lw=.6);ax.set(xlim=(.8,bound),ylim=(.8,bound));ax.set_xticks([1,10,100],['1','10','100']);ax.set_yticks([1,10,100],['1','10','100']);ax.xaxis.set_minor_locator(ticker.NullLocator());ax.yaxis.set_minor_locator(ticker.NullLocator());legend(ax,ncol=1,title='Terminal-family estimate');r.record('all-axis-diagnostics',x)
    r.note(f'All {wide.shape[0]} predefined low-shift genes are shown. Both-route terminal-family estimates exist for {int(valid.eq(2).sum())} genes; missing estimates remain unavailable, not zero. Original native spans reconstruct DWUI.')
    r.panel('B');w=x.pivot(index='gene_id',columns='route',values='terminal_shift_nt').dropna();names=x.groupby('gene_id').gene_symbol.first();s=x.groupby('gene_id').strict_region_pairs.first();ax=r.ax(0,kind='scatter');ax.scatter(w.CIP,w.WOT,s=16,color=GREY,lw=0)
    for gene,row in w[s.reindex(w.index).gt(0)].iterrows():
        ax.scatter(row.CIP,row.WOT,s=24,c=BLUE);ax.annotate(names[gene],(row.CIP,row.WOT),xytext={'Capzb':(-12,34),'Nxn':(-48,-4),'Rad23b':(14,-23),'Rhoa':(28,10)}.get(names[gene],(5,5)),textcoords='offset points',fontsize=7,fontstyle='italic',arrowprops=dict(arrowstyle='-',lw=.5,color=INK,shrinkA=1,shrinkB=3))
    zeros(ax,True);ax.set(xlabel='CIP terminal shift (nt)',ylabel='WOT terminal shift (nt)');r.record('paired-local-shifts',w);r.note(f'Both-route local shifts are available for {len(w)} genes. Strict-interval genes are labelled. Missing local estimates do not establish absence of an interval change.')
    r.panel('C');order=x.loc[x.strict_region_pairs.gt(0),'gene_symbol'].drop_duplicates().sort_values().tolist();z=x[x.gene_symbol.isin(order)].copy();z['row']=[order.index(a) for a in z.gene_symbol];ax=r.ax(0,kind='point_interval',pad=(17,7,1,10));right=r.ax(1,kind='bar',pad=(2,7,1,10));cols=['within_family_DWUI','family_mass_DWUI','nonterminal_DWUI'];labs=['Local shift','Family mass','Outside family'];colorset=[BLUE,'#9A7FB4',ORANGE]
    for i,route in enumerate(['CIP','WOT']):
        zz=z[z.route.eq(route)];y=zz.row+(i-.5)*.28;ax.scatter(zz.terminal_shift_nt,y,s=20,c=STUDY[route],marker=['s','o'][i],label=route)
        for a,yy in zip(zz.itertuples(),y):
            if not np.isfinite(a.terminal_shift_nt):
                both_missing=z.loc[z.gene_symbol.eq(a.gene_symbol),'terminal_shift_nt'].isna().all()
                if not both_missing or i==0:
                    text_y=a.row if both_missing else yy
                    ax.text(.97,text_y,'NA',transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=7);right.text(.97,text_y,'No family',transform=right.get_yaxis_transform(),ha='right',va='center',fontsize=7)
            pos=neg=0
            for k,col in enumerate(cols):
                v=getattr(a,col)*1000;start=pos if v>=0 else neg;right.barh(yy,v,left=start,height=.23,color=colorset[k])
                if v>=0:pos+=v
                else:neg+=v
            right.scatter(a.native_delta_DWUI*1000,yy,marker='D',s=15,facecolor='white',edgecolor=INK,zorder=5)
    for a in [ax,right]:a.axvline(0,c=GREY,ls=':',lw=.6);a.set(yticks=range(len(order)),yticklabels=order,ylim=(len(order)-.4,-.6))
    for lab in ax.get_yticklabels():lab.set_fontstyle('italic')
    ax.set_xlabel('Terminal-family shift (nt)');right.set_xlabel('Contribution to ΔDWUI × 1,000');legend(ax,ncol=2);legend(right,handles=[Rectangle((0,0),1,1,fc=c,label=l) for l,c in zip(labs,colorset)]+[Line2D([],[],ls='',marker='D',mfc='white',mec=INK,label='Native net')],ncol=2)
    r.share_y([ax,right],'Current strict genes; upper Guo2019–CIP and lower Schiebinger2019–2i in every row');r.record('strict-gene-exact-decomposition',z);r.note('All strict genes are shown with route-specific local estimates and decompositions. Native net displacement can combine local position, family mass and outside-family terms. The low-shift selection is not uniformly labelled axis compression; unavailable decomposition terms are omitted rather than replaced by zero.')

def extend(draw):
    result=dict(draw)
    for key,extra in [('D05',extra05),('D18',extra18),('D26',extra26)]:
        old=result[key]
        def call(r,old=old,extra=extra):old(r);extra(r)
        result[key]=call
    result['D27']=D27
    return result
