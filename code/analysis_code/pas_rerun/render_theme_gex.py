from render_theme import *
DRAW={};G='root:notes/2026-09-24-integration-top-labels/'
FAM={'shared_s':'S','shared_r':'R','shared_e':'E','shared_tp':'P','shared_m':'M','shared_k':'K','human_x':'hX','human_n':'hN','human_u':'hU','mouse_b':'mB','mouse_n':'mN','mouse_i':'mI','mouse_y':'mY','mouse_u':'mU'}
U={'GSE118258':'Xing','GSE150311':'Liu','GSE178325':'Guan','GSE221739':'Panariello','gse103221_cip':'CIP','gse103221_osk':'OSK','gse171127_d3_getm':'D3 GETM','gse171127_d3_oskm':'D3 OSKM','gse171127_d6_pooled':'D6','gse171127_d12_pooled':'D12','wot_common_to_2i':'WOT','GSE178325::donor_0605':'Guan0605','GSE178325::donor_0618':'Guan0618'}
def markers(r,sp):
    r.panel('A');cfg=r.json(G+'marker_selection.json');pal=r.json(G+'palette.json')['state_groups'];genes=[g for group in cfg['markers'][sp].values() for g in group];s=r.read(G+'source_data/states.tsv');s=s[s.species.eq(sp)];m=r.read(G+'source_data/marker_by_state.tsv');m=m[m.species.eq(sp)];states=s.state_key.tolist();assert len(states)==(62 if sp=='human' else 49) and len(genes)==(31 if sp=='human' else 30);ax=r.ax(0,kind='matrix',pad=(18,13,2,9));ix={a:i for i,a in enumerate(states)};iy={a:i for i,a in enumerate(genes)}
    px,py,pw,ph=r.pdefs[r.active]['components'][0]['rectangle_mm'];r.fig.text((px+18)/160,1-(py+1)/r.height,'Color  relative expression 0–1    Area  detected-cell fraction',fontsize=7,va='top')
    v=m[m.gene.isin(genes)&m.state_key.isin(states)];ax.scatter(v.state_key.map(ix),v.gene.map(iy),s=18.8*v.detected_fraction,c=v.relative_expression,cmap='Blues',vmin=0,vmax=1,lw=0)
    f=v[v.fallback_small];ax.scatter(f.state_key.map(ix),f.gene.map(iy),s=24.2,facecolor='none',edgecolor=GREY,lw=.4);ax.set(xlim=(-.5,len(s)-.5),ylim=(len(genes)-.5,-2.2));ax.set_yticks(range(len(genes)),genes,fontsize=8,fontstyle='italic');ax.set_xticks(range(len(s)),[str(i+1) for i in range(len(s))],fontsize=7,rotation=90);ax.tick_params(length=0,pad=2)
    if sp=='human':
        from matplotlib.transforms import ScaledTranslation
        for t in ax.get_xticklabels()[1::2]:t.set_transform(t.get_transform()+ScaledTranslation(0,-4/25.4,r.fig.dpi_scale_trans))
    for i,a in enumerate(s.itertuples()):ax.add_patch(Rectangle((i-.5,-2),1,.6,color=pal[a.family],lw=0));ax.text(i,-.75,FAM[a.family].replace('h','').replace('m',''),ha='center',va='center',fontsize=7)
    for pos in s.groupby('unit',sort=False).size().cumsum().to_numpy()[:-1]:ax.axvline(pos-.5,color=GREY,lw=.4)
    for spine in ax.spines.values():spine.set_visible(False)
    
    start=0
    for unit,g in s.groupby('unit',sort=False):
        ax.text(start+(len(g)-1)/2,-2.55,U.get(unit,unit),ha='center',va='bottom',fontsize=7);start+=len(g)
    ax=r.ax(1,kind='bar',pad=(18,1,2,4));ax.bar(range(len(s)),s.n_libraries_ge20,width=.8,color=[pal[a] for a in s.family]);ax.set(xlim=(-.5,len(s)-.5),xlabel='',ylabel='Libraries',xticks=[]);ax.set_yticks([int(s.n_libraries_ge20.max())]);legend(ax,handles=[Line2D([],[],ls='',marker='s',ms=3,color=pal[f],label=FAM[f]) for f in FAM if f in set(s.family)],ncol=7);r.record('marker-values',v);r.record('state-order-key',s.assign(plot_index=range(1,len(s)+1)));r.note('Frozen upstream annotations and native expression, unchanged by PAS reanalysis. Columns follow the linked 1-based state key; vertical rules separate study/route units. Family strip S source/R remodeling/E epithelial/P pluripotency/M mesenchymal; local families remain species-specific. Ring = small-stratum fallback (WOT L11), not ≥20-cell evidence.')
    r.panel('B');x=r.read(G+'source_data/program_by_state_library.tsv.gz');x=x[x.species.eq(sp)&x.n_cells.ge(20)];ax=r.ax(0,kind='distribution',labels=['ECM','Epithelial','Pluripotency']);units=s.unit.drop_duplicates().tolist()
    view=x;assert (x.family==x.state_key.map(s.set_index('state_key').family)).all()
    families=[f for f in FAM if f in set(view.family)]
    for j,prog in enumerate(['ECM','Epithelial','Pluripotency']):
        for k,family in enumerate(families):
            q=view[view.program.eq(prog)&view.family.eq(family)].sort_values(['state_key','library_key']);offset=(k-(len(families)-1)/2)*.055
            ax.scatter(q.mean_log1p_cp10k,np.full(len(q),j+offset),s=3,color=pal[family],alpha=.55,lw=0,rasterized=True,label=FAM[family] if j==0 else None)
    ax.set(yticks=range(3),yticklabels=['ECM','Epithelial','Pluripotency'],ylim=(2.5,-.5),xlabel='Native mean log1p(CP10k)');legend(ax,loc='lower right',bbox_to_anchor=(1,1.02),ncol=3);r.record('all-qualified-library-programs',x)
    ax=r.ax(1,kind='scatter');ax.scatter(s.n_libraries_ge20,s.max_library_fraction,c=[pal[a] for a in s.family],s=10,edgecolor='white',lw=.2);ax.set(xlabel='Libraries with ≥20 cells',ylabel='Largest library fraction',ylim=(0,1.05));legend(ax,handles=[Line2D([],[],ls='',marker='o',ms=3,color=pal[f],label=FAM[f]) for f in FAM if f in set(s.family)],ncol=4);r.record('state-support-and-concentration',s);r.note('Library programme scores preserve native scales. State-family colors identify both program observations and sampling support. Program rows retain all qualified library×state observations; small vertical offsets separate state families, not score values. Maximum library fraction diagnoses concentration, not replication independence.')
def D01(r):markers(r,'human')
def D02(r):markers(r,'mouse')
def D03(r):
    r.panel('A');c=r.read(G+'source_data/display_cells.tsv.gz');pal=r.json(G+'palette.json')['state_groups'];s=r.read(G+'source_data/states.tsv')
    for j,sp in enumerate(['human','mouse']):
        ax=r.ax(0,sub=j,n=2,kind='embedding',pad=(1,1,1,1));q=c[c.species_key.eq(sp)].sample(frac=1,random_state=27);ax.scatter(q.UMAP1,q.UMAP2,c=q.reviewed_group.map(pal),s=.6,lw=0,rasterized=True);ax.set(xlim=(-5.3,15.1),ylim=(-4.4,15.2),aspect='equal');ax.set_anchor('W');ax.text(.01,.98,sp.title(),transform=ax.transAxes,va='top',fontweight='bold');ax.legend(handles=[Line2D([],[],marker='o',ls='',color=pal[f],ms=3,label=FAM[f]) for f in FAM if f in set(q.reviewed_group)],loc='center left',bbox_to_anchor=(1.01,.5),frameon=False,handlelength=.7,borderpad=0,labelspacing=.3);r.record('frozen-display-'+sp,q)
        # Keep the former projection counts in the data audit, without assigning biological meaning to UMAP2 bins.
        edges=np.linspace(-4.4,15.2,29);histrows=[]
        for state in FAM:
            values,_=np.histogram(q.loc[q.reviewed_group.eq(state),'UMAP2'],bins=edges)
            histrows.extend(dict(species=sp,state=state,bin_low=a,bin_high=b,cells=int(n)) for a,b,n in zip(edges[:-1],edges[1:],values))
        r.record('UMAP2-cell-counts-'+sp,pd.DataFrame(histrows),'Archived projection counts; not plotted or interpreted as biological evidence')
    r.note('Frozen joint expression coordinates, separate species views. Colors follow original state-family palette; no expression embedding or annotation is refitted.')
    r.panel('B');e=r.read(G+'source_data/actual_state_pair_edges.tsv.gz');hs=s[s.species.eq('human')];ms=s[s.species.eq('mouse')];m=e.pivot(index='human_state',columns='mouse_state',values='human_row_share').reindex(index=hs.state_key,columns=ms.state_key);ax=r.ax(0,kind='matrix',pad=(23,10,2,11));cm=plt.get_cmap('Blues').copy();cm.set_bad('#DEDEDE');order={f:i for i,f in enumerate(FAM)};hs=hs.assign(_family_order=hs.family.map(order)).sort_values(['_family_order','unit','state_key']);ms=ms.assign(_family_order=ms.family.map(order)).sort_values(['_family_order','unit','state_key']);display=m.reindex(index=hs.state_key,columns=ms.state_key);im=ax.imshow(display,aspect='auto',vmin=0,vmax=1,cmap=cm,interpolation='nearest')
    for frame,which in [(hs,'y'),(ms,'x')]:
        display_family=frame.family.map(lambda f:FAM[f] if FAM[f] in list('SREPM') else 'Other');sizes=display_family.groupby(display_family,sort=False).size();p=np.cumsum(sizes)-sizes/2-.5;labs=sizes.index.tolist()
        for boundary in np.cumsum(sizes)[:-1]-.5:
            (ax.axhline if which=='y' else ax.axvline)(boundary,color=GREY,lw=.35)
        if which=='y':ax.set_yticks(p,labs,fontsize=8)
        else:ax.set_xticks(p,labs,fontsize=7);ax.xaxis.tick_top()
    ax.tick_params(length=0);ax.set_xlabel('');ax.set_ylabel('Human → mouse state families');pos=ax.get_position();cb=r.fig.add_axes([pos.x0,pos.y0-4/r.height,.12,1.4/r.height]);r.fig.colorbar(im,cax=cb,orientation='horizontal',ticks=[0,1]);r.colorbars.append((ax,cb));cb.tick_params(labelsize=7,length=1);cb.text(1.2,.5,'Row weight',transform=cb.transAxes,fontsize=7,va='center');r.record('complete-62x49-state-edges',m)
    p=r.read(G+'source_data/program_by_state_library.tsv.gz');p=p[p.n_cells.ge(20)];r.ecdf(1,[(sp[0].upper()+' '+prog,g.mean_log1p_cp10k) for (sp,prog),g in p.groupby(['species','program'])],'Native program score');r.note('62×49 actual edge matrix uses row-normalized weights. Rows and columns are reordered by frozen state family, then unit/state key; the tick label Other spans the remaining original local families without merging cells or weights; the exported matrix preserves its original order. Gray means no estimated outgoing edge; no cross-species lineage identity is inferred from a link. Native-score ECDFs diagnose scale and composition differences; they do not validate cross-species links. Held-out expression excess supplies the separate link diagnostic.')
    r.panel('C');x=r.read(G+'source_data/S2C_route_expression.tsv');pairs=list(x[['human','mouse']].drop_duplicates().itertuples(index=False,name=None));labs=[U.get(a,a)+' × '+U.get(b,b) for a,b in pairs];ax=r.ax(0,labels=labs,kind='point_interval')
    for i,(h,m) in enumerate(pairs):
        g=x[x.human.eq(h)&x.mouse.eq(m)];o=g[g.setting.eq('Original')].iloc[0];z=g[g.setting.eq('Full variants')&g.status.eq('OK')];ax.hlines(i,o.null_q025-o.null_mean,o.null_q975-o.null_mean,color=GREY,lw=2);ax.scatter(z.excess,i+np.linspace(-.16,.16,len(z)),s=5,c=BLUE,lw=0);ax.scatter(o.excess,i,s=12,marker='D',c=INK)
    ax.set(yticks=range(len(pairs)),yticklabels=labs,ylim=(len(pairs)-.4,-.6),xlabel='Held-out expression correlation excess');ax.axvline(0,color=GREY,ls=':');r.record('complete-held-out-refits',x)
    q=r.read('global/planned_route_contrasts.tsv');routes=q.route.drop_duplicates().tolist();rows=[]
    for ro in routes:
        a=q[q.route.eq(ro)];rows.append([short(ro),'Yes' if ((a.contrast=='P_S')&(a.design=='true_source')).any() else '—','Yes' if ((a.contrast=='P_M')&(a.design=='same_library')).any() else '—'])
    table=r.table(1,['Route','Source P–S','Paired P–M'],rows,widths=[.46,.27,.27]);table.texts[1].set_text('Source\nP–S');table.texts[2].set_text('Paired\nP–M');[t.set_y(.05) for t in table.texts[:3]];r.record('original-distance-support',r.read(G+'source_data/S2A_state_geometry.tsv'));r.record('original-neighbor-support',r.read(G+'source_data/S2B_neighbor_all_scales.tsv'));r.note('Twelve route pairs, each with original point and up to12 full refits; non-evaluable fits are retained in data. Original null ranges are reference diagnostics, not CIs for refits. Source-qualified P–S is a restricted design; Liu remains usable under other modules’ observed-state definitions.')
DRAW.update(D01=D01,D02=D02,D03=D03)
