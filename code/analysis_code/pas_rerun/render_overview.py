"""Preserve frozen GEX panels A/B and rebuild every PAS display in Figure 1."""
import numpy as np
import pandas as pd

from render_common import INK,BLUE,ORANGE,GREY,ROUTES,LABELS,embedding,matrix,zero,join_frozen_scores


def figure1(r):
    from render_common import HOUSE_MAP,STATE_COLORS,HOUSE_INK
    r.main(1);original=r.resource('original_figure1');r.panels.extend(['A','B'])
    cells=pd.read_csv(r.resource('frozen_cell_embedding'),sep='\t')
    scores=r.read('figure_scores/metacell_scores.tsv.gz');r.resource('figure_scores_code_audit')
    cells=join_frozen_scores(cells,scores);assert len(cells)==551201
    r.label('C',0,0)
    for i,sp in enumerate(['human','mouse']):
        ax=r.axes((3+i*40,171,34,30));q=cells[cells.species_key.eq(sp)]
        im=embedding(ax,q,'centered_DWUI',.06);im.set_cmap(HOUSE_MAP)
        r.text(20+i*40,165.4,sp.title(),ha='center')
        for state in 'SREPM':
            col='stateflow_code'
            if col in q:
                z=q[q[col].eq(state)]
                if len(z):ax.text(z.UMAP1.median(),z.UMAP2.median(),state,ha='center',weight='bold',fontsize=8)
    r.text(3,204,'Gene-centered DWUI\nwithin study')
    cb=r.axes((46,205,25,1.4));r.fig.colorbar(im,cax=cb,orientation='horizontal',ticks=[-.06,0,.06],extend='both')
    r.label('D',0,0)
    nodes=r.read('profiles/all_route_node_positions.tsv.gz');nodes=nodes[nodes.setting.eq('supported50_p08')&nodes.depth.eq(20)]
    wide=nodes.pivot(index=['species','route','gene_id'],columns='node',values='position_index').reindex(columns=list('SREPM'))
    parts=[];denominators=[]
    for (sp,route),g in wide.groupby(level=['species','route']):
        observed=[state for state in 'SREPM' if g[state].notna().any()]
        if 'S' not in observed:continue
        measured=g[observed].dropna();mean=measured.sub(measured.S,axis=0).mean()
        parts.append(dict(route=route,**mean.to_dict()))
        denominators.append(dict(species=sp,route=route,common_genes=len(measured),states=''.join(observed)))
    order=[v for v in ROUTES if v in {v['route'] for v in parts}]
    values=pd.DataFrame(parts).set_index('route').reindex(index=order,columns=list('REPM'))
    r.display('D-state-values',values.reset_index());r.display('D-state-denominators',pd.DataFrame(denominators))
    ax=r.axes((116,126.5,41,45));limit=max(.04,np.ceil(np.nanmax(abs(values.to_numpy()))/.02)*.02)
    for i,route in enumerate(order):
        for j,state in enumerate('REPM'):
            ax.barh(i,.86,left=j-.43,height=.68,color='#EDF0F2')
            value=values.loc[route,state]
            if np.isfinite(value):
                ax.plot([j,j],[i-.34,i+.34],color='#94A1AA',lw=.35)
                ax.barh(i,value/limit*.43,left=j,height=.58,color=HOUSE_MAP(.9 if value>0 else .1))
            else:ax.text(j,i,'NA',fontsize=7.2,color=GREY,ha='center',va='center')
    ax.set(xlim=(-.5,3.5),ylim=(len(order)-.5,-.5));ax.set_yticks(range(len(order)),[LABELS[v] for v in order],fontsize=7.2)
    ax.set_xticks(range(4),list('REPM'));ax.xaxis.tick_top();ax.tick_params(length=0)
    for a in ax.spines.values():a.set_visible(False)
    for t in ax.get_xticklabels():t.set_color(STATE_COLORS[t.get_text()])
    r.text(118,114,'Mean ΔDWUI vs observed S',ha='center')
    r.text(119,120,f'Proximal  −{limit:.2f}    0    +{limit:.2f}  Distal',ha='center',fontsize=7.2)
    r.text(119,177,'Documented-source P–S',ha='center')
    source=r.read('global/global_direction_summary.tsv');source=source[source.setting.eq('supported50_p08')&source.contrast.eq('P_S')&source.design.eq('true_source')]
    stable=r.read('global/direction_stability.tsv');stable=stable[stable.contrast.eq('P_S')&stable.design.eq('true_source')]
    order=[v for v in ROUTES if v in set(source.route)];ax=r.axes((121,183,34,16))
    for i,route in enumerate(order):
        a=stable[stable.route.eq(route)].iloc[0];ax.plot([a.mean_min,a.mean_max],[i,i],color=GREY,lw=1.7)
        for j,background in enumerate(['available','fixed_settings_within_route']):
            q=source[source.route.eq(route)&source.background.eq(background)];assert len(q)==1
            ax.plot(q['mean'].iloc[0],i,['o','|'][j],ms=3.3,mec=HOUSE_INK,mfc='white' if j==0 else HOUSE_INK)
    ax.axvline(0,color=GREY,lw=.5,ls=':');ax.set_yticks(range(len(order)),[LABELS[v] for v in order],fontsize=7.2)
    ax.set_ylim(len(order)-.5,-.5);ax.set_xticks([-.04,0,.02]);ax.set_xlabel('Mean ΔDWUI',labelpad=1)
    r.text(81,208,'○ Available   | Fixed genes   — 16-scenario range',fontsize=7.2)
    r.display('D-primary-endpoints',source);r.display('D-sensitivity',stable)
    r.save('Fig/Figure.1.pdf','A/B preserved from the approved apadetect layout. Current gene-centered PAS scores and state means replace C/D. Grey bars in D are ranges across the 16 current settings/backgrounds, not confidence intervals or BAM counting rules.',reuse_figure1=original)
