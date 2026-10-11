"""Native polyApipe versus frozen counts, with complete missing-data denominators."""
import numpy as np
import pandas as pd

from render_common import BLUE,ORANGE,INK,GREY,empty,zero,LABELS


def endpoint(frame):
    return frame[(frame.contrast.eq('P_S')&frame.design.eq('true_source'))|(frame.contrast.eq('P_M')&frame.design.eq('same_library'))]


def supplement5(r):
    r.start('Figure S5',150,210)
    available=r.read('caller_comparison/caller_availability.tsv')
    r.label('A',2,2,'Full-catalogue PAS matching')
    distance=r.read('caller_comparison/bidirectional_PAS_distances.tsv');sids=distance.study_id.drop_duplicates().tolist()
    ax=r.axes((22,17,49,35))
    for j,direction in enumerate(['native_all_to_legacy_primary','legacy_primary_to_native_all']):
        for tau,marker in [(25,'o'),(50,'s'),(100,'^')]:
            q=distance[distance.direction.eq(direction)&distance.tolerance_nt.eq(tau)].set_index('study_id').reindex(sids)
            ax.scatter(q.matched/q.full_catalogue_denominator,np.arange(len(q))+(j-.5)*.3,s=11,marker=marker,
                facecolors=BLUE if j==0 else 'white',edgecolors=BLUE if j==0 else ORANGE,
                label=f'{"N→O" if j==0 else "O→N"} ≤{tau}')
    ax.set_yticks(range(len(sids)),[s.replace('GSE','') for s in sids],fontsize=6.5);ax.invert_yaxis();ax.set_xlim(0,1);ax.set_xlabel('Matched / full catalogue')
    ax.legend(frameon=False,ncol=3,fontsize=6.5,loc='upper left',bbox_to_anchor=(-.25,-.3),columnspacing=.7,handletextpad=.2)
    stats=endpoint(r.read('caller_comparison/effect_agreement_summary.tsv'));stats=stats[stats.contributors.eq('fixed_original_pairs_complete')]
    main=stats[stats.view.eq('native_terminal_sites')&stats.comparison.eq('frozen_main')]
    common=stats[stats.view.eq('matched_sites_common_metacells')&stats.comparison.eq('repro_same_PAS_same_contributors')]
    keys=['study_id','route','contrast','design'];joined=main.merge(common,on=keys,how='outer',suffixes=('_native','_common'),validate='one_to_one')
    labels=[LABELS.get(a.route,a.route).replace('2020','').replace('2019','').replace('2022','').replace('2023','')+' '+a.contrast[-1] for a in joined.itertuples()]
    y=np.arange(len(joined))
    r.label('B',79,2,'Fixed original-gene background')
    ax=r.axes((109,17,34,35))
    ax.barh(y,joined.frozen_main_gene_background_native,color=GREY,label='Old background')
    ax.barh(y-.15,joined.jointly_measurable_genes_native,height=.28,color=BLUE,label='Native fixed pairs')
    ax.barh(y+.15,joined.jointly_measurable_genes_common,height=.28,color=ORANGE,label='Common PAS / MC')
    ax.set_yticks(y,labels,fontsize=6.5);ax.invert_yaxis();ax.set_xlabel('Genes');ax.legend(frameon=False,fontsize=6.5,loc='upper left',bbox_to_anchor=(-.15,-.3))
    for j,(letter,view,title) in enumerate([('C','native','Native axes / fixed pairs'),('D','common','Common PASs, MCs and pairs')]):
        r.label(letter,2+j*77,77,title);ax=r.axes((31+j*77,91,35,35));zero(ax,False)
        for i,a in joined.iterrows():
            n=a.get('jointly_measurable_genes_'+view,np.nan)
            if not np.isfinite(n) or not n:
                ax.text(.02,i,'No common genes',transform=ax.get_yaxis_transform(),fontsize=6.5,color=GREY);continue
            ax.plot(a['spearman_'+view],i-.13,'o',ms=3,color=BLUE,label='Spearman ρ' if i==0 else None)
            ax.plot(a['direction_agree_'+view]/n,i+.13,'s',ms=3,color=ORANGE,label='Direction agreement' if i==0 else None)
        ax.set_yticks(y,labels,fontsize=6.5);ax.invert_yaxis();ax.set_xlim(-1,1);ax.set_xlabel('ρ / agreeing fraction')
        if j==0:ax.legend(frameon=False,ncol=2,fontsize=6.5,loc='upper left',bbox_to_anchor=(-.55,-.3),columnspacing=.7)
    joined.to_csv(r.out/'S5-fixed-pair-native-common-comparison.tsv',sep='\t',index=False)
    r.label('E',2,146,'All old recurrent gene–route records')
    old=r.read('caller_comparison/complete_legacy_recurrent_gene_migration.tsv.gz')
    old=old.drop_duplicates(['species','contrast','gene_id','legacy_threshold','route','study_id'])
    categories=['same_direction_above_threshold','same_direction_below_threshold','zero','reversed','not_measurable']
    table=[]
    for (sp,co,tau),g in old.groupby(['species','contrast','legacy_threshold']):
        key='migration_'+str(float(tau))
        counts={k:int(g[key].eq(k).sum()) for k in categories};counts['not_measurable']+=int(g[key].isna().sum())
        assert sum(counts.values())==len(g)
        table.append(dict(species=sp,contrast=co,threshold=tau,legacy_genes=g.gene_id.nunique(),gene_route_records=len(g),**counts))
    table=pd.DataFrame(table);table.to_csv(r.out/'S5-legacy-roster-outcomes.tsv',sep='\t',index=False)
    ax=r.axes((26,161,43,26));left=np.zeros(len(table));colors=[BLUE,'#97BDD9','#D9DDE1',ORANGE,GREY]
    for name,label,color in zip(categories,['Retained','Below τ','Zero','Reversed','Unavailable'],colors):
        vals=table[name]/table.gene_route_records;ax.barh(range(len(table)),vals,left=left,color=color,label=label);left+=vals
    ax.set_yticks(range(len(table)),[f'{a.species[0].upper()} {a.contrast[-1]} {a.threshold:g}' for a in table.itertuples()],fontsize=6.5);ax.invert_yaxis();ax.set_xlim(0,1);ax.set_xlabel('Fraction of old roster records')
    ax.legend(frameon=False,ncol=3,fontsize=6.5,loc='upper left',bbox_to_anchor=(-.4,-.3),columnspacing=.7,handlelength=.9)
    r.label('F',79,146,'Named examples: measured routes')
    examples=endpoint(r.read('caller_comparison/main_example_caller_support.tsv.gz'))
    examples=examples[examples.contributors.eq('fixed_original_pairs_complete')].copy();examples['identity']=examples.gene_symbol.str.upper()
    names=sorted(examples.identity.unique());support=[]
    views=['native_terminal_sites','matched_sites_common_metacells']
    for name in names:
        for view in views:
            q=examples[examples.identity.eq(name)&examples.view.eq(view)]
            assert not q.duplicated(['study_id','route','contrast','design']).any()
            support.append(dict(gene_symbol=name,view=view,planned_routes=len(q),measured_routes=int(q.delta_primary_axis.notna().sum())))
    support=pd.DataFrame(support);support.to_csv(r.out/'S5-named-example-support.tsv',sep='\t',index=False)
    ax=r.axes((105,161,37,32))
    if names:
        values=support.pivot(index='gene_symbol',columns='view',values='measured_routes').reindex(index=names,columns=views)
        ax.imshow(values.to_numpy(),aspect='auto',cmap='Blues',vmin=0,vmax=max(1,values.to_numpy().max()))
        ax.set_yticks(range(len(names)),names,fontsize=6.5);ax.set_xticks([0,1],['Native','Common'],fontsize=6.5)
        for i in range(len(names)):
            for j in range(2):
                value=values.iloc[i,j]
                ax.text(j,i,str(value),ha='center',va='center',fontsize=6.5,
                    color='white' if value>values.to_numpy().max()*.55 else INK)
    else:empty(ax,'No named gene identities available')
    scape=available[available.caller.eq('SCAPE')]
    r.fig.text(80/150,1-205/210,f'SCAPE: {len(scape)} studies reserved; separate run.',fontsize=6.5)
    r.save('Fig/SuppFig.tool-comparison.reserved.pdf','Current native polyApipe and frozen old-pipeline counts are compared on fixed original endpoint pairs. A uses full catalogues, including unassigned and unmatched sites. B retains old backgrounds. C/D use only jointly measurable genes and retain their denominators in B. E shows every old recurrent gene-route record at its original threshold, including unavailable records; these fractions are not gene-level recurrence retention rates. F counts qualified fixed-pair route/contrast records for each prespecified named identity. All source effects, eligibility, threshold classes and legacy gene rosters are exported. SCAPE remains explicitly reserved for its separate seven-study comparison.')
