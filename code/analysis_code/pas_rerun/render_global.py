"""Fresh global backgrounds, source sensitivity and positional metrics."""
import json
import numpy as np
import pandas as pd

from render_common import BLUE,ORANGE,INK,GREY,LABELS,empty,zero,matrix


def supplement6(r):
    r.start('Figure S6',150,220)
    full=r.read('global/global_direction_summary.tsv');full=full[full.contrast.eq('P_S')&full.design.eq('true_source')]
    main=full[full.setting.eq('supported50_p08')&full.background.eq('available')];routes=main.route.tolist()
    labels=[LABELS.get(a,a) for a in routes]
    r.label('A',2,2,'Documented-source P–S means and response fractions')
    ax=r.axes((44,15,38,20));zero(ax,False)
    for i,a in enumerate(main.itertuples()):
        ax.plot(a.mean,i,'o',ms=3,color=BLUE if a.species=='human' else ORANGE);ax.plot(a.median,i,'D',ms=3,mfc='white',mec=INK)
    ax.set_yticks(range(len(main)),[f'{lab} ({n})' for lab,n in zip(labels,main.n)],fontsize=6.5);ax.invert_yaxis();ax.set_xlabel('Mean / median ΔDWUI')
    ax=r.axes((93,15,48,20))
    for i,a in enumerate(main.to_dict('records')):
        for j,t in enumerate(['002','005']):
            vals=np.array([a['proximal_'+t],a['n']-a['proximal_'+t]-a['distal_'+t],a['distal_'+t]],float)
            if a['n']==0:continue
            vals/=a['n'];left=0
            for v,color in zip(vals,[BLUE,GREY,ORANGE]):ax.barh(i+(j-.5)*.3,v,left=left,height=.25,color=color);left+=v
    ax.set_yticks([]);ax.invert_yaxis();ax.set_xlim(0,1);ax.set_xlabel('Proximal / below / distal')
    ax.text(.5,1.04,'Upper: 0.02; lower: 0.05',transform=ax.transAxes,ha='center',fontsize=6.5)
    r.label('B',2,45,'Available and common-ortholog means')
    ax=r.axes((44,57,42,20));zero(ax,False)
    for i,route in enumerate(routes):
        q=full[full.route.eq(route)];values=q['mean'].dropna()
        if len(values):ax.plot([values.min(),values.max()],[i,i],color=GREY,lw=1.5)
        for background,marker in [('available','o'),('common_one_to_one_ortholog_all_settings','s')]:
            g=q[q.setting.eq('supported50_p08')&q.background.eq(background)]
            if len(g):ax.plot(g['mean'].iloc[0],i,marker,ms=3,mec=INK,mfc=INK if marker=='o' else 'white')
    ax.set_yticks(range(len(routes)),labels,fontsize=6.5);ax.invert_yaxis();ax.set_xlabel('Mean ΔDWUI')
    r.label('C',91,45,'Settings / backgrounds')
    settings=['supported50_p08','supported50_p09','supported50_p08_mc5','expanded_primary_p08']
    backgrounds=['available','fixed_settings_within_route','common_within_species_all_settings','common_one_to_one_ortholog_all_settings']
    q=full.pivot(index='route',columns=['setting','background'],values='mean').reindex(index=routes,columns=pd.MultiIndex.from_product([settings,backgrounds]))
    ax=r.axes((94,57,47,20));im=ax.imshow(np.ma.masked_invalid(q.to_numpy()),cmap='RdBu_r',vmin=-.05,vmax=.05,aspect='auto')
    ax.set_yticks([]);ax.set_xticks(range(len(q.columns)),list('AFSO')*4,fontsize=6.5);ax.set_title('0.8  |  0.9  |  5 MC  |  Expanded',fontsize=6.5)
    for boundary in [3.5,7.5,11.5]:ax.axvline(boundary,color='white',lw=.6)
    cb=r.axes((94,84,47,1.2));r.fig.colorbar(im,cax=cb,orientation='horizontal',ticks=[-.05,0,.05]);cb.tick_params(labelsize=6.5,pad=0)
    q.to_csv(r.out/'S6-setting-background-matrix.tsv',sep='\t')
    r.label('D',2,92,'Source-library sensitivity')
    source=r.read('global/single_source_library_summary.tsv');source=source[source.contrast.eq('P_S')]
    ax=r.axes((44,104,42,20));zero(ax,False)
    for i,route in enumerate(routes):
        q=source[source.route.eq(route)]
        for j,(background,marker) in enumerate([('available','o'),('fixed_genes_all_sources','s'),('matched_target_libraries_all_sources','^')]):
            g=q[q.background.eq(background)];ax.scatter(g['mean'],np.full(len(g),i)+(j-1)*.16,s=9,marker=marker,c=[BLUE if x=='human' else ORANGE for x in g.species])
    ax.set_yticks(range(len(routes)),labels,fontsize=6.5);ax.invert_yaxis();ax.set_xlabel('Mean ΔDWUI per source')
    r.label('E',91,92,'Tandem UTRs')
    tandem=r.read('tandem/tandem_utr_summary.tsv');q=tandem[tandem.setting.eq('supported50_p08')&tandem.background.eq('available')&tandem.contrast.eq('P_S')].set_index('route').reindex(routes)
    ax=r.axes((94,104,47,20));zero(ax,False)
    ax.scatter(q.mean_delta_nt,range(len(q)),s=12,color=INK);ax.scatter(q.median_delta_nt,range(len(q)),s=12,facecolors='white',edgecolors=INK,marker='D');ax.set_yticks([]);ax.invert_yaxis();ax.set_xlabel('Mean / median change (nt)')
    r.label('F',2,138,'Four positional summaries')
    metrics=r.read('metrics/metric_summaries.tsv');q=metrics[metrics.level.eq('route')&metrics.setting.eq('supported50_p08')&metrics.contrast.eq('P_S')&metrics.design.eq('true_source')&metrics.threshold.eq(.02)]
    ax=r.axes((44,150,42,20));zero(ax,False)
    for j,(metric,color,marker) in enumerate([('DWUI',BLUE,'o'),('RWUI',ORANGE,'s'),('PDUI',GREY,'^'),('PPUI','#318F9C','D')]):
        g=q[q.metric.eq(metric)].set_index('route').reindex(routes)
        ax.scatter(g['mean']*(-1 if metric=='PPUI' else 1),np.arange(len(g))+(j-1.5)*.13,s=8,marker=marker,color=color,label='−PPUI' if metric=='PPUI' else metric)
    ax.set_yticks(range(len(routes)),labels,fontsize=6.5);ax.invert_yaxis();ax.set_xlabel('Mean index change');ax.legend(frameon=False,ncol=4,fontsize=6.5,loc='lower center',bbox_to_anchor=(.3,1.03),columnspacing=.6,handletextpad=.3)
    r.label('G',91,138,'Redistribution')
    effects=r.read('state_vectors/route_effects.primary.tsv.gz')
    ax=r.axes((103,150,38,20));summaries=[]
    for co,design,marker,color in [('P_S','true_source','o',BLUE),('P_M','same_library','^',ORANGE)]:
        q=effects[effects.contrast.eq(co)&effects.design.eq(design)].copy()
        q['PAS_count']=q.PAS_ids.map(lambda a:len(json.loads(a)));q=q[q.PAS_count.ge(3)]
        hidden=int((q.delta_DWUI.abs().lt(.02)&q.MPRO.ge(.1)).sum())
        ax.scatter(q.delta_DWUI,q.MPRO,s=1,alpha=.2,color=color,marker=marker,rasterized=True)
        summaries.append(dict(contrast=co,design=design,gene_route_records=len(q),small_net_MPRO01=hidden))
    ax.axhline(.1,color=GREY,ls=':',lw=.6)
    for x in [-.02,.02]:ax.axvline(x,color=GREY,ls=':',lw=.6)
    ax.set(xlabel='ΔDWUI',ylabel='MPRO')
    for i,(a,color) in enumerate(zip(summaries,[BLUE,ORANGE])):
        ax.text(1,1.22-.18*i,f"{a['contrast'].replace('_','–')}: {a['small_net_MPRO01']}/{a['gene_route_records']}",ha='right',va='bottom',transform=ax.transAxes,fontsize=6.5,color=color)
    pd.DataFrame(summaries).to_csv(r.out/'S6-redistribution-denominators.tsv',sep='\t',index=False)
    r.label('H',2,186,'Tracked examples · current PAS vectors')
    axes=r.read('state_vectors/PAS_axes.tsv.gz');axes=axes[axes.setting.eq('supported50_p08')]
    for j,name in enumerate(['COMMD5','CNN3']):
        ax=r.axes((35+j*67,197,43,12));g=effects[effects.gene_symbol.str.upper().eq(name)&effects.route.eq('178325_0809')&effects.contrast.eq('P_S')&effects.design.eq('true_source')]
        if g.empty:empty(ax,name+': not measurable');continue
        a=g.iloc[0];p=axes[axes.study_id.eq(a.study_id)].drop_duplicates('atlas_pas_id').set_index('atlas_pas_id').loc[json.loads(a.PAS_ids)].position
        for field,color,marker,label in [('reference_usage',GREY,'o','S'),('target_usage',BLUE,'s','P')]:ax.plot(p,np.array(json.loads(a[field])),marker=marker,ms=2,color=color,lw=.5,label=label)
        ax.set_title(name,fontsize=7);ax.set_xlim(0,1);ax.set_ylim(0,1);ax.set_ylabel('Usage');ax.set_xlabel('Full-primary position')
        ax.legend(frameon=False,ncol=2,loc='upper right',fontsize=6.5,handlelength=.9,columnspacing=.5)
    r.save('Fig/SuppFig.S5.pdf','New PAS effects define every background, source sensitivity, tandem-UTR subset and metric. Setting ranges are descriptive, not confidence intervals. Composition examples retain their previous identities and use current vectors even if they no longer meet the old selection rationale. No fixed 710-ortholog denominator or old effect size is carried forward.')
