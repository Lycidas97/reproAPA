"""Cis model denominators, coefficients and complete prediction sensitivity."""
import numpy as np
from matplotlib.ticker import MaxNLocator,FormatStrFormatter

from render_common import BLUE,ORANGE,INK,GREY,empty,zero
from render_cross import interval_rows


def supplement18(r):
    r.start('Figure S18',150,214)
    models=r.read('cis/recomputed_models.tsv');primary=models[models.primary]
    r.label('A',2,2,'Model populations and matched S/R/B controls')
    ax=r.axes((36,15,65,25));ax.barh(range(len(primary)),primary.n,color=GREY,label='All model genes');ax.barh(range(len(primary)),primary.positive,color=INK,label='Aligned responses')
    labels=[('Hs→Mm' if a.discovery=='human' else 'Mm→Hs')+' '+a.contrast for a in primary.itertuples()]
    ax.set_yticks(range(len(primary)),labels);ax.invert_yaxis();ax.set_xlabel('Genes');ax.legend(frameon=False,ncol=2,loc='upper left',bbox_to_anchor=(0,-.5))
    for i,a in enumerate(primary.itertuples()):ax.text(a.n,i,f' {a.positive}/{a.n}',va='center',fontsize=6.5)
    counts=r.read('cis/descriptive_class_counts.tsv').set_index(['contrast','discovery'])
    ax=r.axes((113,15,30,25))
    for j,(name,color) in enumerate([('Shared',BLUE),('Recurrent-only',ORANGE),('Matched background',GREY)]):
        values=[counts.loc[(a.contrast,a.discovery),name] for a in primary.itertuples()]
        ax.barh(np.arange(len(primary))+(j-1)*.22,values,height=.2,color=color,label=['S','R','B'][j])
    ax.invert_yaxis();ax.set_yticks([]);ax.set_xlabel('Unique matched genes');ax.legend(frameon=False,ncol=3,loc='upper left',bbox_to_anchor=(-.15,-.5),columnspacing=.6,handlelength=.8)
    r.read('cis/feature_distribution_summary.tsv');r.read('cis/feature_correlations.tsv');r.read('cis/descriptive_matched_features.tsv.gz')
    coef=r.read('cis/recomputed_coefficients.tsv');coef=coef[coef.layer.eq('supported50')&coef.tau.eq(.02)]
    terms=['driver_conserved','frac_conserved_discovery_pas','A_g','phastcons_driver_mean','motif_rank_driver_mean']
    r.label('B',2,62,'Primary standardized coefficients · profile 95% intervals')
    for j,(sp,co) in enumerate([('human','P-S'),('human','P-M'),('mouse','P-S'),('mouse','P-M')]):
        q=coef[coef.discovery.eq(sp)&coef.contrast.eq(co)].set_index('term').reindex(terms).reset_index()
        ax=r.axes((43+j*26,75,20,32));interval_rows(ax,q,['Driver conserved','PAS conservation','Architecture','phastCons','Motif rank'] if j==0 else ['']*5,'beta','ci_low','ci_high')
        ax.set_title(('H→M' if sp=='human' else 'M→H')+' '+co,fontsize=6.5);ax.set_xlabel('β')
        ax.xaxis.set_major_locator(MaxNLocator(2));ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    r.label('C',2,121,'All 32 configurations · paired AUC gain and joint q')
    cv=r.read('cis/recomputed_cv.tsv');cv=cv[cv.model.isin(['M1','All cis'])]
    keys=['layer','contrast','discovery','tau'];average=cv.groupby(keys+['model']).auc.mean().unstack('model')
    average['gain']=average.get('All cis',np.nan)-average.get('M1',np.nan);grid=models.merge(average[['gain']].reset_index(),on=keys,how='left',validate='one_to_one')
    ax=r.axes((39,135,100,22));contexts=[('human','P-S'),('human','P-M'),('mouse','P-S'),('mouse','P-M')]
    for a in grid.itertuples():
        y=contexts.index((a.discovery,a.contrast));x=[.01,.02,.03,.05].index(a.tau)+(0 if a.layer=='supported50' else 4)
        if np.isfinite(a.gain):
            color=BLUE if a.gain<0 else ORANGE;ax.scatter(x,y,s=8+400*abs(a.gain),facecolors=color if a.q_boot_planned_32<.05 else 'white',edgecolors=color,lw=.7)
        else:ax.plot(x,y,'x',ms=3,color=GREY)
    ax.set_xticks(range(8),[f'{x:g}' for x in [.01,.02,.03,.05]*2]);ax.set_yticks(range(4),[('Hs→Mm' if sp=='human' else 'Mm→Hs')+' '+co for sp,co in contexts]);ax.invert_yaxis();ax.axvline(3.5,color=GREY,lw=.5);ax.set_xlabel('Reference-supported                      Expanded primary')
    ax.set_ylim(3.5,-.5)
    r.label('D',2,173,'Prediction gains · fixed-fold gene-block 95% intervals')
    ci=r.read('cis/paired_cv_uncertainty.tsv');families=['M0','PAS','Architecture','Sequence','All cis']
    for j,metric in enumerate(['auc','logloss']):
        ax=r.axes((34+j*64,186,47,18));q=ci[ci.metric.eq(metric)&ci.model.isin(families)]
        for k,(sp,co) in enumerate(contexts):
            f=q[q.discovery.eq(sp)&q.contrast.eq(co)]
            for a in f.itertuples():
                y=families.index(a.model)+(k-1.5)*.18
                # Display benefit consistently: M1−M0 for recurrence, others−M1;
                # lower log loss is a gain, whereas higher AUC is a gain.
                sign=(-1 if a.model=='M0' else 1)*(1 if metric=='auc' else -1)
                if np.isfinite(a.delta):
                    lo,hi=sorted([sign*a.lo,sign*a.hi]);ax.plot([lo,hi],[y,y],color=BLUE if sp=='human' else ORANGE,lw=.6)
                    ax.plot(sign*a.delta,y,'o' if co=='P-S' else 's',ms=2,color=BLUE if sp=='human' else ORANGE)
        zero(ax,False);ax.set_yticks(range(5),['Recurrence','PAS','Architecture','Sequence','Combined cis'] if j==0 else ['']*5,fontsize=6.5);ax.invert_yaxis();ax.set_xlabel('AUC gain' if metric=='auc' else 'Log-loss decrease')
    r.save('Fig/SuppFig.S17.pdf','New PAS sequence, architecture and driver features define every model input. A reports actual model/response denominators and unique genes in the fixed matched Shared, Recurrent-only and Background sets. Matched control weights and feature-wise observed distributions/correlations are exported, with reused genes counted once for population size. B uses primary profile likelihood intervals. C retains all 32 planned positions, with non-estimable predictions marked. D preserves paired folds, shared gene bootstrap draws and reused predictions; intervals condition on fitted models.')
