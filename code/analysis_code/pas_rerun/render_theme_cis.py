from render_theme import *
DRAW={}
CO=[('human','P-S'),('human','P-M'),('mouse','P-S'),('mouse','P-M')]
FEATURES=['driver_conserved','frac_conserved_discovery_pas','A_g','phastcons_driver_mean','motif_rank_driver_mean']
FL=['Driver','PAS fraction','A_g','phastCons','Motif rank']
def counts(r,k):
    q=r.read('cis/descriptive_class_counts.tsv');f=q.set_index(['discovery','contrast']).reindex(CO);a=f[['model_Y1','model_n']].copy();a['model_n']-=a.model_Y1;a.columns=['Y = 1','Y = 0'];a.index=CL;ax=r.bars(k,a,'Model genes');legend(ax,handles=[Line2D([],[],c=BLUE,lw=3,label='Same direction'),Line2D([],[],c=ORANGE,lw=3,label='Other')],ncol=2);r.note('Model Y=1 is validation_same_direction: same effect direction in the validation species relative to discovery. It is not the Shared descriptive class; Y=0 includes all other eligible model genes.');r.record('distinct-descriptive-cohorts',q)
def D19(r):
    r.panel('A');ax=r.ax(0,kind='structure');ax.set(xlim=(0,1),ylim=(0,1));ax._concept_schematic=True
    ax.text(.22,.98,'Architecture',ha='center',va='top',fontsize=8);ax.text(.76,.98,'Sequence',ha='center',va='top',fontsize=8)
    # Definition schematic, not fabricated gene measurements: matched PASs and local windows.
    for yy,sites in [(.70,[.06,.17,.29,.43]),(.48,[.06,.22,.29,.43])]:
        ax.plot([.04,.45],[yy,yy],color=GREY,lw=.7)
        for xx in sites:ax.plot(xx,yy,marker='|',ms=7,color=BLUE if xx in [.06,.29,.43] else GREY)
    for xx in [.06,.29,.43]:ax.plot([xx,xx],[.52,.66],ls=':',color=BLUE,lw=.6)
    ax.text(.04,.34,'Driver 0–2 · PAS fraction',fontsize=7,va='top');ax.text(.04,.16,'A_g = 2m / (n_H + n_M)',fontsize=7,va='top')
    ax.plot([.59,.93],[.70,.70],color=GREY,lw=.8);ax.add_patch(Rectangle((.70,.63),.10,.14,facecolor='#E5DDEB',edgecolor='none'));ax.plot(.75,.70,marker='|',ms=9,c=INK)
    ax.annotate('',xy=(.63,.46),xytext=(.715,.62),arrowprops=dict(arrowstyle='-',color=GREY,lw=.6));ax.annotate('',xy=(.87,.46),xytext=(.785,.62),arrowprops=dict(arrowstyle='-',color=GREY,lw=.6))
    ax.plot([.59,.92],[.43,.43],color=INK,lw=.7)
    for xx in [.61,.65,.69,.73,.77,.81,.85,.89]:ax.plot(xx,.43,marker='|',ms=4,c='#8766A4')
    ax.text(.76,.28,'Local phastCons',ha='center',fontsize=7);ax.text(.76,.11,'Motif rank',ha='center',fontsize=7)
    r.record('definition-0',pd.DataFrame({'label':['Architecture\nDriver · PAS fraction · A_g','Sequence\nphastCons · motif rank'],'value':['','']}),'Concept schematic; no invented gene or sequence observations')
    counts(r,1);r.note('Definition schematic: driver = number of conserved gain/loss PASs (0–2); PAS fraction = conserved discovery PASs / all discovery PASs; A_g = 2 matched PASs / (human PASs + mouse PASs). The schematic has no empirical coordinates. phastCons measures local sequence conservation. Motif rank is the original sequence rank; neither is evidence of selection. Model and descriptive matched cohorts have distinct denominators.')
    r.panel('B');x=r.read('cis/descriptive_matched_features.tsv.gz');x=x[x.layer.eq('supported50')&x.tau.eq(.02)];classes=['Shared','Recurrent-only','Matched background']
    for j,key in enumerate(['A_g','phastcons_driver_mean']):
        ax=r.ax(0,sub=j,n=2,kind='distribution',labels=CL,pad=(25,7,2,11))
        for i,(sp,co) in enumerate(CO):
            for c,cls in enumerate(classes):
                z=x[x.discovery.eq(sp)&x.contrast.eq(co)&x['class'].eq(cls)].dropna(subset=[key]).sort_values(key)
                if z.empty:
                    ax.text(.97,i+(c-1)*.23,'n = 0',transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=7);continue
                v=z[key].to_numpy();w=z.match_weight.to_numpy();y=np.cumsum(w)/w.sum();yy=i+(c-1)*.27
                edges=np.linspace(0,1,13);freq,_=np.histogram(v,edges,weights=w);assert np.isclose(freq.sum(),w.sum())
                ridge=.19*freq/freq.max();ax.fill_between(edges,yy,yy-np.r_[ridge,ridge[-1]],step='post',color=PAL[c],alpha=.32,lw=.45,edgecolor=PAL[c])
                ax.scatter(v,np.full(len(v),yy+.035),s=3,c=PAL[c],alpha=.35,lw=.45,marker='|',rasterized=True)
                qq=np.interp([.25,.5,.75],y,v);ax.plot(qq[[0,2]],[yy]*2,color=PAL[c],lw=1.6);ax.plot(qq[1],yy,'o',ms=2.2,mfc='white',mec=PAL[c],mew=.65)
        ax.set(yticks=range(4),yticklabels=CL,ylim=(3.5,-.55),xlim=(-.03,1.03),xticks=[0,.5,1],xlabel=key if j==0 else 'Local phastCons');ax.legend(handles=[Line2D([],[],c=PAL[i],lw=2,label=c) for i,c in enumerate(['Shared','Recurrent','Background'])],frameon=False,ncol=3,loc='lower left',bbox_to_anchor=(0,1),columnspacing=.4,handlelength=.7)
    r.record('weighted-descriptive-distributions',x);r.note('Weighted empirical step silhouettes use12 equal bins on0–1 and the original match_weight. Each row is normalized to its own peak, so height is not sample size. Rugs retain all finite source values without shifting the measured coordinate; thick segments and hollow points show the same weighted IQR and median as v12. All four discovery/contrast backgrounds, including empty cells, are retained; these are descriptive distributions, not significance tests.')
    r.panel('C');q=r.read('cis/recomputed_coefficients.tsv');q=q[q.layer.eq('supported50')&q.tau.eq(.02)&q.term.isin(FEATURES)].copy();labels=[context(sp,co)+' · '+FL[j] for sp,co in CO for j in range(5)];q['_row']=[CO.index((a.discovery,a.contrast))*5+FEATURES.index(a.term) for a in q.itertuples()];r.dots(0,q,'beta',labels,lo='ci_low',hi='ci_high',xlabel='Adjusted logistic coefficient β');r.note('Intervals are shown only where finite in the source; absent coefficient intervals are not zero-width uncertainty estimates.')
    c=r.read('cis/feature_correlations.tsv');pairs=[(a,b) for i,a in enumerate(FEATURES) for b in FEATURES[i+1:]];rows=[]
    for a,b in pairs:
        z=c[c.feature1.eq(a)&c.feature2.eq(b)].set_index(['discovery','contrast']).reindex(CO);rows.append([FL[FEATURES.index(a)].replace('PAS fraction','PAS%').replace('phastCons','Cons.').replace('Motif rank','Motif')+'–'+FL[FEATURES.index(b)].replace('PAS fraction','PAS%').replace('phastCons','Cons.').replace('Motif rank','Motif'),*['—' if not np.isfinite(v) else f'{v:.2f}' for v in z.spearman]])
    ax=r.ax(1,kind='point_interval',labels=[a[0] for a in rows],pad=(12,8,2,10))
    for j,co in enumerate(CO):
        vals=np.array([float(a[j+1]) if a[j+1]!='—' else np.nan for a in rows]);ax.scatter(vals,np.arange(10)+(j-1.5)*.14,s=13,marker=['o','s','^','D'][j],c=PAL[j],label=context(*co))
    ax.set(yticks=range(10),yticklabels=[a[0] for a in rows],ylim=(9.5,-.5),xlim=(-1,1),xlabel='Feature-pair Spearman ρ');ax.axvline(0,color=GREY,ls=':');legend(ax,loc='lower right',bbox_to_anchor=(1,1.01),ncol=2);r.record('all-feature-correlations',c);r.note('Feature pairs are named directly; PAS%=conserved discovery PAS fraction, Cons.=local phastCons, Motif=motif rank. Four original model contexts are separate; correlations are descriptive, not adjusted effects.')

def D20(r):
    r.panel('A');r.note('M0 measurement; M1 adds recurrence; Architecture, Sequence and All cis are additions to M1. The model definition is in the caption.');counts(r,1);r.note('M1 already includes recurrence score. Architecture and Sequence are separate additions to M1; their gains must not be summed.')
    r.panel('B');q=r.read('cis/paired_cv_uncertainty.tsv');q=q[q.layer.eq('supported50')&q.tau.eq(.02)&q.model.ne('M1')].copy();models=['M0','PAS','Architecture','Sequence','All cis'];names=['M1−M0','PAS−M1','Architecture−M1','Sequence−M1','All cis−M1'];labs=[context(sp,co)+' '+name for sp,co in CO for name in names]
    q['_row']=[CO.index((a.discovery,a.contrast))*5+models.index(a.model) for a in q.itertuples()]
    for j,metric in enumerate(['auc','logloss']):
        z=q[q.metric.eq(metric)].copy();reverse=z.model.eq('M0') ^ (metric=='logloss');old=z.copy();z.loc[reverse,'delta']=-old.loc[reverse,'delta'];z.loc[reverse,'lo']=-old.loc[reverse,'hi'];z.loc[reverse,'hi']=-old.loc[reverse,'lo'];ax=r.dots(0,z,'delta',labs if j==0 else ['']*20,lo='lo',hi='hi',sub=j,n=2,xlabel='AUC gain' if j==0 else 'Log-loss reduction',pad=(32 if j==0 else 2,5,1,11));ax.set_title('Higher = better',fontsize=7,pad=3)
    cv=r.read('cis/recomputed_cv.tsv');cv=cv[cv.layer.eq('supported50')&cv.tau.eq(.02)];w=cv.pivot(index=['discovery','contrast','repeat'],columns='model',values=['auc','logloss']);z=pd.DataFrame({'auc_gain':w[('auc','All cis')]-w[('auc','M1')],'loss_reduction':w[('logloss','M1')]-w[('logloss','All cis')]}).reset_index();ax=r.ax(1,kind='scatter',pad=(13,3,2,10))
    for j,co in enumerate(CO):
        g=z[z.discovery.eq(co[0])&z.contrast.eq(co[1])];ax.scatter(g.auc_gain,g.loss_reduction,s=13,c=PAL[j],marker=['o','s','^','D'][j],label=context(*co),alpha=.75)
    zeros(ax,True);ax.set(xlabel='All-cis AUC gain',ylabel='All-cis log-loss reduction',xticks=[-.04,0,.04]);legend(ax,loc='best',ncol=1);r.record('joint-repeat-gains',z);r.note('Rows identify each model and baseline: recurrence is M1−M0; PAS, Architecture, Sequence and All cis are separate additions to M1. Points and intervals summarize the mean of ten within-repeat paired OOF metric differences. A single gene draw is reused across the ten fixed-fold repeats; models are not refitted. Log-loss signs and interval ends are reversed to display reduction. Joint plot shows the same ten original repeat-level All-cis−M1 gains used in the mean, not independent datasets. All-cis joint q is not assigned to feature-group gains.')
    r.panel('C');x=r.read('cis/recomputed_cv.tsv');x=x[x.model.isin(['M1','All cis'])];g=x.groupby(['layer','contrast','discovery','tau','model'])[['auc','logloss']].mean().unstack('model');r.record('all-32-setting-means',g);cols=['supported50','expanded_primary']
    for j,metric in enumerate(['auc','logloss']):
        ax=r.ax(0,sub=j,n=2,kind='line',pad=(13,13,2,11));dif=g[(metric,'All cis')]-g[(metric,'M1')];dif=dif*(1 if metric=='auc' else -1);z=dif.rename('gain').reset_index()
        for (layer,co,sp),a in z.groupby(['layer','contrast','discovery']):ax.plot(a.tau,a.gain,color=PAL[CO.index((sp,co))],ls='-' if layer=='supported50' else '--',marker='o',ms=2,label=context(sp,co)+(' S' if layer=='supported50' else ' E'))
        ax.axhline(0,color=GREY,lw=.6);ax.set(xlabel='Effect threshold τ',ylabel='AUC gain' if j==0 else 'Log-loss reduction',xticks=[.01,.02,.03,.05]);legend(ax,handles=[Line2D([],[],c=PAL[i],lw=1,label=context(*co)) for i,co in enumerate(CO)]+[Line2D([],[],c=INK,lw=1,label='Supported'),Line2D([],[],c=INK,lw=1,ls='--',label='Expanded')],loc='lower right',bbox_to_anchor=(1,1.01),ncol=2)
    r.note('All 32 cohort×threshold×PAS-layer settings are retained. S = supported50; E = expanded primary. Curves are repeat-mean scores, with no new q values.')
DRAW.update(D19=D19,D20=D20)
