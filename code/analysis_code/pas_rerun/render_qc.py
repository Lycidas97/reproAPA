"""PAS-support panels, separating the retained historical counting audit."""
import numpy as np
import pandas as pd

from downstream_runtime import STUDIES
from render_common import BLUE,ORANGE,INK,GREY,LABELS,empty,zero


def supplement4(r):
    r.start('Figure S4',150,210)
    cat=r.read('pas_qc/catalog_summary.tsv');primary=cat[cat.layer.eq('primary')]
    r.label('A',2,2,'Native PAS annotation support within 50 nt')
    ax=r.axes((15,13,125,27));x=np.arange(len(primary))
    ax.bar(x-.17,primary.supported50_sites/primary.sites,.3,color=BLUE,label='PAS fraction')
    ax.bar(x+.17,primary.supported50_UMI/primary.UMI,.3,color=ORANGE,label='UMI fraction')
    ax.set_xticks(x,primary.study_id.str.replace('GSE',''),rotation=20);ax.set_ylim(0,1);ax.set_ylabel('Supported / primary');ax.legend(frameon=False,ncol=2,loc='lower right',bbox_to_anchor=(1,1.05))
    seq=[];curves=[];reads=[];matched=[]
    for sid in STUDIES:
        folder='WOT_GSE115943' if sid=='GSE115943' else sid
        sp=cat.loc[cat.study_id.eq(sid),'species'].iloc[0]
        seq.append(r.read(f'pas_qc/{folder}/sequence_summary.tsv').assign(study_id=sid,species=sp))
        curves.append(r.read(f'pas_qc/{folder}/sequence_curves.tsv.gz').assign(study_id=sid,species=sp))
        reads.append(r.read(f'pas_qc/{folder}/read_profile_curves.tsv.gz').assign(study_id=sid,species=sp))
        ledger=r.read(f'pas_qc/{folder}/sequence_background_matches.tsv.gz');metrics=r.read(f'pas_qc/{folder}/sequence_site_metrics.tsv.gz')
        ids=set(ledger.loc[ledger.status.eq('matched'),'atlas_pas_id'])
        for name,categories in [('Native candidate',['supported50','unsupported50']),('Matched control',['matched_genomic_background'])]:
            f=metrics[metrics.atlas_pas_id.isin(ids)&metrics.category.isin(categories)]
            assert len(f)==len(ids)
            matched.append(dict(study_id=sid,species=sp,category=name,n=len(f),positive=int(f.upstream_signal_40_10.astype(str).str.lower().eq('true').sum())))
    seq=pd.concat(seq);curves=pd.concat(curves);reads=pd.concat(reads);matched=pd.DataFrame(matched)
    r.label('B',2,51,'Canonical upstream signal: matched windows')
    ax=r.axes((15,63,52,28));matched=matched.groupby(['species','category'])[['n','positive']].sum().reset_index()
    for i,sp in enumerate(['human','mouse']):
        f=matched[matched.species.eq(sp)].set_index('category')
        for j,name in enumerate(['Native candidate','Matched control']):
            row=f.loc[name];ax.bar(i+(j-.5)*.3,row.positive/row.n if row.n else np.nan,.27,color=[BLUE,GREY][j],label=name if i==0 else None)
            ax.text(i+(j-.5)*.3,.98,f'n={row.n}',rotation=90,va='top',ha='center',fontsize=6.5,transform=ax.get_xaxis_transform())
    ax.set_xticks([0,1],['Human','Mouse']);ax.set_ylim(0,1);ax.set_ylabel('Signal-positive fraction');ax.legend(frameon=False,fontsize=6.5,loc='upper center',bbox_to_anchor=(.5,-.22),ncol=2)
    r.label('C',78,51,'Canonical signal position')
    ax=r.axes((91,70,51,21));q=curves[curves.alignment.eq('representative')&curves.feature.eq('canonical_signal_start')]
    colors={'supported50':BLUE,'unsupported50':ORANGE,'native_misprime_flagged_in_median_library':'#9865A3','annotation_excluded_in_median_library':GREY,'matched_genomic_background':INK}
    short={'supported50':'Supported','unsupported50':'Unsupported','native_misprime_flagged_in_median_library':'Native flag','annotation_excluded_in_median_library':'Excluded','matched_genomic_background':'Matched control'}
    for category,color in colors.items():
        f=q[q.category.eq(category)].groupby('position_nt')[['count','denominator']].sum()
        ax.plot(f.index,f['count']/f.denominator,color=color,label=short[category])
    ax.set(xlabel='Signal start relative to PAS (nt)',ylabel='Window frequency',xlim=(-80,50))
    ax.legend(frameon=False,fontsize=6.5,loc='lower center',bbox_to_anchor=(.5,1.02),ncol=3,handlelength=.7,columnspacing=.6)
    r.label('D',2,106,'Read-end QC: selected local windows')
    for i,sp in enumerate(['human','mouse']):
        ax=r.axes((15+i*66,118,52,24))
        for category,color in [('supported50',BLUE),('unsupported50',ORANGE)]:
            f=reads[reads.species.eq(sp)&reads.category.eq(category)].groupby('position_nt')[['normalized_site_sum','selected_sites']].sum()
            ax.plot(f.index,f.normalized_site_sum/f.selected_sites,color=color,label=short[category])
        ax.set(title=sp.title(),xlabel='Read end relative to PAS (nt)',ylabel='Mean site fraction',xlim=(-400,400))
        if i==0:ax.legend(frameon=False)
    r.label('E',2,158,'Historical six-rule audit · original pipeline')
    changes=pd.read_csv(r.resource('historical_counting_steps'),sep='\t')
    fields=[key for key in changes if key.endswith('p90_abs_delta_pa')];assert len(fields)==5
    ax=r.axes((14,171,53,20))
    for _,g in changes.iterrows():
        ax.plot(range(5),g[fields].to_numpy(float),'o-',ms=2,color=GREY)
    ax.set_xticks(range(5),['Dedup','Intervals','Atlas','Primary','All'],rotation=25,ha='right',fontsize=6.5)
    ax.set_ylabel('P90 |ΔDWUI difference|');ax.set_xlabel('Successive historical counting steps')
    r.label('F',78,158,'Native depth sensitivity\nObserved-state profiles')
    nodes=r.read('profiles/all_route_node_positions.tsv.gz');nodes=nodes[nodes.setting.eq('supported50_p08')]
    wide=nodes.pivot(index=['species','route','depth','gene_id'],columns='node',values='position_index').reindex(columns=['S','P','M'])
    rows=[]
    for ref in ['S','M']:
        effect=(wide.P-wide[ref]).dropna().rename('effect').reset_index()
        base=effect[effect.depth.eq(20)].set_index(['species','route','gene_id']).effect
        for depth in [40,80]:
            alt=effect[effect.depth.eq(depth)].set_index(['species','route','gene_id']).effect
            for sp in ['human','mouse']:
                a=base.loc[base.index.get_level_values('species')==sp];b=alt.reindex(a.index);common=b.notna()
                diff=(b[common]-a[common]).abs()
                rows.append(dict(species=sp,reference=ref,depth=depth,base_n=len(a),retained_n=int(common.sum()),
                    fraction=common.mean() if len(a) else np.nan,median=diff.median(),q90=diff.quantile(.9)))
    data=pd.DataFrame(rows);data.to_csv(r.out/'S4-native-profile-depth-summary.tsv',sep='\t',index=False)
    ax=r.axes((91,174,22,17))
    for sp,color in [('human',BLUE),('mouse',ORANGE)]:
        f=data[data.species.eq(sp)];ax.scatter(np.arange(len(f)),f.fraction,color=color,s=10,label=sp)
    ax.set_xticks(range(4),['S40','S80','M40','M80'],rotation=60,ha='right',fontsize=6.5);ax.set_ylim(0,1);ax.set_ylabel('Retained / 20 UMIs')
    ax.legend(frameon=False,fontsize=6.5,ncol=2,loc='upper center',bbox_to_anchor=(.5,1.42),borderaxespad=0,handletextpad=.2,columnspacing=.5)
    ax=r.axes((125,174,18,17))
    for sp,color in [('human',BLUE),('mouse',ORANGE)]:
        f=data[data.species.eq(sp)];ax.scatter(np.arange(len(f)),f['median'],s=9,edgecolors=color,facecolors='white')
        ax.scatter(np.arange(len(f)),f.q90,s=9,color=color)
    ax.set_xticks(range(4),['S40','S80','M40','M80'],rotation=60,ha='right',fontsize=6.5);ax.set_ylabel('|Effect difference|')
    r.save('Fig/SuppFig.S4.pdf','A–D and F use native polyApipe data. E is explicitly historical original-pipeline evidence and is excluded from native tests. Read-end QC does not re-quantify the UMI matrix. F uses observed-state profile contrasts; full effect differences accompany its source table.')
