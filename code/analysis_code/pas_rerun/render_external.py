"""Fresh candidate interval traces and external coordinate intersections."""
import json
import numpy as np
import pandas as pd

from render_common import BLUE,ORANGE,INK,GREY,empty,zero
from render_regulatory import inclusion_example


def focus_features(features,gene):
    """Use the manuscript's named motif examples without selecting effect magnitudes."""
    q=features[features.gene.str.upper().eq(gene.upper())]
    if gene.upper()=='EIF4EBP2':selected=q.kind.eq('shared_region')|(q.kind.eq('site')&q.factor.str.upper().eq('PCBP1'))
    elif gene.upper()=='RAD23B':selected=q.kind.eq('shared_region')|(q.kind.eq('site')&q.factor.eq('AAGTGCT'))
    else:selected=q.kind.eq('shared_region')|q.kind.eq('site')
    return q[selected].drop_duplicates(['gene','kind','chrom','strand','start0','end0'])


def longread_panel(r,ax,gene):
    transcripts=r.read('longread/transcript_audit.tsv');pairs=r.read('longread/compatible_isoform_pairs.tsv')
    selected=focus_features(r.read('examples/external_query_features.tsv'),gene)
    q=pairs[pairs.gene.str.upper().eq(gene.upper())&pairs.feature_id.isin(selected.feature_id)]
    ids=sorted(set(q.included_transcript)|set(q.proximal_transcript)) if len(q) else []
    if not ids:empty(ax,'No current interval has a qualified\nsame-chain inclusion/exclusion pair');return
    t=transcripts[transcripts.transcript_id.isin(ids)].sort_values(['terminal_start0_mm10','terminal_end0_mm10','transcript_id'])
    if len(t)>10:raise ValueError(gene+' has >10 qualifying models; external panel must be reflowed')
    status=r.read('longread/feature_retention.tsv');status=status[status.transcript_id.isin(ids)]
    mapping=r.read('longread/coordinate_liftover.tsv');mapping=mapping[mapping.feature_id.isin(selected.feature_id)&mapping.mapping_status.eq('unique_gapfree_reciprocal')]
    strand=t.strand.iloc[0];assert t.strand.nunique()==1
    anchor=t.terminal_start0_mm10.min() if strand=='+' else t.terminal_end0_mm10.max()
    for i,a in enumerate(t.itertuples()):
        lo,hi=((a.terminal_start0_mm10-anchor,a.terminal_end0_mm10-anchor) if strand=='+' else (anchor-a.terminal_end0_mm10,anchor-a.terminal_start0_mm10))
        ax.plot([lo,hi],[i,i],color=GREY,lw=3)
        ax.text(1.02,i,str(a.n_samples_TPM_ge1),transform=ax.get_yaxis_transform(),va='center',fontsize=6.5)
        for f in mapping.itertuples():
            x=(f.mapping_start0+f.mapping_end0)/2;x=x-anchor if strand=='+' else anchor-x
            row=status[status.transcript_id.eq(a.transcript_id)&status.feature_id.eq(f.feature_id)]
            if not len(row):continue
            included=row.status.iloc[0]=='included';ax.plot(x,i,'s' if f.kind=='shared_region' else 'v',ms=3,mfc=BLUE if included else 'white',mec=BLUE,mew=.5)
    samples=len(r.read('longread/sample_inventory.tsv'))
    ax.set_yticks(range(len(t)),t.transcript_id.str.replace(gene+'.','',regex=False),fontsize=6.5);ax.invert_yaxis();ax.set_xlabel('nt from terminal splice boundary');ax.text(.82,1.05,f'TPM≥1: n/{samples}',transform=ax.transAxes,fontsize=6.5,va='bottom',ha='center')


def inclusion_sensitivity(r,ax,gene):
    data=r.read('examples/inclusion_effect_sensitivity.tsv');data=data[data.species.eq('mouse')&data.gene.eq(gene)&data.contrast.eq('SP')]
    data=data[data.kind.eq('site')&data.label.str.upper().eq('PCBP1')] if gene=='EIF4EBP2' else data[data.kind.eq('interval')]
    if data.empty:empty(ax,'No new qualified interval');return
    combos=list(data[['route','start0','end0']].drop_duplicates().itertuples(index=False,name=None))
    for i,(route,a,b) in enumerate(combos):
        q=data[data.route.eq(route)&data.start0.eq(a)&data.end0.eq(b)]
        for j,row in enumerate(q.itertuples()):
            y=i+(j-.5)*.2
            if np.isfinite(row.equal_library):
                if np.isfinite(row.leave_one_min+row.leave_one_max):ax.plot([100*row.leave_one_min,100*row.leave_one_max],[y,y],color=GREY,lw=.7)
                ax.plot(100*row.equal_library,y,'o' if row.cutoff==20 else 's',ms=3,color=BLUE)
                ax.plot(100*row.time_equal,y,'|',ms=5,color=ORANGE)
    zero(ax,False);ax.set_yticks(range(len(combos)),[a+' '+str(b)+'–'+str(c) for a,b,c in combos],fontsize=6.5);ax.invert_yaxis();ax.set_xlabel('P–S Δ inclusion (pp)')


def supplement22(r):
    r.start('Figure S22',150,235)
    for i,gene in enumerate(['HMGN2','EIF4EBP2']):
        r.label(['A','B'][i],2+i*77,2,{'HMGN2':'Hmgn2','EIF4EBP2':'Eif4ebp2'}[gene]+' · library–state inclusion and influence')
        inclusion_example(r,r.axes((16+i*76,17,52,19)),gene,False,legend=False)
        inclusion_sensitivity(r,r.axes((32+i*76,51,36,20)),gene)
    r.label('C',2,84,'Hmgn2 · annotated PASs and newly modeled sites')
    axes=r.read('examples/annotated_example_PAS.tsv');axes=axes[axes.species.eq('mouse')&axes.gene.eq('HMGN2')]
    features=r.read('examples/inclusion_features.tsv');features=features[features.species.eq('mouse')&features.gene.eq('HMGN2')&features.kind.eq('site')]
    ax=r.axes((27,97,42,17))
    if axes.empty:empty(ax,'No new annotated family')
    else:
        strand=axes.strand.iloc[0];origin=axes.cut_interbase.min() if strand=='+' else axes.cut_interbase.max()
        for i,route in enumerate(['CIP','WOT']):
            q=axes[axes.route.eq(route)];x=q.cut_interbase-origin if strand=='+' else origin-q.cut_interbase
            if len(q):ax.plot([x.min(),x.max()],[i,i],color=GREY,lw=2);ax.plot(x,np.full(len(x),i),'|',ms=7,color=INK)
            q=features[features.route.eq(route)];mid=(q.start0+q.end0)/2;mid=mid-origin if strand=='+' else origin-mid
            ax.plot(mid,np.full(len(mid),i),'v',ms=3,color=BLUE)
        ax.set_yticks([0,1],['Guo–CIP','Schiebinger'],fontsize=6.5);ax.invert_yaxis();ax.set_xlabel('nt from proximal PAS†')
    r.label('D',80,84,'Eif4ebp2 · documented-source influence')
    inf=r.read('examples/fixed_source_DWUI_influence.tsv');inf=inf[inf.omission_kind.eq('library')&inf.state.isin(['E','P'])]
    ax=r.axes((104,97,37,17))
    if inf.empty:empty(ax,'No qualified native state profile')
    else:
        for i,(route,state) in enumerate([(a,b) for a in ['CIP','WOT'] for b in ['E','P']]):
            q=inf[inf.route.eq(route)&inf.state.eq(state)];base=q[q.omitted.eq('none')];values=q.loc[~q.omitted.eq('none'),'delta_DWUI'].dropna()
            if len(values):ax.plot([values.min(),values.max()],[i,i],color=GREY,lw=1)
            if len(base):ax.plot(base.delta_DWUI.iloc[0],i,'o',ms=3,color=INK)
        zero(ax,False);ax.set_yticks(range(4),['CIP E–S','CIP P–S','WOT E–S','WOT P–S'],fontsize=6.5);ax.invert_yaxis();ax.set_xlabel('ΔDWUI; leave-library range')
    for i,gene in enumerate(['Hmgn2','Eif4ebp2']):
        r.label(['E','F'][i],2+i*77,129,gene+' · current external structural pairs')
        longread_panel(r,r.axes((28+i*76,145,38,24)),gene)
    r.label('G',2,183,'All ten previously audited loci · new interval support')
    audit=r.read('examples/ten_locus_migration.tsv');genes=audit[['species','gene']].drop_duplicates()
    assert len(genes)==10,'The historical tracking panel must preserve all ten loci'
    for column in range(2):
        ax=r.axes((25+column*75,196,40,26));rows=[];measured=False
        for i,a in enumerate(genes.iloc[column*5:(column+1)*5].itertuples()):
            q=audit[audit.species.eq(a.species)&audit.gene.eq(a.gene)];rows.append(a.gene)
            for j,(lo,hi) in enumerate([('effect_a_min','effect_a_max'),('effect_b_min','effect_b_max')]):
                if q[lo].notna().any():
                    measured=True;ax.plot([q[lo].min()*100,q[hi].max()*100],[i+(j-.5)*.23]*2,color=BLUE if j==0 else ORANGE,lw=.8)
                else:ax.text(.01,i,'×',transform=ax.get_yaxis_transform(),color=GREY,fontsize=6.5)
            ax.text(1.02,i,f'{int(q.matched_region_pairs.sum())} '+('+' if q.exact_partial_recurrence.any() else '?'),transform=ax.get_yaxis_transform(),fontsize=6.5,va='center')
        zero(ax,False);ax.set_yticks(range(len(rows)),rows,fontsize=6.5);ax.invert_yaxis();ax.set_xlabel('Study A/B interval effects (pp)')
        if not measured:ax.set_xticks([]);ax.set_xlabel('No qualified changed regions')
    r.save('Fig/SuppFig.S21.pdf','Current native PASs determine all modeled intervals and site coordinates. A/B show library-state observed ranges and 20/50-UMI/time-equal/leave-library effects, focusing Eif4ebp2 on every current Pcbp1-compatible site. D keeps documented source identity and reweights after each omitted physical library. E/F show structural pairs for the current shared intervals and named motif examples; all queried motifs and transcript pairs remain in S27. G preserves all ten historical gene identities but no former effect ranges or recurrence labels.')


def supplement23(r):
    r.start('Figure S23',150,236)
    r.label('A',2,2,'Rad23b · all qualified common-library inclusion changes')
    inclusion_example(r,r.axes((34,16,105,22)),'RAD23B',True,legend=False)
    r.label('B',2,55,'New interval / seed queries in compatible external transcripts')
    longread_panel(r,r.axes((38,71,97,40)),'Rad23b')
    r.label('C',2,130,'Requeried archived AGO2 cluster overlap')
    query=r.read('ago2/queries.mm9.tsv');mapping=r.read('ago2/query_mapping.tsv');overlap=r.read('ago2/overlap_summary.tsv')
    selected=focus_features(r.read('examples/external_query_features.tsv'),'Rad23b')
    query=query[query.feature_id.isin(selected.feature_id)];mapping=mapping[mapping.feature_id.isin(selected.feature_id)];overlap=overlap[overlap.feature_id.isin(selected.feature_id)]
    track=r.axes((33,145,106,34));ax=r.axes((51,200,88,23))
    if query.empty:empty(ax,'No currently recurrent, reciprocally mappable Rad23b query')
    else:
        q=mapping[mapping.mapping_status.eq('unique_gapfree_reciprocal')].drop_duplicates('feature_id')
        assert q.chrom.nunique()==1 and q.strand.nunique()==1
        strand=q.strand.iloc[0];anchor=q.start0.min() if strand=='+' else q.end0.max()
        def oriented(a,b):return (a-anchor,b-anchor) if strand=='+' else (anchor-b,anchor-a)
        for a in q.itertuples():
            lo,hi=oriented(a.start0,a.end0);row=0 if a.kind=='shared_region' else 1
            track.plot([lo,hi],[row,row],color=BLUE if row==0 else ORANGE,lw=3)
            if row==0:track.text((lo+hi)/2,-.45,f'{int(a.end0-a.start0)} nt',ha='center',fontsize=6.5)
        hits=r.read('ago2/overlapping_clusters.mm39.tsv')
        labels=['WT_merged','KO_merged','WT1A','WT1B','WT2']
        if len(hits):
            mapped=hits[hits.feature_id.isin(selected.feature_id)&hits.mm39_status.eq('unique_gapfree_reciprocal')].drop_duplicates(['label','mm39_chrom','mm39_start0','mm39_end0'])
            for a in mapped.itertuples():
                assert a.mm39_chrom==q.chrom.iloc[0]
                lo,hi=oriented(a.mm39_start0,a.mm39_end0);row=2+labels.index(a.label)
                track.plot([lo,hi],[row,row],color=GREY,lw=2)
        track.set_yticks(range(7),['Interval','miRNA site']+labels,fontsize=6.5);track.invert_yaxis()
        track.set_xlabel(f'GRCm39 transcript-oriented nt from {q.chrom.iloc[0]}:{int(anchor)}',fontsize=7)
        frame=overlap.pivot(index='feature_id',columns='label',values='overlap_bp_sum')
        frame=frame.reindex(columns=['WT_merged','KO_merged','WT1A','WT1B','WT2'])
        if len(frame)>10:raise ValueError('AGO2 has >10 current queries; reflow with all coordinates')
        values=frame.to_numpy(float);im=ax.imshow(np.ma.masked_invalid(values),cmap='Blues',aspect='auto',vmin=0)
        for i in range(len(frame)):
            for j in range(len(frame.columns)):
                v=values[i,j]
                if np.isfinite(v):ax.text(j,i,str(int(v)),ha='center',va='center',fontsize=6.5,color='white' if v>np.nanmax(values)/2 else INK)
        lookup=query.drop_duplicates('feature_id').set_index('feature_id')
        labels=[str(lookup.loc[a,'kind'])+' '+str(lookup.loc[a,'start0'])+'–'+str(lookup.loc[a,'end0']) for a in frame.index]
        ax.set_yticks(range(len(frame)),labels,fontsize=6.5);ax.set_xticks(range(5),frame.columns,fontsize=6.5);ax.set_title('mm9 query coordinates · summed overlap bp',fontsize=7)
    if query.empty:empty(track,'No mappable coordinate track')
    r.save('Fig/SuppFig.S22.pdf','Interval/seed identities come from current recurrent features; the display focuses on the prespecified AAGTGCT seed and its current shared interval. The full motif query inventory remains in S27/S28. A uses common physical libraries. B preserves same-intron-chain compatibility and TPM≥1 sample-column counts. C shows GRCm39 tracks only for uniquely reciprocal cluster mappings, followed by overlap bases at mm9 query coordinates. These are archived cluster bounds, not read coverage; WT1A/B are technical libraries of WT1. No historical interval length, overlap or seed distance is assumed.')
