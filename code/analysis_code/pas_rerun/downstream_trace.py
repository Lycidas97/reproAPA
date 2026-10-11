"""Trace new cross-species members and retain the old roster as migration data."""
import argparse

import numpy as np
import pandas as pd

from downstream_runtime import Context,Study,SETTINGS
from downstream_cross import CODES


def run(ctx):
    shared=ctx.read(ctx.analysis/'cross_overlap/new_shared_members.tsv.gz')
    shared=shared[shared.layer.eq('supported50')&shared.threshold.eq(.02)].copy()
    shared['contrast']=shared.contrast.map({'PS':'P_S','PM':'P_M'})
    old=ctx.read(ctx.resource('legacy_cross_members'))
    key=['human_gene_id','mouse_gene_id','contrast']
    roster=old[key].merge(shared[key],on=key,how='outer',indicator=True,validate='one_to_one')
    roster['old_member']=roster._merge.ne('right_only');roster['new_member']=roster._merge.ne('left_only')
    roster['migration']=roster._merge.map({'both':'retained','left_only':'not_reselected','right_only':'newly_selected'}).astype(str)
    roster=roster.drop(columns='_merge');ctx.table(roster,'member_migration.tsv')
    published=ctx.read(ctx.resource('legacy_cross_support'))
    ctx.table(published,'historical_published_support.tsv')
    global_effects=ctx.read(ctx.analysis/'global/qualified_route_effects.tsv.gz')
    context=ctx.read(ctx.analysis/'global/global_direction_summary.tsv')
    stability=ctx.read(ctx.analysis/'global/direction_stability.tsv')
    source=ctx.read(ctx.analysis/'global/single_source_library_summary.tsv')
    nodes=ctx.read(ctx.analysis/'cross_overlap/all_mapped_nodes.tsv.gz')
    selected=roster[['human_gene_id','contrast','old_member','new_member']]
    ctx.table(nodes.merge(selected,on='human_gene_id',validate='many_to_many'),'selected_observed_state_nodes.tsv.gz')
    lookup={sid:Study(ctx,sid) for sid,_ in CODES.values()};rows=[];effects=[]
    for member in roster.to_dict('records'):
        for code,(sid,route) in CODES.items():
            sp='mouse' if code in ['W','C'] else 'human';gene=member[sp+'_gene_id'];study=lookup[sid]
            for setting in SETTINGS:
                core=study.core('node_positions',setting);f=core[core.route.eq(route)&core.gene_id.eq(gene)]
                measured=f.set_index('node').position_index.to_dict()
                for state in ['S','R','E','P','M']:
                    rows.append(dict(member,species=sp,study_id=sid,route=route,code=code,gene_id=gene,setting=setting,state=state,
                                     APA=measured.get(state,np.nan),status='measured' if state in measured else 'not_measurable'))
                design='true_source' if member['contrast']=='P_S' else 'same_library'
                q=global_effects[global_effects.species.eq(sp)&global_effects.study_id.eq(sid)&global_effects.route.eq(route)&
                    global_effects.gene_id.eq(gene)&global_effects.setting.eq(setting)&global_effects.contrast.eq(member['contrast'])&global_effects.design.eq(design)]
                assert len(q)<=1
                axis=study.axis(gene,setting) if gene in study.genes else pd.DataFrame()
                effects.append(dict(member,species=sp,study_id=sid,route=route,code=code,gene_id=gene,setting=setting,design=design,
                    delta=q.delta.iloc[0] if len(q) else np.nan,new_endpoint_support=bool(len(q) and q.delta.iloc[0]<=-.02),
                    status='measured' if len(q) else 'gene_not_primary_multiPAS' if gene not in study.genes else
                        'fewer_than_two_retained_PAS' if len(axis)<2 else 'endpoint_or_joint_support_not_qualified'))
    columns=key+['old_member','new_member','migration','species','study_id','route','code','gene_id','setting','design','delta','new_endpoint_support','status']
    effects=pd.DataFrame(effects,columns=columns);ctx.table(effects,'member_endpoint_effects.tsv.gz')
    ctx.table(rows,'member_state_trajectories.tsv.gz',columns=key+['species','study_id','route','gene_id','setting','state','APA','status'])
    join=['species','study_id','route','contrast','design','setting']
    trace=effects.merge(context,on=join,suffixes=('_member','_background'),how='left',validate='many_to_many')
    trace['member_proximal_background_distal']=trace.new_endpoint_support&trace['mean'].gt(0)
    trace['background_measurable']=trace.n.fillna(0).gt(0)
    ctx.table(trace,'member_background_trace.tsv.gz')
    ctx.table(effects.merge(stability,on=[x for x in join if x!='setting'],how='left',validate='many_to_one'),'member_setting_source_stability.tsv.gz')
    ctx.table(source.merge(roster[['contrast']].drop_duplicates(),on='contrast',validate='many_to_one'),'source_background_scenarios.tsv')
    ctx.finish(dict(new_members_determine_current_set=True,historical_roster_retained_separately=True,
        all_four_settings_traced=True,missing_effects_are_NA=True,observed_S_not_conflated_with_documented_source=True,
        background_and_source_sensitivities_joined=True),old_members=len(old),new_members=len(shared))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);a=p.parse_args();run(Context(a.config,'cross_trace',a.out))
