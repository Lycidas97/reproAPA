"""Freeze the A03 training pilot at the hierarchy specified in the original plan."""
import numpy as np
import pandas as pd
from common import *


def allocate(capacity,budget):
    capacity=np.asarray(capacity,dtype=int);allocation=np.zeros(len(capacity),dtype=int)
    remaining=min(int(budget),int(capacity.sum()))
    while remaining:
        available=np.flatnonzero(allocation<capacity)
        share,remainder=divmod(remaining,len(available))
        addition=np.minimum(capacity[available]-allocation[available],share+(np.arange(len(available))<remainder))
        assert addition.sum()>0
        allocation[available]+=addition;remaining-=int(addition.sum())
    assert allocation.sum()==min(budget,capacity.sum()) and (allocation<=capacity).all()
    uncapped=allocation[allocation<capacity]
    assert len(uncapped)<2 or np.ptp(uncapped)<=1
    return allocation


def run():
    out=RUN/'A03/pilot_balance_v2';out.mkdir(exist_ok=True)
    registry_path=RUN/'process_support_v1/library_registry.tsv.gz'
    legacy_path=RUN/'shared/pilot_indices.npy'
    sources=[Path(__file__),RUN/'shared/atlas.parquet',registry_path,legacy_path]
    signatures=[dict(path=str(p),sha256=sha(p)) for p in sources]
    receipt_path=out/'complete.json'
    if receipt_path.exists():
        rec=json.loads(receipt_path.read_text());assert rec['inputs']==signatures
        for item in rec['outputs']:assert sha(item['path'])==item['sha256']
        print('Existing balanced A03 pilot verified.',flush=True);return
    assert not list((RUN/'A03/models').rglob('complete.json')), 'Freeze before any A03 model result.'
    x=obs();registry=pd.read_csv(registry_path,sep='\t').fillna('')
    assert registry.sample_key.is_unique and len(registry)==161
    identities=x[['sample_key','species_key','true_study']].drop_duplicates()
    joined=identities.merge(registry,on='sample_key',validate='one_to_one')
    assert len(joined)==161
    for col in ['species_key','true_study']:assert joined[col+'_x'].equals(joined[col+'_y'])
    frame=x[['atlas_index','cell_id','species_key','true_study','sample_key']].copy()
    for col in ['facet_id','time_label','condition']:frame[col]=frame.sample_key.map(registry.set_index('sample_key')[col])
    assert frame[['facet_id','time_label','condition']].notna().all().all() and frame.time_label.ne('').all()
    rng=np.random.default_rng(SEED);selected=[]
    for species,species_frame in frame.groupby('species_key',sort=True):
        studies=list(species_frame.groupby('true_study',sort=True))
        study_budget=allocate([len(d) for key,d in studies],30000)
        for (study,study_frame),budget in zip(studies,study_budget):
            strata=list(study_frame.groupby(['facet_id','time_label','condition'],sort=True,dropna=False))
            stratum_budget=allocate([len(d) for key,d in strata],int(budget))
            for (stratum,stratum_frame),count in zip(strata,stratum_budget):
                libraries=list(stratum_frame.groupby('sample_key',sort=True))
                library_budget=allocate([len(d) for key,d in libraries],int(count))
                for (library,data),take in zip(libraries,library_budget):
                    selected.extend(rng.choice(data.atlas_index.to_numpy(),int(take),replace=False))
    selected=np.sort(np.asarray(selected,dtype=np.int64));assert len(selected)==60000 and len(np.unique(selected))==60000
    membership=frame.iloc[selected].copy();assert np.array_equal(membership.atlas_index,selected)
    assert membership.groupby('species_key').size().eq(30000).all() and membership.sample_key.nunique()==161
    groups=['species_key','true_study','facet_id','time_label','condition','sample_key']
    counts=frame.groupby(groups,sort=True).size().rename('available_cells').to_frame()
    counts['training_cells']=membership.groupby(groups,sort=True).size().reindex(counts.index,fill_value=0)
    assert (counts.training_cells>0).all()
    np.save(out/'pilot_indices.npy',selected)
    membership.to_csv(out/'training_cells.tsv.gz',sep='\t',index=False)
    counts.reset_index().to_csv(out/'library_allocation.tsv',sep='\t',index=False)
    assert sha(legacy_path)==signatures[-1]['sha256']
    write_json(receipt_path,dict(cells=60000,libraries=161,seed=SEED,
        allocation='Equal species, then capped equal study, facet/time/condition, and library allocation; lexicographic remainder assignment.',
        time_source='Registered time_label resolves actual day or stage; original condition identities are retained.',
        source_states_used_for_selection=False,legacy_pilot_unchanged=True,
        legacy_overlap_cells=int(np.intersect1d(selected,np.load(legacy_path)).size),inputs=signatures,
        outputs=[dict(path=str(out/name),sha256=sha(out/name)) for name in ['pilot_indices.npy','training_cells.tsv.gz','library_allocation.tsv']]))
    print(json.dumps(json.loads(receipt_path.read_text()),indent=2),flush=True)


if __name__=='__main__':run()
