import argparse
import importlib.metadata
import resource
import numpy as np
import pandas as pd
from scipy import sparse
from common import *

p=argparse.ArgumentParser()
for field in ['scope','features','batch','output']:p.add_argument('--'+field,required=True)
p.add_argument('--seed',type=int,required=True)
args=p.parse_args()
started=time.time()
import torch
import anndata as ad
import scvi

assert torch.cuda.is_available() and torch.cuda.device_count()==1
torch.set_num_threads(2);torch.set_float32_matmul_precision('high');scvi.settings.seed=args.seed
x=obs();features=np.load(RUN/'shared/feature_sets.npz')[args.features]
counts=sparse.load_npz(RUN/'shared/ortholog_counts.npz')[:,features].tocsr()
data=ad.AnnData(X=counts,obs=x[['species_key','true_study']].copy())
data.obs_names=x.cell_id
data.obs['species_study']=data.obs.species_key+'|'+data.obs.true_study
batch_key='species_key' if args.batch=='species' else 'species_study'
if args.scope=='pilot':
    pilot=RUN/'A03/pilot_balance_v2';pilot_receipt=json.loads((pilot/'complete.json').read_text())
    pilot_path=pilot/'pilot_indices.npy'
    expected=next(item['sha256'] for item in pilot_receipt['outputs'] if item['path']==str(pilot_path))
    assert sha(pilot_path)==expected and pilot_receipt['cells']==60000
    train=np.load(pilot_path)
    training_membership=dict(path=str(pilot_path),sha256=expected,receipt_sha256=sha(pilot/'complete.json'))
else:
    train=np.arange(len(x))
    training_membership=dict(selection='all atlas cells',atlas_sha256=sha(RUN/'shared/atlas.parquet'))
training=data[train].copy()
scvi.model.SCVI.setup_anndata(training,batch_key=batch_key)
model=scvi.model.SCVI(training,n_latent=30,gene_likelihood='nb')
torch.cuda.reset_peak_memory_stats();training_started=time.time()
model.train(max_epochs=100,early_stopping=True,accelerator='gpu',devices=1,batch_size=2048)
training_seconds=time.time()-training_started
epochs_recorded=len(model.history['elbo_train']);assert 1<=epochs_recorded<=100
out=Path(args.output)
model.save(out/'model',overwrite=False)
encoding_started=time.time()
latent=model.get_latent_representation(adata=data,batch_size=2048).astype(np.float32)
encoding_seconds=time.time()-encoding_started
assert latent.shape==(len(x),30) and np.isfinite(latent).all()
np.save(out/'latent.npy',latent)
for key,value in model.history.items():
    assert np.isfinite(value.to_numpy(dtype=float)).all(), ('Non-finite training history',key)
    value.to_csv(out/(key+'.tsv'),sep='\t')
runtime=dict(total_seconds=time.time()-started,training_seconds=training_seconds,encoding_seconds=encoding_seconds,epochs_recorded=epochs_recorded,
    requested_max_epochs=100,early_stopping_enabled=True,float32_matmul_precision=torch.get_float32_matmul_precision(),
    peak_gpu_allocated_bytes=int(torch.cuda.max_memory_allocated()),peak_gpu_reserved_bytes=int(torch.cuda.max_memory_reserved()),
    peak_process_RSS_KiB=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
write_json(out/'runtime.json',runtime)
if args.scope=='pilot':
    # A planning estimate from the observed pilot; full training and early stopping may scale differently.
    full_model_seconds=training_seconds/epochs_recorded*100*len(x)/len(train)
    write_json(out/'full_matrix_budget_estimate.json',dict(basis='Observed pilot seconds per epoch and proportional cell count; not a runtime guarantee.',
        training_cells=len(train),full_training_cells=len(x),observed_epochs=epochs_recorded,
        projected_full_model_training_seconds_at_100_epochs=full_model_seconds,
        projected_12_model_training_hours_at_100_epochs=12*full_model_seconds/3600,
        excludes='Full-data loading, checkpoint writing and encoding; four pilot arms provide separate estimates.'))
write_json(out/'complete.json',dict(scope=args.scope,features=args.features,batch=batch_key,seed=args.seed,
    training_cells=len(train),encoded_cells=len(x),feature_count=len(features),latent_dimensions=30,
    max_epochs=100,batch_size=2048,latent_sha256=sha(out/'latent.npy'),
    versions={n:importlib.metadata.version(n) for n in ['torch','scvi-tools','anndata','numpy']},
    physical_GPU=os.environ['CUDA_VISIBLE_DEVICES'],full_dataset_refit=args.scope=='full',
    training_membership=training_membership,feature_sets_sha256=sha(RUN/'shared/feature_sets.npz'),
    training_script_sha256=sha(__file__),runtime=runtime))
