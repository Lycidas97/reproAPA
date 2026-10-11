"""Read-only deterministic cross-audit of A1; never repeat its simulations."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.stats import false_discovery_control, fisher_exact

HERE = Path(__file__).resolve().parent
OTHER = HERE.parent/'background'
READS = {}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path, **kwargs):
    READS[str(path)] = sha(path)
    return pd.read_csv(path, sep='\t', float_precision='round_trip', **kwargs)


def main():
    for name in ['PROTOCOL.zh.md', 'analyze.py', 'verify.py', 'test_model.py', 'REPORT.zh.md', 'freeze-receipt.json']:
        READS[str(OTHER/name)] = sha(OTHER/name)
    freeze = json.loads((OTHER/'freeze-receipt.json').read_text())
    assert all(sha(path) == fingerprint for path, fingerprint in freeze['input_sha256'].items())
    assert sha(freeze['global_protocol']) == freeze['global_sha256']
    assert sha(freeze['module_protocol']) == freeze['module_sha256']
    assert sha(OTHER/'analyze.py') == freeze['code_sha256']
    inputs = read(OTHER/'results/background-input.tsv.gz')
    tests = read(OTHER/'results/raw-tests.tsv', dtype={'seed':'string'})
    saved_layers = read(OTHER/'results/layer-null-diagnostics.tsv.gz')
    own_rows = []
    assert tests.groupby('family').size().to_dict() == {'A1-conditional-reference-redesign':8,'A1-direction-composition-redesign':4}
    for (source, species), f in inputs.groupby(['pas_source','species']):
        routes = ['GSE118258'] if species == 'human' else ['CIP','WOT']
        p,m,s = [f[[r+'__'+state+'_DWUI' for r in routes]].to_numpy() for state in ['P','M','S']]
        sign = f.direction.to_numpy()[:,None]
        eligible = f.candidate.to_numpy() & ((p-m)*sign >= .02).all(axis=1)
        dual = eligible & ((p-s)*sign >= .02).all(axis=1)
        np.testing.assert_array_equal(eligible, f.pm_eligible)
        np.testing.assert_array_equal(dual, f.dual)
        for model in ['bounded_MS','S_vector']:
            total_mu=total_var=0.;lower=upper=0;variable_events=fixed_zero=fixed_one=0;max_roundtrip=0.
            for signature,q in f.groupby('stratum'):
                a,b,c = [q[[r+'__'+state+'_DWUI' for r in routes]].to_numpy() for state in ['P','M','S']]
                n=len(q)
                if model == 'S_vector':
                    ss=np.broadcast_to(c[None,:,:],(n,n,len(routes)))
                else:
                    z=np.zeros_like(b)
                    for donor in range(n):
                        for route in range(len(routes)):
                            if c[donor,route]>b[donor,route]:
                                z[donor,route]=(c[donor,route]-b[donor,route])/(1-b[donor,route])
                            elif c[donor,route]<b[donor,route]:
                                z[donor,route]=(c[donor,route]-b[donor,route])/b[donor,route]
                    ss=b[:,None,:]+np.maximum(z,0)[None,:,:]*(1-b[:,None,:])+np.minimum(z,0)[None,:,:]*b[:,None,:]
                    rebuilt=np.zeros_like(ss)
                    up=ss>b[:,None,:];down=ss<b[:,None,:]
                    np.divide(ss-b[:,None,:],1-b[:,None,:],out=rebuilt,where=up)
                    np.divide(ss-b[:,None,:],b[:,None,:],out=rebuilt,where=down)
                    deviation=float(np.abs(rebuilt-z[None,:,:]).max());max_roundtrip=max(max_roundtrip,deviation)
                    assert deviation<1e-10
                assert np.isfinite(ss).all() and (ss>=0).all() and (ss<=1).all()
                ps=a[:,None,:]-ss
                np.testing.assert_allclose(ps,(a-b)[:,None,:]+(b[:,None,:]-ss),atol=1e-12)
                h=(q.pm_eligible.to_numpy()[:,None]&((ps*q.direction.to_numpy()[:,None,None])>=.02).all(axis=2)).astype(int)
                np.testing.assert_array_equal(np.diag(h),q.dual)
                assert not h[~q.pm_eligible.to_numpy()].any()
                row_sum=h.sum(axis=1);col_sum=h.sum(axis=0);total=int(h.sum())
                mu=total/n
                cross=total**2-int(row_sum@row_sum)-int(col_sum@col_sum)+int((h*h).sum())
                second=(h*h).sum()/n+cross/(n*(n-1)) if n>1 else mu*mu
                variance=float(second-mu*mu)
                recorded=saved_layers[(saved_layers.pas_source==source)&(saved_layers.species==species)&(saved_layers.model==model)&(saved_layers.stratum==signature)].iloc[0]
                np.testing.assert_allclose([mu,variance],[recorded['mean'],recorded.variance],atol=1e-10)
                total_mu+=mu;total_var+=variance
                i,j=linear_sum_assignment(h);lower+=int(h[i,j].sum())
                i,j=linear_sum_assignment(h,maximize=True);upper+=int(h[i,j].sum())
                counts=row_sum[q.candidate.to_numpy()]
                variable_events+=int(((counts>0)&(counts<n)).sum());fixed_zero+=int((counts==0).sum());fixed_one+=int((counts==n).sum())
            row=tests[(tests.pas_source==source)&(tests.species==species)&(tests.test==model)].iloc[0]
            np.testing.assert_allclose([total_mu,total_var],[row.analytic_null_mean,row.analytic_null_variance],atol=1e-9)
            assert [lower,upper]==[row.attainable_min,row.attainable_max]
            path=OTHER/'results'/f'null-{source}-{species}-{model}.npz';READS[str(path)]=sha(path)
            draws=np.load(path)['draws'];exceed=int((draws>=row.observed).sum())
            assert len(draws)==20000 and exceed==row.exceedances
            np.testing.assert_allclose((1+exceed)/20001,row.p_value,atol=1e-15)
            own_rows.append(dict(pas_source=source,species=species,model=model,observed=int(row.observed),
                analytic_mean=total_mu,analytic_variance=total_var,attainable_min=lower,attainable_max=upper,
                candidate_event_variable=variable_events,candidate_event_fixed_zero=fixed_zero,candidate_event_fixed_one=fixed_one,
                max_residual_roundtrip_error=max_roundtrip,exceedances=exceed,raw_p=row.p_value,status='PASS'))
        row=tests[(tests.pas_source==source)&(tests.species==species)&(tests.test=='direction_composition')].iloc[0]
        q=f[f.candidate];table=[[int((q.dual&q.direction.eq(-1)).sum()),int((q.dual&q.direction.eq(1)).sum())],
                              [int((~q.dual&q.direction.eq(-1)).sum()),int((~q.dual&q.direction.eq(1)).sum())]]
        assert table==[[row.dual_proximal,row.dual_distal],[row.complement_proximal,row.complement_distal]]
        np.testing.assert_allclose(fisher_exact(table,alternative='two-sided').pvalue,row.p_value,atol=1e-12)
    for family,q in tests.groupby('family'):
        tests.loc[q.index,'q_BH_independent']=false_discovery_control(q.p_value,method='bh')
    root_corrected=read(OTHER.parent/'background-root-correction.tsv').set_index('test_id')
    np.testing.assert_allclose(tests.q_BH_independent,tests.test_id.map(root_corrected.q_BH_family),atol=1e-12)
    tests[['test_id','family','p_value','q_BH_independent']].to_csv(HERE/'background-BH-audit.tsv',sep='\t',index=False,float_format='%.17g')
    pd.DataFrame(own_rows).to_csv(HERE/'background-moments-audit.tsv',sep='\t',index=False,float_format='%.17g')
    assert all(sha(path)==fingerprint for path,fingerprint in READS.items())
    assert all(sha(path)==fingerprint for path,fingerprint in freeze['input_sha256'].items())
    result=dict(status='PASS',finished_utc=datetime.now(timezone.utc).isoformat(),new_simulations=0,
                randomization_rows_checked=8,direction_tables_checked=4,original_inputs_SHA_unchanged=True,
                root_family_BH_matches_independent_scipy=True,
                other_module_readonly_SHA_unchanged=True,read_sha256=READS,code_sha256=sha(__file__))
    (HERE/'background-audit-receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(pd.DataFrame(own_rows).to_string(index=False));print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
