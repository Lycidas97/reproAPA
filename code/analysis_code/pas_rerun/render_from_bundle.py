"""Render an explicit table recipe without importing a statistical adapter.

Recipes contain only table/column names, equality filters, labels and geometry.
They cannot execute code. All estimates/intervals/denominators come from the bundle.
"""
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def render(bundle,recipe,out):
    bundle=Path(bundle);out=Path(out);assert not out.exists(),'Use a fresh render directory'
    marker=json.loads((bundle/'bundle-complete.json').read_text())
    assert marker['status']=='statistics_bundle_complete'
    assert sha(bundle/'table-manifest.json')==marker['table_manifest_sha256']
    assert sha(bundle/'plot-contract.json')==marker['plot_contract_sha256']
    manifest=json.loads((bundle/'table-manifest.json').read_text());tables={r['path']:r for r in manifest['tables']}
    contract=json.loads((bundle/'plot-contract.json').read_text());panels={r['panel']:r for r in contract['panels']}
    spec=json.loads(Path(recipe).read_text());out.mkdir(parents=True);used={}
    for item in spec['panels']:
        panel=item['panel'];assert panel in panels
        path=item['table'];assert path in panels[panel]['tables']
        assert sha(bundle/path)==tables[path]['sha256'];used[path]=tables[path]['sha256']
        f=pd.read_csv(bundle/path,sep='\t')
        for key,values in item.get('where',{}).items():f=f[f[key].isin(values)]
        assert not f.empty,f'{panel}: requested result unavailable'
        label=f[item['label']].astype(str).tolist();x=f[item['value']].to_numpy(float)
        assert len(set(label))==len(label),'Recipe must select a unique display row per label'
        fig,ax=plt.subplots(figsize=(item.get('width',5),max(2,len(f)*.28)))
        y=list(range(len(f)))
        if item['kind']=='forest':
            ax.plot(x,y,'o',color='#0072b2',ms=4)
            if 'low' in item:
                lo=f[item['low']].to_numpy(float);hi=f[item['high']].to_numpy(float)
                ax.hlines(y,lo,hi,color='#263d4a',linewidth=.8)
        elif item['kind']=='bar':ax.barh(y,x,color='#0072b2',height=.65)
        else:raise ValueError('Only forest and bar marks are implemented in this table consumer')
        ax.set_yticks(y,label);ax.invert_yaxis();ax.set_xlabel(item['xlabel']);ax.set_title(panel+' · '+item['title'])
        if 'reference' in item:ax.axvline(item['reference'],color='#888888',lw=.6,ls='--')
        ax.spines[['top','right']].set_visible(False);fig.tight_layout()
        fig.savefig(out/(panel+'.pdf'));fig.savefig(out/(panel+'.png'),dpi=180);plt.close(fig)
    (out/'render-receipt.json').write_text(json.dumps(dict(run_id=marker['run_id'],recipe_sha256=sha(recipe),tables=used,
        outputs={p.name:sha(p) for p in out.iterdir() if p.suffix in ['.pdf','.png']},statistics_recomputed=False,
        publication_layout_verified=False,production_replaced=False),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--bundle',required=True);p.add_argument('--recipe',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();render(a.bundle,a.recipe,a.out)
