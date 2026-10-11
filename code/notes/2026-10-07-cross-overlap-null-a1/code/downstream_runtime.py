"""Explicit data access and provenance shared by the PAS rerun adapters."""
import ast
import gzip
import hashlib
import json
import re
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy import sparse

from polyapipe_prepare import digest, save
from polyapipe_effects import metric_weights

STUDIES = ['GSE103221','GSE118258','GSE150311','GSE171127','GSE178325','GSE221739','GSE115943']
SETTINGS = {'supported50_p08':(True,.8,3,20), 'supported50_p09':(True,.9,3,20),
            'supported50_p08_mc5':(True,.8,5,20),'expanded_primary_p08':(False,.8,3,20)}
KEY = ['species','study_id','route','gene_id','gene_symbol','contrast','design']


def bh(values):
    p=np.asarray(values,float)
    assert np.isfinite(p).all() and ((p>=0)&(p<=1)).all()
    order=np.argsort(p);out=np.empty(len(p))
    out[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    return out


def seed(*keys):
    return 20261003+int(hashlib.sha256('|'.join(map(str,keys)).encode()).hexdigest()[:12],16)


def local_state_labels(frame,metacell_ids):
    labels=frame.set_index('metacell_id').local_state_id
    assert labels.index.is_unique
    labels=labels.reindex(metacell_ids)
    assert labels.notna().all()
    return labels.astype(str).str.rsplit('::',n=1).str[-1].to_numpy()


def resolve_gene_symbols(genes,annotation):
    """Replace ID placeholders using the pinned annotation, preserving genuine names."""
    ids=genes.gene_id.astype(str).str.split('.').str[0]
    symbols=genes.gene_symbol.fillna('').astype(str).copy()
    placeholder=symbols.eq('')|symbols.str.split('.').str[0].eq(ids)
    names=ids.map(annotation)
    symbols.loc[placeholder&names.notna()]=names[placeholder&names.notna()]
    return symbols


class Context:
    def __init__(self, config, node, out):
        self.config_path=Path(config).resolve()
        self.config=json.loads(self.config_path.read_text())
        self.node=node;self.out=Path(out).resolve();self.out.mkdir(parents=True,exist_ok=True)
        assert not (self.out/'result.json').exists(), 'Use a fresh adapter output directory'
        self.project=Path(self.config['project_root'])
        self.native=Path(self.config['native_run'])
        self.gex=Path(self.config['frozen_gex'])
        self.analysis=Path(self.config['analysis_root'])
        self.inputs={str(self.config_path):digest(self.config_path)}
        self.run_id=self.config['run_id']

    def register(self,path):
        p=Path(path).resolve();value=digest(p)
        if str(p) in self.inputs:assert self.inputs[str(p)]==value, f'Input changed: {p}'
        self.inputs[str(p)]=value
        return p

    def read(self,path,**kwargs):
        p=self.register(path)
        return pd.read_csv(p,sep='\t',float_precision='round_trip',**kwargs)

    def resource(self,key):
        item=self.config['resources'][key]
        p=self.register(item['path'])
        assert self.inputs[str(p)]==item['sha256'], f'Resource changed: {key}'
        return p

    def table(self,frame,name,columns=None):
        if not isinstance(frame,pd.DataFrame):frame=pd.DataFrame(frame)
        if not len(frame.columns):
            assert columns, f'Empty output requires explicit schema: {name}'
            frame=pd.DataFrame(columns=columns)
        path=self.out/name;path.parent.mkdir(parents=True,exist_ok=True)
        save(frame,path)
        return path

    def functions(self,relative,names,namespace=None):
        """Use reviewed numerical functions without importing legacy drivers."""
        override=self.config.get('source_overrides',{}).get(relative)
        path=self.register(override['path'] if override else self.project/relative)
        if override:assert self.inputs[str(path)]==override['sha256'], f'Reviewed source changed: {relative}'
        source=path.read_text();tree=ast.parse(source)
        selected=[n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef)) and n.name in names]
        assert {n.name for n in selected}==set(names),(relative,names)
        ns=dict(np=np,pd=pd,json=json,hashlib=hashlib,Path=Path,sparse=sparse,
                __file__=str(path),__name__='pas_rerun_numerics',HERE=self.out,
                save=lambda frame,name:self.table(frame,name),bh=bh)
        ns.update(namespace or {})
        exec(compile(ast.Module(body=selected,type_ignores=[]),str(path),'exec'),ns)
        return SimpleNamespace(**{name:ns[name] for name in names})

    def finish(self,checks,**details):
        assert checks and all(v is True for v in checks.values())
        assert all(digest(p)==value for p,value in self.inputs.items()), 'Inputs changed during adapter execution'
        audit=dict(node=self.node,run_id=self.run_id,inputs=self.inputs,checks=checks,**details)
        (self.out/'audit.json').write_text(json.dumps(audit,indent=2,default=str)+'\n')
        files=[p for p in sorted(self.out.rglob('*')) if p.is_file() and p.name not in ['result.json','command.log']]
        result=dict(node=self.node,run_id=self.run_id,input_fingerprint=self.config['input_fingerprint'],
            status='complete',checks=checks,artifacts=[dict(path=str(p.relative_to(self.out)),sha256=digest(p),
                role='audit' if p.name=='audit.json' else 'source_data') for p in files])
        assert len(files)>1, 'A node needs real source data as well as an audit'
        (self.out/'result.json').write_text(json.dumps(result,indent=2)+'\n')


class Study:
    def __init__(self,ctx,study):
        self.ctx=ctx;self.id=study
        folder='WOT_GSE115943' if study=='GSE115943' else study
        self.root=ctx.native/'studies'/folder
        self.f=ctx.read(self.root/'features.tsv.gz',low_memory=False)
        self.meta=ctx.read(self.root/'selected_metacells.tsv.gz').reset_index(drop=True)
        full=ctx.read(self.root/'metacells.tsv.gz')
        assert np.array_equal(full.metacell_id.iloc[self.meta.matrix_row],self.meta.metacell_id)
        self.counts=sparse.load_npz(ctx.register(self.root/'metacell_by_pas_counts.npz')).tocsr()[self.meta.matrix_row].tocsc()
        assert np.issubdtype(self.counts.dtype,np.integer) and (self.counts.data>=0).all()
        self.design=ctx.read(self.root/'sampling_design.tsv')
        self.pairs=ctx.read(self.root/'endpoint_library_design.tsv')
        primary=self.f[self.f.keep_for_apa_primary.astype(str).str.lower().eq('true') & self.f.primary_gene_id.notna()].copy()
        signed=primary.cut_site*np.where(primary.strand.eq('+'),1,-1)
        lo=signed.groupby(primary.primary_gene_id).transform('min');hi=signed.groupby(primary.primary_gene_id).transform('max')
        primary['position']=(signed-lo)/(hi-lo)
        primary['native_span']=hi-lo
        shape=primary.groupby('primary_gene_id').agg(n=('cut_site','nunique'),c=('chrom','nunique'),s=('strand','nunique'))
        self.primary=primary[primary.primary_gene_id.isin(shape.query('n>=2 and c==1 and s==1').index)]
        self.genes={g:f.copy() for g,f in self.primary.groupby('primary_gene_id',sort=False)}
        self.gex_root=ctx.gex/folder;self.gex_counts=None

    def axis(self,gene,setting='supported50_p08'):
        f=self.genes[gene].copy()
        if SETTINGS[setting][0]:f=f[f.matched_source_count_50nt.ge(1)].copy()
        if f.empty:return f
        w=metric_weights(f)
        for i,m in enumerate(['WPP','RWUI','PDUI','PPUI']):f[m]=w[:,i]
        f['rank']=(f.cut_site*np.where(f.strand.eq('+'),1,-1)).rank(method='first').astype(int)
        f['n']=len(f);f['study_id']=self.id
        return f.reset_index(drop=True)

    def gene_counts(self,axis):
        return self.counts[:,axis.feature_index_1based.to_numpy()-1].toarray()

    def load_gex(self):
        if self.gex_counts is None:
            gm=self.ctx.read(self.gex_root/'metacells.tsv.gz')
            self.gex_genes=self.ctx.read(self.gex_root/'genes.tsv.gz')
            ix=pd.Index(gm.metacell_id).get_indexer(self.meta.metacell_id)
            assert (ix>=0).all()
            self.gex_counts=sparse.load_npz(self.ctx.register(self.gex_root/'metacell_by_gene_counts.npz')).tocsr()[ix].tocsc()
            self.gex_total=gm.total_gex.to_numpy()[ix]
            self.gex_index=pd.Index(self.gex_genes.gene_id)

    def expression(self,gene,rows):
        self.load_gex();name=str(gene).split('.')[0]
        if name not in self.gex_index:return np.full(len(rows),np.nan)
        data=self.gex_counts[:,self.gex_index.get_loc(name)].toarray().ravel()
        return data[rows]/self.gex_total[rows]*1e6

    def gex_symbols(self):
        self.load_gex();species=self.meta.species.iloc[0]
        path=self.ctx.resource(species+'_annotation_gtf')
        if not hasattr(self.ctx,'_annotation_gene_names'):self.ctx._annotation_gene_names={}
        if species not in self.ctx._annotation_gene_names:
            names={}
            with gzip.open(path,'rt') if str(path).endswith('.gz') else path.open() as handle:
                for line in handle:
                    if line.startswith('#'):continue
                    fields=line.rstrip().split('\t')
                    if fields[2]!='gene':continue
                    attr=dict(re.findall(r'(\w+) "([^"]*)"',fields[8]))
                    if 'gene_name' not in attr:continue
                    gene=attr['gene_id'].split('.')[0];symbol=attr['gene_name']
                    assert gene not in names or names[gene]==symbol
                    names[gene]=symbol
            assert names;self.ctx._annotation_gene_names[species]=names
        return resolve_gene_symbols(self.gex_genes,self.ctx._annotation_gene_names[species])

    def expression_summary(self,gene,rows):
        self.load_gex();name=str(gene).split('.')[0]
        if name not in self.gex_index:return dict(GEX=np.nan,GEX_CPM=np.nan,GEX_count=np.nan)
        counts=self.gex_counts[:,self.gex_index.get_loc(name)].toarray().ravel()[rows]
        return dict(GEX=float(np.log2(1+counts/self.gex_total[rows]*1e6).mean()),
            GEX_CPM=float(counts.sum()/self.gex_total[rows].sum()*1e6),GEX_count=int(counts.sum()))

    def state(self,data,route,library,node,setting='supported50_p08',depth=None,exclude=()):
        _,purity,nmc,minimum=SETTINGS[setting];minimum=minimum if depth is None else depth
        groups=self.state_groups(purity)
        all_rows=groups.get((route,library,node),np.array([],dtype=int))
        if exclude:all_rows=all_rows[~np.isin(self._state_ids[all_rows],list(exclude))]
        rows=all_rows[data[all_rows].sum(1)>=minimum]
        c=data[rows];usage=(c/c.sum(1)[:,None]).mean(0) if len(c) else np.zeros(data.shape[1])
        n_cells=int(self._state_cells[rows].sum());group_cells=int(self._state_cells[all_rows].sum())
        return dict(rows=rows,counts=c,usage=usage,site_counts=c.sum(0),site_mc=(c>0).sum(0),
            n_mc=len(rows),n_cells=n_cells,group_mc=len(all_rows),group_cells=group_cells,
            endpoint_eligible=len(rows)>=nmc and n_cells>=30,
            profile_eligible=len(rows)>=nmc and group_cells>=30)

    def state_groups(self,purity):
        # Only frozen metadata membership is cached; gene-depth qualification is recomputed.
        if not hasattr(self,'_state_groups'):
            self._state_groups={};self._state_cells=self.meta.n_cells.to_numpy()
            self._state_ids=self.meta.metacell_id.to_numpy()
        if purity not in self._state_groups:
            self._state_groups[purity]={key:part.index.to_numpy() for key,part in
                self.meta[self.meta.purity.ge(purity)].groupby(['route','library_key','node'],sort=False)}
        return self._state_groups[purity]

    def core(self,kind,setting='supported50_p08'):
        return self.ctx.read(self.root/'effects'/f'{setting}_{kind}.tsv.gz',low_memory=False)

    def profile_states(self,data,setting='supported50_p08',depth=20):
        _,purity,nmc,_=SETTINGS[setting]
        states={}
        for key,rows in self.state_groups(purity).items():
            if len(rows)>=nmc and self._state_cells[rows].sum()>=30:
                states[key]=self.state(data,*key,setting,depth)
        supported=set()
        for route in {k[0] for k in states}:
            values=[v for k,v in states.items() if k[0]==route]
            if ((sum(v['site_counts'] for v in values)>=10)&(sum(v['site_mc'] for v in values)>=2)).sum()>=2:
                supported.add(route)
        return {k:v for k,v in states.items() if k[0] in supported and v['profile_eligible']}
