"""GO BP / KEGG ORA with frozen evaluable backgrounds and complete audit tables."""
from pathlib import Path
from collections import defaultdict
from functools import lru_cache
import hashlib,json
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

P=Path(__file__).resolve().parents[1];SOURCE=P/'sources/figure2-cde-v1';DATA=P/'data/figure2-cde-v1';OUT=P/'outputs/figure2-analysis-v1'

def save(d,name):d.to_csv(OUT/name,sep='\t',index=False,compression='gzip' if name.endswith('.gz') else None)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def bh(p):
    a=np.asarray(p);order=np.argsort(a,kind='stable');q=np.empty(len(a));q[order]=np.minimum(1,np.minimum.accumulate((a[order]*len(a)/np.arange(1,len(a)+1))[::-1])[::-1]);return q

def ontology():
    source=(SOURCE/'go-basic.obo').read_text();terms={};aliases={}
    for block in source.split('\n[Term]\n')[1:]:
        record=defaultdict(list)
        for line in block.split('\n[')[0].splitlines():
            if ': ' in line:
                k,v=line.split(': ',1);record[k].append(v)
        if not record['id'] or record['is_obsolete']==['true']:continue
        term=record['id'][0];parents=[x.split()[0] for x in record['is_a']]
        parents += [x.split()[1] for x in record['relationship'] if x.startswith('part_of ')]
        terms[term]=dict(name=record['name'][0],namespace=record['namespace'][0],parents=parents)
        for a in record['alt_id']:aliases[a]=term
    @lru_cache(None)
    def ancestors(t):
        if t not in terms:return frozenset()
        return frozenset([t]).union(*(ancestors(p) for p in terms[t]['parents']))
    version=next(l for l in source.splitlines() if l.startswith('data-version:'))
    return terms,aliases,ancestors,version

def biomart(path,names):
    assert path.read_bytes().rstrip().endswith(b'[success]'),path
    d=pd.read_csv(path,sep='\t',header=None,names=names,dtype=str,keep_default_na=False)
    return d[~d.gene_id.eq('[success]')]

def main():
    OUT.mkdir(exist_ok=True);terms,aliases,ancestors,version=ontology();all_results=[];audit=[];mapping=[];background_sets={}
    for sp,org in [('mouse','mmu'),('human','hsa')]:
        bg=pd.read_csv(DATA/f'{sp}_evaluable_background.tsv',sep='\t');bgids=set(bg.gene_id);symbols=dict(zip(bg.gene_id,bg.gene_symbol));assert bg.gene_id.is_unique
        ids=biomart(SOURCE/f'{sp}_gene_ids.tsv',['gene_id','symbol','entrez_id']);ids=ids[ids.gene_id.isin(bgids)]
        entrez=ids[ids.entrez_id.ne('')].groupby('gene_id').entrez_id.agg(set).to_dict();rev=defaultdict(set)
        for gene,values in entrez.items():
            for eid in values:rev[eid].add(gene)
        if (SOURCE/f'{sp}_go.tsv').exists():go=biomart(SOURCE/f'{sp}_go.tsv',['gene_id','term_id','term_name','namespace'])
        else:
            go=biomart(SOURCE/f'{sp}_go_ids.tsv',['gene_id','term_id'])
            go['namespace']=go.term_id.map(lambda t:terms.get(aliases.get(t,t),{}).get('namespace',''))
        go=go[go.gene_id.isin(bgids)&go.namespace.eq('biological_process')&go.term_id.ne('')]
        sets=defaultdict(set);direct=defaultdict(set);unresolved=set()
        for row in go.itertuples():
            t=aliases.get(row.term_id,row.term_id)
            if t not in terms:unresolved.add(row.term_id);continue
            direct[row.gene_id].add(t)
            for parent in ancestors(t):
                if terms[parent]['namespace']=='biological_process' and parent!='GO:0008150':sets[parent].add(row.gene_id)
        gnames={t:terms[t]['name'] for t in sets}
        links=pd.read_csv(SOURCE/f'{org}_pathway_genes.tsv',sep='\t',header=None,names=['entrez','pathway'],dtype=str)
        ks=defaultdict(set)
        for row in links.itertuples():ks[row.pathway.replace('path:','')].update(rev.get(row.entrez.split(':',1)[1],set()))
        kn=pd.read_csv(SOURCE/f'{org}_pathway_names.tsv',sep='\t',header=None,names=['term_id','name'],dtype=str)
        knames=dict(zip(kn.term_id.str.replace('path:',''),kn.name.str.replace(r' - (Mus musculus|Homo sapiens).*','',regex=True)))
        for database,termsets,names in [('GO_BP',sets,gnames),('KEGG',ks,knames)]:
            termsets={t:genes for t,genes in termsets.items() if genes};annotated=set().union(*termsets.values())
            assert annotated<=bgids
            eligible={t:g for t,g in termsets.items() if 10<=len(g)<=500};N=len(annotated)
            background_sets[(sp,database)]=eligible
            for gene in sorted(bgids):mapping.append(dict(species=sp,database=database,gene_id=gene,gene_symbol=symbols[gene],annotated=gene in annotated,ncbi_gene_ids=';'.join(sorted(entrez.get(gene,set()))),ncbi_id_count=len(entrez.get(gene,set())),direct_GO_BP_terms=len(direct.get(gene,set()))))
            for cluster,d in bg[bg.selected].groupby('cluster'):
                full=set(d.gene_id);fg=full&annotated;n=len(fg);rows=[]
                for t,genes in sorted(eligible.items()):
                    hit=fg&genes;k=len(hit);K=len(genes);expected=n*K/N
                    rows.append(dict(species=sp,database=database,cluster=cluster,term_id=t,term_name=names[t],hits=k,foreground_size=n,foreground_unmapped=len(full-fg),background_hits=K,background_size=N,evaluable_background=len(bg),expected_hits=expected,fold_enrichment=k/expected,p_value=float(hypergeom.sf(k-1,N,K,n)),hit_genes=';'.join(sorted(hit)),hit_symbols=';'.join(str(symbols[g]) for g in sorted(hit)),background_genes=';'.join(sorted(genes))))
                result=pd.DataFrame(rows);result['q_value']=bh(result.p_value);result['display_eligible']=result.hits.ge(3)&result.fold_enrichment.gt(1)
                all_results.append(result)
                audit.append(dict(species=sp,database=database,cluster=cluster,evaluable_background=len(bg),annotated_background=N,foreground=len(full),annotated_foreground=n,terms_tested=len(result),q05=int(result.q_value.lt(.05).sum()),display_eligible=int(result.display_eligible.sum()),unresolved_GO_IDs=sorted(unresolved) if database=='GO_BP' else []))
    results=pd.concat(all_results,ignore_index=True)
    # Also retain a stricter correction across all four groups within a species/database.
    results['q_four_clusters']=results.groupby(['species','database']).p_value.transform(lambda x:bh(x))
    save(results,'enrichment_all_tests.tsv.gz');save(pd.DataFrame(audit),'enrichment_background_audit.tsv');save(pd.DataFrame(mapping),'gene_annotation_coverage.tsv')
    # Select three nonredundant terms per database by minimum q then p across groups/species.
    # GO uses shared IDs; KEGG compares species-specific prefixes via its 5-digit pathway ID.
    results['shared_term']=results.term_id.where(results.database.eq('GO_BP'),results.term_id.str[-5:])
    selected=[]
    for database in ['GO_BP','KEGG']:
        candidates=results[results.database.eq(database)&results.display_eligible].sort_values(['q_value','p_value','term_id']).drop_duplicates('shared_term')
        chosen=[]
        for row in candidates.itertuples():
            redundant=False
            for previous in chosen:
                # Reject strong annotation-set redundancy in either species (Jaccard > 0.5).
                for sp,org in [('mouse','mmu'),('human','hsa')]:
                    a=row.shared_term if database=='GO_BP' else org+row.shared_term;b=previous if database=='GO_BP' else org+previous
                    ga=background_sets[(sp,database)].get(a,set());gb=background_sets[(sp,database)].get(b,set())
                    if ga and gb and len(ga&gb)/len(ga|gb)>.5:redundant=True
            if redundant:continue
            chosen.append(row.shared_term);selected.append(dict(database=database,shared_term=row.shared_term,term_name=row.term_name,selection_species=row.species,selection_cluster=row.cluster,selection_q=row.q_value,selection_p=row.p_value,display_order=len(selected)))
            if len(chosen)==3:break
    selection=pd.DataFrame(selected);save(selection,'main_enrichment_terms.tsv')
    plot=results.merge(selection[['database','shared_term','display_order']],on=['database','shared_term'],validate='many_to_one');save(plot,'main_enrichment_values.tsv')
    hashes={str(f):sha(f) for f in SOURCE.iterdir() if f.is_file() and not f.name.endswith('.json')}
    hashes.update({str(f):sha(f) for f in DATA.glob('*_evaluable_background.tsv')})
    (OUT/'enrichment_audit.json').write_text(json.dumps(dict(status='PASS',ontology_version=version,background='database-annotated subset of frozen structurally eligible complete-profile genes',statistical_unit='unique Ensembl gene; any mapped NCBIGene annotation; no repeated gene counting',GO_propagation='is_a and part_of; obsolete/unresolved terms excluded; BioMart GO xrefs, all evidence types',BH_family='all 10–500-background-gene terms per species/database/cluster; q_four_clusters additionally available',selection='three terms/database; hits>=3, fold>1; min q then p; reject Jaccard>0.5 to a previously selected term in either species',source_sha256=hashes,script_sha256=sha(Path(__file__))),indent=2)+'\n')
    print(pd.DataFrame(audit).drop(columns=['unresolved_GO_IDs']).to_string(index=False));print(selection.to_string(index=False))

if __name__=='__main__':main()
