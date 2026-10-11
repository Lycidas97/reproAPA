"""Replay frozen recognition-site rules on the existing union models locally."""
from pathlib import Path
from collections import defaultdict, Counter
import ast
import csv
import gzip
import hashlib
import importlib.util
import json
import re
import subprocess
import sys
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
UNION = HERE.parent/'2026-10-01-mpro-biology'
F5 = ROOT/'analysis_code/release_2026-09-24/sources/F5'
RES = HERE/'resources'
ROUTES = {'mouse':['CIP','WOT'], 'human':['GSE118258','GSE150311']}
MR = RES/'mouse_time_aware_integration_pilot_v1'
HR = RES/'human_time_aware_integration_v1'
FIMO = Path('/private/tmp/meme-5.5.2/src/fimo')

def module(name, path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

structure=module('union_structure',UNION/'legacy_structure_consensus.py')
core=module('union_core',F5/'full-cohort-v2/sources/execution-scripts/core.py')
migration=module('union_migration',F5/'annotation-pas-global-v1/migrate.py')
paired=module('union_paired',ROOT/'notes/2026-09-29-core-sensitivity/figure5_paired_pm.py')
ns={'defaultdict':defaultdict}
source=(MR/'regulatory_catalog_v1/scripts/build_catalog.py').read_text()
for node in ast.parse(source).body:
    if isinstance(node,ast.FunctionDef) and node.name in ['rc','canonical_patterns','seed_hits']:
        exec(compile(ast.get_source_segment(source,node),'<frozen seed functions>','exec'),ns)

def bundles():
    raw=json.loads((UNION/'remote_resource_paths.json').read_text())
    return {sp:next(v for k,v in raw.items() if '/'+sp+'_time_aware' in k) for sp in ROUTES}

def prepare(sp):
    native={r:json.loads((UNION/f'native/{sp}/{r}/models.json').read_text()) for r in ROUTES[sp]}
    wanted={m['region_id'] for genes in native.values() for g in genes.values()
            if g['eligible_contrasts'] for m in g['family']['models']}
    regions=pd.read_csv(UNION/f'{sp}_catalog_utrs.tsv.gz',sep='\t').set_index('region_id')
    regions=regions.loc[sorted(wanted)].reset_index().to_dict('records')
    seqs={}; rid=None; pieces=[]
    with gzip.open(RES/f'{sp}_sequences.fa.gz','rt') as f:
        for line in f:
            if line.startswith('>'):
                if rid in wanted: seqs[rid]=''.join(pieces).upper()
                rid=line[1:].split()[0]; pieces=[]
            elif rid in wanted: pieces.append(line.strip())
        if rid in wanted: seqs[rid]=''.join(pieces).upper()
    assert set(seqs)==wanted
    out=HERE/sp; out.mkdir(exist_ok=True)
    with (out/'query.fa').open('w') as f:
        for r in regions:
            rid=r['region_id']; assert len(seqs[rid])==r['end0']-r['start0']
            f.write('>'+rid+'\n'+seqs[rid]+'\n')
    print(sp,'regions',len(wanted),'bases',sum(map(len,seqs.values())),flush=True)
    return native,{r['region_id']:r for r in regions},seqs

def scan(sp):
    native,regions,seqs=prepare(sp); out=HERE/sp; bundle=bundles()[sp]
    bank=(MR/'regulatory_query_v2/resources/transferred_scannable.meme' if sp=='mouse'
          else HR/'regulatory_query_v1/resources/human_direct_scannable.meme')
    bg=MR/'regulatory_catalog_v1/background.txt' if sp=='mouse' else RES/'human_background.txt'
    assert subprocess.check_output([str(FIMO),'--version'],text=True).strip()=='5.5.2'
    argv=[str(FIMO),'--text','--norc','--thresh','0.0001','--bfile',str(bg),str(bank),str(out/'query.fa')]
    if not (out/'fimo.tsv').exists():
        with (out/'fimo.partial.tsv').open('w') as f,(out/'fimo.stderr.log').open('w') as log:
            subprocess.run(argv,stdout=f,stderr=log,check=True)
        (out/'fimo.partial.tsv').rename(out/'fimo.tsv')
    ann=pd.read_csv(UNION/'mouse_mirna_annotations.tsv' if sp=='mouse' else RES/'regulatory_resource_revision_v1/resources/human_mirna_annotations.tsv',sep='\t').fillna('')
    mirnames={seed:sorted({n for v in g.mirbase_exact_aliases for n in v.split(';') if n}) for seed,g in ann.groupby('seed_dna')}
    prov=(pd.read_csv(UNION/'mouse_reference_motif_provenance.tsv.gz',sep='\t') if sp=='mouse'
          else pd.DataFrame(json.loads((RES/'human_motif_provenance.json').read_text())))
    direct_names={mid:sorted(set(g.RBP_Name)) for mid,g in prov[prov.scan_status.eq('scan_ready')].groupby('Motif_ID')}
    transfer_names=defaultdict(set)
    for r in bundle.get('transfer_provenance',[]): transfer_names[r['Motif_ID']].add(r['mouse_symbol'])
    cached={rid:[] for rid in regions}; unmapped=Counter()
    def add(rid,start,end,category,feature,kind):
        r=regions[rid]
        gs,ge=(r['start0']+start,r['start0']+end) if r['strand']=='+' else (r['end0']-end,r['end0']-start)
        names=mirnames.get(feature,[]) if category=='miRNA' else direct_names.get(feature,[]) if category=='RBP' else sorted(transfer_names[feature])
        if not names:
            unmapped[category+':'+feature]+=1
            return  # Frozen pipeline excludes motif/seed links without names.
        cached[rid].append(dict(chrom=r['chrom'],strand=r['strand'],genomic_start0=gs,genomic_end0=ge,
            category=category,feature=feature,names=names,site_type=kind))
    with (out/'fimo.tsv').open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if not r.get('motif_id') or r['motif_id'].startswith('#'): continue
            rid=r['sequence_name']; start,end=int(r['start'])-1,int(r['stop'])
            assert r['strand']=='+' and seqs[rid][start:end]==r['matched_sequence'].upper()
            if sp=='mouse': assert r['motif_id'] in bundle['active_motif_ids']
            add(rid,start,end,'RBP_transfer' if sp=='mouse' else 'RBP',r['motif_id'],'PWM')
    if sp=='mouse':
        direct=pd.read_csv(UNION/'mouse_reference_rbp_hits.tsv.gz',sep='\t')
        for r in direct[direct.region_id.isin(regions)].itertuples():
            assert seqs[r.region_id][r.start0:r.end0]==r.matched_sequence.upper()
            add(r.region_id,r.start0,r.end0,'RBP',r.motif_id,'PWM')
    patterns=ns['canonical_patterns'](bundle['mirna_seeds'])
    matcher=re.compile('(?=('+'|'.join(sorted(patterns,key=lambda x:(-len(x),x)))+'))')
    for rid,seq in seqs.items():
        for start,length,seed,kind in ns['seed_hits'](seq,patterns,matcher): add(rid,start,start+length,'miRNA',seed,kind)
    print(sp,'raw hits',sum(map(len,cached.values())),flush=True)
    coverage=[]; audit=Counter(unmapped_names=sum(unmapped.values())); nrows=0
    dest=HERE/f'full-cohort-v2/data/{sp}';dest.mkdir(parents=True,exist_ok=True)
    with gzip.open(dest/'site_effects.tsv.gz','wt') as f:
        writer=None
        for route,genes in native.items():
            annotated=json.loads((UNION/f'annotation/{sp}/{route}/models.json').read_text())
            for gid,g in genes.items():
                model=annotated[gid]
                if not model['eligible_contrasts']: continue
                hits={}; spans=defaultdict(set)
                for m in g['family']['models']:
                    for h in cached[m['region_id']]:
                        result=structure.classify_site(h,g['family'],g['pas'])
                        if not result['site_status'].startswith('model_consensus'):
                            audit['outside_consensus']+=1;continue
                        key=(h['category'],h['genomic_start0'],h['genomic_end0'])
                        factor=h['feature'] if h['category']=='miRNA' else ';'.join(h['names'])
                        obj=hits.setdefault(key,dict(factors={},types=set()))
                        obj['factors'][factor]=h['names'];obj['types'].add(h['site_type']);spans[key[0],factor].add(key[1:])
                ambiguous={k:core.ambiguous_overlaps(v) for k,v in spans.items()}
                cuts=np.array([p['cut_site'] for p in model['pas']]); fam=model['family']
                stats=migration.state_stats(model)
                pm={t:paired.paired_state_stats(model,t) for t in [20,50]}
                nsites=Counter()
                for (category,start,end),h in hits.items():
                    factors={k:v for k,v in h['factors'].items() if (start,end) not in ambiguous[category,k]}
                    if not factors: audit['all_factor_links_overlap_ambiguous']+=1;continue
                    nsites[category]+=1
                    mask=migration.mask_for(start,end,fam['strand'],cuts)
                    base=dict(species=sp,route=route,gene_id=gid,gene_symbol=g['gene_symbol'],family_id=fam['family_id'],category=category,chrom=fam['chrom'],strand=fam['strand'],start0=start,end0=end,
                        names=';'.join(sorted({n for v in factors.values() for n in v})),feature_groups=';'.join(sorted(factors)),site_types=';'.join(sorted(h['types'])),mask_json=json.dumps(np.column_stack([mask,mask]).tolist(),separators=(',',':')),exact=True,variable=bool(mask.min()!=mask.max()))
                    for contrast in sorted(set(model['eligible_contrasts'])):
                        a,b=contrast;pa,pb=model['profiles'][a],model['profiles'][b]
                        effect=float((np.array(pb['usage'])-np.array(pa['usage']))@mask)
                        row=dict(base,contrast=contrast,pooled_lo=effect,pooled_hi=effect,reference_UMI=pa['total'],target_UMI=pb['total'])
                        for t in [20,50]:
                            if contrast=='MP':
                                effect=float(pm[t]['delta_usage']@mask);na=nb=len(pm[t]['libraries'])
                            else:
                                sa,sb=stats[t,a],stats[t,b];effect=float((sb[0]-sa[0])@mask);na,nb=sa[1],sb[1]
                            row.update({f'eq{t}_lo':effect,f'eq{t}_hi':effect,f'eq{t}_n_reference':na,f'eq{t}_n_target':nb})
                        if writer is None: writer=csv.DictWriter(f,fieldnames=list(row),delimiter='\t');writer.writeheader()
                        writer.writerow(row);nrows+=1
                coverage.append(dict(species=sp,route=route,gene_id=gid,gene_symbol=g['gene_symbol'],**nsites))
            print(sp,route,'site rows',nrows,flush=True)
    pd.DataFrame(coverage).to_csv(dest/'site_coverage.tsv',sep='\t',index=False)
    (out/'scan_manifest.json').write_text(json.dumps(dict(status='PASS',fimo_version='5.5.2',argv=argv,regions=len(regions),site_effect_rows=nrows,audit=dict(audit),sha256={str(p.relative_to(HERE)):digest(p) for p in [bank,bg,out/'query.fa',out/'fimo.tsv',dest/'site_effects.tsv.gz']}),indent=2)+'\n')

if __name__=='__main__': scan(sys.argv[1])
