"""Resolve shared terminal-exon families without identifying full transcripts."""
import hashlib
import json
from collections import defaultdict


def stable_id(prefix, values):
    return prefix+hashlib.sha256(json.dumps(values,sort_keys=True,separators=(',', ':')).encode()).hexdigest()[:20]


def contains_cut(model, pas):
    return (model['gene_id']==pas['gene_id'] and model['chrom']==pas['chrom'] and model['strand']==pas['strand']
            and model['start0'] < pas['member_cut_site_min'] <= pas['member_cut_site_max'] <= model['end0'])


def resolve_families(pas_rows, models):
    if len({p['atlas_pas_id'] for p in pas_rows})!=len(pas_rows):
        raise ValueError('Duplicate PAS would duplicate the denominator')
    groups=defaultdict(list)
    for m in models:
        if m['strand'] not in ['+','-'] or not m['start0']<m['end0']:
            raise ValueError('Invalid terminal model coordinates')
        groups[(m['gene_id'],m['chrom'],m['strand'],m['exon_anchor'])].append(m)
    families=[]
    for key, members in sorted(groups.items()):
        ids=sorted(m['structure_id'] for m in members)
        if len(set(ids))!=len(ids):
            raise ValueError('Duplicate model')
        families.append(dict(family_id=stable_id('FAM_', [key,ids]),gene_id=key[0],chrom=key[1],strand=key[2],
            exon_anchor=key[3],models=members,common_start0=max(m['start0'] for m in members),
            common_end0=min(m['end0'] for m in members)))
    inventory=[]
    for p in pas_rows:
        if not p['member_cut_site_min']<=p['cut_site']<=p['member_cut_site_max']:
            raise ValueError('Representative cut outside cluster')
        hits=[f for f in families if any(contains_cut(m,p) for m in f['models'])]
        r=dict(p,family_id='',family_status='no_annotated_family',supporting_models='',sibling_boundary_models='')
        if len(hits)>1:
            r['family_status']='ambiguous_between_terminal_exon_families'
        elif len(hits)==1:
            family=hits[0]
            if any(m['catalog_status']!='matched' for m in family['models']):
                r['family_status']='incomplete_model_resources'
            else:
                r.update(family_id=family['family_id'],family_status='assigned',
                    supporting_models=';'.join(sorted(m['structure_id'] for m in family['models'] if contains_cut(m,p))),
                    sibling_boundary_models=';'.join(sorted(m['structure_id'] for m in family['models'] if not contains_cut(m,p))))
        inventory.append(r)
    for family in families:
        rows=sorted([p for p in inventory if p['family_id']==family['family_id']],key=lambda p:p['member_cut_site_min'])
        if any(a['member_cut_site_max']>=b['member_cut_site_min'] for a,b in zip(rows,rows[1:])):
            for p in rows:
                p.update(family_id='',family_status='overlapping_pas_clusters')
        family['pas_ids']=[p['atlas_pas_id'] for p in rows if p['family_id']]
    return families,inventory


def genomic_retention(site, pas):
    if site['chrom']!=pas['chrom'] or site['strand']!=pas['strand']:
        raise ValueError('Site/PAS coordinate system mismatch')
    if not 0<=site['genomic_start0']<site['genomic_end0']:
        raise ValueError('Invalid site interval')
    lo,hi=pas['member_cut_site_min'],pas['member_cut_site_max']
    if site['strand']=='+':
        return (int(site['genomic_end0']<=lo),int(site['genomic_end0']<=hi))
    if site['strand']=='-':
        return (int(site['genomic_start0']>=hi-1),int(site['genomic_start0']>=lo-1))
    raise ValueError('Unknown strand')


def classify_site(site, family, pas_rows):
    support=[m['structure_id'] for m in family['models']
        if m['chrom']==site['chrom'] and m['strand']==site['strand']
        and m['start0']<=site['genomic_start0']<site['genomic_end0']<=m['end0']]
    if len(support)!=len(family['models']):
        return dict(site_status='model_dependent_utr' if support else 'outside_all_annotated_utrs',
            supporting_site_models=';'.join(sorted(support)),retention_by_pas={})
    if len(pas_rows)<2 or set(p['atlas_pas_id'] for p in pas_rows)!=set(family['pas_ids']):
        raise ValueError('Complete multi-PAS family axis required')
    mask={p['atlas_pas_id']:genomic_retention(site,p) for p in pas_rows}
    return dict(site_status='model_consensus_cut_uncertain' if (0,1) in mask.values() else 'model_consensus_exact',
        supporting_site_models=';'.join(sorted(support)),retention_by_pas=mask)


def site_key(site, family_id):
    # Do not count the same genomic prediction once per transcript/UTR query.
    return stable_id('SITE_', [family_id]+[site[k] for k in
        ['kind','feature_id','site_type','chrom','strand','genomic_start0','genomic_end0']])
