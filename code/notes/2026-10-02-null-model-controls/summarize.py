"""Produce compact, complete summaries without selecting favorable settings."""
from common import *

def main():
    functional=[];supported=[]
    for db in ['GO','KEGG']:
        f=read(HERE/f'{db.lower()}_ladder.tsv.gz');f['comparison']=f.analysis_id.str.split('|').str[2]
        for keys,g in f.groupby(['comparison','model']):
            functional.append(dict(database=db,comparison=keys[0],model=keys[1],positions=len(g),
                nominal=int(g.p.lt(.05).sum()),comparison_supported=int(g.q_comparison.lt(.05).sum()),
                database_supported=int(g.q_database_model.lt(.05).sum()),all_grid_supported=int(g.q_all_new_positions.lt(.05).sum()),
                minimum_p=g.p.min(),minimum_q=g.q_comparison.min()))
        supported.append(f[f.q_comparison.lt(.05)])
    save(pd.DataFrame(functional),'functional_summary.tsv');save(pd.concat(supported),'functional_supported.tsv.gz')
    p=read(HERE/'functional_power.tsv.gz');p['comparison']=p.analysis_id.str.split('|').str[2]
    rows=[]
    for key,g in p.groupby(['database','comparison','tested_model','odds_ratio','threshold']):
        rows.append(dict(zip(['database','comparison','tested_model','odds_ratio','threshold'],key))|dict(positions=len(g),
            minimum=g.rejection_probability.min(),q25=g.rejection_probability.quantile(.25),median=g.rejection_probability.median(),
            q75=g.rejection_probability.quantile(.75),maximum=g.rejection_probability.max()))
    save(pd.DataFrame(rows),'functional_power_summary.tsv')
    retention=read(HERE/'motif_site_retention.tsv.gz');gkeys=['species','contrast','category','gene_id']
    genes=retention.groupby(gkeys).agg(regions=('sites','size'),sites=('sites','sum'),
        weighted_sites20=('weighted_sites20','sum'),weighted_sites50=('weighted_sites50','sum')).reset_index()
    save(genes,'motif_retention_genes.tsv')
    rows=[]
    for key,g in genes.groupby(gkeys[:-1]):
        rows.append(dict(zip(gkeys[:-1],key))|dict(genes=len(g),regions=g.regions.sum(),sites=g.sites.sum(),genes_with_sites=int(g.sites.gt(0).sum()),
            genes_net_increase20=int(g.weighted_sites20.gt(0).sum()),genes_net_decrease20=int(g.weighted_sites20.lt(0).sum()),
            weighted_sites20=g.weighted_sites20.sum(),weighted_sites50=g.weighted_sites50.sum(),
            mean_weighted_sites20=g.weighted_sites20.mean(),median_weighted_sites20=g.weighted_sites20.median()))
    save(pd.DataFrame(rows),'motif_retention_summary.tsv')
    cross=read(HERE/'cross_species_ladder.tsv');nulls=np.load(HERE/'cross_species_draws.npz');rows=[]
    for row in cross.itertuples():
        null=nulls[f'{row.contrast}|{row.model}'];hits=np.arange(min(row.human_events,row.mouse_events)+1)
        for alpha in [.05,.025]:
            passing=hits[probability(null,hits)<=alpha]
            rows.append(dict(contrast=row.contrast,model=row.model,alpha=alpha,observed=row.observed,expected=row.expected,
                smallest_rejecting_overlap=int(passing.min()) if len(passing) else np.nan,
                excess_needed=float(passing.min()-row.expected) if len(passing) else np.nan))
    save(pd.DataFrame(rows),'cross_species_resolution.tsv')

if __name__=='__main__':main()
