"""Summarize frozen full GEX for descriptive state review; no model fitting."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path('/root/nfsdata/apadetect/analysis_projects')
P = ROOT / 'human_time_aware_integration_v1'
OUT = P / 'results_state_flow_v1'
PANELS = {
    'mesenchymal_ecm': 'COL1A1 COL1A2 COL3A1 COL6A1 COL6A2 COL6A3 DCN LUM FN1 PRRX1 THY1 PDGFRB',
    'epithelial': 'EPCAM CDH1 KRT8 KRT18 KRT19 CLDN6 CLDN7 DSP OCLN SPINT2',
    'pluripotency_support': 'NANOG LIN28A DPPA4 DNMT3B TDGF1 L1TD1 ESRG SALL4 UTF1 GDF3',
    'naive_support': 'KLF17 TFCP2L1 DPPA3 DNMT3L ALPG NLRP7 KHDC1L TRIML2',
    'primed_context': 'ZIC2 OTX2 FGF5 LEFTY1 LEFTY2 NODAL',
    'endoderm_secretory': 'GATA4 GATA6 SOX17 FOXA1 FOXA2 TTR APOA1 AFP FABP1 HNF4A ALB',
    'keratinocyte': 'KRT5 KRT14 KRT6A KRT17 TP63 IVL SFN LY6D SLURP1 SLURP2',
    'neural': 'SOX1 PAX6 NEUROD1 TUBB3 NEFL SOX11 NES MAP2',
    'trophoblast': 'GATA3 TFAP2C GATA2 KRT7 TEAD4 ELF5 HAND1',
    'ifn': 'ISG15 IFIT1 IFIT2 IFIT3 MX1 MX2 STAT1 OAS1 IFI6 IFI27',
    'cycle': 'MKI67 TOP2A CCNB1 CDC20 PCNA MCM2 TYMS RRM2',
    'exogenous_factor_caution': 'POU5F1 SOX2 KLF4 MYC',
}


def main():
    OUT.mkdir(exist_ok=True)
    registry = pd.read_csv(P/'results_v1/prepared/primary_cells.tsv.gz', sep='\t', low_memory=False)
    all_rows, lib_rows, coverage, checks, provenance = [], [], [], [], []
    for study, expected in registry.groupby('study_id', sort=False):
        version = 'results_qc_v3_2_gse221739_v1' if study == 'GSE221739' else 'results_qc_v3_2_studywide_v1'
        src = ROOT/'reprogramming_qc_v3_2_clustering_release_v1'/version/study.lower()/'scopes/studywide'
        cells = pd.read_csv(src/'prepared/cells.tsv.gz', sep='\t')
        features = pd.read_csv(src/'prepared/features.tsv.gz', sep='\t')
        frozen = pd.read_csv(src/'frozen/cell_clusters.tsv.gz', sep='\t')
        assert cells.cell_id.is_unique and set(cells.cell_id) == set(expected.cell_id)
        ann = expected.set_index('cell_id').loc[cells.cell_id]
        assert frozen.set_index('cell_id').loc[cells.cell_id].cluster_id.tolist() == ann.cluster_id.tolist()
        wanted = sorted({g for genes in PANELS.values() for g in genes.split()})
        indices = [i for i, g in enumerate(features.gene_symbol) if g in wanted]
        symbols = features.gene_symbol.iloc[indices].tolist()
        assert len(symbols) == len(set(symbols)), 'Ambiguous gene symbols require explicit aggregation'
        matrix_file = src/'prepared/log1p_gex.npz'
        x = sparse.load_npz(matrix_file)
        assert x.shape == (len(cells), len(features))
        x = x[:, indices].toarray()
        assert np.isfinite(x).all() and (x >= 0).all()
        for panel, genes in PANELS.items():
            for gene in genes.split():
                coverage.append(dict(study_id=study,panel=panel,gene_symbol=gene,found=gene in symbols))
        for cluster in sorted(ann.cluster_id.unique()):
            mask = ann.cluster_id.eq(cluster).to_numpy()
            sub = x[mask]
            for j, gene in enumerate(symbols):
                all_rows.append(dict(study_id=study,cluster_id=cluster,gene_symbol=gene,n_cells=int(mask.sum()),
                    mean_log1p=float(sub[:,j].mean()),detected_fraction=float((sub[:,j]>0).mean())))
        for (cluster, library), idx in ann.reset_index().groupby(['cluster_id','library_id']).groups.items():
            sub = x[np.asarray(idx)]
            for panel, genes in PANELS.items():
                jj=[symbols.index(g) for g in genes.split() if g in symbols]
                lib_rows.append(dict(study_id=study,cluster_id=cluster,library_id=library,panel=panel,
                    n_cells=len(idx),n_genes=len(jj),mean_log1p=float(sub[:,jj].mean()),
                    mean_detected_fraction=float((sub[:,jj]>0).mean())))
        checks.append(dict(study_id=study,n_cells=len(cells),n_clusters=ann.cluster_id.nunique(),
            exact_cell_set=True,exact_frozen_cluster_labels=True,matrix_shape_valid=True))
        for f in [src/'prepared/cells.tsv.gz',src/'prepared/features.tsv.gz',src/'prepared/effective_pipeline.tsv',matrix_file]:
            h=hashlib.sha256()
            with f.open('rb') as handle:
                for chunk in iter(lambda:handle.read(8*1024*1024),b''):h.update(chunk)
            provenance.append(dict(path=str(f),sha256=h.hexdigest(),size_bytes=f.stat().st_size))
        print(study, len(cells), ann.cluster_id.nunique(), 'PASS', flush=True)
    for name, rows in [('state_gene_expression',all_rows),('state_library_panel_expression',lib_rows),
                       ('signature_coverage',coverage),('input_checks',checks),('expression_provenance',provenance)]:
        pd.DataFrame(rows).to_csv(OUT/(name+'.tsv'),sep='\t',index=False)
    (OUT/'descriptive_panels.json').write_text(json.dumps(PANELS,indent=2)+'\n')


if __name__ == '__main__':
    main()
