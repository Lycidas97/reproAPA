"""Coverage-only inspection: never calculate APA usage, shifts or predictors."""
from pathlib import Path
import hashlib
import io
import json

import numpy as np
import pandas as pd
from scipy import sparse

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
INPUT = ROOT / 'notes/2026-10-10-marker-reference-sensitivity/inputs.npz'
ROSTER = ROOT / 'notes/2026-10-10-rna-module-enrichment/target_members.tsv'


def frame(z, key):
    return pd.read_json(io.StringIO(str(z[key])), orient='records')


z = np.load(INPUT, allow_pickle=False)
sid = 'GSE115943'
meta = frame(z, sid + '_meta_json')
features = frame(z, sid + '_features_json')
features['column'] = np.arange(len(features))
supported = features[features.coordinate_eligible & features.matched_source_count_50nt.ge(1)]
counts = sparse.csc_matrix((z[sid + '_pas_data'], z[sid + '_pas_indices'],
                           z[sid + '_pas_indptr']), shape=tuple(z[sid + '_pas_shape']))
core = pd.read_csv(ROSTER, sep='\t').query('contrast=="P_M" and foreground_002')
assert len(core) == 27 and core.mouse_gene_id.is_unique
pre = json.loads((HERE / 'coverage_preflight.json').read_text())
libraries = sorted(pre['all24']['gene_counts_per_library'])
assert len(libraries) == 24
rows = []
for gene in core.itertuples(index=False):
    sites = supported[supported.primary_gene_id.eq(gene.mouse_gene_id)]
    assert len(sites) >= 2
    signed = sites.cut_site.to_numpy() * np.where(sites.strand.eq('+'), 1, -1)
    data = counts[:, sites.column.to_numpy()].toarray()
    for library in libraries:
        groups = []
        record = dict(gene_id=gene.mouse_gene_id, gene_symbol=gene.mouse_gene_symbol,
                      library_key=library, n_axis_PAS=len(sites))
        for node in ['P', 'M']:
            ix = np.flatnonzero(meta.library_key.eq(library) & meta.original_selected &
                                meta.original_node.eq(node) & meta.original_purity.ge(.8))
            assert len(ix) >= 3 and meta.n_cells.iloc[ix].sum() >= 30
            raw = data[ix]
            qualified_ix = ix[raw.sum(1) >= 20]
            qualified = data[qualified_ix]
            record.update({f'pre_gene_n_MC_{node}': len(ix),
                           f'pre_gene_n_cells_{node}': int(meta.n_cells.iloc[ix].sum()),
                           f'raw_UMI_{node}': int(raw.sum()),
                           f'raw_detected_PAS_{node}': int((raw.sum(0) > 0).sum()),
                           f'qualified_n_MC_{node}': len(qualified_ix),
                           f'qualified_n_cells_{node}': int(meta.n_cells.iloc[qualified_ix].sum()),
                           f'qualified_UMI_{node}': int(qualified.sum()),
                           f'qualified_detected_PAS_{node}': int((qualified.sum(0) > 0).sum())})
            groups.append((raw, qualified))
        record['raw_50_UMI_each_state'] = min(record['raw_UMI_P'], record['raw_UMI_M']) >= 50
        joint = ((groups[0][1].sum(0) + groups[1][1].sum(0) >= 10) &
                 ((groups[0][1] > 0).sum(0) + (groups[1][1] > 0).sum(0) >= 2))
        rawjoint = ((groups[0][0].sum(0) + groups[1][0].sum(0) >= 10) &
                    ((groups[0][0] > 0).sum(0) + (groups[1][0] > 0).sum(0) >= 2))
        record['qualified_joint_supported_distinct_PAS'] = int(np.unique(signed[joint]).size)
        record['raw_joint_supported_distinct_PAS'] = int(np.unique(signed[rawjoint]).size)
        qualified_enough = min(record['qualified_n_MC_P'], record['qualified_n_MC_M']) >= 3 and min(record['qualified_n_cells_P'], record['qualified_n_cells_M']) >= 30
        record['native_pass'] = qualified_enough and record['qualified_joint_supported_distinct_PAS'] >= 2
        record['native_failure_reason'] = ('PASS' if record['native_pass'] else
            'fewer_than_3_gene_qualified_MC_or_30_cells' if not qualified_enough else
            'fewer_than_2_joint_supported_distinct_PAS')
        record['raw50_with_joint_pass'] = record['raw_50_UMI_each_state'] and record['raw_joint_supported_distinct_PAS'] >= 2
        rows.append(record)

audit = pd.DataFrame(rows)
audit.to_csv(HERE / 'raw_coverage_audit.tsv', sep='\t', index=False)
summary = dict(stage='coverage-only, no APA proportions or differences or expression predictors calculated',
    input_sha256=hashlib.sha256(INPUT.read_bytes()).hexdigest(),
    roster_sha256=hashlib.sha256(ROSTER.read_bytes()).hexdigest(),
    n_libraries=len(libraries), n_core_genes=len(core), n_gene_library_rows=len(audit),
    failure_reasons=audit.native_failure_reason.value_counts().to_dict(),
    count_each_library=audit.groupby('library_key')[['native_pass','raw_50_UMI_each_state','raw50_with_joint_pass']].sum().astype(int).to_dict('index'),
    common_genes_each_gate={gate:audit.groupby('gene_id')[gate].all().loc[lambda x:x].index.tolist()
        for gate in ['native_pass','raw_50_UMI_each_state','raw50_with_joint_pass']},
    native_counts_reproduce_preflight=all(int(audit[audit.library_key.eq(lib)].native_pass.sum()) == n
        for lib,n in pre['all24']['gene_counts_per_library'].items()),
    source_provenance={k:v for k,v in json.loads(str(z['provenance_json'])).items()
        if 'WOT_GSE115943' in k and ('pas_counts' in k or 'metacells.tsv' in k or 'features.tsv' in k)})
(HERE / 'raw_coverage_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
print(json.dumps({k:summary[k] for k in ['n_gene_library_rows','failure_reasons','common_genes_each_gate','native_counts_reproduce_preflight']},indent=2))
