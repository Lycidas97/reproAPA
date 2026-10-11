"""Versioned bindings from the established supplement design to one verified bundle.

Old reference strings are design identifiers, never filesystem fallbacks.  Only
frozen GEX inputs may remain outside the PAS bundle. Display-only reshapes are
explicit here; tests, candidate selection and uncertainty come from producers.
"""
from pathlib import Path
import json
import pandas as pd
from result_records import ResultBundle, contained
from pipeline import digest

ALIASES = {
    'prepared:E04-LOO-gene-assignments.tsv.gz':'supplement_statistics/LOO-gene-assignments.tsv.gz',
    'prepared:E04-center-display.tsv':'supplement_statistics/LOO-center-display.tsv',
    'prepared:E10-continuous-primary-gene-pairs.tsv.gz':'supplement_statistics/continuous-gene-pairs.tsv.gz',
    'prepared:E17-endpoint-14-term-positions.tsv':'supplement_statistics/endpoint-GO-KEGG-display.tsv',
    'axis:well_six_unit_downsample_summary.tsv':'supplement_statistics/well_six_unit_downsample_summary.tsv',
    'well_inference/gene_metric_tests.tsv.gz':'supplement_statistics/well-metric-tests-audited.tsv.gz',
    'reference_comparison/PAS_support_comparison.tsv':'pas_qc/catalog_summary.tsv',
}
LEGACY_NAMES = {
    'E03b-tandem-comparison.tsv':'panel_statistics/tandem-comparison.tsv',
    'E03c-state-reference-summary.tsv':'geometry/state-reference-summary.tsv',
    'E09-E12-transcript-models.tsv':'geometry/transcript-models.tsv.gz',
    'E09-E12-exon-CDS.tsv':'geometry/exon-CDS.tsv.gz',
    'E14d-axis-decomposition-128.tsv':'geometry/axis-decomposition.tsv',
    'E14d-native-sites.tsv.gz':'geometry/native-sites.tsv.gz',
    'small-position-large-redistribution-full.tsv':'geometry/low-shift-cohort.tsv',
    'Figure3C-counts.tsv':'supplement_statistics/dual-reference-availability.tsv',
    'Figure3C-human-members.tsv':'supplement_statistics/dual-reference-all-measurable.tsv.gz',
    'Figure3C-mouse-members.tsv':'supplement_statistics/dual-reference-all-measurable.tsv.gz',
    'Figure4D-source-bank.tsv':'same_pas/v2_primary_correlation_summary.tsv',
    'Figure4D-coarse-null.tsv':'same_pas/v2_primary_conditional_null_tests.tsv',
    'Figure4D-rank-capacity-null.tsv':'same_pas/capacity_sensitivity_tests.tsv',
    'Figure5C-five-genes.tsv':'geometry/low-shift-cohort.tsv',
    'Figure5D-relationship-counts.tsv':'supplement_statistics/partial-relationship-counts.tsv',
}

class SupplementBundle:
    def __init__(self, records, renderer):
        self.reader=ResultBundle(Path(records).parent,renderer.config)
        self.renderer=renderer
        self.bindings=[]

    def path(self, ref):
        r=self.renderer
        if ref.startswith('root:notes/2026-09-24-integration-top-labels/'):
            path=(Path(__file__).resolve().parents[2]/ref[5:]).resolve()
            assert str(path) in r.pins, 'Unregistered frozen GEX input: '+ref
            assert digest(path)==r.pins[str(path)], 'Frozen GEX input changed: '+ref
            source='frozen_GEX'
            name=ref
        else:
            assert 'v3' not in ref, 'v3 is excluded from this delivery: '+ref
            name=ALIASES.get(ref,ref)
            if ref.startswith('base:'):
                name=LEGACY_NAMES[Path(ref).name]
            for prefix,target in [('axis:','axis_extensions/'),('hyp:A_naive/','naive_state/'),
                                  ('hyp:B_shape/','shape_ortholog/'),('hyp:C_reference/','same_pas/'),
                                  ('hyp:D_direction/','direction_difference/'),('hyp:E_mirna/','external_mirna/')]:
                if name.startswith(prefix):name=target+name[len(prefix):]
            key='tables/'+name
            assert key in self.reader.manifest, 'Unbound supplement input: '+ref+' -> '+key
            path=contained(self.reader.root,key)
            expected=self.reader.manifest[key]['sha256']
            if str(path) not in r.inputs:assert digest(path)==expected, 'Changed table: '+key
            source=self.reader.identity['source_key']
            self.reader.used[key]=expected
        if str(path) not in r.inputs:r.inputs[str(path)]=digest(path)
        h=r.inputs[str(path)]
        binding=dict(figure=r.fid,logical_panel=r.active,design_reference=ref,table=name,
                     path=str(path),sha256=h,source_key=source)
        if binding not in self.bindings:self.bindings.append(binding)
        r.sources.add(ref)
        r.evidence[r.active]['sources'][ref]=binding
        return path

    def read(self, ref, **kwargs):
        path=self.path(ref);r=self.renderer;key=(str(path),str(kwargs))
        if key not in r.cache:
            r.cache[key]=pd.read_csv(path,sep=kwargs.pop('sep','\t'),float_precision='round_trip',**kwargs)
        x=r.cache[key].copy()
        if ref=='reference_comparison/PAS_support_comparison.tsv':
            x=pd.concat([x[['study_id','species','layer','sites',f'supported{n}_sites']].rename(
                columns={'sites':'denominator',f'supported{n}_sites':'supported_v2'}).assign(window_nt=n)
                for n in [25,50,100]],ignore_index=True)
        elif Path(ref).name=='Figure5C-five-genes.tsv':x=x[x.strict_region_pairs.gt(0)].copy()
        elif Path(ref).name=='E03b-tandem-comparison.tsv':
            x=x[x.contrast.eq('P_S')&x.design.eq('true_source')&x.setting.eq('supported50_p08')&x.background.eq('available')].copy()
        elif Path(ref).name=='Figure3C-counts.tsv':
            x=x.groupby('species').agg(candidate_genes=('gene_id','size'),measurable_genes=('matched_measurable','sum'),dual_genes=('dual_qualified_gene','sum')).reset_index()
        elif Path(ref).name.startswith('Figure3C-') and 'members' in ref:
            x=x[x.species.eq('human' if 'human' in ref else 'mouse')].copy()
        return x

    def write_bindings(self):
        out=self.renderer.out
        current=pd.DataFrame(self.bindings)
        previous=out/'supplement-source-bindings.tsv'
        if previous.exists():
            old=pd.read_csv(previous,sep='\t')
            assert set(old.source_key)<={self.reader.identity['source_key'],'frozen_GEX'}, 'Mixed-source output directory'
            current=pd.concat([old[~old.figure.isin(current.figure)],current],ignore_index=True)
        current.to_csv(previous,sep='\t',index=False)
        provenance=self.reader.provenance()
        provenance['tables']={'tables/'+v.table:v.sha256 for v in current.itertuples() if v.source_key!='frozen_GEX'}
        (out/'bundle-provenance.json').write_text(json.dumps(provenance,indent=2))
