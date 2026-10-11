"""Seal the complete H1 delivery without modifying upstream inputs."""
from pathlib import Path
import hashlib,json
import pandas as pd
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];B=ROOT/'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 validation=json.loads((HERE/'validation.json').read_text());assert validation['status']=='PASS'
 tests=(HERE/'tests.log').read_text();assert 'Ran 7 tests' in tests and '\nOK\n' in tests and 'skipped' not in tests.lower()
 sources=[HERE/n for n in ['PROTOCOL.zh.md','extract_remote.py','collect.py','collect_axes.py','analyze.py','verify.py','test_analysis.py','finalize.py','inputs.npz','raw_primary_axis_features.tsv.gz']]
 sources += [ROOT/'notes/2026-10-06-dual-main-plot-data/provenance/polyapipe-pipeline.config.json',ROOT/'notes/2026-09-26-state-comparability-audit/expression_inputs.npz',B/'state_vectors/PAS_axes.primary.tsv.gz',B/'cross_overlap/ortholog_matched_features.tsv.gz']
 for sid in ['GSE103221','WOT_GSE115943','GSE118258','GSE150311','GSE178325','GSE221739']:sources.append(B/'metacell'/sid/'metacell_pas_coverage.tsv.gz')
 own={p.name:sha(p) for p in sorted(HERE.iterdir()) if p.is_file() and p.name!='manifest.json' and '.partial.' not in p.name}
 panels=pd.read_csv(HERE/'panel_eligibility.tsv',sep='\t');a=pd.read_csv(HERE/'context_associations.tsv',sep='\t');lib=pd.read_csv(HERE/'library_associations.tsv',sep='\t')
 assert len(a)==1440 and a.status.eq('available').all();assert len(lib)==2592 and lib.status.eq('descriptive').all()
 manifest=dict(status='complete',analysis='H1 cycle-marker transcription and within-state APA association',created='2026-10-10',setting='v2_primary/supported50_p08',protocol_before_new_matched_results=True,protocol_sha256=sha(HERE/'PROTOCOL.zh.md'),remote_read_only=True,manuscript_changed=False,metacells_are_not_biological_replicates=True,cross_species_inference=False,causal_inference=False,tests_passed=7,tests_skipped=0,validation=validation,metacells_scored=3532,contexts_total=len(panels),contexts_available=int(panels.status.eq('available').sum()),context_association_rows=len(a),library_association_rows=len(lib),input_sha256={str(p.relative_to(ROOT)):sha(p) for p in sources},output_sha256={p.name:sha(p) for p in sorted(HERE.glob('*.tsv*'))},own_file_sha256=own,remote_input_sha256=json.loads((HERE/'remote_input_sha256.json').read_text()))
 (HERE/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps({k:manifest[k] for k in ['status','tests_passed','tests_skipped','metacells_scored','contexts_total','contexts_available','context_association_rows','library_association_rows']},indent=2))
if __name__=='__main__':main()
