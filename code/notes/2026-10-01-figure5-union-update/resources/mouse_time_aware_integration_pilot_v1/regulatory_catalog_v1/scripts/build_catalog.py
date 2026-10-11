"""Effect-independent mouse multi-PAS sequence catalog; no matrix modification."""
import argparse
import csv
import hashlib
import io
import json
import re
import sqlite3
import subprocess
import sys
import time
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

ROOT = Path('/root/nfsdata/apadetect/analysis_projects/mouse_time_aware_integration_pilot_v1')
REF = Path('/root/nfsdata/apadetect/data/genome/mouse/refdata-gex-GRCm39-2024-A')
MIRNA = Path('/root/nfsdata/multica/worktrees/apadetect/DOC-248-reproapa-mirna-resources/analyses/reprogramming-mirna-expression-evidence-v2/tables/mouse_reprogramming_candidate_rpm10.tsv')
FIMO = '/root/miniconda3/envs/apasim/bin/fimo'
STUDIES = ['gse103221_pluripotency_reprogramming_core',
           'gse115943_pluripotency_reprogramming_core',
           'gse171127_pluripotency_vs_trophoblast_bifurcation']


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')


def rc(seq):
    return seq.upper().translate(str.maketrans('ACGTUN', 'TGCAAN'))[::-1]


def genomic_interval(start0, end0, strand, offset0, length):
    assert 0 <= offset0 and offset0 + length <= end0 - start0
    return ((start0 + offset0, start0 + offset0 + length) if strand == '+'
            else (end0 - offset0 - length, end0 - offset0))


def window(cut1, strand, chrom_length):
    # cut_site denotes the 1-based genomic cleavage base in the released atlas.
    c = cut1 - 1
    lo, hi = (c - 200, c + 100) if strand == '+' else (c - 99, c + 201)
    lo, hi = max(0, lo), min(chrom_length, hi)
    return lo, hi, c - lo if strand == '+' else hi - 1 - c


def canonical_patterns(seeds):
    patterns = defaultdict(list)
    for seed in sorted(set(seeds)):
        seed = seed.upper().replace('U', 'T')
        assert len(seed) == 7 and set(seed) <= set('ACGT')
        patterns[rc(seed) + 'A'].append((seed, '8mer'))
        patterns[rc(seed)].append((seed, '7mer-m8'))
        patterns[rc(seed[:6]) + 'A'].append((seed, '7mer-A1'))
    return patterns


def seed_hits(seq, patterns, matcher):
    # Lookahead retains overlapping sites. Suppress nested 7mers for each seed only.
    by_seed_anchor = {}
    rank = {'8mer': 3, '7mer-m8': 2, '7mer-A1': 1}
    for m in matcher.finditer(seq):
        pos = m.start()
        for size in (8, 7):
            word = seq[pos:pos + size]
            if len(word) != size:
                continue
            for seed, kind in patterns.get(word, []):
                anchor = pos - 1 if kind == '7mer-A1' else pos
                key = (seed, anchor)
                candidate = (pos, size, seed, kind)
                if key not in by_seed_anchor or rank[kind] > rank[by_seed_anchor[key][3]]:
                    by_seed_anchor[key] = candidate
    return sorted(by_seed_anchor.values())


def candidate_subset(df):
    good = df.keep_for_sensitivity & df.primary_gene_id.notna() & df.primary_gene_id.ne('')
    selected = df[good].copy()
    sizes = selected.groupby('primary_gene_id').atlas_pas_id.nunique()
    return selected[selected.primary_gene_id.isin(sizes[sizes >= 2].index)].copy()


def prepare(out):
    dbpath = out / 'catalog.sqlite'
    if dbpath.exists():
        raise RuntimeError('Catalog exists; use a new version, not an implicit overwrite')
    con = sqlite3.connect(dbpath)
    sources = [REF / 'fasta/genome.fa', REF / 'fasta/genome.fa.fai', REF / 'genes/genes.gtf', MIRNA,
               ROOT / 'results_tandem_apa_v1/all_reference_terminal_utr_candidates.tsv.gz',
               out / 'resources/CisBP-RNA2_Mus_musculus.zip', Path(__file__)]
    fasta = pysam.FastaFile(str(REF / 'fasta/genome.fa'))
    sizes = dict(zip(fasta.references, fasta.lengths))
    windows = {}
    genes = set()
    totals = []
    for study in STUDIES:
        source = ROOT / 'results_metacell_pas_v7_v1/studies' / study / 'features.tsv.gz'
        sources.append(source)
        df = pd.read_csv(source, sep='\t', low_memory=False)
        assert df.atlas_pas_id.is_unique
        assert set(df.genome_build) == {'GRCm39_2024_A'}
        sel = candidate_subset(df)
        assert (sel.peak_start == sel.cut_site - 1).all() and (sel.peak_end == sel.cut_site).all()
        assert set(sel.strand) <= {'+', '-'}
        sel.insert(0, 'study_id', study)
        ids = []
        for r in sel.itertuples():
            assert 1 <= r.cut_site <= sizes[r.chrom]
            lo, hi, anchor = window(r.cut_site, r.strand, sizes[r.chrom])
            rid = f'W_{r.chrom}_{r.cut_site}_{r.strand}'
            windows[rid] = (rid, 'pas_window', r.chrom, lo, hi, r.strand, anchor, None, 0)
            ids.append(rid)
        sel['window_id'] = ids
        sel.to_sql('pas', con, if_exists='append', index=False)
        gene = sel.groupby('primary_gene_id', sort=True).agg(
            gene_symbol=('primary_gene_symbol', 'first'), biotype=('primary_gene_biotype', 'first'),
            n_pas=('atlas_pas_id', 'nunique'), n_primary_pas=('keep_for_apa_primary', 'sum'),
            n_conflicted=('transcript_model_conflict', 'sum'), n_ambiguous=('multi_gene_ambiguous', 'sum'))
        gene = gene.reset_index().rename(columns={'primary_gene_id': 'gene_id'})
        gene.insert(0, 'study_id', study)
        gene.to_sql('genes', con, if_exists='append', index=False)
        gene.to_csv(out / f'{study}.candidate_genes.tsv.gz', sep='\t', index=False)
        genes.update(gene.gene_id)
        totals.append(dict(study=study, genes=len(gene), pas=len(sel), primary_pas=int(sel.keep_for_apa_primary.sum())))
        print(json.dumps(totals[-1]), flush=True)
    con.execute('CREATE UNIQUE INDEX pas_key ON pas(study_id,atlas_pas_id)')
    con.execute('CREATE INDEX pas_gene ON pas(primary_gene_id,chrom,strand)')
    con.execute('CREATE INDEX pas_window ON pas(window_id)')
    refs = pd.read_csv(sources[4], sep='\t')
    refs = refs[refs.gene_id.isin(genes)]
    keys = ['gene_id', 'chrom', 'strand', 'utr_start_1based', 'utr_end_1based']
    refs = refs.groupby(keys, as_index=False).agg(transcript_ids=('transcript_id', lambda v: ';'.join(sorted(set(v)))))
    regions = list(windows.values())
    utrs = []
    for j, r in enumerate(refs.itertuples()):
        start, end = int(r.utr_start_1based) - 1, int(r.utr_end_1based)
        lo, hi = (start, min(sizes[r.chrom], end + 2000)) if r.strand == '+' else (max(0, start - 2000), end)
        rid = f'U{j:06d}'
        regions.append((rid, 'terminal_utr_plus2kb', r.chrom, lo, hi, r.strand, None, r.gene_id, end - start))
        utrs.append((rid, r.gene_id, r.transcript_ids, start, end, end-start,
                     'annotation_compatible_terminal_utr;extension_transcription_and_neighbor_conflicts_not_validated'))
    cols = ['region_id', 'region_type', 'chrom', 'start0', 'end0', 'strand', 'cut_offset0', 'gene_id', 'annotated_length']
    pd.DataFrame(regions, columns=cols).to_sql('regions', con, index=False)
    pd.DataFrame(utrs, columns=['region_id','gene_id','transcript_ids','annotated_start0','annotated_end0','annotated_length','evidence_scope']).to_sql('utrs', con, index=False)
    con.executescript('''
    CREATE UNIQUE INDEX region_key ON regions(region_id);
    CREATE INDEX region_gene ON regions(gene_id,chrom,strand);
    CREATE INDEX utr_gene ON utrs(gene_id);
    CREATE TABLE pas_utr AS SELECT p.study_id,p.atlas_pas_id,r.region_id,
      CASE WHEN p.strand='+' THEN p.cut_site-r.start0 ELSE r.end0-p.cut_site+1 END AS included_length,
      CASE WHEN p.strand='+' THEN p.member_cut_site_min-r.start0 ELSE r.end0-p.member_cut_site_max+1 END AS min_included_length,
      CASE WHEN p.cut_site BETWEEN u.annotated_start0+1 AND u.annotated_end0 THEN 'annotated_utr' ELSE 'tes_extension_candidate' END AS cut_class
    FROM pas p JOIN regions r ON p.primary_gene_id=r.gene_id AND p.chrom=r.chrom AND p.strand=r.strand
    JOIN utrs u ON u.region_id=r.region_id
    WHERE p.member_cut_site_min>=r.start0+1 AND p.member_cut_site_max<=r.end0;
    CREATE INDEX pas_utr_pas ON pas_utr(study_id,atlas_pas_id);
    CREATE INDEX pas_utr_region ON pas_utr(region_id,min_included_length);
    ''')
    counts = Counter()
    with open(out / 'sequences.fa', 'w') as handle:
        for rid, kind, chrom, lo, hi, strand, anchor, gene, annlen in regions:
            seq = fasta.fetch(chrom, lo, hi).upper()
            if strand == '-':
                seq = rc(seq)
            assert len(seq) == hi-lo
            handle.write(f'>{rid}\n{seq}\n')
            counts.update(seq)
    with open(out / 'background.txt', 'w') as f:
        f.write('\n'.join(f'{b} {counts[b]/sum(counts[x] for x in "ACGT"):.12f}' for b in 'ACGT')+'\n')
    mirna = pd.read_csv(MIRNA, sep='\t')
    assert mirna.analysis_sequence_group_id.is_unique and (mirna.max_rpm > 10).all()
    mirna['seed_dna'] = mirna.seed_2_8.str.replace('U','T')
    mirna.to_sql('mirna_expression', con, index=False)
    motif_registry(out, con)
    con.commit()
    assert con.execute('PRAGMA integrity_check').fetchone()[0] == 'ok'
    con.close()
    manifest = dict(created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        selection='Per-study >=2 assigned sensitivity-inclusive PAS; no APA effect, state, or expression selection',
        studies=totals, unique_genes=len(genes), unique_pas_windows=len(windows), utr_regions=len(utrs),
        sequence_bases=sum(counts.values()), reference='GRCm39_2024_A', human_analysis_started=False,
        rbp_scope='CisBP mouse direct status D only; incomplete RBP universe',
        coordinate_contract='PAS cut_site 1-based base; intervals 0-based half-open; transcript strand only; windows [-200,+100) relative to cleavage base',
        sources=[dict(path=str(f),bytes=f.stat().st_size,sha256=sha(f)) for f in sources])
    dump(out/'prepare.json',manifest)
    print(json.dumps({k:v for k,v in manifest.items() if k!='sources'}), flush=True)


def motif_registry(out, con):
    with zipfile.ZipFile(out / 'resources/CisBP-RNA2_Mus_musculus.zip') as z:
        assert z.testzip() is None
        rows = pd.read_csv(io.BytesIO(z.read('RBP_Information.txt')), sep='\t')
        rows = rows[(rows.RBP_Status == 'D') & (rows.Motif_ID != '.')]
        def status(r):
            name='pwms_all_motifs/'+r.Motif_ID+'.txt'
            if r.Motif_Type=='JPLE':return 'inferred_motif_not_direct_experiment'
            if name not in z.namelist():return 'matrix_missing'
            if not z.read(name).strip():return 'matrix_empty'
            return 'scan_ready'
        rows['scan_status']=[status(r) for r in rows.itertuples()]
        con.execute('DROP TABLE IF EXISTS motif_provenance')
        rows.to_sql('motif_provenance', con, index=False)
        rows.to_csv(out/'mouse_direct_motif_provenance.tsv',sep='\t',index=False)
        with open(out/'mouse_direct.meme','w') as f:
            f.write('MEME version 4\n\nALPHABET= ACGT\n\nstrands: +\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n')
            for mid in sorted(rows.loc[rows.scan_status.eq('scan_ready'),'Motif_ID'].unique()):
                pwm=pd.read_csv(io.BytesIO(z.read('pwms_all_motifs/'+mid+'.txt')),sep='\t')[['A','C','G','U']].to_numpy(float)
                assert len(pwm)>0 and np.isfinite(pwm).all() and (pwm>=0).all()
                assert np.allclose(pwm.sum(1),1,atol=1e-5)
                pwm=pwm/pwm.sum(1)[:,None]
                f.write(f'\nMOTIF {mid}\nletter-probability matrix: alength= 4 w= {len(pwm)} nsites= 20 E= 0\n')
                np.savetxt(f,pwm,fmt='%.10f')


def resume_prepare(out):
    # Resume only the failed resource conversion; candidate tables and FASTA stay intact.
    con=sqlite3.connect(out/'catalog.sqlite')
    for table in ['genes','pas','regions','utrs','pas_utr','mirna_expression']:
        assert con.execute('SELECT count(*) FROM '+table).fetchone()[0]>0
    assert con.execute('SELECT count(*) FROM pas').fetchone()[0]==719229
    motif_registry(out,con)
    fa=pysam.FastaFile(str(out/'sequences.fa'))
    assert len(fa.references)==con.execute('SELECT count(*) FROM regions').fetchone()[0]
    sources=[REF/'fasta/genome.fa',REF/'fasta/genome.fa.fai',REF/'genes/genes.gtf',MIRNA,
        ROOT/'results_tandem_apa_v1/all_reference_terminal_utr_candidates.tsv.gz',
        out/'resources/CisBP-RNA2_Mus_musculus.zip',Path(__file__)]
    sources.extend(ROOT/'results_metacell_pas_v7_v1/studies'/s/'features.tsv.gz' for s in STUDIES)
    manifest=dict(created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
        selection='Per-study >=2 assigned sensitivity-inclusive PAS; no APA effect, state, or expression selection',
        studies=pd.read_sql_query('SELECT study_id,count(*) AS genes,sum(n_pas) AS pas,sum(n_primary_pas) AS primary_pas FROM genes GROUP BY study_id',con).to_dict('records'),
        unique_genes=con.execute('SELECT count(DISTINCT gene_id) FROM genes').fetchone()[0],
        unique_pas_windows=con.execute("SELECT count(*) FROM regions WHERE region_type='pas_window'").fetchone()[0],
        utr_regions=con.execute('SELECT count(*) FROM utrs').fetchone()[0],sequence_bases=sum(fa.lengths),
        reference='GRCm39_2024_A',human_analysis_started=False,
        rbp_scope='6 mouse RBPs / 7 available experimental PWMs; 5 JPLE records excluded, Nova2 empty matrix blocked',
        recovery='Candidate tables and sequence extraction completed in attempt01; only motif conversion repaired',
        coordinate_contract='PAS 1-based cleavage base; BED 0-based half-open; transcript strand only; window [-200,+100)',
        sources=[dict(path=str(f),bytes=f.stat().st_size,sha256=sha(f)) for f in sources])
    con.commit();con.close();dump(out/'prepare.json',manifest)
    print(json.dumps({k:v for k,v in manifest.items() if k!='sources'}),flush=True)


def scan_mirna(out):
    con=sqlite3.connect(out/'catalog.sqlite')
    con.execute('CREATE TABLE mirna_hits(region_id TEXT,seed_dna TEXT,site_type TEXT,start0 INTEGER,end0 INTEGER,genomic_start0 INTEGER,genomic_end0 INTEGER,region_class TEXT)')
    seeds=[r[0] for r in con.execute('SELECT DISTINCT seed_dna FROM mirna_expression')]
    patterns=canonical_patterns(seeds)
    matcher=re.compile('(?=('+'|'.join(sorted(patterns, key=lambda x:(-len(x),x)))+'))')
    fasta=pysam.FastaFile(str(out/'sequences.fa'))
    records=con.execute("SELECT region_id,start0,end0,strand,annotated_length FROM regions WHERE region_type='terminal_utr_plus2kb'").fetchall()
    total=0
    for i,(rid,lo,hi,strand,annlen) in enumerate(records):
        seq=fasta.fetch(rid)
        rows=[]
        for pos,n,seed,kind in seed_hits(seq,patterns,matcher):
            gs,ge=genomic_interval(lo,hi,strand,pos,n)
            cls='annotated_utr' if pos+n<=annlen else 'tes_extension_candidate' if pos>=annlen else 'annotation_boundary_spanning'
            rows.append((rid,seed,kind,pos,pos+n,gs,ge,cls))
        con.executemany('INSERT INTO mirna_hits VALUES(?,?,?,?,?,?,?,?)',rows)
        total+=len(rows)
        if i%1000==0:
            con.commit();print(json.dumps(dict(stage='mirna',regions=i,hits=total)),flush=True)
    con.executescript('CREATE INDEX mirna_region ON mirna_hits(region_id,end0); CREATE INDEX mirna_seed ON mirna_hits(seed_dna);')
    con.commit();con.close()
    dump(out/'mirna_complete.json',dict(utr_regions=len(records),seeds=len(seeds),hits=total,score='sequence_match_only_not_TargetScan_context_score'))


def scan_rbp(out):
    con=sqlite3.connect(out/'catalog.sqlite')
    con.execute('CREATE TABLE rbp_hits(motif_id TEXT,region_id TEXT,start0 INTEGER,end0 INTEGER,score REAL,pvalue REAL,matched_sequence TEXT)')
    command=[FIMO,'--text','--norc','--thresh','0.0001','--bfile',str(out/'background.txt'),str(out/'mouse_direct.meme'),str(out/'sequences.fa')]
    dump(out/'fimo_command.json',dict(argv=command,version=subprocess.check_output([FIMO,'--version'],text=True).strip(),qvalues_computed=False))
    total=0;batch=[]
    with open(out/'fimo.stderr.log','w') as log:
        proc=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=log,text=True)
        reader=csv.DictReader(proc.stdout,delimiter='\t')
        required={'motif_id','sequence_name','start','stop','strand','score','p-value','matched_sequence'}
        if not required.issubset(reader.fieldnames or []):
            proc.terminate();proc.wait();raise RuntimeError(str(reader.fieldnames))
        try:
            for r in reader:
                if not r.get('motif_id') or r['motif_id'].startswith('#'):
                    continue
                assert r['strand']=='+'
                a,b=int(r['start'])-1,int(r['stop'])
                assert b-a==len(r['matched_sequence'])
                batch.append((r['motif_id'],r['sequence_name'],a,b,float(r['score']),float(r['p-value']),r['matched_sequence']))
                if len(batch)>=20000:
                    con.executemany('INSERT INTO rbp_hits VALUES(?,?,?,?,?,?,?)',batch);total+=len(batch);batch=[];con.commit()
            assert proc.wait()==0,'FIMO failed; see stderr'
        except BaseException:
            if proc.poll() is None:proc.terminate();proc.wait()
            raise
    con.executemany('INSERT INTO rbp_hits VALUES(?,?,?,?,?,?,?)',batch);total+=len(batch)
    con.executescript('CREATE INDEX rbp_region ON rbp_hits(region_id,end0); CREATE INDEX rbp_motif ON rbp_hits(motif_id);')
    con.commit();con.close()
    dump(out/'rbp_complete.json',dict(hits=total,scope='mouse_direct_6_RBPs_7_available_experimental_PWMs_no_expression_filter',pvalue='per_match_threshold_not_regulatory_evidence',text_stream_no_top_score_cap=True))


def finalize(out):
    con=sqlite3.connect(out/'catalog.sqlite')
    checks={}
    checks['sqlite_integrity']=con.execute('PRAGMA integrity_check').fetchone()[0]=='ok'
    checks['no_single_pas_genes']=con.execute('SELECT count(*) FROM genes WHERE n_pas<2').fetchone()[0]==0
    checks['all_pas_have_windows']=con.execute('SELECT count(*) FROM pas p LEFT JOIN regions r ON r.region_id=p.window_id WHERE r.region_id IS NULL').fetchone()[0]==0
    for table in ['mirna_hits','rbp_hits']:
        checks[table+'_valid_bounds']=con.execute(f'SELECT count(*) FROM {table} h LEFT JOIN regions r USING(region_id) WHERE r.region_id IS NULL OR h.start0<0 OR h.end0>r.end0-r.start0 OR h.start0>=h.end0').fetchone()[0]==0
    fa=pysam.FastaFile(str(out/'sequences.fa'))
    genome=pysam.FastaFile(str(REF/'fasta/genome.fa'))
    checks['sample_strand_sequences']=True
    for strand in ['+','-']:
        for kind in ['pas_window','terminal_utr_plus2kb']:
            for rid,chrom,lo,hi in con.execute('SELECT region_id,chrom,start0,end0 FROM regions WHERE strand=? AND region_type=? ORDER BY region_id LIMIT 30',(strand,kind)):
                s=genome.fetch(chrom,lo,hi).upper();s=s if strand=='+' else rc(s)
                checks['sample_strand_sequences'] &= s==fa.fetch(rid)
    checks['sample_fimo_sequence']=all(fa.fetch(rid,a,b)==s for rid,a,b,s in con.execute('SELECT region_id,start0,end0,matched_sequence FROM rbp_hits ORDER BY rowid LIMIT 300'))
    con.executescript('''
    CREATE VIEW mirna_pas_containment AS SELECT p.study_id,p.atlas_pas_id,p.primary_gene_id AS gene_id,h.*,
      u.cut_class,CASE WHEN h.end0<=u.min_included_length THEN 1 ELSE 0 END AS fully_included_for_all_cluster_cuts
      FROM pas p JOIN pas_utr u USING(study_id,atlas_pas_id) JOIN mirna_hits h USING(region_id);
    CREATE VIEW rbp_pas_windows AS SELECT p.study_id,p.atlas_pas_id,p.primary_gene_id AS gene_id,p.primary_gene_symbol AS gene_symbol,
      h.motif_id,h.start0-r.cut_offset0 AS relative_start,h.end0-r.cut_offset0 AS relative_end,h.score,h.pvalue
      FROM pas p JOIN regions r ON r.region_id=p.window_id JOIN rbp_hits h ON h.region_id=p.window_id;
    ''')
    summary=dict(status='PASS' if all(checks.values()) else 'FAIL',checks=checks,
        counts={t:con.execute('SELECT count(*) FROM '+t).fetchone()[0] for t in ['genes','pas','regions','utrs','pas_utr','mirna_hits','rbp_hits']},
        unique_candidate_genes=con.execute('SELECT count(DISTINCT gene_id) FROM genes').fetchone()[0],
        genes_with_reference_utr=con.execute('SELECT count(DISTINCT gene_id) FROM utrs').fetchone()[0],
        analysis_effects_used=False,human_started=False,
        limitations=['Mouse direct motif subset, not full RBP coverage','No expression/CLIP validation or enrichment test','TES extension and alternative branch compatibility require downstream review','Transcript-overlapping hits are not independent or additive'])
    pd.read_sql_query('SELECT r.gene_id,r.region_type,count(*) AS motif_hits FROM rbp_hits h JOIN regions r USING(region_id) GROUP BY r.gene_id,r.region_type',con).to_csv(out/'rbp_region_summary.tsv.gz',sep='\t',index=False)
    con.commit();con.close()
    dump(out/'validation.json',summary)
    print(json.dumps(summary),flush=True)
    if not all(checks.values()):raise RuntimeError('Catalog validation failed')
    (out/'release.ok').write_text('MOUSE_EFFECT_INDEPENDENT_CATALOG_DIRECT_MOTIF_V1\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--stage',choices=['prepare','resume-prepare','mirna','rbp','finalize','all'],default='all')
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    if args.stage=='resume-prepare':
        resume_prepare(args.out);sys.exit(0)
    stages={'prepare':prepare,'mirna':scan_mirna,'rbp':scan_rbp,'finalize':finalize}
    for name,func in stages.items():
        if args.stage in ['all',name]:
            print(json.dumps(dict(stage=name,status='started')),flush=True);func(args.out)
