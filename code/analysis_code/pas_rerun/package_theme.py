"""Verify and package the current supplement without modifying production assets."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,html,json,subprocess,shutil
import pandas as pd
import pymupdf as fitz
from result_records import ResultBundle

ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

# Functions reuse data frames across logical blocks without reading them again.
# Record those dependencies explicitly; a file-read log alone is not a panel map.
GEX='root:notes/2026-09-24-integration-top-labels/'
REUSED_INPUTS={
    **{f'D{i:02d}B':('A',[GEX+'source_data/states.tsv',GEX+'palette.json']) for i in [1,2]},
    'D03B':('A',[GEX+'source_data/states.tsv']),
    'D04B':('A',['pas_qc/catalog_summary.tsv']),
    'D07C':('B',['profiles/architecture_gene_metrics.tsv.gz']),
    'D08B':('A',['moran/moran_results.tsv']),
    'D13C':('B',['intermediate/process/same_series_process_classes.tsv.gz']),
    'D14B':('A',['profile_gex/enrichment_all_tests.tsv.gz','profile_gex/rank_enrichment_all_tests.tsv.gz']),
    'D14C':('A',['profile_gex/rank_enrichment_all_tests.tsv.gz']),
    'D18B':('A',['same_pas/v2_primary_correlation_summary.tsv']),
    'D18C':('A',['same_pas/v2_primary_correlation_summary.tsv']),
    'D21B':('A',['site_statistics/all_changed_region_pairs.tsv.gz','site_statistics/region_consistency_statistics.tsv']),
    'D22C':('B',['supplement_statistics/partial-selected-site-pairs.tsv']),
    'D26B':('A',['examples/recurrent_region_pairs.tsv.gz']),
    'D26C':('A',['examples/recurrent_region_pairs.tsv.gz']),
    'D26E':('D',['geometry/axis-decomposition.tsv']),
    'D27B':('A',['geometry/axis-decomposition.tsv']),
    'D27C':('A',['geometry/axis-decomposition.tsv']),
}

def resolve_sources(receipt,binding_index):
    fid=receipt['figure'];resolved={}
    for logical,evidence in receipt['panels'].items():
        sources={}
        for ref,source in evidence['sources'].items():
            bound=binding_index[(fid,logical,ref)]
            assert source['path']==bound['path'] and source['sha256']==bound['sha256']
            sources[ref]=dict(bound,reused_from='')
        if logical in REUSED_INPUTS:
            producer,tables=REUSED_INPUTS[logical];producer=fid+producer
            for table in tables:
                found=[v for (f,p,_),v in binding_index.items() if f==fid and p==producer and v['table']==table]
                assert len(found)==1,(logical,producer,table)
                bound=found[0];sources[bound['design_reference']]=dict(bound,logical_panel=logical,reused_from=producer)
        assert sources,('Panel has no resolved input',logical)
        resolved[logical]=sources
    return resolved

def package(out):
    out=Path(out).resolve();provenance=json.loads((out/'bundle-provenance.json').read_text())
    config=json.loads((out.parent/'render.config.json').read_text());reader=ResultBundle(Path(config['result_records']).parent,config)
    assert provenance['source']==reader.identity
    contract=json.loads((out/'review-contract.json').read_text());captions=json.loads((out/'caption-contract.json').read_text())
    definitions={f['id']:f for f in contract['figures']};assert len(definitions)==27
    binding_table=pd.read_csv(out/'supplement-source-bindings.tsv',sep='\t')
    binding_index={(v['figure'],v['logical_panel'],v['design_reference']):v for v in binding_table.to_dict('records')}
    import render_theme_panels,render_theme_review,inspect
    scripts={};source_files=set()
    for fid,fn in render_theme_panels.DRAW.items():
        scripts[fid]=[f'{inspect.getsourcefile(fn)}:{inspect.getsourcelines(fn)[1]} ({fn.__name__})']
    for fid,fn in [('D05',render_theme_review.extra05),('D18',render_theme_review.extra18),('D26',render_theme_review.extra26),('D27',render_theme_review.D27)]:
        scripts.setdefault(fid,[]).append(f'{inspect.getsourcefile(fn)}:{inspect.getsourcelines(fn)[1]} ({fn.__name__})')
    combined=fitz.open();panels=[];bindings=[];summary=[];all_inputs={};toc=[];sections=[]
    layout_rows=[];resolved_bindings=[]
    for i in range(1,28):
        fid=f'D{i:02d}';r=json.loads((out/'verification'/f'{fid}.json').read_text());pdf=out/'Fig'/f'SuppFig.{fid}.pdf'
        resolved=resolve_sources(r,binding_index)
        resolved_bindings.extend(v for sources in resolved.values() for v in sources.values())
        assert digest(pdf)==r['pdf_sha256'] and not r['geometry_flags'],fid
        assert r['height_mm']<=190 and max(r['frame_counts'].values())<=2,fid
        assert all(v['max_y_edge_difference_mm']<.01 and v['max_y_limit_difference']<1e-8 for v in r['shared_y_checks']),fid
        doc=fitz.open(pdf);assert len(doc)==1 and not doc.is_repaired;page=doc[0]
        sizes=[s['size'] for b in page.get_text('dict')['blocks'] if b['type']==0 for l in b['lines'] for s in l['spans']]
        assert min(sizes)*150/160>=6.49,(fid,min(sizes))
        fonts=subprocess.run(['/opt/homebrew/bin/pdffonts',str(pdf)],capture_output=True,text=True,check=True)
        lines=fonts.stdout.splitlines()[2:];assert lines and all('NimbusSanL' in s and 'yes' in s for s in lines),fid
        (out/'verification'/f'{fid}.fonts.txt').write_text(fonts.stdout)
        for tag,dpi in [('placed',150),('placed-2x',300)]:
            target=out/'verification'/tag;target.mkdir(exist_ok=True)
            page.get_pixmap(matrix=fitz.Matrix(dpi/72*150/160,dpi/72*150/160)).save(target/f'SuppFig.{fid}.png')
        combined.insert_pdf(doc);toc.append([1,f'S{i} — {captions[fid][0]}',i]);doc.close()
        shutil.copy2(pdf,out/'Fig'/f'SuppFig.S{i}.pdf')
        caption='Figure S'+str(i)+'. '+captions[fid][0]+'. '+captions[fid][1]
        (out/'captions'/f'S{i}.en.txt').write_text(caption+'\n')
        for actual in r['panel_mapping']:
            n='S'+str(i)+actual['panel'][-1]
            row=dict(panel=n,design_panel=actual['panel'],figure=f'S{i}',topic=definitions[fid]['title'],
                source_logical_panels=';'.join(actual['source_logical_panels']),components=json.dumps(actual['components']),
                rectangle_mm=json.dumps(actual['rectangle_mm']),shared_y_basis=';'.join(actual['shared_y_basis']),
                drawing_script=';'.join(scripts[fid]),pdf=str(pdf),source_key=reader.identity['source_key'])
            panels.append(row)
            for logical in actual['source_logical_panels']:
                for ref,source in resolved[logical].items():
                    p=Path(source['path']);h=source['sha256']
                    if str(p) not in all_inputs:assert digest(p)==h,str(p);all_inputs[str(p)]=h
                    if source['source_key']!='frozen_GEX':
                        assert source['source_key']==reader.identity['source_key']
                        assert p.is_relative_to(reader.root/'tables')
                        assert reader.manifest[str(p.relative_to(reader.root))]['sha256']==h
                        assert 'v3' not in str(p.relative_to(reader.root))
                    bindings.append(dict(panel=n,logical_panel=logical,design_reference=ref,source_table=source['table'],
                        source_path=str(p),sha256=h,source_key=source['source_key'],drawing_script=row['drawing_script'],
                        reused_from=source['reused_from'],
                        dependency_scope='logical_block; exact displayed subsets in linked panel-data exports',
                        displayed_exports=';'.join(v['path'] for v in r['panels'][logical]['exports'])))
            envelope=json.loads((out/'verification'/f'{fid}.content-envelopes.json').read_text())
            selected=next(v for v in envelope if v['panel']==actual['panel'])
            layout_rows.append(dict(panel=n,**{k:v for k,v in selected.items() if k!='panel'}))
        for logical,ev in r['panels'].items():
            for exported in ev['exports']:
                assert digest(out/exported['path'])==exported['sha256'],exported['path']
        summary.append(dict(figure=f'S{i}',design_figure=fid,height_mm=r['height_mm'],panels=len(r['panel_mapping']),
            pdf_sha256=r['pdf_sha256'],bounds_flags=len(r['geometry_flags']),minimum_placed_font_pt=min(sizes)*150/160))
        sections.append(f'<section id="{fid}"><a id="S{i}"></a><h2>S{i} · {html.escape(definitions[fid]["title"])}</h2><p><a href="Fig/SuppFig.S{i}.pdf">矢量 PDF</a> · <a href="verification/placed-2x/SuppFig.{fid}.png">2× 校样</a> · <a href="captions/S{i}.en.txt">英文图注</a> · <a href="captions/{fid}.zh.txt">统计与来源说明</a></p><img loading="lazy" src="Fig/SuppFig.{fid}.png" alt="Figure S{i}"><details><summary>Figure legend</summary><p>{html.escape(caption)}</p></details></section>')
    combined.set_toc(toc);combined.save(out/'Supplement-polyApipe-figures.pdf',garbage=4,deflate=True);combined.close()
    pd.DataFrame(panels).to_csv(out/'supplement-panel-map.tsv',sep='\t',index=False)
    pd.DataFrame(bindings).drop_duplicates().to_csv(out/'supplement-panel-data-script-map.tsv',sep='\t',index=False)
    pd.DataFrame(resolved_bindings).drop_duplicates().to_csv(out/'supplement-resolved-source-bindings.tsv',sep='\t',index=False)
    pd.DataFrame(summary).to_csv(out/'figure-manifest.tsv',sep='\t',index=False)
    (out/'source-hashes.json').write_text(json.dumps(all_inputs,indent=2))
    (out/'verification'/'all-content-envelopes.json').write_text(json.dumps(layout_rows,indent=2))
    code=list((ROOT/'analysis_code/pas_rerun').glob('render_theme*.py'))+[ROOT/'analysis_code/pas_rerun'/n for n in ['render_candidate.py','render_common.py','package_theme.py','result_records.py']]
    (out/'renderer-code-hashes.json').write_text(json.dumps({str(p):digest(p) for p in code},indent=2))
    before=json.loads((out.parent/'production-before.json').read_text())
    changed=[p for p,h in before.items() if digest(ROOT/p)!=h];assert not changed,changed
    visual_path=out/'verification'/'visual-review.json';visual_status='PENDING'
    if visual_path.exists():
        visual=json.loads(visual_path.read_text());assert visual['source_key']==reader.identity['source_key']
        assert {v['figure'] for v in visual['figures']}==set(definitions)
        for v in visual['figures']:
            fid=v['figure'];receipt=json.loads((out/'verification'/f'{fid}.json').read_text())
            assert v['status']=='PASS' and v['pdf_sha256']==receipt['pdf_sha256'],fid
            assert v['panels_reviewed']==[p['panel'] for p in receipt['panel_mapping']],fid
            for tag,h in v['proof_sha256'].items():assert digest(out/'verification'/tag/f'SuppFig.{fid}.png')==h,fid
        visual_status='PASS'
    verification=dict(status='AUTOMATED_GATES_PASS',checked_at_utc=datetime.now(timezone.utc).isoformat(),source=reader.identity,
        figures=len(summary),panels=len(panels),unique_sources=len(all_inputs),panel_source_bindings=len(pd.DataFrame(bindings).drop_duplicates()),
        maximum_height_mm=max(v['height_mm'] for v in summary),geometry_flags=0,production_changed=changed,
        visual_review_file='verification/visual-review.json',visual_review_status=visual_status,
        statistical_tests_recomputed=False,old_PAS_fallback=False)
    (out/'verification'/'package-audit.json').write_text(json.dumps(verification,indent=2))
    nav=' '.join(f'<a href="#D{i:02d}">S{i}</a>' for i in range(1,28))
    header='<p>polyApipe common · v2 primary · 27 幅 / '+str(len(panels))+' 个子图。已接入当前完整统计包；保留冻结 GEX 框架。</p><p><a href="Supplement-polyApipe-figures.pdf">合并 PDF</a> · <a href="supplement-panel-data-script-map.tsv">子图—数据—脚本映射</a> · <a href="verification/package-audit.json">自动校验</a> · <a href="verification/visual-review.json">目视审查</a> · <a href="../REPORT.zh.md">接入报告</a></p><p>候选交付版；原生产图和正文保持可恢复的原版。示例：S22 为 Snx30，S26 为 Trip12 两个区间。</p>'
    (out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>polyApipe supplement v16</title><style>body{font:15px system-ui;color:#263d4a;max-width:1100px;margin:30px auto;padding:0 20px;background:#f5f7f8}a{color:#0072b2}nav{position:sticky;top:0;background:#f5f7f8;padding:12px;line-height:2;z-index:1}nav a{margin-right:13px}section{background:white;padding:24px;margin:26px 0;border:1px solid #d8e0e5;scroll-margin-top:120px}img{width:100%;max-width:945px;height:auto;display:block}details{margin:18px 0;line-height:1.6}h1{font-size:24px}h2{font-size:19px}</style><h1>polyApipe · Supplement v16</h1>'+header+'<nav>'+nav+'</nav>'+''.join(sections))
    print(json.dumps(verification,indent=2))
