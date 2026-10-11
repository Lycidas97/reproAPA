"""Reseal the v11b source-data and source-code archives (2026-10-11, MIT release).
Adds the in-place main-figure revision scripts and receipts of this directory, the MIT LICENSE and the
third-party notice, refreshes every pinned file to its current bytes, and rebuilds both archives with the
project packager. The data archive binds the new code archive. Pre-reseal archives are kept in
before-archives-mit/. Nothing is overwritten unless both archives pass the packager's member checks."""
import importlib.util, json, shutil, zipfile
from pathlib import Path

RUN = Path(__file__).resolve().parent
ROOT = RUN.parents[1]
spec = importlib.util.spec_from_file_location('packager', ROOT/'github/reproapa-analysis/scripts/package_source_data.py')
pack = importlib.util.module_from_spec(spec); spec.loader.exec_module(pack)
NAMES = ('Supplementary_Source_Code_v11b.zip', 'Supplementary-Source-Data-20261011-v11b.zip')
STAGE = RUN/'repacked-mit'; BACKUP = RUN/'before-archives-mit'; WORK = RUN/'work'
OWNER = 'main-figure caption revision'
P = RUN.relative_to(ROOT).as_posix() + '/'


def artifacts(a):
    return (a, a.with_suffix('.manifest.tsv'), a.with_suffix('.verification.json'), a.with_suffix('.zip.sha256'))


def selection(a):
    with zipfile.ZipFile(a) as z:
        return json.loads(z.read('SOURCE-SELECTION.json'))


def bind(rel):
    p = ROOT/rel
    with zipfile.ZipFile(p) as z: n = len(z.namelist())
    return dict(path=rel, bytes=p.stat().st_size, sha256=pack.digest(p), members=n)


def refresh(sel, additions):
    rows = {r['path']: r for r in sel['files']}
    for rel, role in additions:
        assert (ROOT/rel).is_file(), rel
        rows.setdefault(rel, dict(path=rel, archive_path='project/'+rel, roles=[role], owners=[OWNER]))
    for rel, r in rows.items():
        p = ROOT/rel; r.update(bytes=p.stat().st_size, sha256=pack.digest(p))
    sel['files'] = sorted(rows.values(), key=lambda r: r['path'])
    if 'source_files' in sel['coverage_counts']:
        sel['coverage_counts']['source_files'] = len(rows)


def build(sel, name):
    WORK.mkdir(exist_ok=True)
    pin = WORK/(name+'.selection.json'); pin.write_text(json.dumps(sel, indent=2)+'\n')
    out = STAGE/name; pack.package(ROOT, pin, out)
    for p in artifacts(out): shutil.copy2(p, ROOT/'table'/p.name)


LICENSE_NOTE = ('\nLicense: the authors\' code, configurations, processed source data and\n'
                'documentation in this archive are released under the MIT License (LICENSE).\n'
                'Third-party components keep their own terms (docs/third-party.md).\n')

if __name__ == '__main__':
    assert not BACKUP.exists() and not STAGE.exists(), 'earlier run present'
    code, data = [selection(ROOT/'table'/n) for n in NAMES]
    BACKUP.mkdir()
    for n in NAMES:
        for p in artifacts(ROOT/'table'/n): shutil.copy2(p, BACKUP/p.name)
    workbook = pack.digest(ROOT/'table/Supplementary_Tables_S1-S10.xlsx')
    historical = bind('table/Supplementary-Source-Data-20261009.zip'); s9 = bind('table/Supplementary_Source_Code_S9.zip')
    try:
        code['scope'] += ' Includes the in-place main-figure revision of Figures 1-4 and the MIT license.'
        code['companion'] = historical; code['additional_companions'] = [s9]
        code['archive_readme'] += (
            '\nThe current Figures 1-4 are in-place edits of the v11 production figures\n'
            '(notes/2026-10-11-caption-figure-fixes/). fix_fig1.py-fix_fig4.py apply\n'
            'label, panel-order and colour edits through fixlib.py, which uses the\n'
            'archived tools-local and pdfedit helpers. captions.py and apply_text_v2.py\n'
            'carry the matching caption and citation edits. No values were recomputed.\n'
            + LICENSE_NOTE)
        refresh(code, [(P+f, 'figure_revision_code') for f in
                       ('fixlib.py', 'fix_fig1.py', 'fix_fig2.py', 'fix_fig3.py', 'fix_fig4.py', 'diffcheck.py')]
                      + [(P+'captions.py', 'caption_revision_code'), (P+'apply_text_v2.py', 'caption_revision_code'),
                         (P+'repackage.py', 'packaging_entrypoint'),
                         ('LICENSE', 'license'), ('github/reproapa-analysis/LICENSE', 'license'),
                         ('github/reproapa-analysis/docs/third-party.md', 'third_party_notice')])
        build(code, NAMES[0])

        data['scope'] += ' Main Figures 1-4 after the caption review, and the MIT license.'
        data['companion'] = historical; data['additional_companions'] = [bind('table/'+NAMES[0]), s9]
        data['archive_readme'] += (
            '\nFigures 1-4 were revised in place after the main-figure caption review\n'
            '(Figure 2 rows reordered to follow first citation, Figure 3 d/e swapped,\n'
            'Figure 4 recoloured, label fixes). Values are unchanged. Edit receipts and\n'
            'the record are in notes/2026-10-11-caption-figure-fixes/, and\n'
            'Fig/figure-manifest.tsv pins the production assets.\n' + LICENSE_NOTE)
        refresh(data, [(P+'README.zh.md', 'figure_revision_record')]
                      + [(P+f'out/receipt.fig{i}.json', 'figure_revision_receipt') for i in (1, 2, 3, 4)]
                      + [(P+'out/apply_text_v2.log.json', 'caption_revision_receipt'),
                         ('Fig/figure-manifest.tsv', 'figure_provenance_guide'), ('LICENSE', 'license')])
        assert pack.digest(ROOT/'table/Supplementary_Tables_S1-S10.xlsx') == workbook
        build(data, NAMES[1])
    except BaseException:
        for p in BACKUP.iterdir(): shutil.copy2(p, ROOT/'table'/p.name)
        raise
    assert pack.digest(ROOT/'table/Supplementary_Tables_S1-S10.xlsx') == workbook
    rec = []
    for n in NAMES:
        r = json.loads((ROOT/'table'/n).with_suffix('.verification.json').read_text())
        rec.append(dict(archive=n, sha256=r['sha256'], members=r['members']))
    (WORK/'repackage-result.json').write_text(json.dumps(dict(status='PASS', workbook_sha256=workbook, archives=rec,
                                                               pre_reseal=str(BACKUP.relative_to(ROOT))), indent=2)+'\n')
    print(json.dumps(rec, indent=1))
