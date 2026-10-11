"""Result-record access for the five main-figure builders (records adapter).

Every number a builder draws comes from one ResultBundle through `Run.path/df/rows`, addressed by the
fixed component key `<panel>:<table>` of panel-data-interface.tsv. The builder keeps the locked master
canvas and the production layout; when the new data do not fit that layout (axis range, row capacity,
missing example, changed claim) the build still finishes, and the problem is written to the receipt
as a review item instead of being hidden or asserted against old values.
"""
import argparse, csv, gzip, hashlib, json, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p/'Fig').is_dir() and (p/'main.tex').exists())
sys.path.insert(0, str(ROOT/'analysis_code/pas_rerun'))
from result_records import ResultBundle   # noqa: E402
from pipeline import digest, write         # noqa: E402

ADAPTER = HERE.parent
# Components that the records adapter's contract added after a bundle was built (82- vs 81-component interface).
# Only these keys may be read through the table manifest when the bundle's result-records.json lacks them.
OUT_OF_CONTRACT = {
    'F3C:recurrence/recurrence_members.tsv.gz':
        'added to panel-data-interface by the 2026-10-05 records adapter (82 components); this bundle was built under the 81-component interface',
    'F5F:external_mirna/external_seed_expression.tsv':
        '5F site rug restricted to PSC-enriched seeds (user decision 2026-10-06); the S24 seed-expression table is in the bundle manifest but not bound to F5F',
    'F5F:rbp_expression/state_expression.tsv':
        '5F transferred-motif RBP row limited to RBPs expressed in P of both routes (user decision 2026-10-06); table is in the bundle manifest but not bound to F5F',
}
# Locked master canvases (unchanged files; their SHA-256 is asserted by each builder).
MASTERS = ROOT/'notes/2026-10-03-polyapipe-integrated-reference/narrative-reconstruction/main-figures-v5/masters'
MASTER_F3 = ROOT/'revision-backups/before-polyapipe-candidate-20261003-075521-0700/Fig/Figure.3.pdf'
FIGURE_PANELS = {'F1': ['F1A', 'F1B', 'F1C', 'F1D'], 'F2': ['F2A', 'F2B', 'F2C', 'F2D', 'F2E', 'F2F'],
                 'F3': ['F3A', 'F3B', 'F3C', 'F3D', 'F3E'], 'F4': ['F4A', 'F4B', 'F4C', 'F4D', 'F4E', 'F4F'],
                 'F5': ['F5A', 'F5B', 'F5C', 'F5D', 'F5E', 'F5F']}
# Frozen GEX tables behind the master-retained 1A/1B artwork (identical in every source; checked, not redrawn).
FROZEN_MASTER_INPUTS = {
    'frozen/display_cells.tsv.gz': '10dcd2b4d68a', 'frozen/states.tsv': 'c34f1e61de46',
    'frozen/S2A_state_geometry.tsv': '11c16390251f', 'frozen/actual_state_pair_edges.tsv.gz': '8a6efb325032'}


def _pas_parts(pas_id):
    """(position, strand) of a PAS id. Handles the colon form '...:<strand>:<pos>' and the polyApipe form
    '[polyApipe:<route>:]<chrom>_<pos>_<f|r>' (f = plus, r = minus strand)."""
    parts = str(pas_id).split(':')
    if parts[-1].lstrip('-').isdigit() and len(parts) >= 2 and parts[-2] in ('+', '-'):
        return int(parts[-1]), parts[-2]
    chrom, pos, s = parts[-1].rsplit('_', 2)
    strand = {'f': '+', 'r': '-', '+': '+', '-': '-'}.get(s)
    if strand is None or not pos.isdigit():
        raise ValueError(f'unrecognised PAS id: {pas_id}')
    return int(pos), strand


def pas_pos(pas_id):
    return _pas_parts(pas_id)[0]


def pas_strand(pas_id):
    return _pas_parts(pas_id)[1]


class Run:
    def __init__(self, figure, argv=None):
        p = argparse.ArgumentParser(description=f'{figure} from one result bundle (records adapter)')
        p.add_argument('--records', required=True, help='statistics/bundle/result-records.json of ONE source')
        p.add_argument('--out', required=True, help='new candidate directory (may exist only as an empty dispatcher folder)')
        p.add_argument('--selection', default=str(ADAPTER/'display-selection.json'),
                       help='example-gene slots: pinned genes or rank rule; validated against this bundle')
        p.add_argument('--stage', default='all', choices=['all', 'draft', 'retile', 'finish'],
                       help='all (default) = draw, re-tile and write the receipt; the split stages exist for slow re-tiles (Figure 3)')
        p.add_argument('--parts', default=None, help='with --stage retile: comma-separated panels to re-tile in this call')
        a = p.parse_args(argv)
        self.figure = figure; self.panels = FIGURE_PANELS[figure]; self.stage = a.stage
        self.records = Path(a.records).resolve(); self.out = Path(a.out).resolve()
        self.out.mkdir(parents=True, exist_ok=True)
        self.rb = ResultBundle(self.records.parent)
        if a.stage not in ('retile', 'finish'):
            clash = [q.name for q in self.out.iterdir() if q.name not in ('command.log',)]
            if clash:
                raise SystemExit(f'--out must be a new directory (found {clash[:5]})')
        self.comp = {c['component_id']: c for c in self.rb.record['components'] if c['figure'] == figure}
        self.review = []; self.derived = []; self.subsets = []
        self.selection_path = Path(a.selection).resolve()
        self.selection = json.loads(self.selection_path.read_text()) if self.selection_path.exists() else {}
        self.chosen = {}; self._paths = {}
        if a.stage in ('retile', 'finish'):
            self._resume(a.stage, a.parts)
            raise SystemExit(0)

    # ---------------------------------------------------------------- reading
    def path(self, panel, table):
        key = f'{panel}:{table}'
        if key not in self._paths:              # each component is fingerprint-checked once per build
            if key not in self.rb.components and key in OUT_OF_CONTRACT:
                self._paths[key] = self._manifest_table(panel, table, key)
            else:
                self._paths[key] = self.rb.table(panel, key)
        return self._paths[key]

    def _manifest_table(self, panel, table, key):
        """Explicitly listed component that the bundle's (older) interface does not register: read only after the
        table-manifest SHA (itself bound to bundle-complete.json) and the file digest agree; reported as a review item."""
        rel = 'tables/'+table
        meta = self.rb.manifest.get(rel)
        if meta is None:
            raise SystemExit(f'{key}: not registered and not in the table manifest of this bundle')
        path = (self.rb.root/rel).resolve()
        h = digest(path)
        if h != meta['sha256']:
            raise SystemExit(f'{key}: file digest differs from the table manifest')
        self.rb.used[rel] = h
        self.flag(panel, 'out_of_contract_read', f'{key}: {OUT_OF_CONTRACT[key]}; read via table-manifest SHA {h[:12]}',
                  table=rel, sha256=h)
        return path

    def df(self, panel, table, **kw):
        import pandas as pd
        return pd.read_csv(self.path(panel, table), sep='\t', **kw)

    def rows(self, panel, table, pred=None):
        p = self.path(panel, table)
        with (gzip.open(p, 'rt') if str(p).endswith('.gz') else open(p)) as f:
            return [r for r in csv.DictReader(f, delimiter='\t') if pred is None or pred(r)]

    # ---------------------------------------------------------------- review items
    def flag(self, panel, kind, detail, **data):
        self.review.append(dict(panel=panel, kind=kind, detail=detail, **data))

    def axis(self, panel, name, values, lo, hi):
        """Locked axis [lo, hi]. Out-of-range data are clipped by the caller and reported for a layout decision."""
        vals = [float(v) for v in values if v is not None and v == v]
        if not vals:
            return True
        vmin, vmax = min(vals), max(vals)
        if vmin < lo-1e-12 or vmax > hi+1e-12:
            self.flag(panel, 'axis_range', f'{name}: data {vmin:.4g}..{vmax:.4g} outside locked axis {lo:g}..{hi:g}',
                      data_range=[vmin, vmax], axis=[lo, hi])
            return False
        return True

    def capacity(self, panel, name, n, cap, subset=False):
        """subset=True: the layout shows a ranked display subset by design (recorded, not a review item)."""
        if subset:
            self.subsets.append(dict(panel=panel, what=name, shown=min(n, cap), of=n))
            return min(n, cap)
        if n > cap:
            self.flag(panel, 'row_capacity', f'{name}: {n} items, layout holds {cap}; first {cap} drawn', n=n, capacity=cap)
        return min(n, cap)

    def unavailable(self, panel, what):
        self.flag(panel, 'unavailable', what)

    def claim(self, panel, statement, holds):
        if not holds:
            self.flag(panel, 'claim_changed', f'no longer holds with this source: {statement}')
        return holds

    def derive(self, panel, what):
        """Display transform still done in the renderer (registered gap; no new statistic)."""
        self.derived.append(dict(panel=panel, what=what))

    # ---------------------------------------------------------------- examples
    def example(self, slot, panel, candidates, key='gene_id'):
        """Choose an example gene for `slot` from statistics-side eligible rows `candidates` (list of dicts).
        display-selection.json either pins a gene (must be eligible here) or applies the rank rule."""
        spec = self.selection.get(slot, {'mode': 'rank'})
        ranked = sorted(candidates, key=lambda r: float(r.get('rank', 0) or 0))
        pick = None; how = spec.get('mode', 'rank')
        if how == 'pin':
            pin = spec['gene_id']
            pick = next((r for r in candidates if r.get(key) == pin), None)
            if pick is None:
                self.flag(panel, 'example_ineligible', f'{slot}: pinned {spec.get("gene_symbol", pin)} is not eligible in this source',
                          pinned=pin, eligible=[r.get(key) for r in ranked[:10]])
        else:
            pick = ranked[0] if ranked else None
        if pick is None:
            self.unavailable(panel, f'{slot}: no eligible example in this source')
        self.chosen[f'{slot}@{panel}'] = dict(mode=how, chosen=None if pick is None else {k: pick.get(k) for k in pick},
                                 eligible=len(candidates))
        return pick

    # ---------------------------------------------------------------- finish
    def touch_all(self):
        for cid, c in self.comp.items():
            self.path(c['panel'], cid.split(':', 1)[1])
        if self.figure == 'F1':
            for t, short in FROZEN_MASTER_INPUTS.items():
                h = self.rb.used.get('tables/'+t, '')
                if not h.startswith(short):
                    raise SystemExit(f'{t} differs from the frozen GEX behind the master-retained 1A/1B ({h[:12]} != {short}); '
                                     '1A/1B must be redrawn before this source can use Figure 1')

    def _resume(self, stage, parts):
        st = json.loads((self.out/'_state.json').read_text())
        if st['records_sha256'] != digest(self.records):
            raise SystemExit('result records changed since the draft stage')
        self.review, self.derived, self.chosen, self.subsets = st['review'], st['derived'], st['chosen'], st.get('subsets', [])
        self.selection_path = Path(st['selection']); self._paths = {}
        n = int(self.figure[1])
        sys.path.insert(0, str(HERE))
        from retile import retile
        if stage == 'retile':
            retile(n, d=self.out, parts=parts.split(',') if parts else None, assemble=False)
            print(self.figure, 'retiled parts', parts)
        else:
            self._complete(n, self.out/f'Figure.{n}.draft.pdf', json.loads((self.out/'panel-values.json').read_text()), prebuilt=True)

    def finish(self, fig_number, draft, values):
        """draft: Figure.N.draft.pdf in self.out. Writes the retiled final, previews, values and receipt."""
        self.touch_all()
        if self.stage == 'draft':
            write(self.out/'panel-values.json', values)
            write(self.out/'_state.json', dict(records_sha256=digest(self.records), review=self.review, derived=self.derived,
                                              chosen=self.chosen, subsets=self.subsets, selection=str(self.selection_path)))
            print(self.figure, 'draft written; continue with --stage retile [--parts ...] and --stage finish')
            return None
        return self._complete(fig_number, draft, values)

    def _complete(self, fig_number, draft, values, prebuilt=False):
        self.touch_all()
        sys.path.insert(0, str(HERE))
        from retile import retile, check
        conflict = check(fig_number, d=self.out)
        rep = (retile(fig_number, d=self.out, parts=[] if prebuilt else None, assemble=True)) if not conflict else None
        if conflict:
            self.flag('all', 'letter_slot', f'content inside panel-letter slots {conflict}; final not retiled')
        final = self.out/f'Figure.{fig_number}.pdf'
        arts = [draft, self.out/f'Figure.{fig_number}.draft.png']
        if final.exists():
            subprocess.run(['pdftoppm', '-r', '150', '-png', '-singlefile', str(final), str(self.out/f'Figure.{fig_number}')], check=True)
            arts += [final, self.out/f'Figure.{fig_number}.png']
        write(self.out/'panel-values.json', values)
        arts.append(self.out/'panel-values.json')
        receipt = self.rb.provenance()
        builder = ADAPTER/f'figure{fig_number}'/'build.py'
        receipt.update(figure=self.figure, panels=self.panels, builder=str(builder.relative_to(ROOT)), builder_sha256=digest(builder),
                       selection=dict(file=str(self.selection_path), sha256=digest(self.selection_path) if self.selection_path.exists() else None,
                                      slots=self.chosen),
                       layout='needs_layout_review' if self.review else 'in_place', review=self.review,
                       renderer_display_derivations=self.derived, display_subsets=self.subsets, retile=None if rep is None else dict(new_height=rep['new_height']),
                       publication_review='pending', artifacts=[dict(path=a.name, sha256=digest(a)) for a in arts if a.exists()])
        write(self.out/'render-receipt.json', receipt)
        print(self.figure, receipt['layout'], len(self.review), 'review item(s)')
        return receipt
