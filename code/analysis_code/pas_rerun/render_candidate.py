"""Render named figures into the isolated candidate-artifact directory."""
import argparse

from render_common import Renderer
from render_overview import figure1
from render_profiles import figure2,supplement7,supplement8
from render_qc import supplement4
from render_recurrence import figure3,supplement9,supplement10,supplement11
from render_process import supplement12,supplement13,supplement14
from render_cross import figure4,supplement15,supplement16,supplement17
from render_cis import supplement18
from render_caller import supplement5
from render_global import supplement6
from render_regulatory import figure5,supplement19,supplement20,supplement21
from render_external import supplement22,supplement23

FIGURES={'1':figure1,'2':figure2,'3':figure3,'S4':supplement4,'S7':supplement7,'S8':supplement8,
    'S9':supplement9,'S10':supplement10,'S11':supplement11,'S12':supplement12,'S13':supplement13,'S14':supplement14,
    '4':figure4,'5':figure5,'S5':supplement5,'S6':supplement6,
    'S15':supplement15,'S16':supplement16,'S17':supplement17,'S18':supplement18,
    'S19':supplement19,'S20':supplement20,'S21':supplement21,'S22':supplement22,'S23':supplement23}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    p.add_argument('--figures',required=True,help='Comma-separated figure IDs; missing input receipts fail immediately')
    p.add_argument('--theme-supplement',action='store_true')
    a=p.parse_args()
    if a.theme_supplement:
        from render_theme import run
        run(a.config,a.out,a.figures)
        raise SystemExit(0)
    r=Renderer(a.config,a.out)
    for name in a.figures.split(','):FIGURES[name](r)
