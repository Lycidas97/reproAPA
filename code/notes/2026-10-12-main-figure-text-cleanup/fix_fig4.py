"""Figure 4 (2026-10-12 main-figure text cleanup).
4b  title 'Mouse P–M, 122 genes with a strict recurrent interval' -> 'Mouse P–M, 122 genes' (reference and
    denominator stay; the selection is defined in the caption).
4c  x-axis title and the key below it move down 2 pt: the title stood 0.3 pt under the tick labels.
4e  'ΔI' column head moves up 2 pt (it touched '−23 pp'); the key moves down 2 pt (0.5 pt under the tick labels).
4f  the '0' between '−50' and '+50' on both weighted-site axes is dropped (2.7 pt gaps); the dotted zero line marks 0."""
import json
from fixlib import Job, chars_of, INK
J = Job(4)
t4b = chars_of(J.src, 'Mouse', (18, 106, 45, 117))
J.remove((20.2, 108.3, 209.8, 116.9), '4b title (redrawn shorter)')
J.move_word('ΔI', (119, 353, 131, 364), -2.0, note='4e column head clear of −23 pp')
J.drop_word('0', (291, 402, 300, 412), note='4f miRNA weighted-site axis zero label')
J.drop_word('0', (427, 402, 436, 412), note='4f RBP weighted-site axis zero label')


def draw(t):
    t.text(t4b['origin'], t4b['base'], 'Mouse P–M, 122 genes', t4b['size'], color=t4b['color'])


rec = J.run(draw=draw, shifts=[
    ((340.0, 173.45, 453.5, 191.5), 2.0, '4c x-axis title and key'),
    ((0.0, 399.95, 136.0, 418.5), 2.0, '4e key'),
])
print(json.dumps({k: rec[k] for k in ('removed_kinds', 'words_removed', 'words_added', 'sha256')}, indent=1, ensure_ascii=False))
json.dump(rec, open('out/receipt.fig4.json', 'w'), indent=1, ensure_ascii=False)
