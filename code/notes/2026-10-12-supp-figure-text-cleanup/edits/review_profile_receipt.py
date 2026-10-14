TARGET = 'notes/2026-10-10-supplement-v11/figures/profile.py'
# 2026-10-12 review fix (reviews/2026-10-12-supp-figure-revision-check, item 5): the receipt recorded the a-c source as
# 'Fig/Figure.2.pdf' although F2 is the pinned pre-v11 backup, and the d source (V11_S20_D) as an absolute sandbox path.
# Both are now written relative to the project root. The PDF itself is unchanged.
EDITS = [
 ("rec = dict(figure=F.name, page_pt=[PAGE_W, H], pdf_sha256=sha(out), redrawn=False,\n"
  "           kept=[dict(source='Fig/Figure.2.pdf', source_sha=sha(F2),",
  "def rel(p): return os.path.relpath(Path(p).resolve(), ROOT.resolve())\n"
  "rec = dict(figure=F.name, page_pt=[PAGE_W, H], pdf_sha256=sha(out), redrawn=False,\n"
  "           kept=[dict(source=rel(F2), source_sha=sha(F2),"),
 ("dict(source=str(S7), source_sha=sha(S7),",
  "dict(source=rel(S7), source_sha=sha(S7),"),
]
