# Six-study unified RNA QC production

## Frozen scope

- WOT/GSE115943: 126 libraries, including four DiPSC endpoint libraries.
- GSE103221: 21 libraries across the CIP and OSK routes.
- GSE118258: 6 libraries from D0 through TRA-1-60 positive/negative D16.
- GSE150311: 4 libraries covering pooled D0-D7 intermediates and three endpoint cultures.
- GSE171127: 4 libraries covering pooled D6/D12 and route-specific D3 GETM/OSKM.
- GSE178325: 32 libraries; recomputed in production and checked against the accepted v3.2 membership hashes.

Total: 193 libraries and 778,207 Cell Ranger called cells. GSE122662 is not
included because it is absent from the aligned expression inventory used by the
current WOT object. It must not be reported as processed in this release.

## Execution contract

Every library is one Snakemake job. Basic RNA QC alone determines primary,
sensitivity-only, and technical-exclude membership. Scrublet is diagnostic-only
and cannot remove cells, change tiers, or change a sample gate. The preflight
rule validates all H5 shapes, called-cell counts, reviewed metadata, and exact
shared method parameters before any library job can start.

The production command uses 96 CPU slots and a 512 GB scheduler memory pool.
Each library job requests 12 GB, so memory limits concurrent jobs to at most 42.
Do not use `--keep-going`: a failed library or validator is a stop condition.

The frozen wrapper runs preflight and dry-run before starting production, refuses
an active duplicate PID, and keeps Snakemake's greedy `--resources` argument last:

```bash
bash production/run_96core.sh start
bash production/run_96core.sh status
```

```bash
ROOT=/root/nfsdata/apadetect/analysis_projects/reprogramming_unified_rna_qc_v3.2_release
cd "$ROOT"

SNAKEMAKE=/root/miniconda3/envs/single-cell-qc/bin/snakemake

"$SNAKEMAKE" --snakefile production/Snakefile \
  --config run_mode=resume \
  --cores 1 \
  preflight --printshellcmds

"$SNAKEMAKE" --snakefile production/Snakefile \
  --config run_mode=resume \
  --cores 96 --jobs 96 \
  --rerun-incomplete --latency-wait 120 \
  --printshellcmds --dry-run \
  --resources mem_mb=512000

nohup "$SNAKEMAKE" --snakefile production/Snakefile \
  --config run_mode=resume \
  --cores 96 --jobs 96 \
  --rerun-incomplete --latency-wait 120 \
  --printshellcmds --resources mem_mb=512000 \
  > production_qc_v3.2.log 2>&1 &
```

Production is accepted only when all six `independent_validation.tsv` files
contain no failure and every contact sheet has been visually reviewed. A study's
clustering must be rerun whenever its primary membership hash differs from the
membership used by the previous clustering release.

After visual review, update `production/manual_visual_review.tsv` and run:

```bash
bash production/run_96core.sh acceptance
```

The independent acceptance checker revalidates per-library completion, called-cell
closure, diagnostic-only Scrublet invariants, all three tier membership counts and
hashes, tier disjointness, independent validators, and contact-sheet coverage.
