# CLAUDE.md — Codebase: NGS.analysis v3

Conventions for working on this repository. The root `hpc_work/CLAUDE.md`
(environment, SSH, micromamba envs, SLURM rules, README convention,
shared-tool rules) still applies in full and is loaded automatically. This
file only adds what is specific to this codebase or differs from the root
defaults.

## Project

**What this is**: version 3 of the `NGS.analysis` primary-processing
pipeline (pypiper + looper, submitted to SLURM). It maps and QCs ChIP-seq,
ATAC-seq, CUT&RUN (`CR`), CUT&Tag (`CT`) and RNA-seq, in `genes` or
`repeats` mode. `README.md` describes usage, `CHANGELOG.md` lists every
change from v2.

**Goals of v3** (the design constraints behind most decisions):

1. **Fix the v2 data-losing bugs**: Trimmomatic palindrome mode without
   `keepBothReads` dropped read 2 of every fragment shorter than the read
   length, and `.dedup.unique.bam` still contained duplicates and chrM.
2. **genes mode is modernised**: ENCODE-style `<s>.filt.bam`, CPM bigwig,
   MACS3 peaks + FRiP, NRF/PBC, corrected TSS score.
3. **repeats mode reproduces v2 exactly**, except for the trimming fix.
   Its known problems are documented in
   `docs/repeat_strategy_open_questions.md` and **deliberately not fixed**.
   A future redesign of repeat read assignment is planned but not started.
4. **Statistics have one unit convention**: counts in fragments (PE pairs,
   SE reads), percentages 0–100. Every key is defined in
   `docs/stats_units.md`.

| Version | Code | Branch | Env | Status |
|---|---|---|---|---|
| v2 | `/store24/project24/becgsc_001/coding/NGS.analysis/` | `main` | `ngs.v2` | **production**, used by live projects (`01.Angela`, …) |
| v3 | `/store24/project24/becgsc_001/coding/NGS.analysis.v3/` (this folder) | `v3` | `ngs.v3` | in development / validation |

Both are clones of the same GitHub repo (`GunnarSchotta/NGS.analysis`).
The root CLAUDE.md's "Primary NGS processing" section still describes v2;
do not update it to v3 until the user decides v3 replaces v2.

## Data and output locations

| What | Path (HPC) | Status |
|---|---|---|
| This codebase | `/store24/project24/becgsc_001/coding/NGS.analysis.v3/` | edit here (with approval, see rules) |
| v2 code | `/store24/project24/becgsc_001/coding/NGS.analysis/` | **read-only** for v3 work |
| Validation on real data (T1, T2) | `/store24/project24/becgsc_001/analysis/NGS.analysis.v3.validation/` | write test runs here |
| Simulation benchmark | `/store24/project24/becgsc_001/analysis/NGS.analysis.v3.benchmark/` | write benchmark runs here |
| `ngs.v3` env | `/store24/project24/becgsc_001/micromamba/envs/ngs.v3/` | rebuild only via `env/build_ngs.v3.sh` |
| Genome references | `/store24/project24/becgsc_001/genomes/{mm10,hg38,hs1,mhaESC}/` | shared, read-only except when building a new genome |

Test and benchmark outputs never go into the repo (`example/results/` is
gitignored). The validation and benchmark folders each have a `README.md`
with question, steps, outputs and results; update them when a run changes
a conclusion.

The T2 reference for genes mode is the Edenhofer re-processing
(`analysis/Edenhofer/02.retrim.cutadapt/`). Read it, never write there.

## Repository layout

```
NGS.analysis.v3/
├── NGS.analysis.py                 sample pipeline (pypiper); __version__ here
├── NGS.analysis.collator.py        project pipeline (looper runp) -> NGS.summarizer.R
├── NGS.peaks.py                    peaks vs control (separate looper run)
├── ngs_qc.py, pyTssEnrichment.py   QC helpers (insert size, NRF/PBC, TSS score, plots)
├── generate_report.py              HTML report from pipestat results
├── check_project.py                preflight check, run before looper run
├── NGS.shiny.app.R                 interactive viewer (unchanged from v2)
├── *_pipeline_interface.yaml       looper interfaces (sample / project / peaks)
├── pipestat_*_schema.yaml          pipestat result schemas (sample / project / peaks)
├── NGS.analysis.yaml               pypiper tool config (bare names, resolved in ngs.v3)
├── genomes/<genome>.yaml           genome resources (indices, SAF, TSS, blacklist, mito, canonical, macs_gsize)
├── resources/genomes/              build scripts for T2T genome resources
├── env/                            ngs.v3 build script + pinned package lists
├── docs/                           stats_units, repeat strategy open questions, repeat literature
├── tests/                          T1 / T2 validation scripts, test_schema_keys.py
├── benchmark/                      ART + adapter simulation, evaluate.py, summarize.py
└── example/                        minimal example project (RNA, mm10)
```

## Environment

Everything runs in `ngs.v3` (exact copy of `ngs.v2` — conda, pip and R
packages — plus MACS3 3.0.4 and gffread):

```bash
E=/store24/project24/becgsc_001/micromamba/envs/ngs.v3/bin
export PATH=$E:$PATH          # looper, python3, samtools, bowtie2, STAR, macs3, ...
```

- The interfaces call `$E/python3` by absolute path and find scripts via
  `{looper.piface_dir}`. **No other installation paths are hardcoded** in
  interfaces or pipeline code; keep it that way.
- `NGS.analysis.py` puts the directory of its python first on `PATH`, so
  `NGS.analysis.yaml` uses bare tool names.
- Do not `pip install` / `micromamba install` into `ngs.v3` ad hoc. If a
  package is needed, change `env/build_ngs.v3.sh` and the pinned lists
  (`env/ngs.v3.explicit.txt`, `env/ngs.v3.pip.txt`), and ask first.
  `ngs.v2` is never modified.

## Design invariants

Check a change against these before proposing it:

| Area | Invariant | How it is checked |
|---|---|---|
| repeats mode | Commands for alignment, filtering, bigwig, IAP coverage and featureCounts identical to v2. Only trimming differs (`--legacy-trim` restores v2 trimming, for validation only) | T1 comparison R1 in `tests/compare_T1.py` (command identity after path normalisation) |
| Statistics | Counts in fragments, percentages 0–100, unit stated in the schema | `docs/stats_units.md` |
| pipestat schema | Every key `NGS.analysis.py` reports exists in `pipestat_results_schema.yaml` (pipestat raises `ColumnNotFoundError` otherwise). The full schema applies to every sample, not filtered by protocol | `tests/test_schema_keys.py` |
| Project pipeline | Own `pipeline_name` (`NGS.analysis.project`) and schema; schema top-level key is `project:` | T1 runs |
| Genomes | All per-genome paths come from `genomes/<genome>.yaml`; CLI arguments override. A new genome = YAML + row in `genomes/README.md` (+ build script if derived) | `check_project.py` checks the files exist |
| Version guard | `NGS.analysis.py` refuses an output folder with results from another major version (`NGS.analysis.version` marker). Bump `__version__` and `CHANGELOG.md` together | — |
| NGS.peaks | Separate interface, never listed in the project's `.looper.yaml` (looper would start it with the main pipeline) | — |
| Sample names | No name may be a prefix of another at a `_` boundary (pipestat flag glob) | `check_project.py` |

## Validation workflow

Any change to pipeline behaviour gets validated before it is called done:

1. **Static**: `$E/python3 tests/test_schema_keys.py` (login node, seconds).
2. **T1** (8 samples × 200k fragments, all protocols, PE + SE, mm10 + hg38):
   `tests/T1/make_v3_runs.sh submit`, then `tests/compare_T1.py <validation>/T1`
   (via SLURM). R0 = v2 noise floor, R1 = repeats mode must reproduce v2,
   R2 = trimming effect, R3 = genes mode v2 vs v3.
3. **T2** (full depth) only when a change can affect real-data results:
   `tests/T2/make_T2_runs.sh`, `compare_T2_eden.sh` (genes mode vs
   Edenhofer `02.retrim.cutadapt`), `compare_T2_repeats.py` (repeats mode vs
   existing v2 outputs in `01.Angela`).
4. **Benchmark** (`benchmark/run_benchmark.sh <genome>`) only for changes
   to alignment or trimming settings.

`make_v3_runs.sh` deletes and recreates its run folders; do not point it at
anything other than the validation folder.

## Documentation

- `CHANGELOG.md`: every user-visible change (outputs, statistics, defaults,
  new files) under the current version heading.
- `README.md`: usage and the "What the pipeline does" table.
- `docs/stats_units.md`: every new or changed statistic.
- `docs/repeat_strategy_open_questions.md`: new evidence on repeat
  assignment goes in "Evidence collected so far"; literature in
  `docs/repeat_literature.md` (with DOIs).
- Style of the existing docs: short sentences, plain language, tables for
  parameter lists, numbers with units.

## SLURM template (validation / benchmark jobs)

```bash
#!/bin/bash
#SBATCH --job-name=ngsv3_jobname
#SBATCH --partition=slim16
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/store24/project24/becgsc_001/analysis/NGS.analysis.v3.validation/logs/ngsv3_jobname_%j.out
#SBATCH --error=/store24/project24/becgsc_001/analysis/NGS.analysis.v3.validation/logs/ngsv3_jobname_%j.err
set -euo pipefail
E=/store24/project24/becgsc_001/micromamba/envs/ngs.v3/bin
export PATH=$E:$PATH
CODE=/store24/project24/becgsc_001/coding/NGS.analysis.v3
```

(Use `.../NGS.analysis.v3.benchmark/logs/` for benchmark jobs.) Pipeline
jobs themselves are submitted by looper with the `compute:` block of each
interface (slim16; sample 16 cores / 64 GB / 1 day).

## Git

- Work on branch **`v3`**. Never commit to or merge into `main` without an
  explicit request: `main` is v2 production.
- Run git **on the HPC** (`ssh hpc "cd .../coding/NGS.analysis.v3 && git ..."`).
  The local SSHFS mount reports "dubious ownership", and the push
  credentials (`github-ngs` SSH alias) only exist there.
- Commit messages: `v3: <what>` for code, `docs: <what>` for docs only.

## Rules for Claude (codebase additions)

1. Root rule 10 applies: this is shared software. Get an explicit yes
   before changing pipeline code, interfaces, schemas or the env. Once
   approved, commit on `v3`, push to `origin v3` from the HPC and report
   the commit to the user.
2. Never edit `coding/NGS.analysis/` (v2) or the `ngs.v2` env.
3. Do not change repeats-mode commands, even to fix a documented issue in
   `docs/repeat_strategy_open_questions.md`. Those are decided together in
   the planned repeats redesign; add findings to that document instead.
4. Every new reported statistic: schema entry with unit, entry in
   `docs/stats_units.md`, `tests/test_schema_keys.py` passes.
5. No hardcoded paths to a user's home or to this checkout in pipeline
   code; use `{looper.piface_dir}`, `genomes/<genome>.yaml` or arguments.
6. Run pipelines and tests on the HPC (SLURM for anything beyond
   seconds), with outputs in `analysis/NGS.analysis.v3.{validation,benchmark}/`.
7. Update `CHANGELOG.md`, `README.md` and the validation/benchmark
   READMEs in the same change as the code, not afterwards.
