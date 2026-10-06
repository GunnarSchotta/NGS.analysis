# CLAUDE.md — NGS.analysis v3 development

Conventions for developing NGS.analysis v3. The root `hpc_work/CLAUDE.md`
(SSH, path mapping, SLURM rules, shared-tool rules) still applies and is
loaded automatically. This file covers what v3 development needs: where the
tools and resources are, how to test a change, and which questions are open.

## What v3 is

Version 3 of the pypiper/looper primary-processing pipeline (ChIP, ATAC,
CUT&RUN, CUT&Tag, RNA; `genes` or `repeats` mode). Usage is in
`README.md`, all changes from v2 in `CHANGELOG.md`.

- **genes mode** is modernised: v3 trimming, ENCODE-style `.filt.bam`, CPM
  bigwig, MACS3 peaks, FRiP, NRF/PBC, corrected TSS score.
- **repeats mode** is v2 plus the trimming fix. This is deliberate:
  everything else waits for the repeats redesign (see open questions).
- Statistics are counts in fragments and percentages on a 0–100 scale
  (`docs/stats_units.md`).

| | v2 (production) | v3 (development) |
|---|---|---|
| Code | `coding/NGS.analysis/` (branch `main`) | `coding/NGS.analysis.v3/` (branch `v3`, this folder) |
| Env | `ngs.v2` | `ngs.v3` |

Both are clones of `GunnarSchotta/NGS.analysis`. Never edit v2 or `ngs.v2`
from here.

## Tools: where they are

All paths are on the HPC (`/store24/project24/becgsc_001/...`, shortened
to `.../` below).

**`ngs.v3` env** (`.../micromamba/envs/ngs.v3/bin`): exact copy of `ngs.v2`
plus MACS3 and gffread. Call by full path or put it first on `PATH`; no
`micromamba activate`.

| Tool | Version | Used for |
|---|---|---|
| Trimmomatic | 0.39 | trimming (all modes) |
| bowtie2 | 2.5.4 | genes-mode chromatin alignment |
| STAR | 2.7.11b | repeats mode, RNA |
| RSEM | (ngs.v2) | RNA quantification |
| samtools | 1.22.1 | filtering, stats |
| Picard | 3.4.0 | MarkDuplicates |
| featureCounts (subread) | 2.1.1 | repeat family / element counts |
| bamCoverage (deepTools) | (ngs.v2) | bigwigs |
| MACS3 | 3.0.4 | peaks |
| gffread | | GFF3 → GTF for genome builds |
| looper / pipestat / peppy / pypiper | 2.1.1 / 0.13.1 / 0.40.8 / 0.15.1 | framework |
| python3, Rscript | | pipeline scripts, `NGS.summarizer.R` |

`NGS.analysis.yaml` uses bare tool names: `NGS.analysis.py` puts its own
python's directory first on `PATH`. The interfaces call
`.../envs/ngs.v3/bin/python3` and find scripts via `{looper.piface_dir}`.
Keep pipeline code free of other installation paths.

Changing the env: edit `env/build_ngs.v3.sh` and the pinned lists
(`env/ngs.v3.explicit.txt`, `env/ngs.v3.pip.txt`, R packages in
`env/ngs.v3.R_packages_from_v2.txt`) and ask first. No ad-hoc installs.

**Tools outside `ngs.v3`** (used by tests, benchmark or open questions):

| Tool | Path | Notes |
|---|---|---|
| ART (`art_illumina`) | `.../micromamba/envs/art/bin/` | read simulation for `benchmark/` |
| cutadapt 5.2 | `.../micromamba/envs/bioenv/bin/` | reference trimming (Edenhofer `02.retrim.cutadapt`); not in ngs.v3 |
| GenMap | `.../micromamba/envs/bioenv/bin/genmap` | mappability; index only for mm10 (`.../genomes/mm10/genmap/index`) |
| RepeatMasker | `repeatmasker_env` | `-species` broken, see root CLAUDE.md |
| liftOver | `liftover` env | coordinate conversion between builds |

**Not installed** (candidates from the open questions, test in a new env
first): Allo, SmartMap, T3E, Telescope, TEtranscripts/TElocal, SQuIRE,
SalmonTE, TEspeX.

## Genome resources

The pipeline reads `genomes/<genome>.yaml` (key definitions in
`genomes/README.md`). Data are in `.../genomes/<genome>/`.

| Genome | Assembly | Repeat annotation | Blacklist |
|---|---|---|---|
| `mm10` | GRCm38 (UCSC) | UCSC RepeatMasker (RepBase names) | ENCFF547MET |
| `hg38` | GRCh38 (UCSC) | UCSC RepeatMasker | GRCh38 unified |
| `hs1` | T2T-CHM13v2.0, PAR-masked chrY | UCSC hs1 RepeatMasker | none |
| `mhaESC` | mhaESC v1.1 + mT2T-Y v1.1 (C57BL/6) | Dfam names, incl. chrY | none |

- T2T builds: `resources/genomes/{download_t2t_sources.sh,prep_t2t.py,build_genome.sh}`.
  Each genome folder has a `PREP_REPORT.txt`.
- A new genome needs a YAML, a row in `genomes/README.md` and, if derived,
  a build script. `check_project.py` checks that the files exist.
- IAP gag coverage BEDs exist only for mm10 (`gag.{plus,minus}.15k*.bed`).

## Testing a change

| Step | Command | Where | Output |
|---|---|---|---|
| Schema keys | `$E/python3 tests/test_schema_keys.py` | login node | stdout |
| Preflight | `$E/python3 check_project.py` in a project folder | login node | stdout |
| T1: 8 samples × 200k fragments, all protocols, PE + SE | `tests/T1/make_v3_runs.sh submit`, then `tests/compare_T1.py <T1 dir>` | SLURM | `analysis/NGS.analysis.v3.validation/T1/` |
| T2: full depth (Edenhofer ATAC genes mode, 01.Angela repeats mode) | `tests/T2/make_T2_runs.sh`, `compare_T2_eden.sh`, `compare_T2_repeats.py` | SLURM | `.../validation/T2/` |
| Simulation benchmark | `sbatch benchmark/run_benchmark.sh <genome>`, then `benchmark/summarize.py` | SLURM | `analysis/NGS.analysis.v3.benchmark/` |

(`E=/store24/project24/becgsc_001/micromamba/envs/ngs.v3/bin`.)

- T1 comparisons: R0 = v2 vs v2 noise, R1 = v3 repeats mode with
  `--legacy-trim` must match v2 command for command, R2 = effect of the new
  trimming, R3 = genes mode v2 vs v3.
- `make_v3_runs.sh` deletes and recreates its run folders.
- Results and conclusions are in the `README.md` of the validation and
  benchmark folders. Update them when a run changes a conclusion.
- SLURM logs: `.../analysis/NGS.analysis.v3.{validation,benchmark}/logs/`,
  job names `ngsv3_<name>`, `%j` only (no `%x`).

**Status (2026-10-03):** T1, T2 and the benchmark on all four genomes are
done. v3 repeats mode reproduces v2 (full-depth family/element r ≥ 0.9998).
v3 genes mode matches the Edenhofer cutadapt re-processing within 0.3%.

## Open questions for further development

Full lists: `docs/repeat_strategy_open_questions.md` (known issues,
questions, evidence) and `docs/repeat_literature.md` (papers with DOIs and
what each means for us). Benchmark numbers:
`analysis/NGS.analysis.v3.benchmark/README.md`. New evidence goes into
those files.

### Repeat read assignment (main open topic)

Goal: the most reliable assignment of reads to **individual repeat
copies**. Repeats mode is frozen until these are decided together.

**Known issues in current repeats mode** (documented, not fixed on
purpose):
1. Duplicates are marked but not removed (BAM, bigwig, counts, CR/CT splits).
2. Duplicate detection is unreliable for multimappers, because the random locus depends on the read name.
3. featureCounts `-p` without `--countReadPairs` counts mates, so PE counts are about 2× fragments.
4. Repeat counting is unstranded (`-s 0`), also for stranded RNA.
5. No `-O`: overlapping RepeatMasker entries become `Unassigned_Ambiguity`.
6. STAR `--alignMatesGapMax 350` cuts nucleosomal ATAC/CUT&RUN fragments.
7. RSEM runs on a transcriptome BAM made with `--outSAMmultNmax 1`.
8. STAR `--alignEndsType EndToEnd`: no soft-clipping of residual adapter.
9. Short fragments recovered by the v3 trimming add many multimappers.
10. CR/CT split puts TLEN 0 records in the sub-nucleosomal BAM.
11. IAP normalisation counts STAR pairs, `bedtools coverage` counts mates.

**Evidence so far** (ART simulation from the reference itself, so an
upper bound):
- STAR "unique" (MAPQ 255) is wrong for 3–5% (2×100) and 6–7% (2×50) of
  young L1 fragments (L1MdT/L1MdA, L1HS); IAPEz 1.4–3%; other young ERVs
  and SVA ≤ 0.5%.
- bowtie2 MAPQ ≥ 30 keeps fewer fragments uniquely (10–47% of young
  ERV/L1) but places ≥ 99.6% correctly: a candidate stricter "unique"
  definition for element-level counts.
- Of everything repeats mode places, including random multimapper
  placement, 94–98% is at the correct locus.
- 2×50 roughly halves the unique fraction of young L1/IAP compared with 2×100.

**Decisions to make:**
- **Multimapper allocation**: random-1 (current), fractional 1/n, EM
  (Telescope, TEtranscripts, SQuIRE, SalmonTE), or signal-aware for
  chromatin (Allo: CNN, bowtie2 `-k 25`, trained on TF/ATAC and untested
  on broad marks; SmartMap: Bayesian, PE only). For family-level chromatin
  counts: T3E-style 1/n weighting with an input background.
- **Locus-level RNA**: EM tools rank best (Schwarz 2022), but locus-level
  false positives can outnumber true loci (Savytska 2022), so count filters
  and ideally TSS evidence are needed. Exonised TE fragments inflate
  counts (TEspeX).
- **Deduplication for multimappers**: sequence-level dedup before
  alignment, UMIs, or report with and without duplicates.
- **Per-copy mappability**: GenMap for the read/fragment lengths in use
  (k = 36–150): which families are addressable at element level at all?
- **Read length**: gain of PE and of 2×100 over 2×50 for young L1/IAP/MERVL.

### Reference genome and annotation

- **hs1 vs hg38**: outside satellites, same unique fractions and accuracy
  when reads come from the reference itself. hs1 adds unmappable satellite
  arrays.
- **mhaESC vs mm10**: outside satellites, mhaESC is 2–3 points less unique
  at the same accuracy. Its Dfam annotation splits the youngest L1
  subfamilies; for L1MdTf_I/II and L1MdA_I about 15% of STAR "unique"
  fragments are at the wrong copy.
- **Annotation comes with the reference**: RepBase (mm10) and Dfam
  (mhaESC) family names differ, so per-family results are not portable
  without a name mapping.
- **Alternative mouse T2T**: Francis et al. 2025 (C57BL/6J + CAST/EiJ) is
  not built. Choose mhaESC, Francis, or support both.
- **Not yet measured**: the cross-reference effect, i.e. reads from
  sequence missing in hg38/mm10 forced onto wrong paralogues. Needs a
  simulation from hs1/mhaESC aligned to hg38/mm10.

### Strain background (mouse)

ES lines are often 129 or mixed. SNPs/indels and non-reference TE
insertions (TEs are 75% of SV bases, Ferraj 2023) push reads onto
reference paralogues. Options: strain-specific assemblies (Lilue 2018,
Helmy 2025), masking known polymorphic loci, or at least flagging copies
in known SVs. Mismatch allowance (`--outFilterMismatchNmax 3`) is part of
this.

### Benchmark follow-ups

- Cross-reference simulation (above).
- Injected SNPs and indels (strain divergence).
- Family-level (not only locus-level) correctness of random placement.
- Alternative strategies on the ART set: bowtie2 `-k`, Allo, EM.
- Real data: Setdb1 KO RNA-seq and H3K9me3 ChIP (01.Angela), comparing
  family- and element-level results between strategies.
- Ground truth for copy-level claims: one Nanopore pilot (locus-specific
  methylation, Ewing 2020) before committing to an element-level strategy.

### Smaller open points outside repeats

- Trimmomatic `MINLEN:30` drops < 30 bp fragments that cutadapt `-m 20`
  keeps (about 0.2% of ATAC fragments, relevant for footprinting).
- Input normalisation for family-level ChIP counts is not done in any mode.

## Working rules

1. Root rule 10 applies: ask before changing pipeline code, interfaces,
   schemas or the env. Approved changes are committed on `v3` and pushed
   from the HPC (`ssh hpc "cd .../coding/NGS.analysis.v3 && git ..."`;
   the local mount reports "dubious ownership", and the push key is only
   on the HPC). Never commit to `main`.
2. A new reported statistic needs a schema entry with its unit, an entry
   in `docs/stats_units.md`, and `tests/test_schema_keys.py` must pass.
3. Do not change repeats-mode commands piecemeal. Record findings in
   `docs/repeat_strategy_open_questions.md` instead.
4. Test and benchmark outputs go into
   `analysis/NGS.analysis.v3.{validation,benchmark}/`, never into the repo.
5. Update `CHANGELOG.md`, `README.md` and the relevant docs in the same
   change as the code.
