# Changelog

## 3.1.0 (branch `v3`, 2026-11)

### Mitochondrial reads removed before `<s>.bam` (chromatin: ATAC, CHIP, CR, CT; both modes)
- Reads are still aligned against the full genome, chrM included, so MAPQ decides between chrM and nuclear NUMT copies.
  The aligner output then passes `samtools view -e 'rname != "chrM" && rnext != "chrM"'` before sorting: every record of
  a fragment with the read or its mate on chrM is split off (no orphan mates). This is not PEPATAC-style prealignment.
- Everything computed from `<s>.bam` is therefore nuclear only: Picard `Duplication_pct`, NRF/PBC, and in repeats mode
  `.dedup.unique.bam`, its bigwig (RPKM denominator), peaks, FRiP, insert size, TSS score and featureCounts.
  Genes-mode `.filt.bam` already excluded chrM; there only `Duplication_pct`, NRF and PBC change.
- Why: in ATAC the near-saturated chrM reads dominated the duplication and complexity statistics, and in repeats mode
  chrM peaks survived the blacklist and inflated FRiP. Vincendeau ATAC.HERVKC4 under 3.0.0: reported duplication
  64–70% vs about 43% for nuclear fragments; 42–64% of the in-peak fragments were on chrM peaks.
  Effect of 3.1.0 on the same data: `analysis/NGS.analysis.v3.validation/T3_mito/README.md`.
  Measured (T3_mito, 6 libraries): duplication 64–70% → 44–47%, NRF 0.32–0.38 → 0.55–0.58, FRiP 17–23% → 7–13%, chrM peaks 8–9 → 0 per sample, LTR14 non-duplicate fragments identical.
- New keys: `Mito_filter` (`removed` / `kept`), `Nuclear_fragments`. `Mito_fragments` now counts a fragment if either
  mate is on chrM; `Aligned_fragments` and alignment percentages still include chrM.
- bowtie2 `Aligned_fragments` / `Unique_fragments` / `Aligned_proper_pairs` are counted over `<s>.bam` plus the split-off BAM.
- RNA keeps chrM. `--keep-mito` (hidden; `keep_mito` sample attribute) restores the 3.0 behaviour; T1 R1 uses it with
  `--legacy-trim` so that repeats mode can still be compared with v2 command for command.
- Not changed: the repeats IAP normalisation (STAR log totals, which include chrM; known issue 11).

### Framework and environment
- The version guard compares major.minor: 3.1 refuses an output folder written by 3.0 (pypiper would reuse the old
  chrM-containing BAMs).
- Picard 3.4.0 → 3.5.0 (MarkDuplicates stores the physical location as int: no overflow in optical-duplicate detection).
  samtools stays 1.22.1 (the 1.22.4 backport is not on bioconda; 1.23+ deferred, roadmap S1).

### Documentation
- New `docs/roadmap.md`: prioritised open items, the decision against per-protocol sub-pipelines, package status.
- README: QC/peaks rows corrected (complexity on `<s>.bam`; repeats peaks `--keep-dup 1`), per-sample `genome_config`.
- `docs/repeat_strategy_open_questions.md`: issue 9 measured and closed, NUMT issue 12, T2T evidence moved to "collected".

## 3.0.0 (branch `v3`, 2026-10)

### Trimming (all modes)
- **PE:** Trimmomatic `ILLUMINACLIP:...:2:30:10:1:true`.
  - v2 used palindrome mode without `keepBothReads` and with `minAdapterLength` 8. That dropped read 2 for all fragments shorter than the read length (17–26% of pairs in ATAC data) and lost fragments of 53–59 bp entirely.
  - Validated against cutadapt `-O 1`: 99.88% identical fragments.
- **SE:** simple-clip threshold 10 → 7 (adapter matches of about 12 nt).
- Unpaired outputs are gzipped and deleted at cleanup; the Trimmomatic log is parsed for statistics.

### genes mode (chromatin)
- bowtie2 `--dovetail`, so fully clipped short pairs stay concordant.
- New analysis BAM `<s>.filt.bam`:
  - proper pairs, duplicates removed (`-F 1804`), MAPQ ≥ 30;
  - canonical chromosomes, mitochondria removed, blacklist removed;
  - orphan mates removed (name sort → `fixmate -r` → `-f 2`).
  - Replaces `<s>.dedup.unique.bam`, which in v2 still contained the duplicates and chrM.
- `<s>.filt.bw`: CPM, fragment-extended (v2: RPKM, not extended, from the duplicate-containing BAM).
- CR/CT nuc/subnuc split now from the filtered BAM; TLEN 0 records are no longer counted as sub-nucleosomal.

### repeats mode
- **Unchanged:** alignment, filtering, bigwig, IAP coverage (including its normalisation) and featureCounts commands are identical to v2. Known issues are documented, not fixed: `docs/repeat_strategy_open_questions.md`.
- **Added outputs:** peaks, FRiP, insert size, library complexity and the corrected TSS score.

### Statistics
- All counts in fragments, percentages 0–100, documented in `docs/stats_units.md`.
- New keys replace the ambiguous v2 keys: `Trimmed_fragments`, `Aligned_fragments`, `Unique_fragments`, `Mito_pct`, `Duplication_pct`, `Filtered_fragments`, …
- **Fixed:**
  - STAR "mapped" included "too many loci" reads;
  - mito % mixed records and pairs;
  - `Percent_duplication` was a 0–1 fraction;
  - the STAR log was parsed by line number;
  - `File_mb` rounding.
- **New:** insert size statistics, NRF/PBC1/PBC2, peaks and FRiP, pipeline mode and version.

### QC
- **`pyTssEnrichment.py` rewritten:**
  - insertion sites in genome coordinates for minus-strand TSSs (v2 placed the far end on the wrong side);
  - SE reads handled;
  - robust chunking.
- **TSS score:** normalised by both flanks; peak searched within ±500 bp of the TSS.
- Python plots replace `plot.TSS.enrichment.R` (optigrab no longer needed) and `frag_distribution.R`.

### Peaks
- **MACS3 inside the pipeline** (no control): narrow/broad from `peak_mode`, `auto` from `target`.
- **New `NGS.peaks.py` + `peaks_pipeline_interface.yaml`:** peaks against the `control` sample, run after the main run.

### Benchmark and validation
- `benchmark/`: ART simulation (2×100, 2×50) and an adapter read-through set; `summarize.py` writes the summary tables.
  Results and conclusions: `analysis/NGS.analysis.v3.benchmark/README.md`.
- `tests/`: T1 (200k-fragment smoke tests, v2 vs v3) and T2 (full-depth runs; Edenhofer ATAC vs `02.retrim.cutadapt`).

### Framework
- **Genome resources** in `genomes/<genome>.yaml` instead of per-project path lists. New T2T assemblies: hs1 (T2T-CHM13v2.0) and mhaESC (mhaESC v1.1 + mT2T-Y v1.1, release 2026-01-29, with chrY genes and RepeatMasker; the chrY PAR is detected and hard-masked, the telomere-only chrY end is not). Build scripts in `resources/genomes/`.
- **Full pipestat schema** for every sample. The v2 protocol-filtered schema made pipestat reject keys reported in repeats mode.
- **Project schema fixed:**
  - `projects:` → `project:`, which pipestat requires. Hidden in v2 because the collator used the cached sample schema.
  - The project pipeline has its own name, `NGS.analysis.project`.
- **Version guard:** refuses output folders with results from another major version.
- **`check_project.py` preflight** checks sample-name prefix collisions, files, genome resources and controls.
- **Summarizer:**
  - genome per sample (mixed-genome projects work);
  - `<project>_files.tsv` and `<project>_BAM_files.tsv` restored, with correct absolute paths;
  - file paths kept in the stats summary;
  - configurable IGV path mapping (`igv_path_map`).
- **`generate_report.py`:** units in column headers, newest flag wins, project records read from `NGS.analysis.project`.
- **Environment:** `ngs.v3` = exact copy of `ngs.v2` (conda + pip + R packages) plus MACS3 3.0.5 (with its anndata/zarr/h5py dependencies; zlib 1.3.1 → 1.3.2) and gffread. 3.0.5 rather than 3.0.4: fixes an `hmmratac` crash in `refine_peaks` (MACS #735) and two `--call-summits` bugs (#747, #748); the pipeline's `callpeak` without `--call-summits` is not affected by these. Interfaces use the environment's python and `{looper.piface_dir}` (no hard-coded installation paths).
- Removed `NGS.analysis_output_schema.yaml` (stale), `plot.TSS.enrichment.R` and `frag_distribution.R`.
