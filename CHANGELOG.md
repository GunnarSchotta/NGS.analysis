# Changelog

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

### Framework
- **Genome resources** in `genomes/<genome>.yaml` instead of per-project path lists. New T2T assemblies: hs1 (T2T-CHM13v2.0) and mhaESC (T2T mouse with chrY, PAR masked).
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
- **Environment:** `ngs.v3` = exact copy of `ngs.v2` (conda + pip + R packages) plus MACS3 and gffread. Interfaces use the environment's python and `{looper.piface_dir}` (no hard-coded installation paths).
- Removed `NGS.analysis_output_schema.yaml` (stale), `plot.TSS.enrichment.R` and `frag_distribution.R`.
