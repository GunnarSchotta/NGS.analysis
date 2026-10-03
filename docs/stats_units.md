# Statistics reported by NGS.analysis v3

**Convention:** every count is in **fragments**:
- paired-end (PE): read pairs;
- single-end (SE): reads.

Every percentage is on a 0–100 scale. The schema (`pipestat_results_schema.yaml`) states the unit of each key.

The exceptions are `Raw_reads`, `Fastq_reads` and `PF_reads`. pypiper reports these itself, and they count **reads** (for PE both mates, so 2 × pairs).

## Why v2 numbers are not comparable

| v2 key | problem in v2 | v3 replacement |
|---|---|---|
| `Trimmed_reads`, `Trim_loss_rate` | counted both mates | `Trimmed_fragments`, `Trimmed_pct` (Trimmomatic log) |
| `Mapped_reads` (STAR) | read pairs, including "too many loci" reads that are not in the BAM | `Aligned_fragments` = unique + multi |
| `Mapped_reads` (bowtie2) | all mapped records (both mates) | `Aligned_fragments` = primary R1 mapped (`-f 64 -F 2308`) |
| `Mitochondrial_reads_percentage` | records / pairs (STAR) or records / trimmed reads (bowtie2) | `Mito_pct` = mito fragments / aligned fragments, same for both aligners |
| `Percent_duplication` | 0–1 fraction | `Duplication_pct` (×100) |
| `Read_pair_duplicates` / `Read_duplicates` | different units PE vs SE | `Duplicate_fragments` |
| `Mapped_reads_filtered` | records in a BAM that still contained duplicates | `Filtered_fragments` (+ `Filtered_nondup_fragments` in repeats mode) |

## Definitions

- **Trimming** (`Trim_input_fragments`, `Trimmed_fragments`, `Trim_R1_only_fragments`, `Trim_R2_only_fragments`, `Trim_dropped_fragments`): taken from the Trimmomatic summary line. Only fragments with both mates surviving (PE) are aligned.
- **Alignment:**
  - **bowtie2:**
    - `Aligned_fragments` = primary R1 records that are mapped.
    - `Unique_fragments` = the same at MAPQ ≥ 30.
    - `Aligned_proper_pairs` = properly paired.
    - Percentages are relative to `Trimmed_fragments`.
  - **STAR:**
    - `Unique_fragments` / `Multimapped_fragments` / `Too_many_loci_fragments` are taken from `Log.final.out`, parsed by key.
    - `Aligned_fragments` = unique + multi.
    - `Unmapped_fragments` = input − aligned.
- **Mitochondria:** primary aligned fragments on the genome's `mito` chromosome (genome YAML), as a percentage of aligned fragments.
- **Duplicates:** Picard MarkDuplicates metrics, parsed by header. `Duplication_pct` = Picard `PERCENT_DUPLICATION` × 100.
- **Filtered fragments:** the analysis BAM.
  - **genes mode:** `.filt.bam`. It keeps proper pairs; removes duplicates, MAPQ < 30 reads, non-canonical chromosomes, the mito chromosome and blacklist regions; and drops orphan mates (fixmate). `Filtered_pct` = relative to trimmed fragments.
  - **repeats mode:** `.dedup.unique.bam`, as in v2: unique (MAPQ 255) alignments, with duplicates **marked but kept**. `Filtered_nondup_fragments` excludes them.
- **Insert size:** proper-pair R1 records of the analysis BAM, |TLEN|. Reported as median, mode, and % below 150, 150–299 and ≥ 300 bp. A plot is produced.
- **Library complexity** (ENCODE):
  - Computed on uniquely aligned fragments *before* duplicate removal.
  - Fragment key: PE = chromosome, left position, |TLEN|, strand; SE = chromosome, 5′ end, strand.
  - `NRF` = distinct positions / total fragments.
  - `PBC1` = positions seen once / distinct positions.
  - `PBC2` = positions seen once / positions seen twice.
  - ENCODE guide values: NRF > 0.9; PBC1 > 0.9; PBC2 > 10 for ideal libraries.
- **TSS score** (ATAC):
  - Built from the Tn5 insertion profile (±2 kb, 1-bp resolution; `pyTssEnrichment.py`), normalised by the mean of the outer 100 bp on both sides.
  - The score is the mean of ±50 bp around the maximum within ±500 bp of the TSS.
  - It is not identical to the v2/PEPATAC score, which normalised by the upstream flank only and searched for the maximum over the whole window.
- **Peaks:** MACS3, blacklist removed.
  - `Peaks_n` = number of peaks.
  - `FRiP_pct` = fragments assigned to peaks by featureCounts (fragment-level with `--countReadPairs` for PE, duplicates ignored), as a percentage of usable fragments (assigned + no feature + ambiguous + overlap-length).
