# Repeat analysis: open questions for a later deep-dive

In v3, **repeats mode is unchanged from v2.** It still uses the Teissandier et al. 2019 strategy:
- STAR, reporting one random locus per multimapper;
- a `-q 255` unique BAM, with duplicates marked but not removed;
- featureCounts at family level (all reads) and at element level (unique reads);
- IAP gag coverage.

v3 changes only two things here: the upstream trimming fix and correctly labelled statistics. The goal of a
future revision is the most reliable assignment of reads to **individual repeat copies**. This document
collects what has to be decided first (agreed 2026-10-03).

## Known issues in the current implementation (not fixed on purpose)

1. **Duplicates are not removed** in `.dedup.unique.bam`, its bigwig, the element-level counts or the CR/CT splits. The family-level counts include duplicates as well.
2. **Duplicate detection is unreliable for multimappers.** The random locus choice depends on the read name, so PCR copies of one molecule can land at different copies.
3. **featureCounts counts mates separately** (`-p` without `--countReadPairs` in subread ≥ 2.0.2), so PE counts are about 2× fragments.
4. **Repeat counting is unstranded** (`-s 0`), also for stranded RNA-seq, so sense and antisense transcription are merged.
5. **No `-O`**, so reads overlapping nested or overlapping RepeatMasker entries are `Unassigned_Ambiguity`.
6. **The STAR mate gap is capped at 350 bp** (`--alignMatesGapMax 350`), which cuts nucleosomal ATAC and CUT&RUN fragments.
7. **RSEM** runs on the transcriptome BAM produced with `--outSAMmultNmax 1`.
8. **STAR `--alignEndsType EndToEnd`** gives no soft-clipping, so residual adapter bases (SE, < 12 nt) or protruding mates (`--alignEndsProtrude 0`) prevent alignment.
9. **Short fragments recovered by the v3 trimming** (30–60 bp) add many multimappers. Their effect on family counts must be measured (validation comparison R2).
10. **CR/CT split:** records with TLEN 0 go to the sub-nucleosomal BAM.
11. **IAP normalisation** uses STAR unique + multi + too-many-loci *pairs*, while `bedtools coverage` counts every mate.

## Questions to settle

- **Theoretical mappability per copy and family** for the read and fragment lengths actually used (genmap, k = 36 to 150). Which families are addressable at element level at all?
- **Reference choice:**
  - mm10 vs the T2T mouse (mhaESC + T2T-Y) vs mm39;
  - hg38 vs T2T-CHM13 (hs1).
  - Young families sit in the regions T2T resolves; does element-level assignment improve?
- **Strain and individual polymorphisms:**
  - SNPs and indels relative to the reference (C57BL/6 reference vs 129 / mixed ES lines);
  - polymorphic TE insertions absent from the reference (reads forced onto paralogues);
  - mismatch allowance (`--outFilterMismatchNmax 3`).
  - Consider strain-specific references or variant-aware simulation.
- **Allocation of multimappers:**
  - random (current), fractional (1/n), EM (TEtranscripts, Telescope, SQuIRE, SalmonTE);
  - signal-aware allocation for chromatin data: Allo (Morrissey et al., Genome Res 2024; CNN on the profile of uniquely mapped reads; requires bowtie2 `-k 25`; trained on TF/ATAC/DNase, not on broad marks); SmartMap (PLoS Comput Biol 2021; Bayesian, PE only);
  - family-level enrichment against input with multimapper weighting: T3E (Mobile DNA 2022).
- **Locus-level RNA quantification:** Telescope / TElocal / SQuIRE / SalmonTE. Schwarz et al., Brief Bioinform 2022, rank SalmonTE* and Telescope highest for locus-level. Gene-overlapping TE loci need to be handled.
- **Deduplication strategy for multimappers:** sequence-level deduplication before alignment, UMIs, or reporting both with and without duplicates.
- **Read and fragment length:** PE 2×100 vs 2×50/60; the gain of PE for young L1/IAP/MERVL (Teissandier: +10–30% mapping).

## Evidence to collect

- The ART simulation benchmark (`benchmark/`): per-family mapping % and true-positive rate for the current settings vs alternatives, on mm10/hg38 vs T2T, with injected SNPs and indels.
- Real data: Setdb1 KO RNA-seq and H3K9me3 ChIP (01.Angela), comparing family- and element-level results between strategies.
