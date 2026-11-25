# Repeat analysis: open questions for a later deep-dive

**Since 3.3.0, chromatin samples no longer use repeats mode:** they are aligned once with bowtie2 and get repeat
counts in two MAPQ tiers with per-element uncertainty flags (`docs/repeats.md`). Issues 1, 3, 5, 6 and 11 below are
fixed for chromatin there; 2 (duplicates of multimappers) and the allocation question remain. The text below
describes **RNA repeats mode**, which is unchanged from v2 and still uses the Teissandier et al. 2019 strategy:
- STAR, reporting one random locus per multimapper;
- a `-q 255` unique BAM, with duplicates marked but not removed;
- featureCounts at family level (all reads) and at element level (unique reads);
- IAP gag coverage.

v3.0 changed only two things here: the upstream trimming fix and correctly labelled statistics. v3.1 adds the chrM
split for chromatin samples (fragments touching chrM are removed before `<s>.bam`; `--keep-mito` restores v2). The goal of a
future revision is the most reliable assignment of reads to **individual repeat copies**. This document
collects what has to be decided first (agreed 2026-10-03).

Literature for each question (PubMed search 2026-10-03, with DOIs and what each paper means for us): `repeat_literature.md`.

## Known issues in the current implementation (not fixed on purpose)

1. **Duplicates are not removed** in `.dedup.unique.bam`, its bigwig, the element-level counts or the CR/CT splits. The family-level counts include duplicates as well.
2. **Duplicate detection is unreliable for multimappers.** The random locus choice depends on the read name, so PCR copies of one molecule can land at different copies.
3. **featureCounts counts mates separately** (`-p` without `--countReadPairs` in subread ≥ 2.0.2), so PE counts are about 2× fragments.
4. **Repeat counting is unstranded** (`-s 0`), also for stranded RNA-seq, so sense and antisense transcription are merged.
5. **No `-O`**, so reads overlapping nested or overlapping RepeatMasker entries are `Unassigned_Ambiguity`.
6. **The STAR mate gap is capped at 350 bp** (`--alignMatesGapMax 350`), which cuts nucleosomal ATAC and CUT&RUN fragments.
7. **RSEM** runs on the transcriptome BAM produced with `--outSAMmultNmax 1`.
8. **STAR `--alignEndsType EndToEnd`** gives no soft-clipping, so residual adapter bases (SE, < 12 nt) or protruding mates (`--alignEndsProtrude 0`) prevent alignment.
9. **Short fragments recovered by the v3 trimming** (30–60 bp) add multimappers. *Measured 2026-10-07, see evidence below: acceptable, kept.*
10. **CR/CT split:** records with TLEN 0 go to the sub-nucleosomal BAM.
11. **IAP normalisation** uses STAR unique + multi + too-many-loci *pairs*, while `bedtools coverage` counts every mate. Since 3.1 the coverage is nuclear only (chromatin), but the denominator from `Log.final.out` still includes chrM.
12. **NUMTs under random-1 placement.** A read that maps equally well to chrM and a nuclear NUMT is placed on either copy at random. The 3.1 chrM split removes the chrM-placed half; the NUMT-placed half stays in `<s>.bam` (MAPQ < 255, so not in the unique BAM, but in family-level counts and IAP coverage). Any multimapper strategy (A1) should treat chrM as a competing locus.

## Questions to settle

- **Theoretical mappability per copy and family** for the read and fragment lengths actually used (genmap, k = 36 to 150). Which families are addressable at element level at all?
- **Reference choice:**
  - mouse: mm10 vs a T2T mouse (mhaESC + T2T-Y, built; Francis et al. 2025 C57BL/6J + CAST/EiJ, not built) vs mm39 (conservative update, not built);
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

## Evidence collected so far

The ART benchmark (`analysis/NGS.analysis.v3.benchmark/README.md`) simulates reads from the reference itself, with no SNPs, so its numbers are an upper bound. In the mm10/hg38 runs of 2026-10-03:
- **STAR "unique" (MAPQ 255) is not error-free for the youngest L1.**
  - L1MdT/L1MdA and L1HS: 3–5% (2×100) and 6–7% (2×50) of the unique fragments are at the wrong locus.
  - IAPEz: 1.4–3%. MMERVK10C, ETnERV, SVA, HERVK, LTR5_Hs: ≤ 0.5%.
  - Element-level counts for young L1 therefore carry a few percent misassignment even in the best case.
- **Unique fractions:**
  - STAR MAPQ 255 keeps 29–82% of young-ERV/L1 fragments (2×100);
  - bowtie2 MAPQ ≥ 30 keeps 10–47%, at ≥ 99.6% correct locus;
  - bowtie2 MAPQ ≥ 42 keeps 0.6–9%.

  bowtie2 MAPQ ≥ 30 could serve as an alternative, stricter "unique" definition for element-level counts.
- **Placement accuracy and read length:**
  - Of all fragments placed in repeats mode, including random multimapper placement, 94–98% are at the correct locus.
  - 2×50 roughly halves the unique fraction of young L1/IAP compared with 2×100.
- **Real data, uniqueness rule (01.Angela ATAC, 2×60, 2026-10-07):** STAR MAPQ 255 (`v3rep`) vs bowtie2 MAPQ ≥ 30 (`v3genes`), same trimming, duplicates excluded (`analysis/NGS.analysis.v3.peaks.benchmark/mappability/README.md`).
  - bowtie2 MAPQ ≥ 30 keeps 97.8% of fragments overall, but only 64% on IAPEz-int, 70% on L1Md_T, 76% on L1Md_A and 77–79% on RLTR10C, MMERVK10C-int and ETnERV-int. L1 and ERVK elements < 5% diverged keep 73–82%; elements > 15% diverged lose nothing.
  - The Setdb1 KO d6/d0 log2FC per family is unchanged: r = 0.996 across about 1,250 families; IAPEz-int 2.00 vs 2.05 in ES.
  - Per element it changes: IAPEz-int copies with log2FC > 1 and ≥ 20 fragments drop from 680 to 395 in ES and from 182 to 107 in XEN. The bowtie2 set is almost a subset of the STAR set.
  - Peaks only found with STAR are enriched on young L1 (11–16% vs 2–4% of shared peaks in ES).
  - So the uniqueness rule matters for element-level and young-repeat peak results, not for family-level results.
- **Real data, short fragments from the v3 trimming (issue 9; 01.Angela ATAC, 2026-10-07,** `analysis/NGS.analysis.v3.peaks.benchmark/trimming/README.md`**):**
  - Repeat classes gain ×1.17–1.28 fragments vs ×1.36 genome-wide; young ERV/L1 gain least (IAPEz-int ×1.11, L1Md_T ×1.12).
  - Setdb1 KO d6/d0 log2FC per family: r = 0.998 (ES) and 0.996 (XEN), no family moves by > 0.3. Element log2FC r = 0.98–0.99.
  - Conclusion: the extra multimappers do not change family-level results; the new trimming stays.
- **STAR vs bowtie2 is a threshold question, not an aligner question (2026-11-24).**
  - **Teissandier 2019 itself found no accuracy difference** (mouse chr1, PE 2×100; Table 1): unique STAR 95.4 % mapped / 99.81 % correct vs bowtie2 (AS > XS) 94.6 % / 99.94 %; random-1 STAR 95.93 % vs bowtie2 95.90 % correct. STAR was recommended for speed (5 vs 36 min, bowtie2 2.1.0, 4 threads); locus-level accuracy per family was not compared between aligners.
  - **MAPQ curves on our ART set** (`analysis/NGS.analysis.v3.benchmark/mapq_curve/`, mm10, `benchmark/mapq_curve.py`): bowtie2 MAPQ ≥ 2 reproduces STAR MAPQ 255 in yield and accuracy (all repeats 95.7 % kept / 99.80 % correct vs 95.6 % / 99.87 %; L1MdT 48.5 / 94.2 vs 47.4 / 95.6; IAPEz-int 29.8 / 96.1 vs 29.0 / 98.6; MERVL-int 29.4 / 97.0 vs 29.1 / 97.9; 2×100). STAR's "unique" is a loci count (255 = no second alignment within 1 score point), so it cannot be made stricter; bowtie2's MAPQ grades the score gap and gives a tunable curve: MAPQ ≥ 10 keeps L1MdT 30.9 % at 99.69 % correct, IAPEz-int 13.6 % at 99.63 %, MERVL-int 17.2 % at 100 %. STAR cannot reach ≥ 99.5 % on young L1/IAP at any setting.
  - **Speed on real data:** STAR (5,000-loci settings) took 1–3× (median about 2×) the time of bowtie2 `--very-sensitive` on the same 30 01.Angela libraries, 16 threads each; worst on H3K9me3 ES (9 h vs 4 h). STAR also needs about 30 GB RAM vs a few GB.
- **Strain divergence and non-reference TE insertions (H2, 2026-11-25;** `analysis/NGS.analysis.v3.benchmark/README.md`, `snp/`**):** mm10 chr1 with 0–0.7 % SNPs (+ indels) and 300 non-reference young-TE copies, bowtie2 v3 settings.
  - Young-L1 precision at 129-like divergence (0.2 %, 2×100): MAPQ ≥ 2 88 %, ≥ 10 98.5 %, ≥ 30 98.6 %; 2×50: 83 / 94 / 94 %. MAPQ ≥ 30 adds no precision over ≥ 10 under divergence but loses up to 75 % of young-element fragments. Recommended locus-level cutoff: **MAPQ ≥ 10**.
  - Reads from non-reference copies: 19–21 % (2×100) / 13–14 % (2×50) are placed on a reference copy with MAPQ ≥ 10, 6–9 % with ≥ 30. A MAPQ cutoff cannot remove false loci created by polymorphic insertions (G3, H6).
- **Element calibration: which elements receive wrong fragments (2026-11-25;** `analysis/NGS.analysis.v3.benchmark/element_calib/`**):** the SNP/non-reference TE simulation, scored by the RepeatMasker divergence of the element where a fragment *lands*.
  - bowtie2 gives both mates of a concordant pair the same MAPQ: R1, min and max of the mates give identical tiers.
  - At MAPQ ≥ 10 and 0.2 % strain divergence (2×100), L1/LTR elements < 5 % diverged receive 5.5–11 % wrong fragments; 5–8 %: 0.6 %; 8–12 %: 0.07 %. 80–97 % of the wrong fragments come from non-reference copies. Without strain divergence the same young bins receive 0.1–0.4 %. 2×50: 7–15 % (< 5 %), 1.1 % (5–8 %).
  - → 3.3.0 `young` flag = divergence < 5 %.
- **STAR random-1 placement is not reproducible between runs** with several threads (STAR documents this for `--outMultimapperOrder Random`): in the 3.3.0 T1 run, about half of the 23,401 RNA multimappers were placed differently from the v2 run with identical commands (IAP coverage −17 %, r = 0.68 at 200k fragments). bowtie2's choice is seeded per read and reproducible.
- **T2T, human (hs1, done):** with reads simulated from the reference itself, hs1 and hg38 give the same young-TE unique fractions and TP rates once satellites are excluded. hs1 adds satellite arrays that are not uniquely mappable.
- **T2T, mouse (mhaESC, done):**
  - Excluding satellites, mhaESC is 2–3 points less uniquely mappable than mm10, at unchanged TP.
  - Same-named families (IAPEz, MMERVK10C, RLTR10C, B1/B2) behave the same.
  - Its Dfam annotation splits off the youngest L1 subfamilies; for L1MdTf_I/II and L1MdA_I, about 15% of STAR "unique" fragments are at the wrong copy (2×100).
  - The choice of reference also fixes the repeat annotation (RepBase vs Dfam names), which changes per-family results.

## Evidence still to collect

- **ART benchmark:**
  - a **cross-reference** simulation (reads from hs1 / mhaESC aligned to hg38 / mm10), which would measure the expected real T2T benefit: reads from sequence missing in the old reference being forced onto wrong paralogues;
  - family-level (not only locus-level) correctness of random placement;
  - alternatives (bowtie2 `-k`, Allo, EM).
- Real data: Setdb1 KO RNA-seq and H3K9me3 ChIP (01.Angela), comparing family- and element-level results between strategies.
