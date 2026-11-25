# NGS.analysis v3: development roadmap (as of 2026-11-25)

One line per item: what, where it is documented, kind (**D** decision, **B** bug, **F** feature, **A** analysis),
size (S < 1 day, M a few days + T1 rerun, L 1–3 weeks or a benchmark, XL redesign). Item IDs are stable; refer to
them in commits and READMEs. Repeats-mode details and evidence: `repeat_strategy_open_questions.md`.

## Done in 3.1.0
- **chrM split** for chromatin samples (ATAC/CHIP/CR/CT, both modes): fragments with the read or mate on chrM are
  removed straight after alignment, before `<s>.bam`. Duplication, NRF/PBC, peaks, FRiP, bigwigs and counts are now
  nuclear only; new keys `Nuclear_fragments`, `Mito_filter`. Prompted by the Vincendeau ATAC run (reported duplication
  66% vs about 43% nuclear; 42–64% of in-peak fragments on chrM peaks). Aligned against the full genome (not
  PEPATAC-style prealignment), so MAPQ arbitrates chrM vs nuclear NUMTs.
- Picard 3.5.0 (MarkDuplicates optical-duplicate overflow fix). Version guard now on major.minor.
- Stale docs fixed (README QC row, peaks `--keep-dup`, per-sample `genome_config`, issue 9 measured, T2T evidence).

## Done in 3.2.0
- **D1** (2026-11-24): peaks benchmark conclusions written, `analysis/NGS.analysis.v3.peaks.benchmark/README.md`.
  Main finding: MACS3 differs from MACS2 only through `-g mm` (1.87e9 in MACS2, 2.65e9 in MACS3); MACS3 with
  `-g 1.87e9` reproduces MACS2 exactly.
- **D2** (decided 2026-11-24): (a) keep MACS3's larger genome size, pinned as a number (mm10/mhaESC 2,652,783,500,
  hg38/hs1 2,913,022,398); (b) 01.Angela histone settings: broad `-q 1e-4 --broad-cutoff 0.2`, H3K27ac broad;
  (c) per-sample no-control peaks documented as QC-only for histone marks. MACS3 3.1.0 (identical peaks, 5–14× faster).
- **D3** (decided 2026-11-24): ATAC `--call-summits`; `Peaks_n` counts regions, new `Peak_summits_n`; FRiP per region
  (with one feature per summit record FRiP would drop 32.8 → 25.8 %). BAMPE stays; Tn5 mode not added; HMMRATAC dropped.
- Validation: `analysis/NGS.analysis.v3.validation/T4_macs3.1/`, T1 re-run (R1 still identical), T1 `NGS.peaks` smoke test.

## Done in 3.3.0
- **A0 implemented** (2026-11-25): one bowtie2 mode for chromatin with gene and repeat outputs (`docs/repeats.md`);
  RNA keeps `genes`/`repeats` (repeats frozen). Repeat counts in two tiers (all fragments; locus MAPQ ≥ 10), per-element
  uncertainty markers (`locus_frac`; `young` = < 5 % divergence, calibrated on the SNP/non-reference TE simulation;
  `overlap`; `blacklist`), `Young_locus_pct`. Fixes issues 1, 3, 5, 6, 11 for chromatin (12: NUMT-placed reads
  have low MAPQ and stay out of the locus tier). Gene outputs unchanged (T1 identical to 3.2.0).
- Validation: T1 (`analysis/NGS.analysis.v3.validation/T1/compare/`, R4), T5 Setdb1 KO ATAC at full depth
  (`analysis/NGS.analysis.v3.validation/T5_single_mode/`), MAPQ 10 vs 30 for the filtered BAM
  (`analysis/NGS.analysis.v3.peaks.benchmark/mapq10/`), element calibration (`analysis/NGS.analysis.v3.benchmark/element_calib/`).

## Next: small, independent changes (each needs a decision first)
| ID | Item | Source | Kind | Size |
|---|---|---|---|---|
| D5 | Default control for histone peaks without a matched input. Tested 2026-11-25 (`analysis/NGS.analysis.v3.peaks.benchmark/control_test/`): a fixed real input of another cell type ≈ matched input (Jaccard 0.71–0.76); simulated input recovers 86–96 % of peaks but adds 10–50 % weaker, more promoter-proximal regions; `--nolambda`/no control unsuitable. To decide: per-genome default input in `genomes/<g>.yaml`, simulated input as fallback | control test README | D/F | M |
| D4 | Project-level consensus peak set (`make.common.peakset` logic) so peak numbers are comparable across samples | `analysis/Vincendeau/ATAC.HERVKC4/README.md` | F | M |
| B2 | RNA genes mode uses STAR's default multimapper cap (10): reads with more loci are dropped, which suppressed HERV signal in the HERVKC4 RNA-seq | `analysis/Vincendeau/HERVKC4/locus_quantification/README.md` | D | S |
| B3 | ~~MAPQ ≥ 30 removes transgene/TRE and young-repeat signal from `.filt.bw`; optional all-MAPQ bigwig~~ done in 3.3.0 (`.filt.bw` now MAPQ ≥ 10, plus `.all.bw`) | `analysis/LiLi/playground/05.*`, `06.*` READMEs | F | S |
| B1 | Trimmomatic `MINLEN:30` drops < 30 bp fragments that cutadapt `-m 20` keeps (about 0.2% of ATAC; matters for footprinting) | validation README | D | S |
| B4 | bowtie2 end-to-end: no soft-clips, no junction evidence at transgene borders | LiLi 06 README | D | S (low priority) |
| G4 | No blacklist for hs1 / mhaESC (chrM and artefact peaks survive there); IAP gag BEDs and GenMap index only for mm10 | `genomes/README.md` | F | M each |
| E2 | Repeats-mode CR/CT nuc/subnuc bigwigs are not registered as report objects | `NGS.analysis.py` CR/CT block | B | S |
| C5 | Mito % per fragment is reported, but not the chrM duplication rate (library-quality indicator for ATAC); optional | – | F | S |
| S1 | samtools 1.22.1 → 1.25 (1.22.4 backport not on bioconda; 1.24/1.25 change `fixmate` on supplementary records, `markdup`, `stats`) | `env/build_ngs.v3.sh` | D | S + T1 |

## Repeats-mode redesign (XL; chromatin part started with A0 in 3.3.0, RNA repeats mode stays frozen)

**Decided 2026-11-25 (A0, mapping):** the focus is locus-level analysis.
- Chromatin (ChIP/ATAC/CUT&RUN/CUT&Tag) is aligned with **bowtie2** (v3 genes-mode settings `--very-sensitive -X 2000
  --dovetail`, end-to-end) in **one mapping mode**; STAR is dropped for chromatin. Evidence: Teissandier 2019 found no
  accuracy difference between STAR and bowtie2 (STAR chosen for speed); on our ART set STAR "unique" (MAPQ 255) =
  bowtie2 MAPQ ≥ 2 in yield and accuracy, while bowtie2's graded MAPQ gives stricter locus tiers (≥ 10: ≥ 99.6 %
  correct for young L1/ERV) that STAR cannot; STAR took about 2× bowtie2's time on real 01.Angela data
  (`analysis/NGS.analysis.v3.benchmark/README.md`, MAPQ curves; `repeat_strategy_open_questions.md`).
- RNA stays on **STAR** (splicing).
- Locus-level MAPQ tier: the SNP/indel and non-reference TE simulation (H2, done 2026-11-25) supports **MAPQ ≥ 10**:
  at 129-like divergence it keeps young-L1 precision at 98.5 % (2×100; MAPQ ≥ 2: 88 %), and ≥ 30 adds no precision
  but loses up to 75 % of young-element fragments. Non-reference TE insertions still place 13–21 % of their reads on
  a reference copy at MAPQ ≥ 10; only strain information (G3) or long reads (H6) address that. Final choice of the
  tier: with the module implementation. A family-level sub-analysis (how many copies drive a family effect; overlapping/nested
  annotation) is done after mapping (A11).

| ID | Item | Kind |
|---|---|---|
| A11 | Family-level sub-analysis after mapping: is a family effect carried by the majority of copies or by a few loci, and do those loci overlap other (nested, non-ERV) annotation? Per-sample basis in 3.3.0 (family `n_el_locus`, `top10_locus_pct`, element `overlap` flag); the across-sample analysis (which copies drive a differential effect) is still open | F |
| A12 | Non-reference insertions: flag families/elements with known strain-polymorphic copies (catalogues of polymorphic TE insertions, e.g. Mouse Genomes Project SVs) or detect them from the data (discordant pairs at element ends; input/WGS of the line). 3.3.0 only marks the risk (young flag) | D/F |
| A1 | Multimapper allocation: random-1 (current), 1/n, EM (Telescope, TEtranscripts, SQuIRE, SalmonTE), signal-aware (Allo, SmartMap), T3E for family-level chromatin | D |
| A2 | "Unique" rule for element-level counts: STAR MAPQ 255 vs bowtie2 MAPQ ≥ 30 (real data: family results unchanged, IAPEz element calls 680 → 395) | D |
| A3 | Deduplication: chromatin 3.3.0 removes duplicates (issue 1); still unreliable for multimappers (issue 2: duplicates can land on different copies, slightly inflating the all tier of young families); sequence-level dedup, UMIs, or report both | D/B |
| A4 | featureCounts: `-p` without `--countReadPairs` (3), unstranded `-s 0` (4), no `-O` (5). Chromatin: 3 and 5 fixed in 3.3.0; RNA repeats mode unchanged | B/D |
| A5 | STAR: `--alignMatesGapMax 350` (6), RSEM on `--outSAMmultNmax 1` BAM (7), `EndToEnd` (8), `--outFilterMismatchNmax 3` | B/D |
| A6 | CR/CT TLEN 0 to subnuc (10), IAP normalisation pairs vs mates and chrM in the denominator (11), NUMTs under random placement (12) | B |
| A7 | Locus-level RNA quantification with count filters and TSS evidence | D/F |
| A8 | Input normalisation for family-level ChIP counts | F |
| A9/A10 | Per-copy mappability (GenMap, k = 36–150); gain of PE and 2×100 for young L1/IAP/MERVL | A |
| G1–G3 | Mouse reference (mhaESC, Francis 2025, mm39), RepBase ↔ Dfam names, strain background (129, polymorphic TEs) | D |
| H1–H6 | Benchmarks: cross-reference simulation, SNP/indel injection (H2, done 2026-11-25), family-level correctness of random placement, alternative strategies on the ART set, 01.Angela RNA/H3K9me3 comparison, Nanopore ground truth | A |

The redesign is also when the code gets split into modules (see below): the new repeats strategy is written as its own
module next to the frozen v2 one, selectable per project, and compared on T1/T2 and real data before it replaces it.

## Decision: no per-protocol sub-pipelines (2026-11-24)
Considered: separate looper pipelines per protocol (ATAC, ChIP, CR/CT, RNA). Not done, because:
- **Protocol-specific code is small.** ATAC adds the TSS score (about 17 lines), CR/CT the nuc/subnuc split (about 17),
  RNA RSEM and its bigwig (about 30); ChIP has nothing of its own. About 45% of `NGS.analysis.py` (setup, trimming,
  statistics) is shared by everything.
- **The real forks are aligner and mode** (bowtie2-genes, STAR-repeats, RNA), not protocol. A protocol split would
  duplicate the alignment/filter code across ATAC/ChIP/CR/CT.
- **Framework cost.** `generate_report.py` (`infer_pipeline_name`, flag lookup) and `NGS.summarizer.R` assume one sample
  pipeline name and one schema; separate pipelines need several schemas, flag sets and report paths. Looper 2.1.1 could
  pick an interface per sample (`sample_modifiers.imply` → `pipeline_interfaces`), but this is untested here.
- **Planned instead:** with the repeats redesign, move the code into modules (common: input/trimming/stats;
  chromatin-genes; chromatin-repeats; RNA) behind the same entry point, schema and pipeline name.

## Package status (checked 2026-11-17)
- **PEP stack is current:** looper 2.1.1, pipestat 0.13.1, piper (pypiper) 0.15.1, yacman 1.0.0, ubiquerg 0.9.3 are the
  latest stable releases. peppy 0.50 is pre-release only (it absorbs eido and pephubclient; expect conflicts) – stay
  on 0.40.8. eido 0.2.6 / pephubclient 0.6.0 bring nothing we use.
- **Updated in 3.1.0:** Picard 3.4.0 → 3.5.0.
- **Deferred:** samtools 1.25 (S1); MACS3 3.1.0 (with D2; released 2026-10-21, speed + `--threads`, fixes not relevant to our `--keep-dup` settings); bowtie2 2.5.5 (only `-k` determinism and AVX dispatch;
  relevant for A1 if bowtie2 `-k` is used).
- **Do not upgrade without re-validation:** deepTools 4.0 (Rust rewrite; `--ignoreDuplicates`/`--exactScaling`
  removed, bigwig values not identical), Trimmomatic 0.41 (requires Java 25; 0.40 adds poly-G clipping).
- STAR 2.7.11b and subread 2.1.1 are still the latest releases.
