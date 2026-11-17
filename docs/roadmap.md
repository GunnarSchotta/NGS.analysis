# NGS.analysis v3: development roadmap (as of 2026-11-17)

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

## Next: small, independent changes (each needs a decision first)
| ID | Item | Source | Kind | Size |
|---|---|---|---|---|
| D1 | Write up the peaks benchmark conclusions (numbers in `results/RESULTS.md` since 2026-10-07) | `analysis/NGS.analysis.v3.peaks.benchmark/README.md` | doc | S |
| D2 | Peak defaults for histone marks: H3K27ac narrow vs broad, broad-cutoff 0.1 vs 0.2, q 0.01 vs 1e-4, `--call-summits` (v3 vs 01.Angela Jaccard 0.37 for H3K27ac). Combine with MACS3 3.1.0 (faster; changes control-based results) | peaks benchmark `results/RESULTS.md` | D | M |
| D3 | ATAC peak mode: BAMPE vs Tn5-shifted cut sites vs HMMRATAC (benchmark arms exist) | peaks benchmark README | D | S |
| D4 | Project-level consensus peak set (`make.common.peakset` logic) so peak numbers are comparable across samples | `analysis/Vincendeau/ATAC.HERVKC4/README.md` | F | M |
| B2 | RNA genes mode uses STAR's default multimapper cap (10): reads with more loci are dropped, which suppressed HERV signal in the HERVKC4 RNA-seq | `analysis/Vincendeau/HERVKC4/locus_quantification/README.md` | D | S |
| B3 | MAPQ ≥ 30 removes transgene/TRE and young-repeat signal from `.filt.bw`; optional all-MAPQ bigwig | `analysis/LiLi/playground/05.*`, `06.*` READMEs | F | S |
| B1 | Trimmomatic `MINLEN:30` drops < 30 bp fragments that cutadapt `-m 20` keeps (about 0.2% of ATAC; matters for footprinting) | validation README | D | S |
| B4 | bowtie2 end-to-end: no soft-clips, no junction evidence at transgene borders | LiLi 06 README | D | S (low priority) |
| G4 | No blacklist for hs1 / mhaESC (chrM and artefact peaks survive there); IAP gag BEDs and GenMap index only for mm10 | `genomes/README.md` | F | M each |
| E2 | Repeats-mode CR/CT nuc/subnuc bigwigs are not registered as report objects | `NGS.analysis.py` CR/CT block | B | S |
| C5 | Mito % per fragment is reported, but not the chrM duplication rate (library-quality indicator for ATAC); optional | – | F | S |
| S1 | samtools 1.22.1 → 1.25 (1.22.4 backport not on bioconda; 1.24/1.25 change `fixmate` on supplementary records, `markdup`, `stats`) | `env/build_ngs.v3.sh` | D | S + T1 |

## Repeats-mode redesign (XL, decided together; repeats mode stays frozen until then)
| ID | Item | Kind |
|---|---|---|
| A1 | Multimapper allocation: random-1 (current), 1/n, EM (Telescope, TEtranscripts, SQuIRE, SalmonTE), signal-aware (Allo, SmartMap), T3E for family-level chromatin | D |
| A2 | "Unique" rule for element-level counts: STAR MAPQ 255 vs bowtie2 MAPQ ≥ 30 (real data: family results unchanged, IAPEz element calls 680 → 395) | D |
| A3 | Deduplication: duplicates only marked (issue 1), unreliable for multimappers (issue 2); sequence-level dedup, UMIs, or report both | D/B |
| A4 | featureCounts: `-p` without `--countReadPairs` (3), unstranded `-s 0` (4), no `-O` (5) | B/D |
| A5 | STAR: `--alignMatesGapMax 350` (6), RSEM on `--outSAMmultNmax 1` BAM (7), `EndToEnd` (8), `--outFilterMismatchNmax 3` | B/D |
| A6 | CR/CT TLEN 0 to subnuc (10), IAP normalisation pairs vs mates and chrM in the denominator (11), NUMTs under random placement (12) | B |
| A7 | Locus-level RNA quantification with count filters and TSS evidence | D/F |
| A8 | Input normalisation for family-level ChIP counts | F |
| A9/A10 | Per-copy mappability (GenMap, k = 36–150); gain of PE and 2×100 for young L1/IAP/MERVL | A |
| G1–G3 | Mouse reference (mhaESC, Francis 2025, mm39), RepBase ↔ Dfam names, strain background (129, polymorphic TEs) | D |
| H1–H6 | Benchmarks: cross-reference simulation, SNP/indel injection, family-level correctness of random placement, alternative strategies on the ART set, 01.Angela RNA/H3K9me3 comparison, Nanopore ground truth | A |

The redesign is also when the code gets split into modules (see below): the new repeats strategy is written as its own
module next to the frozen v2 one, selectable per project, and compared on T1/T2 and real data before it replaces it.

## Decision: no per-protocol sub-pipelines (2026-11-17)
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
- **Deferred:** samtools 1.25 (S1); MACS3 3.1.0 (with D2); bowtie2 2.5.5 (only `-k` determinism and AVX dispatch;
  relevant for A1 if bowtie2 `-k` is used).
- **Do not upgrade without re-validation:** deepTools 4.0 (Rust rewrite; `--ignoreDuplicates`/`--exactScaling`
  removed, bigwig values not identical), Trimmomatic 0.41 (requires Java 25; 0.40 adds poly-G clipping).
- STAR 2.7.11b and subread 2.1.1 are still the latest releases.
