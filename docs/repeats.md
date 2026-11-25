# Repeat outputs of chromatin samples (3.3.0+)

ATAC, ChIP, CUT&RUN and CUT&Tag samples are aligned once with bowtie2, and the gene outputs and the repeat outputs both
come from that alignment. There is no separate repeats mode for chromatin any more (RNA keeps `genes`/`repeats`).
The decision and its evidence are in `docs/roadmap.md` (A0) and `analysis/NGS.analysis.v3.benchmark/README.md`.

## Two fragment sets, two tiers

| Set | Definition | Used for |
|---|---|---|
| filtered (`<s>.filt.bam`) | proper pairs, duplicates removed, canonical chromosomes, no chrM, **MAPQ ≥ 10**, blacklist removed | peaks, FRiP, `<s>.filt.bw` (also the locus-plot bigwig), QC, CR/CT split |
| all (`<s>.all.bam`, temporary) | proper pairs, duplicates removed, canonical chromosomes, no chrM, **any MAPQ**, blacklist kept | `<s>.all.bw`, repeat counts, IAP coverage |

bowtie2 places a multimapper at one of its best loci, chosen at random. Fragments of the all set are therefore at the
right copy of their family for most families, but not necessarily at the right copy (locus).

Repeat counts are reported in two **tiers** of the all set:
- **all**: every fragment (MAPQ ≥ 0). Use for **family-level** questions.
- **locus**: MAPQ ≥ 10. Use for **element-level** (individual copy) questions. bowtie2 gives both mates of a
  concordant pair the same MAPQ, so the tier is a property of the fragment.

The locus tier is the same MAPQ cutoff as the filtered set, so `<s>.filt.bw` (CPM, MAPQ ≥ 10) is the bigwig for locus
plots; `<s>.all.bw` shows every fragment. (3.3.0 development had a separate `.locus.bw` while `.filt.bam` was MAPQ 30;
it was dropped when the filtered set moved to MAPQ 10.)

Why MAPQ ≥ 10 (ART simulation of mm10 chr1, 2×100, `analysis/NGS.analysis.v3.benchmark/`):
- STAR's "unique" (MAPQ 255) is the same as bowtie2 MAPQ ≥ 2 in yield and accuracy. It cannot be made stricter.
- With a sample that differs from the reference like a 129 strain from C57BL/6 (0.2 % SNPs + indels, non-reference young TE
  copies), young L1 fragments are at the correct copy in 88 % (MAPQ ≥ 2), 98.5 % (≥ 10) and 98.6 % (≥ 30) of cases.
  MAPQ ≥ 30 adds no precision but keeps up to 75 % fewer young-element fragments.

Why the filtered set (peaks, FRiP) uses MAPQ ≥ 10 too (real data, `analysis/NGS.analysis.v3.peaks.benchmark/mapq10/`,
4 ATAC and 6 histone libraries, 2 replicates per mark): MAPQ 10 adds 3 % of fragments (ATAC, H3K27ac) to 4–9 %
(H3K9me3, H4K20me3). ATAC and H3K27ac peaks are unchanged (Jaccard 0.96–0.97, FRiP ≤ 0.2 points). H3K9me3 and
H4K20me3 domains shift (Jaccard 0.81–0.89), less than the replicate-to-replicate difference at MAPQ 30, and replicate
concordance does not drop. Fragments in young elements are 5–12 % of the extra fragments (0.1–1.1 % of all fragments),
and young elements are flagged. Numbers that change relative to 3.2.0: fragment counts, FRiP, peak numbers, NRF/PBC,
`Unique_*`, TSS score.

## Files (`repeats_<genome>/`)

**`<s>.repeats.family.tsv`**, one row per family of the genome's repeat SAF:

| Column | Meaning |
|---|---|
| `frag_all`, `frag_locus` | fragments overlapping ≥ 1 element of the family (all / locus tier). A fragment touching two families counts for both (featureCounts `-O`). |
| `locus_pct` | `frag_locus / frag_all` × 100 |
| `n_elements` | elements of the family in the annotation |
| `n_el_locus` | elements with ≥ 1 locus-tier fragment |
| `young_locus_pct` | % of the family's element-level locus assignments on young elements |
| `top10_locus_pct` | % of the family's element-level locus assignments in its 10 strongest elements. High = the family signal is carried by a few copies. |

**`<s>.repeats.elements.tsv.gz`**, one row per element with ≥ 1 fragment:

| Column | Meaning |
|---|---|
| `element_id`, `family`, `class`, `repfamily`, `chr`, `start`, `end`, `strand` | from the element SAF and RepeatMasker |
| `div_pct` | divergence of the element from its family consensus (RepeatMasker), % |
| `frag_all`, `frag_locus` | fragments (all / locus tier) |
| `locus_frac` | `frag_locus / frag_all`. Low = the element's signal comes mostly from multimappers placed at random. |
| `flags` | comma list, `.` if none: `young`, `overlap`, `blacklist` (below) |

Counting details:
- fragments (PE: read pairs, `--countReadPairs`);
- duplicates removed (Picard). Duplicates of multimappers can land on different copies and then escape detection, so
  the all tier is slightly inflated for young, high-copy families (roadmap A3).

## Uncertainty markers

The pipeline cannot know the genome of the sample (strain, polymorphic insertions, changes in culture). It marks where
that matters instead:

| Marker | Where | What it means |
|---|---|---|
| `locus_frac` | element, per sample | Share of the element's fragments that are placed with MAPQ ≥ 10. |
| `young` flag | element, static | Interspersed repeat < 5 % diverged from its consensus. In the simulation, such elements received 5–15 % of their MAPQ ≥ 10 fragments from elsewhere when the sample differed from the reference like a 129 strain, mostly from **non-reference copies** (a copy present in the sample but not in the reference sends its reads to its closest reference relative). No MAPQ cutoff removes this. Without strain differences the same elements receive < 0.4 % wrong fragments. Above 5 % divergence: 0.6 % (5–8 %), < 0.2 % (≥ 8 %). |
| `overlap` flag | element, static | Another annotation entry covers ≥ 10 % of the element (nested or overlapping annotation). Fragments there are counted for both, so an effect at such an element may belong to the other entry. mhaESC's Dfam annotation has many overlapping hits (46 % of elements overlap another entry vs 4 % in mm10). |
| `blacklist` flag | element, static | The element overlaps the genome blacklist (artefact-prone signal). The filtered set removes these regions, the repeat counts do not. |
| `Young_locus_pct` | sample statistic | % of all element-level locus assignments that fall on young elements: how much a sample's locus-level result depends on copies where non-reference insertions cannot be excluded. |
| `Repeat_locus_pct` | sample statistic | % of repeat fragments in the locus tier. |

Suggested use:
- **Family level:** `frag_all`; check `top10_locus_pct` and `n_el_locus` to see whether an effect is carried by many
  copies or by a few, and look at those few (overlap flag, genomic context).
- **Element level:** `frag_locus`; treat `young` elements as candidates that need independent support (strain
  information, flanking unique signal, long reads), not as established copies.

Calibration of the young threshold: `analysis/NGS.analysis.v3.benchmark/README.md`, section "Element calibration".

## Static annotation (per genome)

`genomes/<g>.yaml: repeats_annotation` points to `repeat_elements.<g>.tsv.gz`, built once per genome by
`resources/genomes/build_repeat_annotation.py`. It is in the same row order as the element SAF; the pipeline stops
with an error if they differ (rebuild after changing the SAF). Divergence sources: UCSC `rmsk.txt` `milliDiv` (mm10,
hg38), RepeatMasker `.out` (hs1) and RepeatMasker GFF score (mhaESC). Elements missing from the divergence source
(mm10 3.5 %, hg38 2 %: the SAF and the rmsk table come from different RepeatMasker releases) have `div_pct` NA and no
young flag.

## IAP coverage (mm10)

Per-base depth of the all set over the IAP gag reference (every mate counts, as in v2), per million `All_fragments`.
v2 normalised by STAR unique + multi + too-many-loci pairs including chrM, so absolute values differ from v2.
