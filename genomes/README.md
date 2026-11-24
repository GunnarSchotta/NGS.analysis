# Genome resource files

One YAML per genome assembly. The pipeline loads `genomes/<genome>.yaml` (next to `NGS.analysis.py`) for the
`--genome` given by the sample (PEP `imply` on `organism`), unless `--genome-config <file>` is passed.
Explicit command-line paths (`--Bowtie2_index`, `--STAR_genome_index`, ...) override YAML values.

| key | meaning |
|---|---|
| `fasta`, `fai` | genome FASTA and its samtools index |
| `bowtie2_index` | Bowtie2 index prefix |
| `star_genome_index` | STAR index without annotation (repeats mode, non-RNA) |
| `star_rna_index` | STAR index with gene annotation (RNA) |
| `rsem_index` | RSEM reference prefix (optional) |
| `repeats_saf`, `repeats_safid` | RepeatMasker SAF: family-level (GeneID = repName) and element-level (GeneID = repName.chr:start-end) |
| `tss_bed` | TSS positions, BED with strand in column 6 |
| `blacklist` | ENCODE-style blacklist BED, or `null` |
| `mito` | mitochondrial chromosome name |
| `canonical` | regex of chromosomes kept in filtered BAMs (mito is always removed) |
| `macs_gsize` | MACS3 `-g` value: a number (since 3.2.0; a shorthand still works but depends on the MACS version) |
| `iap_plus_bed`, `iap_minus_bed` | IAP gag coverage BEDs (repeats mode; mm10 only), relative to the pipeline folder |

Effective genome size (since 3.2.0 pinned as numbers, because the shorthands changed between MACS versions:
`mm` = 1.87e9 in MACS2, 2,652,783,500 in MACS3; `hs` = 2.7e9 vs 2,913,022,398):
- mm10 = 2,652,783,500 and hg38 = 2,913,022,398, MACS3's values (non-N length of GRCm38 / GRCh38). Same peaks as
  `-g mm` / `-g hs` in MACS3 (validated in `analysis/NGS.analysis.v3.peaks.benchmark/`); compared with the MACS2 value
  they give more, weaker peaks (broad marks +66–72 %, ATAC +8–12 %), which a replicate-consensus filter removes again
  for ATAC.
- mhaESC and hs1 use the mm10 and hg38 numbers. Their FASTAs have more non-N sequence (2,867,426,149 and
  3,114,496,529; `resources/genomes/count_nonN.py`), but the extra T2T sequence is mostly satellite / segmental
  duplication that is largely unmappable with short reads; the same number keeps peak calling comparable to
  mm10/hg38 analyses.

| genome | assembly | notes |
|---|---|---|
| mm10 | GRCm38 (UCSC) | ENCODE blacklist ENCFF547MET; IAP gag BEDs for repeats mode |
| hg38 | GRCh38 (UCSC) | GRCh38 unified blacklist |
| hs1 | T2T-CHM13v2.0 analysis set (Nurk et al. 2022, Science, doi:10.1126/science.abj6987; chrY: Rhie et al. 2023, Nature, doi:10.1038/s41586-023-06457-y) | PAR-masked chrY, rCRS chrM; CAT/Liftoff GENCODE v35; UCSC hs1 RepeatMasker; no blacklist |
| mhaESC | mhaESC v1.1 + mT2T-Y v1.1 (C57BL/6 T2T mouse, release 2026-01; Liu et al. 2024, Science, doi:10.1126/science.adq8191; chrY: Li et al. 2026, Science, doi:10.1126/science.aea2249) | chrY PAR hard-masked by the build; combined gene annotation incl. chrY; RepeatMasker (Dfam names) for all chromosomes incl. chrY (chrY: many cross-species labels for Y amplicon repeats); no blacklist |

T2T resources are built with `resources/genomes/build_genome.sh` (sources: `download_t2t_sources.sh`; details in
`genomes/<g>/PREP_REPORT.txt` in the genome folder).

An independent T2T assembly of C57BL/6J (and CAST/EiJ) exists: Francis et al. 2025, Nat Genet, doi:10.1038/s41588-025-02367-z. It is not built here; whether v3 should use it instead of or in addition to mhaESC is open (see `docs/repeat_literature.md`).
