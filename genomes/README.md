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
| `macs_gsize` | MACS3 `-g` value (shorthand or number) |
| `iap_plus_bed`, `iap_minus_bed` | IAP gag coverage BEDs (repeats mode; mm10 only), relative to the pipeline folder |

Effective genome size: MACS shorthands (`mm` = 1.87e9, `hs` = 2.7e9) for mm10/hg38, as used in earlier
analyses; deepTools 100-bp mappable sizes for the T2T assemblies.
