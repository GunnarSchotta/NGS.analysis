# Repeat analysis: literature review (PubMed search, 2026-10-03)

Background for `repeat_strategy_open_questions.md`. All entries were retrieved from PubMed; the DOI is given for each.

**Search:** PubMed, topic by topic:
1. allocation of multimapped reads in ChIP/ATAC;
2. locus-level TE quantification in RNA-seq and its benchmarks;
3. mappability;
4. T2T references (human and mouse);
5. mouse strain variation and polymorphic TE insertions;
6. long-read approaches;
7. duplicate handling for multimappers.

Not exhaustive. bioRxiv has no keyword search, so preprints only appear when PubMed indexes them.

## 1. Mappability and read length

| Paper | Key point for us |
|---|---|
| Teissandier et al. 2019, Mobile DNA. [10.1186/s13100-019-0192-1](https://doi.org/10.1186/s13100-019-0192-1) | Basis of our repeats mode: simulated reads (mouse and human), aligner comparison, quantification recommendations, mappability per family. |
| Sexton & Han 2019, Mobile DNA. [10.1186/s13100-019-0172-5](https://doi.org/10.1186/s13100-019-0172-5) | Human: 68%, 85% and 88% of all RepeatMasker TEs are uniquely mappable with 2×50, 2×76 and 2×100 PE reads (3 mismatches). Older TEs are addressable at locus level. Matches our ART benchmark. |
| Pockrandt et al. 2020, Bioinformatics (GenMap). [10.1093/bioinformatics/btaa222](https://doi.org/10.1093/bioinformatics/btaa222) | Fast (k, e)-mappability; also across several genomes. Tool for per-copy theoretical mappability (bioconda). |
| Karimzadeh et al. 2018, NAR (Umap/Bismap). [10.1093/nar/gky677](https://doi.org/10.1093/nar/gky677) | Ready-made mappability tracks for hg19/hg38/mm9/mm10. Bismap covers bisulfite data. |
| Cechova 2020, Genes (review). [10.3390/genes12010048](https://doi.org/10.3390/genes12010048) | Review: read length and genome complexity decide multimapping; long reads and T2T as the remedy. |

## 2. Allocation of multimapped reads: chromatin (ChIP, ATAC, CUT&RUN)

| Paper | Approach | Key point for us |
|---|---|---|
| Chung et al. 2011, PLoS Comput Biol (CSEM). [10.1371/journal.pcbi.1002111](https://doi.org/10.1371/journal.pcbi.1002111) | EM, fractional allocation weighted by local read density | Multi-reads add up to 30% depth and new peaks, mostly in segmental duplications. |
| Zhang & Keleş 2014, Bioinformatics (cnvCSEM). [10.1093/bioinformatics/btu402](https://doi.org/10.1093/bioinformatics/btu402) | CSEM initialised with copy-number information | Copy-number variation biases density-based allocation. |
| Zeng et al. 2015, PLoS Comput Biol (Perm-seq). [10.1371/journal.pcbi.1004491](https://doi.org/10.1371/journal.pcbi.1004491) | Allocation with priors from other data (e.g. DNase, histone ChIP) | Priors improve allocation accuracy in highly repetitive regions. |
| Shah & Ruthenburg 2021, PLoS Comput Biol (SmartMap). [10.1371/journal.pcbi.1008926](https://doi.org/10.1371/journal.pcbi.1008926) | Bayesian weights per alignment from read distribution and alignment quality, on top of standard aligners | MNase/ChIP/ATAC: up to 53% more depth, +18% of the genome analysable, >140,000 additional repeats; fast. |
| Almeida da Paz & Taher 2022, Mobile DNA (T3E). [10.1186/s13100-022-00285-z](https://doi.org/10.1186/s13100-022-00285-z) | Family-level enrichment: each mapping weighted by 1/(number of loci), background from the **input** | Discarding multimappers underestimates young families. Random-permutation backgrounds give false positives and negatives. Directly relevant to our family-level ChIP counts. |
| Morrissey et al. 2024, Genome Res (Allo). [10.1101/gr.278638.123](https://doi.org/10.1101/gr.278638.123) | Probabilistic allocation plus a CNN that recognises peak-shaped read distributions; writes a corrected BAM | Thousands of new CTCF peaks; most useful in young TEs, centromeres and segmental duplications. A drop-in BAM fits our pipeline structure. Needs multiple alignments per read (bowtie2 `-k`). Trained on peak-type data (TF/ATAC), not broad marks. |

## 3. TE quantification: RNA-seq

| Paper | Level | Key point for us |
|---|---|---|
| Jin et al. 2015, Bioinformatics (TEtranscripts). [10.1093/bioinformatics/btv422](https://doi.org/10.1093/bioinformatics/btv422) | family (EM) | Standard family-level tool. TElocal (same lab) is the locus-level variant; no PubMed paper. |
| Jeong et al. 2018, Pac Symp Biocomput (SalmonTE) | family | Fast, salmon-based. (No DOI in PubMed.) |
| Yang et al. 2019, NAR (SQuIRE). [10.1093/nar/gky1301](https://doi.org/10.1093/nar/gky1301) | locus (EM) | First locus-specific pipeline; mouse tissues. |
| Bendall et al. 2019, PLoS Comput Biol (Telescope). [10.1371/journal.pcbi.1006453](https://doi.org/10.1371/journal.pcbi.1006453) | locus (Bayesian reassignment) | Locus-level HERV expression. |
| Schwarz et al. 2022, Brief Bioinform. [10.1093/bib/bbab417](https://doi.org/10.1093/bib/bbab417) | benchmark | Locus-level differential expression works well for PE data; SalmonTE (slightly modified) and Telescope perform best. |
| Savytska et al. 2022, Front Genet. [10.3389/fgene.2022.1026847](https://doi.org/10.3389/fgene.2022.1026847) | benchmark | **At locus level, false positives exceed the true active loci** for all tools, including SQuIRE, TElocal, SalmonTE and featureCounts unique/fraction/random. Count filtering and TSS profiling help. A warning for element-level RNA claims. |
| Ansaloni et al. 2022, Bioinformatics (TEspeX). [10.1093/bioinformatics/btac526](https://doi.org/10.1093/bioinformatics/btac526) | consensus | Excludes reads from TE fragments embedded in genes and transcripts (exonised TEs), which otherwise inflate apparent TE expression. |
| Lee et al. 2025, Genome Biol (LocusMasterTE). [10.1186/s13059-025-03522-9](https://doi.org/10.1186/s13059-025-03522-9) | locus (EM + long-read prior) | Long-read RNA-seq as a prior improves short-read locus assignment. Only applicable with matched long reads. |
| Tabaro et al. 2024, Brief Bioinform (3t-seq). [10.1093/bib/bbae467](https://doi.org/10.1093/bib/bbae467) | pipeline | Snakemake pipeline for genes + TEs + tRNAs; a point of comparison for pipeline design. |
| Single-cell: IRescue (Polimeni 2024, NAR, [10.1093/nar/gkae793](https://doi.org/10.1093/nar/gkae793)); MATES (Wang 2024, Nat Commun, [10.1038/s41467-024-53114-7](https://doi.org/10.1038/s41467-024-53114-7)); Stellarscope (Reyes-Gopar 2025, Cell Rep Methods, [10.1016/j.crmeth.2025.101086](https://doi.org/10.1016/j.crmeth.2025.101086)) | sc, family / locus | Not needed for bulk. MATES uses read context (deep learning) for locus allocation and also handles other modalities. |

## 4. Reference genomes (T2T)

| Paper | Key point for us |
|---|---|
| Nurk et al. 2022, Science (T2T-CHM13). [10.1126/science.abj6987](https://doi.org/10.1126/science.abj6987) | About 200 Mb of new sequence: centromeric satellites, segmental duplications, acrocentric arms. |
| Aganezov et al. 2022, Science. [10.1126/science.abl3533](https://doi.org/10.1126/science.abl3533) | CHM13 improves short-read mapping and variant calling for 3,202 samples and removes tens of thousands of spurious variants per sample. That is the **cross-reference** effect our own benchmark has not measured yet. |
| Hoyt et al. 2022, Science. [10.1126/science.abk3112](https://doi.org/10.1126/science.abk3112) | De novo repeat annotation of CHM13, new satellite arrays, transcriptionally active retroelements. Resource for hs1 repeat annotation. |
| Gershman et al. 2022, Science. [10.1126/science.abj5089](https://doi.org/10.1126/science.abj5089) | 166,058 previously unresolved ChIP-seq peaks once short-read data are mapped to CHM13; paralogue-specific regulation. |
| Rhie et al. 2023, Nature (T2T-Y). [10.1038/s41586-023-06457-y](https://doi.org/10.1038/s41586-023-06457-y) | Complete human Y (HG002); more than 30 Mb added to GRCh38-Y. hs1 includes it. |
| Liu et al. 2024, Science (mhaESC). [10.1126/science.adq8191](https://doi.org/10.1126/science.adq8191) | T2T mouse genome from haploid ESCs (C57BL/6); more than 7.7% new sequence (rDNA, pericentromeres, subtelomeres). **Our `mhaESC` genome.** |
| Li et al. 2026, Science (mT2T-Y). [10.1126/science.aea2249](https://doi.org/10.1126/science.aea2249) | Complete mouse Y chromosome (95.21 Mb, C57BL/6), 8.7 Mb new, 142 new genes, PAR recombination loci; combined with mhaESC as "T2T mhaESC+Y", a complete C57BL/6 reference. **The chrY in our mhaESC build (v1.1).** |
| Francis et al. 2025, Nat Genet. [10.1038/s41588-025-02367-z](https://doi.org/10.1038/s41588-025-02367-z) | **An independent T2T assembly for C57BL/6J, plus CAST/EiJ**: 213 Mb new sequence, 517 protein-coding genes, PAR boundary, KZFP loci. An alternative to mhaESC to consider. CAST/EiJ is useful for F1 hybrid (allele-specific) designs. |
| Packiaraj & Thakur 2024, Genome Biol. [10.1186/s13059-024-03184-z](https://doi.org/10.1186/s13059-024-03184-z) | Mouse minor and major satellites are heterogeneous in sequence and organisation; H3K9me3/CENP-A ChIP on long-read B6 assemblies. |

## 5. Strain variation and polymorphic TE insertions (mouse)

| Paper | Key point for us |
|---|---|
| Nellåker et al. 2012, Genome Biol. [10.1186/gb-2012-13-6-r45](https://doi.org/10.1186/gb-2012-13-6-r45) | 103,798 polymorphic TE variants across 17 strains (short reads). |
| Lilue et al. 2018, Nat Genet. [10.1038/s41588-018-0223-8](https://doi.org/10.1038/s41588-018-0223-8) | De novo assemblies of 16 laboratory strains; the most divergent regions are enriched in TEs and recent retrotransposition. Includes the Mouse Genomes Project strains; check that 129S1/SvImJ is among them before using it. |
| Ferraj et al. 2023, Cell Genomics. [10.1016/j.xgen.2023.100291](https://doi.org/10.1016/j.xgen.2023.100291) | Long-read SVs in 20 strains: 413,758 SVs over 13% of the reference. **TEs are 39% of SVs and 75% of altered bases. TE heterogeneity changes chromatin accessibility in mouse ESCs.** |
| Helmy et al. 2025, Cell Genomics. [10.1016/j.xgen.2025.101074](https://doi.org/10.1016/j.xgen.2025.101074) | 17 long-read, annotated strain genomes; strain-specific annotation improves RNA-seq mapping (+5.1% for PWK). |
| López-Cortegano et al. 2025, Genome Res. [10.1101/gr.279982.124](https://doi.org/10.1101/gr.279982.124) | De novo mutations in inbred lines: TE insertions are the second most frequent structural mutation, so even "pure" B6 colonies acquire new TE copies. |

## 6. Long-read alternatives

| Paper | Key point for us |
|---|---|
| Ewing et al. 2020, Mol Cell. [10.1016/j.molcel.2020.10.024](https://doi.org/10.1016/j.molcel.2020.10.024) | Nanopore gives locus-specific CpG methylation of young L1/SVA, including non-reference insertions. A ground truth for copy-level chromatin claims. |
| Smits et al. 2023, Methods Mol Biol. [10.1007/978-1-0716-2883-6_9](https://doi.org/10.1007/978-1-0716-2883-6_9) | ONT protocol for polymorphic TE insertions and their methylation. |

## 7. Duplicates and multimappers

No dedicated method was found in PubMed for PCR-duplicate detection among multimapped reads. Our open issue 2 (a random locus per duplicate copy) remains a design decision. Options are listed in `repeat_strategy_open_questions.md`: deduplication before alignment at the sequence or fragment level, or UMIs.

## What this means for the NGS.analysis repeat strategy

1. **Chromatin, element level.**
   - Signal-aware allocation (Allo, SmartMap) is the current state of the art for peak-type data. Random placement, as in our repeats mode, spreads signal across copies.
   - A candidate is an **Allo or SmartMap branch** built on bowtie2 `-k`, benchmarked with our ART set. Allo's CNN suits TF/ATAC data; its suitability for broad H3K9me3 domains is untested.
2. **Chromatin, family level.**
   - T3E-style 1/n weighting with an input-based background is a principled alternative to random-1 plus featureCounts.
   - It addresses the young-family underestimation that our STAR random-1 avoids, and adds a proper input normalisation that we do not do yet.
3. **RNA, element level.**
   - The EM tools (Telescope, SalmonTE) are best in Schwarz 2022.
   - However, Savytska 2022 shows that **locus-level false discoveries can outnumber true active loci**. Any element-level RNA output needs count filters, and ideally TSS evidence.
   - TEspeX's exonised-fragment issue also applies to our family counts.
4. **Reference.**
   - Two T2T C57BL/6 references exist: mhaESC + mT2T-Y (Westlake; our build) and Francis 2025 (B6J + CAST).
   - Choose one for v3, or support both. The key real-data benefit is the cross-reference effect (Aganezov 2022), which our benchmark has not tested yet.
5. **Strain background.**
   - ES-cell lines are often 129 or mixed background.
   - Non-reference TE insertions (Ferraj 2023: 75% of SV bases are TEs) send reads to reference paralogues and create false copy-level signal.
   - **Options:**
     - strain-specific assemblies (Lilue 2018, Helmy 2025);
     - masking known polymorphic loci;
     - at least flagging copies within known SVs.
6. **Validation.**
   - Nanopore methylation (Ewing 2020) or long-read data on our own lines would give a ground truth for copy-level claims.
   - This would be worth one pilot sample before committing to an element-level strategy.
