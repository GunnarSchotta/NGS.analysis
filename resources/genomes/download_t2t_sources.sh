#!/bin/bash
# Login node only (compute nodes have no internet). Downloads the T2T source files used by build_genome.sh.
#   hs1 (T2T-CHM13v2.0): analysis set with PAR-masked Y and rCRS chrM (T2T consortium, AWS),
#                        CAT/Liftoff GENCODE v35 gene annotation (T2T consortium), UCSC hs1 RepeatMasker (.out).
#   mhaESC: Ensembl GRCm39 release 110 GTF, only to map Liftoff ENSMUST transcript IDs to gene symbols.
set -euo pipefail
G=/store24/project24/becgsc_001/genomes
mkdir -p $G/hs1/source $G/mhaESC/source
cd $G/hs1/source
wget -q -c https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz
wget -q -c https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13.draft_v2.0.gene_annotation.gff3
wget -q -c http://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.repeatMasker.out.gz
gzip -f chm13.draft_v2.0.gene_annotation.gff3
cd $G/mhaESC/source
wget -q -c https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz
md5sum $G/hs1/source/* $G/mhaESC/source/* > $G/T2T_sources.md5
date > $G/T2T_sources.download_done
