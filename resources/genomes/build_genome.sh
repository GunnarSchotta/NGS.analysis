#!/bin/bash
# Build T2T genome resources for NGS.analysis v3.  Usage: build_genome.sh {hs1|mhaESC} submit
# Submits: prep (FASTA, GTF, TSS, RepeatMasker SAFs; prep_t2t.py) -> bowtie2-build | STAR genome (no GTF) |
#          STAR RNA (GTF, sjdbOverhang 100) | RSEM, each after prep. Outputs in genomes/<g>/ (mm10/hg38 layout).
set -euo pipefail
G=$1
E=/store24/project24/becgsc_001/micromamba/envs/ngs.v3/bin
H=/store24/project24/becgsc_001/coding/NGS.analysis.v3/resources/genomes
O=/store24/project24/becgsc_001/genomes/$G
L=/store24/project24/becgsc_001/genomes/build_logs; mkdir -p $L $O
FA=$O/Sequence/WholeGenomeFasta/genome.fa
GTF=$O/Annotation/Genes/genes.gtf
sb() { sbatch --parsable -p slim16 -J "t2t_${G}_$1" -c $2 --mem=$3 -t $4 -o $L/t2t_${G}_$1_%j.out -e $L/t2t_${G}_$1_%j.err ${5:-} --wrap "$6"; }
P=$(sb prep 2 48G 6:00:00 "" "$E/python3 $H/prep_t2t.py $G $O")
D="--dependency=afterok:$P"
sb bowtie2 16 64G 24:00:00 "$D" "mkdir -p $O/Sequence/Bowtie2Index && PATH=$E:\$PATH $E/bowtie2-build --threads 16 $FA $O/Sequence/Bowtie2Index/genome"
sb stargenome 16 64G 24:00:00 "$D" "mkdir -p $O/STARgenome && $E/STAR --runMode genomeGenerate --runThreadN 16 --genomeDir $O/STARgenome --genomeFastaFiles $FA --outFileNamePrefix $O/STARgenome/"
sb starrna 16 64G 24:00:00 "$D" "mkdir -p $O/STAR && $E/STAR --runMode genomeGenerate --runThreadN 16 --genomeDir $O/STAR --genomeFastaFiles $FA --sjdbGTFfile $GTF --sjdbOverhang 100 --outFileNamePrefix $O/STAR/"
sb rsem 8 32G 24:00:00 "$D" "mkdir -p $O/RSEM && PATH=$E:\$PATH $E/rsem-prepare-reference --gtf $GTF -p 8 $FA $O/RSEM/RSEM"
echo "submitted $G (prep job $P)"
