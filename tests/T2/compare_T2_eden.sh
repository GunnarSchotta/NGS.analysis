#!/bin/bash
#SBATCH --job-name=ngsv3_T2_eden_compare
#SBATCH --partition=slim16
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=/store24/project24/becgsc_001/analysis/NGS.analysis.v3.validation/logs/ngsv3_T2_eden_compare_%j.out
#SBATCH --error=/store24/project24/becgsc_001/analysis/NGS.analysis.v3.validation/logs/ngsv3_T2_eden_compare_%j.err
# T2 acceptance, genes mode: v3 <s>.filt.bam vs analysis/Edenhofer/02.retrim.cutadapt clean BAMs (cutadapt -O 1 -m 20,
# bowtie2 --dovetail, -f 2 -F 1804 -q 30, chr1-22/X/Y). The 02 BAMs are passed through the same blacklist + fixmate
# orphan removal as v3, so the remaining difference is Trimmomatic vs cutadapt (and MINLEN 30 vs 20).
# Output: T2/eden_vs_02retrim.tsv (fragments = R1 records; length classes from |TLEN|).
set -euo pipefail
export PATH=/store24/project24/becgsc_001/micromamba/envs/ngs.v3/bin:$PATH
V=/store24/project24/becgsc_001/analysis/NGS.analysis.v3.validation/T2
R=$V/runs/eden_genes/results/samples
O=/store24/project24/becgsc_001/analysis/Edenhofer/02.retrim.cutadapt/bam_clean
BL=/store24/project24/becgsc_001/genomes/hg38/GRCh38_unified_blacklist.bed
T=$V/tmp_eden_compare; mkdir -p $T
out=$V/eden_vs_02retrim.tsv
lens() {  # fragments per |TLEN| class: <30, 30-59, 60-149, >=150
  samtools view -f 64 "$1" | awk '{t=$9<0?-$9:$9; if (t<30) a++; else if (t<60) b++; else if (t<150) c++; else d++}
                                  END{printf "%d\t%d\t%d\t%d", a, b, c, d}'
}
printf "sample\tv3_filt\t02_clean\t02_clean_blfix\tdiff_pct_v3_vs_02blfix\tv3_lt30\tv3_30_59\tv3_60_149\tv3_ge150\t02blfix_lt30\t02blfix_30_59\t02blfix_60_149\t02blfix_ge150\n" > $out
for s in $(tail -n +2 $V/runs/eden_genes/sample.table.csv | cut -d, -f1); do
  v3=$R/$s/aligned_hg38/$s.filt.bam
  [ -s "$v3" ] || { echo "skip $s (no v3 filt.bam)"; continue; }
  [ -s $T/$s.blfix.bam ] || samtools view -u -f 2 -F 1804 $O/$s.clean.bam | bedtools intersect -v -abam stdin -b $BL |
    samtools sort -n -@ 8 -m 1G -T $T/$s.n -O BAM - | samtools fixmate -r - - | samtools view -b -f 2 -F 1804 -o $T/$s.blfix.bam -
  a=$(samtools view -c -f 64 $v3); b=$(samtools view -c -f 64 $O/$s.clean.bam); c=$(samtools view -c -f 64 $T/$s.blfix.bam)
  printf "%s\t%s\t%s\t%s\t%.2f\t%s\t%s\n" $s $a $b $c $(echo "100*($a-$c)/$c" | bc -l) "$(lens $v3)" "$(lens $T/$s.blfix.bam)" >> $out
done
cat $out
rm -rf $T
