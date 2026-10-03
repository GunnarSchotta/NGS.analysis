#!/bin/bash
#SBATCH --job-name=ngsv3_benchmark
#SBATCH --partition=slim16
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=/store24/project24/becgsc_001/analysis/NGS.analysis.v3.benchmark/logs/ngsv3_benchmark_%j.out
#SBATCH --error=/store24/project24/becgsc_001/analysis/NGS.analysis.v3.benchmark/logs/ngsv3_benchmark_%j.err
# Usage: sbatch run_benchmark.sh <genome>   (mm10 | hg38 | hs1 | mhaESC; resources from genomes/<genome>.yaml)
# 1) ART (Teissandier-style): chr1, HiSeq2500 profile, PE 2x100 (frag 200+-20) and 2x50 (frag 150+-20), 5x
#    aligned with: v2 genes (bowtie2 --very-sensitive -X 2000, unique = MAPQ 42), v3 genes (+ --dovetail, MAPQ 30),
#    repeats mode STAR (v2 = v3 settings; unique = MAPQ 255); per-family metrics for fragments overlapping RepeatMasker
# 2) adapter read-through set (simulate_adapters.py, 300k fragments 30-150 bp, 60-bp reads, Nextera):
#    v2 trimming + v2 bowtie2 vs v3 trimming + v3 bowtie2 -> recovery per fragment length
set -euo pipefail
G=$1
E=/store24/project24/becgsc_001/micromamba/envs/ngs.v3/bin
ART=/store24/project24/becgsc_001/micromamba/envs/art/bin/art_illumina
CODE=/store24/project24/becgsc_001/coding/NGS.analysis.v3
B=/store24/project24/becgsc_001/analysis/NGS.analysis.v3.benchmark/$G
export PATH=$E:$PATH
mkdir -p $B && cd $B
y() { $E/python3 -c "import yaml;print(yaml.safe_load(open('$CODE/genomes/$G.yaml'))['$1'])"; }
FA=$(y fasta); BT2=$(y bowtie2_index); STARIDX=$(y star_genome_index); SAF=$(y repeats_saf)
N=16
# rmsk BED (name in col 4, class col 5) from the family SAF (class: name only, SAF has no class)
[ -s rmsk.bed ] || awk 'NR>1 && $2=="chr1"{OFS="\t"; s=$3-1; if (s<0) s=0; print $2,s,$4,$1,$1,$5}' $SAF | sort -k1,1 -k2,2n > rmsk.bed  # chr1 (simulated) only
[ -s chr1.fa ] || samtools faidx $FA chr1 > chr1.fa

align_all() {  # $1 = prefix of R1/R2 fastq
  local P=$1
  [ -s $P.v2genes.bam ] || bowtie2 -p $N --very-sensitive -X 2000 -x $BT2 -1 ${P}_R1.fq -2 ${P}_R2.fq 2> $P.v2genes.log | samtools view -b -o $P.v2genes.bam -
  [ -s $P.v3genes.bam ] || bowtie2 -p $N --very-sensitive -X 2000 --dovetail -x $BT2 -1 ${P}_R1.fq -2 ${P}_R2.fq 2> $P.v3genes.log | samtools view -b -o $P.v3genes.bam -
  if [ ! -s $P.repeats.bam ]; then
    STAR --runThreadN $N --outSAMtype BAM Unsorted --runMode alignReads --outFilterMultimapNmax 5000 \
      --outSAMmultNmax 1 --outFilterMismatchNmax 3 --outMultimapperOrder Random --winAnchorMultimapNmax 5000 \
      --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 350 --seedSearchStartLmax 30 \
      --alignTranscriptsPerReadNmax 30000 --alignWindowsPerReadNmax 30000 --alignTranscriptsPerWindowNmax 300 \
      --seedPerReadNmax 3000 --seedPerWindowNmax 300 --seedNoneLociPerWindow 1000 --genomeDir $STARIDX \
      --readFilesIn ${P}_R1.fq ${P}_R2.fq --outFileNamePrefix $P.star.
    mv $P.star.Aligned.out.bam $P.repeats.bam
  fi
}

for L in 100 50; do
  P=art_${L}; M=$([ $L = 100 ] && echo 200 || echo 150)
  if [ ! -s ${P}_R1.fq ]; then
    $ART -ss HS25 -i chr1.fa -p -l $L -f 5 -m $M -s 20 -sam -na -rs 42 -o ${P}_ > /dev/null
    mv ${P}_1.fq ${P}_R1.fq; mv ${P}_2.fq ${P}_R2.fq; mv ${P}_.sam ${P}.truth.sam
  fi
  align_all $P
  for a in v2genes:q42 v3genes:q30 repeats:q255 repeats_any:q0; do
    lab=${a%%:*}; q=${a##*:}; bam=$P.${lab%_any}.bam
    [ -s eval_${P}_${lab}.tsv ] || $E/python3 $CODE/benchmark/evaluate.py art $P.truth.sam rmsk.bed $bam "$G:$P:$lab:q$q" eval_${P}_${lab}.tsv
  done
done

# adapter read-through set
P=adapt60
[ -s ${P}_R1.fastq ] || $E/python3 $CODE/benchmark/simulate_adapters.py $FA chr1 300000 60 $P
AD=$CODE/NexteraPE-PE.fa
[ -s ${P}.v2_R1.fq ] || trimmomatic PE -phred33 -threads $N ${P}_R1.fastq ${P}_R2.fastq ${P}.v2_R1.fq /dev/null ${P}.v2_R2.fq /dev/null ILLUMINACLIP:$AD:2:30:10 MINLEN:30 2> ${P}.v2.trim.log
[ -s ${P}.v3_R1.fq ] || trimmomatic PE -phred33 -threads $N ${P}_R1.fastq ${P}_R2.fastq ${P}.v3_R1.fq /dev/null ${P}.v3_R2.fq /dev/null ILLUMINACLIP:$AD:2:30:10:1:true MINLEN:30 2> ${P}.v3.trim.log
[ -s ${P}.v2.bam ] || bowtie2 -p $N --very-sensitive -X 2000 -x $BT2 -1 ${P}.v2_R1.fq -2 ${P}.v2_R2.fq 2> ${P}.v2.bt2.log | samtools view -b -o ${P}.v2.bam -
[ -s ${P}.v3.bam ] || bowtie2 -p $N --very-sensitive -X 2000 --dovetail -x $BT2 -1 ${P}.v3_R1.fq -2 ${P}.v3_R2.fq 2> ${P}.v3.bt2.log | samtools view -b -o ${P}.v3.bam -
$E/python3 $CODE/benchmark/evaluate.py adapt ${P}.v2.bam "$G:adapt60:v2trim_v2bt2:q42" eval_adapt_v2.tsv ${P}_R1.fastq
$E/python3 $CODE/benchmark/evaluate.py adapt ${P}.v3.bam "$G:adapt60:v3trim_v3bt2:q30" eval_adapt_v3.tsv ${P}_R1.fastq
echo "[$(date)] benchmark $G done"
