#!/bin/bash
#SBATCH --job-name=ngsv3_snpbench
#SBATCH --partition=slim16
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --array=0-3
#SBATCH --output=/store24/project24/becgsc_001/analysis/NGS.analysis.v3.benchmark/logs/ngsv3_snpbench_%A_%a.out
#SBATCH --error=/store24/project24/becgsc_001/analysis/NGS.analysis.v3.benchmark/logs/ngsv3_snpbench_%A_%a.err
# SNP/indel and non-reference TE insertion benchmark (roadmap H2), mm10 chr1:
#   array task = SNP rate 0 / 0.0005 / 0.002 / 0.007 (indels at 1/10 of it; 300 non-reference young-TE copies for
#   rates > 0) -> ART 2x100 (200+-20) and 2x50 (150+-20), 5x, from the mutated chr1 -> bowtie2 v3 settings against the
#   unmodified whole mm10 -> eval_snp.py (MAPQ curves per family, nonref placements).
# Output: analysis/NGS.analysis.v3.benchmark/snp/rate_<r>/ ; idempotent.
set -euo pipefail
E=/store24/project24/becgsc_001/micromamba/envs/ngs.v3/bin
ART=/store24/project24/becgsc_001/micromamba/envs/art/bin/art_illumina
CODE=/store24/project24/becgsc_001/coding/NGS.analysis.v3
B=/store24/project24/becgsc_001/analysis/NGS.analysis.v3.benchmark
export PATH=$E:$PATH
RATES=(0 0.0005 0.002 0.007)
R=${RATES[$SLURM_ARRAY_TASK_ID]}
NINS=$([ "$R" = 0 ] && echo 0 || echo 300)
D=$B/snp/rate_$R; mkdir -p "$D"; cd "$D"
BT2=$($E/python3 -c "import yaml;print(yaml.safe_load(open('$CODE/genomes/mm10.yaml'))['bowtie2_index'])")
[ -s sample.fa ] || $E/python3 $CODE/benchmark/mutate_genome.py $B/mm10/chr1.fa $B/mm10/rmsk.bed "$R" "$NINS" 42 sample | tee mutate.log
for spec in 100:200 50:150; do
  L=${spec%%:*}; M=${spec##*:}; P=art_$L
  [ -s ${P}_1.fq ] || $ART -ss HS25 -i sample.fa -p -l $L -f 5 -m $M -s 20 -sam -na -rs 42 -o ${P}_ > /dev/null
  [ -s $P.bt2.bam ] || bowtie2 -p 16 --very-sensitive -X 2000 --dovetail -x "$BT2" -1 ${P}_1.fq -2 ${P}_2.fq 2> $P.bt2.log \
      | samtools view -b -o $P.bt2.bam -
  [ -s eval_$P.tsv ] || $E/python3 $CODE/benchmark/eval_snp.py ${P}_.sam sample.events.tsv $B/mm10/rmsk.bed eval_$P.tsv bowtie2_v3genes=$P.bt2.bam
done
rm -f art_*_[12].fq
echo "rate $R done"
