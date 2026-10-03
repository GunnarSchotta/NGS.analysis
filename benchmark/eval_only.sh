#!/bin/bash
# Re-run the ART evaluations of one genome (after fixing rmsk.bed). Usage: sbatch eval_only.sh <genome> <readlen>
#SBATCH --partition=slim16
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=6:00:00
#SBATCH --job-name=ngsv3_bench_eval
#SBATCH --output=/store24/project24/becgsc_001/analysis/NGS.analysis.v3.benchmark/logs/ngsv3_bench_eval_%j.out
#SBATCH --error=/store24/project24/becgsc_001/analysis/NGS.analysis.v3.benchmark/logs/ngsv3_bench_eval_%j.err
set -euo pipefail
G=$1; L=$2; E=/store24/project24/becgsc_001/micromamba/envs/ngs.v3/bin; export PATH=$E:$PATH
CODE=/store24/project24/becgsc_001/coding/NGS.analysis.v3
cd /store24/project24/becgsc_001/analysis/NGS.analysis.v3.benchmark/$G
P=art_$L
for a in v2genes:42 v3genes:30 repeats:255 repeats_any:0; do
  lab=${a%%:*}; q=${a##*:}
  $E/python3 $CODE/benchmark/evaluate.py art $P.truth.sam rmsk.bed $P.${lab%_any}.bam "$G:$P:$lab:q$q" eval_${P}_${lab}.tsv
done
