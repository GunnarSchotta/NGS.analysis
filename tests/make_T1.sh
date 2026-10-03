#!/bin/bash
#SBATCH --job-name=ngsv3_T1_subsample
#SBATCH --partition=slim16
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=6:00:00
#SBATCH --array=1-8
#SBATCH --output=/store24/project24/becgsc_001/analysis/NGS.analysis.v3.validation/logs/ngsv3_T1_subsample_%A_%a.out
#SBATCH --error=/store24/project24/becgsc_001/analysis/NGS.analysis.v3.validation/logs/ngsv3_T1_subsample_%A_%a.err
# T1 smoke-test data: 200k reads/pairs per sample, reservoir sample (seed 42) of one real library per protocol.
set -euo pipefail
OUT=/store24/project24/becgsc_001/analysis/NGS.analysis.v3.validation/T1/fastq
PY=/store24/project24/becgsc_001/micromamba/envs/ngs.v3/bin/python3
S=/store24/project24/becgsc_001/coding/NGS.analysis.v3/tests/subsample_fastq.py
A=/store24/project24/becgsc_001/schottalab/01.Angela/00.raw.data/00.paper
D=/store24/project24/becgsc_001/schottalab/00.external.datasets/data
case $SLURM_ARRAY_TASK_ID in
 1) N=chip_K9;  F="$A/Histones/01.ES.XEN.wt26/read1_GS1807.fastq.gz $A/Histones/01.ES.XEN.wt26/read2_GS1807.fastq.gz";;
 2) N=input_ES; F="$A/Histones/01.ES.XEN.wt26/read1_GS1815.fastq.gz $A/Histones/01.ES.XEN.wt26/read2_GS1815.fastq.gz";;
 3) N=atac_mm;  F="$A/ATACseq/01.ES.XEN.wt26/read1_GS1677.fastq.gz $A/ATACseq/01.ES.XEN.wt26/read2_GS1677.fastq.gz";;
 4) N=cr_gata6; F="$D/Thompson.2022.iXEN/TFs/SRR15295347_GSM5484568_Gata6_0h_ixen_CR_DL_12_rep1_Mus_musculus_OTHER_1.fastq.gz $D/Thompson.2022.iXEN/TFs/SRR15295347_GSM5484568_Gata6_0h_ixen_CR_DL_12_rep1_Mus_musculus_OTHER_2.fastq.gz";;
 5) N=ct_k27ac; F="$D/Thompson.2022.iXEN/Histones/SRR15295333_GSM5484554_DL_214b_CT_H3K27ac_0h_rep1_Mus_musculus_OTHER_1.fastq.gz $D/Thompson.2022.iXEN/Histones/SRR15295333_GSM5484554_DL_214b_CT_H3K27ac_0h_rep1_Mus_musculus_OTHER_2.fastq.gz";;
 6) N=chip_se;  F="$D/Cernilogar.2019/SRR7427636_GSM3223314_H3K9me3.d0.r1_Mus_musculus_ChIP-Seq.fastq.gz";;
 7) N=rna;      F="$A/RNAseq/01.ctrl.Setdb1KO/250827_VL00118_536_AAHCVCKM5.1rep_251002_VL00118_548_AAHGLWLM5.2.3.reps/read1_GS2152.fastq.gz $A/RNAseq/01.ctrl.Setdb1KO/250827_VL00118_536_AAHCVCKM5.1rep_251002_VL00118_548_AAHGLWLM5.2.3.reps/read2_GS2152.fastq.gz";;
 8) N=atac_hs;  F="/store24/project24/becgsc_020/data/Edenhofer_ATAC/read1_GS2071.fastq.gz /store24/project24/becgsc_020/data/Edenhofer_ATAC/read2_GS2071.fastq.gz";;
esac
for f in $F; do [ -s "$f" ] || { echo "missing $f"; exit 1; }; done
"$PY" "$S" 200000 42 "$OUT/$N" $F
