#!/bin/bash
# T2 (full-depth real data):
#  eden_genes   : 4 Edenhofer ATAC libraries, v3 genes mode -> compare with analysis/Edenhofer/02.retrim.cutadapt
#  angela_repeats: 1 RNA (Setdb1 ctrl d0) + 1 H3K9me3 ChIP (wt26 ES), v3 repeats mode -> compare with existing v2 outputs
set -euo pipefail
V=/store24/project24/becgsc_001/analysis/NGS.analysis.v3.validation/T2
V3=/store24/project24/becgsc_001/coding/NGS.analysis.v3
LOOPER=/store24/project24/becgsc_001/micromamba/envs/ngs.v3/bin/looper
mk() {  # name mode sample_table_body sources_dir_R1 sources_dir_R2
  local R=$V/runs/$1; mkdir -p $R
  printf "sample_name,study,cell_type,target,genotype,batch,protocol,organism,fastq1,fastq2,read_type,read1,read2\n%s" "$3" > $R/sample.table.csv
  cat > $R/analysis.configuration.yaml <<YAML
name: T2_$1
pep_version: 2.1.0
sample_table: sample.table.csv
pipeline_mode: $2
sample_modifiers:
  derive:
    attributes: [read1, read2]
    sources:
      R1: "{src1}"
      R2: "{src2}"
  imply:
    - if: {protocol: RNA}
      then: {strandedness: reverse}
    - if: {organism: mouse}
      then: {genome: mm10}
    - if: {organism: human}
      then: {genome: hg38}
YAML
  cat > $R/.looper.yaml <<YAML
pep_config: analysis.configuration.yaml
output_dir: results
pipeline_interfaces:
   - $V3/sample_pipeline_interface.yaml
   - $V3/project_pipeline_interface.yaml
pipestat:
  results_file_path: "{record_identifier}/stats.yaml"
  flag_file_dir: results/flags
YAML
}
E=/store24/project24/becgsc_020/data/Edenhofer_ATAC
mk eden_genes genes "eNSPC_LP_A_rep1,T2,eNSPC,none,wt,1,ATAC,human,$E/read1_GS2071.fastq.gz,$E/read2_GS2071.fastq.gz,paired,R1,R2
eNSPC_HP_A2_rep2,T2,eNSPC,none,wt,1,ATAC,human,$E/read1_GS2076.fastq.gz,$E/read2_GS2076.fastq.gz,paired,R1,R2
NHDF_rep1,T2,iPSC,none,wt,1,ATAC,human,$E/read1_GS2077.fastq.gz,$E/read2_GS2077.fastq.gz,paired,R1,R2
smNPC_rep1,T2,smNPC,none,wt,1,ATAC,human,$E/read1_GS2081.fastq.gz,$E/read2_GS2081.fastq.gz,paired,R1,R2
"
A=/store24/project24/becgsc_001/schottalab/01.Angela/00.raw.data/00.paper
RN=$A/RNAseq/01.ctrl.Setdb1KO/250827_VL00118_536_AAHCVCKM5.1rep_251002_VL00118_548_AAHGLWLM5.2.3.reps
mk angela_repeats repeats "RNA.T263.ES.d0_r1,T2,ESC,none,wt,1,RNA,mouse,$RN/read1_GS2152.fastq.gz,$RN/read2_GS2152.fastq.gz,paired,R1,R2
wt26.ES.H3K9me3_r1,T2,ESC,H3K9me3,wt,1,CHIP,mouse,$A/Histones/01.ES.XEN.wt26/read1_GS1807.fastq.gz,$A/Histones/01.ES.XEN.wt26/read2_GS1807.fastq.gz,paired,R1,R2
"
for r in eden_genes angela_repeats; do
  sed -i 's#"{src1}"#"{fastq1}"#; s#"{src2}"#"{fastq2}"#' $V/runs/$r/analysis.configuration.yaml
  (cd $V/runs/$r && PATH=$(dirname $LOOPER):$PATH $LOOPER run -p slurm 2>&1 | grep -c Submitted)
done
