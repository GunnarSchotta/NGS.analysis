#!/bin/bash
# Step 0 of the v3 plan: v2 baselines on T1 (genes x1, repeats x2 for the random-multimapper noise floor).
# Creates run folders under the validation area and submits with the ngs.v2 env (bare `python` in the v2
# interface templates resolves to the submitting environment).
set -euo pipefail
T=/store24/project24/becgsc_001/coding/NGS.analysis.v3/tests/T1
V=/store24/project24/becgsc_001/analysis/NGS.analysis.v3.validation/T1
V2CODE=/store24/project24/becgsc_001/coding/NGS.analysis
G=/store24/project24/becgsc_001/genomes
for run in v2_genes:genes v2_repeats_A:repeats v2_repeats_B:repeats; do
  name=${run%%:*}; mode=${run##*:}; R=$V/runs/$name
  mkdir -p "$R"; cp "$T/sample.table.csv" "$R/"
  cat > "$R/analysis.configuration.yaml" <<YAML
name: T1_$name
pep_version: 2.1.0
sample_table: sample.table.csv
pipeline_mode: $mode
sample_modifiers:
  derive:
    attributes: [read1, read2]
    sources:
      R1: "$V/fastq/{fastq1}"
      R2: "$V/fastq/{fastq2}"
  imply:
    - if: {protocol: RNA}
      then: {strandedness: reverse}
    - if: {organism: mouse}
      then:
        genome: mm10
        STAR_RNA_index: $G/mm10/STAR
        STAR_genome_index: $G/mm10/STARgenome
        Bowtie2_index: $G/mm10/Sequence/Bowtie2Index/genome
        repeats_SAF: $G/mm10/Annotation/Genes/rmsk.mm10.160322.SAF
        repeats_SAFid: $G/mm10/Annotation/Genes/rmsk.ids.mm10.160322.SAF
        refgene_tss: $G/mm10/Annotation/Genes/mm10_TSS.bed
        genome_index: $G/mm10/Sequence/WholeGenomeFasta/genome.fa.fai
        rsem_index: $G/mm10/RSEM/RSEM
    - if: {organism: human}
      then:
        genome: hg38
        STAR_RNA_index: $G/hg38/STAR
        STAR_genome_index: $G/hg38/STARgenome
        Bowtie2_index: $G/hg38/Sequence/Bowtie2Index/genome
        repeats_SAF: $G/hg38/Annotation/Genes/rmsk.hg38.291122.SAF
        repeats_SAFid: $G/hg38/Annotation/Genes/rmsk.ids.hg38.291122.SAF
        refgene_tss: $G/hg38/Annotation/Genes/hg38_TSS.bed
        genome_index: $G/hg38/Sequence/WholeGenomeFasta/genome.fa.fai
        rsem_index: $G/hg38/RSEM/RSEM
YAML
  cat > "$R/.looper.yaml" <<YAML
pep_config: analysis.configuration.yaml
output_dir: results
pipeline_interfaces:
   - $V2CODE/sample_pipeline_interface.yaml
pipestat:
  results_file_path: "{record_identifier}/stats.yaml"
  flag_file_dir: results/flags
YAML
done
