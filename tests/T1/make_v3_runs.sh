#!/bin/bash
# v3 runs on T1: genes (new trimming), repeats with --legacy-trim (R1: must reproduce v2), repeats (R2).
# Usage: make_v3_runs.sh [submit]
set -euo pipefail
T=/store24/project24/becgsc_001/coding/NGS.analysis.v3/tests/T1
V=/store24/project24/becgsc_001/analysis/NGS.analysis.v3.validation/T1
V3CODE=/store24/project24/becgsc_001/coding/NGS.analysis.v3
LOOPER=/store24/project24/becgsc_001/micromamba/envs/ngs.v3/bin/looper
for run in v3_genes:genes:no v3_repeats_legacy:repeats:yes v3_repeats:repeats:no; do
  IFS=: read name mode legacy <<< "$run"; R=$V/runs/$name
  rm -rf "$R"; mkdir -p "$R"; cp "$T/sample.table.csv" "$R/"
  APPEND=""; [ "$legacy" = yes ] && APPEND=$'  append:\n    legacy_trim: "true"'
  cat > "$R/analysis.configuration.yaml" <<YAML
name: T1_$name
pep_version: 2.1.0
sample_table: sample.table.csv
pipeline_mode: $mode
sample_modifiers:
$APPEND
  derive:
    attributes: [read1, read2]
    sources:
      R1: "$V/fastq/{fastq1}"
      R2: "$V/fastq/{fastq2}"
  imply:
    - if: {protocol: RNA}
      then: {strandedness: reverse}
    - if: {organism: mouse}
      then: {genome: mm10}
    - if: {organism: human}
      then: {genome: hg38}
YAML
  cat > "$R/.looper.yaml" <<YAML
pep_config: analysis.configuration.yaml
output_dir: results
pipeline_interfaces:
   - $V3CODE/sample_pipeline_interface.yaml
   - $V3CODE/project_pipeline_interface.yaml
pipestat:
  results_file_path: "{record_identifier}/stats.yaml"
  flag_file_dir: results/flags
YAML
  if [ "${1:-}" = submit ]; then (cd "$R" && PATH=$(dirname $LOOPER):$PATH $LOOPER run -p slurm 2>&1 | grep -E "Submitted|valid" | tail -2); fi
done
