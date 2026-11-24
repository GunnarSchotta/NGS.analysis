#!/bin/bash
# T1 NGS.peaks smoke test (since 3.2.0): control-based peaks on the finished v3_genes T1 run.
# Adds a control column (input_ES) for chip_K9 (broad), ct_k27ac (broad since 3.2) and atac_mm (narrow + --call-summits)
# in a separate PEP (the main T1 sample table has no control column) and runs NGS.peaks via .looper.peaks.yaml.
# Run after make_v3_runs.sh has finished. Usage: make_peaks_run.sh [submit]
set -euo pipefail
V=/store24/project24/becgsc_001/analysis/NGS.analysis.v3.validation/T1
V3CODE=/store24/project24/becgsc_001/coding/NGS.analysis.v3
LOOPER=/store24/project24/becgsc_001/micromamba/envs/ngs.v3/bin/looper
R=$V/runs/v3_genes
cd "$R"
awk -F, 'BEGIN{OFS=","} NR==1{print $0,"control"; next}
         {c=($1=="chip_K9"||$1=="ct_k27ac"||$1=="atac_mm")?"input_ES":""; print $0,c}' sample.table.csv > sample.table.peaks.csv
sed 's/^sample_table: .*/sample_table: sample.table.peaks.csv/' analysis.configuration.yaml > analysis.configuration.peaks.yaml
cat > .looper.peaks.yaml <<YAML
pep_config: analysis.configuration.peaks.yaml
output_dir: results
pipeline_interfaces:
   - $V3CODE/peaks_pipeline_interface.yaml
pipestat:
  results_file_path: "{record_identifier}/NGS.peaks_stats.yaml"
  flag_file_dir: results/flags_peaks
YAML
if [ "${1:-}" = submit ]; then
  PATH=$(dirname $LOOPER):$PATH $LOOPER run -p slurm -c .looper.peaks.yaml 2>&1 | grep -E "Submitted|valid" | tail -2
fi
