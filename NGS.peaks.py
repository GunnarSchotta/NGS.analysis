#!/usr/bin/env python3
"""
NGS.peaks (NGS.analysis v3) -- MACS3 peak calling of a sample against its control (input / IgG).

Run AFTER the main NGS.analysis run has completed for the sample and its control:
    looper run -p slurm -c .looper.peaks.yaml
(a copy of .looper.yaml whose pipeline_interfaces lists only <v3>/peaks_pipeline_interface.yaml)
Samples without a `control` column value are skipped (nothing is done, status completed).
Inputs (from the main run): genes mode .filt.bam, repeats mode .dedup.unique.bam of sample and control.
Outputs: <sample>/peaks_ctrl_<genome>/<sample>_peaks[_noBL].{narrow,broad}Peak, FRiP; results in
<sample>/NGS.peaks_stats.yaml (own pipeline name and schema: pipestat_peaks_results_schema.yaml).
"""
__version__ = "3.0.0"

from argparse import ArgumentParser
import os
import sys

import yaml
import pypiper

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ.get("PATH", "")
sys.path.insert(0, SCRIPT_DIR)
import ngs_qc as qc  # noqa: E402

BROAD_TARGETS = {"h3k27me3", "h3k9me3", "h3k36me3", "h3k79me2", "h4k20me3", "h3k9me2", "h2ak119ub"}

p = ArgumentParser(description="NGS.peaks " + __version__)
p = pypiper.add_pypiper_args(p, groups=["pypiper", "looper", "ngs"], required=["sample-name", "output-parent", "genome"])
p.add_argument("--genome-config", default=None)
p.add_argument("--pipeline-mode", required=True, choices=["genes", "repeats"])
p.add_argument("--control", default=None)
p.add_argument("--target", default=None)
p.add_argument("--peak-mode", default="auto", choices=["auto", "narrow", "broad", "none"])
p.add_argument("--se-extend", type=int, default=200)
p.add_argument("--pipestat-config", default=None)
args = p.parse_args()
paired = args.single_or_paired.lower() == "paired"

gcfg_path = args.genome_config or os.path.join(SCRIPT_DIR, "genomes", args.genome_assembly + ".yaml")
GCFG = yaml.safe_load(open(gcfg_path)) if os.path.exists(gcfg_path) else {}
BLACKLIST, GSIZE = GCFG.get("blacklist"), str(GCFG.get("macs_gsize", "hs"))

samples_dir = os.path.abspath(args.output_parent)
outfolder = os.path.join(samples_dir, args.sample_name)


def analysis_bam(sample):
    tag = "filt" if args.pipeline_mode == "genes" else "dedup.unique"
    return os.path.join(samples_dir, sample, "aligned_" + args.genome_assembly, f"{sample}.{tag}.bam")


if args.pipestat_config and os.path.exists(args.pipestat_config):
    cfg = yaml.safe_load(open(args.pipestat_config)) or {}
    cfg["results_file_path"] = os.path.join(outfolder, "NGS.peaks_stats.yaml")
    cfg["schema_path"] = os.path.join(SCRIPT_DIR, "pipestat_peaks_results_schema.yaml")
    os.makedirs(outfolder, exist_ok=True)
    args.pipestat_config = os.path.join(outfolder, "NGS.peaks_pipestat_config.yaml")
    yaml.dump(cfg, open(args.pipestat_config, "w"))

pm = pypiper.PipelineManager(name="NGS.peaks", outfolder=outfolder, pipestat_record_identifier=args.sample_name,
                             pipestat_config_file=args.pipestat_config, args=args, version=__version__)

if not args.control or args.control.lower() in ("none", "na", ""):
    pm.info("No control for this sample; nothing to do.")
    pm.stop_pipeline()
    sys.exit()

mode = args.peak_mode
if mode == "auto":
    mode = "broad" if (args.target or "").lower() in BROAD_TARGETS else "narrow"
t_bam, c_bam = analysis_bam(args.sample_name), analysis_bam(args.control)
for b in (t_bam, c_bam):
    if not os.path.exists(b):
        pm.fail_pipeline(IOError("Missing BAM (has the main pipeline finished?): " + b))
pm.report_result("Control", args.control)
pm.report_result("Peak_mode", mode)

folder = os.path.join(outfolder, "peaks_ctrl_" + args.genome_assembly)
os.makedirs(folder, exist_ok=True)
ptype = "broadPeak" if mode == "broad" else "narrowPeak"
peaks = os.path.join(folder, f"{args.sample_name}_peaks.{ptype}")
peaks_nobl = os.path.join(folder, f"{args.sample_name}_peaks_noBL.{ptype}")
keepdup = "all" if args.pipeline_mode == "genes" else "1"
cmd = (f"macs3 callpeak -t {t_bam} -c {c_bam} -n {args.sample_name} --outdir {folder} -g {GSIZE} -q 0.01 "
       f"--keep-dup {keepdup}" + (" -f BAMPE" if paired else f" -f BAM --nomodel --extsize {args.se_extend}") +
       (" --broad --broad-cutoff 0.1" if mode == "broad" else "") + f" 2> {folder}/{args.sample_name}.macs3.log")
pm.run(cmd, peaks, shell=True)
pm.run((f"bedtools intersect -v -a {peaks} -b {BLACKLIST} > {peaks_nobl}") if BLACKLIST else f"cp {peaks} {peaks_nobl}",
       peaks_nobl, shell=True)
n = qc.count_lines(peaks_nobl)
pm.report_result("Peaks_ctrl_n", n)
if n:
    saf = peaks_nobl + ".saf"
    qc.peaks_to_saf(peaks_nobl, saf)
    fo = os.path.join(folder, args.sample_name + ".frip.txt")
    pm.run(f"featureCounts -F SAF -a {saf} -o {fo} -T {pm.cores} --ignoreDup" + (" -p --countReadPairs " if paired else " ") +
           t_bam, fo + ".summary", nofail=True)
    if os.path.exists(fo + ".summary"):
        s = qc.parse_featurecounts_summary(fo + ".summary")
        usable = sum(s.get(k, 0) for k in ["Assigned", "Unassigned_NoFeatures", "Unassigned_Ambiguity",
                                           "Unassigned_Overlapping_Length"])
        pm.report_result("FRiP_ctrl_pct", round(100 * s.get("Assigned", 0) / usable, 2) if usable else 0)
    pm.report_object("Peaks_ctrl", peaks_nobl)
pm.stop_pipeline()
