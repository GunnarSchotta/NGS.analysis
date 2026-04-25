#!/usr/bin/env python3
"""
NGS Analysis Collator - project-level pipeline
"""

__author__ = ["Gunnar Schotta"]
__email__ = "gunnar.schotta@bmc.med.lmu.de"
__version__ = "2.0.0"

from argparse import ArgumentParser
import os
import sys
import yaml
import pypiper
from ubiquerg import VersionInHelpParser

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

def tool_path(tool_name: str) -> str:
    """Return absolute path to a tool next to this script."""
    return os.path.join(SCRIPT_DIR, tool_name)


def _file_result(path: str, title: str) -> dict[str, str]:
    return {"path": path, "title": title}


def _image_result(path: str, thumbnail_path: str, title: str) -> dict[str, str]:
    return {"path": path, "thumbnail_path": thumbnail_path, "title": title}


def collect_project_results(project_name: str, summary_dir: str) -> dict[str, dict[str, str]]:
    """Collect project-level summary artifacts that were created by the R summarizer."""
    summary_path = os.path.abspath(summary_dir)
    candidates = {
        "stat_summary_file": (
            os.path.join(summary_path, f"{project_name}_stats_summary.tsv"),
            "Statistics summary file",
            "file",
        ),
        "feature_counts_file": (
            os.path.join(summary_path, f"{project_name}_fc_summary.rds"),
            "FeatureCounts classes summary",
            "file",
        ),
        "feature_counts_id_file": (
            os.path.join(summary_path, f"{project_name}_fc_id_summary.rds"),
            "FeatureCounts individual elements summary",
            "file",
        ),
        "IAP_coverage_file": (
            os.path.join(summary_path, f"{project_name}_IAP_coverage_summary.rds"),
            "IAP coverage summary",
            "file",
        ),
        "unstranded_gene_counts_file": (
            os.path.join(summary_path, f"{project_name}_unstranded_gene_counts_summary.rds"),
            "Unstranded gene counts summary",
            "file",
        ),
        "sense_gene_counts_file": (
            os.path.join(summary_path, f"{project_name}_sense_gene_counts_summary.rds"),
            "Sense gene counts summary",
            "file",
        ),
        "IAP_coverage_samples": (
            os.path.join(summary_path, f"{project_name}_IAP_coverage_samples.pdf"),
            os.path.join(summary_path, f"{project_name}_IAP_coverage_samples.png"),
            "IAP coverage samples plot",
            "image",
        ),
        "IAP_coverage_merged": (
            os.path.join(summary_path, f"{project_name}_IAP_coverage_merged.pdf"),
            os.path.join(summary_path, f"{project_name}_IAP_coverage_merged.png"),
            "IAP coverage merged plot",
            "image",
        ),
    }

    results: dict[str, dict[str, str]] = {}
    for key, entry in candidates.items():
        if entry[-1] == "file":
            path, title, _ = entry
            if os.path.exists(path):
                results[key] = _file_result(path, title)
        else:
            path, thumbnail_path, title, _ = entry
            if os.path.exists(path) and os.path.exists(thumbnail_path):
                results[key] = _image_result(path, thumbnail_path, title)
    return results

def parse_arguments():
    """Create parser instance and parse command-line arguments."""
    parser = VersionInHelpParser(
        prog="NGS_analysis_collator", description="NGS analysis collator", version=__version__
    )
    # Important: include 'config' so --config (PEP path) is parsed into args.config_file
    parser = pypiper.add_pypiper_args(parser, groups=["pypiper", "looper", "config"])
    parser.add_argument("-n", "--name", help="Project name (record identifier).", type=str, required=False)
    parser.add_argument("--project-name", dest="project_name", help="Project name from looper.", type=str, required=False)
    parser.add_argument("-r", "--results", help="Output results subdirectory path.", type=str, required=False)
    parser.add_argument("--pipestat-config", dest="pipestat_config", default=None, help="Looper-generated pipestat config (YAML).")
    
    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()
    # Looper 2.1.1 + pipestat 0.13.x initialize project-level records with the
    # fallback identifier "project" unless project_name is present in the
    # generated pipestat config. To keep native `looper report` output coherent,
    # align the collator to that same project record instead of creating a
    # second synthetic project-level record.
    project_record_id = "project"

    # Align the collator working directory with the project-level record that
    # looper/pipestat create natively for project pipelines in this stack.
    outfolder = os.path.abspath(os.path.join(args.output_parent, project_record_id))
    os.makedirs(outfolder, exist_ok=True)

    # Resolve {record_identifier} so project results go to their own file.
    # Write inside outfolder to keep results/ clean.
    if args.pipestat_config and os.path.exists(args.pipestat_config):
        with open(args.pipestat_config) as _f:
            _cfg = yaml.safe_load(_f) or {}
        _rfp = _cfg.get("results_file_path", "")
        if _rfp and "{record_identifier}" in str(_rfp):
            _cfg["results_file_path"] = str(_rfp).replace("{record_identifier}", project_record_id)
            _resolved = os.path.join(outfolder, "pipestat_config.yaml")
            with open(_resolved, "w") as _f:
                yaml.dump(_cfg, _f)
            args.pipestat_config = _resolved

    # IMPORTANT: use the unified pipeline name to match your interfaces/schemas
    pm = pypiper.PipelineManager(
        name="NGS.analysis",
        outfolder=outfolder,
        pipestat_record_identifier=project_record_id,
        pipestat_pipeline_type="project",
        pipestat_config_file=args.pipestat_config,
        args=args,
        version=__version__
    )

    # ---- Run the R summarizer as before ----
    # args.config_file is populated by pypiper (since we added 'config' group)
    # args.results can be None; your interface passes -r {looper.output_dir}, adjust if needed
    r_cmd = (
        f"Rscript {tool_path('NGS.summarizer.R')} "
        f"{args.config_file} {outfolder} {args.results or args.output_parent}"
    )
    pm.run(r_cmd, "lock.max")

    # ---- Create a top-level reports/index.html redirect ----
    # pipestat places its HTML at reports/NGS.analysis/index.html; this redirect
    # makes results/reports/index.html work as a direct entry point.
    _reports_dir = os.path.join(args.output_parent, "reports")
    os.makedirs(_reports_dir, exist_ok=True)
    _redirect = os.path.join(_reports_dir, "index.html")
    with open(_redirect, "w") as _rf:
        _rf.write(
            '<!DOCTYPE html>\n'
            '<html>\n'
            '<head>\n'
            '<meta http-equiv="refresh" content="0; url=NGS.analysis/index.html">\n'
            '<title>NGS.analysis Report</title>\n'
            '</head>\n'
            '<body>\n'
            '<p>Redirecting to '
            '<a href="NGS.analysis/index.html">NGS.analysis/index.html</a></p>\n'
            '</body>\n'
            '</html>\n'
        )

    # Report project-level summary artifacts so the synthetic "summary" record
    # becomes a real record page in the HTML report.
    project_summary_dir = os.path.join(outfolder, "summary")
    project_results = collect_project_results(args.name or project_record_id, project_summary_dir)
    if project_results:
        try:
            pm.pipestat.report(record_identifier=project_record_id, values=project_results)
        except Exception as e:
            pm.warn(f"Pipestat project reporting skipped: {e}")

    pm.stop_pipeline()

if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Pipeline aborted.")
        sys.exit(1)
