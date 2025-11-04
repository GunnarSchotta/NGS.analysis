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

def parse_arguments():
    """Create parser instance and parse command-line arguments."""
    parser = VersionInHelpParser(
        prog="NGS_analysis_collator", description="NGS analysis collator", version=__version__
    )
    # Important: include 'config' so --config (PEP path) is parsed into args.config_file
    parser = pypiper.add_pypiper_args(parser, groups=["pypiper", "looper", "config"])
    parser.add_argument("-n", "--name", help="Project name (record identifier).", type=str, required=False)
    parser.add_argument("-r", "--results", help="Output results subdirectory path.", type=str, required=False)
    #parser.add_argument("--pipestat-config", default=None, help="Looper-generated pipestat config (YAML).")
    
    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()

    # Where the collator writes its own summary artifacts
    outfolder = os.path.abspath(os.path.join(args.output_parent, "summary"))
    os.makedirs(outfolder, exist_ok=True)

    # IMPORTANT: use the unified pipeline name to match your interfaces/schemas
    pm = pypiper.PipelineManager(
        name="NGS.analysis",
        outfolder=outfolder,
        pipestat_record_identifier="summary",
        pipestat_pipeline_type="project",
        args=args,              # ensure pypiper sees --config and other CLI args
        version=__version__
    )

    # ---- Run the R summarizer as before ----
    # args.config_file is populated by pypiper (since we added 'config' group)
    # args.results can be None; your interface passes -r {looper.output_dir}, adjust if needed
    r_cmd = (
        f"Rscript {tool_path('NGS.summarizer.R')} "
        f"{args.config_file} {args.output_parent} {args.results or args.output_parent}"
    )
    pm.run(r_cmd, "lock.max")

    # ---- (Optional) Report project-level outputs to pipestat ----
    # Uncomment/add keys that you have declared in your *project* output schema
    # e.g., in NGS.analysis_output_schema.yaml:
    #   project_summary_dir: {type: string}
    #   collator_status:     {type: string}
    #
    # try:
    #     pm.pipestat.report(values={
    #         "project_summary_dir": outfolder,
    #         "collator_status": "completed"
    #     })
    # except Exception as e:
    #     pm.warn(f"Pipestat reporting skipped: {e}")

    pm.stop_pipeline()

if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Pipeline aborted.")
        sys.exit(1)
