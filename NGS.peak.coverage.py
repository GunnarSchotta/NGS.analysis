#!/usr/bin/env python3
"""
NGS Peak Coverage
"""

__author__ = ["Gunnar Schotta"]
__email__ = "gunnar.schotta@bmc.med.lmu.de"
__version__ = "0.0.1"

from argparse import ArgumentParser
import os
import sys
import pypiper
from ubiquerg import VersionInHelpParser

#run with -C config file -O parent.output -m coverage|density -p peak.file

def parse_arguments():
    """
    Creat parser instance and parse command-line arguments passed to the pipeline

    :return argparse.Namespace: parsed arguments namespace
    """
    parser = VersionInHelpParser(prog="NGS_analysis_collator",
        description='NGS Peak Coverage' , version=__version__)
    parser = pypiper.add_pypiper_args(parser, groups=['pypiper', 'looper'])
    parser.add_argument("-m", "--mode",
                        help="Run mode (coverage|density).", type=str)
    parser.add_argument("-p", "--peak.file",
                        help="Peak file in BED format.", type=str)
    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()
    outfolder = os.path.abspath(os.path.join(args.output_parent, "peak.coverage"))

    pm = pypiper.PipelineManager(name="NGS_Peak_Coverage", outfolder=outfolder,
                                 args=args, version=__version__)

    cmd = 
#    cmd = (f"Rscript {tool_path('NGS.summarizer.R')} "
#           f"{args.config_file} {args.output_parent} {args.results}")

#    pm.run(cmd, "lock.max")
    pm.stop_pipeline()


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Pipeline aborted.")
        sys.exit(1)
