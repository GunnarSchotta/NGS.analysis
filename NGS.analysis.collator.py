#!/usr/bin/env python3
"""
NGS Analysis Collator - project-level pipeline
"""

__author__ = ["Gunnar Schotta"]
__email__ = "gunnar.schotta@bmc.med.lmu.de"
__version__ = "1.1.2"

from argparse import ArgumentParser
import os
import sys
import pypiper
from ubiquerg import VersionInHelpParser

def tool_path(tool_name):
    """
    Return the path to a tool used by this pipeline.

    :param str tool_name: name of the tool (e.g., a script file_name)
    :return str: real, absolute path to tool (expansion and symlink resolution)
    """

    return os.path.join(os.path.dirname(__file__), tool_name)


def parse_arguments():
    """
    Creat parser instance and parse command-line arguments passed to the pipeline

    :return argparse.Namespace: parsed arguments namespace
    """
    parser = VersionInHelpParser(prog="NGS_analysis_collator",
        description='NGS analysis collator' , version=__version__)
    parser = pypiper.add_pypiper_args(parser, groups=['pypiper', 'looper'])
    parser.add_argument("-n", "--name",
                        help="Name of the project to use.", type=str)
    parser.add_argument("-r", "--results",
                        help="Output results sub directory path.", type=str)
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()
    outfolder = os.path.abspath(os.path.join(args.output_parent, "summary"))

    pm = pypiper.PipelineManager(name="NGS_analysis_collator", outfolder=outfolder,
                                 args=args, version=__version__)

    cmd = (f"Rscript {tool_path('NGS.summarizer.R')} "
           f"{args.config_file} {args.output_parent} {args.results}")

    pm.run(cmd, "lock.max")
    pm.stop_pipeline()


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Pipeline aborted.")
        sys.exit(1)
