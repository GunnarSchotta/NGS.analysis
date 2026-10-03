#!/usr/bin/env python3
"""Static check: every result key reported by NGS.analysis.py (report / report_object / dict keys of ngs_qc
helpers) exists in pipestat_results_schema.yaml. pipestat raises ColumnNotFoundError otherwise."""
import os
import re
import sys

import yaml

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
schema = yaml.safe_load(open(os.path.join(HERE, "pipestat_results_schema.yaml")))["samples"]
code = open(os.path.join(HERE, "NGS.analysis.py")).read()
qc = open(os.path.join(HERE, "ngs_qc.py")).read()
keys = set(re.findall(r'report\("([^"]+)"', code)) | set(re.findall(r'report_object\("([^"]+)"', code))
keys |= set(re.findall(r'fastqc\([^,]+, "([^"]+)"', code))
keys |= {"BigWig_filtered_nuc", "BigWig_filtered_subnuc"}                     # report_object("BigWig_filtered_" + part)
keys |= set(re.findall(r'"(Insert_[A-Za-z0-9_]+)":', qc)) | {"NRF", "PBC1", "PBC2"}
keys.discard("BigWig_filtered_")
missing = sorted(k for k in keys if k not in schema)
print(f"{len(keys)} reported keys checked; missing in schema: {missing}")
sys.exit(1 if missing else 0)
