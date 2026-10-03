#!/usr/bin/env python3
"""Preflight check for an NGS.analysis v3 project -- run in the project folder BEFORE `looper run`.

Usage:  python check_project.py [analysis.configuration.yaml]
Checks
  - sample names: unique; no name is a prefix of another at an '_' boundary (pipestat finds status flags with
    the glob {pipeline}_{sample}_*.flag, so 'A' would also match the flags of 'A_rep2'); not 'project'
  - required columns; protocol / read_type values; FASTQ files exist
  - genome resources: genomes/<genome>.yaml exists and the files it names exist
  - controls (column 'control', optional) refer to existing samples with the same genome and read type
Exit status 1 if any ERROR was found.
"""
import os
import sys

import peppy
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
PROTOCOLS = {"CHIP", "CT", "CR", "RNA", "ATAC"}
GENOME_FILES = ["fai", "bowtie2_index", "star_genome_index", "star_rna_index", "repeats_saf", "repeats_safid",
                "tss_bed", "blacklist"]


def main(cfg="analysis.configuration.yaml"):
    errors, warnings = [], []
    prj = peppy.Project(cfg)
    names = [s.sample_name for s in prj.samples]
    seen = set()
    for n in names:
        if n in seen:
            errors.append(f"duplicate sample_name '{n}'")
        seen.add(n)
        if n == "project":
            errors.append("sample_name 'project' is reserved (project-level record)")
    for a in names:
        for b in names:
            if a != b and b.startswith(a + "_"):
                errors.append(f"sample_name '{a}' is a prefix of '{b}' at an '_' boundary: pipestat flag lookup "
                              f"'{{pipeline}}_{a}_*.flag' would also match '{b}'. Rename one of them.")
    mode = prj.config.get("pipeline_mode")
    if mode not in ("genes", "repeats"):
        errors.append(f"pipeline_mode must be genes or repeats, found {mode!r}")
    genomes_checked = {}
    for s in prj.samples:
        n = s.sample_name
        for col in ["protocol", "read_type", "genome", "read1"]:
            if not s.get(col):
                errors.append(f"{n}: missing '{col}'")
        if s.get("protocol") not in PROTOCOLS:
            errors.append(f"{n}: protocol {s.get('protocol')!r} not in {sorted(PROTOCOLS)}")
        if s.get("read_type") not in ("paired", "single"):
            errors.append(f"{n}: read_type must be paired or single")
        reads = [s.get("read1")] + ([s.get("read2")] if s.get("read_type") == "paired" else [])
        for r in reads:
            for f in (r if isinstance(r, list) else [r]):
                if f and not os.path.isfile(f):
                    errors.append(f"{n}: FASTQ not found: {f}")
        g = s.get("genome")
        gcfg = s.get("genome_config") or os.path.join(HERE, "genomes", f"{g}.yaml")
        if g and gcfg not in genomes_checked:
            if not os.path.exists(gcfg):
                errors.append(f"{n}: no genome resource file {gcfg}")
                genomes_checked[gcfg] = None
            else:
                y = yaml.safe_load(open(gcfg))
                genomes_checked[gcfg] = y
                for k in GENOME_FILES:
                    v = y.get(k)
                    if v is None:
                        if k == "blacklist":
                            warnings.append(f"genome {g}: no blacklist")
                        continue
                    if not (os.path.exists(v) or any(os.path.exists(v + e) for e in (".1.bt2", ".grp", "/SA"))):
                        (errors if k in ("fai", "bowtie2_index") else warnings).append(f"genome {g}: {k} missing: {v}")
    by_name = {s.sample_name: s for s in prj.samples}
    for s in prj.samples:
        c = s.get("control")
        if c and str(c).lower() not in ("none", "na", ""):
            if c not in by_name:
                errors.append(f"{s.sample_name}: control '{c}' is not a sample of this project")
            else:
                for col in ("genome", "read_type"):
                    if by_name[c].get(col) != s.get(col):
                        errors.append(f"{s.sample_name}: control '{c}' differs in {col}")
    for w in warnings:
        print("WARNING:", w)
    for e in errors:
        print("ERROR:", e)
    print(f"{len(names)} samples, {len(errors)} errors, {len(warnings)} warnings")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:2]))
