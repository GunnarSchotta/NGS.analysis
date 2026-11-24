#!/usr/bin/env python3
"""Count non-N bases of a genome FASTA (reference for `macs_gsize` in genomes/<g>.yaml).

MACS3's shortcuts are the non-N length of the assembly (mm = 2,652,783,500 GRCm38, hs = 2,913,022,398 GRCh38,
MACS3/Utilities/Constants.py). The pipeline pins those numbers for mm10/hg38 and reuses them for mhaESC/hs1, whose
extra non-N sequence is mostly unmappable satellite (genomes/README.md); this script documents the actual counts.

Usage: count_nonN.py genome.fa [canonical_regex]
Prints total length, total non-N, and non-N over the canonical chromosomes (regex, default chr([0-9]+|X|Y)).
"""
import re
import sys

fa = sys.argv[1]
canon = re.compile(sys.argv[2] if len(sys.argv) > 2 else r"chr([0-9]+|X|Y)")
tot = nonn = canon_nonn = 0
name, is_canon = None, False
with open(fa) as f:
    for line in f:
        if line.startswith(">"):
            name = line[1:].split()[0]
            is_canon = bool(canon.fullmatch(name))
            continue
        s = line.strip()
        n = len(s) - s.count("N") - s.count("n")
        tot += len(s)
        nonn += n
        if is_canon:
            canon_nonn += n
print(f"{fa}\ttotal={tot}\tnonN={nonn}\tnonN_canonical={canon_nonn}")
