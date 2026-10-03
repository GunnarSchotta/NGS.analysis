#!/usr/bin/env python3
"""Simulate PE reads from short fragments WITH adapter read-through (ART cannot simulate fragments shorter than the
read length). Usage: simulate_adapters.py genome.fa chrom n_fragments read_len out_prefix [protocol]
Fragment length uniform 30..150 bp from random positions on <chrom> (no N); Nextera (ATAC/CT) adapters;
0.5% substitution errors; quality 'I'. Read names encode the truth: <chrom>:<start0>:<len>:<i>."""
import random
import sys

import pysam

fa, chrom, n, rl, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
AD1, AD2 = "CTGTCTCTTATACACATCTCCGAGCCCACGAGACTAAGGCGAATCTCGTATGCCGTCTTCTGCTTG", \
           "CTGTCTCTTATACACATCTGACGCTGCCGACGATCTCGTGTAGATCTCGGTGGTCGCCGTATCATT"
rng = random.Random(42)
seq = pysam.FastaFile(fa).fetch(chrom).upper()
comp = str.maketrans("ACGTN", "TGCAN")


def err(s):
    return "".join(rng.choice("ACGT".replace(c, "")) if rng.random() < 0.005 and c in "ACGT" else c for c in s)


with open(out + "_R1.fastq", "w") as o1, open(out + "_R2.fastq", "w") as o2:
    i = 0
    while i < n:
        L = rng.randint(30, 150)
        s = rng.randint(0, len(seq) - L - 1)
        frag = seq[s:s + L]
        if "N" in frag:
            continue
        r1 = (frag + AD1)[:rl]
        r2 = (frag.translate(comp)[::-1] + AD2)[:rl]
        name = f"{chrom}:{s}:{L}:{i}"
        o1.write(f"@{name}/1\n{err(r1)}\n+\n{'I' * rl}\n")
        o2.write(f"@{name}/2\n{err(r2)}\n+\n{'I' * rl}\n")
        i += 1
