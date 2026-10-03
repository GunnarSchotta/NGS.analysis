#!/usr/bin/env python
# pyTssEnrichment.py  (NGS.analysis v3)
#
# Aggregate Tn5 insertion profile around TSSs (single-base BED positions with strand).
# Derived from PEPATAC's pyTssEnrichment.py (J. Buenrostro, V. Reuter, J. Smith); rewritten for v3:
#   - insertion sites are computed in GENOME coordinates for every TSS (the PEPATAC version placed the
#     far fragment end of minus-strand TSSs at pos-5-ilen, i.e. on the wrong side); orientation of minus-strand
#     TSSs is handled once, when the profile is filled
#   - paired-end: each fragment once (from its leftmost, forward-strand mate): ends pos+4 and pos+TLEN-5
#   - single-end (TLEN 0): the read's 5' insertion only (forward: pos+4; reverse: reference_end-5), both strands
#   - MAPQ >= 30, no secondary/supplementary/duplicate/QC-fail records
#   - robust chunking when there are fewer TSSs than threads
# Output (-v -z): one line per position from -u to +d (u+d values) = summed insertion counts.
# The TSS score is computed by the pipeline (ngs_qc.tss_score).

import sys
from multiprocessing import Pool
from optparse import OptionParser

import numpy as np
import pysam

opts = OptionParser(usage="usage: %prog -a reads.bam -b tss.bed -o out.txt [options]")
opts.add_option("-a", help="coordinate-sorted, indexed BAM")
opts.add_option("-b", help="BED of TSS positions (strand in column -s)")
opts.add_option("-o", help="output file")
opts.add_option("-u", default="2000", help="bases upstream (default 2000)")
opts.add_option("-d", default="2000", help="bases downstream (default 2000)")
opts.add_option("-p", default="ends", help="kept for compatibility; only 'ends' is supported")
opts.add_option("-c", default="8", help="threads")
opts.add_option("-s", default="6", help="strand column (1-based), default 6")
opts.add_option("-q", default="30", help="minimum MAPQ (default 30)")
opts.add_option("-z", action="store_true", default=False, help="plain-text output (compatibility)")
opts.add_option("-v", action="store_true", default=False, help="profile output (compatibility)")
options, _ = opts.parse_args()
if not (options.a and options.b and options.o):
    opts.print_help()
    sys.exit(1)

UP, DOWN, MINQ, SCOL = int(options.u), int(options.d), int(options.q), int(options.s) - 1
COLS = UP + DOWN


def load_tss(path):
    out = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 3 or line.startswith(("#", "track")):
            continue
        strand = f[SCOL] if len(f) > SCOL else "+"
        out.append((f[0], (int(f[1]) + int(f[2])) // 2, strand))
    return out


def insertion_sites(r):
    """Tn5 insertion positions (genome coordinates) contributed by this record."""
    if r.is_paired:
        if r.is_reverse or not r.is_proper_pair or r.template_length <= 0:
            return ()
        return (r.reference_start + 4, r.reference_start + r.template_length - 5)
    return (r.reference_end - 5,) if r.is_reverse else (r.reference_start + 4,)


def profile_chunk(tss_chunk):
    prof = np.zeros(COLS)
    bam = pysam.AlignmentFile(options.a, "rb")
    refs = set(bam.references)
    for chrom, center, strand in tss_chunk:
        if chrom not in refs:
            continue
        s_int, e_int = (center - DOWN, center + UP) if strand == "-" else (center - UP, center + DOWN)
        for r in bam.fetch(chrom, max(0, s_int - 1000), e_int + 1000):
            if (r.mapping_quality < MINQ or r.is_secondary or r.is_supplementary or r.is_duplicate
                    or r.is_qcfail or r.is_unmapped):
                continue
            for pos in insertion_sites(r):
                if s_int <= pos < e_int:
                    b = pos - s_int
                    prof[COLS - 1 - b if strand == "-" else b] += 1
    bam.close()
    return prof


if __name__ == "__main__":
    tss = load_tss(options.b)
    n = max(1, int(options.c))
    chunks = [tss[i::n] for i in range(n)]
    chunks = [c for c in chunks if c]
    with Pool(processes=len(chunks) or 1) as pool:
        total = np.sum(pool.map(profile_chunk, chunks), axis=0) if chunks else np.zeros(COLS)
    np.savetxt(options.o, total, fmt="%d")
