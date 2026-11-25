#!/usr/bin/env python3
"""MAPQ threshold curves for the ART benchmark: how many fragments does each aligner keep, and how many of those
are at the true locus, for every MAPQ cutoff? (Decides whether STAR "unique" beats bowtie2 at an equal error rate.)

Usage: mapq_curve.py <truth.sam> <rmsk.bed> <out.tsv> <label>=<bam> [<label>=<bam> ...]

Origin family and true position as in evaluate.py (RepeatMasker element with the largest overlap of the true
fragment; correct = R1 primary on the true chromosome, start within 5 bp). Output, one row per
(label, family, class, mapq): n_family (all simulated fragments of the family), n_at_mapq (R1 primary alignments
with exactly this MAPQ), correct_at_mapq. Cumulative "kept at MAPQ >= t" and its precision are computed downstream.
"""
import collections
import os
import subprocess
import sys

import pysam

truth_sam, rmsk, out = sys.argv[1:4]
bams = dict(a.split("=", 1) for a in sys.argv[4:])

frag_bed = out + ".frag.bed"
with pysam.AlignmentFile(truth_sam, check_sq=False) as t, open(frag_bed, "w") as o:
    for r in t:
        r1 = r.is_read1 if r.is_paired else r.query_name.endswith("/1")
        if r1 and not r.is_unmapped:
            s = min(r.reference_start, r.next_reference_start) if r.is_paired else r.reference_start
            e = s + abs(r.template_length) if r.template_length else r.reference_end
            o.write(f"{r.reference_name}\t{s}\t{e}\t{r.query_name.split('/')[0]}\t{r.reference_start}\n")
best = {}
p = subprocess.Popen(["bedtools", "intersect", "-wo", "-a", frag_bed, "-b", rmsk], stdout=subprocess.PIPE, text=True)
for line in p.stdout:
    c = line.rstrip("\n").split("\t")
    ov, name = int(c[-1]), c[3]
    if name not in best or ov > best[name][0]:
        best[name] = (ov, c[8], c[9], c[0], int(c[4]))
if p.wait():
    sys.exit("bedtools intersect failed")
os.remove(frag_bed)
n_fam = collections.Counter(v[1:3] for v in best.values())

with open(out, "w") as fh:
    fh.write("label\tfamily\tclass\tmapq\tn_family\tn_at_mapq\tcorrect_at_mapq\n")
    for label, bam in bams.items():
        cnt = collections.defaultdict(lambda: [0, 0])          # (family, class, mapq) -> [n, correct]
        with pysam.AlignmentFile(bam) as b:
            for r in b.fetch(until_eof=True):
                if not r.is_read1 or r.is_secondary or r.is_supplementary or r.is_unmapped:
                    continue
                v = best.get(r.query_name.split("/")[0])
                if v is None:
                    continue
                k = (v[1], v[2], r.mapping_quality)
                cnt[k][0] += 1
                if r.reference_name == v[3] and abs(r.reference_start - v[4]) <= 5:
                    cnt[k][1] += 1
        for (fam, cls, q), (n, ok) in sorted(cnt.items()):
            fh.write(f"{label}\t{fam}\t{cls}\t{q}\t{n_fam[(fam, cls)]}\t{n}\t{ok}\n")
        print(f"{label}: done", flush=True)
