#!/usr/bin/env python3
"""Score alignments of reads simulated from a mutated 'sample genome' (mutate_genome.py) against the unmodified reference.

Usage: eval_snp.py <truth.sam> <events.tsv> <rmsk.bed> <out.tsv> <label>=<bam> [...]

- Truth positions (ART -sam, in mutated coordinates) are converted to reference coordinates with the events table.
- Fragments whose R1 starts inside a non-reference TE insertion are labelled family "nonref:<family>", class "nonref":
  they have no true locus in the reference, so every placement is a false locus (correct_at_mapq = 0).
- Other fragments: origin family = RepeatMasker element with the largest overlap of the true fragment (reference
  coordinates); correct = R1 primary on the true chromosome, start within 5 bp of the converted true start.
Output columns as mapq_curve.py: label, family, class, mapq, n_family, n_at_mapq, correct_at_mapq.
"""
import bisect
import collections
import os
import subprocess
import sys

import pysam

truth_sam, events_tsv, rmsk, out = sys.argv[1:5]
bams = dict(a.split("=", 1) for a in sys.argv[5:])

thr, dl, te = [], [], []          # cumulative shift thresholds (mutated coords), shift sizes, TE intervals
for line in list(open(events_tsv))[1:]:
    mpos, rpos, kind, length, fam = line.rstrip("\n").split("\t")
    mpos, length = int(mpos), int(length)
    if kind in ("INS", "TE"):
        thr.append(mpos + length); dl.append(length)
        if kind == "TE":
            te.append((mpos, mpos + length, fam))
    elif kind == "DEL":
        thr.append(mpos); dl.append(length)
order = sorted(range(len(thr)), key=thr.__getitem__)
thr = [thr[i] for i in order]
cum, s = [], 0
for i in order:
    s += dl[i]
    cum.append(s)
te.sort()
te_starts = [t[0] for t in te]


def to_ref(m):
    k = bisect.bisect_right(thr, m)
    return m - (cum[k - 1] if k else 0)


def te_at(m):
    k = bisect.bisect_right(te_starts, m) - 1
    return te[k][2] if k >= 0 and te[k][0] <= m < te[k][1] else None


frag_bed = out + ".frag.bed"
best = {}                          # name -> (overlap, family, class, true chrom, true ref R1 start)
with pysam.AlignmentFile(truth_sam, check_sq=False) as t, open(frag_bed, "w") as o:
    for r in t:
        r1 = r.is_read1 if r.is_paired else r.query_name.endswith("/1")
        if not r1 or r.is_unmapped:
            continue
        name = r.query_name.split("/")[0]
        fam = te_at(r.reference_start)
        if fam:
            best[name] = (10**9, "nonref:" + fam, "nonref", None, None)
            continue
        s = min(r.reference_start, r.next_reference_start) if r.is_paired else r.reference_start
        e = s + abs(r.template_length) if r.template_length else r.reference_end
        rs, re_ = to_ref(s), max(to_ref(e), to_ref(s) + 1)
        o.write(f"{r.reference_name}\t{rs}\t{re_}\t{name}\t{to_ref(r.reference_start)}\n")
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
        cnt = collections.defaultdict(lambda: [0, 0])
        with pysam.AlignmentFile(bam) as b:
            for r in b.fetch(until_eof=True):
                if not r.is_read1 or r.is_secondary or r.is_supplementary or r.is_unmapped:
                    continue
                v = best.get(r.query_name.split("/")[0])
                if v is None:
                    continue
                k = (v[1], v[2], r.mapping_quality)
                cnt[k][0] += 1
                if v[3] is not None and r.reference_name == v[3] and abs(r.reference_start - v[4]) <= 5:
                    cnt[k][1] += 1
        for (fam, cls, q), (n, ok) in sorted(cnt.items()):
            fh.write(f"{label}\t{fam}\t{cls}\t{q}\t{n_fam[(fam, cls)]}\t{n}\t{ok}\n")
        print(f"{label}: done", flush=True)
