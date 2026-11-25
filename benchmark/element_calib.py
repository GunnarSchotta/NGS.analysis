#!/usr/bin/env python3
"""Calibrate per-element uncertainty flags on the ART / SNP simulation (bowtie2 BAMs from run_snp_benchmark.sh).

Question: for an element that RECEIVES fragments at a given MAPQ cutoff, how many of them are not from that
element (misplaced reference reads or reads from non-reference TE copies), depending on the element's divergence
from its consensus (RepeatMasker milliDiv)? And which fragment-MAPQ definition to use for pairs: R1 only, min of
the mates (both mates pass; = NGS.analysis filtered BAM), or max (one anchoring mate is enough)?

Usage: element_calib.py <truth.sam> <events.tsv|-> <rmsk_div.bed> <bam> <out_prefix>
  rmsk_div.bed: chrom start end name class family milliDiv (from UCSC rmsk.txt)
  bam: bowtie2 output in read-name order (mates adjacent), as written by run_snp_benchmark.sh

Fragments: concordant pairs only (flag 2, both primary), like the pipeline. A fragment is correct when R1 is on the
true chromosome within 5 bp of its true start (eval_snp.py rule); fragments from non-reference insertions are
always wrong. Landing element = RepeatMasker entry with the largest overlap of the placed fragment.

Outputs:
  <out>.landing.tsv  sem, mapq_ge, group, div_bin, n_placed, n_wrong, n_wrong_nonref
  <out>.origin.tsv   sem, mapq_ge, group, div_bin, n_sim, n_kept, n_correct   (by the true element)
"""
import bisect
import collections
import os
import subprocess
import sys

import pysam

truth_sam, events_tsv, rmsk, bam, out = sys.argv[1:6]
CUTS = [0, 2, 5, 10, 20, 30]
BINS = [0, 1, 2, 3, 5, 8, 12, 20, 1000]          # % divergence


def div_bin(md):
    d = md / 10.0
    k = bisect.bisect_right(BINS, d) - 1
    return f"{BINS[k]}-{BINS[k + 1]}" if BINS[k + 1] < 1000 else f">={BINS[k]}"


def group(cls, fam, name):
    if cls == "LINE" and fam == "L1":
        return "L1"
    if cls == "LTR":
        return {"ERVK": "ERVK", "ERV1": "ERV1", "ERVL": "ERVL", "ERVL-MaLR": "ERVL-MaLR"}.get(fam, "LTR other")
    if cls == "SINE":
        return "SINE"
    if cls in ("Satellite", "Simple_repeat", "Low_complexity"):
        return "simple/satellite"
    return "other"


# ---- coordinate conversion (as eval_snp.py)
thr, dl, te = [], [], []
if events_tsv != "-":
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


# ---- truth: chrom, ref R1 start, origin element (group, div bin) or nonref
truth = {}                                   # name -> [chrom|None, r1start, origin_key]
tbed = out + ".truth.bed"
with pysam.AlignmentFile(truth_sam, check_sq=False) as t, open(tbed, "w") as o:
    for r in t:
        r1 = r.is_read1 if r.is_paired else r.query_name.endswith("/1")
        if not r1 or r.is_unmapped:
            continue
        name = r.query_name.split("/")[0]
        if te_at(r.reference_start):
            truth[name] = [None, None, ("nonref", "nonref")]
            continue
        s0 = min(r.reference_start, r.next_reference_start)
        e0 = s0 + abs(r.template_length) if r.template_length else r.reference_end
        rs, re_ = to_ref(s0), max(to_ref(e0), to_ref(s0) + 1)
        truth[name] = [r.reference_name, to_ref(r.reference_start), ("none", "none")]
        o.write(f"{r.reference_name}\t{rs}\t{re_}\t{name}\n")


def best_overlap(bed, cb):
    best = {}
    p = subprocess.Popen(["bedtools", "intersect", "-wo", "-a", bed, "-b", rmsk], stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        c = line.rstrip("\n").split("\t")
        ov, name = int(c[-1]), c[3]
        if name not in best or ov > best[name][0]:
            best[name] = (ov, (group(c[-4], c[-3], c[-5]), div_bin(int(c[-2]))))
    if p.wait():
        sys.exit("bedtools intersect failed")
    for name, (ov, key) in best.items():
        cb(name, key)


def set_origin(name, key):
    truth[name][2] = key


best_overlap(tbed, set_origin)
os.remove(tbed)

# ---- placed fragments
SEMS = ("r1", "min", "max")
placed = {}                                  # name -> (q_r1, q_min, q_max, correct)
pbed = out + ".placed.bed"
with pysam.AlignmentFile(bam) as b, open(pbed, "w") as o:
    pend = {}
    for r in b.fetch(until_eof=True):
        if r.is_secondary or r.is_supplementary or r.is_unmapped or not r.is_proper_pair:
            continue
        name = r.query_name.split("/")[0]
        if name not in pend:
            pend[name] = r
            continue
        m = pend.pop(name)
        r1, r2 = (r, m) if r.is_read1 else (m, r)
        tv = truth.get(name)
        if tv is None:
            continue
        ok = tv[0] is not None and r1.reference_name == tv[0] and abs(r1.reference_start - tv[1]) <= 5
        q1, q2 = r1.mapping_quality, r2.mapping_quality
        placed[name] = (q1, min(q1, q2), max(q1, q2), ok)
        s0 = min(r1.reference_start, r2.reference_start)
        e0 = max(r1.reference_end, r2.reference_end)
        o.write(f"{r1.reference_name}\t{s0}\t{e0}\t{name}\n")
landing = {}


def set_landing(name, key):
    landing[name] = key


best_overlap(pbed, set_landing)
os.remove(pbed)

# ---- tables
lt = collections.defaultdict(lambda: [0, 0, 0])
for name, (q1, qmin, qmax, ok) in placed.items():
    key = landing.get(name, ("none", "none"))
    nonref = truth[name][2][0] == "nonref"
    for sem, q in zip(SEMS, (q1, qmin, qmax)):
        for c in CUTS:
            if q >= c:
                v = lt[(sem, c) + key]
                v[0] += 1
                v[1] += not ok
                v[2] += nonref
with open(out + ".landing.tsv", "w") as fh:
    fh.write("sem\tmapq_ge\tgroup\tdiv_bin\tn_placed\tn_wrong\tn_wrong_nonref\n")
    for k, v in sorted(lt.items()):
        fh.write("\t".join(map(str, k + tuple(v))) + "\n")

ot = collections.defaultdict(lambda: [0, 0, 0])
nsim = collections.Counter(v[2] for v in truth.values())
for name, tv in truth.items():
    p = placed.get(name)
    if p is None:
        continue
    for sem, q in zip(SEMS, p[:3]):
        for c in CUTS:
            if q >= c:
                v = ot[(sem, c) + tv[2]]
                v[0] += 1
                v[1] += p[3]
with open(out + ".origin.tsv", "w") as fh:
    fh.write("sem\tmapq_ge\tgroup\tdiv_bin\tn_sim\tn_kept\tn_correct\n")
    for k, v in sorted(ot.items()):
        fh.write("\t".join(map(str, k + (nsim[k[2:]],) + tuple(v[:2]))) + "\n")
print("done", len(truth), "simulated,", len(placed), "concordant pairs placed")
