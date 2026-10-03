#!/usr/bin/env python3
"""Evaluate benchmark alignments against the simulated truth.
  ART sets   : evaluate.py art   <truth.sam> <rmsk.bed> <bam> <label> <out.tsv>
               truth = ART -sam alignment of R1 (true position); fragment origin family = RepeatMasker element with
               the largest overlap of the true fragment. Per family: n fragments, % aligned (primary R1), % aligned
               uniquely (label-specific MAPQ, encoded as label suffix ':q<mapq>'), TP rate of unique alignments
               (same chrom, |pos - true pos| <= 5).
  adapter set: evaluate.py adapt <bam> <label> <out.tsv> <R1.fastq>
               per true fragment length: % recovered = R1 proper pair, MAPQ>=q, correct start (+-2) and length (+-2).
"""
import collections
import os
import subprocess
import sys

import pandas as pd
import pysam


def mapq_of(label):
    return int(label.split(":q")[-1].lstrip("q")) if ":q" in label else 30


if sys.argv[1] == "art":
    truth_sam, rmsk, bam, label, out = sys.argv[2:7]
    q = mapq_of(label)
    # true fragment intervals from R1 truth (flag 64) -> BED, overlap with rmsk (largest overlap)
    frag_bed = out + ".frag.bed"
    with pysam.AlignmentFile(truth_sam, check_sq=False) as t, open(frag_bed, "w") as o:
        for r in t:
            r1 = r.is_read1 if r.is_paired else r.query_name.endswith("/1")
            if r1 and not r.is_unmapped:
                s = min(r.reference_start, r.next_reference_start) if r.is_paired else r.reference_start
                e = s + abs(r.template_length) if r.template_length else r.reference_end
                o.write(f"{r.reference_name}\t{s}\t{e}\t{r.query_name.split('/')[0]}\t{r.reference_start}\n")
    best = {}   # fragment -> (overlap, family, class, true chrom, true R1 start); streamed (whole output = tens of GB)
    p = subprocess.Popen(["bedtools", "intersect", "-wo", "-a", frag_bed, "-b", rmsk], stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        c = line.rstrip("\n").split("\t")
        ov, name = int(c[-1]), c[3]
        if name not in best or ov > best[name][0]:
            best[name] = (ov, c[8], c[9], c[0], int(c[4]))
    if p.wait():
        sys.exit("bedtools intersect failed")
    os.remove(frag_bed)
    stats = collections.defaultdict(lambda: [0, 0, 0, 0])   # n, aligned, unique, unique_correct
    for name, v in best.items():
        stats[v[1:3]][0] += 1
    with pysam.AlignmentFile(bam) as b:
        for r in b.fetch(until_eof=True):
            if not r.is_read1 or r.is_secondary or r.is_supplementary:
                continue
            name = r.query_name.split("/")[0]
            if name not in best:
                continue
            k = best[name][1:3]
            if r.is_unmapped:
                continue
            stats[k][1] += 1
            if r.mapping_quality >= q:
                stats[k][2] += 1
                tc, tp = best[name][3:5]
                if r.reference_name == tc and abs(r.reference_start - tp) <= 5:
                    stats[k][3] += 1
    rows = [{"label": label, "family": k[0], "class": k[1], "n": v[0], "aligned_pct": 100 * v[1] / v[0],
             "unique_pct": 100 * v[2] / v[0], "unique_TP_pct": (100 * v[3] / v[2]) if v[2] else None}
            for k, v in stats.items()]
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
else:
    bam, label, out, fq1 = sys.argv[2:6]
    n_by_len = collections.Counter(int(l.split(":")[2]) for l in open(fq1) if l.startswith("@"))
    q = mapq_of(label)
    tot, ok = collections.Counter(), collections.Counter()
    seen = set()
    with pysam.AlignmentFile(bam) as b:
        for r in b.fetch(until_eof=True):
            if not r.is_read1 or r.is_secondary or r.is_supplementary:
                continue
            chrom, s, L, i = r.query_name.split("/")[0].split(":")
            s, L = int(s), int(L)
            seen.add(i)
            if (r.is_proper_pair and not r.is_unmapped and r.mapping_quality >= q and r.reference_name == chrom
                    and abs(min(r.reference_start, r.next_reference_start) - s) <= 2 and abs(abs(r.template_length) - L) <= 2):
                ok[L] += 1
    d = pd.DataFrame({"label": label, "frag_len": sorted(n_by_len)})
    d["n"] = d.frag_len.map(n_by_len)
    d["recovered"] = d.frag_len.map(lambda L: ok.get(L, 0))
    d["recovered_pct"] = 100 * d.recovered / d.n
    d.to_csv(out, sep="\t", index=False)
