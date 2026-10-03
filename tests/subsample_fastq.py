#!/usr/bin/env python3
"""Reservoir-sample N reads (SE) or read pairs (PE, R1/R2 in lockstep) from gzipped FASTQ. Deterministic (seed).
Usage: subsample_fastq.py N SEED OUT_PREFIX R1.fastq.gz [R2.fastq.gz]
Writes OUT_PREFIX_R1.fastq.gz [and OUT_PREFIX_R2.fastq.gz]."""
import gzip, random, sys, itertools

def records(path):
    with gzip.open(path, "rt") as fh:
        while True:
            rec = list(itertools.islice(fh, 4))
            if len(rec) < 4:
                return
            yield "".join(rec)

n, seed, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
files = sys.argv[4:]
rng = random.Random(seed)
res = []
streams = zip(*[records(f) for f in files])
for i, rec in enumerate(streams):
    if i < n:
        res.append(rec)
    else:
        j = rng.randint(0, i)
        if j < n:
            res[j] = rec
for k in range(len(files)):
    with gzip.open(f"{out}_R{k + 1}.fastq.gz", "wt", compresslevel=4) as o:
        o.writelines(r[k] for r in res)
print(f"kept {len(res)} of {i + 1} records")
