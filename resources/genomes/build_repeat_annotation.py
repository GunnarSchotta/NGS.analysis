#!/usr/bin/env python3
"""Static per-element repeat annotation for NGS.analysis (one file per genome; genomes/<g>.yaml: repeats_annotation).

For every element of the genome's element SAF (repeats_safid; the IDs the pipeline counts):
  element_id family class repfamily chr start end strand length div_pct n_overlap overlap_families blacklist
- div_pct: RepeatMasker divergence from the family consensus in % (UCSC rmsk milliDiv / 10; RepeatMasker .out
  "perc div."; RepeatMasker GFF score). NA if the element is not in the divergence source (annotation versions differ).
- n_overlap / overlap_max_frac / overlap_families: other entries of the same SAF overlapping the element by >= 10 bp
  (nested or overlapping annotation; < 10 bp is RepeatMasker boundary jitter), the largest fraction of the element
  covered by one of them, and up to 5 of their family names (comma-separated, sorted, "." if none). With
  featureCounts -O a fragment in such an overlap is counted for every entry it touches.
- blacklist: 1 if the element overlaps the genome blacklist (0 if there is none).
The pipeline derives its flags (young, nested, blacklist) from these columns; thresholds live in NGS.analysis.py.

Usage: build_repeat_annotation.py <genome> <format: ucsc|out|gff> <source>[,<source>...] [<out.tsv.gz>]
  genome: name of genomes/<genome>.yaml (repeats_safid, blacklist are read from it)
  ucsc: UCSC rmsk.txt (bin swScore milliDiv ... genoName genoStart genoEnd ... repName repClass repFamily ...)
  out:  RepeatMasker .out(.gz);  gff: RepeatMasker GFF(.gz) with Target "Motif:<name>" (mhaESC: Chr01 -> chr1)
Run on the HPC (bedtools from ngs.v3 on PATH), about 10 min and 16 GB for mm10.
"""
import collections
import gzip
import os
import re
import subprocess
import sys
import tempfile

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
MIN_OV = 10
genome, fmt, sources = sys.argv[1], sys.argv[2], sys.argv[3].split(",")
cfg = yaml.safe_load(open(os.path.join(HERE, "..", "..", "genomes", genome + ".yaml")))
safid = cfg["repeats_safid"]
blacklist = cfg.get("blacklist")
out = sys.argv[4] if len(sys.argv) > 4 else os.path.join(os.path.dirname(safid), f"repeat_elements.{genome}.tsv.gz")


def zopen(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def rename_mha(c):          # as prep_t2t.py: Chr01 -> chr1, ChrX -> chrX
    m = re.match(r"Chr0*(\w+)$", c)
    return "chr" + m.group(1) if m else c


# ---- divergence source: (chr, saf_start, end) -> (name, class, repfamily, div_pct)
src = {}
for path in sources:
    for line in zopen(path):
        if fmt == "ucsc":
            f = line.rstrip("\n").split("\t")
            if f[0].startswith("#"):
                continue
            src[(f[5], int(f[6]), int(f[7]))] = (f[10], f[11], f[12], int(f[2]) / 10)
        elif fmt == "out":
            f = line.split()
            if len(f) < 15 or not f[0].isdigit():
                continue
            cls, _, fam = f[10].partition("/")
            src[(f[4], int(f[5]), int(f[6]))] = (f[9], cls, fam or cls, float(f[1]))
        elif fmt == "gff":
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            m = re.search(r'Motif:([^"]+)"', c[8])
            if not m:
                continue
            cls, _, fam = c[2].partition("/")
            src[(rename_mha(c[0]), int(c[3]), int(c[4]))] = (m.group(1), cls, fam or cls,
                                                            float(c[5]) if c[5] != "." else None)
        else:
            sys.exit("format must be ucsc, out or gff")
print(f"divergence source: {len(src):,} entries", flush=True)

# ---- elements of the SAF
rows = []
tmp = tempfile.mkdtemp(dir=os.path.dirname(out))
bed = os.path.join(tmp, "el.bed")
with zopen(safid) as fh, open(bed, "w") as b:
    for line in fh:
        eid, ch, st, en, strand = line.rstrip("\n").split("\t")[:5]
        if eid == "GeneID":
            continue
        st, en = int(st), int(en)
        k = eid.rfind(ch + ":")
        fam = eid[:k - 1] if k > 0 else eid
        s = src.get((ch, st, en))
        div = None if s is None else s[3]
        cls, rfam = (s[1], s[2]) if s else ("NA", "NA")
        rows.append([eid, fam, cls, rfam, ch, st, en, strand, en - st + 1, div])
        b.write(f"{ch}\t{max(0, st - 1)}\t{en}\t{len(rows) - 1}\t{fam}\n")
del src
matched = sum(r[9] is not None for r in rows)
print(f"elements: {len(rows):,}; with divergence: {matched:,} ({100 * matched / len(rows):.1f} %)", flush=True)

subprocess.run(f"sort -k1,1 -k2,2n -S 4G -T {tmp} {bed} -o {bed}", shell=True, check=True)
nov = [0] * len(rows)
cov = [0.0] * len(rows)
ofam = collections.defaultdict(set)
p = subprocess.Popen(["bedtools", "intersect", "-sorted", "-wo", "-a", bed, "-b", bed],
                     stdout=subprocess.PIPE, text=True)
for line in p.stdout:
    c = line.rstrip("\n").split("\t")
    i, j, ov = int(c[3]), int(c[8]), int(c[10])
    if i != j and ov >= MIN_OV:
        nov[i] += 1
        cov[i] = max(cov[i], ov / rows[i][8])
        if len(ofam[i]) < 5:
            ofam[i].add(c[9])
if p.wait():
    sys.exit("bedtools intersect (overlaps) failed")
bl = set()
if blacklist and os.path.exists(blacklist):
    p = subprocess.Popen(["bedtools", "intersect", "-u", "-a", bed, "-b", blacklist], stdout=subprocess.PIPE, text=True)
    bl = {int(line.split("\t")[3]) for line in p.stdout}
    if p.wait():
        sys.exit("bedtools intersect (blacklist) failed")
subprocess.run(["rm", "-rf", tmp])

with gzip.open(out + ".tmp", "wt") as o:
    o.write(f"# build_repeat_annotation.py {genome}: elements {safid}; divergence {fmt} {','.join(sources)}; "
            f"blacklist {blacklist or 'none'}\n")
    o.write("element_id\tfamily\tclass\trepfamily\tchr\tstart\tend\tstrand\tlength\tdiv_pct\tn_overlap\t"
            "overlap_max_frac\toverlap_families\tblacklist\n")
    for i, r in enumerate(rows):
        o.write("\t".join(map(str, r[:9])) + f"\t{'NA' if r[9] is None else round(r[9], 1)}\t{nov[i]}\t"
                f"{round(cov[i], 3)}\t{','.join(sorted(ofam[i])) or '.'}\t{int(i in bl)}\n")
os.replace(out + ".tmp", out)
print(f"written {out}; nested/overlapping: {sum(n > 0 for n in nov):,}; blacklist: {len(bl):,}")
