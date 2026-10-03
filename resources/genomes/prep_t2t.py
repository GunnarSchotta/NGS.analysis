#!/usr/bin/env python3
"""Prepare T2T genome resources for NGS.analysis v3 (run by build_genome.sh, step 'prep').

Usage: prep_t2t.py {hs1|mhaESC} OUTDIR
Writes (layout as for mm10/hg38):
  Sequence/WholeGenomeFasta/genome.fa (+ .fai, chrom.sizes)
  Annotation/Genes/genes.gtf, <g>_TSS.bed, rmsk.<g>.SAF, rmsk.ids.<g>.SAF, rmsk.<g>.bed
  PREP_REPORT.txt (what was done, incl. Y-PAR check for mhaESC)

hs1    : T2T-CHM13v2.0 analysis set (PAR-masked chrY, rCRS chrM); CAT/Liftoff GENCODE v35 annotation;
         UCSC hs1 RepeatMasker (.out). Element IDs: repName_chr:start-end (hg38 convention).
mhaESC : mhaESC v1.1 + mT2T-Y v1.1 (release 2026-01-29; Chr01 -> chr1 ...); combined gene annotation
         mhaESC_v1.1_with_mT2T-Y_v1.1.260129.gff3 (gene_name attribute); RepeatMasker GFF = mouse.241018.repeats.gff
         (chr1-19, X, M) + mT2T-Y_v1.1.repeats.gff (chrY). Element IDs:
         repName.chr:start-end (mm10 convention). chrY pseudoautosomal region: detected by k-mer identity with the
         distal chrX and hard-masked (N) on chrY if present, so that PAR reads are not multimappers.
"""
import gzip
import os
import random
import re
import subprocess
import sys

G, OUT = sys.argv[1], sys.argv[2]
SRC = "/store24/project24/becgsc_001/genomes"
GFFREAD = "/store24/project24/becgsc_001/micromamba/envs/ngs.v3/bin/gffread"
SAMTOOLS = "/store24/project24/becgsc_001/micromamba/envs/ngs.v3/bin/samtools"
seqdir = os.path.join(OUT, "Sequence", "WholeGenomeFasta")
anndir = os.path.join(OUT, "Annotation", "Genes")
os.makedirs(seqdir, exist_ok=True)
os.makedirs(anndir, exist_ok=True)
fa = os.path.join(seqdir, "genome.fa")
report = []


def zopen(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def rename_mha(name):
    n = re.sub(r"^Chr", "", name)
    n = re.sub(r"^0+", "", n)
    return "chr" + n


def read_fasta(path, rename=lambda x: x):
    seqs, name, buf = {}, None, []
    for line in zopen(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = rename(line[1:].split()[0]), []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def write_fasta(seqs, path, width=60):
    with open(path, "w") as o:
        for n, s in seqs.items():
            o.write(">" + n + "\n")
            for i in range(0, len(s), width):
                o.write(s[i:i + width] + "\n")


def find_par(seqs, k=50, win=10000, region=3_000_000, frac=0.8, samples=40):
    """Return (start, end) of a chrY end region identical to the distal chrX (PAR), or None."""
    x, y = seqs["chrX"].upper(), seqs["chrY"].upper()
    rng = random.Random(1)
    xkm = set()
    for xs in (x[-region:], x[:region]):
        xkm |= {xs[i:i + k] for i in range(0, len(xs) - k, 1)}

    def is_par(ws):
        ks = [ws[j:j + k] for j in (rng.randrange(0, len(ws) - k) for _ in range(samples))]
        ks = [s for s in ks if "N" not in s]
        return len(ks) >= samples // 2 and sum(s in xkm for s in ks) / len(ks) >= frac

    hits = []
    for end_side in ("q", "p"):
        pos = len(y) if end_side == "q" else 0
        while True:
            ws = y[pos - win:pos] if end_side == "q" else y[pos:pos + win]
            if len(ws) < win or not is_par(ws):
                break
            pos = pos - win if end_side == "q" else pos + win
        if end_side == "q" and pos < len(y):
            hits.append((pos, len(y)))
        if end_side == "p" and pos > 0:
            hits.append((0, pos))
    # a chromosome end that matches chrX only through its telomere repeat (TTAGGG)n is not a PAR
    telo = lambda s: 6 * (s.count("TTAGGG") + s.count("CCCTAA")) / max(len(s), 1)
    return [h for h in hits if telo(y[h[0]:h[1]]) < 0.5], [h for h in hits if telo(y[h[0]:h[1]]) >= 0.5]


# ----------------------------------------------------------------------------- genome FASTA
if G == "hs1":
    seqs = read_fasta(os.path.join(SRC, "hs1/source/chm13v2.0_maskedY_rCRS.fa.gz"))
    report.append("FASTA: chm13v2.0_maskedY_rCRS (chrY PARs hard-masked by the T2T consortium; chrM = rCRS)")
else:
    seqs = read_fasta(os.path.join(SRC, "mhaESC/source/release_Y1.1/mhaESC_v1.1_with_mT2T-Y_v1.1.251107.fasta.gz"), rename_mha)
    report.append("FASTA: mhaESC_v1.1_with_mT2T-Y_v1.1.251107 (Chr01 -> chr1, ChrX -> chrX, ChrY -> chrY, ChrM -> chrM)")
    par, telo_only = find_par(seqs)
    if telo_only:
        report.append(f"chrY end region(s) matching chrX only via telomere repeat, not masked: {telo_only}")
    if par:
        y = list(seqs["chrY"])
        for s, e in par:
            y[s:e] = "N" * (e - s)
        seqs["chrY"] = "".join(y)
        report.append(f"chrY PAR detected (k-mer identity with distal chrX) and hard-masked: {par}")
    else:
        report.append("chrY PAR check: no chrY end region identical to chrX found; chrY not masked")
write_fasta(seqs, fa)
subprocess.check_call([SAMTOOLS, "faidx", fa])
with open(fa + ".fai") as f, open(os.path.join(OUT, "chrom.sizes"), "w") as o:
    for line in f:
        o.write("\t".join(line.split("\t")[:2]) + "\n")
report.append("chromosomes: " + " ".join(seqs))
del seqs

# ----------------------------------------------------------------------------- gene annotation -> GTF, TSS
gtf = os.path.join(anndir, "genes.gtf")
tss_bed = os.path.join(anndir, f"{G}_TSS.bed")
if G == "hs1":
    gff = os.path.join(SRC, "hs1/source/chm13.draft_v2.0.gene_annotation.gff3.gz")
    tmp_gff = os.path.join(anndir, "tmp.gff3")
    with zopen(gff) as f, open(tmp_gff, "w") as o:
        for line in f:
            if line.startswith("#") or line.split("\t")[2] in ("gene", "transcript", "exon", "CDS"):
                o.write(line)
    name_of = None
else:
    gff = os.path.join(SRC, "mhaESC/source/release_Y1.1/mhaESC_v1.1_with_mT2T-Y_v1.1.260129.gff3.gz")
    tmp_gff = os.path.join(anndir, "tmp.gff3")
    with zopen(gff) as f, open(tmp_gff, "w") as o:
        for line in f:
            if line.startswith("#"):
                o.write(line)
                continue
            c = line.split("\t")
            c[0] = rename_mha(c[0])
            o.write("\t".join(c))
    name_of = None   # the 2026 combined annotation carries gene_name on every transcript
subprocess.check_call([GFFREAD, tmp_gff, "-T", "-o", gtf])
tss = set()
for line in open(tmp_gff):
    if line.startswith("#"):
        continue
    c = line.rstrip("\n").split("\t")
    if c[2] not in ("transcript", "mRNA"):
        continue
    attrs = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
    if name_of is None:
        name = attrs.get("gene_name") or attrs.get("source_gene_common_name") or attrs.get("Parent", "NA")
    else:
        tid = re.sub(r"\.\d+(_\d+)?$", "", attrs.get("ID", ""))
        name = name_of.get(tid, attrs.get("Parent", "NA"))
    pos = int(c[3]) - 1 if c[6] == "+" else int(c[4]) - 1
    tss.add((c[0], pos, name, c[6]))
with open(tss_bed, "w") as o:
    for ch, pos, name, st in sorted(tss):
        o.write(f"{ch}\t{pos}\t{pos}\t{name}\t.\t{st}\n")
os.remove(tmp_gff)
report.append(f"GTF: gffread -T from {os.path.basename(gff)}; TSS BED: {len(tss)} transcript start sites "
              f"(0-length intervals, strand col 6)")

# ----------------------------------------------------------------------------- RepeatMasker -> SAF / BED
saf = os.path.join(anndir, f"rmsk.{G}.SAF")
safid = os.path.join(anndir, f"rmsk.ids.{G}.SAF")
bed = os.path.join(anndir, f"rmsk.{G}.bed")
n = 0
with open(saf, "w") as s, open(safid, "w") as si, open(bed, "w") as b:
    for x in (s, si):
        x.write("GeneID\tChr\tStart\tEnd\tStrand\n")
    if G == "hs1":
        src = os.path.join(SRC, "hs1/source/hs1.repeatMasker.out.gz")
        for line in zopen(src):
            f = line.split()
            if len(f) < 15 or not f[0].isdigit():
                continue
            ch, st, en, strand, rep, cls = f[4], int(f[5]), int(f[6]), "-" if f[8] == "C" else "+", f[9], f[10]
            eid = f"{rep}_{ch}:{st}-{en}"
            s.write(f"{rep}\t{ch}\t{st}\t{en}\t{strand}\n")
            si.write(f"{eid}\t{ch}\t{st}\t{en}\t{strand}\n")
            b.write(f"{ch}\t{st - 1}\t{en}\t{rep}\t{cls}\t{strand}\n")
            n += 1
    else:
        srcs = [os.path.join(SRC, "mhaESC/mouse.241018.repeats.gff.gz"),
                os.path.join(SRC, "mhaESC/source/release_Y1.1/mT2T-Y_v1.1.repeats.gff.gz")]
        src = " + ".join(os.path.basename(x) for x in srcs)
        for line in (l for x in srcs for l in zopen(x)):
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            m = re.search(r'Motif:([^"]+)"', c[8])
            if not m:
                continue
            ch, st, en, strand, rep, cls = rename_mha(c[0]), int(c[3]), int(c[4]), c[6], m.group(1), c[2]
            eid = f"{rep}.{ch}:{st}-{en}"
            s.write(f"{rep}\t{ch}\t{st}\t{en}\t{strand}\n")
            si.write(f"{eid}\t{ch}\t{st}\t{en}\t{strand}\n")
            b.write(f"{ch}\t{st - 1}\t{en}\t{rep}\t{cls}\t{strand}\n")
            n += 1
report.append(f"RepeatMasker: {n} elements from {os.path.basename(src)} -> rmsk.{G}.SAF (family), "
              f"rmsk.ids.{G}.SAF (element), rmsk.{G}.bed (name, class/family, strand)")
with open(os.path.join(OUT, "PREP_REPORT.txt"), "w") as o:
    o.write("\n".join(report) + "\n")
print("\n".join(report))
