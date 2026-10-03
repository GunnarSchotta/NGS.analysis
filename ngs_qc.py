#!/usr/bin/env python3
"""QC and log-parsing helpers for NGS.analysis v3 (imported by NGS.analysis.py).

Unit convention (docs/stats_units.md): every count reported by v3 is in FRAGMENTS
(paired-end: read pairs; single-end: reads) and every percentage is on a 0-100 scale.
"""
import gzip
import os
import re
import subprocess
from collections import Counter

import numpy as np


# --------------------------------------------------------------------------- logs
def parse_trimmomatic_log(path, paired_end):
    """Parse the Trimmomatic summary line (stderr). Returns dict of fragment counts."""
    txt = open(path).read()
    if paired_end:
        m = re.search(r"Input Read Pairs: (\d+) Both Surviving: (\d+) .*?Forward Only Surviving: (\d+) .*?"
                      r"Reverse Only Surviving: (\d+) .*?Dropped: (\d+)", txt)
        if not m:
            raise ValueError("Trimmomatic PE summary not found in " + path)
        i, b, f, r, d = map(int, m.groups())
        return {"input": i, "both": b, "forward_only": f, "reverse_only": r, "dropped": d}
    m = re.search(r"Input Reads: (\d+) Surviving: (\d+) .*?Dropped: (\d+)", txt)
    if not m:
        raise ValueError("Trimmomatic SE summary not found in " + path)
    i, s, d = map(int, m.groups())
    return {"input": i, "both": s, "forward_only": 0, "reverse_only": 0, "dropped": d}


def parse_star_log(path):
    """Parse STAR Log.final.out by key (not by line number). Counts are reads (SE) or pairs (PE)."""
    vals = {}
    for line in open(path):
        if "|" in line:
            k, v = line.split("|", 1)
            vals[k.strip()] = v.strip()

    def num(key):
        return int(vals[key])

    return {
        "input": num("Number of input reads"),
        "unique": num("Uniquely mapped reads number"),
        "multi": num("Number of reads mapped to multiple loci"),
        "too_many_loci": num("Number of reads mapped to too many loci"),
        "unmapped_mismatches": num("Number of reads unmapped: too many mismatches"),
        "unmapped_short": num("Number of reads unmapped: too short"),
        "unmapped_other": num("Number of reads unmapped: other"),
        "chimeric": num("Number of chimeric reads"),
    }


def parse_picard_dup_metrics(path):
    """Parse Picard MarkDuplicates metrics by header (all library rows summed)."""
    lines = open(path).read().splitlines()
    hi = next(i for i, l in enumerate(lines) if l.startswith("LIBRARY\t"))
    head = lines[hi].split("\t")
    rows = []
    for l in lines[hi + 1:]:
        if not l.strip():
            break
        rows.append(dict(zip(head, l.split("\t"))))

    def s(key):
        return sum(int(float(r[key] or 0)) for r in rows)

    out = {k: s(k) for k in ["UNPAIRED_READS_EXAMINED", "READ_PAIRS_EXAMINED", "UNPAIRED_READ_DUPLICATES",
                             "READ_PAIR_DUPLICATES", "READ_PAIR_OPTICAL_DUPLICATES"]}
    # Picard PERCENT_DUPLICATION = (unpaired dups + 2*pair dups) / (unpaired + 2*pairs); recomputed over all rows
    denom = out["UNPAIRED_READS_EXAMINED"] + 2 * out["READ_PAIRS_EXAMINED"]
    out["PERCENT_DUPLICATION"] = 100.0 * (out["UNPAIRED_READ_DUPLICATES"] + 2 * out["READ_PAIR_DUPLICATES"]) / denom \
        if denom else 0.0
    return out


def parse_featurecounts_summary(path):
    """featureCounts *.summary -> {status: count} (single sample column)."""
    out = {}
    for l in open(path).read().splitlines()[1:]:
        k, v = l.split("\t")[:2]
        out[k] = int(v)
    return out


# --------------------------------------------------------------------------- BAM counting
def count_fragments(samtools, bam, paired_end, flags_req=0, flags_excl=0, region=None, mapq=0):
    """Count fragments in a BAM. PE: R1 records only (adds -f 64); SE: all records (after -F/-f)."""
    req = flags_req | (64 if paired_end else 0)
    cmd = [samtools, "view", "-c", "-f", str(req), "-F", str(flags_excl)]
    if mapq:
        cmd += ["-q", str(mapq)]
    cmd.append(bam)
    if region:
        cmd.append(region)
    return int(subprocess.check_output(cmd).decode().strip())


def canonical_bed(fai, canonical_regex, mito, out_bed):
    """BED of canonical chromosomes (regex on names), mitochondrial chromosome excluded."""
    rx = re.compile(r"^(" + canonical_regex + r")$")
    with open(fai) as f, open(out_bed, "w") as o:
        for line in f:
            c, ln = line.split("\t")[:2]
            if rx.match(c) and c != mito:
                o.write(f"{c}\t0\t{ln}\n")
    return out_bed


def chrom_present(fai, chrom):
    return any(l.split("\t")[0] == chrom for l in open(fai))


# --------------------------------------------------------------------------- insert size
def insert_size_qc(samtools, bam, out_prefix, max_len=1000):
    """Fragment-length histogram from proper-pair R1 records (|TLEN|), stats + plot.
    Writes <out_prefix>.fraglen.txt (len, count), .fraglen.pdf/.png. Returns stats dict."""
    p = subprocess.Popen([samtools, "view", "-f", "67", "-F", "3340", bam], stdout=subprocess.PIPE, text=True)
    cnt = Counter()
    for line in p.stdout:
        t = abs(int(line.split("\t", 9)[8]))
        if t > 0:
            cnt[t] += 1
    p.wait()
    n = sum(cnt.values())
    with open(out_prefix + ".fraglen.txt", "w") as o:
        for k in sorted(cnt):
            o.write(f"{k}\t{cnt[k]}\n")
    if n == 0:
        return {}
    lens = np.array(sorted(cnt))
    cum = np.cumsum([cnt[k] for k in lens])
    median = int(lens[np.searchsorted(cum, n / 2)])
    mode = int(max(cnt, key=cnt.get))
    stats = {"Insert_size_median": median, "Insert_size_mode": mode,
             "Insert_lt150_pct": round(100 * sum(v for k, v in cnt.items() if k < 150) / n, 2),
             "Insert_150_300_pct": round(100 * sum(v for k, v in cnt.items() if 150 <= k < 300) / n, 2),
             "Insert_ge300_pct": round(100 * sum(v for k, v in cnt.items() if k >= 300) / n, 2)}
    _plot_fraglen(cnt, n, os.path.basename(out_prefix), out_prefix, max_len)
    return stats


def _plot_fraglen(cnt, n, name, out_prefix, max_len):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x = np.arange(1, max_len + 1)
    y = np.array([cnt.get(i, 0) for i in x]) / n
    fig, ax = plt.subplots(figsize=(4, 3.2))
    ax.plot(x, y, lw=0.8, color="#2166ac")
    ax.set_xlabel("fragment length (bp)")
    ax.set_ylabel("fraction of fragments")
    ax.set_title(name, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_prefix + ".fraglen.pdf")
    fig.savefig(out_prefix + ".fraglen.png", dpi=110)
    plt.close(fig)


# --------------------------------------------------------------------------- library complexity
def library_complexity(bam, paired_end, mapq, tmpdir, threads=4):
    """ENCODE library complexity on the duplicate-MARKED (not removed) BAM, uniquely mapped fragments.
    PE: fragment key = (chrom, leftmost pos, TLEN, R1 strand) from proper-pair R1 records.
    SE: key = (chrom, 5' end, strand).
    NRF = distinct/total; PBC1 = M1/distinct; PBC2 = M1/M2 (M1/M2: positions seen exactly once/twice)."""
    import pysam
    sort = subprocess.Popen(["sort", "-S", "2G", "--parallel=" + str(threads), "-T", tmpdir],
                            stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
    uniq = subprocess.Popen(["uniq", "-c"], stdin=sort.stdout, stdout=subprocess.PIPE, text=True)
    with pysam.AlignmentFile(bam, "rb") as b:
        for r in b.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_qcfail or r.mapping_quality < mapq:
                continue
            if paired_end:
                if not (r.is_read1 and r.is_proper_pair) or r.template_length == 0:
                    continue
                left = min(r.reference_start, r.next_reference_start)
                key = f"{r.reference_name}\t{left}\t{abs(r.template_length)}\t{int(r.is_reverse)}\n"
            else:
                five = r.reference_end if r.is_reverse else r.reference_start
                key = f"{r.reference_name}\t{five}\t{int(r.is_reverse)}\n"
            sort.stdin.write(key)
    sort.stdin.close()
    m1 = m2 = d = t = 0
    for line in uniq.stdout:
        c = int(line.split()[0])
        t += c
        d += 1
        m1 += c == 1
        m2 += c == 2
    uniq.wait()
    sort.wait()
    if t == 0:
        return {}
    return {"NRF": round(d / t, 4), "PBC1": round(m1 / d, 4), "PBC2": round(m1 / m2, 4) if m2 else None}


# --------------------------------------------------------------------------- TSS enrichment
def tss_score(profile_file, flank_bins=100, search_halfwidth=500, smooth_halfwidth=50):
    """TSS enrichment score from a per-bp insertion profile centred on TSSs (pyTssEnrichment -v output).
    Normalise by the mean of the outermost `flank_bins` bins on BOTH sides; the score is the mean of the
    normalised profile in +-smooth_halfwidth around the maximum found within +-search_halfwidth of the TSS.
    Returns (score, normalised_profile) or (0, None) if the flanks are empty."""
    v = np.loadtxt(profile_file, dtype=float).ravel()
    if v.size < 2 * flank_bins + 1:
        return 0.0, None
    flank = np.r_[v[:flank_bins], v[-flank_bins:]].mean()
    if flank <= 0:
        return 0.0, None
    norm = v / flank
    c = v.size // 2
    lo, hi = max(0, c - search_halfwidth), min(v.size, c + search_halfwidth + 1)
    peak = lo + int(np.argmax(norm[lo:hi]))
    s = norm[max(0, peak - smooth_halfwidth):min(v.size, peak + smooth_halfwidth + 1)].mean()
    return round(float(s), 2), norm


def plot_tss(norm, score, name, out_prefix, cutoff=6.0):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x = np.arange(norm.size) - norm.size // 2
    k = 25
    sm = np.convolve(norm, np.ones(k) / k, mode="same")
    fig, ax = plt.subplots(figsize=(4, 3.6))
    ax.plot(x, sm, color="#1a9850" if score >= cutoff else "#d73027", lw=1)
    ax.set_xlabel("distance from TSS (bp)")
    ax.set_ylabel("TSS enrichment (flank-normalised)")
    ax.set_title(f"{name}: TSS score {score}", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_prefix + "_TSS_enrichment.pdf")
    fig.savefig(out_prefix + "_TSS_enrichment.png", dpi=110)
    plt.close(fig)


# --------------------------------------------------------------------------- peaks
def peaks_to_saf(peak_file, saf):
    with open(peak_file) as f, open(saf, "w") as o:
        o.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for i, line in enumerate(f):
            c = line.split("\t")
            o.write(f"peak{i}\t{c[0]}\t{int(c[1]) + 1}\t{c[2]}\t.\n")


def count_lines(path):
    if not os.path.exists(path):
        return 0
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        return sum(1 for _ in f)
