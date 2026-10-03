#!/usr/bin/env python3
"""T1 validation comparisons (plan: R0-R3). Usage: compare_T1.py <validation T1 dir>  -> writes compare/*.tsv + report.

R0  v2_repeats_A vs v2_repeats_B          noise floor of STAR random multimapper placement (v2 vs itself)
R1  v2_repeats_A vs v3_repeats_legacy      v3 repeats mode with v2 trimming must reproduce v2 (within R0)
     - repeat-relevant commands identical after path normalisation
     - featureCounts family (.fc.txt), element (.fc.id.txt), IAP coverage
R2  v3_repeats_legacy vs v3_repeats        effect of the new trimming alone on repeats outputs
R3  v2_genes vs v3_genes                   genes mode: fragments, bigwig correlation (1 kb bins)
"""
import os
import re
import subprocess
import sys

import numpy as np
import pandas as pd

V = sys.argv[1]
R = os.path.join(V, "runs")
OUT = os.path.join(V, "compare")
os.makedirs(OUT, exist_ok=True)
SAMPLES = pd.read_csv(os.path.join(R, "v2_genes", "sample.table.csv"))
names = SAMPLES.sample_name.tolist()
lines = []


def log(s=""):
    print(s)
    lines.append(s)


def sdir(run, s):
    return os.path.join(R, run, "results", "samples", s)


def read_fc(run, s, level):
    f = os.path.join(sdir(run, s), "feature_counts", s + (".fc.txt" if level == "family" else ".fc.id.txt"))
    if not os.path.exists(f):
        return None
    if level == "family":
        d = pd.read_csv(f, sep=" ", skiprows=2, header=None, names=["id", "len", "n"])
    else:
        d = pd.read_csv(f, sep="\t", skiprows=2, header=None, names=["id", "chr", "s", "e", "st", "len", "n"])
    return d.set_index("id")["n"]


def read_iap(run, s):
    f = os.path.join(sdir(run, s), "IAP_coverage", s + ".IAP.norm.coverage.txt")
    return np.loadtxt(f) if os.path.exists(f) else None


def diff_stats(a, b, min_count=100):
    j = pd.concat([a, b], axis=1, keys=["a", "b"]).fillna(0)
    tot_a, tot_b = j.a.sum(), j.b.sum()
    big = j[(j.a + j.b) / 2 >= min_count]
    lfc = np.log2((big.b + 1) / (big.a + 1))
    return {"total_a": int(tot_a), "total_b": int(tot_b), "total_rel_diff_pct": round(100 * (tot_b - tot_a) / tot_a, 3) if tot_a else np.nan,
            "n_features_ge%d" % min_count: len(big), "max_abs_log2fc": round(float(lfc.abs().max()), 4) if len(big) else 0,
            "n_abs_log2fc_gt_0.1": int((lfc.abs() > 0.1).sum()),
            "pearson_log": round(float(np.corrcoef(np.log1p(j.a), np.log1p(j.b))[0, 1]), 6)}


def compare_repeats(run_a, run_b, label):
    rows = []
    for s in names:
        for level in ("family", "element"):
            a, b = read_fc(run_a, s, level), read_fc(run_b, s, level)
            if a is None or b is None:
                continue
            rows.append({"comparison": label, "sample": s, "level": level, **diff_stats(a, b)})
        ia, ib = read_iap(run_a, s), read_iap(run_b, s)
        if ia is not None and ib is not None:
            rows.append({"comparison": label, "sample": s, "level": "IAP_coverage", "total_a": round(ia.sum(), 1),
                         "total_b": round(ib.sum(), 1), "total_rel_diff_pct": round(100 * (ib.sum() - ia.sum()) / ia.sum(), 3),
                         "pearson_log": round(float(np.corrcoef(ia, ib)[0, 1]), 6)})
    return pd.DataFrame(rows)


def norm_cmds(run, s):
    """Repeat-relevant commands from NGS.analysis_commands.sh, paths normalised."""
    f = os.path.join(sdir(run, s), "NGS.analysis_commands.sh")
    keep = ("STAR ", " sort ", " index ", "MarkDuplicates", "view -b -q 255", "--normalizeUsing RPKM",
            "coverage -g", "featureCounts -M", "featureCounts -F SAF -T 1", "awk -vN=15413", "awk '{print $1,$6,$7}'")
    out = []
    for line in open(f):
        line = line.strip()
        if not line or line.startswith("#") or not any(k in line for k in keep):
            continue
        line = re.sub(r"/store24/\S*/(runs/[^/]+)/", "RUN/", line)
        line = re.sub(r"/store24/project24/becgsc_001/micromamba/envs/ngs\.v\d/bin/", "", line)
        line = re.sub(r"/store24/project24/becgsc_001/coding/NGS\.analysis(\.v3)?/", "CODE/", line)
        out.append(line)
    return out


log("# T1 validation comparisons\n")
rep_tsv = os.path.join(OUT, "repeats_R0_R1_R2.tsv")
if os.path.exists(rep_tsv):
    allrep = pd.read_csv(rep_tsv, sep="\t")
else:
    r0 = compare_repeats("v2_repeats_A", "v2_repeats_B", "R0 v2 vs v2")
    r1 = compare_repeats("v2_repeats_A", "v3_repeats_legacy", "R1 v2 vs v3-legacy-trim")
    r2 = compare_repeats("v3_repeats_legacy", "v3_repeats", "R2 v3 legacy vs new trim")
    allrep = pd.concat([r0, r1, r2])
    allrep.to_csv(rep_tsv, sep="\t", index=False)
cols = ["comparison", "sample", "level", "total_rel_diff_pct", "n_abs_log2fc_gt_0.1", "max_abs_log2fc", "pearson_log"]
log("## Repeats mode (featureCounts family/element, IAP)\n")
log(allrep[cols].to_string(index=False))

log("\n## R1 command identity (repeat-relevant commands, paths normalised)\n")
for s in names:
    try:
        a, b = norm_cmds("v2_repeats_A", s), norm_cmds("v3_repeats_legacy", s)
    except FileNotFoundError:
        continue
    diff = [x for x in a if x not in b] + ["+" + x for x in b if x not in a]
    log(f"{s}: {len(a)} vs {len(b)} commands; {'IDENTICAL' if not diff else 'DIFFERENT'}")
    for d in diff[:6]:
        log("    " + d[:250])

log("\n## R3 genes mode: v2 vs v3\n")
rows = []
for s in names:
    import yaml
    st = {}
    for run in ("v2_genes", "v3_genes"):
        y = yaml.safe_load(open(os.path.join(sdir(run, s), "stats.yaml")))
        st[run] = y[next(iter(y))]["sample"][s]
    a, b = st["v2_genes"], st["v3_genes"]
    row = {"sample": s, "v2_Mapped_reads_filtered": a.get("Mapped_reads_filtered"), "v3_Filtered_fragments": b.get("Filtered_fragments"),
           "v2_Trimmed_reads": a.get("Trimmed_reads"), "v3_Trimmed_fragments": b.get("Trimmed_fragments"),
           "v3_Duplication_pct": b.get("Duplication_pct"), "v3_FRiP_pct": b.get("FRiP_pct"), "v3_Peaks_n": b.get("Peaks_n"),
           "v2_TSS_score": a.get("TSS_score"), "v3_TSS_score": b.get("TSS_score")}
    bw2 = os.path.join(sdir("v2_genes", s), "aligned_" + b.get("Genome", "mm10"), s + (".bw" if b.get("Protocol") == "RNA" else ".dedup.unique.bw"))
    bw3 = os.path.join(sdir("v3_genes", s), "aligned_" + b.get("Genome", "mm10"), s + (".bw" if b.get("Protocol") == "RNA" else ".filt.bw"))
    if os.path.exists(bw2) and os.path.exists(bw3):
        npz = os.path.join(OUT, f"bw_{s}.npz")
        raw = os.path.join(OUT, f"bw_{s}.tsv")
        if not os.path.exists(raw):
            subprocess.run(["multiBigwigSummary", "bins", "-b", bw2, bw3, "-bs", "1000", "-p", "8", "-o", npz,
                            "--outRawCounts", raw, "--chromosomesToSkip", "chrM", "chrY"], check=True)
        d = pd.read_csv(raw, sep="\t").iloc[:, 3:5].fillna(0)
        d = d[(d.iloc[:, 0] > 0) | (d.iloc[:, 1] > 0)]
        row["bigwig_1kb_spearman"] = round(float(d.iloc[:, 0].corr(d.iloc[:, 1], method="spearman")), 4)
        row["bigwig_1kb_pearson_log"] = round(float(np.corrcoef(np.log1p(d.iloc[:, 0]), np.log1p(d.iloc[:, 1]))[0, 1]), 4)
    rows.append(row)
r3 = pd.DataFrame(rows)
r3.to_csv(os.path.join(OUT, "genes_R3.tsv"), sep="\t", index=False)
log(r3.to_string(index=False))
open(os.path.join(OUT, "compare_report.txt"), "w").write("\n".join(lines) + "\n")
