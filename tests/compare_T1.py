#!/usr/bin/env python3
"""T1 validation comparisons (plan: R0-R3). Usage: compare_T1.py <validation T1 dir>  -> writes compare/*.tsv + report.

R0  v2_repeats_A vs v2_repeats_B          noise floor of STAR random multimapper placement (v2 vs itself)
R1  v2_repeats_A vs v3_repeats_legacy      v3 repeats mode with v2 trimming and --keep-mito must reproduce v2 (within R0)
     - repeat-relevant commands identical after path normalisation
     - featureCounts family (.fc.txt), element (.fc.id.txt), IAP coverage
R2  v3_repeats_legacy vs v3_repeats        effect of the new trimming and (3.1+) the chrM split on repeats outputs
R3  v2_genes vs v3_genes                   genes mode: fragments, bigwig correlation (1 kb bins)
R4  (3.3.0, one chromatin mode) a: v3_genes_3.2.0 vs v3_genes   gene outputs (stats, peaks, filtered bigwig); identical
                                   only while FILT_MAPQ = 30; with FILT_MAPQ = 10 (final 3.3.0) the filtered set,
                                   Filtered_fragments, NRF/PBC, TSS_score, peaks and FRiP differ by design
                                b: v2_repeats_A vs v3_genes     chromatin repeat counts, bowtie2 vs v2 STAR:
                                   family frag_all vs v2 .fc.txt (v2 counts mates, keeps duplicates), element
                                   frag_locus vs v2 .fc.id.txt (STAR unique); Repeat_* / Young_locus_pct stats
R5  (3.3.1) v3_genes_3.3.0 vs v3_genes         B1 (Trimmomatic MINLEN 30 -> 20) and B2 (RNA genes-mode STAR cap 10 -> 100):
                                   trimming / alignment / filtered counts, peaks and FRiP per sample; only fragments
                                   of 20-29 bp after trimming (B1) and RNA reads with 11-100 loci (B2) may differ.
                                   v3_repeats_legacy (--legacy-trim) keeps MINLEN 30: R1 stays identical to v2.
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

# ---------------------------------------------------------------- R4 (3.3.0)
def stats(run, s):
    import yaml
    f = os.path.join(sdir(run, s), "stats.yaml")
    if not os.path.exists(f):
        return {}
    y = yaml.safe_load(open(f))
    return y[next(iter(y))]["sample"][s]


def file_md5(f):
    import hashlib
    return hashlib.md5(open(f, "rb").read()).hexdigest() if os.path.exists(f) else None


if os.path.isdir(os.path.join(R, "v3_genes_3.2.0")):
    log("\n## R4a gene outputs: 3.2.0 vs 3.3.0 (v3_genes)\n")
    keys = ["Filtered_fragments", "Duplication_pct", "NRF", "PBC1", "Insert_size_median", "TSS_score", "Peaks_n",
            "Peak_summits_n", "FRiP_pct"]
    rows = []
    for s in names:
        a, b = stats("v3_genes_3.2.0", s), stats("v3_genes", s)
        if b.get("Protocol") == "RNA":
            continue
        row = {"sample": s}
        row["stats_identical"] = all(a.get(k) == b.get(k) for k in keys)
        row["differing"] = ",".join(f"{k}:{a.get(k)}->{b.get(k)}" for k in keys if a.get(k) != b.get(k)) or "-"
        g = b.get("Genome", "mm10")
        pk = [f for f in os.listdir(os.path.join(sdir("v3_genes", s), "peaks_" + g))
              if "_peaks_noBL." in f and f.endswith("Peak")] if os.path.isdir(os.path.join(sdir("v3_genes", s), "peaks_" + g)) else []
        row["peaks_identical"] = all(file_md5(os.path.join(sdir("v3_genes_3.2.0", s), "peaks_" + g, f)) ==
                                     file_md5(os.path.join(sdir("v3_genes", s), "peaks_" + g, f)) for f in pk) if pk else None
        rows.append(row)
    r4a = pd.DataFrame(rows)
    r4a.to_csv(os.path.join(OUT, "R4a_genes_3.2_vs_3.3.tsv"), sep="\t", index=False)
    log(r4a.to_string(index=False))

log("\n## R4b chromatin repeats: v2 STAR repeats mode (v2_repeats_A) vs 3.3.0 bowtie2 (v3_genes)\n")
rows = []
for s in names:
    st = stats("v3_genes", s)
    if not st or st.get("Protocol") == "RNA":
        continue
    g = st.get("Genome", "mm10")
    d = os.path.join(sdir("v3_genes", s), "repeats_" + g)
    fam_f, el_f = os.path.join(d, s + ".repeats.family.tsv"), os.path.join(d, s + ".repeats.elements.tsv.gz")
    row = {"sample": s, "All_fragments": st.get("All_fragments"), "Repeat_pct": st.get("Repeat_pct"),
           "Repeat_locus_pct": st.get("Repeat_locus_pct"), "Young_locus_pct": st.get("Young_locus_pct"),
           "tables": os.path.exists(fam_f) and os.path.exists(el_f)}
    if row["tables"]:
        fam = pd.read_csv(fam_f, sep="\t").set_index("family")
        el = pd.read_csv(el_f, sep="\t").set_index("element_id")
        row["elements_with_frag"] = len(el)
        row["young_elements_pct"] = round(100 * el["flags"].str.contains("young").mean(), 2)
        v2f, v2e = read_fc("v2_repeats_A", s, "family"), read_fc("v2_repeats_A", s, "element")
        if v2f is not None:
            j = pd.concat([v2f, fam.frag_all], axis=1, keys=["v2", "v3"]).fillna(0)
            j = j[(j.v2 + j.v3) > 0]
            row["family_v3_over_v2_total"] = round(j.v3.sum() / j.v2.sum(), 3)
            big = j[(j.v2 + j.v3) / 2 >= 50]
            row["family_pearson_log_ge50"] = round(float(np.corrcoef(np.log1p(big.v2), np.log1p(big.v3))[0, 1]), 4)
        if v2e is not None:
            j = pd.concat([v2e, el.frag_locus], axis=1, keys=["v2", "v3"]).fillna(0)
            j = j[(j.v2 + j.v3) > 0]
            row["element_v3locus_over_v2unique_total"] = round(j.v3.sum() / j.v2.sum(), 3)
            row["element_pearson_log"] = round(float(np.corrcoef(np.log1p(j.v2), np.log1p(j.v3))[0, 1]), 4)
    rows.append(row)
r4b = pd.DataFrame(rows)
r4b.to_csv(os.path.join(OUT, "R4b_repeats_v2_vs_3.3.tsv"), sep="\t", index=False)
log(r4b.to_string(index=False))
if os.path.isdir(os.path.join(R, "v3_genes_3.3.0")):
    log("\n## R5 3.3.0 vs 3.3.1 (v3_genes): B1 MINLEN 20, B2 RNA multimapper cap 100\n")
    keys = ["Trim_input_fragments", "Trimmed_fragments", "Trim_dropped_fragments", "Aligned_fragments", "Unique_fragments",
            "Multimapped_fragments", "Too_many_loci_fragments", "Filtered_fragments", "All_fragments", "Duplication_pct",
            "NRF", "PBC1", "Insert_size_median", "TSS_score", "Peaks_n", "FRiP_pct", "Repeat_pct", "Young_locus_pct"]
    rows = []
    for s in names:
        a, b = stats("v3_genes_3.3.0", s), stats("v3_genes", s)
        if not a or not b:
            continue
        row = {"sample": s, "protocol": b.get("Protocol")}
        for k in keys:
            if a.get(k) != b.get(k):
                row[k] = f"{a.get(k)}->{b.get(k)}"
        row["unchanged"] = all(a.get(k) == b.get(k) for k in keys)
        rows.append(row)
    r5 = pd.DataFrame(rows).fillna("")
    r5.to_csv(os.path.join(OUT, "R5_genes_3.3.0_vs_3.3.1.tsv"), sep="\t", index=False)
    log(r5.T.to_string())

open(os.path.join(OUT, "compare_report.txt"), "w").write("\n".join(lines) + "\n")
