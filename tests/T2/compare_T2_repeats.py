#!/usr/bin/env python3
"""T2 repeats mode: v3 (new trimming, repeats commands identical to v2) vs the existing v2 outputs of the same samples
(01.Angela 00.paper). Differences = trimming change + random multimapper placement (noise floor: T1 R0).
Usage: compare_T2_repeats.py      writes T2/angela_repeats_v2_vs_v3.tsv (+ _families.tsv) and prints a summary
Per sample: STAR input / unique / multi / too-many-loci (Log.final.out, in pairs), featureCounts assigned (family -M,
element unique), family counts (families with >= 100 v2 counts: Pearson r of log2, % families |log2FC| > 0.1,
median log2FC), element counts (elements with >= 20 v2 counts: same), IAP normalised coverage (Pearson r).

3.3.0: the RNA sample is still v2-style STAR (repeats mode frozen), so the comparison above holds for it unchanged.
Chromatin samples (the H3K9me3 ChIP) are now bowtie2 in ONE mode: no STAR log and no feature_counts/ folder; the
repeat tables are repeats_mm10/<s>.repeats.{family.tsv,elements.tsv.gz}. For these the v2 counts are compared with
frag_all (family) and with both element tiers (frag_locus MAPQ >= 10, frag_all). Not identical by design: v2 counts
mates and keeps duplicates, v3 counts fragments after duplicate removal, and STAR "unique" differs from the tiers.
"""
import os

import numpy as np
import pandas as pd

V2 = "/store24/project24/becgsc_001/schottalab/01.Angela/01.primary.processing/00.paper"
V3 = "/store24/project24/becgsc_001/analysis/NGS.analysis.v3.validation/T2/runs/angela_repeats/results/samples"
OUT = "/store24/project24/becgsc_001/analysis/NGS.analysis.v3.validation/T2/angela_repeats_v2_vs_v3"
SAMPLES = {"RNA.T263.ES.d0_r1": f"{V2}/RNAseq/01.ctrl.Setdb1KO/results/samples/RNA.T263.ES.d0_r1",
           "wt26.ES.H3K9me3_r1": f"{V2}/Histones/01.ES.XEN.wt26/results/samples/wt26.ES.H3K9me3_r1"}
STAR_KEYS = {"Number of input reads": "input", "Uniquely mapped reads number": "unique",
             "Number of reads mapped to multiple loci": "multi", "Number of reads mapped to too many loci": "too_many"}


def star(d, s):
    r = {}
    for line in open(f"{d}/aligned_mm10/{s}.Log.final.out"):
        k, _, v = line.partition("|")
        if k.strip() in STAR_KEYS:
            r[STAR_KEYS[k.strip()]] = int(v.strip())
    return r


def fc(path):   # fc.txt: "Geneid Length count" (space, from awk); fc.id.txt: featureCounts table (tab)
    fam = path.endswith(".fc.txt")
    x = pd.read_csv(path, sep=" " if fam else "\t", comment="#", usecols=[0, 2] if fam else [0, 6], dtype={0: str})
    x.columns = ["id", "count"]
    return x.groupby("id")["count"].sum()


def assigned(path):
    return int(pd.read_csv(path, sep="\t", index_col=0).loc["Assigned"].iloc[0])


def compare(a, b, min_count):
    j = pd.concat([a.rename("v2"), b.rename("v3")], axis=1).fillna(0)
    j = j[j.v2 >= min_count]
    l2 = np.log2(j.v3 + 1) - np.log2(j.v2 + 1)
    return j, l2, {"n": len(j), "pearson_log2": round(np.corrcoef(np.log2(j.v2 + 1), np.log2(j.v3 + 1))[0, 1], 5),
                   "pct_abs_log2FC_gt_0.1": round(100 * (l2.abs() > 0.1).mean(), 2),
                   "median_log2FC": round(l2.median(), 4)}


def chromatin_row(s, d2, d3, fams):
    """3.3.0 chromatin sample (bowtie2, tiers all / locus) vs the v2 STAR repeats-mode outputs."""
    import yaml
    y = yaml.safe_load(open(f"{d3}/stats.yaml"))
    st = y[next(iter(y))]["sample"][s]
    fam = pd.read_csv(f"{d3}/repeats_mm10/{s}.repeats.family.tsv", sep="\t").set_index("family")
    el = pd.read_csv(f"{d3}/repeats_mm10/{s}.repeats.elements.tsv.gz", sep="\t").set_index("element_id")
    r = {"sample": s, "v3_mapper": "bowtie2 (3.3.0 chromatin)"}
    for k, n in star(d2, s).items():
        r[f"v2_STAR_{k}"] = n
    r["v2_fc_family_assigned"] = assigned(f"{d2}/feature_counts/{s}.fc.tmp.txt.summary")
    r["v2_fc_element_assigned"] = assigned(f"{d2}/feature_counts/{s}.fc.id.txt.summary")
    for k in ("All_fragments", "Filtered_fragments", "Repeat_pct", "Repeat_locus_pct", "Young_locus_pct"):
        r[f"v3_{k}"] = st.get(k)
    v2f, v2e = fc(f"{d2}/feature_counts/{s}.fc.txt"), fc(f"{d2}/feature_counts/{s}.fc.id.txt")
    r["family_v3all_over_v2_total"] = round(float(fam.frag_all.sum() / v2f.sum()), 3)
    j, l2, m = compare(v2f, fam.frag_all, 100)
    r.update({f"family_{k}": v for k, v in m.items()})
    fams.append(j.assign(log2FC=l2.round(4), sample=s))
    for tier, col in (("locus", "frag_locus"), ("all", "frag_all")):
        _, _, m = compare(v2e, el[col], 20)
        r.update({f"element_{tier}_{k}": v for k, v in m.items()})
    r["element_v3locus_over_v2unique_total"] = round(float(el.frag_locus.sum() / v2e.sum()), 3)
    i2, i3 = f"{d2}/IAP_coverage/{s}.IAP.norm.coverage.txt", f"{d3}/IAP_coverage/{s}.IAP.norm.coverage.txt"
    if os.path.exists(i2) and os.path.exists(i3):
        a, b = (pd.read_csv(f, sep="\t", header=None).iloc[:, -1] for f in (i2, i3))
        if len(a) == len(b):
            r["IAP_cov_pearson"] = round(np.corrcoef(a, b)[0, 1], 5)
            r["IAP_cov_sum_change_pct"] = round(100 * (b.sum() - a.sum()) / a.sum(), 2)
    return r


rows, fams = [], []
for s, d2 in SAMPLES.items():
    d3 = f"{V3}/{s}"
    if os.path.exists(f"{d3}/repeats_mm10/{s}.repeats.family.tsv"):
        rows.append(chromatin_row(s, d2, d3, fams))
        continue
    if not os.path.exists(f"{d3}/feature_counts/{s}.fc.txt"):
        print("skip", s, "(v3 not finished)")
        continue
    r = {"sample": s}
    for v, d in (("v2", d2), ("v3", d3)):
        for k, n in star(d, s).items():
            r[f"{v}_STAR_{k}"] = n
        r[f"{v}_fc_family_assigned"] = assigned(f"{d}/feature_counts/{s}.fc.tmp.txt.summary")
        r[f"{v}_fc_element_assigned"] = assigned(f"{d}/feature_counts/{s}.fc.id.txt.summary")
    for v in ("STAR_input", "STAR_unique", "STAR_multi", "fc_family_assigned", "fc_element_assigned"):
        r[f"{v}_change_pct"] = round(100 * (r[f"v3_{v}"] - r[f"v2_{v}"]) / r[f"v2_{v}"], 2)
    j, l2, m = compare(fc(f"{d2}/feature_counts/{s}.fc.txt"), fc(f"{d3}/feature_counts/{s}.fc.txt"), 100)
    r.update({f"family_{k}": v for k, v in m.items()})
    fams.append(j.assign(log2FC=l2.round(4), sample=s))
    _, _, m = compare(fc(f"{d2}/feature_counts/{s}.fc.id.txt"), fc(f"{d3}/feature_counts/{s}.fc.id.txt"), 20)
    r.update({f"element_{k}": v for k, v in m.items()})
    i2, i3 = (f"{d}/IAP_coverage/{s}.IAP.norm.coverage.txt" for d in (d2, d3))
    if os.path.exists(i2) and os.path.exists(i3):
        a, b = (pd.read_csv(f, sep="\t", header=None).iloc[:, -1] for f in (i2, i3))
        if len(a) == len(b):
            r["IAP_cov_pearson"] = round(np.corrcoef(a, b)[0, 1], 5)
            r["IAP_cov_sum_change_pct"] = round(100 * (b.sum() - a.sum()) / a.sum(), 2)
    rows.append(r)

if rows:
    t = pd.DataFrame(rows)
    t.to_csv(OUT + ".tsv", sep="\t", index=False)
    pd.concat(fams).to_csv(OUT + "_families.tsv", sep="\t")
    print(t.set_index("sample").T.to_string())
    for f in fams:
        print("\n", f["sample"].iloc[0], "largest family changes (v2 >= 100):")
        print(f.reindex(f.log2FC.abs().sort_values(ascending=False).index).head(10)[["v2", "v3", "log2FC"]].to_string())
