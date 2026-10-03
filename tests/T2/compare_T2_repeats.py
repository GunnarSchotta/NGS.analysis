#!/usr/bin/env python3
"""T2 repeats mode: v3 (new trimming, repeats commands identical to v2) vs the existing v2 outputs of the same samples
(01.Angela 00.paper). Differences = trimming change + random multimapper placement (noise floor: T1 R0).
Usage: compare_T2_repeats.py      writes T2/angela_repeats_v2_vs_v3.tsv (+ _families.tsv) and prints a summary
Per sample: STAR input / unique / multi / too-many-loci (Log.final.out, in pairs), featureCounts assigned (family -M,
element unique), family counts (families with >= 100 v2 counts: Pearson r of log2, % families |log2FC| > 0.1,
median log2FC), element counts (elements with >= 20 v2 counts: same), IAP normalised coverage (Pearson r).
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


rows, fams = [], []
for s, d2 in SAMPLES.items():
    d3 = f"{V3}/{s}"
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
