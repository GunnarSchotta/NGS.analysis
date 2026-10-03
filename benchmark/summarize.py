#!/usr/bin/env python3
"""Summarise the benchmark evaluations (run_benchmark.sh outputs) into tables.
Usage: summarize.py <benchmark_dir>     reads <dir>/<genome>/eval_*.tsv, writes <dir>/summary_*.tsv + summary_tables.md
  summary_overall.tsv : genome x read set x arm: aligned %, unique %, TP % of unique (all repeat-overlapping fragments)
  summary_groups.tsv  : the same for selected TE groups. Mouse names are grouped across the RepBase names of the mm10
                        UCSC track (L1Md_T, L1Md_A) and the Dfam subfamily names of mhaESC (L1MdTf_I-III, L1MdA_I-VII)
  summary_adapter.tsv : adapter read-through set, recovered % per fragment-length bin
Arms: v2genes (bowtie2, unique = MAPQ>=42), v3genes (bowtie2 --dovetail, MAPQ>=30), repeats (STAR random-1, MAPQ 255),
repeats_any (same STAR BAM, any MAPQ: TP % = correct locus of all placed fragments, incl. random multimapper placement).
"""
import glob
import os
import re
import sys

import pandas as pd

GROUPS = {
    "mouse": [("L1MdT", r"^(L1Md_T|L1MdTf_.*)$"), ("L1MdA", r"^(L1Md_A|L1MdA_.*)$"), ("L1MdGf", r"^(L1Md_Gf|L1MdGf_.*)$"),
              ("IAPEz-int", r"^IAPEz-int$"), ("MMERVK10C-int", r"^MMERVK10C-int$"), ("ETnERV-int", r"^ETnERV-int$"),
              ("RLTR10C", r"^RLTR10C$"), ("B2_Mm1a", r"^B2_Mm1a$"), ("B1_Mm", r"^B1_Mm$")],
    "human": [(f, "^" + re.escape(f) + "$") for f in
              ["L1HS", "L1PA2", "L1PA3", "AluY", "AluYa5", "SVA_D", "SVA_F", "HERVK-int", "LTR5_Hs"]],
}
SPECIES = {"mm10": "mouse", "mhaESC": "mouse", "hg38": "human", "hs1": "human"}
ARMS = ["v2genes", "v3genes", "repeats", "repeats_any"]


def agg(d):
    al = (d.aligned_pct * d.n / 100).sum()
    un = (d.unique_pct * d.n / 100).sum()
    tp = (d.unique_TP_pct.fillna(0) * d.unique_pct * d.n / 1e4).sum()
    n = d.n.sum()
    return {"n": int(n), "aligned_pct": round(100 * al / n, 2), "unique_pct": round(100 * un / n, 2),
            "unique_TP_pct": round(100 * tp / un, 2) if un else None}


def md(df):
    cols = list(df.columns)
    out = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
    out += ["| " + " | ".join("" if pd.isna(v) else str(v) for v in r) + " |" for r in df.itertuples(index=False)]
    return "\n".join(out)


B = sys.argv[1]
ov, gr, ad = [], [], []
for g in SPECIES:
    fs = glob.glob(os.path.join(B, g, "eval_art_*.tsv"))
    if fs:
        d = pd.concat([pd.read_csv(f, sep="\t") for f in fs])
        lab = d.label.str.split(":")
        d["set"], d["arm"] = lab.str[1].str.replace("art_", "2x"), lab.str[2]
        for (st, arm), x in d.groupby(["set", "arm"]):
            ov.append({"genome": g, "set": st, "arm": arm, **agg(x)})
            for name, rx in GROUPS[SPECIES[g]]:
                y = x[x.family.str.match(rx)]
                if len(y):
                    gr.append({"genome": g, "group": name, "set": st, "arm": arm, **agg(y)})
    for f in sorted(glob.glob(os.path.join(B, g, "eval_adapt_v*.tsv"))):
        d = pd.read_csv(f, sep="\t")
        d["bin"] = pd.cut(d.frag_len, [29, 52, 59, 79, 99, 150], labels=["30-52", "53-59", "60-79", "80-99", "100-150"])
        for b, x in d.groupby("bin", observed=True):
            ad.append({"genome": g, "arm": d.label.iloc[0].split(":")[2], "frag_len": b, "n": int(x.n.sum()),
                       "recovered_pct": round(100 * x.recovered.sum() / x.n.sum(), 1)})

arm_rank = {a: i for i, a in enumerate(ARMS)}
sections = []
if ov:
    rank = {**arm_rank, **{gg: i for i, gg in enumerate(SPECIES)}}
    o = pd.DataFrame(ov).sort_values(["genome", "set", "arm"], key=lambda s: s.map(rank) if s.name != "set" else s)
    o.to_csv(os.path.join(B, "summary_overall.tsv"), sep="\t", index=False)
    sections.append("## All repeat-overlapping fragments\n\n" + md(o))
if gr:
    t = pd.DataFrame(gr)
    t.to_csv(os.path.join(B, "summary_groups.tsv"), sep="\t", index=False)
    t = t[t.arm != "repeats_any"].copy()
    t["unique / TP"] = t.unique_pct.round(1).astype(str) + " / " + t.unique_TP_pct.round(1).astype(str)
    for sp in ("mouse", "human"):
        x = t[t.genome.map(SPECIES) == sp]
        if len(x):
            w = x.pivot_table(index=["group", "set"], columns=["genome", "arm"], values="unique / TP", aggfunc="first")
            w = w[[c for c in [(gg, a) for gg in SPECIES if SPECIES[gg] == sp for a in ARMS] if c in w.columns]]
            order = [n for n, _ in GROUPS[sp]]
            w = w.reindex(sorted(w.index, key=lambda i: (order.index(i[0]), i[1] != "2x100")))
            w.columns = [f"{gg} {a}" for gg, a in w.columns]
            sections.append(f"## {sp}: unique % / TP % of unique, per TE group\n\n" + md(w.reset_index()))
if ad:
    a = pd.DataFrame(ad)
    a.to_csv(os.path.join(B, "summary_adapter.tsv"), sep="\t", index=False)
    w = a.pivot_table(index=["genome", "arm"], columns="frag_len", values="recovered_pct")
    w = w[["30-52", "53-59", "60-79", "80-99", "100-150"]]
    w = w.reindex(sorted(w.index, key=lambda i: (list(SPECIES).index(i[0]), i[1]))).reset_index()
    sections.append("## Adapter read-through: recovered % per fragment length (bp)\n\n" + md(w))
open(os.path.join(B, "summary_tables.md"), "w").write("\n\n".join(sections) + "\n")
print("\n\n".join(sections))
