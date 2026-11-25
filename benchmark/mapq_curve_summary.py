#!/usr/bin/env python3
"""Summarise mapq_curve.py output: for every aligner and MAPQ cutoff, per family group,
kept % (fragments with MAPQ >= cutoff / all simulated fragments of the group), precision % (at the true locus /
kept) and correct % (at the true locus with MAPQ >= cutoff / all simulated), i.e. the usable locus-level signal.

Usage: mapq_curve_summary.py <curve.tsv> [<curve.tsv> ...] > summary.md   (also writes <curve>.summary.tsv)
"""
import re
import sys

import pandas as pd

GROUPS = [("all repeats", r".*"), ("L1MdT", r"^L1Md_T$"), ("L1MdA", r"^L1Md_A$"), ("L1MdGf", r"^L1Md_Gf$"),
          ("IAPEz-int", r"^IAPEz-int$"), ("MMERVK10C-int", r"^MMERVK10C-int$"), ("ETnERV-int", r"^ETnERV-int$"),
          ("MERVL-int", r"^MERVL-int$"), ("RLTR10C", r"^RLTR10C$"), ("B1_Mm", r"^B1_Mm$"),
          # eval_snp.py: fragments from non-reference TE insertions (no true locus: every placement is wrong)
          ("nonref young TE (all)", r"^nonref:"), ("nonref L1Md_T/A", r"^nonref:L1Md_[TA]$"),
          ("nonref ERV (IAPEz, MERVL, MMERVK10C, ETnERV)", r"^nonref:(IAPEz|MERVL|MMERVK10C|ETnERV)-int$")]
CUTS = {"bowtie2_v3genes": [0, 1, 2, 3, 5, 10, 20, 30, 40, 42], "STAR_repeats": [0, 1, 3, 255]}

out = []
for path in sys.argv[1:]:
    d = pd.read_csv(path, sep="\t")
    rows = []
    for label, dl in d.groupby("label"):
        for gname, rx in GROUPS:
            g = dl[dl.family.str.match(rx)]
            if not gname.startswith("nonref"):
                g = g[~g.family.str.startswith("nonref:")]
            if g.empty:
                continue
            n = g.drop_duplicates("family").n_family.sum()
            for c in CUTS.get(label, sorted(g.mapq.unique())):
                k = g[g.mapq >= c]
                kept, ok = k.n_at_mapq.sum(), k.correct_at_mapq.sum()
                rows.append(dict(set=path.rsplit("/", 1)[-1].replace(".tsv", ""), aligner=label, group=gname,
                                 n=int(n), mapq_ge=c, kept_pct=round(100 * kept / n, 2),
                                 precision_pct=round(100 * ok / kept, 2) if kept else None,
                                 correct_pct=round(100 * ok / n, 2)))
    df = pd.DataFrame(rows)
    df.to_csv(path.replace(".tsv", ".summary.tsv"), sep="\t", index=False)
    out.append(df)
df = pd.concat(out)
for (s, grp), g in df.groupby(["set", "group"], sort=False):
    print(f"\n### {s}: {grp} (n = {g.n.iloc[0]:,})\n")
    print("| aligner | MAPQ >= | kept % | precision % | correct % |\n|---|---|---|---|---|")
    for r in g.itertuples():
        print(f"| {r.aligner} | {r.mapq_ge} | {r.kept_pct} | {r.precision_pct} | {r.correct_pct} |")
