#!/usr/bin/env python3
"""Make a 'sample genome' from one reference chromosome: random SNPs, small indels and non-reference TE insertions.

Usage: mutate_genome.py <chrom.fa> <rmsk.bed> <snp_rate> <n_insertions> <seed> <out_prefix>

- SNPs: uniform, rate per bp (bases N untouched).
- Indels: rate = snp_rate / 10, 50 % insertions / 50 % deletions, length 1-10 bp (geometric, p = 0.4).
- Non-reference TE insertions: copies of full-length young elements of this chromosome (families in YOUNG, at least
  MINLEN bp), 0.5 % extra substitutions, inserted at random positions outside existing young elements. Reads from them
  have no true locus in the reference.
Writes <out>.fa (sequence named as the input chromosome) and <out>.events.tsv:
  mut_pos (0-based start in the mutated sequence), ref_pos (0-based reference position), kind (SNP/INS/DEL/TE),
  length (bases added (+) or removed (-) relative to the reference), family (TE only).
Coordinate conversion mutated -> reference is done by mutated_to_ref() in eval_snp.py from these events.
"""
import random
import sys

YOUNG = {"L1Md_T": 5000, "L1Md_A": 5000, "IAPEz-int": 4000, "MMERVK10C-int": 4000, "ETnERV-int": 3000,
         "MERVL-int": 3000}
DIVERGE = 0.005
BASES = "ACGT"

fa, rmsk, snp_rate, n_ins, seed, out = sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
rng = random.Random(seed)
name, seq = None, []
for line in open(fa):
    if line.startswith(">"):
        name = line[1:].split()[0]
    else:
        seq.append(line.strip().upper())
ref = "".join(seq)
L = len(ref)

# young full-length elements on this chromosome (sources for insertions; their spans are not used as insertion sites)
young = [(int(c[1]), int(c[2]), c[3]) for c in (l.split("\t") for l in open(rmsk))
         if c[0] == name and c[3] in YOUNG and int(c[2]) - int(c[1]) >= YOUNG[c[3]]]
blocked = sorted((s, e) for s, e, _ in young)


def mutate_copy(s):
    s = list(s)
    for i in range(len(s)):
        if s[i] in BASES and rng.random() < DIVERGE:
            s[i] = rng.choice([b for b in BASES if b != s[i]])
    return "".join(s)


def in_blocked(p):
    lo, hi = 0, len(blocked)
    while lo < hi:
        m = (lo + hi) // 2
        if blocked[m][1] <= p:
            lo = m + 1
        else:
            hi = m
    return lo < len(blocked) and blocked[lo][0] <= p


# events at reference positions: dict pos -> event (one per position)
events = {}
indel_rate = snp_rate / 10
n_snp, n_indel = int(L * snp_rate), int(L * indel_rate)
for _ in range(n_snp):
    p = rng.randrange(L)
    if ref[p] in BASES and p not in events:
        events[p] = ("SNP", rng.choice([b for b in BASES if b != ref[p]]))
for _ in range(n_indel):
    p = rng.randrange(L - 20)
    if p in events or ref[p] not in BASES:
        continue
    k = 1
    while k < 10 and rng.random() > 0.4:
        k += 1
    if rng.random() < 0.5:
        events[p] = ("INS", "".join(rng.choice(BASES) for _ in range(k)))
    elif "N" not in ref[p:p + k]:
        events[p] = ("DEL", k)
fam_count = {}
for _ in range(n_ins):
    if not young:
        break
    s, e, fam = rng.choice(young)
    for _try in range(100):
        p = rng.randrange(1000, L - 1000)
        if p not in events and ref[p] in BASES and not in_blocked(p) and "N" not in ref[p - 500:p + 500]:
            break
    events[p] = ("TE", mutate_copy(ref[s:e]), fam)
    fam_count[fam] = fam_count.get(fam, 0) + 1

# build mutated sequence and event table
outseq, rows, i, mpos = [], [], 0, 0
for p in sorted(events):
    if p < i:                       # inside a previous deletion
        continue
    outseq.append(ref[i:p])
    mpos += p - i
    ev = events[p]
    if ev[0] == "SNP":
        outseq.append(ev[1]); rows.append((mpos, p, "SNP", 0, "")); mpos += 1; i = p + 1
    elif ev[0] == "INS":
        outseq.append(ref[p] + ev[1]); rows.append((mpos + 1, p + 1, "INS", len(ev[1]), "")); mpos += 1 + len(ev[1]); i = p + 1
    elif ev[0] == "DEL":
        outseq.append(ref[p]); rows.append((mpos + 1, p + 1, "DEL", -ev[1], "")); mpos += 1; i = p + 1 + ev[1]
    else:                           # TE insertion after base p
        outseq.append(ref[p] + ev[1]); rows.append((mpos + 1, p + 1, "TE", len(ev[1]), ev[2])); mpos += 1 + len(ev[1]); i = p + 1
outseq.append(ref[i:])
mut = "".join(outseq)
with open(out + ".fa", "w") as o:
    o.write(f">{name}\n")
    for j in range(0, len(mut), 60):
        o.write(mut[j:j + 60] + "\n")
with open(out + ".events.tsv", "w") as o:
    o.write("mut_pos\tref_pos\tkind\tlength\tfamily\n")
    for r in rows:
        o.write("\t".join(map(str, r)) + "\n")
kinds = {}
for r in rows:
    kinds[r[2]] = kinds.get(r[2], 0) + 1
print(f"{name}: ref {L:,} bp, mutated {len(mut):,} bp; events {kinds}; TE copies {fam_count}")
