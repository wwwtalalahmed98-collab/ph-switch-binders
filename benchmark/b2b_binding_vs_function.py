"""Does measured binding predict measured FUNCTION?

Bits to Binders: 12,000 CD20-targeting minibinders, each built into a chimeric
antigen receptor and screened in a pooled T-cell proliferation assay. A subset
also went to Adaptyv for SPR. So for ~100 designs we have both the field's
usual answer (does it bind?) and the answer that matters (does the cell kill?).

Usage: python -I benchmark/b2b_binding_vs_function.py <12k_all_results.csv> <b2b_summary.csv>
"""
import csv, sys, collections
import numpy as np
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

res = list(csv.DictReader(open(sys.argv[1], newline="", encoding="utf-8", errors="replace")))
spr = list(csv.DictReader(open(sys.argv[2], newline="", encoding="utf-8", errors="replace")))

print("results rows", len(res), " spr rows", len(spr))

# SPR is per-replicate; collapse to one call per sequence
by_seq = collections.defaultdict(list)
for r in spr:
    by_seq[r["sequence"].strip()].append(r)
print("unique SPR sequences:", len(by_seq))

def call(reps):
    b = [x["binding"] for x in reps]
    if "true" in b: return "binds"
    if all(x == "unknown" for x in b): return "no_expression"
    return "no_binding"

spr_call = {s: call(v) for s, v in by_seq.items()}
print(collections.Counter(spr_call.values()))

idx = {r["sequence"].strip(): r for r in res}
matched = [(s, c, idx[s]) for s, c in spr_call.items() if s in idx]
print("matched to 12k table:", len(matched), "of", len(spr_call))

tab = collections.Counter((c, r["fourclass"]) for _, c, r in matched)
classes = ["Up", "Not Sig", "Down", "Not Recovered"]
print("\n          SPR call  |  " + "  ".join("%13s" % c for c in classes) + "   total")
for c in ["binds", "no_binding", "no_expression"]:
    row = [tab[(c, k)] for k in classes]
    print("%18s  |  " % c + "  ".join("%13d" % v for v in row) + "   %5d" % sum(row))
