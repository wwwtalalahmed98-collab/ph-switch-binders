"""Which of 433 computed metrics predicts each experimental gate?

Three gates in order: the DNA could be synthesised, the construct was recovered
in the screen, the CAR drove CD20-specific T-cell proliferation. A metric that
only predicts the first two is predicting manufacturability, not function.

Usage: python -I benchmark/b2b_gates.py <12k_all_metrics.csv> <12k_all_results.csv>
"""
import csv, sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

met = list(csv.DictReader(open(sys.argv[1], newline="", encoding="utf-8", errors="replace")))
res = {r["global_id"]: r for r in
       csv.DictReader(open(sys.argv[2], newline="", encoding="utf-8", errors="replace"))}

def num(v):
    try:
        if v in (None, "", "NA", "nan", "None"): return np.nan
        return float(v)
    except Exception: return np.nan

cols = [c for c in met[0] if c not in ("team","global_id","sequence","dna_sequence",
                                       "dssp","longest_dup_substr","longest_dup_substr_dna")]

def auc(x, y):
    ok = ~np.isnan(x)
    x, y = x[ok], y[ok]
    if y.sum() < 10 or (1-y).sum() < 10: return np.nan, 0
    order = np.argsort(x); r = np.empty(len(x)); r[order] = np.arange(1, len(x)+1)
    # average ranks for ties
    s = np.argsort(x); xs = x[s]; i = 0
    while i < len(xs):
        j = i
        while j+1 < len(xs) and xs[j+1] == xs[i]: j += 1
        if j > i: r[s[i:j+1]] = (i+j+2)/2
        i = j+1
    n1, n0 = y.sum(), (1-y).sum()
    return (r[y==1].sum() - n1*(n1+1)/2) / (n1*n0), int(ok.sum())

gates = {}
ids = [m["global_id"] for m in met]
rr = [res.get(i, {}) for i in ids]
gates["DNA synthesised"]      = np.array([1.0 if r.get("leah_12k_twist_dna_detected")=="True" else 0.0 for r in rr])
gates["recovered in screen"]  = np.array([1.0 if r.get("leah_12k_detected")=="True" else 0.0 for r in rr])
func = np.array([r.get("fourclass","") for r in rr])
mask_fn = (func=="Up") | (func=="Not Sig") | (func=="Down")
gates["T-cell function (Up)"] = np.array([1.0 if f=="Up" else 0.0 for f in func])

for name, y in gates.items():
    if name == "T-cell function (Up)":
        sub = mask_fn
    else:
        sub = np.ones(len(y), bool)
    print("\n=== %s ===  positives %d / %d" % (name, int(y[sub].sum()), int(sub.sum())))
    out = []
    for c in cols:
        x = np.array([num(m[c]) for m in met])
        a, n = auc(x[sub], y[sub])
        if not np.isnan(a): out.append((abs(a-0.5), a, c, n))
    out.sort(reverse=True)
    for _, a, c, n in out[:8]:
        print("   %-52s AUC %.3f   (n=%d)" % (c[:52], a, n))
