"""Is the benchmark number everybody reports measuring the wrong thing?

Published AUCs for binder-selection metrics are POOLED across targets. But in
real use you only ever rank designs WITHIN one target - you have one target and
a pile of candidates for it. If targets differ in their average score, pooling
inflates the number without the metric being able to rank anything.

This computes, for the same metric and the same data:
    pooled AUC        - what papers report
    within-target AUC - what a designer actually experiences

Usage: python -I benchmark/pooled_vs_within.py anthropic <design_summary.csv>
       python -I benchmark/pooled_vs_within.py dtu       <final_dataset.csv>
"""
import csv
import sys

import numpy as np
from sklearn.metrics import roc_auc_score

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
csv.field_size_limit(2 ** 31 - 1)

MODE, PATH = sys.argv[1], sys.argv[2]


def num(v):
    try:
        if v in (None, "", "NA", "nan", "not_tested", "inconclusive"):
            return None
        return float(v)
    except Exception:
        return None


rows = list(csv.DictReader(open(PATH, newline="", encoding="utf-8",
                                errors="replace")))

if MODE == "anthropic":
    metrics = [c for c in rows[0] if c.startswith("ipsae_min_")
               or c.startswith("sc_dockq_")]
    FLIP = set()

    def lab(r):
        b = (r.get("binder_final") or "").strip().lower()
        return 1 if b == "true" else (0 if b == "false" else None)

    def tgt(r):
        return r.get("target")
else:
    metrics = ["af3_ipSAE_min", "boltz1_ipSAE_min", "colab_ipSAE_min",
               "af3_ipae", "boltz1_ipae", "af2_pae_interaction",
               "af3_iptm_model_0", "boltz1_iptm_model_0", "af3_dockQ",
               "boltz1_complex_iplddt_model_0", "af2_plddt_binder"]
    FLIP = {"af3_ipae", "boltz1_ipae", "af2_pae_interaction"}

    def lab(r):
        b = (r.get("binder") or "").strip()
        return 1 if b == "True" else (0 if b == "False" else None)

    def tgt(r):
        return r.get("target_id")


data = []
for r in rows:
    L = lab(r)
    if L is None:
        continue
    data.append((tgt(r), L, r))

print("dataset: %s   rows: %d   binders: %d (%.1f%%)   targets: %d"
      % (MODE, len(data), sum(d[1] for d in data),
         100.0 * sum(d[1] for d in data) / len(data),
         len(set(d[0] for d in data))))
print("\n%-32s %8s %10s %9s %8s" % ("metric", "pooled", "within-tgt", "drop",
                                    "targets"))
print("-" * 72)

summary = []
for m in metrics:
    ys, xs, gs = [], [], []
    for t, L, r in data:
        v = num(r.get(m))
        if v is None:
            continue
        xs.append(-v if m in FLIP else v)
        ys.append(L)
        gs.append(t)
    if len(ys) < 200:
        continue
    xs, ys, gs = np.array(xs), np.array(ys), np.array(gs)
    try:
        pooled = roc_auc_score(ys, xs)
    except ValueError:
        continue
    per, wts = [], []
    for t in set(gs):
        sel = gs == t
        if ys[sel].sum() < 5 or ys[sel].sum() == sel.sum():
            continue
        per.append(roc_auc_score(ys[sel], xs[sel]))
        wts.append(sel.sum())
    if not per:
        continue
    within = float(np.average(per, weights=wts))
    summary.append((pooled, within, m, len(per), per))
    print("%-32s %8.3f %10.3f %9.3f %8d" % (m, pooled, within, pooled - within,
                                            len(per)))

if summary:
    print("-" * 72)
    p = np.mean([s[0] for s in summary])
    w = np.mean([s[1] for s in summary])
    print("%-32s %8.3f %10.3f %9.3f" % ("MEAN ACROSS METRICS", p, w, p - w))
    best = max(summary)
    print("\nbest pooled metric: %s" % best[2])
    print("  its per-target AUCs: %s"
          % ", ".join("%.2f" % v for v in sorted(best[4])))
    below = sum(1 for v in best[4] if v < 0.6)
    print("  targets where it is below 0.60 (near useless): %d of %d"
          % (below, len(best[4])))
    worse = sum(1 for v in best[4] if v < 0.5)
    print("  targets where it points BACKWARDS (<0.50): %d of %d"
          % (worse, len(best[4])))
