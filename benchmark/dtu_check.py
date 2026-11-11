"""Does the DTU meta-analysis conclusion survive its own target imbalance?

2,123 of 3,676 rows are a single target (FGFR2). This recomputes the headline
metric comparison with and without it, and per target.
argv[1] = final_dataset.csv
"""
import csv
import sys
from collections import Counter, defaultdict

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
csv.field_size_limit(2 ** 31 - 1)
rows = list(csv.DictReader(open(sys.argv[1], newline="", encoding="utf-8",
                                errors="replace")))


def num(v):
    try:
        if v in (None, "", "NA", "nan"):
            return None
        return float(v)
    except Exception:
        return None


def auc(pairs):
    pairs = sorted(pairs)
    n = len(pairs)
    r = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and pairs[j + 1][0] == pairs[i][0]:
            j += 1
        a = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            r[k] = a
        i = j + 1
    pos = [r[k] for k in range(n) if pairs[k][1] == 1]
    p, q = len(pos), n - len(pos)
    if p == 0 or q == 0:
        return None
    return (sum(pos) - p * (p + 1) / 2.0) / (p * q)


METRICS = ["af3_ipSAE_min", "boltz1_ipSAE_min", "colab_ipSAE_min",
           "af3_ipae", "boltz1_ipae", "colab_ipae", "af2_pae_interaction",
           "af3_iptm_model_0", "boltz1_iptm_model_0", "colab_iptm_model_0",
           "af2_dockQ", "af3_dockQ", "boltz1_dockQ",
           "boltz1_complex_iplddt_model_0", "af2_plddt_binder"]
# lower-is-better metrics must be flipped so AUC is comparable
FLIP = {"af3_ipae", "boltz1_ipae", "colab_ipae", "af2_pae_interaction"}

lab, tgt = {}, {}
for i, r in enumerate(rows):
    b = (r.get("binder") or "").strip()
    if b in ("True", "False"):
        lab[i] = 1 if b == "True" else 0
        tgt[i] = r.get("target_id") or "(none)"

n1 = sum(lab.values())
print("labelled rows: %d   binders: %d (%.1f%%)" % (len(lab), n1,
                                                    100.0 * n1 / len(lab)))
tc = Counter(tgt[i] for i in lab)
print("targets:", tc.most_common())
print()


def table(idx, title):
    nb = sum(lab[i] for i in idx)
    print("=== %s  (n=%d, binders=%d, %.1f%%) ===" % (title, len(idx), nb,
                                                      100.0 * nb / max(1, len(idx))))
    out = []
    for m in METRICS:
        if m not in rows[0]:
            continue
        pr = []
        for i in idx:
            v = num(rows[i].get(m))
            if v is not None:
                pr.append((-v if m in FLIP else v, lab[i]))
        if len(pr) < 100:
            continue
        a = auc(pr)
        if a is not None:
            out.append((a, m, len(pr)))
    for a, m, n in sorted(out, reverse=True):
        print("   %-32s %6.3f   n=%d" % (m, a, n))
    print()


allidx = list(lab)
table(allidx, "ALL rows, as the paper reports")
table([i for i in allidx if tgt[i] == "FGFR2"], "FGFR2 only (the dominant target)")
table([i for i in allidx if tgt[i] != "FGFR2"], "EVERYTHING EXCEPT FGFR2")

print("=== per target, af3_ipSAE_min only ===")
print("%-18s %7s %8s %8s" % ("target", "n", "binders", "AUC"))
for t, n in tc.most_common():
    idx = [i for i in lab if tgt[i] == t]
    nb = sum(lab[i] for i in idx)
    if nb < 5 or nb == len(idx):
        continue
    pr = [(num(rows[i].get("af3_ipSAE_min")), lab[i]) for i in idx
          if num(rows[i].get("af3_ipSAE_min")) is not None]
    if len(pr) < 30:
        continue
    a = auc(pr)
    print("%-18s %7d %8d %8s" % (t, len(idx), nb,
                                 "%.3f" % a if a is not None else "-"))
