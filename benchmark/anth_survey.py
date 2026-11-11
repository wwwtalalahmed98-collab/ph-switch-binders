"""Head-to-head: ten structure predictors scored on the SAME designs with dual-CRO labels.

Also counts the off-target counter-screens, which are rare in public data.
argv[1] = design_summary.csv
"""
import csv
import math
import sys
from collections import Counter, defaultdict

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
rows = list(csv.DictReader(open(sys.argv[1], newline="", encoding="utf-8")))


def num(v):
    try:
        if v is None or v == "":
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


def rank(xs):
    o = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and xs[o[j + 1]] == xs[o[i]]:
            j += 1
        a = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            r[o[k]] = a
        i = j + 1
    return r


def spearman(xs, ys):
    n = len(xs)
    if n < 20:
        return None
    rx, ry = rank(xs), rank(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    sxx = sum((a - mx) ** 2 for a in rx)
    syy = sum((b - my) ** 2 for b in ry)
    if sxx == 0 or syy == 0:
        return None
    return sxy / math.sqrt(sxx * syy)


print("rows:", len(rows))
print("design_model:", dict(Counter(r["design_model"] for r in rows)))
print("generator:", dict(Counter(r["generator"] for r in rows).most_common(8)))
print("binder_final:", dict(Counter(r["binder_final"] for r in rows)))
print("vendor_agreement:", dict(Counter(r["vendor_agreement"] for r in rows)))
tgt = Counter(r["target"] for r in rows)
print("targets (%d):" % len(tgt), tgt.most_common())
print()

lab = {}
for r in rows:
    b = (r["binder_final"] or "").strip().lower()
    if b in ("true", "yes", "1", "binder"):
        lab[r["uuid"]] = 1
    elif b in ("false", "no", "0", "non-binder", "nonbinder"):
        lab[r["uuid"]] = 0
print("usable binary labels: %d  binders: %d (%.1f%%)"
      % (len(lab), sum(lab.values()), 100.0 * sum(lab.values()) / max(1, len(lab))))

preds = [c[len("ipsae_min_"):] for c in rows[0] if c.startswith("ipsae_min_")]
print("\n=== ipSAE_min head-to-head, identical designs, integrated dual-CRO label ===")
print("%-10s %8s %8s %10s" % ("predictor", "AUC", "n", "coverage"))
out = []
for p in preds:
    pr = [(num(r["ipsae_min_" + p]), lab[r["uuid"]]) for r in rows
          if r["uuid"] in lab and num(r["ipsae_min_" + p]) is not None]
    a = auc(pr)
    if a is not None:
        out.append((a, p, len(pr)))
for a, p, n in sorted(out, reverse=True):
    print("%-10s %8.3f %8d %9.0f%%" % (p, a, n, 100.0 * n / len(lab)))

print("\n=== sc_dockq head-to-head ===")
out2 = []
for p in preds:
    c = "sc_dockq_" + p
    if c not in rows[0]:
        continue
    pr = [(num(r[c]), lab[r["uuid"]]) for r in rows
          if r["uuid"] in lab and num(r[c]) is not None]
    a = auc(pr)
    if a is not None:
        out2.append((a, p, len(pr)))
for a, p, n in sorted(out2, reverse=True):
    print("%-10s %8.3f %8d" % (p, a, n))

kd = {}
for r in rows:
    v = num(r["kd_nM_final"])
    if v and v > 0 and lab.get(r["uuid"]) == 1:
        kd[r["uuid"]] = math.log10(v)
print("\n=== vs affinity: Spearman rho of ipSAE_min against log10 KD, binders only ===")
print("designs with a final KD and a binder call: %d" % len(kd))
out3 = []
for p in preds:
    pr = [(num(r["ipsae_min_" + p]), kd[r["uuid"]]) for r in rows
          if r["uuid"] in kd and num(r["ipsae_min_" + p]) is not None]
    if len(pr) < 30:
        continue
    rho = spearman([a for a, _ in pr], [b for _, b in pr])
    if rho is not None:
        out3.append((abs(rho), rho, p, len(pr)))
for _, rho, p, n in sorted(out3, reverse=True):
    print("%-10s rho %7.3f  (n=%d)" % (p, rho, n))

print("\n=== counter-screens: measured off-target behaviour ===")
for c in ("twist_mouse_binding", "twist_cyno_binding", "twist_clec12a_binding",
          "twist_gdf11_binding", "twist_cas9_apo_binding", "mouse_binding_final",
          "cyno_binding_final"):
    if c in rows[0]:
        cc = Counter((r[c] or "").strip() or "(blank)" for r in rows)
        print("  %-24s %s" % (c, dict(cc.most_common(5))))

on = [r for r in rows if lab.get(r["uuid"]) == 1]
clec = [r for r in on if (r["twist_clec12a_binding"] or "").strip().lower()
        in ("true", "yes", "1")]
print("\n  confirmed binders that ALSO bound the unrelated control CLEC12A: %d of %d"
      % (len(clec), len(on)))
