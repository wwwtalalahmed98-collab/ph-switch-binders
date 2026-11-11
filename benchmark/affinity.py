"""Among designs that DO bind, does any computational metric predict HOW tightly?

Binary binding is the easy question. Affinity is the one the field says is open.
Spearman rank correlation of each metric against log10(KD), binders only.
argv[1] = proteinbase csv.
"""
import csv
import json
import math
import sys
from collections import Counter, defaultdict

csv.field_size_limit(min(sys.maxsize, 2**31 - 1))
path = sys.argv[1]

comp = defaultdict(dict)
kd, target_of = {}, {}


def num(v):
    try:
        return float(v)
    except Exception:
        return None


with open(path, newline="", encoding="utf-8") as fh:
    for row in csv.DictReader(fh):
        pid = list(row.values())[0]
        try:
            evs = json.loads(row.get("evaluations") or "[]")
        except Exception:
            continue
        kds = []
        for ev in evs:
            t, m, v = ev.get("type"), ev.get("metric"), ev.get("value")
            if t == "computational":
                x = num(v)
                if x is not None:
                    comp[m][pid] = x
            elif t == "experimental" and m == "kd":
                x = num(v)
                if x is not None and x > 0:
                    kds.append(x)
            if ev.get("target"):
                target_of[pid] = ev["target"]
        if kds:
            kds.sort()
            kd[pid] = math.log10(kds[len(kds) // 2])


def rank(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            r[order[k]] = avg
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


print("designs with a measured KD: %d" % len(kd))
vals = sorted(kd.values())
print("log10 KD (nM units as stored): min %.2f  median %.2f  max %.2f"
      % (vals[0], vals[len(vals) // 2], vals[-1]))
tc = Counter(target_of[p] for p in kd if p in target_of)
print("targets represented: %d; top:" % len(tc), tc.most_common(6))
print()
print("=== Spearman rho vs log10(KD), binders only ===")
print("   (negative rho = higher metric goes with tighter binding, which is what")
print("    you would want from a score that is supposed to rank affinity)")
print("%-42s %8s %8s" % ("metric", "rho", "n"))
out = []
for m, d in comp.items():
    pairs = [(v, kd[p]) for p, v in d.items() if p in kd]
    if len(pairs) < 100:
        continue
    r = spearman([a for a, _ in pairs], [b for _, b in pairs])
    if r is not None:
        out.append((abs(r), r, m, len(pairs)))
out.sort(reverse=True)
for _, r, m, n in out[:22]:
    print("%-42s %8.3f %8d" % (m, r, n))

print()
print("=== within the three best-represented targets ===")
for tgt, _ in tc.most_common(3):
    sub = {p for p in kd if target_of.get(p) == tgt}
    if len(sub) < 40:
        continue
    print("\n-- %s (n=%d with KD)" % (tgt, len(sub)))
    rows = []
    for m, d in comp.items():
        pairs = [(v, kd[p]) for p, v in d.items() if p in sub]
        if len(pairs) < 40:
            continue
        r = spearman([a for a, _ in pairs], [b for _, b in pairs])
        if r is not None:
            rows.append((abs(r), r, m, len(pairs)))
    rows.sort(reverse=True)
    for _, r, m, n in rows[:5]:
        print("     %-38s rho %6.3f  (n=%d)" % (m, r, n))
