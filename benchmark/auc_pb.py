"""Which computational metric predicts the lab result? Read-only. argv[1] = csv path.

For every numeric computational metric, join it to the experimental `binding`
label on the same protein and compute a rank-based ROC AUC.
AUC 0.50 = coin flip. Below 0.50 = predicts the wrong way.
"""
import csv
import json
import sys
from collections import Counter, defaultdict

csv.field_size_limit(min(sys.maxsize, 2**31 - 1))
path = sys.argv[1]

comp = defaultdict(dict)        # metric -> {pid: value}
label = {}                      # pid -> 0/1
kd = {}                         # pid -> median kd (nM)
target_of = {}                  # pid -> target slug
strength = {}
raw_binding = Counter()


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
        pos = tot = 0
        kds = []
        for ev in evs:
            t, m, v = ev.get("type"), ev.get("metric"), ev.get("value")
            if t == "computational":
                x = num(v)
                if x is not None:
                    comp[m][pid] = x
            elif t == "experimental":
                if m == "binding":
                    raw_binding[str(v)] += 1
                    s = str(v).strip().lower()
                    tot += 1
                    if s in ("true", "1", "yes"):
                        pos += 1
                elif m == "kd":
                    x = num(v)
                    if x is not None:
                        kds.append(x)
                elif m == "binding_strength":
                    strength[pid] = str(v).strip().lower()
                if ev.get("target"):
                    target_of[pid] = ev["target"]
        if tot:
            label[pid] = 1 if pos * 2 > tot else 0
        if kds:
            kds.sort()
            kd[pid] = kds[len(kds) // 2]

print("raw `binding` values seen:", dict(raw_binding.most_common(6)))
n1 = sum(label.values())
print("labelled proteins: %d   binders: %d (%.1f%%)   non-binders: %d"
      % (len(label), n1, 100.0 * n1 / max(1, len(label)), len(label) - n1))
print("binding_strength labels:", dict(Counter(strength.values()).most_common(8)))
print("proteins with a KD:", len(kd))
print()


def auc(pairs):
    """pairs = [(value, label)]. Rank-based AUC with tie handling."""
    pairs = sorted(pairs)
    n = len(pairs)
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and pairs[j + 1][0] == pairs[i][0]:
            j += 1
        r = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[k] = r
        i = j + 1
    pos = [ranks[k] for k in range(n) if pairs[k][1] == 1]
    np_, nn = len(pos), n - len(pos)
    if np_ == 0 or nn == 0:
        return None
    return (sum(pos) - np_ * (np_ + 1) / 2.0) / (np_ * nn)


rows = []
for m, d in comp.items():
    pairs = [(v, label[p]) for p, v in d.items() if p in label]
    if len(pairs) < 200:
        continue
    npos = sum(l for _, l in pairs)
    if npos < 30 or npos == len(pairs):
        continue
    a = auc(pairs)
    if a is None:
        continue
    rows.append((abs(a - 0.5), a, m, len(pairs), npos))

rows.sort(reverse=True)
print("=== ROC AUC of each computational metric against the measured binding outcome ===")
print("%-42s %6s %8s %8s" % ("metric", "AUC", "n", "binders"))
for _, a, m, n, npos in rows[:28]:
    flag = "  <-- predicts BACKWARDS" if a < 0.5 else ""
    print("%-42s %6.3f %8d %8d%s" % (m, a, n, npos, flag))

print()
print("=== the same, restricted to the single biggest target ===")
tc = Counter(target_of[p] for p in label if p in target_of)
for tgt, _ in tc.most_common(3):
    sub = {p for p in label if target_of.get(p) == tgt}
    npos_t = sum(label[p] for p in sub)
    print("\n-- %s: %d labelled, %d binders (%.1f%%)"
          % (tgt, len(sub), npos_t, 100.0 * npos_t / max(1, len(sub))))
    out = []
    for m, d in comp.items():
        pairs = [(v, label[p]) for p, v in d.items() if p in sub]
        npos = sum(l for _, l in pairs)
        if len(pairs) < 100 or npos < 15 or npos == len(pairs):
            continue
        a = auc(pairs)
        if a is not None:
            out.append((abs(a - 0.5), a, m, len(pairs)))
    out.sort(reverse=True)
    for _, a, m, n in out[:6]:
        print("     %-38s AUC %5.3f  (n=%d)" % (m, a, n))
