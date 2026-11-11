"""Do interface metrics add anything ONCE the protein is known to express?

A design that never expressed cannot bind, so a metric that merely predicts
expression will look like it predicts binding. This restricts the comparison to
proteins that DID express, which is the only fair test of an interface metric.
argv[1] = csv path.
"""
import csv
import json
import sys
from collections import Counter, defaultdict

csv.field_size_limit(min(sys.maxsize, 2**31 - 1))
path = sys.argv[1]

comp = defaultdict(dict)
label, expressed, target_of = {}, {}, {}


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
        bpos = btot = epos = etot = 0
        for ev in evs:
            t, m, v = ev.get("type"), ev.get("metric"), ev.get("value")
            if t == "computational":
                x = num(v)
                if x is not None:
                    comp[m][pid] = x
            elif t == "experimental":
                s = str(v).strip().lower()
                if m == "binding":
                    btot += 1
                    bpos += s in ("true", "1", "yes")
                elif m == "expressed":
                    etot += 1
                    epos += s in ("true", "1", "yes")
                if ev.get("target"):
                    target_of[pid] = ev["target"]
        if btot:
            label[pid] = 1 if bpos * 2 > btot else 0
        if etot:
            expressed[pid] = 1 if epos * 2 > etot else 0


def auc(pairs):
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
    a, b = len(pos), n - len(pos)
    if a == 0 or b == 0:
        return None
    return (sum(pos) - a * (a + 1) / 2.0) / (a * b)


exp_ok = {p for p, v in expressed.items() if v == 1}
print("proteins with an expression result:", len(expressed))
print("  of which expressed:", len(exp_ok))
nb = sum(label[p] for p in exp_ok if p in label)
nl = sum(1 for p in exp_ok if p in label)
print("expressed AND labelled for binding: %d, of which binders %d (%.1f%%)"
      % (nl, nb, 100.0 * nb / max(1, nl)))
print()

print("=== AUC against EXPRESSION (all proteins) vs AUC against BINDING (expressed only) ===")
print("%-42s %9s %9s %9s" % ("metric", "->express", "->bind|exp", "n(exp)"))
out = []
for m, d in comp.items():
    pe = [(v, expressed[p]) for p, v in d.items() if p in expressed]
    pb = [(v, label[p]) for p, v in d.items() if p in exp_ok and p in label]
    if len(pe) < 300 or len(pb) < 200:
        continue
    npos = sum(l for _, l in pb)
    if npos < 30 or npos == len(pb):
        continue
    ae, ab = auc(pe), auc(pb)
    if ae is None or ab is None:
        continue
    out.append((abs(ab - 0.5), m, ae, ab, len(pb)))
out.sort(reverse=True)
for _, m, ae, ab, n in out[:24]:
    print("%-42s %9.3f %9.3f %9d" % (m, ae, ab, n))

print()
print("=== per-target: does ipSAE separate binders among EXPRESSED designs? ===")
tc = Counter(target_of[p] for p in exp_ok if p in label and p in target_of)
print("%-32s %6s %7s %9s %9s %9s" % ("target", "n", "binders", "ipsae", "iptm", "plddt"))
for tgt, n in tc.most_common(10):
    sub = {p for p in exp_ok if p in label and target_of.get(p) == tgt}
    npos = sum(label[p] for p in sub)
    if npos < 10 or npos == len(sub):
        continue
    cells = []
    for m in ("boltz2_ipsae", "boltz2_iptm", "boltz2_complex_plddt"):
        pairs = [(v, label[p]) for p, v in comp.get(m, {}).items() if p in sub]
        a = auc(pairs) if len(pairs) >= 40 else None
        cells.append("%9.3f" % a if a is not None else "        -")
    print("%-32s %6d %7d %s" % (tgt[:32], len(sub), npos, "".join(cells)))
