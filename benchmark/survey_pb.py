"""Survey the Proteinbase export. Read-only. Path passed as argv[1]."""
import csv
import json
import sys
from collections import Counter

csv.field_size_limit(min(sys.maxsize, 2**31 - 1))

path = sys.argv[1]
rows = 0
cols = None
design = Counter()
ev_type = Counter()
metric_by_type = Counter()
method_by_metric = Counter()
targets = Counter()
has_exp = 0
has_comp = 0
has_both = 0
kd_vals = []
protein_metrics = Counter()
exp_metric_targets = Counter()

with open(path, newline="", encoding="utf-8") as fh:
    r = csv.DictReader(fh)
    cols = r.fieldnames
    for row in r:
        rows += 1
        design[(row.get("designMethod") or "(blank)")] += 1
        try:
            evs = json.loads(row.get("evaluations") or "[]")
        except Exception:
            evs = []
        e = c = False
        seen = set()
        for ev in evs:
            t = ev.get("type", "")
            m = ev.get("metric", "")
            ev_type[t] += 1
            metric_by_type[(t, m)] += 1
            method_by_metric[(m, ev.get("method", ""))] += 1
            if ev.get("target"):
                targets[ev["target"]] += 1
            if t == "experimental":
                e = True
                exp_metric_targets[m] += 1
                if m == "kd":
                    try:
                        kd_vals.append(float(ev.get("value")))
                    except Exception:
                        pass
            elif t == "computational":
                c = True
            seen.add((t, m))
        for s in seen:
            protein_metrics[s] += 1
        has_exp += e
        has_comp += c
        has_both += (e and c)

print("columns:", cols)
print("rows (proteins):", rows)
print()
print("proteins with >=1 experimental measurement:", has_exp)
print("proteins with >=1 computational metric:   ", has_comp)
print("proteins with BOTH:                       ", has_both)
print()
print("evaluation rows by type:", dict(ev_type))
print()
print("-- top experimental metrics (evaluation count) --")
for (t, m), n in metric_by_type.most_common():
    if t == "experimental":
        print("  %-28s %7d   (on %d proteins)" % (m, n, protein_metrics[(t, m)]))
print()
print("-- top computational metrics (on N proteins) --")
comp = [((t, m), n) for (t, m), n in metric_by_type.items() if t == "computational"]
comp.sort(key=lambda x: -protein_metrics[x[0]])
for (t, m), n in comp[:30]:
    print("  %-28s %7d   (on %d proteins)" % (m, n, protein_metrics[(t, m)]))
print()
print("-- curated metrics --")
for (t, m), n in metric_by_type.most_common():
    if t == "curated":
        print("  %-28s %7d   (on %d proteins)" % (m, n, protein_metrics[(t, m)]))
print()
print("-- targets (evaluation count), top 25 --")
for k, n in targets.most_common(25):
    print("  %-40s %7d" % (k, n))
print("  distinct targets:", len(targets))
print()
print("-- design methods, top 25 --")
for k, n in design.most_common(25):
    print("  %-40s %7d" % (k, n))
print("  distinct design methods:", len(design))
print()
if kd_vals:
    kd_vals.sort()
    n = len(kd_vals)
    print("kd measurements: %d  min=%.3g  median=%.3g  max=%.3g"
          % (n, kd_vals[0], kd_vals[n // 2], kd_vals[-1]))
