"""Affinity is two numbers, not one. Is either half predictable?

KD = koff / kon. A design's affinity can be poor because it never finds the
target (slow kon) or because it falls off (fast koff). These are different
physics and the field almost always reports only KD.

If one half is predictable and the other is not, that is a real finding and it
says where any future model should aim.

Usage: python -I benchmark/kinetics.py <anthropic design_summary.csv>
"""
import csv
import math
import sys

import numpy as np
from scipy.stats import spearmanr

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
rows = list(csv.DictReader(open(sys.argv[1], newline="", encoding="utf-8")))


def num(v):
    try:
        if v in (None, "", "NA", "nan", "not_tested", "inconclusive"):
            return None
        return float(v)
    except Exception:
        return None


preds = sorted(c[len("ipsae_min_"):] for c in rows[0] if c.startswith("ipsae_min_"))

for vendor in ("adaptyv", "twist"):
    kdc, konc, koffc = ("%s_kd_nM" % vendor, "%s_kon_per_Ms" % vendor,
                        "%s_koff_per_s" % vendor)
    if kdc not in rows[0]:
        continue
    recs = []
    for r in rows:
        if (r.get("binder_final") or "").strip().lower() != "true":
            continue
        kd, kon, koff = num(r.get(kdc)), num(r.get(konc)), num(r.get(koffc))
        if kd and kd > 0 and kon and kon > 0 and koff and koff > 0:
            recs.append((r, math.log10(kd), math.log10(kon), math.log10(koff)))
    if len(recs) < 40:
        print("\n%s: only %d binders with full kinetics - skipped"
              % (vendor, len(recs)))
        continue
    print("\n=== %s: %d confirmed binders with kon, koff and KD ===" % (vendor,
                                                                        len(recs)))
    lkd = np.array([x[1] for x in recs])
    lkon = np.array([x[2] for x in recs])
    lkoff = np.array([x[3] for x in recs])
    print("   spread  log10 KD %.2f   log10 kon %.2f   log10 koff %.2f  (SD)"
          % (lkd.std(), lkon.std(), lkoff.std()))
    r_on = spearmanr(lkon, lkd).statistic
    r_off = spearmanr(lkoff, lkd).statistic
    print("   which half drives affinity here?  kon vs KD rho %+.2f"
          "   koff vs KD rho %+.2f" % (r_on, r_off))
    print("\n   %-12s %10s %10s %10s" % ("predictor", "vs KD", "vs kon", "vs koff"))
    best = []
    for p in preds:
        v = np.array([num(x[0].get("ipsae_min_" + p)) for x in recs], dtype=float)
        ok = ~np.isnan(v)
        if ok.sum() < 40:
            continue
        a = spearmanr(v[ok], lkd[ok]).statistic
        b = spearmanr(v[ok], lkon[ok]).statistic
        c = spearmanr(v[ok], lkoff[ok]).statistic
        best.append((abs(c), p, a, b, c))
        print("   %-12s %10.3f %10.3f %10.3f" % (p, a, b, c))
    if best:
        best.sort(reverse=True)
        _, p, a, b, c = best[0]
        print("\n   strongest koff signal: %s at rho %+.3f "
              "(vs %+.3f for KD)" % (p, c, a))
