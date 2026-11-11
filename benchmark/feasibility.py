"""FEASIBILITY TEST, not a product.

Two questions, answered with leave-one-target-out cross-validation so that no
model can cheat by memorising a target:

  A. Can a combination of available scores beat the BEST SINGLE score at
     separating binders from non-binders?
  B. Can anything predict HOW TIGHTLY a binder binds (log10 KD)?

Leave-one-target-out is the honest protocol here: in real use you point a
pipeline at a target nobody has data for. Random k-fold would leak.

Usage: python -I benchmark/feasibility.py <anthropic design_summary.csv>
"""
import csv
import math
import sys
import warnings

import numpy as np
from scipy.stats import spearmanr
from sklearn.ensemble import GradientBoostingClassifier, GradientBoostingRegressor
from sklearn.linear_model import LogisticRegression, RidgeCV
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

KD = "AVILMFWYC"
NEG, POS = "DE", "KR"


def num(v):
    try:
        if v in (None, "", "NA", "nan", "not_tested", "inconclusive"):
            return None
        return float(v)
    except Exception:
        return None


def seqfeat(s):
    s = (s or "").upper()
    n = max(1, len(s))
    f = [len(s),
         sum(c in KD for c in s) / n,
         sum(c in NEG for c in s) / n,
         sum(c in POS for c in s) / n,
         (sum(c in POS for c in s) - sum(c in NEG for c in s)) / n,
         s.count("G") / n, s.count("P") / n, s.count("C") / n]
    best = run = 0
    for c in s:
        run = run + 1 if c in KD else 0
        best = max(best, run)
    f.append(best)
    return f


rows = list(csv.DictReader(open(sys.argv[1], newline="", encoding="utf-8")))
preds = sorted(c[len("ipsae_min_"):] for c in rows[0] if c.startswith("ipsae_min_"))
cols = (["ipsae_min_" + p for p in preds] + ["sc_dockq_" + p for p in preds]
        + ["binder_length", "epitope_n_residues"])

X, y, kd, g = [], [], [], []
for r in rows:
    b = (r.get("binder_final") or "").strip().lower()
    if b not in ("true", "false"):
        continue
    vals = [num(r.get(c)) for c in cols]
    if sum(v is None for v in vals) > 4:
        continue
    med = 0.0
    vals = [med if v is None else v for v in vals]
    X.append(vals + seqfeat(r.get("sequence")))
    y.append(1 if b == "true" else 0)
    kd.append(num(r.get("kd_nM_final")))
    g.append(r.get("target"))

X = np.array(X, float)
y = np.array(y, int)
g = np.array(g)
names = cols + ["len", "hydroph", "neg", "pos", "netq", "gly", "pro", "cys", "hrun"]
print("designs: %d   binders: %d (%.1f%%)   targets: %d"
      % (len(y), y.sum(), 100.0 * y.mean(), len(set(g))))
print("features: %d" % X.shape[1])
print("protocol: leave-one-target-out\n")

# ---------- A. binary ----------
print("=== A. binder vs non-binder ===")
best_single, best_name = 0.5, None
for j, nm in enumerate(names):
    if not nm.startswith("ipsae_min_"):
        continue
    a = roc_auc_score(y, X[:, j])
    if a > best_single:
        best_single, best_name = a, nm
print("best single score (pooled, optimistic): %-22s AUC %.3f" % (best_name, best_single))

for label, mk in [("logistic regression", lambda: make_pipeline(
                       StandardScaler(), LogisticRegression(max_iter=4000, C=0.3))),
                  ("gradient boosting", lambda: GradientBoostingeClassifier()
                   if False else GradientBoostingClassifier(random_state=0,
                                                            n_estimators=150,
                                                            max_depth=2))]:
    oof = np.full(len(y), np.nan)
    for t in sorted(set(g)):
        tr, te = g != t, g == t
        if y[tr].sum() < 10 or y[tr].sum() == tr.sum():
            continue
        m = mk()
        m.fit(X[tr], y[tr])
        oof[te] = m.predict_proba(X[te])[:, 1]
    ok = ~np.isnan(oof)
    print("%-22s held-out AUC %.3f   (n=%d)" % (label, roc_auc_score(y[ok], oof[ok]),
                                                ok.sum()))

# same single score evaluated honestly, pooled over held-out targets only
print("\n   per-target AUC of the best single score, where it is computable:")
for t in sorted(set(g)):
    m = g == t
    if y[m].sum() < 5 or y[m].sum() == m.sum():
        continue
    j = names.index(best_name)
    print("     %-16s n=%4d binders=%3d  AUC %.3f"
          % (t, m.sum(), y[m].sum(), roc_auc_score(y[m], X[m, j])))

# ---------- B. affinity ----------
print("\n=== B. affinity among confirmed binders ===")
sel = np.array([(yy == 1 and k is not None and k > 0) for yy, k in zip(y, kd)])
Xa = X[sel]
ya = np.array([math.log10(k) for k, s in zip(kd, sel) if s])
ga = g[sel]
print("binders with a KD: %d across %d targets" % (len(ya), len(set(ga))))
bs, bn = 0.0, None
for j, nm in enumerate(names):
    r = spearmanr(Xa[:, j], ya).statistic
    if abs(r) > abs(bs):
        bs, bn = r, nm
print("best single feature (pooled, optimistic): %-22s rho %+.3f" % (bn, bs))

for label, mk in [("ridge regression", lambda: make_pipeline(StandardScaler(),
                                                             RidgeCV())),
                  ("gradient boosting", lambda: GradientBoostingRegressor(
                      random_state=0, n_estimators=200, max_depth=2))]:
    oof = np.full(len(ya), np.nan)
    for t in sorted(set(ga)):
        tr, te = ga != t, ga == t
        if tr.sum() < 40 or te.sum() < 5:
            continue
        m = mk()
        m.fit(Xa[tr], ya[tr])
        oof[te] = m.predict(Xa[te])
    ok = ~np.isnan(oof)
    r = spearmanr(oof[ok], ya[ok]).statistic
    print("%-22s held-out rho %+.3f   (n=%d)" % (label, r, ok.sum()))

print("\n   within-target rho (does it rank designs for ONE target?):")
for label, mk in [("ridge", lambda: make_pipeline(StandardScaler(), RidgeCV()))]:
    for t in sorted(set(ga)):
        tr, te = ga != t, ga == t
        if tr.sum() < 40 or te.sum() < 10:
            continue
        m = mk()
        m.fit(Xa[tr], ya[tr])
        r = spearmanr(m.predict(Xa[te]), ya[te]).statistic
        print("     %-16s n=%3d   rho %+.3f" % (t, te.sum(), r))
