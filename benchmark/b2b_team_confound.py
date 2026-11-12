"""Is the composition signal biology, or a proxy for which team submitted it?

29 teams, hit rates from under 1% to 48%. Any feature that happens to describe
one successful team's house style will look predictive when the teams are
pooled. The honest number is the within-team one.

Usage: python -I benchmark/b2b_team_confound.py <12k_all_metrics.csv> <12k_all_results.csv>
"""
import csv, sys, collections
import numpy as np
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
met = list(csv.DictReader(open(sys.argv[1], newline="", encoding="utf-8", errors="replace")))
res = {r["global_id"]: r for r in csv.DictReader(open(sys.argv[2], newline="", encoding="utf-8", errors="replace"))}

def num(v):
    try:
        if v in (None,"","NA","nan","None"): return np.nan
        return float(v)
    except Exception: return np.nan

def auc(x, y):
    ok = ~np.isnan(x); x, y = x[ok], y[ok]
    if y.sum() < 5 or (1-y).sum() < 5: return np.nan
    s = np.argsort(x); r = np.empty(len(x)); r[s] = np.arange(1, len(x)+1)
    xs = x[s]; i = 0
    while i < len(xs):
        j = i
        while j+1 < len(xs) and xs[j+1] == xs[i]: j += 1
        if j > i: r[s[i:j+1]] = (i+j+2)/2
        i = j+1
    n1, n0 = y.sum(), (1-y).sum()
    return (r[y==1].sum() - n1*(n1+1)/2)/(n1*n0)

rr = [res.get(m["global_id"], {}) for m in met]
four = np.array([r.get("fourclass","") for r in rr])
keep = (four=="Up")|(four=="Not Sig")|(four=="Down")
y = np.array([1.0 if f=="Up" else 0.0 for f in four])
team = np.array([m["team"] for m in met])

feats = ["ratio_KE","boltz_A_BC_max_ipSAE","chai1_A_BC_max_ipSAE","boltz_hs_iplddt",
         "boltz_rosetta_A_BC_dG_separated","amina_solubility","rnn_solubility",
         "esm2_binder_LL_car","boltz_sap_score","boltz_plddt_binder_none",
         "chai1_iptm_B","boltz_A_B_ipSAE","msa_depth","closest_ab_pident"]

print("Pooled AUC for predicting T-cell function (Up), among %d recovered designs\n" % keep.sum())
for f in feats:
    x = np.array([num(m[f]) for m in met])
    print("   %-34s %.3f" % (f, auc(x[keep], y[keep])))

print("\nHow concentrated are the 707 hits across teams?")
hits = collections.Counter(team[(four=="Up")])
tot  = collections.Counter(team[keep])
rows = sorted(tot, key=lambda t: -hits[t]/max(tot[t],1))
print("   %-26s %6s %7s %7s" % ("team","hits","tested","rate"))
for t in rows[:8]:
    print("   %-26s %6d %7d %6.1f%%" % (t[:26], hits[t], tot[t], 100*hits[t]/tot[t]))
print("   ... %d teams total" % len(tot))

print("\nLeave-one-team-out: does ratio_KE still rank within a single team?")
for f in ["ratio_KE","boltz_A_BC_max_ipSAE"]:
    x = np.array([num(m[f]) for m in met])
    per = []
    for t in tot:
        s = keep & (team==t)
        a = auc(x[s], y[s])
        if not np.isnan(a): per.append(a)
    per = np.array(per)
    print("   %-24s mean within-team AUC %.3f over %d teams (range %.2f-%.2f)"
          % (f, per.mean(), len(per), per.min(), per.max()))
