# Build the 20-design Track-3 set from the mouse-aware switch run (tnf_switch5_repair_fix).
#
# The switch5 run offers up to four histidine placements per parent and keeps only those
# that hold BOTH species. This picks one placement per parent and ten parents, then
# verifies every switch against its own control before writing anything.
#
# Selection rule, applied in order (objectives are ranked switch > mouse > affinity):
#   1. prefer two histidines over one - one His leaves 50% of molecules uncharged at
#      pH 6.0 and still binding, two leaves 25%
#   2. prefer more histidines facing a CONSERVED anchor (present in human and mouse),
#      because a His facing a human-only cation cannot switch in mouse
#   3. prefer the higher worse-species i_ptm, min(human, mouse)
# Then take the ten highest-ranked parents.
import csv, os, sys

SRC = sys.argv[1]
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "submission")
TOTAL_SLOTS = 20
N_SWITCH = 12
MIN_MOUSE_ONSITE = 5
# A control's only job is to make a pH response attributable. A twin that carries a
# native histidine of its own cannot serve as a pH-inert baseline, so it does not earn
# a slot - the slot goes to another switch instead, which can at least win objective 1.
# Controls therefore go to the highest-ranked switches whose twin is histidine-free.

ctl, var = {}, {}
for row in open(SRC, encoding="utf-8"):
    f = row.rstrip("\n").split("\t")
    if not f or not f[0]:
        continue
    if f[0] == "C":
        ctl[f[1]] = dict(name=f[1], length=int(f[2]), repair_level=f[3], repair=f[4],
                         seq=f[5], h=float(f[6]), hp=float(f[7]), m=float(f[8]),
                         mp=float(f[9]), mo=int(f[10]), free_cys=int(f[11]),
                         parent_h=float(f[12]), parent_m=float(f[13]),
                         parent_m_known=bool(int(f[14])), rule_best=f[15])
    elif f[0] == "V":
        var.setdefault(f[1], []).append(dict(
            parent=f[1], kind=f[2], his=f[3], seq=f[4], h=float(f[5]), hp=float(f[6]),
            hr=float(f[7]), m=float(f[8]), mp=float(f[9]), mr=float(f[10]),
            mo=int(f[11]), n_his=int(f[12]), n_anchor=int(f[13]),
            switch_score=float(f[14]), free_cys=int(f[15]), n_runs=int(f[16]),
            partners=f[17]))

# ---- verification: a switch is only interpretable if it differs from its twin at the
# ---- histidine positions and nowhere else. Anything that fails here is dropped, loudly.
AA = set("ACDEFGHIKLMNPQRSTVWY")
problems = []

def verify(c, v):
    bad = []
    if len(v["seq"]) != len(c["seq"]):
        bad.append("length %d vs control %d" % (len(v["seq"]), len(c["seq"])))
        return bad
    if len(c["seq"]) != c["length"]:
        bad.append("control length %d != declared %d" % (len(c["seq"]), c["length"]))
    if set(v["seq"]) - AA:
        bad.append("non-standard residues %s" % sorted(set(v["seq"]) - AA))
    subs = v["his"].split(",")
    want = {}
    for s in subs:
        frm, pos, to = s[0], int(s[1:-1]), s[-1]
        if to != "H":
            bad.append("%s does not substitute to H" % s)
        if c["seq"][pos - 1] != frm:
            bad.append("%s: control has %s at %d" % (s, c["seq"][pos - 1], pos))
        if v["seq"][pos - 1] != "H":
            bad.append("%s: switch has %s at %d, not H" % (s, v["seq"][pos - 1], pos))
        want[pos - 1] = "H"
    diff = [i for i in range(len(v["seq"])) if v["seq"][i] != c["seq"][i]]
    if sorted(diff) != sorted(want):
        bad.append("differs at %s, expected only %s"
                   % ([i + 1 for i in diff], sorted(i + 1 for i in want)))
    if len(subs) != v["n_his"]:
        bad.append("n_his %d but %d substitutions named" % (v["n_his"], len(subs)))
    if v["free_cys"]:
        bad.append("%d free cysteine(s) remain" % v["free_cys"])
    if v["mo"] < MIN_MOUSE_ONSITE:
        bad.append("only %d mouse-site contacts" % v["mo"])
    return bad

clean = {}
for p, vs in var.items():
    c = ctl[p]
    for v in vs:
        bad = verify(c, v)
        if bad:
            problems.append((p, v["his"], bad))
        else:
            clean.setdefault(p, []).append(v)

# ---- pick one placement per parent, then the ten best parents
def rank_variant(v):
    return (v["n_his"], v["n_anchor"], min(v["h"], v["m"]))

chosen = {}
for p, vs in clean.items():
    chosen[p] = max(vs, key=rank_variant)

order = sorted(chosen.values(), key=rank_variant, reverse=True)
picked, bench = order[:N_SWITCH], order[N_SWITCH:]

n_controls = TOTAL_SLOTS - len(picked)
clean_twin = [v for v in picked if "H" not in ctl[v["parent"]]["seq"]]
with_controls = {v["parent"] for v in clean_twin[:n_controls]}
if len(clean_twin) < n_controls:
    raise SystemExit("only %d histidine-free twins available, %d control slots to fill - "
                     "decide explicitly rather than silently pairing a compromised twin"
                     % (len(clean_twin), n_controls))

# ---- write
os.makedirs(OUT, exist_ok=True)
des_rows, met_rows = [], []
for i, v in enumerate(picked, 1):
    c = ctl[v["parent"]]
    sw = "tnf_sw_%02d" % i
    ct = "tnf_ct_%02d" % i if v["parent"] in with_controls else ""
    des_rows.append((sw, v["seq"], "single_chain"))
    met_rows.append(dict(name=sw, role="switch", paired_with=ct, internal_id=v["parent"],
                         length=len(v["seq"]), n_his=v["n_his"], his_substitutions=v["his"],
                         his_partner=v["partners"], n_anchor_his=v["n_anchor"],
                         switch_score=v["switch_score"], repair_level=c["repair_level"],
                         repair=c["repair"], human_i_ptm=v["h"], human_plddt=v["hp"],
                         human_retained=v["hr"], mouse_i_ptm=v["m"], mouse_plddt=v["mp"],
                         mouse_retained=v["mr"], mouse_on_site=v["mo"],
                         parent_human_i_ptm=c["parent_h"], parent_mouse_i_ptm=c["parent_m"],
                         control_has_native_his=int("H" in c["seq"])))
for i, v in enumerate(picked, 1):
    if v["parent"] not in with_controls:
        continue
    c = ctl[v["parent"]]
    sw, ct = "tnf_sw_%02d" % i, "tnf_ct_%02d" % i
    des_rows.append((ct, c["seq"], "single_chain"))
    met_rows.append(dict(name=ct, role="control", paired_with=sw, internal_id=c["name"],
                         length=len(c["seq"]), n_his=0, his_substitutions="",
                         his_partner="", n_anchor_his=0, switch_score=0.0,
                         repair_level=c["repair_level"], repair=c["repair"],
                         human_i_ptm=c["h"], human_plddt=c["hp"], human_retained=1.0,
                         mouse_i_ptm=c["m"], mouse_plddt=c["mp"],
                         mouse_retained=round(c["m"] / c["parent_m"], 3) if c["parent_m"] else None,
                         mouse_on_site=c["mo"], parent_human_i_ptm=c["parent_h"],
                         parent_mouse_i_ptm=c["parent_m"],
                         control_has_native_his=int("H" in c["seq"])))

with open(os.path.join(OUT, "challenge2_proposed_v3.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["name", "sequence", "molecule_class"])
    w.writerows(des_rows)
with open(os.path.join(OUT, "challenge2_metrics_v3.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(met_rows[0].keys()))
    w.writeheader()
    w.writerows(met_rows)

# ---- report
print("VERIFICATION")
if problems:
    for p, his, bad in problems:
        print("  DROPPED %-14s %-12s %s" % (p, his, "; ".join(bad)))
else:
    print("  every variant differs from its control at the named histidines and nowhere else")
print("  %d parents with at least one verified placement; %d verified placements total"
      % (len(clean), sum(len(v) for v in clean.values())))
print("\nSELECTED (%d switch + %d control = %d designs)"
      % (len(picked), n_controls, len(des_rows)))
print("%-5s %-14s %-12s %3s %3s %6s %6s %6s %5s  %s"
      % ("slot", "lineage", "histidines", "nH", "anc", "human", "mouse", "min", "site", "twin"))
for i, v in enumerate(picked, 1):
    c = ctl[v["parent"]]
    twin = ("ct_%02d" % i if v["parent"] in with_controls
            else "none - twin has native H")
    print("sw_%02d %-14s %-12s %3d %3d %6.3f %6.3f %6.3f %5d  %s"
          % (i, v["parent"], v["his"], v["n_his"], v["n_anchor"], v["h"], v["m"],
             min(v["h"], v["m"]), v["mo"], twin))
print("\nNOT SELECTED (verified, but outside the switch slots)")
for v in bench:
    print("      %-14s %-12s %3d %3d %6.3f %6.3f %6.3f %5d"
          % (v["parent"], v["his"], v["n_his"], v["n_anchor"], v["h"], v["m"],
             min(v["h"], v["m"]), v["mo"]))
nh2 = sum(1 for v in picked if v["n_his"] == 2)
hisctl = sum(1 for v in des_rows
             if v[0].startswith("tnf_ct_") and "H" in v[1])
print("\n%d of %d switches carry two histidines; %d of %d controls carry a native histidine"
      % (nh2, len(picked), hisctl, n_controls))
print("lowest mouse i_ptm in the set: %.3f; lowest mouse-site contacts: %d"
      % (min(v["m"] for v in picked), min(v["mo"] for v in picked)))
