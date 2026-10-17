# Rebuild the 20-design Track-3 set on burial, the variable that actually predicts the
# charge effect (r = -0.58) rather than the Cbeta-to-cation distance we originally
# designed with (r = +0.24).
#
# Selection rule, applied in order and stated so it can be argued with:
#
#   1. CONSERVED ANCHOR FIRST. Keep only placements whose histidines face one of the two
#      cations present in both species (K166 chain A, R108 chain B). A histidine facing a
#      human-only cation cannot push off the mouse target, and mouse is ranked above
#      affinity. If a lineage has no anchor-facing placement, fall back to all of them.
#   2. BURIAL ON BINDING, worst histidine first. A histidine that stays solvent-exposed
#      when bound protonates at no cost and cannot contribute to release. This is
#      geometry from the predicted complex, not a model's opinion.
#   3. Then rank lineages by that burial and take the top twelve.
#
# This deliberately prefers single histidines over pairs. The Henderson-Hasselbalch
# argument for two (one His leaves 50% of molecules binding at pH 6.0, two leaves 25%)
# assumes BOTH are doing work. Every pair we have has a partner buried by <= 0.10 and
# one buried by 0.000, so the effective count is one and the pH 7.4 cost is two.
import csv, os, sys

SWITCH5, BURIAL = sys.argv[1], sys.argv[2]
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "submission")
TOTAL_SLOTS, N_SWITCH = 20, 12
MIN_MOUSE_ONSITE = 5

ctl, var = {}, {}
for row in open(SWITCH5, encoding="utf-8"):
    f = row.rstrip("\n").split("\t")
    if not f or not f[0]:
        continue
    if f[0] == "C":
        ctl[f[1]] = dict(name=f[1], length=int(f[2]), repair_level=f[3], repair=f[4],
                         seq=f[5], h=float(f[6]), hp=float(f[7]), m=float(f[8]),
                         mp=float(f[9]), mo=int(f[10]), parent_h=float(f[12]),
                         parent_m=float(f[13]))
    elif f[0] == "V":
        var[(f[1], f[3])] = dict(parent=f[1], kind=f[2], his=f[3], seq=f[4],
                                 h=float(f[5]), hp=float(f[6]), hr=float(f[7]),
                                 m=float(f[8]), mp=float(f[9]), mr=float(f[10]),
                                 mo=int(f[11]), n_his=int(f[12]), n_anchor=int(f[13]),
                                 switch_score=float(f[14]), partners=f[17])

n_join = 0
for row in open(BURIAL, encoding="utf-8"):
    f = row.rstrip("\n").split("\t")
    if len(f) < 7:
        continue
    key = (f[0], f[1])
    if key not in var:
        raise SystemExit("burial row %s has no matching switch5 placement - the two "
                         "files disagree, refusing to build from them" % (key,))
    var[key].update(burial=float(f[2]), burial_sum=float(f[3]), charge=float(f[4]),
                    burial_human=float(f[5]), burial_mouse=float(f[6]))
    n_join += 1
missing = [k for k, v in var.items() if "burial" not in v]
if missing:
    raise SystemExit("%d placements have no burial measurement: %s" % (len(missing), missing))

# ---- verification: a switch must differ from its twin ONLY at the histidines --------
AA = set("ACDEFGHIKLMNPQRSTVWY")
problems = []
def verify(c, v):
    bad = []
    if len(v["seq"]) != len(c["seq"]):
        return ["length %d vs control %d" % (len(v["seq"]), len(c["seq"]))]
    if set(v["seq"]) - AA:
        bad.append("non-standard residues")
    want = {}
    for s in v["his"].split(","):
        frm, pos = s[0], int(s[1:-1])
        if s[-1] != "H": bad.append("%s does not substitute to H" % s)
        if c["seq"][pos-1] != frm: bad.append("%s: control has %s at %d" % (s, c["seq"][pos-1], pos))
        if v["seq"][pos-1] != "H": bad.append("%s: switch has %s at %d" % (s, v["seq"][pos-1], pos))
        want[pos-1] = "H"
    diff = [i for i in range(len(v["seq"])) if v["seq"][i] != c["seq"][i]]
    if sorted(diff) != sorted(want):
        bad.append("differs at %s, expected only %s"
                   % ([i+1 for i in diff], sorted(i+1 for i in want)))
    if v["mo"] < MIN_MOUSE_ONSITE:
        bad.append("only %d mouse-site contacts" % v["mo"])
    return bad

clean = {}
for key, v in var.items():
    bad = verify(ctl[v["parent"]], v)
    if bad: problems.append((key, bad))
    else: clean.setdefault(v["parent"], []).append(v)

# ---- the rule ----------------------------------------------------------------------
chosen, notes = {}, {}
for p, vs in clean.items():
    anchored = [v for v in vs if v["n_anchor"] >= 1]
    pool, why = (anchored, "anchor-facing") if anchored else (vs, "NO anchor-facing option")
    best = max(pool, key=lambda v: (v["burial"], min(v["h"], v["m"])))
    chosen[p] = best
    notes[p] = why
    alt = max(vs, key=lambda v: v["burial"])
    if alt is not best:
        notes[p] += "; %s has more burial (%.3f) but faces a non-conserved cation" % (
            alt["his"], alt["burial"])

order = sorted(chosen.values(), key=lambda v: (v["burial"], min(v["h"], v["m"])),
               reverse=True)
picked, bench = order[:N_SWITCH], order[N_SWITCH:]
n_controls = TOTAL_SLOTS - len(picked)
clean_twin = [v for v in picked if "H" not in ctl[v["parent"]]["seq"]]
if len(clean_twin) < n_controls:
    raise SystemExit("only %d histidine-free twins among the picks, need %d"
                     % (len(clean_twin), n_controls))
with_ctl = {v["parent"] for v in clean_twin[:n_controls]}

# ---- write --------------------------------------------------------------------------
os.makedirs(OUT, exist_ok=True)
des, met = [], []
for i, v in enumerate(picked, 1):
    c = ctl[v["parent"]]
    sw = "tnf_sw_%02d" % i
    ct = "tnf_ct_%02d" % i if v["parent"] in with_ctl else ""
    des.append((sw, v["seq"], "single_chain"))
    met.append(dict(name=sw, role="switch", paired_with=ct, internal_id=v["parent"],
                    length=len(v["seq"]), n_his=v["n_his"], his_substitutions=v["his"],
                    his_partner=v["partners"], n_anchor_his=v["n_anchor"],
                    burial_on_binding=v["burial"], charge_lys_minus_gln=v["charge"],
                    repair_level=c["repair_level"], repair=c["repair"],
                    human_i_ptm=v["h"], human_plddt=v["hp"], human_retained=v["hr"],
                    mouse_i_ptm=v["m"], mouse_plddt=v["mp"], mouse_retained=v["mr"],
                    mouse_on_site=v["mo"], parent_human_i_ptm=c["parent_h"],
                    parent_mouse_i_ptm=c["parent_m"],
                    control_has_native_his=int("H" in c["seq"])))
for i, v in enumerate(picked, 1):
    if v["parent"] not in with_ctl: continue
    c = ctl[v["parent"]]
    sw, ct = "tnf_sw_%02d" % i, "tnf_ct_%02d" % i
    des.append((ct, c["seq"], "single_chain"))
    met.append(dict(name=ct, role="control", paired_with=sw, internal_id=c["name"],
                    length=len(c["seq"]), n_his=0, his_substitutions="", his_partner="",
                    n_anchor_his=0, burial_on_binding="", charge_lys_minus_gln="",
                    repair_level=c["repair_level"], repair=c["repair"],
                    human_i_ptm=c["h"], human_plddt=c["hp"], human_retained=1.0,
                    mouse_i_ptm=c["m"], mouse_plddt=c["mp"],
                    mouse_retained=round(c["m"]/c["parent_m"], 3) if c["parent_m"] else None,
                    mouse_on_site=c["mo"], parent_human_i_ptm=c["parent_h"],
                    parent_mouse_i_ptm=c["parent_m"],
                    control_has_native_his=int("H" in c["seq"])))

with open(os.path.join(OUT, "challenge2_proposed_v4.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh); w.writerow(["name", "sequence", "molecule_class"]); w.writerows(des)
with open(os.path.join(OUT, "challenge2_metrics_v4.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(met[0].keys())); w.writeheader(); w.writerows(met)

# ---- report -------------------------------------------------------------------------
print("joined %d burial measurements onto %d verified placements" % (n_join, len(var)))
if problems:
    for k, bad in problems: print("  DROPPED %s %s" % (k, "; ".join(bad)))
else:
    print("every placement differs from its control only at the named histidines")
print("\nSELECTED (%d switch + %d control)" % (len(picked), n_controls))
print("%-5s %-14s %-12s %3s %3s %7s %8s %7s %7s %5s %s"
      % ("slot", "lineage", "histidines", "nH", "anc", "burial", "charge", "human",
         "mouse", "site", "twin"))
for i, v in enumerate(picked, 1):
    print("sw_%02d %-14s %-12s %3d %3d %7.3f %8.3f %7.3f %7.3f %5d %s"
          % (i, v["parent"], v["his"], v["n_his"], v["n_anchor"], v["burial"],
             v["charge"], v["h"], v["m"], v["mo"],
             "ct_%02d" % i if v["parent"] in with_ctl else "-"))
print("\nNOT SELECTED")
for v in bench:
    print("      %-14s %-12s %7.3f  %s" % (v["parent"], v["his"], v["burial"],
                                           "burial below the cut"))
print("\nrule notes per lineage:")
for p in sorted(notes):
    if "NO anchor" in notes[p] or "more burial" in notes[p]:
        print("  %-14s %s" % (p, notes[p]))
dead = [v for v in picked if v["burial"] < 0.02]
print("\n%d of %d chosen placements have a histidine buried by <0.02 (was 4 of 12)"
      % (len(dead), len(picked)))
print("%d of %d switches carry two histidines" % (sum(1 for v in picked if v["n_his"] == 2), len(picked)))
print("lowest mouse i_ptm %.3f, lowest mouse-site contacts %d, lowest burial %.3f"
      % (min(v["m"] for v in picked), min(v["mo"] for v in picked),
         min(v["burial"] for v in picked)))
