"""Choose which BoltzGen designs are worth an AlphaFold screen, and in what order.

WHY THIS STEP EXISTS. BoltzGen returns a few hundred designs in no particular order,
and the three-state AF2 screen can afford about 24 of them (roughly 5 minutes each on
a T4, times three states, times two models). Feeding it the first 24 rows would screen
an arbitrary subset. This picks them.

WHAT IT RANKS ON, and what it does not. We ran BoltzGen's generation steps only: its
own quality filter lives behind ESMFold2 scoring, which we skipped. So there is no
model-provided quality score here, and nothing in this file predicts binding. All it
does is the buildability screen this project has used since challenge 1 - can this be
expressed, is it plausibly a fold - plus a diversity cap so the 24 are not 24 copies
of one idea. Binding is what the AF2 screen afterwards is for.

THE DIVERSITY CAP IS THE POINT. Twenty-four near-identical sequences is one bet, not
twenty-four. A design is only taken if it is below MAX_IDENTITY to everything already
taken. Because BoltzGen designs vary in length, identity is measured on the shorter of
the two over its best ungapped offset, which is cruder than an alignment and is
deliberately conservative: it over-estimates similarity, so it rejects more than a
proper alignment would.

Run:  python pmhc/prep_boltzgen.py [boltzgen_designs.csv]
"""
import csv
import json
import os
import random
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "pmhc"))
import screen  # noqa: E402

RESULTS = os.path.join(ROOT, "pmhc", "results")

N_PICK = 24
MAX_IDENTITY = 0.60
TIER_ORDER = {"good": 0, "marginal": 1, "fail": 2, "fold-breaker": 3}


def best_offset_identity(a, b):
    """Highest ungapped positional identity over any offset, scored on the shorter one.

    The shorter sequence is slid across the longer one through EVERY register,
    including ones where it hangs off either end; positions with no partner simply
    count as mismatches. The earlier version only slid from offset 0 rightwards, so a
    design carrying one extra residue at the front - a bare leading M, say - was
    thrown a register out and scored 0.20 against a chain it was 90% identical to.
    Nine of the BoltzGen designs were mis-scored exactly that way.
    """
    if not a or not b:
        return 0.0
    if len(a) > len(b):
        a, b = b, a
    n, m, best = len(a), len(b), 0
    for off in range(-(n - 1), m):
        same = 0
        for i, ch in enumerate(a):
            j = off + i
            if 0 <= j < m and ch == b[j]:
                same += 1
        best = max(best, same)
    return round(best / float(n), 4)


# Immunoglobulin / T-cell-receptor variable-domain framework motifs. A TCR is the
# natural binder of a pMHC and the official prescreen for this target holds 144 of
# them, so a design that is one is a de novo rule risk, whatever it scores.
IG_MOTIFS = ("WYRQ", "WYQQ", "WYKQ", "WYKY", "YRQK", "GAGTKLTV", "GQGTKLTV",
             "GGGTKLTV", "GPGTKLTV", "GAGTKVTV", "GGGTQVTV", "YYCA", "YYCV",
             "YYCQ", "TLTCT", "TLSCR", "TLTCR", "VTLTC", "TITCT")


def ig_hits(seq):
    return sum(1 for m in IG_MOTIFS if m in seq)


def target_chains():
    """Sequences of the chains BoltzGen was given, read from the structure itself."""
    t3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
          "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
          "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
          "TYR": "Y", "VAL": "V"}
    path = os.path.join(ROOT, "data", "pmhc", "states", "target_state_g12d.pdb")
    if not os.path.exists(path):
        raise SystemExit(
            "no target structure at %s - run build_states.py and then "
            "build_counter_states.py first." % path)
    seqs, seen = {}, set()
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")):
            key = (line[21], line[22:26])
            if key in seen:
                continue
            seen.add(key)
            seqs.setdefault(line[21], []).append(
                t3.get(line[17:20].strip().upper(), "X"))
    return {c: "".join(v) for c, v in seqs.items()}


def main(argv):
    path = (argv[1] if len(argv) > 1
            else os.path.join(RESULTS, "boltzgen_designs.csv"))
    if not os.path.exists(path):
        raise SystemExit(
            "no BoltzGen designs at %s\n"
            "Download boltzgen_designs.csv from the Kaggle run output into\n"
            "pmhc/results/ and run this again." % path)
    rows = list(csv.DictReader(open(path)))
    if not rows:
        raise SystemExit("%s is empty" % path)
    scol = next((c for c in rows[0] if c.lower() in ("sequence", "seq")), None)
    ncol = next((c for c in rows[0] if c.lower() in ("name", "id")), None)
    if scol is None:
        raise SystemExit("no sequence column in %s (%s)" % (path, list(rows[0])))

    seen, designs = set(), []
    for i, r in enumerate(rows, 1):
        s = (r[scol] or "").strip().upper()
        if not s or s in seen:
            continue
        seen.add(s)
        d = screen.report(s)
        d["seq"] = s
        d["name"] = (r[ncol].strip() if ncol and r.get(ncol) else "bg_%03d" % i)
        designs.append(d)

    print("unique designs: %d of %d rows" % (len(designs), len(rows)))

    # --- Is it a binder, or is it one of the things we handed the model? ---
    # Run 4 returned 172 copies of beta-2 microglobulin. Run 5 put B2M into the
    # target to close that hole; this measures how much of run 5 is still an echo,
    # and separately how much of it is a T-cell receptor, which is the OTHER obvious
    # answer for a pMHC and a de novo rule risk.
    tgt = target_chains()
    print()
    print("--- what the model was given, and what it gave back ---")
    for c in sorted(tgt):
        print("    target chain %s  %3d aa" % (c, len(tgt[c])))
    # A flat identity threshold is wrong here. A 60-residue design slid across a
    # 275-residue chain gets ~216 chances to line up, so it scores 0.3-0.4 on chance
    # alone, while a 100-residue one scores lower by chance. The null is measured per
    # design instead: shuffle it, which keeps its length and composition and destroys
    # its order, and see what identity that gets. Real echoes sit far above it.
    rnd = random.Random(20261108)
    for d in designs:
        d["id_target"] = max(best_offset_identity(d["seq"], t) for t in tgt.values())
        nulls = []
        for _ in range(5):
            sh = list(d["seq"])
            rnd.shuffle(sh)
            sh = "".join(sh)
            nulls.append(max(best_offset_identity(sh, t) for t in tgt.values()))
        d["id_null"] = round(sum(nulls) / len(nulls), 4)
        d["id_excess"] = round(d["id_target"] - d["id_null"], 4)
        d["ig"] = ig_hits(d["seq"])
    echo = [d for d in designs if d["id_excess"] > 0.12 or d["id_target"] > 0.6]
    ex = sorted(d["id_excess"] for d in designs)
    print("    chance identity to a target, per design   : median %.2f"
          % sorted(d["id_null"] for d in designs)[len(designs) // 2])
    print("    excess over that chance level             : median %+.2f  max %+.2f"
          % (ex[len(ex) // 2], ex[-1]))
    igs = [d for d in designs if d["ig"] >= 2 and d not in echo]
    print("    more than 60%% identical to a target chain : %d of %d"
          % (len(echo), len(designs)))
    print("    immunoglobulin / TCR framework motifs     : %d of %d"
          % (len(igs), len(designs)))
    keep = {id(d) for d in echo} | {id(d) for d in igs}
    designs = [d for d in designs if id(d) not in keep]
    print("    left after removing both                  : %d" % len(designs))
    if not designs:
        print()
        print("NOTHING SURVIVES. Every design is either a copy of the target or an")
        print("antibody/TCR framework. Report that, do not lower the bar.")
        return 1
    tiers = {}
    for d in designs:
        tiers[d["tier"]] = tiers.get(d["tier"], 0) + 1
    print()
    print("--- buildability screen (expression and fold plausibility only) ---")
    for t in ("good", "marginal", "fail", "fold-breaker"):
        n = tiers.get(t, 0)
        print("    %-14s %3d  (%4.1f%%)" % (t, n, 100.0 * n / len(designs)))
    lens = sorted(d["length"] for d in designs)
    kds = sorted(d["kd"] for d in designs)
    print("    length         min %d  median %d  max %d"
          % (lens[0], lens[len(lens) // 2], lens[-1]))
    print("    Kyte-Doolittle min %+.2f  median %+.2f  max %+.2f"
          % (kds[0], kds[len(kds) // 2], kds[-1]))
    ncys = sum(1 for d in designs if "C" in d["seq"])
    print("    contain cysteine %d of %d" % (ncys, len(designs)))

    pool = [d for d in designs if d["tier"] in ("good", "marginal")]
    print()
    print("eligible (good or marginal): %d" % len(pool))
    if not pool:
        print()
        print("NOTHING PASSES THE BUILDABILITY SCREEN. That is a result about the")
        print("generated set, not a bug: report it rather than lowering the bar.")
        return 1
    pool.sort(key=lambda d: (TIER_ORDER[d["tier"]], d["distance"]))

    picked = []
    rejected_similar = 0
    for d in pool:
        if any(best_offset_identity(d["seq"], p["seq"]) > MAX_IDENTITY
               for p in picked):
            rejected_similar += 1
            continue
        picked.append(d)
        if len(picked) == N_PICK:
            break

    print("dropped for being >%.0f%% identical to one already picked: %d"
          % (100 * MAX_IDENTITY, rejected_similar))
    print("selected: %d" % len(picked))
    if len(picked) < N_PICK:
        print("NOTE: fewer than %d distinct enough designs were available." % N_PICK)

    if len(picked) > 1:
        pair = max(best_offset_identity(a["seq"], b["seq"])
                   for i, a in enumerate(picked) for b in picked[i + 1:])
        print("highest identity between any two selected: %.2f" % pair)

    print()
    print("%-22s %4s %-9s %6s %6s %6s %6s"
          % ("name", "len", "tier", "kd", "chg", "dist", "excess"))
    for d in picked:
        print("%-22s %4d %-9s %+6.2f %6.2f %6.2f %6.2f"
              % (d["name"][:22], d["length"], d["tier"], d["kd"],
                 d["charged_frac"], d["distance"], d["id_excess"]))

    out_csv = os.path.join(RESULTS, "boltzgen_selected.csv")
    with open(out_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "seq"])
        for d in picked:
            w.writerow([d["name"], d["seq"]])
    json.dump(picked, open(os.path.join(RESULTS, "boltzgen_selected.json"), "w"),
              indent=1)
    print()
    print("wrote %s" % out_csv)
    print("Upload that as a Kaggle dataset and run the three-state AF2 screen on it.")
    print()
    print("Nothing above predicts binding. It is a buildability screen plus a")
    print("diversity cap, used only to choose which designs are worth GPU time.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
