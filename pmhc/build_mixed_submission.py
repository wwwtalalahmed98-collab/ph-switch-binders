"""Assemble the final 20 for Challenge 5: ten scaffold-docked, ten BoltzGen.

WHY A MIX, AND WHY THIS IS NOT A RANKING. As of 8 Nov the three-state AlphaFold
screen is known to be blind: a T-cell receptor that provably binds this exact pMHC
scored 0.158 against a shuffled sequence's 0.134, a gap of +0.024 on a measurement
whose run-to-run noise is 0.022. So nothing in this file predicts binding, and the
order of the twenty carries no claim about affinity. What it does carry:

  - expression and fold plausibility (pmhc/screen.py, used since challenge 1)
  - diversity, so twenty designs are twenty bets rather than one
  - de novo standing, which is the one criterion still measurable today

THE TWO HALVES DIFFER IN EXACTLY ONE WAY THAT MATTERS.

  scaffold-docked  Backbones from the de novo scaffold set published with Liu et al.
                   (Science 389:386-391, 2025), built for OTHER pMHC targets, docked
                   here and fully resurfaced by ProteinMPNN. Sequence novelty is
                   clean - the official prescreen returns prior-art similarity 0.0 -
                   but the FOLD is reused, so ProteinTyper may cap these at novelty
                   level 2, below the level-3 floor. That is unmeasured, not
                   disproven, and it is the reason for not submitting twenty of them.

  boltzgen         Backbones generated from nothing by BoltzGen against this target,
                   so novel in fold as well as sequence. The same run produced 59 of
                   120 designs that were just copies of a target chain and 18 that
                   were antibody or TCR frameworks; those are excluded, and what is
                   left carries prior-art, design and target similarity all 0.0 from
                   the official prescreen.

Submitting ten of each hedges the fold-novelty risk, matches the organisers' stated
preference for diversity across method families, and costs nothing, since a single
submission is capped at ten designs regardless.

TIE-BREAK, AND WHAT IT IS WORTH. Within the scaffold-docked half, designs are ordered
by delta_vs_shuffled from the ProteinMPNN specificity screen: how much more the pocket
prefers Asp at peptide position 6 than a shuffled control does. It is a sequence
model's opinion and it has never been checked against a control the way the AlphaFold
screen just was. It breaks ties; it does not claim anything.

Run:  python pmhc/build_mixed_submission.py
"""
import csv
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "pmhc"))
import screen  # noqa: E402
from prep_boltzgen import best_offset_identity  # noqa: E402
from qc_adaptyv import liabilities  # noqa: E402

SUB = os.path.join(ROOT, "submission")
RESULTS = os.path.join(ROOT, "pmhc", "results")

N_EACH = 10
MAX_IDENTITY = 0.60      # between any two designs in the final twenty
MAX_PER_PARENT = 3       # scaffold-docked half only
MOLECULE_CLASS = "miniprotein"   # what ProteinTyper called both sets


def load(path, scol_names, ncol_names):
    if not os.path.exists(path):
        raise SystemExit("missing %s" % path)
    rows = list(csv.DictReader(open(path)))
    if not rows:
        raise SystemExit("%s is empty" % path)
    cols = list(rows[0].keys())
    scol = next((c for c in cols if c.lower() in scol_names), None)
    ncol = next((c for c in cols if c.lower() in ncol_names), None)
    if scol is None:
        raise SystemExit("no sequence column in %s (%s)" % (path, cols))
    out = []
    for i, r in enumerate(rows, 1):
        s = (r[scol] or "").strip().upper()
        if not s:
            continue
        r["seq"] = s
        r["name"] = (r[ncol].strip() if ncol and r.get(ncol) else "row_%03d" % i)
        out.append(r)
    return out


def pick(pool, n, per_parent=None, chosen=()):
    """Take n, refusing anything too close to something already taken."""
    taken, dropped_sim, dropped_parent = [], 0, 0
    counts = {}
    for d in pool:
        if per_parent is not None:
            p = d.get("parent", "")
            if counts.get(p, 0) >= per_parent:
                dropped_parent += 1
                continue
        if any(best_offset_identity(d["seq"], o["seq"]) > MAX_IDENTITY
               for o in list(chosen) + taken):
            dropped_sim += 1
            continue
        taken.append(d)
        if per_parent is not None:
            counts[d.get("parent", "")] = counts.get(d.get("parent", ""), 0) + 1
        if len(taken) == n:
            break
    return taken, dropped_sim, dropped_parent


def report(tag, pool, kept):
    tiers = {}
    for d in pool:
        tiers[d["tier_calc"]] = tiers.get(d["tier_calc"], 0) + 1
    print()
    print("--- %s ---" % tag)
    print("    pool %d; tiers %s" % (len(pool), tiers))
    print("    eligible (good or marginal): %d" % len(kept))


def main():
    print("=" * 72)
    print("FINAL TWENTY - ten scaffold-docked, ten BoltzGen")
    print("=" * 72)
    print("No design here is backed by a working prediction of binding. The AF2")
    print("screen was shown to be blind on 8 Nov (a known binder of this target")
    print("scored +0.024 over a shuffle, against 0.022 of noise), so this file")
    print("is a buildability, diversity and novelty selection and nothing more.")

    scaf = load(os.path.join(SUB, "challenge5_candidates.csv"),
                ("seq", "sequence"), ("name", "id"))
    bolt = load(os.path.join(RESULTS, "boltzgen_selected.csv"),
                ("seq", "sequence"), ("name", "id"))

    for d in scaf + bolt:
        d.update(screen.report(d["seq"]))
        d["tier_calc"] = d["tier"]

    # A control is a known binder of this target. It must never reach a submission.
    bad = [d for d in scaf + bolt if d["name"].startswith("CTRL_")]
    if bad:
        raise SystemExit("refusing to build from a set containing controls: %s"
                         % [d["name"] for d in bad])

    # The organisers' own expression filters, which this project had never run until
    # 8 Nov. The dominant one is the instability index: 12 of the first 20 we built
    # failed it. Anything at redesign level is now excluded outright, and 'clean'
    # beats 'flag' in the ordering.
    for d in scaf + bolt:
        d["qc"] = liabilities(d["seq"])
    def qc_ok(d):
        return d["qc"]["verdict"] != "redesign"
    for tag, pool in (("scaffold-docked", scaf), ("boltzgen", bolt)):
        t = {}
        for d in pool:
            t[d["qc"]["verdict"]] = t.get(d["qc"]["verdict"], 0) + 1
        print()
        print("    %s expression QC: %s" % (tag, t))

    scaf_ok = [d for d in scaf
               if d["tier_calc"] in ("good", "marginal") and qc_ok(d)]
    bolt_ok = [d for d in bolt
               if d["tier_calc"] in ("good", "marginal") and qc_ok(d)]
    report("scaffold-docked", scaf, scaf_ok)
    report("boltzgen", bolt, bolt_ok)

    def fnum(d, k):
        try:
            return float(d.get(k) or 0.0)
        except ValueError:
            return 0.0

    # Scaffold half: buildability first, then the unvalidated specificity tie-break.
    scaf_ok.sort(key=lambda d: (0 if d["qc"]["verdict"] == "clean" else 1,
                                0 if d["tier_calc"] == "good" else 1,
                                -fnum(d, "delta_vs_shuffled"),
                                -fnum(d, "pep_fraction")))
    bolt_ok.sort(key=lambda d: (0 if d["qc"]["verdict"] == "clean" else 1,
                                0 if d["tier_calc"] == "good" else 1,
                                d["qc"]["instability"]))
    # BoltzGen half: the file already arrives ordered by tier then distance.

    a, sa, pa = pick(scaf_ok, N_EACH, per_parent=MAX_PER_PARENT)
    b, sb, _ = pick(bolt_ok, N_EACH, chosen=a)
    print()
    print("    scaffold-docked taken %d  (dropped %d too similar, %d over the"
          " per-parent cap of %d)" % (len(a), sa, pa, MAX_PER_PARENT))
    print("    boltzgen        taken %d  (dropped %d too similar to something"
          " already taken)" % (len(b), sb))
    if len(a) < N_EACH or len(b) < N_EACH:
        print()
        print("    WARNING: one half is short. The twenty below is not twenty.")

    final = []
    for i, d in enumerate(a + b, 1):
        final.append({"name": "krasG12D_%02d" % i,
                      "sequence": d["seq"],
                      "molecule_class": MOLECULE_CLASS,
                      "source": ("scaffold-docked" if d in a else "boltzgen"),
                      "internal_name": d["name"],
                      "parent": d.get("parent", ""),
                      "length": d["length"],
                      "tier": d["tier_calc"],
                      "kd_hydropathy": d["kd"],
                      "charged_frac": d["charged_frac"],
                      "instability": d["qc"]["instability"],
                      "gravy": d["qc"]["gravy"],
                      "qc_verdict": d["qc"]["verdict"],
                      "parent_identity": d.get("parent_identity", "")})

    seqs = [d["sequence"] for d in final]
    if len(set(seqs)) != len(seqs):
        raise SystemExit("duplicate sequences in the final set")
    for d in final:
        if not (10 <= len(d["sequence"]) <= 250):
            raise SystemExit("%s is %d aa, outside 10-250"
                             % (d["name"], len(d["sequence"])))
    worst = max(best_offset_identity(x, y)
                for i, x in enumerate(seqs) for y in seqs[i + 1:])
    print()
    print("    final set: %d designs, lengths %d-%d, highest identity between any"
          " two %.2f" % (len(final), min(len(s) for s in seqs),
                         max(len(s) for s in seqs), worst))

    print()
    hdr = "%-14s %-16s %4s %-9s %7s %7s %-6s %s"
    print(hdr % ("submit name", "source", "len", "tier", "instab", "GRAVY",
                 "qc", "internal"))
    for d in final:
        print(hdr % (d["name"], d["source"], d["length"], d["tier"],
                     "%.1f" % d["instability"], "%+.3f" % d["gravy"],
                     d["qc_verdict"], d["internal_name"][:22]))

    out = os.path.join(SUB, "challenge5_final20.csv")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "sequence", "molecule_class"])
        for d in final:
            w.writerow([d["name"], d["sequence"], d["molecule_class"]])
    det = os.path.join(SUB, "challenge5_final20_detail.tsv")
    with open(det, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        keys = ["name", "source", "internal_name", "parent", "length", "tier",
                "kd_hydropathy", "charged_frac", "instability", "gravy",
                "qc_verdict", "parent_identity"]
        w.writerow(keys)
        for d in final:
            w.writerow([d[k] for k in keys])
    json.dump(final, open(os.path.join(RESULTS, "final20.json"), "w"), indent=1)
    print()
    print("wrote %s" % out)
    print("wrote %s" % det)
    print()
    print("The order of these twenty is not a confidence ranking. Nothing in this")
    print("project currently predicts binding for either half.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
