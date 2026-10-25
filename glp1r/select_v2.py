"""Select the 20 GLP1R designs, across THREE populations.

READ THIS FIRST - HOW MUCH TO TRUST THIS FILE
---------------------------------------------
build_submission.py was written before any GLP1R design existed, which is why its
thresholds are worth something. THIS file is not that. The two-domain population did not
exist when that rule was written, and these thresholds were set on 25 Oct AFTER seeing
all 14 designs from the sealed-target run. A rule written after seeing the data is a rule
that can be tuned to pick what you already liked, and this one should be read with that
discount applied. It is recorded as an amendment rather than folded in silently.

What is NOT tuned after the fact: the measurements themselves (contact counts, depth,
polarity) are geometric facts computed by the same code that measured every earlier
design, and the decision to drop the threader class rests on a controlled experiment
whose terms were fixed in advance.

THE THREE POPULATIONS
---------------------
1. LANDING-PAD BINDERS - 12 slots. Bind the extracellular groove. Mechanistically these
   are ANTAGONISTS; exendin(9-39) is this molecule class and is a real drug. Rule
   unchanged from build_submission.py: i_ptm >= 0.50, on_site >= 8 of 22, hydropathy
   <= 0.0, rank by i_ptm then on_site.

2. TWO-DOMAIN DESIGNS - 5 slots. Grip the groove AND reach into the channel, N-terminus
   first, which is how the natural hormone works and the only shape that can activate.
   pass: ecd_on >= 4, tm_on >= 5, threaded (N-terminus is the deep end)
   rank: tm_on, then ecd_on

3. DUAL-SITE REVERSED - 3 slots. Touch both contact sets extensively but go in
   C-terminus first. By mechanism these are blockers, not agonists. They are included as
   binders, not as agonists, and the methods file says so.
   pass: ecd_on >= 5, tm_on >= 15, NOT threaded
   rank: tm_on + ecd_on

UNFILLED SLOTS ROLL OVER TO POPULATION 1. The two-domain pool came up one short of
its 5 slots. The choice was: leave a slot empty, or pad the two-domain pool with a design
that failed its own filter. Neither helps. Instead the spare goes to the landing-pad pool,
which has 25 designs passing a rule written before any design existed. This costs us one
agonist-shaped slot and buys one more design screened by the only rule here with
pre-registered thresholds. It is stated as a rule rather than done by editing N_ECD.

THREADERS GET ZERO SLOTS, and that is the one decision here with pre-registered
evidence behind it. The sealed-target run was set up as a controlled experiment: same
job, same method, one variable changed (the hole in the target closed). Threaders failed
again - 4 of 7 touched ZERO channel contacts. The earlier failure was the method, not
just our broken crop. We do not spend slots on a class we have now failed with twice.

NO CONFIDENCE THRESHOLD ON POPULATIONS 2 AND 3. This is deliberate and it is the single
biggest departure from how everything else in this project was selected. Across the 14
sealed-target designs, predicted interface confidence correlates -0.651 with touching the
hormone's channel contacts, and -0.804 among the two-domain designs alone. The previous
run showed -0.963 on 11 designs against a different target. Two targets, two design
tasks, same inverted sign. The highest-confidence design in the set touches ONE of 30
channel residues; the lowest touches 19. An i_ptm floor would delete the designs we want
and keep the ones we do not. Direct geometric measurement of the thing we actually care
about beats a proxy that has now pointed the wrong way twice.

BUILDABILITY GATE, added 25 Oct, AFTER the designs existed. Adaptyv expresses these in
cells, and two sequence features reliably stop that working: many free cysteines (they
pair with each other and the chain mis-folds) and long single-residue runs (they
aggregate). This was not screened in this challenge or the previous two, and 13 of the
first 20 selected designs tripped it - one groove design carried ELEVEN consecutive
leucines, another EIGHT cysteines.

The gate is deliberately loose: at most 6 cysteines and no run longer than 6. That removes
the designs that are obviously unbuildable while still filling all 20 slots from the 25
groove designs that passed the pre-registered rule. A strict version (3 Cys, run 4) would
leave only 4 groove designs and an 11-design submission.

Being honest about what this is: a filter invented after seeing the data, which is the same
weakness as the two-domain and reversed thresholds. It is more defensible than those only
because it screens whether a molecule can be MADE, not whether it does what we are trying
to predict. It is applied to every population, including the agonist candidates - none of
which it removes.

WHAT IS NOT SCREENED. Activation. Nothing available to us tests whether any design
switches the receptor on. We can build the shape that permits it; we cannot test it.
Foldability is also NOT filtered: the best-placed design (both_007) has helicity 0.14 and
pLDDT 0.458 and may simply not fold. That risk is reported per design rather than used to
cut, because cutting on it would remove the design that best matches the native geometry.
"""
import csv
import difflib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.abspath(os.path.join(HERE, "..", "submission"))
N_ECD, N_TWO, N_REV = 12, 5, 3
MAX_IDENT = 0.80

MAX_CYS, MAX_RUN = 6, 6
ECD = dict(i_ptm=0.50, on_site=8, kd=0.0)
TWO = dict(ecd=4, tm=5, kd=0.5)
REV = dict(ecd=5, tm=15)


def buildable(seq):
    """None if fine, else why not. See the BUILDABILITY GATE note in the header."""
    cys = seq.count("C")
    run, best = 1, 1
    for i in range(1, len(seq)):
        run = run + 1 if seq[i] == seq[i - 1] else 1
        best = max(best, run)
    run = best
    why = []
    if cys > MAX_CYS:
        why.append("%d cysteines > %d" % (cys, MAX_CYS))
    if run > MAX_RUN:
        why.append("run of %d identical > %d" % (run, MAX_RUN))
    return "; ".join(why) or None


def ident(a, b):
    return difflib.SequenceMatcher(None, a, b).ratio()


def n_ecd(r):
    """on-site groove contacts, whichever schema the record came from."""
    if "ecd_on" in r:
        return int(r["ecd_on"])
    return len(r.get("on_site") or [])


def n_tm(r):
    if "tm_on" in r:
        return int(r["tm_on"])
    return int(r.get("tm_contacts") or 0)


def dedupe(cands):
    kept, dropped = [], []
    for c in cands:
        twin = next((k for k in kept
                     if ident(k["sequence"], c["sequence"]) >= MAX_IDENT), None)
        (dropped.append((c, twin["name"])) if twin else kept.append(c))
    return kept, dropped


def screen(recs, pop):
    passed, failed = [], []
    for r in recs:
        why = []
        bad_build = buildable(r["sequence"])
        if bad_build:
            why.append(bad_build)
        kd = r.get("hydropathy")
        ip = r.get("true_iptm")
        if pop == "ecd":
            if ip is None or ip < ECD["i_ptm"]:
                why.append("i_ptm %.3f < %.2f" % (ip or 0, ECD["i_ptm"]))
            if n_ecd(r) < ECD["on_site"]:
                why.append("groove %d < %d" % (n_ecd(r), ECD["on_site"]))
            if kd is None or kd > ECD["kd"]:
                why.append("hydropathy %+.2f > %+.2f" % (kd or 0, ECD["kd"]))
        elif pop == "two":
            if not r.get("threaded"):
                why.append("N-terminus is not the deep end")
            if n_ecd(r) < TWO["ecd"]:
                why.append("groove %d < %d" % (n_ecd(r), TWO["ecd"]))
            if n_tm(r) < TWO["tm"]:
                why.append("channel %d < %d" % (n_tm(r), TWO["tm"]))
            if kd is None or kd > TWO["kd"]:
                why.append("hydropathy %+.2f > %+.2f" % (kd or 0, TWO["kd"]))
        else:
            if r.get("threaded"):
                why.append("N-terminus IS the deep end (belongs in the two-domain pool)")
            if n_ecd(r) < REV["ecd"]:
                why.append("groove %d < %d" % (n_ecd(r), REV["ecd"]))
            if n_tm(r) < REV["tm"]:
                why.append("channel %d < %d" % (n_tm(r), REV["tm"]))
        (failed if why else passed).append((r, why))
    return [p[0] for p in passed], failed


POPS = (
    ("ecd", "landing-pad binders", N_ECD,
     lambda r: (r.get("true_iptm") or 0, n_ecd(r)), "glp1r_ecd_%02d"),
    ("two", "two-domain designs", N_TWO,
     lambda r: (n_tm(r), n_ecd(r)), "glp1r_dual_%02d"),
    ("rev", "dual-site, reversed", N_REV,
     lambda r: (n_tm(r) + n_ecd(r),), "glp1r_rev_%02d"),
)


def main(paths):
    recs = []
    for p in paths:
        got = json.load(open(p))
        recs.extend(got)
        print("read %3d designs from %s" % (len(got), os.path.basename(p)))
    for r in recs:
        r.setdefault("on_site", [])

    groove_pool = [r for r in recs if r.get("mode") == "ecd"]
    channel_pool = [r for r in recs if r.get("mode") in ("both", "thread2", "thread")]
    print("\n%d groove-targeted, %d channel-targeted\n" % (len(groove_pool),
                                                           len(channel_pool)))

    # Populations 2 and 3 are screened first so any slots they cannot fill roll over to
    # population 1, whose rule is the pre-registered one. See the header.
    chosen = []
    shortfall = 0
    for pop, label, want, rank, tag in reversed(POPS):
        pool = groove_pool if pop == "ecd" else channel_pool
        if pop == "ecd":
            want += shortfall
            if shortfall:
                print("=== %d slot(s) rolled over from the channel populations ==="
                      % shortfall)
        ok, bad = screen(pool, pop)
        ok.sort(key=rank, reverse=True)
        kept, dropped = dedupe(ok)
        take = kept[:want]
        print("=== %s: %d of %d passed, %d near-duplicates dropped, taking %d of %d ==="
              % (label, len(ok), len(pool), len(dropped), len(take), want))
        if len(take) < want:
            print("  SHORT by %d. Not padding with designs that failed our own filters."
                  % (want - len(take)))
            if pop != "ecd":
                shortfall += want - len(take)
        for r in take:
            print("  %-18s len %3d  i_ptm %.3f  pLDDT %.3f  kd %+5.2f  hlx %.2f  "
                  "groove %2d  channel %2d  deep %5.1f"
                  % (r["name"], r["length"], r.get("true_iptm") or 0,
                     r.get("plddt") or 0, r.get("hydropathy") or 0,
                     r.get("helicity") or 0, n_ecd(r), n_tm(r), r.get("deepest") or 0))
        for r, tw in dropped:
            print("  dup  %-18s ~ %s" % (r["name"], tw))
        for r, why in bad[:6]:
            print("  out  %-18s %s" % (r["name"], "; ".join(why)))
        print()
        chosen.extend((pop, tag, r) for r in take)

    if not chosen:
        raise SystemExit("nothing passed. Do not submit; read the failures above.")

    os.makedirs(OUT, exist_ok=True)
    prop = os.path.join(OUT, "challenge3_proposed.csv")
    mets = os.path.join(OUT, "challenge3_metrics.csv")
    with open(prop, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "sequence", "molecule_class"])
        counts = {}
        for pop, tag, r in chosen:
            counts[pop] = counts.get(pop, 0) + 1
            r["sub_name"] = tag % counts[pop]
            # Four formats are allowed: single chain protein, nanobody, scFv, Fab.
            # Everything here is one chain of 30-70 aa, inside the stated 10-250 range.
            w.writerow([r["sub_name"], r["sequence"], "single_chain"])
    cols = ["sub_name", "name", "mode", "length", "true_iptm", "design_iptm", "plddt",
            "hydropathy", "helicity", "deepest", "nterm_depth", "threaded"]
    with open(mets, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(cols + ["groove_contacts", "channel_contacts", "population"])
        for pop, tag, r in chosen:
            w.writerow([r.get(c, "") for c in cols] + [n_ecd(r), n_tm(r), pop])

    print("wrote %d designs" % len(chosen))
    print("  %s" % prop)
    print("  %s" % mets)
    print("\nACTIVATION IS UNSCREENED for every design above. Nothing we have tests it.")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit("usage: select_v2.py <designs.json> [more.json ...]")
    main(sys.argv[1:])
