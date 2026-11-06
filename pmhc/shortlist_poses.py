"""Choose which docked poses go forward to sequence redesign, and say why.

The filter is written here, before any sequences exist, so that the choice cannot be
tuned after seeing which poses happen to give nice-looking designs. That has been the
failure mode in this project more than once.

THE RULE, fixed in advance:

  1. No backbone clash with the target. Sequence design can replace a side chain; it
     cannot move a backbone. A pose with a backbone clash is unusable, full stop.
  2. At least MIN_MUT heavy-atom contacts to the Asp6 mutation. Specificity for G12D
     over wild type has to come from somewhere, and the only thing that distinguishes
     them is that one side chain. A pose that does not touch it cannot read it.
  3. At least MIN_PEP_FRACTION of the interface on the peptide rather than the MHC.
     This is the Liu et al. criterion: binders with extensive MHC contact bind every
     peptide the MHC presents, which is the opposite of what we are scored on.

Everything passing goes forward. ProteinMPNN is CPU-cheap, so there is no reason to cut
further at this stage - and cutting further would mean choosing on criteria we have not
justified.

Run:  python pmhc/shortlist_poses.py
"""
import json
import os
import sys
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "pmhc", "results")
POSES = os.path.join(RESULTS, "docked_poses.json")

MIN_MUT = 4              # heavy-atom contacts to Asp6
MIN_PEP_FRACTION = 0.20  # share of the interface that is peptide, not MHC
MUT_POS = 6


def main():
    if not os.path.exists(POSES):
        raise SystemExit("run pmhc/dock_scaffolds.py first (missing %s)" % POSES)
    rows = json.load(open(POSES))
    n0 = len(rows)

    # Apply the three criteria one at a time so the attrition is visible.
    step1 = [r for r in rows if r["bb_clashes"] == 0]
    step2 = [r for r in step1 if r["n_mut_contacts"] >= MIN_MUT]
    step3 = [r for r in step2 if r["pep_fraction"] >= MIN_PEP_FRACTION]

    print("--- shortlist, criteria applied in order ---")
    print("    docked poses                                 : %d" % n0)
    print("    ... no backbone clash                        : %d" % len(step1))
    print("    ... and >= %d contacts to Asp6                : %d"
          % (MIN_MUT, len(step2)))
    print("    ... and >= %.2f of interface on the peptide   : %d"
          % (MIN_PEP_FRACTION, len(step3)))
    print()

    if not step3:
        raise SystemExit("nothing survived the pose filter")

    keep = sorted(step3, key=lambda r: (-r["n_mut_contacts"], -r["pep_fraction"]))

    # Length spread matters: the Liu paper found success rates are length dependent,
    # and a shortlist that is all one size samples one geometry.
    lens = [r["binder_len"] for r in keep]
    print("--- what the shortlist looks like ---")
    print("    binder lengths: min %d median %d max %d"
          % (min(lens), sorted(lens)[len(lens) // 2], max(lens)))
    print("    length buckets: %s"
          % sorted(Counter((l // 20) * 20 for l in lens).items()))
    print()

    # Which peptide positions does the shortlist read? If every pose reads the same
    # two residues we have one binding mode, not a diverse set.
    pos_hits = Counter()
    for r in keep:
        for p, n in r["per_peptide_residue"].items():
            if n > 0:
                pos_hits[int(p)] += 1
    print("    peptide positions contacted (of %d poses):" % len(keep))
    pep = "VVVGADGVGK"
    for p in range(1, len(pep) + 1):
        mark = "  <-- the mutation" if p == MUT_POS else ""
        bar = "#" * int(30 * pos_hits.get(p, 0) / len(keep))
        print("      %2d %s  %3d  %-30s%s"
              % (p, pep[p - 1], pos_hits.get(p, 0), bar, mark))
    print()

    print("--- shortlist (all %d going forward to redesign) ---" % len(keep))
    print("    %-10s %5s %4s %5s %5s %6s %6s"
          % ("name", "len", "mut", "pep", "mhc", "pepfr", "rmsd"))
    for r in keep:
        print("    %-10s %5d %4d %5d %5d %6.2f %6.2f"
              % (r["name"], r["binder_len"], r["n_mut_contacts"],
                 r["n_pep_contacts"], r["n_mhc_contacts"], r["pep_fraction"],
                 r["align_rmsd"]))

    out = os.path.join(RESULTS, "shortlist_poses.json")
    json.dump({
        "criteria": {"bb_clashes": 0, "min_mut_contacts": MIN_MUT,
                     "min_pep_fraction": MIN_PEP_FRACTION},
        "attrition": {"docked": n0, "no_bb_clash": len(step1),
                      "and_mut_contact": len(step2), "and_pep_centric": len(step3)},
        "poses": keep,
    }, open(out, "w"), indent=1)
    print()
    print("wrote %s" % out)
    print()
    print("NOT DONE: these still carry their original scaffold sequences. The next")
    print("step replaces every amino acid in the binder with ProteinMPNN.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
