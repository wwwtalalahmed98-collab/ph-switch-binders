"""Merge the two redesign passes into the candidate set the specificity screen will see.

There are two runs to draw from:

  mpnn_all.json        pass 1, plain ProteinMPNN. Good sequences, but 50-60% identical
                       to the scaffold each backbone came from.
  mpnn_divergent.json  pass 2, with a penalty on the parent's own residue at every
                       designed position. Lower identity, at some cost in sequence
                       quality.

Both are legitimate; they trade different things. This script applies the de novo
identity cap first - that is a hard rule, not a preference - and then takes the best
remaining sequences per pose, so the specificity screen spends its time on candidates
that could actually be submitted.

WHY THE CAP COMES FIRST. If the cap were applied last, the shortlist would fill with
high-identity sequences that score well and then get thrown away, and the screen would
have spent its budget on them. Worse, it would be tempting to raise the cap to keep
them. Ordering the filters this way removes that temptation.

Run:  python pmhc/merge_candidates.py [--cap 0.45] [--per-pose 4] [--max 140]
"""
import json
import os
import sys
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "pmhc", "results")

DEFAULT_CAP = 0.45
DEFAULT_PER_POSE = 4
DEFAULT_MAX = 140          # keeps the specificity screen inside ~40 minutes


def load(path, source):
    if not os.path.exists(path):
        print("    %s: not present" % os.path.basename(path))
        return []
    rows = json.load(open(path))
    for r in rows:
        r["source"] = source
    print("    %-24s %5d sequences" % (os.path.basename(path), len(rows)))
    return rows


def main(argv):
    cap = DEFAULT_CAP
    per_pose = DEFAULT_PER_POSE
    cap_total = DEFAULT_MAX
    if "--cap" in argv:
        cap = float(argv[argv.index("--cap") + 1])
    if "--per-pose" in argv:
        per_pose = int(argv[argv.index("--per-pose") + 1])
    if "--max" in argv:
        cap_total = int(argv[argv.index("--max") + 1])

    print("--- inputs ---")
    rows = []
    rows += load(os.path.join(RESULTS, "mpnn_all.json"), "plain")
    rows += load(os.path.join(RESULTS, "mpnn_divergent.json"), "divergent")
    if not rows:
        raise SystemExit("no redesign output found - run pmhc/run_mpnn.py first")

    print()
    print("--- filters, in this order ---")
    n0 = len(rows)
    rows = [r for r in rows if r["tier"] in ("good", "marginal")]
    n1 = len(rows)
    rows = [r for r in rows if "C" not in r["seq"]]
    n2 = len(rows)
    rows = [r for r in rows if r["parent_identity"] < cap]
    n3 = len(rows)
    # identical sequences can arise from both passes
    seen = set()
    uniq = []
    for r in rows:
        if r["seq"] in seen:
            continue
        seen.add(r["seq"])
        uniq.append(r)
    rows = uniq
    print("    all sequences                         : %d" % n0)
    print("    tier good or marginal                 : %d" % n1)
    print("    no cysteine                           : %d" % n2)
    print("    parent identity below %.2f             : %d" % (cap, n3))
    print("    distinct sequences                    : %d" % len(rows))
    if not rows:
        raise SystemExit("nothing survived the cap of %.2f - either lower the bar "
                         "deliberately and say so, or report that the de novo "
                         "requirement cannot be met this way" % cap)

    by_src = Counter(r["source"] for r in rows)
    print("    contributed by: %s" % dict(by_src))
    print("    distinct poses: %d" % len({r["parent"] for r in rows}))

    rank = {"good": 0, "marginal": 1}
    by_pose = {}
    for r in rows:
        by_pose.setdefault(r["parent"], []).append(r)
    picked = []
    for parent, rs in by_pose.items():
        rs.sort(key=lambda r: (rank[r["tier"]], r["mpnn_score"]))
        picked.extend(rs[:per_pose])

    # Order the overall shortlist by the interface criteria, then trim to budget.
    picked.sort(key=lambda r: (rank[r["tier"]], -r["n_mut_contacts"],
                               -r["pep_fraction"], r["mpnn_score"]))
    if len(picked) > cap_total:
        print()
        print("    trimming %d -> %d to keep the specificity screen inside its"
              " time budget" % (len(picked), cap_total))
        picked = picked[:cap_total]

    for i, r in enumerate(picked, 1):
        r["design"] = "g12d_%s_%s_%03d" % (r["parent"], r["source"][:4], i)

    print()
    print("--- shortlist for the specificity screen ---")
    print("    candidates      : %d" % len(picked))
    print("    distinct poses  : %d" % len({r["parent"] for r in picked}))
    print("    by source       : %s" % dict(Counter(r["source"] for r in picked)))
    print("    by tier         : %s" % dict(Counter(r["tier"] for r in picked)))
    ids = sorted(r["parent_identity"] for r in picked)
    print("    identity        : min %.2f median %.2f max %.2f"
          % (ids[0], ids[len(ids) // 2], ids[-1]))
    lens = sorted(r["length"] for r in picked)
    print("    length          : min %d median %d max %d"
          % (lens[0], lens[len(lens) // 2], lens[-1]))

    out = os.path.join(RESULTS, "mpnn_candidates.json")
    json.dump(picked, open(out, "w"), indent=1)
    print()
    print("wrote %s" % out)
    print()
    print("NOT DONE: no specificity measured on these yet, and no structure")
    print("predicted. The identity cap reduces the de novo exposure; it does not")
    print("remove it, because the backbones are still published scaffolds.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
