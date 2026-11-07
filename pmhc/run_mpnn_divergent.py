"""A second redesign pass that is pushed AWAY from the parent scaffold sequence.

WHY. The first pass produced good sequences but at 50-60% identity to the scaffold
each backbone came from (min 0.23, median 0.55, max 0.78). That is not a bug: given a
fixed backbone, ProteinMPNN recovers a largely native-like sequence, because the
backbone is most of what determines the sequence. But Track 3 says "de novo only;
modifying an existing binder is not allowed", and a sequence 55% identical to a
published designed binder is uncomfortably close to the thing that rule prohibits.

Raising the sampling temperature barely helps - at T=0.3 the median was still 0.53 -
because temperature adds noise everywhere rather than specifically discouraging the
native residue.

WHAT THIS DOES. It applies a per-position penalty to the parent scaffold's own amino
acid at every designed position, so ProteinMPNN has to find a different residue unless
the backbone really insists. ColabDesign exposes the bias array that makes this
possible: prep_inputs builds `_inputs["bias"]` with shape (L, 21), and already uses it
to pin fixed positions and to implement rm_aa, so subtracting from the native entry at
each designed position is using the mechanism as intended rather than reaching around
it.

THE TRADE, STATED PLAINLY. Pushing away from the native sequence costs sequence
quality: the further from what the backbone wants, the worse ProteinMPNN's own score
and the more likely the design does not fold. That is why this is a SECOND pass rather
than a replacement. The final selection draws from both runs, and the trade is visible
because both the identity and the ProteinMPNN score are reported for every sequence.

WHAT IT DOES NOT FIX. Lowering identity from 0.55 to 0.35 does not settle the rule
question. The backbone still came from a published scaffold. Only changing the backbone
itself - partial diffusion, which the source paper uses and which we have not got
running - would do that. This reduces the exposure; it does not remove it, and the
submission says so.

Run:  python pmhc/run_mpnn_divergent.py [--penalty 2.0] [--poses 30]
"""
import json
import os
import sys
import time
from collections import Counter

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "pmhc"))
import screen  # noqa: E402

RESULTS = os.path.join(ROOT, "pmhc", "results")
DOCKED = os.path.join(ROOT, "data", "pmhc", "docked")

ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"
N_SEQ = 16
TEMPERATURES = (0.2, 0.3)
PEP = "VVVGADGVGK"
TARGET_LEN = 189
SEED = 37


def main(argv):
    penalty = 2.0
    n_poses = 30
    if "--penalty" in argv:
        penalty = float(argv[argv.index("--penalty") + 1])
    if "--poses" in argv:
        n_poses = int(argv[argv.index("--poses") + 1])

    sl = os.path.join(RESULTS, "shortlist_poses.json")
    poses = json.load(open(sl))["poses"]
    # Best poses by the interface criteria, so the limited budget goes where the
    # geometry is already right.
    poses = sorted(poses, key=lambda p: (-p["n_mut_contacts"], -p["pep_fraction"]))
    poses = poses[:n_poses]

    print("--- divergent redesign pass ---")
    print("    poses            : %d (best by Asp6 contact, then peptide share)"
          % len(poses))
    print("    native penalty   : %.2f applied to the parent residue at every"
          % penalty)
    print("                       designed position")
    print("    temperatures     : %s, %d sequences each" % (list(TEMPERATURES), N_SEQ))
    print("    expected total   : %d sequences"
          % (len(poses) * len(TEMPERATURES) * N_SEQ))
    print()

    from colabdesign.mpnn import mk_mpnn_model
    model = mk_mpnn_model(model_name="v_48_020", seed=SEED)

    rows = []
    n_seen = n_cys = n_bad = 0
    t0 = time.time()

    for i, p in enumerate(poses, 1):
        pdb = os.path.join(DOCKED, "%s_docked.pdb" % p["name"])
        if not os.path.exists(pdb):
            continue
        native = p["scaffold_seq"]
        model.prep_inputs(pdb_filename=pdb, chain="A,B", fix_pos="A", rm_aa="C")

        # Penalise the parent's own residue at each designed position. Designed
        # positions are chain B, i.e. the tail of the (189 + len) array.
        bias = np.asarray(model._inputs["bias"]).copy()
        L = bias.shape[0]
        start = L - len(native)
        if start != TARGET_LEN:
            raise SystemExit("expected %d target residues before the binder, got %d"
                             % (TARGET_LEN, start))
        n_applied = 0
        for k, aa in enumerate(native):
            j = ALPHABET.find(aa)
            if j >= 0:
                bias[start + k, j] -= penalty
                n_applied += 1
        if n_applied != len(native):
            raise SystemExit("penalised %d of %d positions" % (n_applied, len(native)))
        model._inputs["bias"] = bias
        if i == 1:
            print("    verified: penalty applied at %d of %d binder positions, "
                  "offset %d" % (n_applied, len(native), start))
            print()

        for temp in TEMPERATURES:
            out = model.sample(num=N_SEQ, batch=1, temperature=temp)
            for k in range(len(out["seq"])):
                full = str(out["seq"][k])
                n_seen += 1
                if "/" not in full:
                    n_bad += 1
                    continue
                tgt, seq = full.split("/")[0], full.split("/")[-1]
                if len(tgt) != TARGET_LEN or not tgt.endswith(PEP):
                    n_bad += 1
                    continue
                if "C" in seq:
                    n_cys += 1
                    continue
                r = screen.report(seq)
                r["name"] = "%s_docked" % p["name"]
                r["parent"] = p["name"]
                r["seq"] = seq
                r["temperature"] = temp
                r["mpnn_score"] = round(float(out["score"][k]), 4)
                r["parent_identity"] = screen.identity(seq, native)
                r["native_penalty"] = penalty
                for key in ("n_mut_contacts", "n_pep_contacts", "n_mhc_contacts",
                            "pep_fraction", "binder_len", "align_rmsd"):
                    r[key] = p[key]
                rows.append(r)
        if i % 10 == 0 or i == len(poses):
            print("  %d/%d poses  (%.0f s)" % (i, len(poses), time.time() - t0),
                  flush=True)

    print()
    print("--- what came back ---")
    print("    sequences generated  : %d" % n_seen)
    print("    target corrupted     : %d" % n_bad)
    print("    dropped for cysteine : %d" % n_cys)
    print("    screened             : %d" % len(rows))
    if not rows:
        raise SystemExit("nothing came back")

    tiers = Counter(r["tier"] for r in rows)
    print()
    print("--- buildability tier (denominator %d) ---" % len(rows))
    for t in ("good", "marginal", "fail", "fold-breaker"):
        n = tiers.get(t, 0)
        print("    %-14s %5d  (%4.1f%%)" % (t, n, 100.0 * n / len(rows)))

    ok = [r for r in rows if r["tier"] in ("good", "marginal")]
    ids = sorted(r["parent_identity"] for r in ok)
    print()
    print("--- did the penalty actually move identity? ---")
    if ids:
        print("    usable sequences : %d" % len(ok))
        print("    identity         : min %.2f  median %.2f  max %.2f"
              % (ids[0], ids[len(ids) // 2], ids[-1]))
        print("    first pass was   : min 0.23  median 0.55  max 0.78")
        for cap in (0.30, 0.35, 0.40, 0.45):
            s = [r for r in ok if r["parent_identity"] < cap]
            print("    below %.2f        : %4d sequences, %d poses"
                  % (cap, len(s), len({r["parent"] for r in s})))
    else:
        print("    NONE - the penalty destroyed sequence quality at this strength.")
        print("    Report that as the result; do not lower the screen to rescue it.")

    sc = [r["mpnn_score"] for r in ok]
    if sc:
        sc.sort()
        print()
        print("    the cost: ProteinMPNN score min %.3f median %.3f max %.3f"
              % (sc[0], sc[len(sc) // 2], sc[-1]))
        print("    (first pass best scores were around 0.83-1.10; higher is worse)")

    out = os.path.join(RESULTS, "mpnn_divergent.json")
    json.dump(rows, open(out, "w"), indent=1)
    print()
    print("wrote %s (%d sequences)" % (out, len(rows)))
    print()
    print("NOT DONE: no specificity measured and no structure predicted for any of")
    print("these. Lowering identity does not settle the de novo question either -")
    print("the backbone still came from a published scaffold.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
