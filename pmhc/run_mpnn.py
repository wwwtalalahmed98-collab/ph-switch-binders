"""Step 3: redesign the docked scaffold binders with ProteinMPNN, locally on CPU.

This runs on a laptop. The only step in this challenge that genuinely needs a GPU is the
three-state structure prediction at the end, so there is no reason to spend a Kaggle
session - or depend on browser automation - for sequence design.

WHAT IT DOES. For each shortlisted docked pose: hold chain A completely fixed (the
179-residue MHC platform followed by the 10-residue KRAS G12D peptide) and rewrite every
residue of chain B, the binder. The backbone never moves; only the sequence changes.
This is step 3 of the pipeline in Liu et al., Science 389:386-391 (2025).

VERIFIED, NOT ASSUMED. A probe run on one pose confirmed the three things this script
depends on, all of which are the kind of layout assumption that has cost this project
three rounds of measurement:
  - the model reports lengths [189, 71], so chain A is the target and chain B the binder
  - fix_pos='A' fixes positions 0-188, i.e. the whole target including the peptide
  - the returned sequence is "chainA/chainB", and chain A comes back byte-identical
    ending in VVVGADGVGK, so splitting on '/' and taking the last field is the binder

CYSTEINE IS CHECKED, NOT TRUSTED. On the previous run this same library printed
`cysteines removed: True` and then emitted sequences containing cysteine; I reported it
as a working feature without looking at the output. Here rm_aa is still requested, and
then every sequence is inspected and the result of that inspection is printed.

Run:  python pmhc/run_mpnn.py
"""
import json
import os
import sys
import time
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "pmhc"))
import screen  # noqa: E402

RESULTS = os.path.join(ROOT, "pmhc", "results")
DOCKED = os.path.join(ROOT, "data", "pmhc", "docked")

N_SEQ = 16                      # sequences per pose per temperature
TEMPERATURES = (0.1, 0.2, 0.3)
KEEP_PER_POSE = 6               # so the final set is not one backbone repeated
PEP = "VVVGADGVGK"
TARGET_LEN = 189                # 179 MHC platform + 10 peptide


def main():
    sl = os.path.join(RESULTS, "shortlist_poses.json")
    if not os.path.exists(sl):
        raise SystemExit("run pmhc/shortlist_poses.py first")
    poses = json.load(open(sl))["poses"]
    print("shortlisted poses to redesign: %d" % len(poses))
    print()
    print("--- settings ---")
    print("    designed chain : B (binder); chain A held fixed")
    print("    temperatures   : %s, %d sequences each = %d per pose"
          % (list(TEMPERATURES), N_SEQ, N_SEQ * len(TEMPERATURES)))
    print("    weights        : v_48_020")
    print("    expected total : %d sequences"
          % (N_SEQ * len(TEMPERATURES) * len(poses)))
    print()

    from colabdesign.mpnn import mk_mpnn_model
    model = mk_mpnn_model(model_name="v_48_020")

    rows = []
    n_seen = n_cys = 0
    n_target_corrupt = 0
    t0 = time.time()

    for i, p in enumerate(poses, 1):
        pdb = os.path.join(DOCKED, "%s_docked.pdb" % p["name"])
        if not os.path.exists(pdb):
            print("  missing %s - skipped" % pdb)
            continue
        model.prep_inputs(pdb_filename=pdb, chain="A,B", fix_pos="A", rm_aa="C")
        for temp in TEMPERATURES:
            out = model.sample(num=N_SEQ, batch=1, temperature=temp)
            for k in range(len(out["seq"])):
                full = str(out["seq"][k])
                n_seen += 1
                if "/" not in full:
                    n_target_corrupt += 1
                    continue
                tgt, seq = full.split("/")[0], full.split("/")[-1]
                # The target must come back untouched. If it does not, fix_pos did
                # not hold and every number downstream would be meaningless.
                if len(tgt) != TARGET_LEN or not tgt.endswith(PEP):
                    n_target_corrupt += 1
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
                r["seqid_to_native"] = round(float(out["seqid"][k]), 4)
                r["parent_identity"] = screen.identity(seq, p["scaffold_seq"])
                for key in ("n_mut_contacts", "n_pep_contacts", "n_mhc_contacts",
                            "pep_fraction", "binder_len", "align_rmsd"):
                    r[key] = p[key]
                rows.append(r)
        if i % 10 == 0 or i == len(poses):
            print("  %d/%d poses redesigned  (%.0f s elapsed)"
                  % (i, len(poses), time.time() - t0), flush=True)

    print()
    print("--- what came back ---")
    print("    sequences generated      : %d" % n_seen)
    print("    target chain corrupted   : %d  %s"
          % (n_target_corrupt,
             "(fix_pos held)" if n_target_corrupt == 0 else "<-- fix_pos FAILED"))
    print("    dropped for cysteine     : %d  -> rm_aa='C' %s"
          % (n_cys, "WORKED" if n_cys == 0
             else "WAS IGNORED (checked, not trusted)"))
    print("    sequences screened       : %d" % len(rows))
    if not rows:
        raise SystemExit("nothing survived collection")

    tiers = Counter(r["tier"] for r in rows)
    print()
    print("--- buildability tier (denominator %d) ---" % len(rows))
    for t in ("good", "marginal", "fail", "fold-breaker"):
        n = tiers.get(t, 0)
        print("    %-14s %5d  (%4.1f%%)" % (t, n, 100.0 * n / len(rows)))
    print()
    print("    The previous approach produced 720 sequences across 15 hallucinated")
    print("    backbones with ZERO good and ZERO marginal.")

    ident = sorted(r["parent_identity"] for r in rows)
    print()
    print("--- identity to the parent scaffold (for the de novo declaration) ---")
    print("    min %.2f  median %.2f  max %.2f"
          % (ident[0], ident[len(ident) // 2], ident[-1]))

    # selection
    print()
    print("--- selection rule, stated before the numbers are read ---")
    print("    1. tier 'good' before 'marginal'; 'fail'/'fold-breaker' discarded")
    print("    2. within a tier, lower ProteinMPNN score")
    print("    3. at most %d per pose" % KEEP_PER_POSE)
    rank = {"good": 0, "marginal": 1}
    usable = [r for r in rows if r["tier"] in rank]
    if not usable:
        raise SystemExit("nothing reached 'good' or 'marginal' - the scaffold route "
                         "has failed the same test the hallucination route did")
    by_pose = {}
    for r in usable:
        by_pose.setdefault(r["parent"], []).append(r)
    picked = []
    for parent, rs in by_pose.items():
        rs.sort(key=lambda r: (rank[r["tier"]], r["mpnn_score"], r["distance"]))
        picked.extend(rs[:KEEP_PER_POSE])
    picked.sort(key=lambda r: (rank[r["tier"]], -r["n_mut_contacts"],
                               -r["pep_fraction"], r["mpnn_score"]))
    for i, r in enumerate(picked, 1):
        r["design"] = "g12d_%s_%02d" % (r["parent"], i)
    print()
    print("    poses contributing         : %d of %d" % (len(by_pose), len(poses)))
    print("    candidates carried forward : %d" % len(picked))

    print()
    print("--- top 30 candidates ---")
    print("%-24s %-8s %4s %4s %6s %6s %7s"
          % ("design", "tier", "len", "mut", "pepfr", "ident", "mpnn"))
    for r in picked[:30]:
        print("%-24s %-8s %4d %4d %6.2f %6.2f %7.3f"
              % (r["design"], r["tier"], r["length"], r["n_mut_contacts"],
                 r["pep_fraction"], r["parent_identity"], r["mpnn_score"]))

    json.dump(picked, open(os.path.join(RESULTS, "mpnn_candidates.json"), "w"),
              indent=1)
    json.dump(rows, open(os.path.join(RESULTS, "mpnn_all.json"), "w"), indent=1)
    print()
    print("wrote pmhc/results/mpnn_candidates.json (%d) and mpnn_all.json (%d)"
          % (len(picked), len(rows)))
    print()
    print("NOT DONE AND NOT VERIFIED:")
    print("  - No structure prediction has been run on any of these sequences.")
    print("  - No specificity measured. Whether a design prefers Asp6 over Gly6 is")
    print("    the next step and it is what the competition ranks first.")
    print("  - A buildability tier says 'could fold, could be expressed'. It says")
    print("    nothing about binding.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
