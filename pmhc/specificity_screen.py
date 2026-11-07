"""Does the design actually READ the G12D mutation? The screen we have never had.

WHY THIS IS THE MOST IMPORTANT SCRIPT IN THE CHALLENGE. Track 3 ranks entries on
KD(wild type) / KD(G12D) first, KD(apo) / KD(G12D) second, and raw affinity only as a
tie-break. In other words we are scored on discrimination, not strength. Across four
challenges this project has never once measured discrimination. Everything we have
reported - pLDDT, i_pTM, interface PAE - answers "does it bind", not "does it bind THIS
rather than that".

THE METHOD, from Liu et al., Science 389:386-391 (2025):

    "We then used ProteinMPNN to predict the sequence of the target peptide in the
     context of the designed binder and selected those designs for which the
     confidence of the on-target peptide sequence was higher than off-target
     sequences"

Turned around for our target: hold the binder and the MHC fixed, rub out the peptide's
identity, and ask ProteinMPNN what peptide belongs in that pocket. A binder that truly
grips the mutated aspartate will push the model towards Asp at peptide position 6. A
binder merely resting on the groove will not care which residue is there.

The number we take is

    margin = log P(Asp at peptide position 6) - log P(Gly at peptide position 6)

Asp is G12D, Gly is wild type. A positive margin means the pocket prefers the mutant.

THE CONTROL THAT MAKES IT MEAN ANYTHING. The MHC groove has its own preferences at
position 6, with no binder present at all. Without measuring that, a large margin might
be the MHC talking, not our design. So the same measurement is run on the bare target
and every design is reported as a DELTA against that baseline. The baseline is printed
first, before any design scores, so it cannot be chosen after the fact.

WHAT THIS IS NOT. It is a sequence model's opinion about a static backbone, not a
binding free energy. It does not predict KD, and a positive margin is not a promise of
experimental specificity. It is the cheapest available estimate of the quantity we are
ranked on, and it costs no GPU time.

Run:  python pmhc/specificity_screen.py [--limit N]
"""
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "pmhc"))

RESULTS = os.path.join(ROOT, "pmhc", "results")
DOCKED = os.path.join(ROOT, "data", "pmhc", "docked")
WORK = os.path.join(ROOT, "data", "pmhc", "spec")

# ColabDesign's ProteinMPNN alphabet.
ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"
AA_D = ALPHABET.index("D")
AA_G = ALPHABET.index("G")

PEP = "VVVGADGVGK"
N_MHC = 179                      # chain A residues 1..179 are the MHC platform
PEP_START = N_MHC                # 0-based index of peptide residue 1
MUT_POS = 6                      # 1-based position within the peptide
MUT_IDX = PEP_START + MUT_POS - 1   # 0-based index of the Asp, = 184

BACKBONE = ("N", "CA", "C", "O")
T1 = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
      "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
      "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
      "Y": "TYR", "V": "VAL"}


def thread_sequence(pose_pdb, binder_seq, out_pdb, drop_binder=False):
    """Write a backbone-only PDB carrying `binder_seq` on chain B.

    Backbone-only is correct, not a shortcut: ProteinMPNN reads N, CA, C, O and a CB it
    computes itself from those three, so side-chain coordinates are never used. Writing
    only the backbone means we do not have to invent side-chain positions for a
    sequence that has never been packed.

    drop_binder=True writes the target alone, which is the no-binder control.
    """
    lines = []
    serial = 0
    bi = 0
    for line in open(pose_pdb):
        if not line.startswith("ATOM"):
            continue
        name = line[12:16].strip()
        if name not in BACKBONE:
            continue
        chain = line[21]
        if chain == "B":
            if drop_binder:
                continue
            # residue number within chain B, 1-based in our docked files
            resnum = int(line[22:26])
            if resnum - 1 >= len(binder_seq):
                raise SystemExit("binder sequence shorter than the pose backbone "
                                 "(%d residues vs %d)" % (len(binder_seq), resnum))
            resname = T1[binder_seq[resnum - 1]]
            line = line[:17] + resname + line[20:]
        serial += 1
        lines.append(line[:6] + "%5d" % serial + line[11:])
        bi += 1
    with open(out_pdb, "w") as f:
        f.writelines(lines)
        f.write("END\n")
    return out_pdb


def _logprobs(model, pdb, chain, fix_pos, n_samples):
    """Mean log-probability per position over n_samples decodings.

    fix_pos names the positions to DESIGN (inverse=True inverts the set, verified in
    colabdesign/mpnn/model.py line 76-79). Fixed positions get a +1e7 one-hot bias, so
    only the designed positions carry a meaningful distribution - which is why every
    number read out of this is taken at a designed position.
    """
    model.prep_inputs(pdb_filename=pdb, chain=chain, fix_pos=fix_pos, inverse=True)
    out = model.sample(num=n_samples, batch=1, temperature=0.1)
    logits = np.asarray(out["logits"])          # (n, L, 21)
    lp = logits - np.log(np.exp(logits).sum(-1, keepdims=True))
    return lp.mean(axis=0)


def mutation_logprobs(model, pdb, chain="A,B", n_samples=4):
    """P(amino acid at peptide position 6 | everything else), the primary measure.

    Only residue MUT_IDX is left free. The other nine peptide residues, the whole MHC
    and the whole binder are all context. That makes this a clean conditional
    probability for the one position that distinguishes G12D from wild type, rather
    than a number that depends on what the model happened to sample at the other
    peptide positions first.
    """
    return _logprobs(model, pdb, chain, "A%d" % (MUT_IDX + 1), n_samples)


def peptide_logprobs(model, pdb, chain="A,B", n_samples=8):
    """P at every peptide position with the whole peptide free. Secondary sanity
    number: how well the pocket recovers the real peptide as a whole."""
    return _logprobs(model, pdb, chain,
                     "A%d-%d" % (PEP_START + 1, PEP_START + len(PEP)), n_samples)


def check_indexing(model, pdb):
    """Prove MUT_IDX points at the aspartate before trusting any number from it.

    If this index were wrong we would be reading some other residue's preferences and
    the whole screen would produce confident, meaningless numbers - the exact failure
    this project has hit three times.
    """
    model.prep_inputs(pdb_filename=pdb, chain="A,B", fix_pos="A")
    out = model.sample(num=1, batch=1, temperature=0.1)
    full = str(out["seq"][0])
    tgt = full.split("/")[0]
    if len(tgt) != N_MHC + len(PEP):
        raise SystemExit("target chain is %d residues, expected %d"
                         % (len(tgt), N_MHC + len(PEP)))
    if not tgt.endswith(PEP):
        raise SystemExit("target chain does not end in %s - layout changed" % PEP)
    if tgt[MUT_IDX] != "D":
        raise SystemExit("index %d is %r, expected the mutation 'D'"
                         % (MUT_IDX, tgt[MUT_IDX]))
    print("    index check: position %d of the target chain is %r (the G12D "
          "aspartate)" % (MUT_IDX, tgt[MUT_IDX]))
    print("    peptide read back: %s" % tgt[PEP_START:])
    return True


def main(argv):
    limit = None
    if "--limit" in argv:
        limit = int(argv[argv.index("--limit") + 1])

    cand_path = os.path.join(RESULTS, "mpnn_candidates.json")
    if not os.path.exists(cand_path):
        raise SystemExit("run pmhc/run_mpnn.py first (missing %s)" % cand_path)
    cands = json.load(open(cand_path))
    if limit:
        cands = cands[:limit]
    os.makedirs(WORK, exist_ok=True)

    print("candidates to screen: %d" % len(cands))
    print()
    from colabdesign.mpnn import mk_mpnn_model
    model = mk_mpnn_model(model_name="v_48_020")

    print("--- verifying the layout before measuring anything ---")
    first_pose = os.path.join(DOCKED, "%s_docked.pdb" % cands[0]["parent"])
    check_indexing(model, first_pose)
    print()

    # ---- the control, computed and printed BEFORE any design score ----
    print("--- baseline: the bare target, no binder present ---")
    bare = thread_sequence(first_pose, "", os.path.join(WORK, "_bare.pdb"),
                           drop_binder=True)
    lp0 = mutation_logprobs(model, bare, chain="A")
    base_d, base_g = float(lp0[MUT_IDX, AA_D]), float(lp0[MUT_IDX, AA_G])
    baseline = base_d - base_g
    print("    log P(Asp6) = %+.3f   log P(Gly6) = %+.3f" % (base_d, base_g))
    print("    baseline margin (MHC groove alone) = %+.3f" % baseline)
    print()
    print("    Every design below is reported as a DELTA against this number. A")
    print("    design that adds nothing scores about zero no matter how large its")
    print("    own margin looks.")
    print()

    rows = []
    t0 = time.time()
    for i, c in enumerate(cands, 1):
        pose = os.path.join(DOCKED, "%s_docked.pdb" % c["parent"])
        if not os.path.exists(pose):
            continue
        p = thread_sequence(pose, c["seq"],
                            os.path.join(WORK, "%s.pdb" % c["design"]))
        lp = mutation_logprobs(model, p)
        d, g = float(lp[MUT_IDX, AA_D]), float(lp[MUT_IDX, AA_G])
        r = dict(c)
        r["logp_asp6"] = round(d, 4)
        r["logp_gly6"] = round(g, 4)
        r["margin"] = round(d - g, 4)
        r["delta_vs_baseline"] = round((d - g) - baseline, 4)
        rows.append(r)
        if i % 10 == 0 or i == len(cands):
            print("  %d/%d screened  (%.0f s)" % (i, len(cands), time.time() - t0),
                  flush=True)

    if not rows:
        raise SystemExit("nothing screened")

    deltas = sorted(r["delta_vs_baseline"] for r in rows)
    print()
    print("--- did the screen separate anything? (denominator %d) ---" % len(rows))
    print("    delta vs baseline: min %+.3f  median %+.3f  max %+.3f"
          % (deltas[0], deltas[len(deltas) // 2], deltas[-1]))
    pos = sum(1 for d in deltas if d > 0)
    print("    designs that improve on the bare groove: %d of %d" % (pos, len(rows)))
    spread = deltas[-1] - deltas[0]
    print("    spread: %.3f" % spread)
    if spread < 0.2:
        print()
        print("    WARNING: the scores are essentially flat. This screen is then not")
        print("    discriminating between our designs and must not be used to rank")
        print("    them. Report it as a null result rather than reading noise.")

    rows.sort(key=lambda r: -r["delta_vs_baseline"])
    print()
    print("--- top 25 by specificity margin over the bare groove ---")
    print("%-24s %8s %8s %8s %8s %5s"
          % ("design", "logP(D)", "logP(G)", "margin", "delta", "mut"))
    for r in rows[:25]:
        print("%-24s %+8.3f %+8.3f %+8.3f %+8.3f %5d"
              % (r["design"], r["logp_asp6"], r["logp_gly6"], r["margin"],
                 r["delta_vs_baseline"], r["n_mut_contacts"]))

    out = os.path.join(RESULTS, "specificity.json")
    json.dump({"baseline_margin": baseline,
               "baseline_logp_asp6": base_d,
               "baseline_logp_gly6": base_g,
               "designs": rows}, open(out, "w"), indent=1)
    print()
    print("wrote %s" % out)
    print()
    print("WHAT THIS DOES NOT TELL US: it is a sequence model's opinion about a fixed")
    print("backbone. It is not a binding free energy, it does not predict KD, and a")
    print("positive margin is not a promise of experimental specificity. No structure")
    print("has been predicted for any of these sequences yet.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
