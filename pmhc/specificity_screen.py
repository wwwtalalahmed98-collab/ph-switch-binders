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


def make_wt_pose(pose_pdb, out_pdb):
    """The same pose with the peptide's Asp6 renamed to Gly: the wild-type state.

    This is honest only because the files involved are backbone-only as far as
    ProteinMPNN is concerned - it reads N, CA, C, O and computes CB itself
    (colabdesign/mpnn/model.py line 56). Renaming the residue therefore produces
    exactly the wild-type input the model would see, with no invented coordinates.

    The limitation is the one already recorded for the counter-states: wild type and
    G12D share a single backbone here, the one crystallised with the aspartate. That
    makes the comparison cleaner than reality, so a margin measured this way should be
    read as an upper bound on how easy the discrimination is.
    """
    target_resnum = PEP_START + MUT_POS          # 185 in chain A, 1-based
    n_renamed = 0
    out = []
    for line in open(pose_pdb):
        if line.startswith("ATOM") and line[21] == "A" \
                and int(line[22:26]) == target_resnum:
            if line[17:20].strip() != "ASP":
                raise SystemExit("residue %d of chain A is %s, expected ASP"
                                 % (target_resnum, line[17:20]))
            if line[12:16].strip() not in BACKBONE:
                continue                          # drop the aspartate side chain
            line = line[:17] + "GLY" + line[20:]
            n_renamed += 1
        out.append(line)
    if n_renamed == 0:
        raise SystemExit("no atoms renamed at chain A residue %d" % target_resnum)
    open(out_pdb, "w").writelines(out)
    return out_pdb, n_renamed


def binder_fit(model, pose_pdb, binder_seq):
    """ProteinMPNN score of the binder sequence given this target background.

    Score is a mean negative log-likelihood, so LOWER is a better fit. Comparing the
    same binder sequence on the G12D background and the wild-type background asks the
    complementary question to the peptide-recovery test: not "which peptide does this
    pocket want" but "does this binder sequence sit more comfortably against the
    mutant than against wild type". Two independent reads on the same claim are worth
    more than one, especially when neither is a binding free energy.
    """
    model.prep_inputs(pdb_filename=pose_pdb, chain="A,B", fix_pos="A")
    out = model.score(seq=binder_seq)
    return float(np.asarray(out["score"]).reshape(-1)[0])


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

    # One wild-type background per parent pose, built once and reused.
    wt_poses = {}
    for c in cands:
        if c["parent"] in wt_poses:
            continue
        pose = os.path.join(DOCKED, "%s_docked.pdb" % c["parent"])
        if not os.path.exists(pose):
            continue
        wt, _ = make_wt_pose(pose, os.path.join(WORK, "%s_wt.pdb" % c["parent"]))
        wt_poses[c["parent"]] = wt
    print("--- wild-type backgrounds built: %d (Asp6 -> Gly6 on the same backbone) ---"
          % len(wt_poses))
    print()

    # Work parent by parent. The binder-fit measure only needs the model prepared on
    # the parent pose, not on each design, so preparing once per parent and scoring
    # every design of that parent against it turns 2 preparations per design into 2
    # per parent. With ~6 designs per parent that is most of the run time.
    by_parent = {}
    for c in cands:
        by_parent.setdefault(c["parent"], []).append(c)

    rows = []
    done = 0
    t0 = time.time()
    for parent, group in by_parent.items():
        pose = os.path.join(DOCKED, "%s_docked.pdb" % parent)
        if not os.path.exists(pose) or parent not in wt_poses:
            continue

        # (a) binder fit, both backgrounds, one preparation each
        fits = {}
        try:
            model.prep_inputs(pdb_filename=pose, chain="A,B", fix_pos="A")
            mut_scores = {c["design"]: float(
                np.asarray(model.score(seq=c["seq"])["score"]).reshape(-1)[0])
                for c in group}
            model.prep_inputs(pdb_filename=wt_poses[parent], chain="A,B",
                              fix_pos="A")
            wt_scores = {c["design"]: float(
                np.asarray(model.score(seq=c["seq"])["score"]).reshape(-1)[0])
                for c in group}
            fits = {k: (mut_scores[k], wt_scores[k]) for k in mut_scores}
        except Exception as e:                     # noqa: BLE001
            if not rows:
                print("    binder-fit measure unavailable: %s: %s"
                      % (type(e).__name__, e))

        # (b) peptide recovery, which does need the design threaded onto the backbone
        for c in group:
            p = thread_sequence(pose, c["seq"],
                                os.path.join(WORK, "%s.pdb" % c["design"]))
            lp = mutation_logprobs(model, p)
            d, g = float(lp[MUT_IDX, AA_D]), float(lp[MUT_IDX, AA_G])
            r = dict(c)
            r["logp_asp6"] = round(d, 4)
            r["logp_gly6"] = round(g, 4)
            r["margin"] = round(d - g, 4)
            r["delta_vs_baseline"] = round((d - g) - baseline, 4)
            if c["design"] in fits:
                s_mut, s_wt = fits[c["design"]]
                r["binder_score_g12d"] = round(s_mut, 4)
                r["binder_score_wt"] = round(s_wt, 4)
                # score is a negative log-likelihood, so lower is a better fit and a
                # POSITIVE difference means the binder suits the mutant background.
                r["binder_fit_delta"] = round(s_wt - s_mut, 4)
            else:
                r["binder_fit_delta"] = None
            rows.append(r)
            done += 1
            if done % 20 == 0 or done == len(cands):
                print("  %d/%d screened  (%.0f s)" % (done, len(cands),
                                                      time.time() - t0), flush=True)

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

    # The second, independent measure, and whether the two agree.
    fits = [r["binder_fit_delta"] for r in rows if r.get("binder_fit_delta") is not None]
    if fits:
        fs = sorted(fits)
        print()
        print("--- second measure: binder fit on mutant vs wild-type background ---")
        print("    delta (positive = fits the mutant better):"
              " min %+.3f  median %+.3f  max %+.3f"
              % (fs[0], fs[len(fs) // 2], fs[-1]))
        print("    designs fitting the mutant better: %d of %d"
              % (sum(1 for x in fits if x > 0), len(fits)))
        both = [r for r in rows if r.get("binder_fit_delta") is not None]
        if len(both) > 2:
            a = np.array([r["delta_vs_baseline"] for r in both])
            b = np.array([r["binder_fit_delta"] for r in both])
            if a.std() > 1e-9 and b.std() > 1e-9:
                cc = float(np.corrcoef(a, b)[0, 1])
                print("    correlation between the two measures: %+.2f" % cc)
                print("    They ask different questions (which peptide does the")
                print("    pocket want, vs which background suits this binder), so")
                print("    agreement is evidence and disagreement is a caution, not")
                print("    a bug.")

    rows.sort(key=lambda r: -r["delta_vs_baseline"])
    print()
    print("--- top 25 by specificity margin over the bare groove ---")
    print("%-24s %8s %8s %8s %8s %8s %5s"
          % ("design", "logP(D)", "logP(G)", "margin", "delta", "fitdiff", "mut"))
    for r in rows[:25]:
        fd = r.get("binder_fit_delta")
        print("%-24s %+8.3f %+8.3f %+8.3f %+8.3f %8s %5d"
              % (r["design"], r["logp_asp6"], r["logp_gly6"], r["margin"],
                 r["delta_vs_baseline"],
                 ("%+.3f" % fd) if fd is not None else "-", r["n_mut_contacts"]))

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
