"""Choose where to put histidines on a binder so that it lets go in acid.

STAGE A of the pH-switch pass, and the only stage that is free. GPU time is the binding
constraint this week, so this narrows the search on CPU before anything is predicted:
a 70-residue binder has 70 positions, and predicting a histidine at every one of them,
then every useful pair, is thousands of AlphaFold runs we cannot afford. This picks the
handful worth spending predictions on, and says why for each.

THE MECHANISM, and the direction it has to point
------------------------------------------------
Histidine is the only residue that changes charge between pH 7.4 and 6.0. We want binding
LOST in acid, so each histidine must sit where GAINING a positive charge COSTS binding
energy. Two environments do that:

  * buried at the interface against neutral or oily surface - burying a charge is expensive,
    so protonation at pH 6.0 pushes the complex apart;
  * facing one of the target's own positive residues (Arg, Lys) - protonation creates
    like-charge repulsion.

One environment does the OPPOSITE and is the trap: a histidine facing a target Asp or Glu
forms a salt bridge when protonated, which makes binding TIGHTER in acid. That is the FcRn
mechanism and it is what Challenge 1 wanted. Here it is backwards, and it is easy to design
by accident because "histidine near a charged residue" sounds right either way. Positions
like that are reported as WRONG DIRECTION rather than quietly scored low.

HOW MANY, AND HOW FAR APART
---------------------------
The organizers state one titratable site changes affinity at most tenfold per pH unit. From
7.4 to 6.0 is 1.4 units, so one histidine buys at most ~25x and "no detectable binding"
needs ~100x. Two working histidines are the minimum.

But on TNF-alpha we measured that in every histidine pair, one partner did no work. Adjacent
histidines shift each other's pKa, so a pair placed close together silences one of its
members. Hence MIN_SPACING: combinations are only proposed between sites far enough apart to
titrate independently.

WHAT THIS DOES NOT DO. It does not predict the switch. It proposes positions. The prediction
is tfr1/pka_switch.py, run on the structures of the mutants after they are folded.
"""
import argparse
import json
import os

import numpy as np
from Bio import PDB
from Bio.PDB.SASA import ShrakeRupley

CONTACT = 4.5          # A, binder atom to target atom
NEIGHBOUR = 8.0        # A, for reading the chemical neighbourhood across the interface
SALT_BRIDGE = 4.0      # A, close enough that an opposite charge pair is a real interaction
MIN_SPACING = 10.0     # A between CB atoms of two histidines, so they titrate separately
MIN_BURIAL = 15.0      # A^2 of side-chain area buried on binding
TOP_N = 6              # positions proposed per design

THREE2ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
             "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
             "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
             "TYR": "Y", "VAL": "V"}
# Positions we will not touch: Gly and Pro carry backbone roles a histidine cannot take
# over, and an existing His is already what we are trying to create.
SKIP = {"G", "P", "H"}
# Mutating a residue that is itself making a specific interaction removes that interaction.
# A histidine added at the cost of an existing salt bridge is not a switch, it is a broken
# binder: this bit the first test run, which proposed replacing two of barstar's aspartates
# with histidine - the exact residues that hold it onto barnase.
CHARGED = set("DEKR")


def _resolve(model, binder_chain, target_chain):
    ids = [c.id for c in model]
    if binder_chain and target_chain:
        return binder_chain, target_chain
    if len(ids) != 2:
        raise SystemExit("expected two chains, found %s - name them explicitly" % ids)
    a, b = ids
    na = sum(1 for r in model[a] if r.id[0] == " ")
    nb = sum(1 for r in model[b] if r.id[0] == " ")
    return (a, b) if na < nb else (b, a)      # the binder is the shorter chain


def analyse(pdb_path, binder_chain=None, target_chain=None, quiet=False):
    st = PDB.PDBParser(QUIET=True).get_structure("x", pdb_path)
    model = st[0]
    bid, tid = _resolve(model, binder_chain, target_chain)

    # side-chain area buried on binding: complex versus binder alone
    sr = ShrakeRupley()
    sr.compute(model, level="R")
    bound = {r.id[1]: r.sasa for r in model[bid] if r.id[0] == " "}

    free_st = PDB.PDBParser(QUIET=True).get_structure("f", pdb_path)
    for ch in list(free_st[0]):
        if ch.id != bid:
            free_st[0].detach_child(ch.id)
    sr.compute(free_st[0], level="R")
    free = {r.id[1]: r.sasa for r in free_st[0][bid] if r.id[0] == " "}

    tatoms, tinfo = [], []
    for r in model[tid]:
        if r.id[0] != " ":
            continue
        aa = THREE2ONE.get(r.get_resname())
        for a in r:
            if a.element != "H":
                tatoms.append(a.get_coord())
                tinfo.append((aa, r.id[1]))
    tatoms = np.array(tatoms)

    rows = []
    for r in model[bid]:
        if r.id[0] != " " or r.get_resname() not in THREE2ONE:
            continue
        aa = THREE2ONE[r.get_resname()]
        n = r.id[1]
        ratoms = np.array([a.get_coord() for a in r if a.element != "H"])
        d = np.linalg.norm(ratoms[:, None, :] - tatoms[None, :, :], axis=2)
        if d.min() > CONTACT:
            continue                                   # not at the interface at all
        burial = free.get(n, 0.0) - bound.get(n, 0.0)
        near = {tinfo[j] for j in np.where((d < NEIGHBOUR).any(axis=0))[0]}
        basic = sorted(x for x in near if x[0] in "KR")
        acid = sorted(x for x in near if x[0] in "DE")

        # is the wild-type residue itself salt-bridged across the interface?
        bridge = []
        if aa in CHARGED:
            want = "KR" if aa in "DE" else "DE"
            close = np.where((d < SALT_BRIDGE).any(axis=0))[0]
            bridge = sorted({tinfo[j] for j in close if tinfo[j][0] in want})

        reasons, verdict = [], "candidate"
        if bridge:
            verdict = "risky"
            reasons.append("%s is salt-bridged to %s - replacing it removes that contact,"
                           " which costs affinity rather than adding a switch"
                           % (aa, ", ".join("%s%d" % x for x in bridge)))
        if aa in SKIP:
            verdict = "skip"
            reasons.append("%s carries a backbone role or is already His" % aa)
        if burial < MIN_BURIAL:
            verdict = "skip" if verdict == "skip" else "weak"
            reasons.append("buries only %.0f A^2" % burial)
        if acid and not basic:
            verdict = "WRONG DIRECTION"
            reasons.append("faces %s - protonation would form a salt bridge and bind"
                           " TIGHTER in acid" % ", ".join("%s%d" % x for x in acid))
        elif basic:
            reasons.append("faces %s - protonation repels" % ", ".join("%s%d" % x for x in basic))
        if burial >= MIN_BURIAL and not acid and not basic:
            reasons.append("buried against neutral surface")

        cb = r["CB"].get_coord() if "CB" in r else r["CA"].get_coord()
        rows.append(dict(pos=n, wt=aa, burial=round(float(burial), 1),
                         basic=["%s%d" % x for x in basic],
                         acid=["%s%d" % x for x in acid],
                         salt_bridge=["%s%d" % x for x in bridge],
                         verdict=verdict, why="; ".join(reasons),
                         cb=[round(float(v), 2) for v in cb]))

    good = [r for r in rows if r["verdict"] == "candidate"]
    # prefer neutral wild-type residues: swapping one of those for histidine changes the
    # chemistry at that position without deleting an interaction that was already there
    good.sort(key=lambda r: (r["wt"] not in CHARGED, len(r["basic"]) > 0, r["burial"]),
              reverse=True)
    picks = good[:TOP_N]

    # pairs far enough apart to titrate independently
    pairs = []
    for i in range(len(picks)):
        for j in range(i + 1, len(picks)):
            sep = float(np.linalg.norm(np.array(picks[i]["cb"]) - np.array(picks[j]["cb"])))
            if sep >= MIN_SPACING:
                pairs.append(dict(a=picks[i]["pos"], b=picks[j]["pos"],
                                  separation=round(sep, 1),
                                  score=round(picks[i]["burial"] + picks[j]["burial"], 1)))
    pairs.sort(key=lambda p: p["score"], reverse=True)

    if not quiet:
        print("%s: binder chain %s (%d interface residues), target chain %s"
              % (os.path.basename(pdb_path), bid, len(rows), tid))
        print("  %-6s %4s %8s  %s" % ("pos", "wt", "burial", "assessment"))
        for r in sorted(rows, key=lambda r: -r["burial"])[:14]:
            print("  %-6d %4s %8.1f  [%s] %s" % (r["pos"], r["wt"], r["burial"],
                                                 r["verdict"], r["why"]))
        print("  -> proposing %d single placements: %s"
              % (len(picks), ", ".join("%s%dH" % (p["wt"], p["pos"]) for p in picks)))
        print("  -> %d pairs at least %.0f A apart: %s"
              % (len(pairs), MIN_SPACING,
                 ", ".join("%d+%d (%.0f A)" % (p["a"], p["b"], p["separation"])
                           for p in pairs[:5])))
        wrong = [r for r in rows if r["verdict"] == "WRONG DIRECTION"]
        if wrong:
            print("  -> %d position(s) would switch the WRONG WAY and are excluded: %s"
                  % (len(wrong), ", ".join("%s%d" % (r["wt"], r["pos"]) for r in wrong)))
    return dict(pdb=os.path.basename(pdb_path), binder_chain=bid, target_chain=tid,
                interface=rows, singles=picks, pairs=pairs[:8])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("pdb", nargs="+", help="complex structure(s)")
    ap.add_argument("--binder", help="binder chain id")
    ap.add_argument("--target", help="target chain id")
    ap.add_argument("--out", help="write proposals as JSON")
    a = ap.parse_args()
    out = [analyse(p, a.binder, a.target) for p in a.pdb]
    if a.out:
        json.dump(out, open(a.out, "w"), indent=0)
        print("\nwrote %s" % a.out)


if __name__ == "__main__":
    main()
