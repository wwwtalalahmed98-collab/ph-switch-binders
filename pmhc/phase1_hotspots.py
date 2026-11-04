"""Choose the hotspots: the spots we instruct the design to touch.

THE QUESTION THIS ANSWERS. Our entire specificity story is "we contact Asp6", the one
residue that differs between the mutant and the wild type. If Asp6 points down into the
groove rather than up out of it, no binder can read it from outside and the premise of
the whole entry is wrong. That has to be measured before a single GPU hour is spent, not
discovered afterwards.

THE ANSWER KEY WE WERE GIVEN. 7OW6 is not just the target - it is the target with a real
T-cell receptor bound to it, and that receptor does exactly what we are being asked to
do: bind this complex and discriminate the mutation. Its footprint is free information
about where a successful binder sits and which parts of the groove it reads.

Using the footprint is legitimate; it is public structural information about the target,
the same as the coordinates themselves. Copying any part of the receptor's SEQUENCE is
not - the challenge requires de novo designs and bans modifying an existing binder. This
script therefore reports contacts and positions only. No receptor sequence is written to
any output, and none is carried into the design stage.

WHY NOT JUST USE ALL 18 GROOVE RESIDUES. ColabDesign averages its contact term over the
hotspot positions, keeping the best-satisfied num_pos of them. Spread the instruction
across 18 positions and a small binder cannot satisfy them, the average stops responding
to anything it can actually do, and the design drifts. On TfR1 that term was averaged
over 13 hotspots with a cutoff that passed almost any distance, and round 1 aimed at
nothing. The set here is deliberately small and weighted onto the peptide.
"""
import json
import os

import numpy as np
from Bio import PDB
from Bio.PDB.SASA import ShrakeRupley

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "data", "pmhc")
OUT = os.path.join(D, "states")

CONTACT = 4.5            # A, heavy atom: what counts as touching
MUT_POS = 6
PEP_G12D = "VVVGADGVGK"

T3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
      "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
      "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
      "TYR": "Y", "VAL": "V"}
# Tien et al. maximum accessible surface, for turning A^2 into a fraction
MAXASA = {"A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "Q": 225, "E": 223,
          "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240,
          "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174}


def residues(chain):
    return [r for r in chain if r.id[0] == " " and r.get_resname() in T3]


def heavy(r):
    return np.array([a.coord for a in r if a.element != "H"])


def strip(model, keep_chains):
    import copy
    sub = copy.deepcopy(model)
    for ch in list(sub):
        if ch.id not in keep_chains:
            sub.detach_child(ch.id)
    return sub


def main():
    man = json.load(open(os.path.join(OUT, "states.json")))
    m6 = PDB.MMCIFParser(QUIET=True).get_structure("7OW6",
                                                   os.path.join(D, "7OW6.cif"))[0]
    hv_ok = set(man["kept"]["A"])
    hv = [r for r in residues(m6["A"]) if r.id[1] in hv_ok]
    pep = residues(m6["C"])
    if "".join(T3[r.get_resname()] for r in pep) != PEP_G12D:
        raise SystemExit("peptide is not %s - rebuild the states first" % PEP_G12D)

    # ---- 1. is Asp6 reachable at all? exposure of the peptide in the bare complex ----
    bare = strip(m6, {"A", "B", "C"})       # TCR removed: this is what a binder sees
    ShrakeRupley().compute(bare, level="R")
    exposure = {r.id[1]: r.sasa for r in bare["C"]}
    print("peptide exposure in the complex, with the receptor removed")
    print("  %-4s %-3s %8s %8s" % ("pos", "aa", "A^2", "of max"))
    for r in pep:
        aa = T3[r.get_resname()]
        frac = exposure[r.id[1]] / MAXASA[aa]
        mark = "   <- the mutation" if r.id[1] == MUT_POS else ""
        print("  %-4d %-3s %8.1f %7.0f%%%s" % (r.id[1], aa, exposure[r.id[1]],
                                               100 * frac, mark))
    mut_frac = exposure[MUT_POS] / MAXASA[T3[pep[MUT_POS - 1].get_resname()]]
    if mut_frac < 0.05:
        raise SystemExit("Asp%d is %.0f%% exposed - it cannot be read from outside, and "
                         "the premise of this entry is wrong" % (MUT_POS, 100 * mut_frac))
    print("\nAsp%d is %.0f%% exposed. A binder can reach it." % (MUT_POS, 100 * mut_frac))

    # ---- 2. what does the real receptor touch? positions only, never sequence ----
    tcr = [r for cid in ("D", "E") for r in residues(m6[cid])]
    tcr_atoms = np.concatenate([heavy(r) for r in tcr])
    foot_pep, foot_hv = [], []
    for r in pep:
        d = np.linalg.norm(heavy(r)[:, None, :] - tcr_atoms[None, :, :], axis=-1).min()
        if d <= CONTACT:
            foot_pep.append((r.id[1], T3[r.get_resname()], round(float(d), 2)))
    for r in hv:
        d = np.linalg.norm(heavy(r)[:, None, :] - tcr_atoms[None, :, :], axis=-1).min()
        if d <= CONTACT:
            foot_hv.append((r.id[1], T3[r.get_resname()], round(float(d), 2)))
    print("\nthe bound receptor's footprint at %.1f A:" % CONTACT)
    print("  peptide:      %s" % ", ".join("%s%d" % (a, n) for n, a, _ in foot_pep))
    print("  heavy chain:  %s" % ", ".join("%s%d" % (a, n) for n, a, _ in foot_hv))
    touches_mut = any(n == MUT_POS for n, _, _ in foot_pep)
    print("  touches the mutation: %s" % ("yes" if touches_mut else "NO"))

    # ---- 3. the hotspot set ----
    # The receptor reads peptide 4-9, and all six are exposed enough to use. But using
    # all six would be a trap. Asp6 is only 11% exposed while Gly7 next door is 73%, so
    # a contact term spread over the six would be satisfied most cheaply by piling onto
    # Gly7 and never touching the mutation at all - and a binder that does that reads
    # nothing the wild type does not also have. It would bind both and score zero.
    #
    # So the set is the mutation and its immediate flanks only. Every member is adjacent
    # to Asp6, so a design that satisfies them sits over the mutation rather than beside
    # it. This is the same lesson as TfR1 round 1 - a contact term averaged over too many
    # positions stops pointing anywhere - applied before the run instead of after it.
    usable = [n for n, _, _ in foot_pep
              if exposure[n] / MAXASA[T3[pep[n - 1].get_resname()]] >= 0.05]
    if MUT_POS not in usable:
        raise SystemExit("Asp%d is not in the usable peptide set %s" % (MUT_POS, usable))
    pep_hot = [n for n in (MUT_POS - 1, MUT_POS, MUT_POS + 1) if n in usable]
    escape = max((n for n in usable if n not in pep_hot),
                 key=lambda n: exposure[n], default=None)
    if escape is not None:
        print("\nthe escape route this set closes: position %d is %.0f%% exposed "
              "against Asp%d's %.0f%%"
              % (escape, 100 * exposure[escape] / MAXASA[T3[pep[escape - 1].get_resname()]],
                 MUT_POS, 100 * mut_frac))
    # A few groove residues, taken only where the receptor also touches them, closest to
    # the mutation first. These anchor the binder's register; they are not the point.
    mut_atoms = heavy(pep[MUT_POS - 1])
    near_mut = []
    for n, aa, _ in foot_hv:
        r = next(x for x in hv if x.id[1] == n)
        d = np.linalg.norm(heavy(r)[:, None, :] - mut_atoms[None, :, :], axis=-1).min()
        near_mut.append((round(float(d), 2), n, aa))
    near_mut.sort()
    groove_hot = [n for _, n, _ in near_mut[:4]]

    print("\nchosen hotspots: peptide positions %s, num_pos %d"
          % (pep_hot, len(pep_hot)))
    print("  the mutation is %d; the other two are its flanks" % MUT_POS)
    print("  groove residues %s are recorded for register, NOT used as hotspots"
          % groove_hot)
    print("\n  groove picks, by distance to the mutation:")
    for d, n, aa in near_mut[:6]:
        print("    %s%-4d %5.2f A%s" % (aa, n, d, "" if n in groove_hot else "   (not used)"))

    man["hotspots_chosen"] = dict(
        peptide=pep_hot, num_pos=len(pep_hot), i_con_cutoff=8.0,
        register_only=groove_hot, contact=CONTACT,
        escape_position=escape,
        mutation_exposed_fraction=round(float(mut_frac), 3),
        receptor_footprint=dict(peptide=[n for n, _, _ in foot_pep],
                                heavy=[n for n, _, _ in foot_hv]),
        note="positions only; no receptor sequence is used anywhere in this entry",
    )
    json.dump(man, open(os.path.join(OUT, "states.json"), "w"), indent=1)
    print("\nupdated %s" % os.path.join(OUT, "states.json"))


if __name__ == "__main__":
    main()
