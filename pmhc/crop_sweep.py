"""Compare candidate crops of the pMHC target: size against the cut face each one opens.

The full complex is 384 residues. On TfR1 a 108-residue target gave about 36 designs in
4 h 20 m on two GPUs; 384 is more than three times that and AlphaFold does not scale
kindly, so generating against the full complex would buy very few designs in the time
left. Cropping is therefore not optional - the question is only where to cut.

The cost of cutting is buried hydrophobic surface that becomes exposed and is not exposed
in the real protein. On challenge 3 the optimiser found exactly such a face and dived in.
This prints, for each candidate crop, how much new hydrophobic surface it creates, scaled
against the peptide's own surface so the number means something.

Note the one piece of luck in this challenge: a binder that grips the cut face would grip
it identically in all three states, so the counter-screen would reject it. The cut face
wastes design budget rather than corrupting the result. That is a reason to tolerate a
larger face here than on an affinity-scored challenge, not a reason to ignore it.
"""
import copy
import json
import os

import numpy as np
from Bio import PDB
from Bio.PDB.SASA import ShrakeRupley

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "data", "pmhc")

T3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
      "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
      "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
      "TYR": "Y", "VAL": "V"}
HYDROPHOBIC = set("AVLIMFWCY")
MIN_SEGMENT = 3   # drop dangling 1-2 residue fragments; they are noise, as on TfR1


def residues(chain):
    return [r for r in chain if r.id[0] == " " and r.get_resname() in T3]


def sasa(model, keep):
    sub = copy.deepcopy(model)
    for ch in list(sub):
        if ch.id not in keep:
            sub.detach_child(ch.id)
            continue
        allowed = set(keep[ch.id])
        for r in list(ch):
            if r.id[0] != " " or r.id[1] not in allowed:
                ch.detach_child(r.id)
    ShrakeRupley().compute(sub, level="R")
    return {(ch.id, r.id[1]): r.sasa for ch in sub for r in ch}


def despur(nums):
    """Drop runs shorter than MIN_SEGMENT - isolated residues float free in the crop."""
    nums = sorted(nums)
    out, run = [], []
    for n in nums:
        if run and n == run[-1] + 1:
            run.append(n)
        else:
            if len(run) >= MIN_SEGMENT:
                out += run
            run = [n]
    if len(run) >= MIN_SEGMENT:
        out += run
    return out


def main():
    man = json.load(open(os.path.join(D, "states", "states.json")))
    m6 = PDB.MMCIFParser(QUIET=True).get_structure("7OW6",
                                                   os.path.join(D, "7OW6.cif"))[0]
    hv, b2m, pep = residues(m6["A"]), residues(m6["B"]), residues(m6["C"])
    hv_ok, b2m_ok = set(man["kept"]["A"]), set(man["kept"]["B"])
    hv = [r for r in hv if r.id[1] in hv_ok]
    b2m = [r for r in b2m if r.id[1] in b2m_ok]
    pep_nums = [r.id[1] for r in pep]
    # keyed by (chain, resnum): a crop that removes a3 can leave b2m exposed too, and
    # counting only the heavy chain would hide exactly that.
    aa_of = {("A", r.id[1]): T3[r.get_resname()] for r in hv}
    aa_of.update({("B", r.id[1]): T3[r.get_resname()] for r in b2m})
    aa_of.update({("C", r.id[1]): T3[r.get_resname()] for r in pep})

    full = {"A": [r.id[1] for r in hv], "B": [r.id[1] for r in b2m], "C": pep_nums}
    ref = sasa(m6, full)
    pep_area = sum(ref[k] for k in ref if k[0] == "C")

    pep_atoms = np.array([a.coord for r in pep for a in r if a.element != "H"])

    def within(radius):
        out = []
        for r in hv:
            xyz = np.array([a.coord for a in r if a.element != "H"])
            d = np.linalg.norm(xyz[:, None, :] - pep_atoms[None, :, :], axis=-1).min()
            if d <= radius:
                out.append(r.id[1])
        return despur(out)

    cands = [("full complex", full)]
    cands.append(("platform 1-180", {"A": [n for n in full["A"] if n <= 180],
                                     "C": pep_nums}))
    cands.append(("platform + b2m", {"A": [n for n in full["A"] if n <= 180],
                                     "B": full["B"], "C": pep_nums}))
    for radius in (14.0, 16.0, 18.0, 20.0, 24.0):
        cands.append(("%.0f A around peptide" % radius,
                      {"A": within(radius), "C": pep_nums}))

    print("%-24s %6s %10s %10s" % ("crop", "res", "new phob", "vs peptide"))
    rows = []
    for name, keep in cands:
        s = sasa(m6, keep)
        new = sum(s[k] - ref[k] for k in s
                  if s[k] > ref[k] + 1.0 and aa_of.get(k) in HYDROPHOBIC)
        n = sum(len(v) for v in keep.values())
        rows.append(dict(name=name, residues=n, new_hydrophobic=round(float(new), 1),
                         ratio=round(float(new / pep_area), 2),
                         keep={k: list(v) for k, v in keep.items()}))
        print("%-24s %6d %9.0f %9.2fx" % (name, n, new, new / pep_area))

    print("\npeptide surface in the full complex: %.0f A^2" % pep_area)
    out = os.path.join(D, "states", "crop_sweep.json")
    json.dump(dict(peptide_surface=round(float(pep_area), 1), crops=rows),
              open(out, "w"), indent=1)
    print("wrote %s" % out)


if __name__ == "__main__":
    main()
