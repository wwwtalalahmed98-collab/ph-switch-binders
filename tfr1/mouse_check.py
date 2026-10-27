"""Does each candidate epitope survive the species swap?

EPITOPE_CHOICE.md ranked patches by human/mouse sequence identity. Identity is necessary
but not sufficient: identical residues can sit on a different backbone, point a different
way, or be buried in one species and exposed in the other. This script checks the shape,
not just the letters.

Mouse structure is the AlphaFold DB model for Q62351. That is a prediction, not an
experiment, so every conclusion here is only as good as that model - which is why per-patch
mean pLDDT is reported alongside every number and a low-confidence patch is called out.

Comparison is like for like: human SASA is recomputed on the human MONOMER here, because
the mouse model is a monomer. (EPITOPE_CHOICE.md used dimer-context SASA, which is the
right number for "can a binder reach it" but the wrong one for comparing against a
monomeric model.)
"""
import io
import json
import os

import numpy as np
import requests
from Bio import PDB
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
from Bio.PDB.SASA import ShrakeRupley

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "data", "tfr1")
MOUSE_PDB = os.path.join(D, "mouse_AF_Q62351.pdb")
ECTO = (121, 760)

THREE2ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
             "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
             "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
             "TYR": "Y", "VAL": "V"}
MAXASA = {"A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "Q": 225, "E": 223,
          "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240,
          "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174}


def fasta(p):
    return "".join(l.strip() for l in io.open(p) if not l.startswith(">"))


def hm_map(h, m):
    al = PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score = -11, -1
    a = al.align(h, m)[0]
    out = {}
    for (hs, he), (ms, me) in zip(a.aligned[0], a.aligned[1]):
        for k in range(he - hs):
            out[hs + k + 1] = ms + k + 1
    return out


def load(path, parser, chain=None):
    st = parser.get_structure("x", path)[0]
    ch = [c for c in st if chain is None or c.id == chain]
    return st, ch


def residues(chain):
    out = {}
    for r in chain:
        if r.id[0] == " " and r.get_resname() in THREE2ONE:
            out[r.id[1]] = r
    return out


def main():
    if not os.path.exists(MOUSE_PDB):
        u = requests.get("https://alphafold.ebi.ac.uk/api/prediction/Q62351",
                         timeout=40).json()[0]["pdbUrl"]
        io.open(MOUSE_PDB, "w", newline="\n").write(requests.get(u, timeout=120).text)
        print("downloaded %s" % os.path.basename(MOUSE_PDB))

    human = fasta(os.path.join(D, "P02786.fasta"))
    mouse = fasta(os.path.join(D, "Q62351.fasta"))
    h2m = hm_map(human, mouse)

    hst, hch = load(os.path.join(D, "11QC.cif"), PDB.MMCIFParser(QUIET=True), "A")
    mst, mch = load(MOUSE_PDB, PDB.PDBParser(QUIET=True))
    hres, mres = residues(hch[0]), residues(mch[0])

    # monomer SASA for both, so the two are comparable
    for st in (hst, mst):
        ShrakeRupley().compute(st, level="R")

    def rel(r):
        return r.sasa / MAXASA[THREE2ONE[r.get_resname()]]

    plddt = {n: float(np.mean([a.get_bfactor() for a in r])) for n, r in mres.items()}

    # global superposition on all ectodomain CA pairs present in both
    pairs = [(hres[n], mres[h2m[n]]) for n in sorted(hres)
             if ECTO[0] <= n <= ECTO[1] and n in h2m and h2m[n] in mres
             and "CA" in hres[n] and "CA" in mres[h2m[n]]]
    sup = PDB.Superimposer()
    sup.set_atoms([p[0]["CA"] for p in pairs], [p[1]["CA"] for p in pairs])
    print("global superposition over %d CA pairs: RMSD %.2f A" % (len(pairs), sup.rms))
    sup.apply([a for r in mch[0] for a in r])

    sites = {}
    # Prefer the combined funnel (conservation + glycans + transferrin site) when it exists;
    # ranking on conservation alone is what sent the earlier passes to sites that later
    # turned out to be glycan-shadowed or on top of the transferrin footprint.
    sl = os.path.join(D, "tfr1_epitope_shortlist.json")
    if os.path.exists(sl):
        for p in json.load(open(sl))["sites"][:10]:
            sites["%s%d" % (p["centre_aa"], p["centre"])] = p["members"]
    else:
        for p in json.load(open(os.path.join(D, "tfr1_patches_r14.json")))[:8]:
            sites["%s%d" % (p["centre_aa"], p["centre"])] = p["members"]
    sites["APICAL TIP (recommended)"] = json.load(
        open(os.path.join(D, "tfr1_patches.json")))["apical"]

    print("\n%-26s %4s %5s %6s %6s %6s %6s  %s"
          % ("patch", "n", "cons", "RMSD", "pLDDT", "hSASA", "mSASA", "substitutions"))
    rows = []
    for name, mem in sites.items():
        mem = [n for n in mem if n in hres and n in h2m and h2m[n] in mres]
        if len(mem) < 4:
            continue
        ha = [hres[n]["CA"] for n in mem if "CA" in hres[n]]
        ma = [mres[h2m[n]]["CA"] for n in mem if "CA" in mres[h2m[n]]]
        # local superposition: is the patch the same SHAPE, independent of placement
        loc = PDB.Superimposer()
        loc.set_atoms(ha, ma)
        subs = [(n, THREE2ONE[hres[n].get_resname()],
                 THREE2ONE[mres[h2m[n]].get_resname()]) for n in mem
                if THREE2ONE[hres[n].get_resname()] != THREE2ONE[mres[h2m[n]].get_resname()]]
        row = dict(patch=name, n=len(mem),
                   cons=100.0 * (len(mem) - len(subs)) / len(mem),
                   rmsd=float(loc.rms),
                   plddt=float(np.mean([plddt[h2m[n]] for n in mem])),
                   hsasa=float(np.mean([rel(hres[n]) for n in mem])),
                   msasa=float(np.mean([rel(mres[h2m[n]]) for n in mem])),
                   subs=["%s%d%s" % (a, n, b) for n, a, b in subs],
                   members=mem)
        rows.append(row)
        print("%-26s %4d %4.0f%% %6.2f %6.1f %6.2f %6.2f  %s"
              % (name, row["n"], row["cons"], row["rmsd"], row["plddt"],
                 row["hsasa"], row["msasa"], ", ".join(row["subs"])[:60]))

    json.dump(rows, open(os.path.join(D, "tfr1_mouse_check.json"), "w"), indent=0)
    print("\nwrote data/tfr1/tfr1_mouse_check.json")
    print("RMSD is after superposing that patch alone, so it measures shape, not position.")
    print("pLDDT below ~70 means the mouse model is unsure there and the row is weak.")


if __name__ == "__main__":
    main()
