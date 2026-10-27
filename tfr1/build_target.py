"""Build the TfR1 design target around the chosen epitope, and prove the crop is honest.

Designs are generated against a cropped target because the full ectodomain dimer is far too
large to run AlphaFold on repeatedly. Cropping is also how Challenge 3 went wrong: cutting a
receptor down created a cavity that did not exist in the real protein, and the optimiser
found it and dived in. At the depth the designs reached, the intact receptor had 211 atoms
and the crop had 4.

TfR1 is soluble and globular rather than a membrane bundle, so there is no channel to fake,
but the same class of error is still available: slicing through the protein exposes interior
surface that is buried in reality, and a designed binder will happily bind to it.

So this does not just cut. It measures. For every residue kept, it compares solvent exposure
in the crop against the same residue in the full dimer, and reports the total NEW hydrophobic
surface the cut created. If that number is large, the crop is lying and the radius should go
up. The number is printed whether or not anyone looks at it, and written to the map file.

Numbering is UniProt throughout, checked residue by residue against P02786 before anything is
written.
"""
import io
import json
import os

import numpy as np
from Bio import PDB
from Bio.PDB.SASA import ShrakeRupley

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "data", "tfr1")
RADIUS = 20.0          # A around the epitope centroid; generous on purpose
MIN_SEGMENT = 3        # drop dangling 1-2 residue fragments, they are noise
SITE = "E559"

THREE2ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
             "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
             "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
             "TYR": "Y", "VAL": "V"}
HYDROPHOBIC = set("AVLIMFWCY")
MAXASA = {"A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "Q": 225, "E": 223,
          "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240,
          "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174}


class _Sel(PDB.Select):
    def __init__(self, chain, keep):
        self.chain, self.keep = chain, set(keep)

    def accept_model(self, m):
        return m.id == 0

    def accept_chain(self, c):
        return c.id == self.chain

    def accept_residue(self, r):
        return r.id[0] == " " and r.id[1] in self.keep

    def accept_atom(self, a):
        return (not a.is_disordered()) or a.get_altloc() in (" ", "A")


def main():
    human = "".join(l.strip() for l in io.open(os.path.join(D, "P02786.fasta"))
                    if not l.startswith(">"))
    short = json.load(open(os.path.join(D, "tfr1_epitope_shortlist.json")))
    site = next(s for s in short["sites"]
                if "%s%d" % (s["centre_aa"], s["centre"]) == SITE)
    epitope = sorted(site["members"])
    print("epitope %s: %d residues  %s" % (SITE, len(epitope), epitope))

    st = PDB.MMCIFParser(QUIET=True).get_structure("x", os.path.join(D, "11QC.cif"))[0]
    chains = [c for c in st if sum(1 for r in c if r.id[0] == " ") > 100]
    ref = chains[0].id

    # --- numbering check against UniProt, before any geometry is used ---
    bad = 0
    for c in chains:
        for r in c:
            if r.id[0] != " " or r.get_resname() not in THREE2ONE:
                continue
            n = r.id[1]
            if 1 <= n <= len(human) and THREE2ONE[r.get_resname()] != human[n - 1]:
                bad += 1
    print("numbering check against P02786: %d mismatches" % bad)
    if bad:
        raise SystemExit("PDB numbering disagrees with UniProt. Stopping.")

    # exposure in the FULL dimer, which is the real molecule
    ShrakeRupley().compute(st, level="R")
    full = {}
    for r in st[ref]:
        if r.id[0] == " " and r.get_resname() in THREE2ONE:
            full[r.id[1]] = r.sasa

    pos = {r.id[1]: r["CA"].get_coord() for r in st[ref]
           if r.id[0] == " " and "CA" in r and r.get_resname() in THREE2ONE}
    centre = np.mean([pos[n] for n in epitope if n in pos], axis=0)

    keep = sorted(n for n, xyz in pos.items()
                  if np.linalg.norm(xyz - centre) <= RADIUS)
    # drop isolated fragments: a 1-2 residue island has no structure to hold it
    # Build a NEW list when a run ends. Appending `cur` and then clearing it stores a
    # reference and then empties what was stored, so every run comes back identical and
    # the segment filter silently deletes almost everything.
    runs, cur = [], [keep[0]]
    for a, b in zip(keep, keep[1:]):
        if b - a <= 2:
            cur.append(b)
        else:
            runs.append(cur)
            cur = [b]
    runs.append(cur)
    keep = sorted(n for run in runs if len(run) >= MIN_SEGMENT for n in run)
    missing = [n for n in epitope if n not in keep]
    if missing:
        raise SystemExit("crop lost epitope residues %s - raise RADIUS" % missing)
    print("crop: %d residues within %.0f A, in %d segments"
          % (len(keep), RADIUS, sum(1 for r in runs if len(r) >= MIN_SEGMENT)))

    out = os.path.join(D, "TfR1_%s_target.pdb" % SITE)
    io_ = PDB.PDBIO()
    io_.set_structure(st)
    io_.save(out, _Sel(ref, keep))

    # --- did the cut invent surface? compare crop against the full molecule ---
    cr = PDB.PDBParser(QUIET=True).get_structure("c", out)
    ShrakeRupley().compute(cr, level="R")
    new_total, new_phobic, worst = 0.0, 0.0, []
    for r in cr[0][ref]:
        n = r.id[1]
        if n not in full:
            continue
        aa = THREE2ONE[r.get_resname()]
        d = r.sasa - full[n]
        if d > 1.0:
            new_total += d
            if aa in HYDROPHOBIC:
                new_phobic += d
            worst.append((d, aa, n))
    worst.sort(reverse=True)
    print("\nsurface the cut created, against the intact dimer:")
    print("  newly exposed area   : %.0f A^2" % new_total)
    print("  of which hydrophobic : %.0f A^2 (%.0f%%)"
          % (new_phobic, 100 * new_phobic / max(new_total, 1)))
    print("  biggest offenders    : %s"
          % ", ".join("%s%d +%.0f" % (a, n, d) for d, a, n in worst[:6]))
    ep_new = sum(d for d, a, n in worst if n in epitope)
    print("  on the epitope itself: %.0f A^2 (want this near zero - the epitope should look"
          " the same cropped as it does whole)" % ep_new)

    m = dict(site=SITE, pdb="11QC", chain=ref, radius=RADIUS,
             epitope=epitope, kept=keep, n_kept=len(keep),
             new_surface_A2=round(new_total, 1),
             new_hydrophobic_A2=round(new_phobic, 1),
             new_surface_on_epitope_A2=round(ep_new, 1),
             sequence="".join(THREE2ONE[r.get_resname()] for r in cr[0][ref]
                              if r.id[0] == " "))
    json.dump(m, open(os.path.join(D, "TfR1_%s_map.json" % SITE), "w"), indent=0)
    print("\nwrote %s" % out)
    print("wrote %s" % os.path.join(D, "TfR1_%s_map.json" % SITE))


if __name__ == "__main__":
    main()
