"""Cross-check our ECD epitope against 3IOL - the structure Adaptyv's own
design-a-protein.com hands people as the GLP-1R target.

Ours came from 6X18 (full-length receptor, cryo-EM, GLP-1 bound). Theirs is
3IOL (isolated ECD, crystal, GLP-1 bound). If the two agree on which receptor
residues the hormone grips, our landing-pad binders are aimed at the epitope
the official on-ramp points at. If they disagree, we are designing against a
conformation the lab assay may not present.
"""
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1 as T31

P = PDBParser(QUIET=True)
CUT, ECD_END = 4.5, 145


def contacts(path, rec_chain, pep_chain, limit=None):
    m = P.get_structure("x", path)[0]
    rec = [r for r in m[rec_chain] if r.id[0] == " "]
    pep = [r for r in m[pep_chain] if r.id[0] == " "]
    pa = np.array([a.get_coord() for r in pep for a in r])
    out = {}
    for r in rec:
        if limit and r.id[1] > limit:
            continue
        ra = np.array([a.get_coord() for a in r])
        d = np.min(np.linalg.norm(ra[:, None, :] - pa[None, :, :], axis=-1))
        if d <= CUT:
            out[r.id[1]] = T31.get(r.get_resname(), "X")
    return out, len(rec), len(pep)


ours, nr6, np6 = contacts("data/glp1r/6X18.pdb", "R", "P", limit=ECD_END)
theirs, nr3, np3 = contacts("data/glp1r/3IOL.pdb", "A", "B")

print("6X18 (ours):   receptor %d res, hormone %d res -> %d ECD contacts"
      % (nr6, np6, len(ours)))
print("3IOL (theirs): receptor %d res, hormone %d res -> %d ECD contacts"
      % (nr3, np3, len(theirs)))

both = sorted(set(ours) & set(theirs))
only_o = sorted(set(ours) - set(theirs))
only_t = sorted(set(theirs) - set(ours))

# identity check: does the residue number mean the same amino acid in both?
mismatch = [(k, ours[k], theirs[k]) for k in both if ours[k] != theirs[k]]

print("\nshared     %3d  %s" % (len(both), ", ".join("%s%d" % (ours[k], k) for k in both)))
print("6X18 only  %3d  %s" % (len(only_o), ", ".join("%s%d" % (ours[k], k) for k in only_o)))
print("3IOL only  %3d  %s" % (len(only_t), ", ".join("%s%d" % (theirs[k], k) for k in only_t)))
print("\nresidue-identity mismatches at shared positions: %d %s"
      % (len(mismatch), mismatch if mismatch else "(numbering agrees)"))

j = len(both) / len(set(ours) | set(theirs))
print("\nJaccard overlap %.2f   |   %d of our %d ECD contacts (%.0f%%) are in 3IOL's set"
      % (j, len(both), len(ours), 100.0 * len(both) / len(ours)))
