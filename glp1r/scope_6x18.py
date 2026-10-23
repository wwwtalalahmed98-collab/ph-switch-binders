# Scope the GLP1R agonist problem from the active-state structure 6X18.
#
# The question this answers is not "can we bind the receptor" - we know a soluble
# extracellular domain is bindable, that is the same problem as challenges 1 and 2. The
# question is whether the ACTIVATION half is reachable at all, because the winner has to
# be an agonist and extracellular-domain binding alone gives an antagonist
# (exendin 9-39 is exactly that drug).
#
# So: measure how GLP-1 itself does it. Where does each part of the hormone touch the
# receptor, and how deep does its N-terminus sit? That depth is the bar a designed
# agonist has to clear, and it becomes the geometric proxy we select on - declared up
# front and labelled unvalidated, rather than discovered after submitting, which is how
# the TNF switch rule went wrong.
import sys
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1 as T31

PDB = sys.argv[1]
mdl = PDBParser(QUIET=True).get_structure("x", PDB)[0]
rec = [r for r in mdl["R"] if r.id[0] == " "]
pep = [r for r in mdl["P"] if r.id[0] == " "]

# The extracellular domain of GLP1R is roughly UniProt 24-145; the first transmembrane
# helix starts around 146. Split on that rather than guessing from coordinates.
ECD_END = 145
ecd = [r for r in rec if r.id[1] <= ECD_END]
tm = [r for r in rec if r.id[1] > ECD_END]
print("receptor chain R: %d residues, %d-%d" % (len(rec), rec[0].id[1], rec[-1].id[1]))
print("  extracellular domain (<=%d): %d residues" % (ECD_END, len(ecd)))
print("  transmembrane and beyond   : %d residues" % len(tm))
print("GLP-1 chain P: %d residues, %d-%d, %s"
      % (len(pep), pep[0].id[1], pep[-1].id[1],
         "".join(T31.get(r.get_resname(), "X") for r in pep)))

def ca(r):
    return r["CA"].get_coord() if "CA" in r else None

# Depth axis: from the extracellular domain centroid into the transmembrane bundle.
# Everything is projected onto it, so "depth" is one number per residue, in angstroms,
# increasing as you go further into the membrane.
ecd_c = np.mean([ca(r) for r in ecd if ca(r) is not None], axis=0)
tm_c = np.mean([ca(r) for r in tm if ca(r) is not None], axis=0)
axis = tm_c - ecd_c
axis = axis / np.linalg.norm(axis)
def depth(xyz):
    return float(np.dot(xyz - ecd_c, axis))

print("\n==== GLP-1 DEPTH PROFILE (0 = extracellular domain centroid) ====")
print("%4s %3s %8s  %s" % ("resi", "aa", "depth A", "receptor residues within 4.5 A"))
rec_atoms = {}
for r in rec:
    rec_atoms[r.id[1]] = np.array([a.get_coord() for a in r])

prof = []
for r in pep:
    c = ca(r)
    if c is None:
        continue
    d = depth(c)
    pa = np.array([a.get_coord() for a in r])
    near = []
    for rid, ra in rec_atoms.items():
        if np.min(np.linalg.norm(pa[:, None, :] - ra[None, :, :], axis=-1)) <= 4.5:
            near.append(rid)
    near.sort()
    prof.append((r.id[1], T31.get(r.get_resname(), "X"), d, near))
    tag = "ECD" if all(x <= ECD_END for x in near) and near else (
          "TM " if near and all(x > ECD_END for x in near) else
          "both" if near else "-")
    print("%4d %3s %8.1f  %-4s %s" % (r.id[1], T31.get(r.get_resname(), "X"), d, tag,
                                      ",".join(str(x) for x in near[:12])))

deepest = max(prof, key=lambda p: p[2])
shallow = min(prof, key=lambda p: p[2])
print("\ndeepest  GLP-1 residue: %s%d at %.1f A" % (deepest[1], deepest[0], deepest[2]))
print("shallowest             : %s%d at %.1f A" % (shallow[1], shallow[0], shallow[2]))
print("span the hormone covers: %.1f A" % (deepest[2] - shallow[2]))

# How the contacts divide. This is the number that matters: the share of the hormone's
# grip that comes from the extracellular domain versus the transmembrane bundle.
ecd_contacts = sorted({x for _, _, _, n in prof for x in n if x <= ECD_END})
tm_contacts = sorted({x for _, _, _, n in prof for x in n if x > ECD_END})
print("\nreceptor residues contacted by GLP-1:")
print("  extracellular domain: %d  %s" % (len(ecd_contacts), ecd_contacts))
print("  transmembrane bundle: %d  %s" % (len(tm_contacts), tm_contacts))

n_term = [p for p in prof if p[0] <= 15]
print("\nthe activation half (GLP-1 residues 7-15):")
for resi, aa, d, near in n_term:
    tmn = [x for x in near if x > ECD_END]
    print("  %s%-3d depth %6.1f A, %d transmembrane contacts" % (aa, resi, d, len(tmn)))
print("\nIf a designed agonist cannot put residues near %.0f A depth with contacts in the"
      % max(p[2] for p in n_term))
print("transmembrane bundle, it is an extracellular-domain binder and will not activate.")
