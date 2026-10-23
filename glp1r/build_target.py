# Build a croppable GLP1R target from the active-state structure 6X18.
#
# Two competing costs. A tight shell around GLP-1 keeps the residue count affordable but
# shatters the pocket into dozens of fragments, which is not a structure any more. A
# loose shell keeps the pocket intact but the target gets long, and AfDesign cost rises
# steeply with target length - the TNF challenge ran a 122-residue target at about 43 s
# per iteration.
#
# Compromise: take a shell, then bridge short gaps so the kept residues form a smaller
# number of real segments. Report the trade so the choice is made on numbers.
import sys
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import protein_letters_3to1 as T31

PDB, OUT = sys.argv[1], sys.argv[2]
RADIUS = float(sys.argv[3]) if len(sys.argv) > 3 else 12.0
MAXGAP = int(sys.argv[4]) if len(sys.argv) > 4 else 4
ECD_END = 145

mdl = PDBParser(QUIET=True).get_structure("x", PDB)[0]
rec = [r for r in mdl["R"] if r.id[0] == " "]
pep = [r for r in mdl["P"] if r.id[0] == " "]
patoms = np.array([a.get_coord() for r in pep for a in r])

mind = {}
for r in rec:
    ra = np.array([a.get_coord() for a in r])
    mind[r.id[1]] = float(np.min(np.linalg.norm(ra[:, None, :] - patoms[None, :, :], axis=-1)))

keep = sorted(k for k, v in mind.items() if v <= RADIUS)
present = {r.id[1] for r in rec}

# bridge gaps of up to MAXGAP residues, but only where the receptor actually has those
# residues modelled - never invent a connection across a disordered stretch
filled = set(keep)
for a, b in zip(keep, keep[1:]):
    if 1 < b - a <= MAXGAP + 1 and all(x in present for x in range(a + 1, b)):
        filled.update(range(a + 1, b))
keep2 = sorted(filled)

def segments(ks):
    out, start, prev = [], None, None
    for k in ks:
        if start is None:
            start = prev = k; continue
        if k == prev + 1:
            prev = k; continue
        out.append((start, prev)); start = prev = k
    if start is not None:
        out.append((start, prev))
    return out

segs = segments(keep2)
ecd = [k for k in keep2 if k <= ECD_END]
tm = [k for k in keep2 if k > ECD_END]
print("radius %.0f A, gaps up to %d bridged" % (RADIUS, MAXGAP))
print("  kept %d residues (%d extracellular, %d transmembrane) in %d segments"
      % (len(keep2), len(ecd), len(tm), len(segs)))
print("  before bridging: %d residues in %d segments" % (len(keep), len(segments(keep))))

# the GLP-1 contacts must survive the crop, or the pocket we are designing into is gone
contacts = set()
for r in pep:
    ra = np.array([a.get_coord() for a in r])
    for rr in rec:
        rra = np.array([a.get_coord() for a in rr])
        if np.min(np.linalg.norm(ra[:, None, :] - rra[None, :, :], axis=-1)) <= 4.5:
            contacts.add(rr.id[1])
lost = sorted(contacts - set(keep2))
print("  GLP-1 contact residues retained: %d of %d%s"
      % (len(contacts) - len(lost), len(contacts),
         "" if not lost else "  LOST: %s" % lost))
if lost:
    raise SystemExit("crop drops contacts the hormone actually uses - widen the radius")

class Keep(Select):
    def accept_chain(self, c): return c.id == "R"
    def accept_residue(self, r): return r.id[0] == " " and r.id[1] in filled
    def accept_atom(self, a): return a.element != "H"

io = PDBIO(); io.set_structure(mdl)
io.save(OUT, Keep())
print("  wrote %s" % OUT)

# hotspots: the transmembrane contacts are what an AGONIST has to reach. Quote them
# separately from the extracellular ones so the two design jobs can be aimed differently.
tm_contacts = sorted(x for x in contacts if x > ECD_END)
ecd_contacts = sorted(x for x in contacts if x <= ECD_END)
print("\nhotspots for the ACTIVATION job (transmembrane, %d):" % len(tm_contacts))
print("  " + ",".join("R%d" % x for x in tm_contacts))
print("hotspots for the BINDING job (extracellular domain, %d):" % len(ecd_contacts))
print("  " + ",".join("R%d" % x for x in ecd_contacts))
seq = "".join(T31.get(r.get_resname(), "X") for r in rec if r.id[1] in filled)
print("\ncropped target length seen by the model: %d" % len(seq))
