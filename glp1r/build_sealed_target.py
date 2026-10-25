"""Rebuild the sealed GLP-1R target locally, identically to the Kaggle cell.

cellGLP4_sealed_setup.py writes this file inside a Kaggle session, so it has never
existed on disk here. RFD4-Proteina runs on hosted hardware and needs the target
uploaded from this machine, so the same crop has to be reproducible locally.

Same constants, same order of operations, same assertions as the notebook cell. The
printed residue count and depth-band table must match the Kaggle log line for line -
if they do not, the two targets are not the same molecule and nothing downstream is
comparable.
"""
import json
import os

import numpy as np
from Bio.PDB import PDBIO, PDBParser, Select
from Bio.PDB.Polypeptide import protein_letters_3to1 as T31

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "glp1r", "6X18.pdb")
OUTDIR = os.path.join(ROOT, "data", "glp1r")

RADIUS, MAXGAP, ECD_END = 10.0, 4, 145
SEAL_R, SEAL_LO, SEAL_HI = 18.0, 28.0, 62.0

mdl = PDBParser(QUIET=True).get_structure("x", RAW)[0]
rec = [r for r in mdl["R"] if r.id[0] == " "]
pep = [r for r in mdl["P"] if r.id[0] == " "]
native_seq = "".join(T31.get(r.get_resname(), "X") for r in pep)
if native_seq[0] != "H" or pep[0].id[1] != 7:
    raise SystemExit("chain P does not start at His7 - the depth reference would be wrong")

patoms = np.array([a.get_coord() for r in pep for a in r])
present = {r.id[1] for r in rec}
mind, ratoms = {}, {}
for r in rec:
    ra = np.array([a.get_coord() for a in r])
    ratoms[r.id[1]] = ra
    mind[r.id[1]] = float(np.min(np.linalg.norm(ra[:, None, :] - patoms[None, :, :], axis=-1)))

base = sorted(k for k, v in mind.items() if v <= RADIUS)
filled = set(base)
for a, b in zip(base, base[1:]):
    if 1 < b - a <= MAXGAP + 1 and all(x in present for x in range(a + 1, b)):
        filled.update(range(a + 1, b))


def ca(r):
    return r["CA"].get_coord() if "CA" in r else None


prelim = [r for r in rec if r.id[1] in filled]
ecd_c = np.mean([ca(r) for r in prelim if r.id[1] <= ECD_END and ca(r) is not None], axis=0)
tm_c = np.mean([ca(r) for r in prelim if r.id[1] > ECD_END and ca(r) is not None], axis=0)
AXIS = (tm_c - ecd_c) / np.linalg.norm(tm_c - ecd_c)
OFF = float(ecd_c @ AXIS)


def depth(x):
    return float(np.dot(x - ecd_c, AXIS))


def radial(x):
    v = x - ecd_c
    return float(np.linalg.norm(v - np.dot(v, AXIS) * AXIS))


for r in rec:
    c = ca(r)
    if c is not None and SEAL_LO <= depth(c) <= SEAL_HI and radial(c) <= SEAL_R:
        filled.add(r.id[1])
ks = sorted(filled)
for a, b in zip(ks, ks[1:]):
    if 1 < b - a <= MAXGAP + 1 and all(x in present for x in range(a + 1, b)):
        filled.update(range(a + 1, b))

contacts = set()
for r in pep:
    ra = np.array([a.get_coord() for a in r])
    for rid, rra in ratoms.items():
        if np.min(np.linalg.norm(ra[:, None, :] - rra[None, :, :], axis=-1)) <= 4.5:
            contacts.add(rid)
lost = sorted(contacts - filled)
if lost:
    raise SystemExit("sealed crop drops GLP-1 contacts %s" % lost)


class Keep(Select):
    def accept_chain(self, c):
        return c.id == "R"

    def accept_residue(self, r):
        return r.id[0] == " " and r.id[1] in filled

    def accept_atom(self, a):
        return a.element != "H"


TARGET = os.path.join(OUTDIR, "GLP1R_sealed.pdb")
pdbio = PDBIO()
pdbio.set_structure(mdl)
pdbio.save(TARGET, Keep())
kept_res = [r for r in rec if r.id[1] in filled]
print("sealed target: %d residues, all %d GLP-1 contacts retained"
      % (len(kept_res), len(contacts)))


def band(residues, lo, hi):
    xyz = np.array([a.get_coord() for r in residues for a in r])
    d = xyz @ AXIS - OFF
    return int(((d >= lo) & (d < hi)).sum())


print("\ndepth band    sealed target    intact receptor    kept")
for lo, hi in ((33.2, 40), (40, 45), (45, 50), (50, 60)):
    nc, nf = band(kept_res, lo, hi), band(rec, lo, hi)
    print("%5.1f-%4.1f A %13d %18d %7.0f%%"
          % (lo, hi, nc, nf, 100.0 * nc / nf if nf else 0.0))

TM_C = sorted(x for x in contacts if x > ECD_END)
ECD_C_ = sorted(x for x in contacts if x <= ECD_END)
native_depth = {r.id[1]: round(depth(ca(r)), 2) for r in pep if ca(r) is not None}
cmap = [dict(uniprot=r.id[1], aa=T31.get(r.get_resname(), "X"),
             depth=round(depth(ca(r)), 2) if ca(r) is not None else None,
             region="ecd" if r.id[1] <= ECD_END else "tm",
             contact=r.id[1] in contacts)
        for r in kept_res]
json.dump(dict(order=cmap, ecd_end=ECD_END, radius=RADIUS, seal_r=SEAL_R,
               ecd_centroid=[float(x) for x in ecd_c], axis=[float(x) for x in AXIS],
               native_glp1=dict(seq=native_seq, start=pep[0].id[1], depth=native_depth),
               tm_contacts=TM_C, ecd_contacts=ECD_C_),
          open(os.path.join(OUTDIR, "glp1r_sealed_map.json"), "w"), indent=1)

print("\nhotspots: %d activation (channel), %d binding (groove), %d combined"
      % (len(TM_C), len(ECD_C_), len(contacts)))
print("\nchannel : %s" % ",".join(str(x) for x in TM_C))
print("groove  : %s" % ",".join(str(x) for x in ECD_C_))
print("\nwrote %s" % TARGET)
