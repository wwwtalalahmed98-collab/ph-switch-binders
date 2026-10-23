# How big a target do we have to carry? AfDesign cost scales steeply with target
# length - the TNF challenge used a 122-residue crop at about 43 s per iteration, and
# the whole GLP1R receptor in 6X18 is 384 residues, which would be far too slow.
#
# So: crop to a shell around the GLP-1 peptide itself, which is the binding channel we
# care about, and see how the residue count falls with radius. This decides whether the
# agonist route is affordable at all.
import sys
import numpy as np
from Bio.PDB import PDBParser

mdl = PDBParser(QUIET=True).get_structure("x", sys.argv[1])[0]
rec = [r for r in mdl["R"] if r.id[0] == " "]
pep = [r for r in mdl["P"] if r.id[0] == " "]
patoms = np.array([a.get_coord() for r in pep for a in r])
ECD_END = 145

mind = {}
for r in rec:
    ra = np.array([a.get_coord() for a in r])
    mind[r.id[1]] = float(np.min(np.linalg.norm(ra[:, None, :] - patoms[None, :, :], axis=-1)))

print("%6s %6s %6s %6s  %s" % ("radius", "total", "ECD", "TM", "longest contiguous runs"))
for rad in (8, 10, 12, 14, 16, 18, 20):
    keep = sorted(k for k, v in mind.items() if v <= rad)
    ecd = [k for k in keep if k <= ECD_END]
    tm = [k for k in keep if k > ECD_END]
    runs, start, prev = [], None, None
    for k in keep:
        if start is None:
            start = prev = k; continue
        if k == prev + 1:
            prev = k; continue
        runs.append((start, prev)); start = prev = k
    if start is not None:
        runs.append((start, prev))
    runs.sort(key=lambda ab: ab[1] - ab[0], reverse=True)
    top = " ".join("%d-%d" % ab for ab in runs[:4])
    print("%6d %6d %6d %6d  %d segments, longest: %s"
          % (rad, len(keep), len(ecd), len(tm), len(runs), top))

print("\nreference: the TNF challenge target was 122 residues at ~43 s per iteration;")
print("cost rises steeply with length, so a 350-residue target is not affordable.")
