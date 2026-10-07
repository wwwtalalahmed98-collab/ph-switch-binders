"""Calibrate the acid-OFF switch geometry from real structures.

Challenge 1 asked where nature PUTS a histidine next to a carboxylate, because the
protonated imidazolium forms a salt bridge there and binding switches on in acid.

Challenge 2 needs the opposite, so the question inverts: where does nature AVOID
putting a histidine next to an arginine or lysine? A histidine packed hard against a
cation is fine while it is neutral at pH 7.4 and costly once protonated at pH 6.0.
The distances that are rare in real structures are exactly the distances at which
that cost bites - so the rarity itself is the measurement.
"""
import glob, math, statistics as st
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1

HIS_N = ("ND1", "NE2")
ACID_O = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}
CAT_N = {"ARG": ("NE", "NH1", "NH2"), "LYS": ("NZ",)}
P = PDBParser(QUIET=True)


def sidechain_pairs(paths):
    """Closest sidechain nitrogen-to-partner distances for every His in the set."""
    acid, cat, cb_acid, cb_cat = [], [], [], []
    for p in paths:
        try:
            mdl = P.get_structure("x", p)[0]
        except Exception:
            continue
        res = [r for ch in mdl for r in ch if r.id[0] == " "]
        his = [r for r in res if r.get_resname() == "HIS"]
        acids = [r for r in res if r.get_resname() in ACID_O]
        cats = [r for r in res if r.get_resname() in CAT_N]
        for h in his:
            hn = [h[a].get_coord() for a in HIS_N if a in h]
            if not hn:
                continue
            cb = h["CB"].get_coord() if "CB" in h else None
            for partners, names, nlist, cblist in (
                    (acids, ACID_O, acid, cb_acid), (cats, CAT_N, cat, cb_cat)):
                best, best_cb = 1e9, None
                for q in partners:
                    if q is h:
                        continue
                    qa = [q[a].get_coord() for a in names[q.get_resname()] if a in q]
                    if not qa:
                        continue
                    d = min(np.linalg.norm(x - y) for x in hn for y in qa)
                    if d < best:
                        best = d
                        if cb is not None:
                            best_cb = min(np.linalg.norm(cb - y) for y in qa)
                if best < 12.0:
                    nlist.append(float(best))
                    if best_cb is not None:
                        cblist.append(float(best_cb))
    return acid, cat, cb_acid, cb_cat


paths = sorted(glob.glob("data/calib/*.pdb") + glob.glob("data/structures/*.pdb")
               + glob.glob("tnf/data/*.pdb"))
print("calibration set: %d structures" % len(paths))
acid, cat, cb_acid, cb_cat = sidechain_pairs(paths)
print("His-carboxylate closest-approach samples: %d" % len(acid))
print("His-cation      closest-approach samples: %d" % len(cat))


def pct(v, q):
    v = sorted(v)
    return v[max(0, min(len(v) - 1, int(q * len(v))))]


print("\nclosest sidechain N->partner distance (A)")
print("%-18s %6s %6s %6s %6s %6s" % ("", "p1", "p5", "p25", "median", "n<4.0A"))
for label, v in (("His-carboxylate", acid), ("His-cation", cat)):
    print("%-18s %6.2f %6.2f %6.2f %6.2f %6d"
          % (label, pct(v, .01), pct(v, .05), pct(v, .25), st.median(v),
             sum(1 for x in v if x < 4.0)))

close_a = 100 * sum(1 for x in acid if x < 4.0) / len(acid)
close_c = 100 * sum(1 for x in cat if x < 4.0) / len(cat)
print("\nwithin 4.0 A: %.1f%% of His-carboxylate pairs, %.1f%% of His-cation pairs"
      % (close_a, close_c))
print("ratio: nature puts a His that close to an acid %.1fx more often than to a cation"
      % (close_a / close_c if close_c else float("inf")))

print("\nCbeta(His) -> partner distance when the sidechains are in contact")
for label, v in (("to carboxylate", cb_acid), ("to cation", cb_cat)):
    sel = [x for x in v if x < 9.0]
    print("  %-16s p5 %.2f  median %.2f  p95 %.2f  (n=%d)"
          % (label, pct(sel, .05), st.median(sel), pct(sel, .95), len(sel)))
