"""Final epitope funnel for TfR1: conservation, glycans and the transferrin site together.

The earlier scripts applied these filters one at a time, and each pass overturned the one
before. Ranking by sequence identity promoted a patch whose centre residue loses a charge in
the mouse; the structural check fixed that and promoted the helical-domain face; the glycan
check then found N727 sitting inside that face. Each check was right and each was blind to
the next one, which is why the answer kept moving.

This applies all of them to every candidate patch at once, so the shortlist is stable rather
than provisional.

FILTERS, in the order a design actually has to survive them:
  1. exposed        - mean relative SASA over the patch, on the assembled dimer
  2. transferrin    - no patch residue in, and the centroid clear of, the 1SUV/1DE4 footprint
  3. glycan         - centroid beyond REACH of any N-glycosylation site in either species
  4. conservation   - fraction of patch residues identical human/mouse

Distances are advisory, not physics: no glycan geometry is modelled and the mouse structure
is a prediction. What this gives is a ranked shortlist with the reason each candidate failed,
not a guarantee.
"""
import io
import json
import os
import re

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "data", "tfr1")
RADIUS = 14.0        # patch radius, about a miniprotein footprint
MIN_MEMBERS = 12
MIN_SASA = 0.20
TF_CLEAR = 12.0      # centroid must be at least this far from the transferrin footprint
GLYCAN_REACH = 20.0
ECTO = (121, 760)


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


def sequons(seq):
    return {m.start() + 1 for m in re.finditer(r"(?=N[^P][ST])", seq)}


def main():
    cons = json.load(open(os.path.join(D, "tfr1_conservation.json")))
    res = {r["uniprot"]: r for r in cons["residues"]}
    pat = json.load(open(os.path.join(D, "tfr1_patches.json")))
    tf = [n for n in pat["transferrin_site"] if n in res]
    tfx = np.array([res[n]["xyz"] for n in tf])

    human, mouse = fasta(os.path.join(D, "P02786.fasta")), fasta(os.path.join(D, "Q62351.fasta"))
    h2m = hm_map(human, mouse)
    m2h = {v: k for k, v in h2m.items()}
    hs = {int(n) for n in sequons(human) if ECTO[0] <= n <= ECTO[1]}
    ms = {int(m2h[n]) for n in sequons(mouse)
          if n in m2h and ECTO[0] <= m2h[n] <= ECTO[1]}
    gly = sorted(int(n) for n in (hs | ms) if n in res)
    glyx = np.array([res[n]["xyz"] for n in gly])
    print("glycosylation sites in either species, resolved: %s"
          % ", ".join("N%d%s" % (n, "" if n in hs and n in ms else
                                 ("[h]" if n in hs else "[m]")) for n in gly))
    print("transferrin/HFE footprint: %d residues\n" % len(tf))

    exp = [r for r in res.values() if r["rel_sasa"] >= MIN_SASA]
    xyz = np.array([r["xyz"] for r in exp])
    anchor = np.array(res[min(res)]["xyz"])

    rows, rejected = [], {"small": 0, "transferrin": 0, "glycan": 0}
    for i, centre in enumerate(exp):
        d = np.linalg.norm(xyz - xyz[i], axis=1)
        mem = [exp[j] for j in np.where(d <= RADIUS)[0]]
        if len(mem) < MIN_MEMBERS:
            rejected["small"] += 1
            continue
        nums = [m["uniprot"] for m in mem]
        c = np.mean([m["xyz"] for m in mem], axis=0)
        dtf = float(np.linalg.norm(tfx - c, axis=1).min())
        # a patch that CONTAINS footprint residues is out regardless of its centroid
        overlap = sorted(set(nums) & set(tf))
        if overlap or dtf < TF_CLEAR:
            rejected["transferrin"] += 1
            continue
        dg = float(np.linalg.norm(glyx - c, axis=1).min())
        kg = gly[int(np.argmin(np.linalg.norm(glyx - c, axis=1)))]
        if dg < GLYCAN_REACH:
            rejected["glycan"] += 1
            continue
        nc = sum(m["conserved"] for m in mem)
        rows.append(dict(centre=centre["uniprot"], centre_aa=centre["aa"], n=len(mem),
                         conserved=nc, frac=nc / float(len(mem)),
                         mean_sasa=float(np.mean([m["rel_sasa"] for m in mem])),
                         dist_tf=dtf, dist_glycan=dg, nearest_glycan=int(kg),
                         height=float(np.linalg.norm(c - anchor)),
                         members=sorted(nums)))

    print("candidate patches: %d too small, %d hit the transferrin site, %d glycan-shadowed,"
          " %d survive" % (rejected["small"], rejected["transferrin"], rejected["glycan"],
                           len(rows)))
    if not rows:
        raise SystemExit("nothing survives all three filters. Loosen one and say which.")

    rows.sort(key=lambda r: (r["frac"], r["n"]), reverse=True)
    print("\n%-9s %4s %6s %6s %6s %7s %7s  %s"
          % ("centre", "size", "cons", "sasa", "dTf", "dGlycan", "height", "span"))
    for r in rows[:14]:
        print("%-9s %4d %5.0f%% %6.2f %6.0f %7.0f %7.0f  %d-%d"
              % ("%s%d" % (r["centre_aa"], r["centre"]), r["n"], 100 * r["frac"],
                 r["mean_sasa"], r["dist_tf"], r["dist_glycan"], r["height"],
                 min(r["members"]), max(r["members"])))

    # collapse to distinct sites: patches sharing most members are the same place
    sites, seen = [], []
    for r in rows:
        s = set(r["members"])
        if any(len(s & t) / float(min(len(s), len(t))) > 0.5 for t in seen):
            continue
        seen.append(s)
        sites.append(r)
    print("\n%d distinct sites after collapsing overlapping patches:" % len(sites))
    for r in sites[:8]:
        print("  %s%-6d %2d residues, %.0f%% conserved, %.0f A from Tf site,"
              " %.0f A from N%d" % (r["centre_aa"], r["centre"], r["n"], 100 * r["frac"],
                                    r["dist_tf"], r["dist_glycan"], r["nearest_glycan"]))

    json.dump(dict(filters=dict(radius=RADIUS, min_members=MIN_MEMBERS,
                                min_sasa=MIN_SASA, tf_clear=TF_CLEAR,
                                glycan_reach=GLYCAN_REACH),
                   glycosylation=gly, sites=sites[:12], all_passing=rows[:40]),
              open(os.path.join(D, "tfr1_epitope_shortlist.json"), "w"), indent=0)
    print("\nwrote data/tfr1/tfr1_epitope_shortlist.json")


if __name__ == "__main__":
    main()
