"""Phase 1b: cetuximab footprint conservation + ranked epitope shortlist.

Tests the central hypothesis of the plan: that the cetuximab epitope on domain
III (the epitope the challenge page recommends) carries human/mouse
substitutions that explain cetuximab's lack of mouse cross-reactivity, and is
therefore the wrong surface to target for Objective 2.

Then ranks alternative surface patches by cross-species conservation, acidic
(Asp/Glu) content for acid-ON histidine pairing, and functional relevance.
"""
import json
import os

import numpy as np
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import protein_letters_3to1

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results")
OFFSET = 24  # PDB number + 24 = UniProt number (verified 99.7% in phase1_conservation)

MAX_SASA = {
    "A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "E": 223, "Q": 225,
    "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240,
    "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174,
}
DOMAINS = [("I (L1)", 25, 188), ("II (CR1)", 189, 334),
           ("III (L2)", 335, 505), ("IV (CR2)", 506, 645)]
# Domain II dimerization arm (mature 242-259 -> UniProt +24).
DIMER_ARM = range(266, 284)


def load(path):
    return str(next(SeqIO.parse(path, "fasta")).seq)


def identity_map(human, mouse):
    al = Align.PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score, al.mode = -11, -1, "global"
    aln = al.align(human, mouse)[0]
    out = {}
    for (hs, he), (ms, _) in zip(*aln.aligned):
        for off in range(he - hs):
            out[hs + off + 1] = (human[hs + off], mouse[ms + off],
                                 human[hs + off] == mouse[ms + off])
    for i, aa in enumerate(human, 1):
        out.setdefault(i, (aa, "-", False))
    return out


def domain_of(p):
    for n, lo, hi in DOMAINS:
        if lo <= p <= hi:
            return n
    return "other"


def main():
    human = load(os.path.join(ROOT, "data", "human_EGFR_P00533.fasta"))
    mouse = load(os.path.join(ROOT, "data", "mouse_Egfr_Q01279.fasta"))
    ident = identity_map(human, mouse)

    st = PDBParser(QUIET=True).get_structure("x",
                                             os.path.join(ROOT, "data", "structures", "6ARU.pdb"))
    model = st[0]
    # SASA must be computed on EGFR alone: with the Fab present, the cetuximab
    # face is artificially buried and would be excluded from patch analysis.
    ShrakeRupley().compute(model["A"], level="R")

    def residues(chain_id):
        out = []
        for r in model[chain_id]:
            if r.id[0] != " ":
                continue
            try:
                aa = protein_letters_3to1[r.get_resname()]
            except KeyError:
                continue
            out.append((r.id[1], aa, r))
        return out

    egfr = residues("A")
    fab = residues("B") + residues("C")
    fab_atoms = np.array([a.get_coord() for _, _, r in fab for a in r])

    # --- Cetuximab footprint: EGFR residues with any atom within 5 A of Fab ---
    footprint = []
    for num, aa, r in egfr:
        atoms = np.array([a.get_coord() for a in r])
        d = np.linalg.norm(atoms[:, None, :] - fab_atoms[None, :, :], axis=2)
        if d.min() <= 5.0:
            up = num + OFFSET
            h, m, same = ident[up]
            footprint.append(dict(uniprot=up, aa=aa, mouse_aa=m, identical=same,
                                  domain=domain_of(up), min_dist=round(float(d.min()), 2)))

    n_same = sum(1 for f in footprint if f["identical"])
    subs = [f for f in footprint if not f["identical"]]
    print("=" * 74)
    print("CETUXIMAB FOOTPRINT ON HUMAN EGFR (6ARU, any atom within 5 A of Fab)")
    print("=" * 74)
    print(f"Footprint size: {len(footprint)} residues")
    print(f"Identical in mouse: {n_same}/{len(footprint)} = {n_same/len(footprint):.1%}")
    print(f"Substituted in mouse: {len(subs)}")
    for f in subs:
        print(f"   {f['aa']}{f['uniprot']}{f['mouse_aa']}  "
              f"(domain {f['domain']}, {f['min_dist']} A from Fab)")
    print("\nFull footprint:",
          ", ".join(f"{f['aa']}{f['uniprot']}" for f in footprint))

    # --- Surface patches ---
    rows = []
    for num, aa, r in egfr:
        up = num + OFFSET
        if up not in ident or "CA" not in r:
            continue
        h, m, same = ident[up]
        rows.append(dict(uniprot=up, aa=aa, mouse_aa=m, identical=same,
                         rel_sasa=r.sasa / MAX_SASA.get(aa, 200),
                         domain=domain_of(up), xyz=r["CA"].get_coord()))

    xyz = np.array([r["xyz"] for r in rows])
    exposed = np.array([r["rel_sasa"] >= 0.20 for r in rows])
    same = np.array([r["identical"] for r in rows])
    acid = np.array([r["aa"] in ("D", "E") for r in rows])
    fp_set = {f["uniprot"] for f in footprint}

    patches = []
    for i, row in enumerate(rows):
        if not exposed[i]:
            continue
        nb = (np.linalg.norm(xyz - xyz[i], axis=1) <= 13.0) & exposed
        n = int(nb.sum())
        if n < 12:
            continue
        cons = float(same[nb].mean())
        n_acid_c = int((acid & nb & same).sum())
        members = [rows[j]["uniprot"] for j in np.where(nb)[0]]
        patches.append(dict(
            center=row["uniprot"], center_aa=row["aa"], domain=row["domain"],
            n=n, conservation=round(cons, 3), n_acidic_conserved=n_acid_c,
            mean_exposure=round(float(np.mean([rows[j]["rel_sasa"]
                                               for j in np.where(nb)[0]])), 3),
            overlaps_cetuximab=len(set(members) & fp_set),
            in_dimer_arm=len(set(members) & set(DIMER_ARM)),
            mismatches=[f"{rows[j]['aa']}{rows[j]['uniprot']}{rows[j]['mouse_aa']}"
                        for j in np.where(nb & ~same)[0]],
            members=members))

    # Rank: fully conserved first, then carboxylate capacity, then patch size.
    for p in patches:
        p["score"] = round(p["conservation"] ** 6 * (p["n_acidic_conserved"] + 1)
                           * np.log1p(p["n"]), 2)
    patches.sort(key=lambda p: -p["score"])

    os.makedirs(OUT, exist_ok=True)
    json.dump(footprint, open(os.path.join(OUT, "cetuximab_footprint.json"), "w"), indent=1)
    json.dump(patches, open(os.path.join(OUT, "patches_ranked.json"), "w"), indent=1,
              default=float)

    print("\n" + "=" * 74)
    print("TOP CANDIDATE EPITOPE PATCHES (13 A radius, >=12 exposed residues)")
    print("=" * 74)
    print(f"{'center':>7} {'aa':>3} {'domain':<10} {'n':>3} {'cons':>7} {'acidC':>6} "
          f"{'exp':>5} {'cetux':>6} {'arm':>4} {'score':>6}")
    for p in patches[:25]:
        print(f"{p['center']:>7} {p['center_aa']:>3} {p['domain']:<10} {p['n']:>3} "
              f"{p['conservation']:>7.1%} {p['n_acidic_conserved']:>6} "
              f"{p['mean_exposure']:>5.2f} {p['overlaps_cetuximab']:>6} "
              f"{p['in_dimer_arm']:>4} {p['score']:>6.2f}")

    print("\nFully conserved patches (100% human/mouse identity):")
    full = [p for p in patches if p["conservation"] == 1.0]
    print(f"  {len(full)} of {len(patches)} patches")
    by_up = {r["uniprot"]: r for r in rows}
    for p in full[:10]:
        carbox = [f"{by_up[m]['aa']}{m}(rSASA {by_up[m]['rel_sasa']:.2f})"
                  for m in sorted(p["members"])
                  if by_up[m]["aa"] in ("D", "E") and by_up[m]["identical"]]
        print(f"   center {p['center_aa']}{p['center']} ({p['domain']}) n={p['n']} "
              f"cetux_overlap={p['overlaps_cetuximab']} arm={p['in_dimer_arm']}")
        print(f"      residues: {','.join(str(m) for m in sorted(p['members']))}")
        print(f"      conserved carboxylates (His-pairing partners): "
              f"{', '.join(carbox) if carbox else 'none'}")


if __name__ == "__main__":
    main()
