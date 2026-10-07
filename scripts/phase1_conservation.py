"""Phase 1: human/mouse EGFR conservation map and epitope shortlist.

Aligns human (P00533) and mouse (Q01279) EGFR, maps per-residue identity onto
the human ECD structure (PDB 6ARU chain A), computes solvent accessibility, and
ranks candidate surface patches by cross-species conservation and carboxylate
density (Asp/Glu available for histidine pairing in an acid-ON pH switch).
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
HUMAN = os.path.join(ROOT, "data", "human_EGFR_P00533.fasta")
MOUSE = os.path.join(ROOT, "data", "mouse_Egfr_Q01279.fasta")
PDB = os.path.join(ROOT, "data", "structures", "6ARU.pdb")
OUT = os.path.join(ROOT, "results")

# Maximum SASA per residue (Tien et al. 2013, theoretical) for relative SASA.
MAX_SASA = {
    "A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "E": 223, "Q": 225,
    "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240,
    "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174,
}

# UniProt P00533 feature boundaries (mature numbering = UniProt - 24).
DOMAINS = [
    ("I  (L1)", 25, 188),
    ("II (CR1/dimerization arm)", 189, 334),
    ("III (L2, cetuximab epitope)", 335, 505),
    ("IV (CR2)", 506, 645),
]


def load_seq(path):
    return str(next(SeqIO.parse(path, "fasta")).seq)


def align_pair(a, b):
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "global"
    return aligner.align(a, b)[0]


def identity_by_human_index(alignment, human, mouse):
    """Return dict: 1-based human UniProt position -> (human_aa, mouse_aa, identical)."""
    out = {}
    for (h_start, h_end), (m_start, m_end) in zip(*alignment.aligned):
        for off in range(h_end - h_start):
            hi, mi = h_start + off, m_start + off
            out[hi + 1] = (human[hi], mouse[mi], human[hi] == mouse[mi])
    # Positions in gaps are unaligned -> treat as non-identical, mouse unknown.
    for i, aa in enumerate(human, start=1):
        out.setdefault(i, (aa, "-", False))
    return out


def chain_residues(structure, chain_id="A"):
    chain = structure[0][chain_id]
    res = []
    for r in chain:
        if r.id[0] != " ":
            continue
        try:
            aa = protein_letters_3to1[r.get_resname()]
        except KeyError:
            continue
        res.append((r.id[1], aa, r))
    return res


def find_offset(pdb_res, human_seq):
    """Empirically find offset so that pdb_number + offset == UniProt number."""
    best, best_score = None, -1
    numbers = [n for n, _, _ in pdb_res]
    aas = {n: a for n, a, _ in pdb_res}
    for offset in range(-5, 60):
        score = sum(
            1
            for n in numbers
            if 1 <= n + offset <= len(human_seq) and human_seq[n + offset - 1] == aas[n]
        )
        if score > best_score:
            best, best_score = offset, score
    return best, best_score / len(numbers)


def domain_of(uniprot_pos):
    for name, lo, hi in DOMAINS:
        if lo <= uniprot_pos <= hi:
            return name
    return "other"


def main():
    human, mouse = load_seq(HUMAN), load_seq(MOUSE)
    aln = align_pair(human, mouse)
    ident_map = identity_by_human_index(aln, human, mouse)

    ecd = [p for p in range(25, 646)]
    ecd_ident = sum(1 for p in ecd if ident_map[p][2])
    print(f"Human/mouse EGFR full-length identity: "
          f"{sum(1 for p in range(1, 1211) if ident_map[p][2]) / 1210:.1%}")
    print(f"Extracellular region (25-645) identity: {ecd_ident / len(ecd):.1%}")
    for name, lo, hi in DOMAINS:
        n = sum(1 for p in range(lo, hi + 1) if ident_map[p][2])
        print(f"  Domain {name:<30} {n / (hi - lo + 1):6.1%}  ({hi - lo + 1} aa)")

    structure = PDBParser(QUIET=True).get_structure("egfr", PDB)
    ShrakeRupley().compute(structure[0], level="R")
    pdb_res = chain_residues(structure, "A")
    offset, frac = find_offset(pdb_res, human)
    print(f"\n6ARU chain A: {len(pdb_res)} residues resolved; "
          f"PDB numbering + {offset} = UniProt numbering (match {frac:.1%})")

    rows = []
    for num, aa, res in pdb_res:
        up = num + offset
        if up not in ident_map:
            continue
        h_aa, m_aa, same = ident_map[up]
        rel = res.sasa / MAX_SASA.get(aa, 200)
        ca = res["CA"].get_coord() if "CA" in res else None
        if ca is None:
            continue
        rows.append(dict(pdb=num, uniprot=up, aa=aa, mouse_aa=m_aa, identical=same,
                         sasa=round(res.sasa, 1), rel_sasa=round(rel, 3),
                         domain=domain_of(up), xyz=ca.tolist()))

    coords = np.array([r["xyz"] for r in rows])
    exposed = np.array([r["rel_sasa"] >= 0.20 for r in rows])
    identical = np.array([r["identical"] for r in rows])
    acidic = np.array([r["aa"] in ("D", "E") for r in rows])

    # Patch = all residues with CA within 11 A of a surface-exposed center.
    patches = []
    for i, row in enumerate(rows):
        if not exposed[i]:
            continue
        d = np.linalg.norm(coords - coords[i], axis=1)
        nb = (d <= 11.0) & exposed
        n = int(nb.sum())
        if n < 8:
            continue
        cons = float(identical[nb].mean())
        n_acid = int((acidic & nb).sum())
        n_acid_cons = int((acidic & nb & identical).sum())
        patches.append(dict(
            center_uniprot=row["uniprot"], center_pdb=row["pdb"], center_aa=row["aa"],
            domain=row["domain"], n_residues=n, conservation=round(cons, 3),
            n_acidic=n_acid, n_acidic_conserved=n_acid_cons,
            mismatches=[f"{rows[j]['aa']}{rows[j]['uniprot']}{rows[j]['mouse_aa']}"
                        for j in np.where(nb & ~identical)[0]],
            members=[int(rows[j]["uniprot"]) for j in np.where(nb)[0]],
        ))

    # Score: fully conserved patches with carboxylates for His pairing.
    for p in patches:
        p["score"] = round(p["conservation"] ** 3 * (1 + p["n_acidic_conserved"]), 3)
    patches.sort(key=lambda p: (-p["score"], -p["n_acidic_conserved"]))

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "residues.json"), "w") as fh:
        json.dump(rows, fh, indent=1)
    with open(os.path.join(OUT, "patches.json"), "w") as fh:
        json.dump(patches, fh, indent=1)

    print(f"\n{len(patches)} candidate surface patches evaluated.")
    print("\nTop 20 patches by conservation x conserved-carboxylate content:")
    print(f"{'ctr(UP)':>8} {'aa':>3} {'dom':<30} {'n':>3} {'cons':>6} "
          f"{'acid':>5} {'acidC':>6} {'score':>6}")
    for p in patches[:20]:
        print(f"{p['center_uniprot']:>8} {p['center_aa']:>3} {p['domain']:<30} "
              f"{p['n_residues']:>3} {p['conservation']:>6.1%} {p['n_acidic']:>5} "
              f"{p['n_acidic_conserved']:>6} {p['score']:>6.2f}")

    dom3 = [p for p in patches if p["domain"].startswith("III")]
    print(f"\nTop 10 within domain III (recommended functional epitope), "
          f"{len(dom3)} patches:")
    for p in dom3[:10]:
        print(f"{p['center_uniprot']:>8} {p['center_aa']:>3} n={p['n_residues']:>3} "
              f"cons={p['conservation']:>6.1%} acidC={p['n_acidic_conserved']:>2} "
              f"mismatch={','.join(p['mismatches'][:6]) or 'none'}")


if __name__ == "__main__":
    main()
