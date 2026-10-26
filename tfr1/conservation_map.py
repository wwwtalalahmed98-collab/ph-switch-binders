"""Map human/mouse TfR1 identity onto the ectodomain surface, and rank conserved patches.

WHY THIS EXISTS
---------------
Challenge 4 wants one binder that works on human AND mouse TfR1. The organizers recommend
the apical domain, and also state that the apical patch differs at six of nine positions
while the two ectodomains are 78% identical overall. So the recommended site is the worst
place to meet the cross-species requirement. This script finds the places that are not.

A patch is scored on three things, all of which have to hold at once:
  - conserved: human and mouse residues identical across the patch
  - exposed: actually reachable by a binder on the assembled ectodomain
  - off-limits sites excluded: we do not want the transferrin site

NUMBERING. Everything is reported in UniProt numbering for human TfR1 (P02786). The PDB
author numbering for 11QC is checked against the UniProt sequence residue by residue and
the script refuses to continue if they disagree. Three numbering bugs in this project all
came from two systems silently differing.
"""
import io
import json
import os
import sys

import requests
from Bio import PDB
from Bio.Align import PairwiseAligner
from Bio.PDB.SASA import ShrakeRupley

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "data", "tfr1")
PDB_ID = "11QC"
ECTO = (121, 760)          # the construct Adaptyv supplies, UniProt numbering
PATCH_RADIUS = 10.0        # A, about the footprint of a miniprotein interface
MIN_REL_SASA = 0.20        # a residue counts as exposed above this fraction

THREE2ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
             "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
             "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
             "TYR": "Y", "VAL": "V"}

# Tien et al. 2013, theoretical maximum SASA per residue, used to make exposure a
# fraction rather than a raw area that depends on residue size.
MAXASA = {"A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "Q": 225, "E": 223,
          "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240,
          "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174}


def fasta(path):
    return "".join(l.strip() for l in io.open(path) if not l.startswith(">"))


def fetch_pdb():
    p = os.path.join(D, PDB_ID + ".cif")
    if not os.path.exists(p):
        r = requests.get("https://files.rcsb.org/download/%s.cif" % PDB_ID, timeout=60)
        r.raise_for_status()
        io.open(p, "w", encoding="utf-8", newline="\n").write(r.text)
    return p


def align_species(h, m):
    """Per-position identity, indexed by human UniProt residue number."""
    al = PairwiseAligner()
    al.open_gap_score, al.extend_gap_score = -11, -1
    al.substitution_matrix = __import__(
        "Bio.Align.substitution_matrices", fromlist=["load"]).load("BLOSUM62")
    a = al.align(h, m)[0]
    same = {}
    for (hs, he), (ms, me) in zip(a.aligned[0], a.aligned[1]):
        for k in range(he - hs):
            same[hs + k + 1] = (h[hs + k] == m[ms + k])
    ident = sum(same.values()) / float(len(same))
    return same, ident, len(same)


def main():
    human = fasta(os.path.join(D, "P02786.fasta"))
    mouse = fasta(os.path.join(D, "Q62351.fasta"))
    same, ident, naln = align_species(human, mouse)
    print("human %d aa, mouse %d aa" % (len(human), len(mouse)))
    print("aligned %d positions, %.1f%% identical overall" % (naln, 100 * ident))

    ecto = [n for n in same if ECTO[0] <= n <= ECTO[1]]
    ei = sum(same[n] for n in ecto) / float(len(ecto))
    print("ectodomain %d-%d: %d positions, %.1f%% identical"
          % (ECTO[0], ECTO[1], len(ecto), 100 * ei))
    print("  (the challenge page states 78%% for the ectodomain - if the number above is"
          " far from that, stop and work out why before trusting anything downstream)")

    st = PDB.MMCIFParser(QUIET=True).get_structure(PDB_ID, fetch_pdb())[0]
    chains = [c for c in st if sum(1 for r in c if r.id[0] == " ") > 100]
    print("\n%s: using chains %s" % (PDB_ID, [c.id for c in chains]))

    # --- numbering check against UniProt, before any geometry is trusted ---
    bad, seen = [], 0
    for c in chains:
        for r in c:
            if r.id[0] != " " or r.get_resname() not in THREE2ONE:
                continue
            n = r.id[1]
            if not (1 <= n <= len(human)):
                continue
            seen += 1
            if THREE2ONE[r.get_resname()] != human[n - 1]:
                bad.append((c.id, n, r.get_resname(), human[n - 1]))
    print("numbering check: %d residues compared, %d mismatches" % (seen, len(bad)))
    if bad:
        for b in bad[:10]:
            print("   chain %s %s%d is %s in the PDB, %s in UniProt" % (b[0], b[2], b[1],
                                                                        b[2], b[3]))
        raise SystemExit("PDB author numbering does not match UniProt. Stopping.")

    ShrakeRupley().compute(st, level="R")

    # Keep ONE subunit's geometry. SASA was computed on the whole dimer above, so the
    # subunit-subunit interface is correctly counted as buried, but every coordinate
    # below comes from the same chain. Mixing the two copies puts residues 50-100 A
    # apart into the same "patch" and makes every spatial result nonsense - this was a
    # real bug here, caught because the patch counts came out absurdly low.
    ref = chains[0].id
    res = {}
    for c in chains:
        if c.id != ref:
            continue
        for r in c:
            if r.id[0] != " " or r.get_resname() not in THREE2ONE:
                continue
            n = r.id[1]
            if n not in same:
                continue
            aa = THREE2ONE[r.get_resname()]
            try:
                ca = r["CA"].get_coord()
            except KeyError:
                continue
            res[n] = dict(uniprot=n, aa=aa, chain=c.id,
                          rel_sasa=round(float(r.sasa / MAXASA[aa]), 3),
                          conserved=bool(same[n]),
                          xyz=[round(float(x), 2) for x in ca])
    print("geometry taken from chain %s only; SASA computed on the full dimer" % ref)

    exposed = [r for r in res.values() if r["rel_sasa"] >= MIN_REL_SASA]
    print("resolved %d residues, %d exposed (rel SASA >= %.2f)"
          % (len(res), len(exposed), MIN_REL_SASA))

    out = os.path.join(D, "tfr1_conservation.json")
    json.dump(dict(pdb=PDB_ID, ectodomain=list(ECTO), overall_identity=round(ident, 4),
                   ectodomain_identity=round(ei, 4), residues=list(res.values())),
              open(out, "w"), indent=0)
    print("wrote %s" % out)


if __name__ == "__main__":
    main()
