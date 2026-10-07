"""TNF-alpha Challenge 2, Phase 1: where to bind.

Objectives 2 (mouse cross-reactivity) and 3 (acid-OFF pH switch) are decided by the
choice of epitope before any backbone exists, exactly as in Challenge 1. Two numbers
decide it here:
  - human/mouse identity, locally at the candidate patch, not globally. The challenge
    page states the soluble domains are only 79% identical - nine points lower than
    EGFR's ectodomain - so conservation is a tighter constraint this time.
  - cationic residues (Arg/Lys) at the interface. The acid-OFF switch needs a
    histidine that is neutral at pH 7.4 and protonated at pH 6.0, placed where that
    new positive charge collides with an existing one.
"""
import json, pathlib, sys
import numpy as np
from Bio import Align
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import protein_letters_3to1

D = pathlib.Path("tnf/data")
OUT = pathlib.Path("tnf/results"); OUT.mkdir(parents=True, exist_ok=True)
MAXS = {"A":129,"R":274,"N":195,"D":193,"C":167,"E":223,"Q":225,"G":104,"H":224,
        "I":197,"L":201,"K":236,"M":224,"F":240,"P":159,"S":155,"T":172,"W":285,
        "Y":263,"V":174}


def fasta(p):
    lines = pathlib.Path(p).read_text().split("\n")
    return "".join(l.strip() for l in lines[1:] if l and not l.startswith(">"))


def align(a, b):
    al = Align.PairwiseAligner()
    al.open_gap_score, al.extend_gap_score = -11, -1
    al.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
    return al.align(a, b)[0]


def chain_seq(chain):
    out = []
    for r in chain:
        aa = protein_letters_3to1.get(r.get_resname())
        if r.id[0] == " " and aa:
            out.append((r.id[1], aa, r))
    return out


def contacts(struct_path, target_chains, partner_chains, cutoff=5.0):
    """Target residues with any atom within `cutoff` of any partner chain atom."""
    mdl = PDBParser(QUIET=True).get_structure("x", struct_path)[0]
    present = {c.id for c in mdl}
    tc = [c for c in target_chains if c in present]
    pc = [c for c in partner_chains if c in present]
    if not tc or not pc:
        return None
    pa = np.array([a.get_coord() for c in pc for r in mdl[c] for a in r])
    hits = {}
    for c in tc:
        for num, aa, r in chain_seq(mdl[c]):
            ra = np.array([a.get_coord() for a in r])
            d = np.min(np.linalg.norm(ra[:, None, :] - pa[None, :, :], axis=-1))
            if d <= cutoff:
                hits.setdefault(num, (aa, round(float(d), 2)))
    return hits


human = fasta(D / "P01375.fasta")
mouse = fasta(D / "P06804.fasta")
sol_h = human[76:233]                      # UniProt 77-233, the assayed construct

# ---- 1. human/mouse identity across the soluble domain -------------------------
aln = align(human, mouse)
hmap = {}                                   # human UniProt index (1-based) -> mouse aa
hi = mi = 0
for (hs, he), (ms, me) in zip(aln.aligned[0], aln.aligned[1]):
    for k in range(he - hs):
        hmap[hs + k + 1] = mouse[ms + k]
ident = [(i, human[i - 1], hmap.get(i)) for i in range(77, 234)]
same = sum(1 for _, a, b in ident if a == b)
print("soluble domain 77-233: %d/%d identical = %.1f%%  (page states 79%%)"
      % (same, len(ident), 100 * same / len(ident)))

# ---- 2. map onto 1TNF and get exposure -----------------------------------------
mdl = PDBParser(QUIET=True).get_structure("t", D / "1TNF.pdb")[0]
chains = [c.id for c in mdl]
print("1TNF chains:", chains)
seqA = chain_seq(mdl["A"])
print("chain A: %d residues, numbering %d-%d" % (len(seqA), seqA[0][0], seqA[-1][0]))

# PDB numbering of 1TNF is mature TNF (1-157) = UniProt - 76. Verify.
ok = sum(1 for num, aa, _ in seqA if 1 <= num <= 157 and sol_h[num - 1] == aa)
print("PDB->UniProt offset +76 matches %d/%d resolved residues" % (ok, len(seqA)))

ShrakeRupley().compute(mdl, level="R")      # trimer: exposure in the assembled form
rsasa = {}
for num, aa, r in chain_seq(mdl["A"]):
    rsasa[num] = r.sasa / MAXS.get(aa, 200)

# ---- 3. functional epitopes -----------------------------------------------------
sets = {}
for label, path, tgt, partner in [
        ("TNFR2 (3ALQ)", D / "3ALQ.pdb", ["A", "B", "C"], ["D", "E", "F", "R"]),
        ("adalimumab (3WD5)", D / "3WD5.pdb", ["A", "B", "C"], ["H", "L"]),
        ("infliximab (4G3Y)", D / "4G3Y.pdb", ["A", "B", "C"], ["H", "L"]),
]:
    h = contacts(path, tgt, partner)
    if h is None:
        print("%-20s chains not found as guessed" % label); continue
    sets[label] = h
    print("%-20s %d contact residues" % (label, len(h)))

json.dump({k: {str(n): v for n, v in s.items()} for k, s in sets.items()},
          open(OUT / "epitopes_raw.json", "w"), indent=1)
json.dump({"identity": [[i, a, b] for i, a, b in ident],
           "rsasa": {str(k): round(v, 3) for k, v in rsasa.items()}},
          open(OUT / "conservation.json", "w"), indent=1)
print("\nwrote tnf/results/epitopes_raw.json and conservation.json")
