"""Phase 1c: functional-site footprints and final epitope shortlist.

Computes, empirically from structures, (a) the EGF ligand footprint on EGFR
(1IVO, EGF-bound dimerized conformation) and (b) the receptor-receptor
dimerization interface (1IVO chain A vs chain B), then intersects both with the
cross-species conserved patches from phase1_epitopes.py.

An epitope is only therapeutically credible if blocking it plausibly blocks
function, so overlap with the ligand site or the dimer interface is the
functional-relevance criterion the competition asks for.
"""
import json
import os

import numpy as np
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results")
OFFSET_6ARU = 24


def load(p):
    return str(next(SeqIO.parse(p, "fasta")).seq)


def residues(model, chain_id):
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


def find_offset(res, seq):
    best, score = 0, -1
    aas = {n: a for n, a, _ in res}
    for off in range(-5, 60):
        s = sum(1 for n in aas
                if 1 <= n + off <= len(seq) and seq[n + off - 1] == aas[n])
        if s > score:
            best, score = off, s
    return best, score / len(res)


def contacts(res_a, atoms_b, cutoff=5.0):
    out = {}
    for num, aa, r in res_a:
        at = np.array([a.get_coord() for a in r])
        d = np.linalg.norm(at[:, None, :] - atoms_b[None, :, :], axis=2).min()
        if d <= cutoff:
            out[num] = (aa, round(float(d), 2))
    return out


def main():
    human = load(os.path.join(ROOT, "data", "human_EGFR_P00533.fasta"))
    st = PDBParser(QUIET=True).get_structure("x",
                                             os.path.join(ROOT, "data", "structures", "1IVO.pdb"))
    m = st[0]
    egfr_a = residues(m, "A")
    off, frac = find_offset(egfr_a, human)
    print(f"1IVO chain A: {len(egfr_a)} residues, PDB + {off} = UniProt "
          f"(match {frac:.1%})")

    egf_atoms = np.array([a.get_coord() for _, _, r in residues(m, "C") for a in r])
    egfr_b_atoms = np.array([a.get_coord() for _, _, r in residues(m, "B") for a in r])

    lig = contacts(egfr_a, egf_atoms)
    dim = contacts(egfr_a, egfr_b_atoms)
    lig_up = {n + off for n in lig}
    dim_up = {n + off for n in dim}

    print(f"\nEGF ligand footprint on EGFR: {len(lig)} residues")
    print("  " + ", ".join(f"{a}{n + off}" for n, (a, _) in sorted(lig.items())))
    print(f"\nReceptor dimerization interface (chain A vs B): {len(dim)} residues")
    print("  " + ", ".join(f"{a}{n + off}" for n, (a, _) in sorted(dim.items())))

    patches = json.load(open(os.path.join(OUT, "patches_ranked.json")))
    cetux = {f["uniprot"] for f in
             json.load(open(os.path.join(OUT, "cetuximab_footprint.json")))}

    for p in patches:
        mem = set(p["members"])
        p["overlap_ligand"] = len(mem & lig_up)
        p["overlap_dimer"] = len(mem & dim_up)
        p["functional"] = p["overlap_ligand"] + p["overlap_dimer"]

    full = [p for p in patches if p["conservation"] == 1.0]
    full.sort(key=lambda p: (-(p["functional"] > 0), -p["n_acidic_conserved"],
                             -p["functional"], -p["n"]))

    print("\n" + "=" * 78)
    print("CONSERVED PATCHES RANKED BY FUNCTIONAL RELEVANCE + His-PAIRING CAPACITY")
    print("=" * 78)
    print(f"{'center':>7} {'dom':<10} {'n':>3} {'acidC':>6} {'EGFsite':>8} "
          f"{'dimer':>6} {'cetux':>6}")
    for p in full[:15]:
        print(f"{p['center_aa']}{p['center']:<6} {p['domain']:<10} {p['n']:>3} "
              f"{p['n_acidic_conserved']:>6} {p['overlap_ligand']:>8} "
              f"{p['overlap_dimer']:>6} {p['overlaps_cetuximab']:>6}")

    cet_lig = len(cetux & lig_up)
    print(f"\nSanity check — cetuximab footprint overlaps the EGF site at "
          f"{cet_lig}/{len(cetux)} residues (cetuximab is a ligand-blocking antibody).")
    acid_cet = [f"{f['aa']}{f['uniprot']}" for f in
                json.load(open(os.path.join(OUT, "cetuximab_footprint.json")))
                if f["aa"] in ("D", "E")]
    print(f"Carboxylates in the cetuximab footprint (His-pairing partners): "
          f"{', '.join(acid_cet) if acid_cet else 'none'}  "
          f"-> {len(acid_cet)}/{len(cetux)} residues")

    json.dump(full, open(os.path.join(OUT, "epitope_shortlist.json"), "w"),
              indent=1, default=float)
    print(f"\nWrote {len(full)} fully conserved patches to results/epitope_shortlist.json")


if __name__ == "__main__":
    main()
