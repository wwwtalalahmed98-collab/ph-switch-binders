"""Phase 2 prep: build BindCraft target constructs and verify epitope accessibility.

Strips the cetuximab Fab from 6ARU, writes the full tethered ECD and domain-trimmed
design constructs, and checks whether each candidate epitope is occluded by the rest
of the receptor. An epitope that is exposed in the trimmed construct but buried in the
full ECD would produce binders that cannot engage the protein Adaptyv actually assays.
"""
import json
import os

from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import protein_letters_3to1

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STRUCT = os.path.join(ROOT, "data", "structures")
OUT = os.path.join(ROOT, "targets")
OFFSET = 24  # PDB + 24 = UniProt

MAX_SASA = {
    "A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "E": 223, "Q": 225,
    "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240,
    "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174,
}

# Epitopes from Phase 1, in UniProt numbering.
EPITOPES = {
    "A_domainI_ligand": dict(
        hotspots=[44, 45, 46, 49, 52, 74, 75],
        carboxylates=[45, 46, 75, 97],
        construct=(25, 188),   # domain I (L1)
        note="domain I ligand-binding face; blocks EGF binding"),
    "B_domainII_dimerarm": dict(
        hotspots=[274, 275, 277, 278, 279, 281, 282, 283],
        carboxylates=[278, 282, 321, 602],
        construct=(189, 334),  # domain II (CR1)
        note="dimerization arm; blocks receptor dimerization"),
    "C_domainII_K294": dict(
        hotspots=[282, 283, 293, 294, 296, 319, 320, 321],
        carboxylates=[257, 282, 319, 320, 321],
        construct=(189, 334),
        note="backup; highest carboxylate density"),
}


class ChainA(Select):
    def __init__(self, lo=None, hi=None):
        self.lo, self.hi = lo, hi

    def accept_chain(self, chain):
        return chain.id == "A"

    def accept_residue(self, residue):
        if residue.id[0] != " ":
            return False
        if self.lo is None:
            return True
        return self.lo <= residue.id[1] <= self.hi

    def accept_atom(self, atom):
        return atom.element != "H" and not atom.is_disordered() or atom.get_altloc() in (" ", "A")


def sasa_map(structure_entity):
    ShrakeRupley().compute(structure_entity, level="R")
    out = {}
    for r in structure_entity.get_residues():
        if r.id[0] != " ":
            continue
        try:
            aa = protein_letters_3to1[r.get_resname()]
        except KeyError:
            continue
        out[r.id[1]] = (aa, r.sasa / MAX_SASA.get(aa, 200))
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    parser = PDBParser(QUIET=True)
    io = PDBIO()

    # Full tethered ECD, Fab removed.
    st = parser.get_structure("x", os.path.join(STRUCT, "6ARU.pdb"))
    io.set_structure(st)
    full_path = os.path.join(OUT, "EGFR_ECD_tethered.pdb")
    io.save(full_path, ChainA())
    full_sasa = sasa_map(parser.get_structure("f", full_path)[0])
    print(f"Wrote {os.path.basename(full_path)} "
          f"({len(full_sasa)} residues, cetuximab Fab removed)")

    report = {}
    for name, ep in EPITOPES.items():
        lo, hi = ep["construct"]
        pdb_lo, pdb_hi = lo - OFFSET, hi - OFFSET
        st2 = parser.get_structure("x", os.path.join(STRUCT, "6ARU.pdb"))
        io.set_structure(st2)
        path = os.path.join(OUT, f"EGFR_{name}.pdb")
        io.save(path, ChainA(pdb_lo, pdb_hi))
        trim_sasa = sasa_map(parser.get_structure("t", path)[0])

        print(f"\n=== Epitope {name} ===")
        print(f"  {ep['note']}")
        print(f"  construct UniProt {lo}-{hi} (PDB {pdb_lo}-{pdb_hi}), "
              f"{len(trim_sasa)} residues -> {os.path.basename(path)}")
        print(f"  {'residue':>9} {'rSASA full':>11} {'rSASA trim':>11} {'occluded':>9}")
        occluded = []
        for up in ep["hotspots"]:
            p = up - OFFSET
            aa_f, rf = full_sasa.get(p, ("?", float("nan")))
            aa_t, rt = trim_sasa.get(p, ("?", float("nan")))
            flag = "YES" if (rt - rf) > 0.25 else ""
            if flag:
                occluded.append(up)
            print(f"  {aa_f}{up:>8} {rf:>11.2f} {rt:>11.2f} {flag:>9}")
        buried_full = [up for up in ep["hotspots"]
                       if full_sasa.get(up - OFFSET, ("?", 0))[1] < 0.15]
        print(f"  hotspots buried in the full ECD (rSASA < 0.15): "
              f"{buried_full if buried_full else 'none'}")
        print(f"  hotspots gaining >0.25 rSASA on trimming (trim artefact): "
              f"{occluded if occluded else 'none'}")

        report[name] = dict(
            construct_pdb=os.path.basename(path),
            construct_uniprot=[lo, hi], construct_pdb_range=[pdb_lo, pdb_hi],
            hotspots_uniprot=ep["hotspots"],
            hotspots_pdb=[h - OFFSET for h in ep["hotspots"]],
            carboxylates_uniprot=ep["carboxylates"],
            carboxylates_pdb=[c - OFFSET for c in ep["carboxylates"]],
            buried_in_full_ecd=buried_full, trim_artefact=occluded,
            note=ep["note"])

    json.dump(report, open(os.path.join(OUT, "targets.json"), "w"), indent=1)
    print(f"\nWrote {os.path.join('targets', 'targets.json')}")
    print("\nBindCraft hotspot strings (PDB numbering, chain A):")
    for name, r in report.items():
        print(f"  {name}: {','.join('A' + str(h) for h in r['hotspots_pdb'])}")


if __name__ == "__main__":
    main()
