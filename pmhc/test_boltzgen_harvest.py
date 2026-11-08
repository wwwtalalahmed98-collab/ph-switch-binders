"""Prove the BoltzGen harvester can read a BoltzGen-style mmCIF before we pay for one.

WHY THIS EXISTS. A completed BoltzGen run produced all 48 inverse-folded structures
and the harvester parsed 0 of them, because it picked the designed chain by chain
LETTER and BoltzGen does not carry the ids from the input yaml through to its output.
Seventy minutes of GPU bought one line saying "skipped 48". This test costs no GPU and
fails in a second.

WHAT IT PINS:
  - the mmCIF reader returns one sequence per chain, from CA atoms, in file order
  - the designed chain is found by SEQUENCE, whatever letter it was given
  - it is still found when the letters are the ones we would have guessed wrong
  - an ambiguous file is refused rather than guessed at
  - a target chain alone is never mistaken for a design

Run:  python pmhc/test_boltzgen_harvest.py
"""
import importlib.util
import os
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "kaggle", "cellPMHCk_boltzgen.py")

T1 = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
      "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
      "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
      "Y": "TYR", "V": "VAL"}

HEADER = """data_design
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
"""


def write_cif(path, chains):
    """chains: {chain_id: one_letter_sequence}. Writes N, CA, C per residue."""
    n = 0
    with open(path, "w") as fh:
        fh.write(HEADER)
        for cid, seq in chains.items():
            for i, aa in enumerate(seq, 1):
                for at, el in (("N", "N"), ("CA", "C"), ("C", "C")):
                    n += 1
                    fh.write("ATOM %d %s %s %s %s %d %d %s "
                             "%.3f %.3f %.3f 1.00 50.00\n"
                             % (n, el, at, T1[aa], cid, i, i, cid,
                                i * 1.0, 0.0, 0.0))
        fh.write("#\n")


def load_module():
    spec = importlib.util.spec_from_file_location("bg", SRC)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main():
    bg = load_module()
    tmp = tempfile.mkdtemp(prefix="bgtest_")
    fails = []

    mhc = "AG" * 137 + "A"          # 275 aa, stands in for the HLA heavy chain
    pep = "VVVGADGVGK"              # 10 aa, the real peptide
    design = "EKLAREAVKQ" * 8       # 80 aa
    targets = {mhc, pep}

    # 1. reader returns one sequence per chain
    p = os.path.join(tmp, "a.cif")
    write_cif(p, {"A": mhc, "B": pep, "C": design})
    ch = bg.seqs_from_cif(p)
    if {k: len(v) for k, v in ch.items()} != {"A": 275, "B": 10, "C": 80}:
        fails.append("reader returned %s" % {k: len(v) for k, v in ch.items()})
    if ch.get("B") != pep:
        fails.append("peptide read back as %r" % ch.get("B"))

    # 2. the design is found when it sits on chain C -- the letter the OLD code
    #    treated as the peptide. This is the exact case that failed.
    cid, seq = bg.pick_design(ch, targets)
    if seq != design:
        fails.append("design on chain C not found (got %r on %s)"
                     % (seq and seq[:12], cid))

    # 3. and when the labels are shuffled any other way
    for layout in ({"A": design, "B": mhc, "C": pep},
                   {"X": pep, "Y": design, "Z": mhc},
                   {"B": mhc, "D": design, "E": pep}):
        q = os.path.join(tmp, "l.cif")
        write_cif(q, layout)
        cid, seq = bg.pick_design(bg.seqs_from_cif(q), targets)
        if seq != design:
            fails.append("layout %s -> got %r" % (list(layout), seq and seq[:12]))

    # 4. a file with two non-target chains of designable length is AMBIGUOUS and
    #    must be refused, not guessed
    q = os.path.join(tmp, "amb.cif")
    write_cif(q, {"A": mhc, "B": design, "C": "KLKLKLKLKL" * 7})
    cid, seq = bg.pick_design(bg.seqs_from_cif(q), targets)
    if cid is not None:
        fails.append("ambiguous file was not refused (picked %s)" % cid)

    # 5. target chains alone yield nothing
    q = os.path.join(tmp, "tgt.cif")
    write_cif(q, {"A": mhc, "B": pep})
    cid, seq = bg.pick_design(bg.seqs_from_cif(q), targets)
    if cid is not None:
        fails.append("a target-only file produced a 'design' on %s" % cid)

    # 6. an unparsable file returns {} rather than raising
    q = os.path.join(tmp, "junk.cif")
    open(q, "w").write("not a cif at all\n")
    if bg.seqs_from_cif(q) != {}:
        fails.append("junk file did not return {}")

    print("checks run: 6 groups")
    if fails:
        print()
        print("FAIL")
        for f in fails:
            print("  - %s" % f)
        return 1
    print("PASS - the designed chain is found by sequence on every layout tried,")
    print("       ambiguity is refused, and junk does not raise.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
