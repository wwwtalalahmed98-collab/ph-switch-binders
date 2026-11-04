"""Prove the contact measurement works on the chain layout save_pdb actually produces.

Version 2 of the probe spent 32 minutes of GPU time, completed a design, and then died
in measure() because save_pdb merges the two target chains into one. Each of these
failures costs half an hour, so the fix is tested here against a synthetic structure
that reproduces the exact layout from the failure - chains {'A': 189, 'B': 77} - before
any more GPU time is spent.

Three cases, and the test is only useful because two of them are meant to be ZERO:

  touching   a binder placed against Asp6          -> on_mutation > 0
  missing    a binder placed against Gly7 instead  -> on_mutation == 0, on_peptide > 0
  far away   a binder 50 A from everything         -> all zero

The middle case is the one that matters. It is the failure mode this whole challenge
turns on: a design that contacts the peptide but not the mutation reads nothing the
wild type lacks. If measure() cannot tell that case from the first, every number we
report downstream is worthless.
"""
import io
import os
import sys

import numpy as np
from Bio import PDB
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TARGET = os.path.join(ROOT, "data", "pmhc", "states", "target_design.pdb")
SCRIPT = os.path.join(ROOT, "kaggle", "cellPMHCa_probe.py")
OUT = os.path.join(ROOT, "data", "pmhc", "states")


def load_probe():
    """Import the probe's functions without running main()."""
    src = io.open(SCRIPT, encoding="utf-8").read()
    marker = "\nmain()\n"
    if not src.endswith(marker):
        raise SystemExit("probe script no longer ends with a bare main() call; "
                         "this test strips that line and must be updated")
    ns = {"__name__": "probe_under_test"}
    exec(compile(src[: -len(marker)], SCRIPT, "exec"), ns)
    return ns


def build(path, anchor_atoms, offset, n_binder=77):
    """Merged target chain A (groove + peptide), plus a binder chain B near `anchor`."""
    p = PDB.PDBParser(QUIET=True)
    m = p.get_structure("t", TARGET)[0]
    T3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
          "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
          "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
          "TYR": "Y", "VAL": "V"}

    st = Structure("s")
    mod = Model(0)
    st.add(mod)
    merged = Chain("A")
    mod.add(merged)

    n = 0
    for cid in ("A", "C"):                      # groove first, then peptide
        for r in m[cid]:
            if r.id[0] != " " or r.get_resname() not in T3:
                continue
            n += 1
            nr = Residue((" ", n, " "), r.get_resname(), "")  # renumbered from 1
            for a in r:
                if a.element == "H":
                    continue
                nr.add(Atom(a.get_name(), a.coord, 0.0, 1.0, " ",
                            a.get_name(), a.serial_number, a.element))
            merged.add(nr)

    binder = Chain("B")
    mod.add(binder)
    base = anchor_atoms.mean(axis=0) + offset
    for i in range(n_binder):
        r = Residue((" ", i + 1, " "), "ALA", "")
        r.add(Atom("CA", base + np.array([i * 0.6, 0.0, 0.0], dtype="f4"),
                   0.0, 1.0, " ", "CA", i + 1, "C"))
        binder.add(r)

    io_ = PDB.PDBIO()
    io_.set_structure(st)
    io_.save(path)
    return path


def main():
    ns = load_probe()
    measure, EXPECTED_A, PEP = ns["measure"], ns["EXPECTED_A"], ns["PEP_G12D"]

    m = PDB.PDBParser(QUIET=True).get_structure("t", TARGET)[0]
    pep = [r for r in m["C"] if r.id[0] == " "]
    asp6 = np.array([a.coord for a in pep[5] if a.element != "H"])
    gly7 = np.array([a.coord for a in pep[6] if a.element != "H"])

    cases = [
        ("touching Asp6", asp6, np.array([0.0, 0.0, 3.0], dtype="f4")),
        ("touching Gly7 only", gly7, np.array([0.0, 0.0, 3.0], dtype="f4")),
        ("far away", asp6, np.array([50.0, 50.0, 50.0], dtype="f4")),
    ]
    results = {}
    for label, anchor, off in cases:
        path = os.path.join(OUT, "_test_%s.pdb" % label.split()[0])
        build(path, anchor, off)
        chains = {c.id: sum(1 for r in c if r.id[0] == " ")
                  for c in PDB.PDBParser(QUIET=True).get_structure("x", path)[0]}
        r = measure(path, EXPECTED_A, PEP)
        results[label] = r
        print("%-20s chains %s -> on_mutation %d  on_peptide %d  on_groove %d"
              % (label, chains, r["on_mutation"], r["on_peptide"], r["on_groove"]))
        os.remove(path)

    print()
    ok = True
    if results["touching Asp6"]["on_mutation"] == 0:
        print("FAIL: a binder against Asp6 reports no mutation contact")
        ok = False
    if results["touching Gly7 only"]["on_mutation"] != 0:
        print("FAIL: a binder that avoids Asp6 still reports a mutation contact - "
              "the metric cannot tell the two apart, which is the whole challenge")
        ok = False
    if results["touching Gly7 only"]["on_peptide"] == 0:
        print("FAIL: the Gly7 binder reports no peptide contact either, so it was "
              "simply placed wrong and this case proves nothing")
        ok = False
    if any(results["far away"][k] for k in ("on_mutation", "on_peptide", "on_groove")):
        print("FAIL: a binder 50 A away reports contacts")
        ok = False
    print("PASS - the chains {'A': 189, 'B': 77} layout is handled, and a design that "
          "touches the peptide\n       but not the mutation is correctly scored zero "
          "on the mutation." if ok else "TEST FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
