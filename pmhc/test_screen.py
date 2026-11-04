"""Check the counter-screen's contact logic against the chain layout it will really see.

The screen scores a three-chain target (HLA + b2m + peptide) and save_pdb merges those
into one chain, so the peptide ends up as the LAST ten residues of a 384-residue chain.
The apo state has no peptide at all and must not be mistaken for a holo one.

Both of those are assumptions about byte layout, which is the category of assumption that
has already cost this project three rounds of measurement. They are checked here rather
than on a GPU.
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
OUT = os.path.join(ROOT, "data", "pmhc", "states")
SCRIPT = os.path.join(ROOT, "kaggle", "cellPMHCc_screen.py")
T3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
      "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
      "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
      "TYR": "Y", "VAL": "V"}


def load_screen():
    src = io.open(SCRIPT, encoding="utf-8").read()
    marker = "\nmain()\n"
    if not src.endswith(marker):
        raise SystemExit("screen script no longer ends with a bare main(); update this test")
    head = src[: -len(marker)]
    # Remove only the two blocks that need a Kaggle box - the colabdesign install and
    # the AlphaFold weight download. Exact text, because the imports this test needs sit
    # BETWEEN them and a blind cut takes those out too (which is how the first version
    # of this loader failed with NameError: PDB).
    blocks = [
        'try:\n    import colabdesign  # noqa: F401\nexcept ImportError:\n'
        '    sh("pip -q install git+https://github.com/sokrypton/ColabDesign.git'
        '@v1.1.1 biopython")\n',
        'if not os.path.exists("%s/params/params_model_1_ptm.npz" % W):\n'
        '    os.makedirs("%s/params" % W, exist_ok=True)\n'
        '    sh("wget -q -O %s/af2.tar https://storage.googleapis.com/alphafold/"\n'
        '       "alphafold_params_2022-12-06.tar" % W)\n'
        '    sh("tar -xf %s/af2.tar -C %s/params" % (W, W))\n'
        '    os.remove("%s/af2.tar" % W)\n'
        'if len(glob.glob("%s/params/params_model_*.npz" % W)) == 0:\n'
        '    raise SystemExit("no AlphaFold parameters - refusing to start")\n',
    ]
    for b in blocks:
        if b not in head:
            raise SystemExit("screen script changed; this test strips a block that is "
                             "no longer present:\n%s" % b[:90])
        head = head.replace(b, "", 1)
    ns = {"__name__": "screen_under_test"}
    exec(compile(head, SCRIPT, "exec"), ns)
    return ns


def merged(state_pdb, path, anchor, offset, n_binder=80):
    """One merged target chain, as save_pdb writes it, plus a binder chain."""
    m = PDB.PDBParser(QUIET=True).get_structure("t", state_pdb)[0]
    st = Structure("s")
    mod = Model(0)
    st.add(mod)
    tgt = Chain("A")
    mod.add(tgt)
    n = 0
    for cid in ("A", "B", "C"):
        if cid not in m:
            continue
        for r in m[cid]:
            if r.id[0] != " " or r.get_resname() not in T3:
                continue
            n += 1
            nr = Residue((" ", n, " "), r.get_resname(), "")
            for a in r:
                if a.element == "H":
                    continue
                nr.add(Atom(a.get_name(), a.coord, 0.0, 1.0, " ",
                            a.get_name(), a.serial_number, a.element))
            tgt.add(nr)
    b = Chain("B")
    mod.add(b)
    base = anchor.mean(axis=0) + offset
    for i in range(n_binder):
        r = Residue((" ", i + 1, " "), "ALA", "")
        r.add(Atom("CA", base + np.array([i * 0.6, 0.0, 0.0], dtype="f4"),
                   0.0, 1.0, " ", "CA", i + 1, "C"))
        b.add(r)
    w = PDB.PDBIO()
    w.set_structure(st)
    w.save(path)
    return path, n


def main():
    ns = load_screen()
    contacts, PEP_G12D, PEP_WT = ns["contacts"], ns["PEP_G12D"], ns["PEP_WT"]

    g12d = os.path.join(OUT, "target_state_g12d.pdb")
    wt = os.path.join(OUT, "target_state_wt.pdb")
    apo = os.path.join(OUT, "target_state_apo.pdb")
    for f in (g12d, wt, apo):
        if not os.path.exists(f):
            raise SystemExit("run pmhc/build_counter_states.py first (missing %s)" % f)

    m = PDB.PDBParser(QUIET=True).get_structure("g", g12d)[0]
    pep = [r for r in m["C"] if r.id[0] == " "]
    asp6 = np.array([a.coord for a in pep[5] if a.element != "H"])

    near = np.array([0.0, 0.0, 3.0], dtype="f4")
    far = np.array([60.0, 60.0, 60.0], dtype="f4")
    ok = True

    p, n = merged(g12d, os.path.join(OUT, "_t_g12d.pdb"), asp6, near)
    r = contacts(p, "g12d", PEP_G12D)
    print("g12d, binder on Asp6   target %d res -> on_mutation %d on_peptide %d on_groove %d"
          % (n, r["on_mutation"], r["on_peptide"], r["on_groove"]))
    if r["on_mutation"] == 0:
        print("FAIL: no mutation contact detected in the merged 3-chain layout")
        ok = False
    os.remove(p)

    p, n = merged(wt, os.path.join(OUT, "_t_wt.pdb"), asp6, near)
    r = contacts(p, "wt", PEP_WT)
    print("wt,   binder on pos 6   target %d res -> on_mutation %d on_peptide %d on_groove %d"
          % (n, r["on_mutation"], r["on_peptide"], r["on_groove"]))
    if r["on_mutation"] == 0:
        print("FAIL: wt state did not resolve its peptide (it ends in Gly, not Asp)")
        ok = False
    os.remove(p)

    p, n = merged(apo, os.path.join(OUT, "_t_apo.pdb"), asp6, near)
    r = contacts(p, "apo", None)
    print("apo,  binder in groove  target %d res -> on_mutation %d on_peptide %d on_groove %d"
          % (n, r["on_mutation"], r["on_peptide"], r["on_groove"]))
    if r["on_peptide"] != 0 or r["on_mutation"] != 0:
        print("FAIL: apo state reported peptide contacts - it has no peptide")
        ok = False
    if r["on_groove"] == 0:
        print("FAIL: apo state reported no groove contact either, so the binder was "
              "placed wrong and this case proves nothing")
        ok = False
    os.remove(p)

    p, n = merged(g12d, os.path.join(OUT, "_t_far.pdb"), asp6, far)
    r = contacts(p, "g12d", PEP_G12D)
    print("g12d, binder 60 A away  target %d res -> on_mutation %d on_peptide %d on_groove %d"
          % (n, r["on_mutation"], r["on_peptide"], r["on_groove"]))
    if any(r.values()):
        print("FAIL: a binder 60 A away reports contacts")
        ok = False
    os.remove(p)

    # the wrong peptide must be rejected, not silently tolerated
    p, _ = merged(g12d, os.path.join(OUT, "_t_bad.pdb"), asp6, near)
    try:
        contacts(p, "g12d", PEP_WT)
        print("FAIL: a g12d structure was accepted as wild type - the sequence check "
              "does not actually guard anything")
        ok = False
    except RuntimeError:
        print("mismatched peptide correctly rejected")
    os.remove(p)

    print("\n%s" % ("PASS - merged 3-chain targets, the apo state, and a wrong peptide "
                    "are all handled" if ok else "TEST FAILED"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
