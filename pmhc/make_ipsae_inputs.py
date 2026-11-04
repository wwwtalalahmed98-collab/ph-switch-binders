"""Convert our saved PAE into the file ipsae.py actually reads, and check it lines up.

WHY THIS EXISTS. ipSAE is the interface-quality score the competition encourages, and
this project has reported "cannot compute ipSAE" in four consecutive challenges. The
reason was never subtle: the score needs the PAE matrix and we never saved it. We save
it now - as .npy, which ipsae.py does not accept. It takes AF2 input as a .json holding

    pae    an N x N matrix      (key 'pae', or 'predicted_aligned_error')
    plddt  a length-N array     (key 'plddt'; zeros are substituted if absent)

alongside the .pdb, where **N is the number of CA atoms in that PDB**. Nothing validates
that. If N disagrees, ipsae.py indexes a matrix of the wrong size and returns numbers
that look fine and mean nothing - which is the exact class of silent failure that cost
this project two rounds of TfR1 measurement.

So this script refuses to write a file whose PAE size does not equal the CA count, and
says what it found.
"""
import json
import os
import sys

import numpy as np
from Bio import PDB


def ca_count(pdb_path):
    """Residues ipsae.py will index: one CA per standard residue, in file order."""
    m = PDB.PDBParser(QUIET=True).get_structure("s", pdb_path)[0]
    n = 0
    for ch in m:
        for r in ch:
            if r.id[0] == " " and "CA" in r:
                n += 1
    return n


def convert(pdb_path, pae_path, plddt_path=None, out_path=None, iptm=None, ptm=None):
    pae = np.asarray(np.load(pae_path), dtype=float)
    if pae.ndim == 3 and pae.shape[0] == 1:
        pae = pae[0]
    if pae.ndim != 2 or pae.shape[0] != pae.shape[1]:
        raise SystemExit("%s: PAE is shape %s, expected a square matrix"
                         % (os.path.basename(pae_path), pae.shape))

    n_ca = ca_count(pdb_path)
    if pae.shape[0] != n_ca:
        raise SystemExit(
            "%s: PAE is %dx%d but the structure has %d CA atoms. ipsae.py does not "
            "check this and would silently score the wrong residues."
            % (os.path.basename(pdb_path), pae.shape[0], pae.shape[1], n_ca))

    if plddt_path and os.path.exists(plddt_path):
        plddt = np.asarray(np.load(plddt_path), dtype=float).reshape(-1)
        if plddt.size != n_ca:
            raise SystemExit("%s: pLDDT has %d entries, structure has %d CA atoms"
                             % (os.path.basename(plddt_path), plddt.size, n_ca))
        # ColabDesign reports pLDDT on 0-1; ipsae.py and AlphaFold use 0-100.
        if float(plddt.max()) <= 1.0:
            plddt = plddt * 100.0
    else:
        plddt = np.zeros(n_ca)

    data = {"pae": pae.tolist(), "plddt": plddt.tolist()}
    if iptm is not None:
        data["iptm"] = float(iptm)
    if ptm is not None:
        data["ptm"] = float(ptm)
    out_path = out_path or pdb_path.replace(".pdb", "_pae.json")
    json.dump(data, open(out_path, "w"))
    return out_path, n_ca, bool(plddt.any())


def main(argv):
    if len(argv) < 3:
        print("usage: make_ipsae_inputs.py <design.pdb> <design_pae.npy> "
              "[design_plddt.npy]")
        return 2
    pdb_path, pae_path = argv[1], argv[2]
    plddt_path = argv[3] if len(argv) > 3 else None
    out, n, has_plddt = convert(pdb_path, pae_path, plddt_path)
    print("wrote %s  (%d residues, pLDDT %s)"
          % (out, n, "included" if has_plddt else "ZEROS - ipSAE will still run but "
             "its pLDDT columns are meaningless"))
    print("\nthen:  python ipsae.py %s %s 10 10" % (out, pdb_path))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
