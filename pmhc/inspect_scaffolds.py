"""Inspect the Baker-lab pMHC scaffold set before we commit to using it.

WHY THIS EXISTS. The published pipeline (Liu et al., Science 389:386, 2025) ships 148
scaffolds that were shown experimentally to arc over an MHC-I groove. Their docking
script works by superposing each scaffold's OWN pMHC (which it calls chain B) onto the
new target, then keeping the scaffold's binder (chain A). That only works if the chain
layout is what the script assumes, and the paper does not state it.

So before building anything on top of these files, read them and report what is actually
in them: how many chains, which chain is the binder, how long each one is, and whether
the chain B really looks like an MHC heavy chain plus a peptide. Everything downstream
depends on these being what we think they are, and this project has already lost three
rounds of work to byte-layout assumptions that were never checked.

Run:  python pmhc/inspect_scaffolds.py <path-to-scaffolds-dir>
"""
import collections
import glob
import os
import sys

T3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
      "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
      "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
      "TYR": "Y", "VAL": "V"}

# The first 25 residues of the HLA class I heavy chain alpha1 domain are near-identical
# across alleles. If a scaffold's chain B starts like this, it is an MHC heavy chain.
HLA_A_START = "GSHSMRYF"


def read_ca(path):
    """CA atoms per chain, in file order: (resnum, one-letter)."""
    chains = collections.OrderedDict()
    for line in open(path):
        if not line.startswith("ATOM"):
            continue
        if line[12:16].strip() != "CA":
            continue
        res = line[17:20].strip()
        if res not in T3:
            continue
        chains.setdefault(line[21], []).append((line[22:27].strip(), T3[res]))
    return chains


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    d = argv[1]
    files = sorted(glob.glob(os.path.join(d, "scaf*.pdb")))
    if not files:
        raise SystemExit("no scaf*.pdb found in %s" % d)

    rows = []
    for f in files:
        ch = read_ca(f)
        seqs = {c: "".join(x[1] for x in v) for c, v in ch.items()}
        rows.append((os.path.basename(f), seqs))

    print("scaffold files found: %d" % len(rows))
    print()
    print("--- three examples, in full ---")
    for name, seqs in rows[:3]:
        print("%s  chains %s" % (name, {c: len(s) for c, s in seqs.items()}))
        for c, s in seqs.items():
            print("    chain %s (%3d): %s" % (c, len(s), s))
        print()

    # Which chains exist, and is the layout uniform?
    layouts = collections.Counter(tuple(sorted(s.keys())) for _, s in rows)
    print("--- chain layouts across all %d scaffolds ---" % len(rows))
    for lay, n in layouts.most_common():
        print("    chains %s : %d scaffolds" % (list(lay), n))
    print()

    # Binder = chain A by the script's assumption. Check lengths are design-like.
    a_len = [len(s.get("A", "")) for _, s in rows]
    b_len = [len(s.get("B", "")) for _, s in rows]
    print("--- chain A (the binder, per align_chainB.py) ---")
    print("    lengths: min %d  median %d  max %d"
          % (min(a_len), sorted(a_len)[len(a_len) // 2], max(a_len)))
    print("    distribution: %s" % sorted(collections.Counter(a_len).items()))
    print()
    print("--- chain B (should be MHC heavy chain + peptide) ---")
    print("    lengths: %s" % sorted(collections.Counter(b_len).items()))

    # Does chain B actually look like an HLA heavy chain?
    hla_like = sum(1 for _, s in rows if s.get("B", "").startswith(HLA_A_START))
    print("    chain B starting with the HLA-A alpha1 motif %r: %d of %d"
          % (HLA_A_START, hla_like, len(rows)))
    if hla_like == 0:
        print("    WARNING: no chain B starts with that motif. Either the scaffolds")
        print("    are numbered/trimmed differently, or chain B is not the heavy chain.")
        print("    Printing the first 20 residues of chain B for three scaffolds:")
        for name, s in rows[:3]:
            print("        %s: %s" % (name, s.get("B", "")[:20]))

    # Cysteines matter: Adaptyv expresses these, and free cysteines are a liability.
    cys = sum(1 for _, s in rows if "C" in s.get("A", ""))
    print()
    print("--- composition notes on the binder chain ---")
    print("    scaffolds whose binder contains cysteine: %d of %d" % (cys, len(rows)))
    print("    (the scaffold sequence is only a starting point - ProteinMPNN will")
    print("     redesign it, so this is informational, not a filter)")

    over = [n for n, s in rows if len(s.get("A", "")) > 250]
    print()
    print("--- competition limits (Track 3: 10-250 aa) ---")
    print("    binders longer than 250 aa: %d %s" % (len(over), over[:10]))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
