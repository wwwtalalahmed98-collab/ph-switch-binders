"""Is any candidate TfR1 epitope shadowed by a sugar?

EPITOPE_CHOICE.md ranked patches on conservation and shape and flagged glycans as the most
likely way that analysis is wrong. TfR1 is glycosylated, 11QC is a crystal structure with no
sugar modelled, and a patch that looks open in the coordinates can be unreachable in the
assay. This checks it.

Two hazards, not one:
  1. a glycan near a patch in BOTH species blocks the binder outright;
  2. a glycan near a patch in only ONE species breaks cross-reactivity specifically, which
     is objective 2 of the challenge. That one is invisible to a conservation score, because
     the sequon can appear or vanish without the patch residues themselves changing.

Glycosylation sites come from UniProt's own annotations for each species, cross-checked
against N-X-S/T sequons computed from the sequence (X != Pro). Where the two disagree the
script says so rather than silently preferring one.

REACH. An N-linked glycan is not a point. A complex biantennary glycan extends roughly
20-25 A from the attachment nitrogen and sweeps a cone, so the numbers below use a 20 A
advisory radius and report the actual distance rather than a yes/no, because the real
exclusion depends on the glycan's orientation, which nothing here models.
"""
import io
import json
import os
import re

import numpy as np
import requests
from Bio import PDB
from Bio.Align import PairwiseAligner, substitution_matrices

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "data", "tfr1")
REACH = 20.0
ECTO = (121, 760)

THREE2ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
             "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
             "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
             "TYR": "Y", "VAL": "V"}


def fasta(p):
    return "".join(l.strip() for l in io.open(p) if not l.startswith(">"))


def annotated_sites(acc):
    """N-linked glycosylation positions UniProt actually annotates."""
    j = requests.get("https://rest.uniprot.org/uniprotkb/%s.json" % acc, timeout=40).json()
    out = {}
    for f in j.get("features", []):
        if f.get("type") != "Glycosylation":
            continue
        desc = (f.get("description") or "")
        pos = f["location"]["start"]["value"]
        out[pos] = desc
    return out


def sequons(seq):
    """N-X-S/T with X != Pro. Returns the Asn position, 1-based."""
    return {m.start() + 1 for m in re.finditer(r"(?=N[^P][ST])", seq)}


def hm_map(h, m):
    al = PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score = -11, -1
    a = al.align(h, m)[0]
    out = {}
    for (hs, he), (ms, me) in zip(a.aligned[0], a.aligned[1]):
        for k in range(he - hs):
            out[hs + k + 1] = ms + k + 1
    return out


def main():
    human = fasta(os.path.join(D, "P02786.fasta"))
    mouse = fasta(os.path.join(D, "Q62351.fasta"))
    h2m = hm_map(human, mouse)
    m2h = {v: k for k, v in h2m.items()}

    ha, ma = annotated_sites("P02786"), annotated_sites("Q62351")
    hs, ms = sequons(human), sequons(mouse)
    print("UniProt-annotated N-glycosylation sites")
    print("  human: %s" % ", ".join("N%d" % n for n in sorted(ha)))
    print("  mouse: %s" % ", ".join("N%d" % n for n in sorted(ma)))
    print("sequons computed from sequence (N-X-S/T, X!=Pro)")
    print("  human: %d, mouse: %d" % (len(hs), len(ms)))
    miss = sorted(n for n in ha if n not in hs)
    if miss:
        print("  NOTE: annotated but not a sequon in our scan: %s"
              % ", ".join("N%d" % n for n in miss))

    # Species comparison on the human numbering, restricted to the supplied ectodomain.
    heq = {n for n in hs if ECTO[0] <= n <= ECTO[1]}
    meq = {m2h[n] for n in ms if n in m2h and ECTO[0] <= m2h[n] <= ECTO[1]}
    both, honly, monly = heq & meq, heq - meq, meq - heq
    print("\nectodomain sequons, in human numbering:")
    print("  in both species : %d  %s" % (len(both), ", ".join("N%d" % n for n in sorted(both))))
    print("  human only      : %d  %s" % (len(honly), ", ".join("N%d" % n for n in sorted(honly))))
    print("  mouse only      : %d  %s" % (len(monly), ", ".join("N%d" % n for n in sorted(monly))))
    print("  (a species-specific glycan is a cross-reactivity hazard a conservation score"
          " cannot see)")

    st = PDB.MMCIFParser(QUIET=True).get_structure("x", os.path.join(D, "11QC.cif"))[0]
    chain = [c for c in st if sum(1 for r in c if r.id[0] == " ") > 100][0]
    pos = {}
    for r in chain:
        if r.id[0] == " " and r.get_resname() in THREE2ONE and "CA" in r:
            pos[r.id[1]] = r["CA"].get_coord()
    print("\nchain %s, %d residues with coordinates" % (chain.id, len(pos)))

    def site_xyz(ns):
        return {n: pos[n] for n in ns if n in pos}

    gb, gh, gm = site_xyz(both), site_xyz(honly), site_xyz(monly)
    print("resolved in the structure: %d shared, %d human-only, %d mouse-only"
          % (len(gb), len(gh), len(gm)))

    rows = json.load(open(os.path.join(D, "tfr1_mouse_check.json")))
    print("\n%-26s %4s %9s %9s %9s  %s"
          % ("patch", "n", "shared", "human-only", "mouse-only", "verdict"))
    out = []
    for r in rows:
        mem = [n for n in r["members"] if n in pos]
        if not mem:
            continue
        c = np.mean([pos[n] for n in mem], axis=0)

        def near(g):
            if not g:
                return None, None
            d = {k: float(np.linalg.norm(v - c)) for k, v in g.items()}
            k = min(d, key=d.get)
            return k, d[k]
        kb, db = near(gb)
        kh, dh = near(gh)
        km, dm = near(gm)
        flags = []
        if db is not None and db < REACH:
            flags.append("blocked both species (N%d, %.0f A)" % (kb, db))
        if dh is not None and dh < REACH:
            flags.append("human-only glycan N%d at %.0f A" % (kh, dh))
        if dm is not None and dm < REACH:
            flags.append("mouse-only glycan N%d at %.0f A" % (km, dm))
        verdict = "; ".join(flags) if flags else "clear"
        out.append(dict(patch=r["patch"], n=len(mem), shared=db, human_only=dh,
                        mouse_only=dm, verdict=verdict))
        print("%-26s %4d %9s %9s %9s  %s"
              % (r["patch"], len(mem),
                 "%.0f" % db if db is not None else "-",
                 "%.0f" % dh if dh is not None else "-",
                 "%.0f" % dm if dm is not None else "-", verdict))

    json.dump(out, open(os.path.join(D, "tfr1_glycan_check.json"), "w"), indent=0)
    print("\nDistances are centre-of-patch to glycosylation-site CA, in angstrom.")
    print("Under %.0f A the glycan may reach the patch; this models no glycan geometry,"
          " so treat it as advisory." % REACH)
    print("wrote data/tfr1/tfr1_glycan_check.json")


if __name__ == "__main__":
    main()
