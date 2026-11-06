"""Dock the published pMHC scaffold set onto our KRAS-G12D / HLA-A*11:01 target.

WHY THIS EXISTS - AND WHY WE ABANDONED THE PREVIOUS APPROACH.

Three generation runs (phases 2, 2B, 3) tried to invent binder backbones from scratch
by back-propagating through AlphaFold2 (AfDesign "hallucination"). They produced 62
candidate structures and exactly zero that survived our own buildability screen. The
diagnosis was unambiguous: when we asked ProteinMPNN to write sequences for those
backbones, it returned glycine/proline soup on all 15 we tried, 720 sequences, zero
survivors. ProteinMPNN only does that when the backbone is not a real protein fold -
it places Gly and Pro where the local geometry demands them, so Gly/Pro soup is a
readout saying "these dihedrals do not belong to a foldable chain".

Liu et al. (Science 389:386-391, 2025) solved exactly our problem - de novo binders to
peptide-MHC-I with specificity for a mutated peptide residue - and say plainly why the
from-scratch route is so wasteful: "only a small fraction of trajectories yield
backbones that interact extensively with the peptide while largely avoiding the MHC."
Their fix was to stop inventing backbones and reuse a set of scaffolds already shown to
arc over an MHC-I groove, docking them onto each new target. They published the
scaffolds and the pipeline (Zenodo 10.5281/zenodo.15169815).

This script is step 2b-i of their pipeline, reimplemented so it runs on CPU here with no
PyRosetta and no PyMOL. Their align_chainB.py uses PyMOL's `align`, which does its own
sequence alignment internally; we do the same thing explicitly with Biopython so that the
alignment quality is a number we can read rather than a side effect we have to trust.

HOW THE DOCKING WORKS. Each scaffold file holds two chains: chain A is the designed
binder, chain B is the pMHC it was originally designed against (MHC platform plus its
peptide, concatenated and renumbered as one chain). MHC-I platforms are structurally
near-identical across alleles, so superposing a scaffold's own MHC onto ours carries its
binder into our target's frame in approximately the pose it was designed to adopt. We
then keep our target and the moved binder, and throw the scaffold's original MHC away.

WHAT THIS DOES NOT DO. It does not optimise anything. A docked pose is a hypothesis
about where a binder could sit, nothing more; most will clash or miss the peptide, which
is why every pose is measured and filtered here rather than assumed to be usable. The
sequence is still the scaffold's original sequence at this stage - ProteinMPNN redesigns
it in the next step. Nothing here has been near a GPU or a structure predictor.

Run:  python pmhc/dock_scaffolds.py <scaffolds-dir>
"""
import glob
import json
import os
import sys

import numpy as np
from Bio import PDB
from Bio.Align import PairwiseAligner
from Bio.PDB.Superimposer import Superimposer

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STATES = os.path.join(ROOT, "data", "pmhc", "states")
OUT = os.path.join(ROOT, "data", "pmhc", "docked")
RESULTS = os.path.join(ROOT, "pmhc", "results")

TARGET = os.path.join(STATES, "target_design.pdb")
PEPTIDE = "VVVGADGVGK"          # KRAS G12D 10-mer, P01116 residues 7-16
MUT_POS = 6                      # 1-based position of the Asp that WT has as Gly

# Geometry thresholds. Heavy-atom distances, in angstroms.
CONTACT = 4.5                    # a contact, as used throughout this project
CLASH = 2.2                      # closer than this is a steric clash, not a contact
MIN_ALIGN_RES = 120              # fewer matched MHC residues than this and the
                                 # superposition is not describing the same fold
MAX_ALIGN_RMSD = 2.5             # platforms are conserved; worse than this is suspect

T3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
      "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
      "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
      "TYR": "Y", "VAL": "V"}


def residues(chain):
    """Standard amino acid residues with a CA, in file order."""
    return [r for r in chain
            if r.id[0] == " " and r.get_resname() in T3 and "CA" in r]


def seq_of(res_list):
    return "".join(T3[r.get_resname()] for r in res_list)


def load_target():
    """Our target: chain A = MHC platform (179 aa), chain C = peptide (10 aa)."""
    m = PDB.PDBParser(QUIET=True).get_structure("t", TARGET)[0]
    heavy = residues(m["A"])
    pep = residues(m["C"])
    if seq_of(pep) != PEPTIDE:
        raise SystemExit("target peptide is %r, expected %r - wrong target file"
                         % (seq_of(pep), PEPTIDE))
    return m, heavy, pep


PEP_LEN_RANGE = (8, 11)          # MHC-I peptides; anything else means we mis-split


def split_scaffold_chain_b(sc_b, target_heavy_seq):
    """Separate a scaffold's chain B into (its MHC residues, its peptide residues).

    The scaffold stores MHC and peptide concatenated in one chain with no marker
    between them. Rather than assuming the peptide is "the last 9 or 10 residues" -
    the class of assumption that has already cost this project three rounds of
    measurement - we find the MHC by aligning against our own MHC sequence and treat
    whatever trails off the end of that alignment as the peptide.

    This uses LOCAL alignment, not global. A global aligner, told to align our
    179-residue MHC against the scaffold's 183-191 residue MHC+peptide, happily
    stretches the alignment over the peptide too - it has no reason not to, since
    extending costs less than the terminal gap. That is exactly what happened on the
    first run of this script: 143 of 152 scaffolds came back as "no peptide found",
    because the peptide had been absorbed into the MHC alignment. Local alignment
    stops at the end of the genuinely similar region, which is the MHC.

    Returns (mhc_residues, peptide_residues, pairs) where pairs is the residue
    correspondence used for superposition.
    """
    sc_seq = seq_of(sc_b)
    al = PairwiseAligner(mode="local", open_gap_score=-10, extend_gap_score=-0.5,
                         match_score=2, mismatch_score=-1)
    aln = al.align(target_heavy_seq, sc_seq)[0]
    t_blocks, s_blocks = aln.aligned
    pairs = []
    for (t0, t1), (s0, s1) in zip(t_blocks, s_blocks):
        for k in range(t1 - t0):
            pairs.append((t0 + k, s0 + k))
    if not pairs:
        return [], [], []
    last_sc_mhc = max(p[1] for p in pairs)
    mhc = sc_b[: last_sc_mhc + 1]
    pep = sc_b[last_sc_mhc + 1:]
    return mhc, pep, pairs


BACKBONE = ("N", "CA", "C", "O")


def atom_coords(res_list, backbone_only=False):
    out = []
    for r in res_list:
        for a in r:
            if a.element == "H":
                continue
            if backbone_only and a.get_name() not in BACKBONE:
                continue
            out.append(a.coord)
    return np.asarray(out, dtype=float)


def min_dists(a, b):
    """For each atom in a, distance to the nearest atom in b."""
    if len(a) == 0 or len(b) == 0:
        return np.array([])
    d = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1)
    return d.min(axis=1)


def measure(binder_res, heavy, pep):
    """Interface description of one docked pose. Pure geometry, no prediction."""
    b = atom_coords(binder_res)
    pep_atoms = atom_coords(pep)
    heavy_atoms = atom_coords(heavy)
    mut_atoms = atom_coords([pep[MUT_POS - 1]])

    d_pep = min_dists(b, pep_atoms)
    d_mhc = min_dists(b, heavy_atoms)
    d_mut = min_dists(b, mut_atoms)

    n_pep = int((d_pep <= CONTACT).sum())
    n_mhc = int((d_mhc <= CONTACT).sum())
    n_mut = int((d_mut <= CONTACT).sum())
    clashes = int(((d_pep < CLASH).sum()) + ((d_mhc < CLASH).sum()))

    # Backbone-only clashes are the number that actually decides whether a pose is
    # usable. The all-atom count above includes side chains of the scaffold's
    # ORIGINAL sequence, every one of which ProteinMPNN is about to replace, so a
    # side-chain clash is not evidence against the pose. A backbone clash is, because
    # sequence design cannot move a backbone.
    bb = atom_coords(binder_res, backbone_only=True)
    tgt_bb = np.vstack([atom_coords(heavy, backbone_only=True),
                        atom_coords(pep, backbone_only=True)])
    d_bb = min_dists(bb, tgt_bb)
    bb_clashes = int((d_bb < CLASH).sum())

    # Per-peptide-residue contact, so we can see which part of the peptide is read.
    per_res = {}
    for i, r in enumerate(pep, start=1):
        dd = min_dists(b, atom_coords([r]))
        per_res[i] = int((dd <= CONTACT).sum())

    total = n_pep + n_mhc
    return {
        "binder_len": len(binder_res),
        "n_pep_contacts": n_pep,
        "n_mhc_contacts": n_mhc,
        "n_mut_contacts": n_mut,
        "pep_fraction": round(n_pep / total, 4) if total else 0.0,
        "clashes": clashes,
        "bb_clashes": bb_clashes,
        "per_peptide_residue": per_res,
    }


def fresh_residue(r, new_number):
    """A standalone copy of residue `r`, renumbered, with no parent.

    Built atom by atom on purpose. The obvious `copy.deepcopy(r)` is a trap: a
    Biopython residue holds a reference to its parent chain, so deepcopy duplicates
    the entire structure graph behind it - for every residue. Writing 279 residues
    per pose that way took ~28 seconds a scaffold, which would have made this script
    a 70-minute job instead of a two-minute one. Constructing residues directly also
    avoids the "id already used" renumbering that the deepcopy version triggered.
    """
    nr = PDB.Residue.Residue((" ", new_number, " "), r.get_resname(), "")
    for a in r:
        if a.element == "H":
            continue
        nr.add(PDB.Atom.Atom(a.get_name(), a.coord, a.get_bfactor(),
                             a.get_occupancy(), " ", a.get_name(),
                             a.serial_number, a.element))
    return nr


def write_pose(path, heavy, pep, binder_res):
    """Write a 2-chain file: A = target (MHC then peptide, renumbered 1..189),
    B = docked binder (renumbered 1..n). ProteinMPNN then fixes A and designs B,
    and the peptide occupies A residues len(heavy)+1 .. len(heavy)+len(pep)."""
    st = PDB.Structure.Structure("d")
    mod = PDB.Model.Model(0)
    st.add(mod)

    ca = PDB.Chain.Chain("A")
    mod.add(ca)
    for n, r in enumerate(list(heavy) + list(pep), start=1):
        ca.add(fresh_residue(r, n))

    cb = PDB.Chain.Chain("B")
    mod.add(cb)
    for i, r in enumerate(binder_res, start=1):
        cb.add(fresh_residue(r, i))

    io = PDB.PDBIO()
    io.set_structure(st)
    io.save(path)


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    scaf_dir = argv[1]
    files = sorted(glob.glob(os.path.join(scaf_dir, "scaf*.pdb")))
    if not files:
        raise SystemExit("no scaf*.pdb in %s" % scaf_dir)

    os.makedirs(OUT, exist_ok=True)
    os.makedirs(RESULTS, exist_ok=True)

    tgt_model, heavy, pep = load_target()
    heavy_seq = seq_of(heavy)
    print("target: MHC platform %d aa, peptide %s (mutation at position %d = %s)"
          % (len(heavy), seq_of(pep), MUT_POS, seq_of(pep)[MUT_POS - 1]))
    print("scaffolds to dock: %d" % len(files))
    print()

    parser = PDB.PDBParser(QUIET=True)
    rows = []
    rejected = {"no_alignment": 0, "short_alignment": 0, "bad_rmsd": 0, "no_peptide": 0}

    for f in files:
        name = os.path.basename(f)[:-4]
        m = parser.get_structure("s", f)[0]
        if "A" not in m or "B" not in m:
            rejected["no_alignment"] += 1
            continue
        binder = residues(m["A"])
        sc_b = residues(m["B"])

        sc_mhc, sc_pep, pairs = split_scaffold_chain_b(sc_b, heavy_seq)
        if not pairs:
            rejected["no_alignment"] += 1
            continue
        if not (PEP_LEN_RANGE[0] <= len(sc_pep) <= PEP_LEN_RANGE[1]):
            rejected["no_peptide"] += 1
            continue
        # Keep only pairs that fall inside the scaffold's MHC part.
        pairs = [(ti, si) for ti, si in pairs if si < len(sc_mhc)]
        if len(pairs) < MIN_ALIGN_RES:
            rejected["short_alignment"] += 1
            continue

        fixed = [heavy[ti]["CA"] for ti, si in pairs]
        moving = [sc_b[si]["CA"] for ti, si in pairs]
        sup = Superimposer()
        sup.set_atoms(fixed, moving)
        if sup.rms > MAX_ALIGN_RMSD:
            rejected["bad_rmsd"] += 1
            continue

        # Move the whole scaffold, then measure only its binder against our target.
        sup.apply(list(m.get_atoms()))

        met = measure(binder, heavy, pep)
        met["name"] = name
        met["align_rmsd"] = round(float(sup.rms), 3)
        met["align_residues"] = len(pairs)
        met["scaffold_peptide"] = seq_of(sc_pep)
        met["scaffold_seq"] = seq_of(binder)
        rows.append(met)

        write_pose(os.path.join(OUT, "%s_docked.pdb" % name), heavy, pep, binder)

    print("docked and measured: %d of %d" % (len(rows), len(files)))
    print("rejected: %s" % rejected)
    print()

    if not rows:
        raise SystemExit("nothing docked - the superposition assumption is wrong")

    rmsds = sorted(r["align_rmsd"] for r in rows)
    print("--- superposition quality (MHC platform onto ours) ---")
    print("    RMSD over matched CA: min %.2f  median %.2f  max %.2f A"
          % (rmsds[0], rmsds[len(rmsds) // 2], rmsds[-1]))
    print()

    bb_free = [r for r in rows if r["bb_clashes"] == 0]
    touch_pep = [r for r in rows if r["n_pep_contacts"] > 0]
    touch_mut = [r for r in rows if r["n_mut_contacts"] > 0]
    print("--- interface, before any filtering (denominator %d docked poses) ---"
          % len(rows))
    print("    backbone-clash-free             : %d" % len(bb_free))
    print("    touching the peptide at all     : %d" % len(touch_pep))
    print("    touching the Asp6 mutation      : %d" % len(touch_mut))
    pf = sorted(r["pep_fraction"] for r in rows)
    print("    peptide fraction of interface   : min %.2f median %.2f max %.2f"
          % (pf[0], pf[len(pf) // 2], pf[-1]))
    print()
    print("    NOTE: 'clash' below is all-atom and includes the scaffold's original")
    print("    side chains, which ProteinMPNN replaces. 'bbcl' is backbone-only and")
    print("    is the count that sequence design cannot fix.")
    print()

    rows.sort(key=lambda r: (r["bb_clashes"], -r["n_mut_contacts"],
                             -r["pep_fraction"], -r["n_pep_contacts"]))
    print("--- top 20 poses: no backbone clash first, then mutation contact ---")
    print("    %-10s %5s %4s %5s %5s %6s %5s %6s %6s"
          % ("name", "len", "mut", "pep", "mhc", "pepfr", "bbcl", "clash", "rmsd"))
    for r in rows[:20]:
        print("    %-10s %5d %4d %5d %5d %6.2f %5d %6d %6.2f"
              % (r["name"], r["binder_len"], r["n_mut_contacts"],
                 r["n_pep_contacts"], r["n_mhc_contacts"], r["pep_fraction"],
                 r["bb_clashes"], r["clashes"], r["align_rmsd"]))

    out_json = os.path.join(RESULTS, "docked_poses.json")
    json.dump(rows, open(out_json, "w"), indent=1)
    print()
    print("wrote %d pose files to %s" % (len(rows), OUT))
    print("wrote measurements to %s" % out_json)
    print()
    print("NOT DONE YET: these are raw docked poses with their ORIGINAL scaffold")
    print("sequences. Nothing has been redesigned, predicted, or scored for binding.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
