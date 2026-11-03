"""Build the Challenge 5 design target and define the three states it must tell apart.

THE TARGET. KRAS G12D peptide VVVGADGVGK (P01116 residues 7-16) in the groove of
HLA-A*11:01, with beta-2-microglobulin underneath. PDB 7OW6, chains A (heavy chain),
B (b2m), C (peptide). Chains D and E are a TCR and are discarded.

THE THREE STATES. The competition scores specificity before affinity, so the design is
measured against three things, not one:

    g12d  HLA + b2m + VVVGADGVGK   <- must bind
    wt    HLA + b2m + VVVGAGGVGK   <- must not bind  (one residue different)
    apo   HLA + b2m                <- must not bind  (no peptide at all)

The only difference between g12d and wt is Gly to Asp at peptide position 6. Everything
else is identical, so the background MUST be identical too. If the g12d background came
from 7OW6 and the wt background from 7OW5, any difference we measured could come from
two crystals rather than from the mutation. So the states share one background, taken
from 7OW6, and 7OW5 is used only to confirm the wild-type sequence is what we think it
is and to tell us which residues both structures actually model.

WHY THIS SCRIPT MEASURES INSTEAD OF ASSUMING. On challenge 3 a crop opened a cavity that
does not exist in the real protein and the optimiser dived into it. The groove here is
formed by the a1/a2 platform, roughly the first 180 residues; a3 and b2m sit underneath.
Cropping to the platform would halve the size of every prediction, which matters when
every design has to be folded three times. But it also slices the protein, so this script
computes how much buried hydrophobic surface that cut would expose and prints the number
whether or not anyone reads it. The decision is recorded, not assumed.

NUMBERING. save_pdb renumbers residues from 1, which silently destroyed two rounds of
measurement on TfR1. Nothing downstream may identify a chain by residue number. The
manifest records, for every chain, the ordered list of original numbers, so position i in
the crop maps back to kept[i] - and downstream code identifies the target by reproducing
its sequence residue for residue.
"""
import copy
import json
import os

import numpy as np
from Bio import PDB
from Bio.PDB.SASA import ShrakeRupley

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "data", "pmhc")
OUT = os.path.join(D, "states")

MUT_POS = 6              # peptide position carrying Gly12Asp
PEP_G12D = "VVVGADGVGK"
PEP_WT = "VVVGAGGVGK"
PLATFORM_LAST = 180      # a1/a2 groove platform, heavy chain 1..180
HOTSPOT_RADIUS = 8.0     # A around the mutated residue, for groove hotspots
DESIGN_RADIUS = 20.0     # A around the peptide: the design target. See crop_sweep.py.
MIN_SEGMENT = 3          # drop dangling 1-2 residue fragments, they are noise

T3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
      "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
      "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
      "TYR": "Y", "VAL": "V"}
HYDROPHOBIC = set("AVLIMFWCY")


def residues(chain):
    """Ordered standard residues, as a writer would emit them."""
    return [r for r in chain if r.id[0] == " " and r.get_resname() in T3]


def seq(rs):
    return "".join(T3[r.get_resname()] for r in rs)


class _Sel(PDB.Select):
    def __init__(self, keep):
        self.keep = {c: set(n) for c, n in keep.items()}

    def accept_model(self, m):
        return m.id == 0

    def accept_chain(self, c):
        return c.id in self.keep

    def accept_residue(self, r):
        return r.id[0] == " " and r.id[1] in self.keep[r.get_parent().id]

    def accept_atom(self, a):
        return (not a.is_disordered()) or a.get_altloc() in (" ", "A")


def sasa_by_residue(model, keep):
    """Per-residue solvent exposure for a sub-structure, keyed (chain, resnum)."""
    sub = copy.deepcopy(model)
    for ch in list(sub):
        if ch.id not in keep:
            sub.detach_child(ch.id)
            continue
        allowed = set(keep[ch.id])
        for r in list(ch):
            if r.id[0] != " " or r.id[1] not in allowed:
                ch.detach_child(r.id)
    ShrakeRupley().compute(sub, level="R")
    return {(ch.id, r.id[1]): r.sasa for ch in sub for r in ch}


def main():
    if not os.path.isdir(OUT):
        os.makedirs(OUT)
    p = PDB.MMCIFParser(QUIET=True)
    m6 = p.get_structure("7OW6", os.path.join(D, "7OW6.cif"))[0]
    m5 = p.get_structure("7OW5", os.path.join(D, "7OW5.cif"))[0]

    hv, b2m, pep = residues(m6["A"]), residues(m6["B"]), residues(m6["C"])

    # --- verify before anything is written ---
    if seq(pep) != PEP_G12D:
        raise SystemExit("7OW6 chain C is %s, expected %s" % (seq(pep), PEP_G12D))
    if seq(residues(m5["C"])) != PEP_WT:
        raise SystemExit("7OW5 chain C is %s, expected %s"
                         % (seq(residues(m5["C"])), PEP_WT))
    diff = [i + 1 for i, (a, b) in enumerate(zip(PEP_G12D, PEP_WT)) if a != b]
    if diff != [MUT_POS]:
        raise SystemExit("peptides differ at %s, expected only %d" % (diff, MUT_POS))
    print("peptide check: 7OW6 %s, 7OW5 %s, differ only at position %d"
          % (seq(pep), seq(residues(m5["C"])), MUT_POS))

    # 7OW6 chain A starts with a GS cloning remnant that 7OW5 lacks. Use only the heavy
    # chain and b2m residues present in BOTH structures, so the background is the same
    # one regardless of which crystal a number was read from.
    common_hv = sorted({r.id[1] for r in hv} & {r.id[1] for r in residues(m5["A"])})
    common_b2m = sorted({r.id[1] for r in b2m} & {r.id[1] for r in residues(m5["B"])})
    dropped = [r.id[1] for r in hv if r.id[1] not in set(common_hv)]
    print("heavy chain %d res, %d shared with 7OW5%s"
          % (len(hv), len(common_hv),
             "; dropped %s" % dropped if dropped else ""))
    print("b2m %d res, %d shared" % (len(b2m), len(common_b2m)))

    keep_full = {"A": common_hv, "B": common_b2m, "C": [r.id[1] for r in pep]}
    keep_plat = {"A": [n for n in common_hv if n <= PLATFORM_LAST],
                 "C": [r.id[1] for r in pep]}

    # --- does cropping to the platform open a face that is buried in reality? ---
    full = sasa_by_residue(m6, keep_full)
    plat = sasa_by_residue(m6, keep_plat)
    aa_of = {r.id[1]: T3[r.get_resname()] for r in hv}
    new_phobic = sum(plat[k] - full[k] for k in plat
                     if k[0] == "A" and plat[k] > full[k] + 1.0
                     and aa_of.get(k[1]) in HYDROPHOBIC)
    pep_area = sum(full[k] for k in full if k[0] == "C")
    print("\nplatform crop (heavy chain 1-%d, no b2m):" % PLATFORM_LAST)
    print("  %d of %d heavy-chain residues kept" % (len(keep_plat["A"]), len(common_hv)))
    print("  new exposed hydrophobic surface %.0f A^2" % new_phobic)
    print("  for scale, the peptide's own surface in the complex is %.0f A^2" % pep_area)
    print("  ratio %.2fx" % (new_phobic / pep_area if pep_area else float("nan")))

    # --- hotspots: the mutated residue, and the groove around it ---
    mut = pep[MUT_POS - 1]
    if T3[mut.get_resname()] != "D":
        raise SystemExit("peptide position %d is %s, not Asp"
                         % (MUT_POS, mut.get_resname()))
    mut_atoms = np.array([a.coord for a in mut if a.element != "H"])
    groove = []
    allowed_hv = set(common_hv)
    for r in hv:
        if r.id[1] not in allowed_hv:
            continue
        d = np.linalg.norm(
            np.array([a.coord for a in r if a.element != "H"])[:, None, :]
            - mut_atoms[None, :, :], axis=-1).min()
        if d <= HOTSPOT_RADIUS:
            groove.append((r.id[1], T3[r.get_resname()], round(float(d), 2)))
    print("\nhotspots - peptide Asp%d plus %d heavy-chain residues within %.0f A:"
          % (MUT_POS, len(groove), HOTSPOT_RADIUS))
    for n, aa, d in groove:
        print("  A %s%-4d %.2f A" % (aa, n, d))

    # --- the crop we actually design against ---
    # crop_sweep.py compares the candidates. No crop is clean: every one opens about
    # 250 A^2 of hydrophobic surface that is buried in the real complex, roughly the
    # size of the peptide's own footprint. Keeping b2m while dropping a3 is WORSE, not
    # better - it exposes b2m's underside - which only became visible once the audit
    # counted every chain instead of the heavy chain alone.
    #
    # We accept a cut face here for a reason that does not apply on an affinity-scored
    # challenge: the face is identical in all three states, so a binder that grips it
    # binds all three and the counter-screen throws it out. It costs design budget, not
    # correctness. The alternative is the 384-residue full complex, which at three
    # folds per design would buy very few designs in the days remaining.
    pep_atoms = np.array([a.coord for r in pep for a in r if a.element != "H"])
    near = []
    for r in hv:
        if r.id[1] not in allowed_hv:
            continue
        xyz = np.array([a.coord for a in r if a.element != "H"])
        if np.linalg.norm(xyz[:, None, :] - pep_atoms[None, :, :],
                          axis=-1).min() <= DESIGN_RADIUS:
            near.append(r.id[1])
    runs, cur = [], []
    for n in sorted(near):
        if cur and n == cur[-1] + 1:
            cur.append(n)
        else:
            if len(cur) >= MIN_SEGMENT:
                runs += cur
            cur = [n]
    if len(cur) >= MIN_SEGMENT:
        runs += cur
    keep_design = {"A": runs, "C": keep_full["C"]}
    print("\ndesign target: %d A around the peptide, %d heavy-chain residues + peptide"
          % (DESIGN_RADIUS, len(runs)))

    # --- write the design target and the manifest ---
    io = PDB.PDBIO()
    io.set_structure(m6)
    for tag, keep in (("target_full", keep_full), ("target_platform", keep_plat),
                      ("target_design", keep_design)):
        path = os.path.join(OUT, "%s.pdb" % tag)
        io.save(path, _Sel(keep))
        print("\nwrote %s" % path)

    hv_seq = seq([r for r in hv if r.id[1] in allowed_hv])
    b2m_seq = seq([r for r in b2m if r.id[1] in set(common_b2m)])
    manifest = dict(
        source=dict(design_target="7OW6", wild_type_reference="7OW5"),
        chains=dict(heavy="A", b2m="B", peptide="C"),
        kept=dict(A=common_hv, B=common_b2m, C=keep_full["C"]),
        sequences=dict(heavy=hv_seq, b2m=b2m_seq),
        states=dict(
            g12d=dict(peptide=PEP_G12D, bind=True, note="the design target"),
            wt=dict(peptide=PEP_WT, bind=False,
                    note="same background, Asp6 to Gly only"),
            apo=dict(peptide=None, bind=False,
                     note="heavy chain and b2m, no peptide"),
        ),
        mutation=dict(peptide_position=MUT_POS, wild_type="G", mutant="D"),
        hotspots=dict(peptide=[MUT_POS], heavy=[n for n, _, _ in groove],
                      radius=HOTSPOT_RADIUS),
        crop_audit=dict(platform_last=PLATFORM_LAST,
                        platform_new_hydrophobic_A2=round(float(new_phobic), 1),
                        peptide_surface_A2=round(float(pep_area), 1)),
        design_target=dict(file="target_design.pdb", radius=DESIGN_RADIUS,
                           chains=["A", "C"], kept_A=runs,
                           note="generation runs against this; the three-state "
                                "counter-screen runs against the full complex"),
    )
    mpath = os.path.join(OUT, "states.json")
    json.dump(manifest, open(mpath, "w"), indent=1)
    print("wrote %s" % mpath)
    print("\nheavy %d aa, b2m %d aa, peptide %d aa; full complex %d aa"
          % (len(hv_seq), len(b2m_seq), len(PEP_G12D),
             len(hv_seq) + len(b2m_seq) + len(PEP_G12D)))


if __name__ == "__main__":
    main()
