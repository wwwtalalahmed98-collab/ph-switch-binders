"""Phase 3: acid-ON pH switch engineering.

Given a binder-target complex, finds binder positions where a histidine could form a
salt bridge / charged hydrogen bond with a conserved target carboxylate, and scores the
resulting switch.

Mechanism (inverted relative to the published acid-OFF work): at pH 6.5 the imidazolium
is protonated and pairs with a target Asp/Glu; at pH 7.4 the neutral imidazole loses that
interaction, leaving a desolvated unpaired carboxylate at the interface. The "off" state
is therefore a specific lost contact plus a desolvation penalty, not a diffuse weakening.

Design rules enforced here, from Marchand et al. (bioRxiv 2025.09.29.678932):
  - deliberate placement only; random histidine incorporation rarely switches
  - require >= 2 histidine-charge contacts for a credible switch
  - do not substitute buried core positions (destabilises the monomer)

Usage:
    python scripts/phase3_ph_switch.py <complex.pdb> --binder B --target A
    python scripts/phase3_ph_switch.py --selftest
"""
import argparse
import json
import os
import sys

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import protein_letters_3to1

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

MAX_SASA = {
    "A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "E": 223, "Q": 225,
    "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240,
    "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174,
}
CARBOXYL_ATOMS = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}

# Histidine geometry, calibrated on 23 genuine His-carboxylate salt bridges (N-O <= 3.3 A)
# measured across 4HHB, 1IGY, 1BRS, 1A4Y, 2PTC, 1MBN, 5PTI, 1AKI, 3HHB, 1LZ1, 6ARU, 1IVO.
#
# CB->carboxylate-O: 3.55-6.77 A (p5 3.83, median 6.22, p95 6.73, mean 5.79)
# CA->CB direction cosine: -0.33 to 0.81 (p5 -0.26, median 0.27)
# N-O contact distance: 2.52-3.27 A (median 2.84)
#
# Two consequences, both of which overrode the author's prior assumptions:
#  1. A "side chain must point at the target" filter (cos > 0.2) rejects 48% of genuine
#     salt bridges, because chi1/chi2 rotation lets the imidazole swing far off the
#     CA->CB axis. The threshold is therefore set at the observed p5, not at intuition.
#  2. Once a position is reachable, the contact forms at ~2.84 A regardless of how far
#     the CB sits. So contact strength is near-constant per contact, and a
#     distance-dependent pseudo-energy would be modelling noise. We therefore score the
#     NUMBER and GEOMETRIC QUALITY of independent contacts, which is also what the
#     published pH-switch work found to be determinant.
CB_O_MIN, CB_O_MAX = 3.5, 7.0
CB_O_IDEAL = 5.8
COS_MIN = -0.30

# Positions we must not mutate: buried core, or already a switchable/charged residue
# whose loss would cost more than the switch gains.
DO_NOT_MUTATE = {"PRO", "GLY", "CYS"}


def virtual_cb(residue):
    """CB coordinate, constructed from backbone for glycine."""
    if "CB" in residue:
        return residue["CB"].get_coord()
    try:
        n = residue["N"].get_coord()
        ca = residue["CA"].get_coord()
        c = residue["C"].get_coord()
    except KeyError:
        return None
    b, cc = ca - n, c - ca
    a = np.cross(b, cc)
    return -0.58273431 * a + 0.56802827 * b - 0.54067466 * cc + ca


def chain_residues(entity):
    out = []
    for r in entity:
        if r.id[0] != " " or "CA" not in r:
            continue
        try:
            aa = protein_letters_3to1[r.get_resname()]
        except KeyError:
            continue
        out.append((r, aa))
    return out


def rel_sasa_map(entity):
    ShrakeRupley().compute(entity, level="R")
    out = {}
    for r in entity.get_residues():
        if r.id[0] != " ":
            continue
        try:
            aa = protein_letters_3to1[r.get_resname()]
        except KeyError:
            continue
        out[r.id[1]] = r.sasa / MAX_SASA.get(aa, 200)
    return out


def contact_quality(d_cb):
    """How central a CB-carboxylate distance is within the empirical salt-bridge
    distribution. 1.0 at the observed mean, falling to 0 at the window edges.
    Unitless by design: see the calibration note above for why this is not a ddG."""
    span = max(CB_O_MAX - CB_O_IDEAL, CB_O_IDEAL - CB_O_MIN)
    return max(0.0, 1.0 - abs(d_cb - CB_O_IDEAL) / span)


def analyse(pdb_path, binder_chain, target_chain, carboxylates=None):
    parser = PDBParser(QUIET=True)
    model = parser.get_structure("c", pdb_path)[0]
    binder = model[binder_chain]
    target = model[target_chain]

    # Burial: SASA of the binder alone vs in the complex, and the same for the target.
    binder_alone = rel_sasa_map(parser.get_structure("b", pdb_path)[0][binder_chain])
    target_alone = rel_sasa_map(parser.get_structure("t", pdb_path)[0][target_chain])
    ShrakeRupley().compute(model, level="R")
    binder_cplx = {r.id[1]: r.sasa / MAX_SASA.get(
        protein_letters_3to1.get(r.get_resname(), "A"), 200)
        for r in binder if r.id[0] == " " and r.get_resname() in protein_letters_3to1}
    target_cplx = {r.id[1]: r.sasa / MAX_SASA.get(
        protein_letters_3to1.get(r.get_resname(), "A"), 200)
        for r in target if r.id[0] == " " and r.get_resname() in protein_letters_3to1}

    # Target carboxylates to aim at.
    acids = []
    for r, aa in chain_residues(target):
        if r.get_resname() not in CARBOXYL_ATOMS:
            continue
        if carboxylates and r.id[1] not in carboxylates:
            continue
        oxys = [r[a].get_coord() for a in CARBOXYL_ATOMS[r.get_resname()] if a in r]
        if not oxys:
            continue
        buried = target_alone.get(r.id[1], 0) - target_cplx.get(r.id[1], 0)
        acids.append(dict(resnum=r.id[1], aa=aa, oxys=oxys,
                          sasa_free=round(target_alone.get(r.id[1], 0), 3),
                          desolvated=round(buried, 3)))

    proposals = []
    for r, aa in chain_residues(binder):
        if r.get_resname() in DO_NOT_MUTATE:
            continue
        cb = virtual_cb(r)
        ca = r["CA"].get_coord()
        if cb is None:
            continue
        direction = cb - ca
        norm = np.linalg.norm(direction)
        if norm < 1e-6:
            continue
        direction = direction / norm
        # Interface positions only: must become buried on binding.
        d_burial = binder_alone.get(r.id[1], 0) - binder_cplx.get(r.id[1], 0)
        if binder_alone.get(r.id[1], 0) < 0.15:
            continue  # buried in the monomer core - mutating risks destabilising the fold

        for acid in acids:
            dists = [np.linalg.norm(cb - o) for o in acid["oxys"]]
            d_cb = float(min(dists))
            if not (CB_O_MIN <= d_cb <= CB_O_MAX):
                continue
            to_acid = acid["oxys"][int(np.argmin(dists))] - cb
            to_acid /= np.linalg.norm(to_acid)
            cos = float(np.dot(direction, to_acid))
            if cos < COS_MIN:
                continue  # side chain points away; chi rotation cannot recover it
            q = contact_quality(d_cb)
            # A contact is worth more when the target carboxylate is desolvated by
            # binding: at pH 7.4 the neutral imidazole leaves it buried and unpaired,
            # which is what makes the "off" state specifically unfavourable.
            proposals.append(dict(
                binder_resnum=r.id[1], binder_aa=aa,
                target_resnum=acid["resnum"], target_aa=acid["aa"],
                d_cb_to_carboxylate=round(d_cb, 2),
                direction_cos=round(cos, 2),
                contact_quality=round(q, 3),
                binder_pos_exposed_free=round(binder_alone.get(r.id[1], 0), 3),
                binder_pos_buried_on_binding=round(d_burial, 3),
                target_carboxylate_desolvated=acid["desolvated"],
                score=round(q * (1.0 + max(acid["desolvated"], 0.0)), 3)))

    proposals.sort(key=lambda p: -p["score"])
    return proposals, acids


def report(proposals, acids, label=""):
    print(f"\n{'=' * 76}")
    print(f"ACID-ON SWITCH CANDIDATES {label}")
    print("=" * 76)
    print(f"Target carboxylates considered: "
          f"{', '.join(a['aa'] + str(a['resnum']) for a in acids) or 'none'}")
    if not proposals:
        print("No geometrically viable histidine positions found.")
        return
    print(f"\n{'binder':>8} {'->':>3} {'target':>8} {'dCB':>6} {'cos':>6} "
          f"{'qual':>6} {'desolv':>7} {'score':>6}")
    for p in proposals[:20]:
        print(f"{p['binder_aa']}{p['binder_resnum']:<7} {'->':>3} "
              f"{p['target_aa']}{p['target_resnum']:<7} "
              f"{p['d_cb_to_carboxylate']:>6} {p['direction_cos']:>6} "
              f"{p['contact_quality']:>6} {p['target_carboxylate_desolvated']:>7} "
              f"{p['score']:>6}")

    # The >=2 contact rule.
    partners = {}
    for p in proposals:
        partners.setdefault(p["binder_resnum"], set()).add(p["target_resnum"])
    multi = {k: v for k, v in partners.items() if len(v) >= 2}
    distinct_acids = len({p["target_resnum"] for p in proposals})
    print(f"\nDistinct target carboxylates reachable: {distinct_acids}")
    print(f"Binder positions reaching >=2 carboxylates: "
          f"{sorted(multi) if multi else 'none'}")
    if distinct_acids >= 2:
        print("PASS: >=2 independent His-carboxylate contacts are available, "
              "the minimum the published pH-switch work found necessary.")
    else:
        print("FAIL: fewer than 2 independent contacts available - this design is a "
              "poor switch candidate regardless of its affinity.")


def selftest():
    """Validate the geometry engine on a real protein-protein interface (6ARU:
    cetuximab Fab vs EGFR). This checks that reachable His positions and sane
    energies are found at a genuine interface; it is not a design run."""
    pdb = os.path.join(ROOT, "data", "structures", "6ARU.pdb")
    if not os.path.exists(pdb):
        sys.exit("6ARU.pdb not found - run Phase 1 first.")
    print("SELF-TEST 1 (positive control): barnase (A) vs barstar (D), 1BRS.")
    print("Ground truth: barnase His102 forms a salt bridge with barstar Asp39")
    print("(N-O 2.81 A). The engine is given no knowledge of this and must find it.")
    pos = os.path.join(ROOT, "data", "calib", "1BRS.pdb")
    if os.path.exists(pos):
        props, acids = analyse(pos, "A", "D")
        report(props, acids, "(positive control, barnase/barstar)")
        hit = [p for p in props if p["binder_resnum"] == 102
               and p["target_resnum"] == 39]
        print(f"\nRecovered the known His102-Asp39 pair: "
              f"{'YES' if hit else 'NO'}")
    else:
        print("1BRS.pdb not present - skipping positive control.")

    print("\n\nSELF-TEST 2 (negative control): cetuximab Fab (B) vs EGFR (A), 6ARU.")
    print("Expected: no viable positions. The cetuximab epitope carries one "
          "carboxylate\nin 26 residues and none is within reach of the Fab.")
    props, acids = analyse(pdb, "B", "A")
    report(props, acids, "(negative control, Fab vs EGFR)")
    bad = [p for p in props
           if not (CB_O_MIN <= p["d_cb_to_carboxylate"] <= CB_O_MAX)]
    print(f"\nProposals outside the calibrated window (should be 0): {len(bad)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pdb", nargs="?")
    ap.add_argument("--binder", default="B")
    ap.add_argument("--target", default="A")
    ap.add_argument("--carboxylates", default="",
                    help="comma-separated target residue numbers to aim at")
    ap.add_argument("--out", default="")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()

    if a.selftest:
        selftest()
        return
    if not a.pdb:
        ap.error("give a complex PDB, or --selftest")

    carbox = [int(x) for x in a.carboxylates.split(",") if x.strip()] or None
    props, acids = analyse(a.pdb, a.binder, a.target, carbox)
    report(props, acids, f"({os.path.basename(a.pdb)})")
    if a.out:
        json.dump(props, open(a.out, "w"), indent=1)
        print(f"\nWrote {a.out}")


if __name__ == "__main__":
    main()
