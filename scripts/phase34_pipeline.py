"""Phase 3+4 driver: raw design output -> submission-ready design records.

Takes a directory of RFdiffusion/BindCraft complex PDBs (target chain A, binder chain B)
plus the ProteinMPNN FASTA, and for each design:

  1. finds binder positions that could reach a conserved target carboxylate (Phase 3
     geometry, calibrated on real salt bridges),
  2. writes the acid-ON switch variant by substituting those positions to histidine,
  3. keeps the unmutated sequence as a pH-insensitive control at the same epitope,
  4. emits designs.json for scripts/phase5_submission.py.

Histidine placement is per-backbone (it depends on CB geometry), so every sequence from a
given backbone inherits that backbone's switch positions.

Usage:
    python scripts/phase34_pipeline.py --pdbs out/ --fasta designs.fasta \\
        --epitope A --out results/designs.json
"""
import argparse
import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from phase3_ph_switch import analyse  # noqa: E402

# Conserved target carboxylates per epitope, in PDB numbering (UniProt - 24).
EPITOPE_CARBOXYLATES = {
    "A": {"pdb": [21, 22, 51, 73], "uniprot": "E45,D46,D75,E97",
          "name": "domain I ligand-binding face"},
    "B": {"pdb": [254, 258, 297, 578], "uniprot": "D278,E282,D321,E602",
          "name": "domain II dimerisation arm"},
}
MAX_HIS = 4          # cap substitutions; more risks destroying affinity outright
MIN_CONTACTS = 2     # published threshold for a credible switch


def read_fasta(path):
    out, name, buf = {}, None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        elif line:
            buf.append(line)
    if name:
        out[name] = "".join(buf)
    return out


def backbone_key(seq_name):
    """ProteinMPNN names look like <backbone>_s3; strip the sample suffix."""
    return seq_name.rsplit("_s", 1)[0]


# Preference for what we substitute. ProteinMPNN places Lys/Arg against a carboxylate-rich
# epitope for electrostatic complementarity; converting those to His turns a pH-INDEPENDENT
# salt bridge into a pH-DEPENDENT one, which is precisely the acid-ON switch. Polar residues
# are neutral choices. Hydrophobics are penalised: substituting one costs interface packing
# and buys only the charge, and aromatics may be doing real binding work.
RESIDUE_PREFERENCE = {
    "K": 1.5, "R": 1.5,                                  # permanent -> pH-dependent
    "Q": 1.2, "N": 1.2, "S": 1.15, "T": 1.15, "A": 1.0,  # cheap to replace
    "D": 1.1, "E": 1.1,                                  # removes like-charge repulsion
    "H": 0.0,                                            # already histidine, no-op
    "M": 0.7, "L": 0.6, "I": 0.6, "V": 0.6,              # packing cost
    "F": 0.4, "Y": 0.4, "W": 0.3,                        # likely load-bearing
}


def pick_positions(proposals, sequence=None, offset=1):
    """Best histidine positions: highest scoring, one per binder position, preferring
    positions that reach distinct carboxylates and that are cheap to substitute."""
    def adjusted(p):
        s = p["score"]
        if sequence:
            i = p["binder_resnum"] - offset
            if 0 <= i < len(sequence):
                s *= RESIDUE_PREFERENCE.get(sequence[i], 0.9)
        return s

    chosen, used_pos, used_acid = [], set(), set()
    for p in sorted(proposals, key=lambda x: -adjusted(x)):
        if sequence and adjusted(p) <= 0:
            continue
        if p["binder_resnum"] in used_pos:
            continue
        if p["target_resnum"] in used_acid and len(chosen) >= 2:
            continue
        chosen.append(p)
        used_pos.add(p["binder_resnum"])
        used_acid.add(p["target_resnum"])
        if len(chosen) >= MAX_HIS:
            break
    return chosen


def apply_his(seq, positions, offset=1):
    """Substitute 1-based binder positions to histidine."""
    s = list(seq)
    applied = []
    for pos in positions:
        i = pos - offset
        if 0 <= i < len(s):
            if s[i] != "H":
                applied.append(f"{s[i]}{pos}H")
            s[i] = "H"
    return "".join(s), applied


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdbs", required=True, help="directory of complex PDBs")
    ap.add_argument("--fasta", required=True, help="ProteinMPNN sequences")
    ap.add_argument("--epitope", default="A", choices=sorted(EPITOPE_CARBOXYLATES))
    ap.add_argument("--binder-chain", default="B")
    ap.add_argument("--target-chain", default="A")
    ap.add_argument("--generator", default="RFdiffusion")
    ap.add_argument("--controls-per-backbone", type=int, default=1)
    ap.add_argument("--selfconsistency", default="",
                    help="selfconsistency.json from the ESMFold cell")
    ap.add_argument("--designable-only", action="store_true",
                    help="drop sequences that failed the pLDDT/scRMSD filter")
    ap.add_argument("--out", default="results/designs.json")
    a = ap.parse_args()

    ep = EPITOPE_CARBOXYLATES[a.epitope]
    seqs = read_fasta(a.fasta)

    # Quality metrics from the ESMFold designability test, keyed by sequence name.
    quality = {}
    if a.selfconsistency and os.path.exists(a.selfconsistency):
        for r in json.load(open(a.selfconsistency)):
            quality[r["name"]] = r
        n_ok = sum(1 for r in quality.values() if r.get("designable"))
        print(f"loaded quality metrics for {len(quality)} sequences "
              f"({n_ok} designable)")
        if a.designable_only:
            before = len(seqs)
            seqs = {k: v for k, v in seqs.items()
                    if quality.get(k, {}).get("designable")}
            print(f"designable filter: {before} -> {len(seqs)} sequences")
    elif a.selfconsistency:
        print(f"WARNING: {a.selfconsistency} not found - ranking will have no "
              f"structural quality signal")
    by_backbone = {}
    for name, s in seqs.items():
        by_backbone.setdefault(backbone_key(name), []).append((name, s))

    pdbs = sorted(glob.glob(os.path.join(a.pdbs, "*.pdb")))
    print(f"{len(pdbs)} complex PDBs, {len(seqs)} sequences across "
          f"{len(by_backbone)} backbones")
    print(f"epitope {a.epitope} ({ep['name']}); target carboxylates {ep['uniprot']}\n")

    records, skipped = [], []
    for pdb in pdbs:
        stem = os.path.splitext(os.path.basename(pdb))[0]
        matches = [k for k in by_backbone if k == stem or k.startswith(stem)]
        if not matches:
            continue
        try:
            proposals, _ = analyse(pdb, a.binder_chain, a.target_chain, ep["pdb"])
        except Exception as e:
            skipped.append(f"{stem}: {e}")
            continue

        # Residue preference depends on the sequence, so pick positions per backbone
        # using its first sequence as representative (MPNN samples agree closely at
        # interface positions at low sampling temperature).
        rep_seq = sorted(by_backbone[matches[0]])[0][1]
        chosen = pick_positions(proposals, rep_seq)
        n_acids = len({p["target_resnum"] for p in chosen})
        switch_score = round(sum(p["score"] for p in chosen), 3)
        positions = [p["binder_resnum"] for p in chosen]

        verdict = "switch" if n_acids >= MIN_CONTACTS else "control-only"
        print(f"{stem}: {len(proposals)} candidate contacts, "
              f"{n_acids} distinct carboxylates reachable -> {verdict}")

        for idx, (name, seq) in enumerate(sorted(by_backbone[matches[0]])):
            q = quality.get(name, {})
            metrics = dict(plddt=q.get("plddt"), sc_rmsd=q.get("sc_rmsd"),
                           mpnn_score=q.get("mpnn_score"),
                           designable=q.get("designable"))
            if n_acids >= MIN_CONTACTS:
                mutated, applied = apply_his(seq, positions)
                records.append(dict(
                    name=f"{name}_sw", sequence=mutated, role="switch",
                    epitope=f"{a.epitope}_{ep['name'].replace(' ', '_')}",
                    n_his_contacts=n_acids, switch_score=switch_score,
                    his_substitutions=",".join(applied) or "none",
                    epitope_conservation=1.0, generator=a.generator,
                    source_backbone=stem, **metrics))
            if idx < a.controls_per_backbone:
                records.append(dict(
                    name=f"{name}_ctl", sequence=seq, role="control",
                    epitope=f"{a.epitope}_{ep['name'].replace(' ', '_')}",
                    n_his_contacts=0, switch_score=0.0,
                    his_substitutions="none", epitope_conservation=1.0,
                    generator=a.generator, source_backbone=stem, **metrics))

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    json.dump(records, open(a.out, "w"), indent=1)

    n_sw = sum(1 for r in records if r["role"] == "switch")
    n_ct = sum(1 for r in records if r["role"] == "control")
    print(f"\n{len(records)} design records ({n_sw} switch, {n_ct} control) -> {a.out}")
    if skipped:
        print(f"{len(skipped)} PDB(s) skipped:")
        for s in skipped[:5]:
            print("   ", s)
    print("\nnext: novelty screen, then phase5_submission.py")
    print(f"  python scripts/phase4_novelty.py <fasta of final designs> "
          f"--out results/novelty.json")
    print(f"  python scripts/phase5_submission.py {a.out} --out submission/")


if __name__ == "__main__":
    main()
