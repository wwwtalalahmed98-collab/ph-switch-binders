"""One-shot processor for the Kaggle output zip.

Unpacks challenge1_outputs.zip, locates the design PDBs, the FASTA and (optionally)
selfconsistency.json, then runs the Phase 3 switch analysis and writes designs.json
ready for the novelty screen and submission assembler.

Usage:
    python scripts/process_kaggle_output.py incoming/challenge1_outputs.zip \\
        --selfconsistency incoming/selfconsistency.json --binder-chain B
"""
import argparse
import glob
import json
import os
import subprocess
import sys
import zipfile
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def unpack(zip_path, dest):
    os.makedirs(dest, exist_ok=True)
    with zipfile.ZipFile(zip_path) as z:
        names = z.namelist()
        z.extractall(dest)
    print(f"unpacked {len(names)} entries -> {dest}")
    return names


def find(dest, pattern):
    return sorted(glob.glob(os.path.join(dest, "**", pattern), recursive=True))


def detect_binder_chain(pdb, lo=10, hi=250):
    """For a two-chain design complex the target is simply the larger chain and the
    binder the smaller. Matching on an expected target length is brittle: if the
    construct is trimmed differently both chains can fall in the binder range."""
    c = Counter()
    for line in open(pdb):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            c[line[21]] += 1
    if len(c) < 2:
        return None, None, dict(c)
    ordered = sorted(c.items(), key=lambda kv: -kv[1])
    target, binder = ordered[0][0], ordered[-1][0]
    if not (lo <= c[binder] <= hi):
        return None, None, dict(c)
    if len(c) > 2:
        print(f"WARNING: {len(c)} chains present {dict(c)}; "
              f"assuming largest is target and smallest is binder")
    return binder, target, dict(c)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("zip_path")
    ap.add_argument("--selfconsistency", default="")
    ap.add_argument("--binder-chain", default="")
    ap.add_argument("--target-chain", default="")
    ap.add_argument("--epitope", default="A")
    ap.add_argument("--generator", default="RFdiffusion")
    ap.add_argument("--dest", default=os.path.join(ROOT, "incoming", "unpacked"))
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "designs.json"))
    a = ap.parse_args()

    if not os.path.exists(a.zip_path):
        sys.exit(f"not found: {a.zip_path}")

    unpack(a.zip_path, a.dest)
    pdbs = [p for p in find(a.dest, "*.pdb") if "traj" not in p.lower()]
    fastas = find(a.dest, "*.fasta") + find(a.dest, "*.fa")
    print(f"found {len(pdbs)} PDBs, {len(fastas)} FASTA file(s)")
    if not pdbs:
        sys.exit("no PDBs in the zip")
    if not fastas:
        sys.exit("no FASTA in the zip")

    pdb_dir = os.path.dirname(pdbs[0])
    fasta = max(fastas, key=os.path.getsize)
    print(f"pdb dir: {pdb_dir}\nfasta:   {fasta}")

    bc, tc, counts = detect_binder_chain(pdbs[0])
    print(f"chains in {os.path.basename(pdbs[0])}: {counts}")
    binder = a.binder_chain or bc
    target = a.target_chain or tc
    if not binder or not target:
        sys.exit(f"could not identify chains automatically (got binder={bc}, "
                 f"target={tc}); pass --binder-chain / --target-chain")
    print(f"using binder chain {binder}, target chain {target}")

    cmd = [sys.executable, os.path.join(ROOT, "scripts", "phase34_pipeline.py"),
           "--pdbs", pdb_dir, "--fasta", fasta,
           "--epitope", a.epitope, "--binder-chain", binder,
           "--target-chain", target, "--generator", a.generator,
           "--out", a.out]
    if a.selfconsistency and os.path.exists(a.selfconsistency):
        cmd += ["--selfconsistency", a.selfconsistency]
    print("\n" + " ".join(cmd) + "\n")
    subprocess.run(cmd, check=True)

    records = json.load(open(a.out))
    sw = [r for r in records if r["role"] == "switch"]
    print(f"\n=== ready for novelty screen ===")
    print(f"{len(records)} records, {len(sw)} switch designs")
    fa = os.path.splitext(a.out)[0] + "_for_blast.fasta"
    with open(fa, "w") as fh:
        for r in records:
            fh.write(f">{r['name']}\n{r['sequence']}\n")
    print(f"wrote {fa}")
    print(f"\n  python scripts/phase4_novelty.py {fa} --out results/novelty.json")
    print(f"  python scripts/phase5_submission.py {a.out} --out submission/")


if __name__ == "__main__":
    main()
