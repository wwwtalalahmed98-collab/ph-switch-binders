# Drop-in replacement for cell 14 (ProteinMPNN). Paste over the existing cell.
# Auto-detects which chain is the diffused binder instead of assuming "B", which is
# the most common way this step fails.

import glob
import os
import subprocess
from collections import Counter

TARGET_LEN = 161          # EGFR domain I construct
BINDER_RANGE = (40, 110)  # designed binder, sampled 55-95

pdbs = sorted(glob.glob("/kaggle/working/out/egfr_a*.pdb"))
print(f"{len(pdbs)} backbones from RFdiffusion")
if not pdbs:
    raise SystemExit("No output PDBs - did cell 9 finish?")


def chain_lengths(path):
    seen = Counter()
    for line in open(path):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            seen[line[21]] += 1
    return seen


counts = chain_lengths(pdbs[0])
print("chains in", os.path.basename(pdbs[0]), "->", dict(counts))

binder = [c for c, n in counts.items() if BINDER_RANGE[0] <= n <= BINDER_RANGE[1]]
target = [c for c, n in counts.items() if abs(n - TARGET_LEN) < 15]
if not binder:
    raise SystemExit(f"Could not identify the binder chain from {dict(counts)} - "
                     f"check the output PDB manually.")
BINDER_CHAIN = binder[0]
print(f"binder chain = {BINDER_CHAIN} ({counts[BINDER_CHAIN]} residues); "
      f"target chain = {target[0] if target else '?'}")

os.makedirs("/kaggle/working/seqs", exist_ok=True)
ok = 0
for p in pdbs:
    cmd = ["python", "/kaggle/working/ProteinMPNN/protein_mpnn_run.py",
           "--pdb_path", p, "--pdb_path_chains", BINDER_CHAIN,
           "--out_folder", "/kaggle/working/seqs",
           "--num_seq_per_target", "8", "--sampling_temp", "0.1",
           "--omit_AAs", "CX", "--use_soluble_model", "--seed", "37"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(os.path.basename(p), "FAILED:")
        print(r.stderr[-600:])
        break
    ok += 1

print(f"\nsequence design done for {ok}/{len(pdbs)} backbones")
print("NOTE: record the binder chain id -", BINDER_CHAIN,
      "- you will need it for the local pipeline (--binder-chain).")
