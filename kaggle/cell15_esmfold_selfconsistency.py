# Add as a NEW cell after cell 14 (ProteinMPNN), before the collection cell.
#
# Self-consistency / designability test: fold each designed sequence with ESMFold and
# compare it to the RFdiffusion backbone it was designed onto. This is the standard
# quality filter for de novo binders, and without it we would have no signal separating
# a well-folded design from a bad one.
#
#   scRMSD < 2.0 A and pLDDT > 80  ->  designable
#
# Binders are 55-95 aa so this is fast on a P100 (a few seconds each).

import glob
import json
import os
from collections import Counter

import numpy as np
import torch
from transformers import AutoTokenizer, EsmForProteinFolding

OUT_JSON = "/kaggle/working/selfconsistency.json"
PLDDT_MIN, RMSD_MAX = 80.0, 2.0


def ca_coords(path, chain=None):
    xyz = []
    for line in open(path):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            if chain and line[21] != chain:
                continue
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(xyz)


def kabsch_rmsd(P, Q):
    """RMSD after optimal superposition."""
    if len(P) != len(Q) or len(P) == 0:
        return float("nan")
    P = P - P.mean(0)
    Q = Q - Q.mean(0)
    V, S, W = np.linalg.svd(P.T @ Q)
    if (np.linalg.det(V) * np.linalg.det(W)) < 0:
        V[:, -1] = -V[:, -1]
    P = P @ (V @ W)
    return float(np.sqrt(((P - Q) ** 2).sum() / len(P)))


def binder_chain_of(path, target_len=161, lo=40, hi=110):
    c = Counter()
    for line in open(path):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            c[line[21]] += 1
    for ch, n in c.items():
        if lo <= n <= hi:
            return ch
    return None


print("loading ESMFold (first run downloads ~2.6 GB) ...")
tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1",
                                             low_cpu_mem_usage=True)
model = model.cuda().eval()
model.esm = model.esm.half()
torch.backends.cuda.matmul.allow_tf32 = True
print("ESMFold ready")

# Map each MPNN sequence back to its parent backbone.
records = []
for fa in sorted(glob.glob("/kaggle/working/seqs/seqs/*.fa")):
    stem = os.path.basename(fa).replace(".fa", "")
    parent = f"/kaggle/working/out/{stem}.pdb"
    if not os.path.exists(parent):
        print("no parent PDB for", stem)
        continue
    ch = binder_chain_of(parent)
    ref = ca_coords(parent, ch)

    entries, hdr = [], None
    for line in open(fa):
        line = line.strip()
        if line.startswith(">"):
            hdr = line
        elif line:
            entries.append((hdr, line))

    for i, (h, seq) in enumerate(entries[1:], start=1):   # entry 0 = input backbone
        name = f"{stem}_s{i}"
        with torch.no_grad():
            t = tok([seq], return_tensors="pt", add_special_tokens=False)
            t = {k: v.cuda() for k, v in t.items()}
            out = model(**t)
        pdb_str = model.output_to_pdb(out)[0]
        fold_path = f"/kaggle/working/folded/{name}.pdb"
        os.makedirs("/kaggle/working/folded", exist_ok=True)
        open(fold_path, "w").write(pdb_str)

        plddt = float(out["plddt"][0, :, 1].mean() * 100)
        rmsd = kabsch_rmsd(ca_coords(fold_path), ref)
        # ProteinMPNN puts its global score in the FASTA header
        mpnn = None
        for field in (h or "").split(","):
            if "global_score" in field:
                try:
                    mpnn = float(field.split("=")[1])
                except Exception:
                    pass
        passed = (plddt >= PLDDT_MIN) and (rmsd <= RMSD_MAX)
        records.append(dict(name=name, backbone=stem, sequence=seq,
                            length=len(seq), plddt=round(plddt, 1),
                            sc_rmsd=round(rmsd, 2), mpnn_score=mpnn,
                            designable=passed))
        print(f"{name}: pLDDT {plddt:5.1f}  scRMSD {rmsd:5.2f}  "
              f"{'PASS' if passed else 'fail'}")

json.dump(records, open(OUT_JSON, "w"), indent=1)
good = [r for r in records if r["designable"]]
print(f"\n{len(good)}/{len(records)} designable "
      f"(pLDDT >= {PLDDT_MIN}, scRMSD <= {RMSD_MAX}) -> {OUT_JSON}")
if records and not good:
    best = sorted(records, key=lambda r: (r["sc_rmsd"], -r["plddt"]))[:5]
    print("\nNothing passed. Best five, so the bar can be relaxed deliberately "
          "rather than silently:")
    for r in best:
        print(f"   {r['name']}: pLDDT {r['plddt']}, scRMSD {r['sc_rmsd']}")
