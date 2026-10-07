# === CHALLENGE 1 : BACKBONE + SEQUENCE GENERATION (AfDesign binder hallucination) ===
# Replaces RFdiffusion. Same target, same epitope hotspots, same length range.
NUM = 40         # aspirational; the time guard below decides the real count
TIME_BUDGET_MIN = 420   # hard stop at 7h so the 12h batch ALWAYS completes and saves
LEN_MIN, LEN_MAX = 55, 95      # keeps every design in the minibinder category (40-100)

import os, glob, random, time
import numpy as np
os.makedirs("/kaggle/working/out", exist_ok=True)
t0 = time.time()

TARGET = "/kaggle/working/targets/EGFR_A_domainI.pdb"
# Epitope A: F44,E45,D46,L49,Q52,Y74,D75 in UniProt = these in PDB numbering
HOTSPOT = "20,21,22,25,28,50,51"

from colabdesign import mk_afdesign_model, clear_mem

random.seed(37)
made = []
for i in range(NUM):
    if (time.time() - t0) / 60 > TIME_BUDGET_MIN:
        print("time budget reached; stopping with %d designs" % len(made), flush=True)
        break
    L = random.randint(LEN_MIN, LEN_MAX)
    name = "egfr_a_%d" % i
    out_pdb = "/kaggle/working/out/%s.pdb" % name
    if os.path.exists(out_pdb):
        made.append(out_pdb); continue
    try:
        clear_mem()
        model = mk_afdesign_model(protocol="binder", data_dir="/kaggle/working")
        model.prep_inputs(pdb_filename=TARGET, chain="A",
                          binder_len=L, hotspot=HOTSPOT)
        model.design_3stage(soft_iters=50, temp_iters=25, hard_iters=5)
        model.save_pdb(out_pdb)
        seq = model.get_seqs()[0]
        loss = float(model.aux["log"].get("loss", float("nan")))
        plddt = float(model.aux["log"].get("plddt", float("nan")))
        ptm = float(model.aux["log"].get("i_ptm", model.aux["log"].get("ptm", float("nan"))))
        print("%s len=%d loss=%.3f plddt=%.3f i_ptm=%.3f  [%.1f min]"
              % (name, L, loss, plddt, ptm, (time.time() - t0) / 60), flush=True)
        with open("/kaggle/working/out/%s.seq" % name, "w") as fh:
            fh.write(seq)
        made.append(out_pdb)
    except Exception as e:
        print("FAILED", name, type(e).__name__, str(e)[:300], flush=True)

print("\ngenerated %d/%d backbones in %.1f min"
      % (len(made), NUM, (time.time() - t0) / 60))
if made:
    # sanity: the saved complex must contain target + binder as separate chains
    from collections import Counter
    c = Counter()
    for line in open(made[0]):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            c[line[21]] += 1
    print("chains in first design:", dict(c))
