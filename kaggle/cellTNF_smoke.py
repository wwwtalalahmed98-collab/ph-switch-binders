# === CHALLENGE 2 (TNF-alpha) : GENERATION SMOKE TEST ===
# Two short trajectories, not a production run. It answers four questions before six
# GPU-hours are committed to them:
#   1. does AfDesign accept a two-chain target and chain-prefixed hotspots at all?
#   2. does the saved complex come back with three chains (A + B target, binder)?
#   3. does independent re-prediction work on a two-chain target?
#   4. how long does one iteration actually take at 304 + 70 residues, so the real
#      run can be sized from measurement instead of from a scaling guess?
import json, os, random, time
import numpy as np

NUM = 2
SOFT, TEMP, HARD = 20, 10, 5      # deliberately short: mechanics, not quality
LEN_MIN, LEN_MAX = 65, 75
TARGET = "/kaggle/working/targets/TNF_AB.pdb"
HOTSPOT = "A86,A87,A88,A90,A91,A92,B29,B32,B33,B34,B35,B36"

t0 = time.time()
os.makedirs("/kaggle/working/tnf_smoke", exist_ok=True)
from colabdesign import mk_afdesign_model, clear_mem

rng = random.Random(7)
rows = []
for i in range(NUM):
    name = "tnf_smoke_%d" % i
    L = rng.randint(LEN_MIN, LEN_MAX)
    try:
        t_prep = time.time()
        clear_mem()
        model = mk_afdesign_model(protocol="binder", data_dir="/kaggle/working")
        model.prep_inputs(pdb_filename=TARGET, chain="A,B", binder_len=L, hotspot=HOTSPOT)
        prep_s = time.time() - t_prep
        if i == 0:
            print("default loss weights:", dict(model.opt.get("weights", {})), flush=True)

        t_des = time.time()
        model.design_3stage(soft_iters=SOFT, temp_iters=TEMP, hard_iters=HARD)
        des_s = time.time() - t_des
        seq = model.get_seqs()[0]
        lg = model.aux["log"]
        design_iptm = float(lg.get("i_ptm", lg.get("ptm", 0.0)))

        t_pred = time.time()
        try:
            model.predict(seq=seq, num_recycles=3, verbose=False)
        except TypeError:
            model.predict(seq=seq, verbose=False)
        pred_s = time.time() - t_pred
        lg = model.aux["log"]
        true_iptm = float(lg.get("i_ptm", lg.get("ptm", 0.0)))
        plddt = float(lg.get("plddt", 0.0))

        out = "/kaggle/working/tnf_smoke/%s.pdb" % name
        model.save_pdb(out)
        ch = {}
        for line in open(out):
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                ch[line[21]] = ch.get(line[21], 0) + 1

        iters = SOFT + TEMP + HARD
        print("%s len=%d  prep %.0fs  design %.0fs (%d iters, %.1f s/iter)  "
              "predict3rc %.0fs  design_iptm=%.3f TRUE_iptm=%.3f plddt=%.3f"
              % (name, L, prep_s, des_s, iters, des_s / iters, pred_s,
                 design_iptm, true_iptm, plddt), flush=True)
        print("    chains in saved complex: %s" % ch, flush=True)
        rows.append(dict(name=name, length=L, s_per_iter=des_s / iters,
                         predict_s=pred_s, design_iptm=round(design_iptm, 3),
                         true_iptm=round(true_iptm, 3), plddt=round(plddt, 3),
                         chains=ch, sequence=seq))
    except Exception as e:
        print("FAILED", name, type(e).__name__, str(e)[:400], flush=True)

json.dump(rows, open("/kaggle/working/tnf_smoke.json", "w"), indent=1)

print("\n==== SMOKE TEST VERDICT ====")
if not rows:
    print("no trajectory completed - two-chain targets are not working as written")
else:
    spi = sum(r["s_per_iter"] for r in rows) / len(rows)
    pre = sum(r["predict_s"] for r in rows) / len(rows)
    chains_ok = all(len(r["chains"]) == 3 for r in rows)
    print("completed %d/%d trajectories" % (len(rows), NUM))
    print("saved complex has 3 chains (A + B target, binder): %s" % chains_ok)
    print("%.1f s per design iteration, %.0f s per 3-recycle prediction" % (spi, pre))
    for iters in (150, 120, 90):
        per = (iters * spi + pre + 40) / 3600.0
        print("   at %3d iterations: %.2f h per design -> %d designs in 6 h on 2 GPUs"
              % (iters, per, int(2 * 6 / per)))
    print("\nChallenge 1 ran at ~6 s/iter on a 161+60 target; anything near 16 s/iter")
    print("here is the expected cost of 304+70 residues, not a fault.")
