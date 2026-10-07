# === CHALLENGE 1 : ROUND 2 GENERATION — two GPUs, deeper optimisation, honest scoring ===
#
# Why round 1 failed (results/ROUND2_DIAGNOSIS.md): design_3stage(50,25,5) is ~a third
# of the ColabDesign default. With only 5 hard-stage iterations the sequence never
# became discrete, so the i_ptm AfDesign logged was the optimiser's own objective on a
# partly-continuous sequence - 0.725 for egfr_a_15, which re-predicts at 0.100.
#
# Three changes:
#   1. Two worker processes, one per T4. Kaggle bills session wall-clock, not GPU-hours,
#      so the second GPU is free throughput.
#   2. 150 iterations with a 20-iteration hard stage (round 1 used 5, which is why
#      its sequences never became discrete).
#   3. Every design is scored by an independent discrete re-prediction at 3 recycles
#      before it is written. Nothing enters the pool on the optimiser's own number again.

import os, subprocess, sys, textwrap, time

TIME_BUDGET_MIN = 360        # 6 h; Kaggle batch limit is 12 h, weekly quota is the real cap
LEN_MIN, LEN_MAX = 52, 72    # compact minibinders; still inside the 40-100 aa category
SOFT, TEMP, HARD = 90, 40, 20
TARGET = "/kaggle/working/targets/EGFR_A_domainI.pdb"
HOTSPOT = "20,21,22,25,28,50,51"

WORKER = r'''
import json, os, random, sys, time
GPU, SHARD, NSHARD, BUDGET = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
os.environ["CUDA_VISIBLE_DEVICES"] = GPU          # must precede the jax import
SOFT, TEMP, HARD = int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7])
LEN_MIN, LEN_MAX = int(sys.argv[8]), int(sys.argv[9])
TARGET, HOTSPOT = sys.argv[10], sys.argv[11]

t0 = time.time()
def log(m): print("[w%d] %s [%.1f min]" % (SHARD, m, (time.time()-t0)/60), flush=True)

from colabdesign import mk_afdesign_model, clear_mem

os.makedirs("/kaggle/working/out2", exist_ok=True)
rec_path = "/kaggle/working/out2/worker%d.jsonl" % SHARD
rng = random.Random(1000 + SHARD)

def iptm_of(model):
    lg = model.aux["log"]
    return float(lg.get("i_ptm", lg.get("ptm", 0.0)))

printed_weights = False
i = 0
while (time.time() - t0) / 60 < BUDGET:
    i += 1
    name = "r2_w%d_%03d" % (SHARD, i)
    L = rng.randint(LEN_MIN, LEN_MAX)
    try:
        clear_mem()
        model = mk_afdesign_model(protocol="binder", data_dir="/kaggle/working")
        model.prep_inputs(pdb_filename=TARGET, chain="A", binder_len=L, hotspot=HOTSPOT)

        if not printed_weights:
            log("default loss weights: %s" % dict(model.opt.get("weights", {})))
            printed_weights = True
        # Weights are left at ColabDesign's binder defaults. Round 1 showed the
        # optimiser does drive its own objective up; depth, not weighting, is the
        # variable under test, and changing both at once would confound the result.

        # design_3stage is the only staged entry point that exists in this
        # ColabDesign build - there is no design_temp, which killed the one
        # trajectory of the first attempt that passed the abort gate.
        mid = float("nan")
        model.design_3stage(soft_iters=SOFT, temp_iters=TEMP, hard_iters=HARD)

        seq = model.get_seqs()[0]
        design_iptm = iptm_of(model)

        # --- the honest number: re-predict the discrete sequence at 3 recycles ---
        try:
            model.predict(seq=seq, num_recycles=3, verbose=False)
        except TypeError:
            model.predict(seq=seq, verbose=False)
        lg = model.aux["log"]
        true_iptm = float(lg.get("i_ptm", lg.get("ptm", 0.0)))
        true_plddt = float(lg.get("plddt", 0.0))

        model.save_pdb("/kaggle/working/out2/%s.pdb" % name)
        with open(rec_path, "a") as fh:
            fh.write(json.dumps(dict(name=name, length=L, sequence=seq,
                                     soft_iptm=mid, design_iptm=round(design_iptm, 3),
                                     true_iptm=round(true_iptm, 3),
                                     plddt=round(true_plddt, 3))) + "\n")
        log("%s len=%d design_iptm=%.3f TRUE_iptm=%.3f plddt=%.3f"
            % (name, L, design_iptm, true_iptm, true_plddt))
    except Exception as e:
        log("%s FAILED %s %s" % (name, type(e).__name__, str(e)[:200]))
log("worker done after %d trajectories" % i)
'''

os.makedirs("/kaggle/working/out2", exist_ok=True)
with open("/kaggle/working/r2_worker.py", "w") as fh:
    fh.write(textwrap.dedent(WORKER))

args = [str(SOFT), str(TEMP), str(HARD), str(LEN_MIN), str(LEN_MAX),
        TARGET, HOTSPOT]
procs = []
for shard, gpu in enumerate(["0", "1"]):
    p = subprocess.Popen([sys.executable, "-u", "/kaggle/working/r2_worker.py",
                          gpu, str(shard), "2", str(TIME_BUDGET_MIN)] + args,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                         bufsize=1)
    procs.append(p)
    print("launched worker %d on GPU %s" % (shard, gpu), flush=True)

# Interleave both workers' output so the log stays readable and progress is visible.
import threading
def pump(p, tag):
    for line in p.stdout:
        print(line.rstrip(), flush=True)
threads = [threading.Thread(target=pump, args=(p, i), daemon=True)
           for i, p in enumerate(procs)]
[t.start() for t in threads]
for p in procs:
    p.wait()
[t.join(timeout=10) for t in threads]

import glob, json
recs = []
for f in sorted(glob.glob("/kaggle/working/out2/worker*.jsonl")):
    for line in open(f):
        line = line.strip()
        if line:
            recs.append(json.loads(line))
recs.sort(key=lambda r: -r["true_iptm"])
json.dump(recs, open("/kaggle/working/round2.json", "w"), indent=1)
print("\n==== ROUND 2: %d designs ====" % len(recs))
for r in recs[:25]:
    print("%-12s len=%2d TRUE_iptm=%.3f (design said %.3f) plddt=%.3f"
          % (r["name"], r["length"], r["true_iptm"], r["design_iptm"], r["plddt"]))
good = [r for r in recs if r["true_iptm"] >= 0.5]
ok = [r for r in recs if r["true_iptm"] >= 0.35]
print("\n%d at i_ptm>=0.50, %d at i_ptm>=0.35, best %.3f"
      % (len(good), len(ok), recs[0]["true_iptm"] if recs else 0.0))
