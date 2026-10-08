# === CHALLENGE 2 (TNF-alpha) : GENERATION, ROUND 2 ===
# Same shape as the EGFR round that worked: two workers, one per T4, a wall-clock guard
# so the batch always completes and saves, results appended to JSONL as each design
# finishes, and every design scored by an INDEPENDENT re-prediction of the discrete
# sequence at 3 recycles rather than by the optimiser's own number.
#
# Added here, because EGFR had to learn it the hard way: each design also reports which
# target residues it actually contacts, in UniProt numbering, mapped through the
# construct map. A confident interface is not necessarily the intended one, and finding
# that out at the end cost a round last time.
import os, subprocess, sys, textwrap, time

TIME_BUDGET_MIN = 165          # 2.75 h, leaving the rest of the quota for switch + confirm
LEN_MIN, LEN_MAX = 55, 80
SOFT, TEMP, HARD = 90, 40, 20
TARGET = "/kaggle/working/targets/TNF_site_crop.pdb"
HOTSPOT = "A86,A87,A88,A90,A91,A92,B29,B32,B33,B34,B35,B36"

WORKER = r'''
import json, os, random, sys, time
GPU, SHARD, BUDGET = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
os.environ["CUDA_VISIBLE_DEVICES"] = GPU          # must precede the jax import
SOFT, TEMP, HARD = int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
LEN_MIN, LEN_MAX = int(sys.argv[7]), int(sys.argv[8])
TARGET, HOTSPOT = sys.argv[9], sys.argv[10]

t0 = time.time()
def log(m): print("[w%d] %s [%.1f min]" % (SHARD, m, (time.time()-t0)/60), flush=True)

import numpy as np
from colabdesign import mk_afdesign_model, clear_mem

CMAP = json.load(open("/kaggle/working/construct_map.json"))["order"]
SITE = {95,105,106,107,108,109,110,111,112,155,159,161,162,163,164,165,166,167,168,
        201,203,207,220,221,223}

os.makedirs("/kaggle/working/tnf_out2", exist_ok=True)
rec_path = "/kaggle/working/tnf_out2/worker%d.jsonl" % SHARD
rng = random.Random(3000 + SHARD)   # NEW SEED: 2000+shard would regenerate round 1 exactly

def epitope_of(pdb_path):
    """Target residues within 5 A of the binder, in UniProt numbering.

    AfDesign concatenates the target chains into one chain and renumbers from 1, so
    output residue n is CMAP[n-1]. Mapping through the map rather than assuming the
    numbering is the whole point.
    """
    from Bio.PDB import PDBParser
    from collections import Counter
    mdl = PDBParser(QUIET=True).get_structure("c", pdb_path)[0]
    sizes = Counter({ch.id: sum(1 for r in ch if r.id[0]==" ") for ch in mdl})
    tgt = max(sizes, key=lambda k: sizes[k]); bnd = min(sizes, key=lambda k: sizes[k])
    tres = [r for r in mdl[tgt] if r.id[0]==" "]
    # the file carries two models, so the chain holds each residue twice
    n = len(CMAP)
    if len(tres) >= 2*n: tres = tres[:n]
    if len(tres) != n:
        return None, "target chain %d residues, map has %d" % (len(tres), n)
    ba = np.array([a.get_coord() for r in mdl[bnd] if r.id[0]==" " for a in r])
    hits = []
    for i, r in enumerate(tres):
        ra = np.array([a.get_coord() for a in r])
        if np.min(np.linalg.norm(ra[:,None,:]-ba[None,:,:], axis=-1)) <= 5.0:
            hits.append(CMAP[i]["uniprot"])
    return sorted(set(hits)), None

printed = False
i = 0
while (time.time() - t0) / 60 < BUDGET:
    i += 1
    name = "tnf2_w%d_%03d" % (SHARD, i)
    L = rng.randint(LEN_MIN, LEN_MAX)
    try:
        clear_mem()
        model = mk_afdesign_model(protocol="binder", data_dir="/kaggle/working")
        model.prep_inputs(pdb_filename=TARGET, chain="A,B", binder_len=L, hotspot=HOTSPOT)
        if not printed:
            log("target length %d, loss weights %s"
                % (model._target_len, dict(model.opt.get("weights", {}))))
            printed = True

        model.design_3stage(soft_iters=SOFT, temp_iters=TEMP, hard_iters=HARD)
        seq = model.get_seqs()[0]
        lg = model.aux["log"]
        design_iptm = float(lg.get("i_ptm", lg.get("ptm", 0.0)))

        try:
            model.predict(seq=seq, num_recycles=3, verbose=False)
        except TypeError:
            model.predict(seq=seq, verbose=False)
        lg = model.aux["log"]
        true_iptm = float(lg.get("i_ptm", lg.get("ptm", 0.0)))
        plddt = float(lg.get("plddt", 0.0))

        out = "/kaggle/working/tnf_out2/%s.pdb" % name
        model.save_pdb(out)
        epi, err = epitope_of(out)
        on_site = sorted(set(epi) & SITE) if epi else []

        with open(rec_path, "a") as fh:
            fh.write(json.dumps(dict(name=name, length=L, sequence=seq,
                                     design_iptm=round(design_iptm,3),
                                     true_iptm=round(true_iptm,3),
                                     plddt=round(plddt,3),
                                     epitope=epi, on_site=on_site,
                                     epitope_error=err)) + "\n")
        log("%s len=%d design_iptm=%.3f TRUE_iptm=%.3f plddt=%.3f  epitope %s, %d on site %s"
            % (name, L, design_iptm, true_iptm, plddt,
               ("%d res" % len(epi)) if epi else ("ERR " + str(err)),
               len(on_site), on_site[:8]))
    except Exception as e:
        log("%s FAILED %s %s" % (name, type(e).__name__, str(e)[:200]))
log("worker done after %d trajectories" % i)
'''

os.makedirs("/kaggle/working/tnf_out2", exist_ok=True)
with open("/kaggle/working/tnf_worker2.py", "w") as fh:
    fh.write(textwrap.dedent(WORKER))

args = [str(SOFT), str(TEMP), str(HARD), str(LEN_MIN), str(LEN_MAX), TARGET, HOTSPOT]
procs = []
for shard, gpu in enumerate(["0", "1"]):
    p = subprocess.Popen([sys.executable, "-u", "/kaggle/working/tnf_worker2.py",
                          gpu, str(shard), str(TIME_BUDGET_MIN)] + args,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         text=True, bufsize=1)
    procs.append(p)
    print("launched worker %d on GPU %s" % (shard, gpu), flush=True)

import threading
def pump(p):
    for line in p.stdout:
        print(line.rstrip(), flush=True)
threads = [threading.Thread(target=pump, args=(p,), daemon=True) for p in procs]
[t.start() for t in threads]
for p in procs:
    p.wait()
[t.join(timeout=10) for t in threads]

import glob, json
recs = []
for f in sorted(glob.glob("/kaggle/working/tnf_out2/worker*.jsonl")):
    for line in open(f):
        line = line.strip()
        if line:
            recs.append(json.loads(line))
recs.sort(key=lambda r: -r["true_iptm"])
json.dump(recs, open("/kaggle/working/tnf_round2.json", "w"), indent=1)

print("\n==== TNF GENERATION: %d designs ====" % len(recs))
print("%-13s %4s %7s %8s %7s %7s" % ("name","len","i_ptm","(design)","pLDDT","on site"))
for r in recs[:30]:
    print("%-13s %4d %7.3f %8.3f %7.3f %7d"
          % (r["name"], r["length"], r["true_iptm"], r["design_iptm"],
             r["plddt"], len(r["on_site"])))
good = [r for r in recs if r["true_iptm"] >= 0.50]
ok = [r for r in recs if r["true_iptm"] >= 0.35]
onsite = [r for r in recs if len(r["on_site"]) >= 5]
print("\n%d at i_ptm>=0.50, %d at >=0.35, best %.3f; %d contact 5+ site residues"
      % (len(good), len(ok), recs[0]["true_iptm"] if recs else 0.0, len(onsite)))
