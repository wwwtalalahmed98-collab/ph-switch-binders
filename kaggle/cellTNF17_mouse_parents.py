# === CHALLENGE 2 (TNF-alpha) : MOUSE SCREEN OF THE PARENT DESIGNS ===
#
# The measured mouse run showed that epitope conservation does not predict mouse binding
# (r = +0.272 over the 11 designs that had both numbers) and that six of the ten
# submitted switches fail objective 2, which is ranked above affinity. Only three
# unsubmitted switches clear both bars, so the usable pool is too small.
#
# Every switch so far was chosen for human binding and checked against mouse afterwards.
# That is backwards. Mouse cross-reactivity is a property of the parent binder's epitope,
# not of the two histidines placed on top of it. This screens the PARENT designs - all of
# them, from all three rounds, before any histidine is placed - against mouse, so the
# switch stage can be spent only on lineages that already bind both species.
#
# Nothing is designed here and nothing is re-scored against human; the human numbers are
# read from the round files as they were recorded.
import glob, json, os, subprocess, sys, textwrap, time

TIME_BUDGET_MIN = 150
MIN_IPTM = 0.50        # the same human bar the switch stages used
MIN_ONSITE = 5

WORKER = r"""
import glob, json, os, sys, time
from collections import Counter
GPU, SHARD, BUDGET = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
MIN_IPTM, MIN_ONSITE = float(sys.argv[4]), int(sys.argv[5])
os.environ["CUDA_VISIBLE_DEVICES"] = GPU
t0 = time.time()
def log(m): print("[w%d] %s [%.1f min]" % (SHARD, m, (time.time()-t0)/60), flush=True)

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1

TARGET = "/kaggle/working/targets/TNF_mouse_site_crop.pdb"
CM = json.load(open("/kaggle/working/construct_map_mouse.json"))
CMAP = CM["order"]
N = len(CMAP)
MSITE = set(CM["mouse_site"])
M2H = {int(k): v for k, v in CM["mouse_to_human_site"].items()}
HOTSPOT = open("/kaggle/working/mouse_hotspot.txt").read().strip()

# The smoke test recorded every chain written twice. A chain overlapping itself reads as
# fully buried and would silently change every contact count, so this asserts the size.
def trim(chain, n, what):
    res = [r for r in chain if r.id[0] == " "]
    if len(res) > n and len(res) % n == 0:
        for r in res[n:]:
            chain.detach_child(r.id)
        res = res[:n]
    if len(res) != n:
        raise ValueError("%s chain has %d residues, expected %d" % (what, len(res), n))
    return res

def epitope(pdb_path, blen):
    mdl = PDBParser(QUIET=True).get_structure("c", pdb_path)[0]
    sizes = Counter({ch.id: sum(1 for r in ch if r.id[0]==" ") for ch in mdl})
    tgt, bnd = max(sizes, key=lambda k: sizes[k]), min(sizes, key=lambda k: sizes[k])
    tres = trim(mdl[tgt], N, "target")
    bres = trim(mdl[bnd], blen, "binder")
    seen = "".join(protein_letters_3to1.get(r.get_resname(), "X") for r in tres)
    if seen != "".join(o["aa"] for o in CMAP):
        raise ValueError("mouse target sequence does not match construct_map_mouse")
    batoms = np.array([a.get_coord() for r in bres for a in r])
    hits = []
    for i, r in enumerate(tres):
        ra = np.array([a.get_coord() for a in r])
        if np.min(np.linalg.norm(ra[:,None,:]-batoms[None,:,:], axis=-1)) <= 5.0:
            hits.append(CMAP[i]["uniprot"])
    return sorted(set(hits))

from colabdesign import mk_afdesign_model, clear_mem
_c = {}
def score(seq, out_pdb):
    if _c.get("L") != len(seq):
        clear_mem()
        af = mk_afdesign_model(protocol="binder", data_dir="/kaggle/working")
        af.prep_inputs(pdb_filename=TARGET, chain="A,B", binder_len=len(seq),
                       hotspot=HOTSPOT)
        _c["L"], _c["af"] = len(seq), af
    af = _c["af"]
    try:
        af.predict(seq=seq, num_recycles=3, verbose=False)
    except TypeError:
        af.predict(seq=seq, verbose=False)
    lg = af.aux["log"]
    af.save_pdb(out_pdb)
    return (round(float(lg.get("i_ptm", lg.get("ptm", 0.0))), 3),
            round(float(lg.get("plddt", 0.0)), 3))

# every parent design from every round file that is attached
recs = {}
for f in sorted(glob.glob("/kaggle/input/**/tnf_round*.json", recursive=True)):
    rnd = os.path.basename(f).replace("tnf_round", "").replace(".json", "")
    for r in json.load(open(f)):
        r["round"] = rnd
        recs.setdefault(r["name"], r)
if not recs:
    raise SystemExit("no tnf_round*.json found - attach the tnf_generate outputs")
kept = [r for r in recs.values()
        if r.get("true_iptm", 0) >= MIN_IPTM and len(r.get("on_site", [])) >= MIN_ONSITE]
# Sorting by LENGTH, not by score: each new binder length forces a model rebuild and
# recompile, which dominates the cost when only predicting.
kept.sort(key=lambda r: (r["length"], r["name"]))
mine = [r for i, r in enumerate(kept) if i % 2 == SHARD]
if SHARD == 0:
    log("%d parent designs across attached rounds, %d pass human i_ptm >= %.2f with %d+ on-site"
        % (len(recs), len(kept), MIN_IPTM, MIN_ONSITE))
log("this worker screens %d against mouse TNF-alpha" % len(mine))

os.makedirs("/kaggle/working/mouse_parent_pdb", exist_ok=True)
rec_path = "/kaggle/working/mouse_parent_w%d.jsonl" % SHARD
for d in mine:
    if (time.time()-t0)/60 > BUDGET:
        log("time budget reached"); break
    try:
        out = "/kaggle/working/mouse_parent_pdb/%s.pdb" % d["name"]
        ip, pl = score(d["sequence"], out)
        epi = epitope(out, d["length"])
        onsite = sorted(set(epi) & MSITE)
        hu = d["true_iptm"]
        row = dict(name=d["name"], round=d["round"], length=d["length"],
                   human_i_ptm=hu, human_on_site=len(d.get("on_site", [])),
                   mouse_i_ptm=ip, mouse_plddt=pl,
                   mouse_epitope=epi, mouse_on_site=onsite,
                   on_site_human_equiv=sorted(M2H[m] for m in onsite if m in M2H),
                   ratio=round(ip / hu, 3) if hu else None)
        with open(rec_path, "a") as fh:
            fh.write(json.dumps(row) + chr(10))
        log("%-16s human %.3f -> mouse %.3f (%3.0f%%)  epitope %2d, %2d on mouse site"
            % (d["name"], hu, ip, 100*(ip/hu if hu else 0), len(epi), len(onsite)))
    except Exception as e:
        log("%s FAILED %s %s" % (d["name"], type(e).__name__, str(e)[:200]))
log("worker done")
"""

with open("/kaggle/working/mouse_parent_worker.py", "w") as fh:
    fh.write(textwrap.dedent(WORKER))

procs = []
for shard, gpu in enumerate(["0", "1"]):
    p = subprocess.Popen([sys.executable, "-u", "/kaggle/working/mouse_parent_worker.py",
                          gpu, str(shard), str(TIME_BUDGET_MIN),
                          str(MIN_IPTM), str(MIN_ONSITE)],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         text=True, bufsize=1)
    procs.append(p)
    print("launched parent-screen worker %d on GPU %s" % (shard, gpu), flush=True)

import threading
def pump(p):
    for line in p.stdout:
        print(line.rstrip(), flush=True)
threads = [threading.Thread(target=pump, args=(p,), daemon=True) for p in procs]
[t.start() for t in threads]
for p in procs:
    p.wait()
[t.join(timeout=10) for t in threads]

rows = []
for f in sorted(glob.glob("/kaggle/working/mouse_parent_w*.jsonl")):
    for line in open(f):
        line = line.strip()
        if line:
            rows.append(json.loads(line))
rows.sort(key=lambda r: -r["mouse_i_ptm"])
json.dump(rows, open("/kaggle/working/tnf_mouse_parents.json", "w"), indent=1)

good = [r for r in rows if r["mouse_i_ptm"] >= 0.45 and len(r["mouse_on_site"]) >= 5]
print("\n==== PARENT DESIGNS AGAINST MOUSE ====")
print("%-18s %3s %4s %8s %8s %7s %8s"
      % ("design", "rnd", "len", "human", "mouse", "ratio", "on site"))
for r in rows:
    print("%-18s %3s %4d %8.3f %8.3f %6.0f%% %8d"
          % (r["name"], r["round"], r["length"], r["human_i_ptm"], r["mouse_i_ptm"],
             100*(r["ratio"] or 0), len(r["mouse_on_site"])))
print("\n%d parents screened; %d bind mouse at i_ptm >= 0.45 with 5+ on-site contacts"
      % (len(rows), len(good)))
print("these are the lineages worth spending the switch stage on:")
for r in good:
    print("  %-18s human %.3f  mouse %.3f  (%.0f%%)  %d mouse-site contacts"
          % (r["name"], r["human_i_ptm"], r["mouse_i_ptm"], 100*r["ratio"],
             len(r["mouse_on_site"])))
print("\n===== COPY FROM HERE =====")
for r in rows:
    print(json.dumps(r, separators=(",", ":")))
print("===== COPY TO HERE =====")
