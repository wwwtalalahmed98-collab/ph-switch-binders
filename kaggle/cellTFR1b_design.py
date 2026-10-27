"""Kaggle cell 2 of 2 - TfR1 Challenge 4. Generate binders at two epitopes, one per GPU.

Run cellTFR1a_setup.py first; this reads the targets it wrote and will refuse to start
without them.

WHAT THIS DOES AND DOES NOT DO. This generates and scores BINDERS. It does not design the
pH switch. That is deliberate and follows what worked on TNF-alpha: generating a good binder
and then scanning histidine placements onto it gave 35 measured placements across 14
lineages, where trying to optimise binding and pH response at once gave neither. The switch
is added in a later pass and scored with tfr1/pka_switch.py, which was validated against
known pH-dependent complexes first.

CYSTEINE IS EXCLUDED. Adaptyv expresses these cell-free, where free cysteines mispair and
the chain misfolds. On the previous challenge 13 of 20 selected designs carried a
buildability liability that was only noticed on the final day. Here it is prevented at
generation rather than filtered afterwards.

WHAT IS MEASURED, per design: contacts against the chosen epitope and against everything
else. The second number matters because the target is a crop, and a crop exposes interior
surface that does not exist in the real protein - about 3400 A^2 of it here, 90% of it more
than 15 A from the epitope. A design that binds that face will show many off-site contacts
and few on-site ones, and is rejected on that, not on confidence.
"""
import os
import subprocess
import sys
import textwrap
import time

W = "/kaggle/working"
TIME_BUDGET_MIN = 330           # 5.5 h. Sized against the GPU quota ACTUALLY left
                                # this week (12 h, not 30), because two further GPU
                                # passes still have to run on the same budget:
                                # histidine scanning for the pH switch, and mouse
                                # cross-reactivity prediction. Generation does not
                                # get to eat the whole allowance.
LEN_MIN, LEN_MAX = 55, 80

if not os.path.exists("%s/tfr1_targets.json" % W):
    raise SystemExit("run cellTFR1a_setup.py first - no targets found")

WORKER = r'''
import json, os, random, sys, time

TAG = sys.argv[1]; GPU = sys.argv[2]; BUDGET = float(sys.argv[3])
LEN_MIN = int(sys.argv[4]); LEN_MAX = int(sys.argv[5])
os.environ["CUDA_VISIBLE_DEVICES"] = GPU            # must precede the jax import
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ["XLA_PYTHON_CLIENT_ALLOCATOR"] = "platform"
t0 = time.time()
def log(m): print("[%s] %s [%.1f min]" % (TAG, m, (time.time()-t0)/60), flush=True)

import numpy as np
from colabdesign import mk_afdesign_model, clear_mem

W = "/kaggle/working"
MAPS = json.load(open("%s/tfr1_targets.json" % W))
M = MAPS[TAG]
TARGET, CHAIN = M["target"], M["chain"]
EPITOPE = set(M["epitope"])
RM_AA_OK = json.load(open("%s/tfr1_gate.json" % W)).get("rm_aa_ok", False)

KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5,
          L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)
def hydropathy(s): return sum(KD[c] for c in s) / len(s)

def longest_run(s):
    best = cur = 1
    for i in range(1, len(s)):
        cur = cur + 1 if s[i] == s[i-1] else 1
        best = max(best, cur)
    return best

def dedupe(chain):
    out, seen = [], set()
    for r in chain:
        if r.id[0] == " " and r.id[1] not in seen:
            seen.add(r.id[1]); out.append(r)
    return out

def measure(pdb_path, blen):
    from Bio.PDB import PDBParser
    mdl = PDBParser(QUIET=True).get_structure("c", pdb_path)[0]
    kept = {ch.id: dedupe(ch) for ch in mdl}
    # identify the target by its residue NUMBERS, never by "the bigger chain"
    tid = None
    for cid, rs in kept.items():
        if set(r.id[1] for r in rs) >= EPITOPE:
            tid = cid; break
    if tid is None:
        raise ValueError("no chain carries the epitope numbering; %s"
                         % {k: len(v) for k, v in kept.items()})
    bid = [c for c in kept if c != tid]
    if len(bid) != 1:
        raise ValueError("expected one binder chain, got %s" % bid)
    bres = kept[bid[0]]
    if len(bres) != blen:
        raise ValueError("binder chain %d residues, expected %d" % (len(bres), blen))

    ba = np.array([a.get_coord() for r in bres for a in r])
    on = off = 0
    touched = []
    for r in kept[tid]:
        ra = np.array([a.get_coord() for a in r])
        if np.min(np.linalg.norm(ra[:, None, :] - ba[None, :, :], axis=-1)) <= 4.5:
            if r.id[1] in EPITOPE:
                on += 1; touched.append(r.id[1])
            else:
                off += 1
    rec = dict(on_site=on, off_site=off, n_epitope=len(EPITOPE),
               on_fraction=round(on / float(max(on + off, 1)), 3),
               touched=sorted(touched))
    cas = [r["CA"].get_coord() for r in bres if "CA" in r]
    if len(cas) > 4:
        d = [float(np.linalg.norm(cas[i] - cas[i+4])) for i in range(len(cas)-4)]
        rec["helicity"] = round(sum(1 for x in d if 5.5 <= x <= 7.0) / len(d), 3)
    return rec

outdir = "%s/tfr1_%s" % (W, TAG)
os.makedirs(outdir, exist_ok=True)
rec_path = "%s/designs.jsonl" % outdir
rng = random.Random(7700 if TAG == "e559" else 7900)

made = i = 0
while (time.time() - t0) / 60 < BUDGET:
    i += 1
    name = "tfr1_%s_%03d" % (TAG, i)
    L = rng.randint(LEN_MIN, LEN_MAX)
    try:
        clear_mem()
        model = mk_afdesign_model(protocol="binder", data_dir=W)
        kw = dict(pdb_filename=TARGET, chain=CHAIN, binder_len=L,
                  hotspot=M["hotspot"])
        if RM_AA_OK:
            kw["rm_aa"] = "C"
        model.prep_inputs(**kw)
        model.design_3stage(soft_iters=60, temp_iters=30, hard_iters=10)
        seq = model.get_seqs()[0]
        dip = float(model.aux["log"].get("i_ptm", 0.0))
        try:
            model.predict(seq=seq, num_recycles=3, verbose=False)
        except TypeError:
            model.predict(seq=seq, verbose=False)
        lg = model.aux["log"]
        r = dict(name=name, site=TAG, length=L, sequence=seq,
                 design_iptm=round(dip, 3),
                 true_iptm=round(float(lg.get("i_ptm", lg.get("ptm", 0.0))), 3),
                 plddt=round(float(lg.get("plddt", 0.0)), 3),
                 hydropathy=round(hydropathy(seq), 3),
                 n_cys=seq.count("C"), longest_run=longest_run(seq))
        out = "%s/%s.pdb" % (outdir, name)
        model.save_pdb(out)
        try:
            r.update(measure(out, L))
        except Exception as e:
            r["measure_error"] = "%s %s" % (type(e).__name__, str(e)[:160])
        with open(rec_path, "a") as fh:
            fh.write(json.dumps(r) + "\n")
        made += 1
        log("%s len=%d i_ptm=%.3f pLDDT=%.3f kd=%+.2f hlx=%.2f | on-site %d/%d, off-site %d"
            " (%.0f%% on) | Cys %d run %d"
            % (name, L, r["true_iptm"], r["plddt"], r["hydropathy"], r.get("helicity", 0.0),
               r.get("on_site", 0), len(EPITOPE), r.get("off_site", 0),
               100 * r.get("on_fraction", 0), r["n_cys"], r["longest_run"]))
    except Exception as e:
        log("%s FAILED %s %s" % (name, type(e).__name__, str(e)[:200]))

log("done: %d designs from %d attempts" % (made, i))
if made == 0:
    raise SystemExit("worker %s produced nothing - read the failures above" % TAG)
'''

with open("%s/tfr1_worker.py" % W, "w") as fh:
    fh.write(textwrap.dedent(WORKER))

print("GPU state at worker launch (index, used, free):")
print(subprocess.run("nvidia-smi --query-gpu=index,memory.used,memory.free "
                     "--format=csv,noheader", shell=True,
                     capture_output=True, text=True).stdout.strip(), flush=True)

procs = []
for gpu, tag in enumerate(("e559", "l537")):
    os.makedirs("%s/tfr1_%s" % (W, tag), exist_ok=True)
    p = subprocess.Popen([sys.executable, "%s/tfr1_worker.py" % W, tag, str(gpu),
                          str(TIME_BUDGET_MIN), str(LEN_MIN), str(LEN_MAX)],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                         bufsize=1)
    procs.append(p)
    print("launched %s on GPU %d (pid %d)" % (tag, gpu, p.pid), flush=True)

import threading


def pump(p):
    for line in p.stdout:
        sys.stdout.write(line)
        sys.stdout.flush()


threads = [threading.Thread(target=pump, args=(p,), daemon=True) for p in procs]
for t in threads:
    t.start()
for p in procs:
    p.wait()
for t in threads:
    t.join(timeout=10)

import glob
import json

recs = []
for f in sorted(glob.glob("%s/tfr1_*/designs.jsonl" % W)):
    for line in open(f):
        line = line.strip()
        if line:
            recs.append(json.loads(line))

print("\n%d designs total" % len(recs))
if recs:
    by_site = {}
    for r in recs:
        by_site.setdefault(r["site"], []).append(r)
    for site, rs in sorted(by_site.items()):
        good = [r for r in rs if r.get("on_site", 0) >= 4]
        print("  %s: %d designs, %d with >=4 on-site contacts, best i_ptm %.3f"
              % (site, len(rs), len(good), max(r["true_iptm"] for r in rs)))
    with open("%s/tfr1_designs.json" % W, "w") as fh:
        json.dump(recs, fh, indent=0)
    print("\n===== COPY FROM HERE =====")
    for r in recs:
        print(json.dumps(r))
    print("===== COPY TO HERE =====")
