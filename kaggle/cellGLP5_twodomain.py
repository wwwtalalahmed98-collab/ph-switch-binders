# === CHALLENGE 3 (GLP1R) : TWO-DOMAIN AGONISTS, AND THREADERS RE-RUN ON A SEALED TARGET ===
#
# GPU 0, "both"   - the design we never tried. The natural hormone gets its GRIP from the
#                   extracellular groove and its ACTIVATION from threading the channel:
#                   one molecule doing two jobs. We built either/or - short threaders aimed
#                   only at the channel, folded binders aimed only at the groove - and that
#                   separation is why our binders cannot activate and our threaders cannot
#                   grip. Nothing required it; it was my choice and it was wrong. These are
#                   45-70 residues aimed at BOTH contact sets at once.
#
# GPU 1, "thread2" - the previous threading job, unchanged, on the sealed target. Same
#                   design task, one variable altered. If the failure was the target, this
#                   succeeds; if it was the method, it fails the same way. Either result is
#                   worth having, which is what makes it worth a GPU.
#
# What is screened, and what is not. Designs are filtered on reaching the hormone's own
# contact residues and on not diving past where the target stops representing the real
# receptor. ACTIVATION IS NOT SCREENED. Nothing available to us tests whether a design
# switches the receptor on; we can only build the shape that makes it possible. The methods
# file says so in those words.
import os, subprocess, sys, textwrap, time

TIME_BUDGET_MIN = 480          # 8 h. The quota reads 23 h remaining, not the 16 I
                               # assumed: Kaggle bills a dual-T4 session at one quota-hour
                               # per wall-clock hour, not two. Wall clock, not quota, is the
                               # constraint now, and 8 h still lands well before the deadline.
SOFT, TEMP, HARD = 90, 40, 20

WORKER = r'''
import json, os, random, sys, time
GPU, MODE, BUDGET = sys.argv[1], sys.argv[2], float(sys.argv[3])
os.environ["CUDA_VISIBLE_DEVICES"] = GPU          # must precede the jax import
# The sealed target is 283 residues against the old crop's 192, and AF2's memory grows
# with the square of the total chain length. The first attempt died with
# RESOURCE_EXHAUSTED on every design. Defaulting to the platform allocator stops jax
# grabbing one large block up front and then fragmenting; it is slower and it fits.
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ["XLA_PYTHON_CLIENT_ALLOCATOR"] = "platform"
SOFT, TEMP, HARD = int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])

t0 = time.time()
def log(m): print("[%s] %s [%.1f min]" % (MODE, m, (time.time()-t0)/60), flush=True)

import numpy as np
from colabdesign import mk_afdesign_model, clear_mem

M = json.load(open("/kaggle/working/glp1r_sealed_map.json"))
ORDER = M["order"]
AXIS = np.array(M["axis"]); OFF = float(np.array(M["ecd_centroid"]) @ AXIS)
WANT_TM, WANT_ECD = set(M["tm_contacts"]), set(M["ecd_contacts"])
TARGET = "/kaggle/working/targets/GLP1R_sealed.pdb"
TSEQ = "".join(d["aa"] for d in ORDER)

if MODE == "both":
    HOTSPOT = open("/kaggle/working/hot_both.txt").read().strip()
    # Was 45-70, and every draw above ~55 ran the T4 out of memory on the sealed
    # target. 45-70 was over-generous anyway: GLP-1 itself grips the groove AND
    # reaches 33 A down the channel in 31 residues, so 40-50 is already longer than
    # the molecule we are imitating. Nothing is given up by capping it here.
    LEN_MIN, LEN_MAX = 40, 50
else:
    HOTSPOT = open("/kaggle/working/hot_tm.txt").read().strip()
    LEN_MIN, LEN_MAX = 30, 45
log("target %d residues, %d hotspots, binder %d-%d"
    % (len(ORDER), len(HOTSPOT.split(",")), LEN_MIN, LEN_MAX))

KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5,
          L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)
def hydropathy(s):
    v = [KD[c] for c in s if c in KD]
    return sum(v) / len(v) if v else 0.0

def dedupe(chain):
    res, seen = [], set()
    for r in list(chain):
        if r.id[0] != " " or r.id[1] in seen:
            chain.detach_child(r.id); continue
        seen.add(r.id[1]); res.append(r)
    return res

def measure(pdb_path, blen):
    from Bio.PDB import PDBParser
    from Bio.PDB.Polypeptide import protein_letters_3to1 as T31
    mdl = PDBParser(QUIET=True).get_structure("c", pdb_path)[0]
    kept = {ch.id: dedupe(ch) for ch in mdl}
    # identify the receptor by sequence, never by "the bigger chain"
    rid = None
    for cid, rs in kept.items():
        if len(rs) == len(ORDER):
            if "".join(T31.get(r.get_resname(), "X") for r in rs) == TSEQ:
                rid = cid; break
    if rid is None:
        raise ValueError("no chain matches the target sequence; %s"
                         % {k: len(v) for k, v in kept.items()})
    bid = [c for c in kept if c != rid]
    if len(bid) != 1:
        raise ValueError("expected one binder chain, got %s" % bid)
    bres = kept[bid[0]]
    if len(bres) != blen:
        raise ValueError("binder chain %d residues, expected %d" % (len(bres), blen))

    ba = np.array([a.get_coord() for r in bres for a in r])
    tm_on = ecd_on = tm_all = ecd_all = 0
    for k, r in enumerate(kept[rid]):
        ra = np.array([a.get_coord() for a in r])
        if np.min(np.linalg.norm(ra[:, None, :] - ba[None, :, :], axis=-1)) <= 4.5:
            u = ORDER[k]["uniprot"]
            if ORDER[k]["region"] == "tm":
                tm_all += 1
                if u in WANT_TM: tm_on += 1
            else:
                ecd_all += 1
                if u in WANT_ECD: ecd_on += 1

    proj = ba @ AXIS - OFF
    nt = np.array([a.get_coord() for a in bres[0]]) @ AXIS - OFF
    ct = np.array([a.get_coord() for a in bres[-1]]) @ AXIS - OFF
    rec = dict(tm_contacts=tm_all, ecd_contacts=ecd_all, tm_on=tm_on, ecd_on=ecd_on,
               deepest=round(float(np.max(proj)), 1),
               nterm_depth=round(float(np.max(nt)), 1),
               cterm_depth=round(float(np.max(ct)), 1))
    rec["threaded"] = bool(rec["nterm_depth"] >= rec["deepest"] - 5.0
                           and rec["nterm_depth"] > rec["cterm_depth"])
    cas = [r["CA"].get_coord() for r in bres if "CA" in r]
    if len(cas) > 4:
        d = [float(np.linalg.norm(cas[i] - cas[i + 4])) for i in range(len(cas) - 4)]
        rec["helicity"] = round(sum(1 for x in d if 5.5 <= x <= 7.0) / len(d), 3)
    return rec

outdir = "/kaggle/working/glp1r_%s" % MODE
os.makedirs(outdir, exist_ok=True)
rec_path = "%s/designs.jsonl" % outdir
rng = random.Random(11000 if MODE == "both" else 11500)

made = i = 0
while (time.time() - t0) / 60 < BUDGET:
    i += 1
    name = "glp1r_%s_%03d" % (MODE, i)
    L = rng.randint(LEN_MIN, LEN_MAX)
    try:
        clear_mem()
        model = mk_afdesign_model(protocol="binder", data_dir="/kaggle/working")
        model.prep_inputs(pdb_filename=TARGET, chain="R", binder_len=L, hotspot=HOTSPOT)
        model.design_3stage(soft_iters=SOFT, temp_iters=TEMP, hard_iters=HARD)
        seq = model.get_seqs()[0]
        dip = float(model.aux["log"].get("i_ptm", 0.0))
        try:
            model.predict(seq=seq, num_recycles=3, verbose=False)
        except TypeError:
            model.predict(seq=seq, verbose=False)
        lg = model.aux["log"]
        r = dict(name=name, mode=MODE, length=L, sequence=seq,
                 design_iptm=round(dip, 3),
                 true_iptm=round(float(lg.get("i_ptm", lg.get("ptm", 0.0))), 3),
                 plddt=round(float(lg.get("plddt", 0.0)), 3),
                 hydropathy=round(hydropathy(seq), 3))
        out = "%s/%s.pdb" % (outdir, name)
        model.save_pdb(out)
        try:
            r.update(measure(out, L))
        except Exception as e:
            r["measure_error"] = "%s %s" % (type(e).__name__, str(e)[:160])
        with open(rec_path, "a") as fh:
            fh.write(json.dumps(r) + "\n")
        made += 1
        log("%s len=%d i_ptm=%.3f pLDDT=%.3f kd=%+.2f hlx=%.2f | grip %d ECD / reach %d TM"
            " | deep=%.1f nterm=%.1f %s"
            % (name, L, r["true_iptm"], r["plddt"], r["hydropathy"],
               r.get("helicity", 0.0), r.get("ecd_on", 0), r.get("tm_on", 0),
               r.get("deepest", 0.0), r.get("nterm_depth", 0.0),
               "THREADED" if r.get("threaded") else "backwards"))
    except Exception as e:
        log("%s FAILED %s %s" % (name, type(e).__name__, str(e)[:200]))

log("done: %d designs from %d attempts" % (made, i))
if made == 0:
    raise SystemExit("worker %s produced nothing - read the failures above" % MODE)
'''

for d in ("both", "thread2"):
    os.makedirs("/kaggle/working/glp1r_%s" % d, exist_ok=True)
with open("/kaggle/working/glp1r_worker2.py", "w") as fh:
    fh.write(textwrap.dedent(WORKER))

procs = []
for gpu, mode in (("0", "both"), ("1", "thread2")):
    p = subprocess.Popen([sys.executable, "-u", "/kaggle/working/glp1r_worker2.py",
                          gpu, mode, str(TIME_BUDGET_MIN),
                          str(SOFT), str(TEMP), str(HARD)],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         text=True, bufsize=1)
    procs.append(p)
    print("launched %s job on GPU %s" % (mode, gpu), flush=True)

import threading
def pump(p):
    for line in p.stdout:
        print(line.rstrip(), flush=True)
th = [threading.Thread(target=pump, args=(p,), daemon=True) for p in procs]
[x.start() for x in th]
for p in procs:
    p.wait()
[x.join(timeout=10) for x in th]

import glob, json
recs = []
for f in sorted(glob.glob("/kaggle/working/glp1r_*/designs.jsonl")):
    for line in open(f):
        line = line.strip()
        if line:
            recs.append(json.loads(line))
json.dump(recs, open("/kaggle/working/glp1r_designs2.json", "w"), indent=1)
if not recs:
    raise SystemExit("no designs at all from either job")

for mode in ("both", "thread2"):
    sub = sorted([r for r in recs if r["mode"] == mode],
                 key=lambda r: -(r.get("tm_on", 0)))
    print("\n==== %s: %d designs ====" % (mode.upper(), len(sub)))
    print("%-20s %4s %7s %6s %6s %6s %7s %7s" % ("name", "len", "i_ptm", "kd",
                                                 "ECDon", "TMon", "deepest", "nterm"))
    for r in sub[:25]:
        print("%-20s %4d %7.3f %+6.2f %6d %6d %7.1f %7.1f %s"
              % (r["name"], r["length"], r["true_iptm"], r["hydropathy"],
                 r.get("ecd_on", 0), r.get("tm_on", 0), r.get("deepest", 0),
                 r.get("nterm_depth", 0), "THREADED" if r.get("threaded") else "back"))
    if sub:
        both_ok = [r for r in sub if r.get("ecd_on", 0) >= 4 and r.get("tm_on", 0) >= 5
                   and r.get("threaded") and r.get("deepest", 99) <= 45]
        print("%d of %d grip AND reach (>=4 ECD, >=5 TM, threaded, within 45 A)"
              % (len(both_ok), len(sub)))

print("\n===== COPY FROM HERE =====")
for r in recs:
    print(json.dumps(r, separators=(",", ":")))
print("===== COPY TO HERE =====")
