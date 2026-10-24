# === CHALLENGE 3 (GLP1R) : GENERATION, TWO JOBS AT ONCE ===
#
# Two different design problems run in parallel, one per T4, because they need different
# targets and there is no point rebuilding a model between them.
#
#   GPU 0, "ecd"    - bind the extracellular landing pad. 76-residue target, so cheap.
#                     This is the job our pipeline already does well, and the one the lab
#                     can measure reliably. It yields an ANTAGONIST by mechanism.
#   GPU 1, "thread" - thread the helical bundle, which is what actually activates the
#                     receptor. 192-residue target, long thin binder. A declared long shot.
#
# WHAT WE CAN AND CANNOT SCREEN FOR. Yesterday's control run (cellGLP2_control.py) tested
# whether AlphaFold can tell the real hormone from the truncated drug that blocks the
# receptor. It cannot: 34.1 A against 33.3 A, inside our noise, and a scrambled nonsense
# sequence reached 37.0 A - deeper than either. So depth MAGNITUDE is not a selection rule
# and nothing here pretends otherwise.
#
# What the same run did separate is depth POLARITY - which end of the peptide points down
# the channel. Both real peptides put their N-terminus at the deep end to within 0.1 A;
# the scrambled control had its N-terminus at -7.4 A while its deepest atom sat at 37.0 A,
# i.e. sprawled in backwards. That is a real "is this threaded at all" filter, and it is
# the only thing we screen threaders on besides confidence. Activation itself is UNSCREENED
# and the methods file has to say so.
#
# Also measured for every design, because the competition records it and we only found out
# late: mean Kyte-Doolittle hydropathy. Our TNF designs were greasy enough to have failed
# that screen. Measuring it during generation rather than after means we can select on it.
import os, subprocess, sys, textwrap, time

TIME_BUDGET_MIN = 420          # 7 h wall clock; Kaggle sessions cap at 12
SOFT, TEMP, HARD = 90, 40, 20

WORKER = r'''
import json, os, random, sys, time
GPU, MODE, BUDGET = sys.argv[1], sys.argv[2], float(sys.argv[3])
os.environ["CUDA_VISIBLE_DEVICES"] = GPU          # must precede the jax import
SOFT, TEMP, HARD = int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])

t0 = time.time()
def log(m): print("[%s] %s [%.1f min]" % (MODE, m, (time.time()-t0)/60), flush=True)

import numpy as np
from colabdesign import mk_afdesign_model, clear_mem

FULL = json.load(open("/kaggle/working/glp1r_map.json"))
AXIS = np.array(FULL["axis"]); ECD_C = np.array(FULL["ecd_centroid"])
AXIS_OFF = float(ECD_C @ AXIS)

if MODE == "ecd":
    M = json.load(open("/kaggle/working/glp1r_ecd_map.json"))
    TARGET = "/kaggle/working/targets/GLP1R_ecd_crop.pdb"
    HOTSPOT = open("/kaggle/working/glp1r_ecd_hotspot.txt").read().strip()
    LEN_MIN, LEN_MAX = 55, 80
    WANT = set(M["contacts"])
else:
    M = FULL
    TARGET = "/kaggle/working/targets/GLP1R_crop.pdb"
    HOTSPOT = open("/kaggle/working/glp1r_tm_hotspot.txt").read().strip()
    LEN_MIN, LEN_MAX = 30, 45
    WANT = set(FULL["tm_contacts"])

ORDER = M["order"]
TSEQ = "".join(d["aa"] for d in ORDER)
log("target %s, %d residues, %d hotspots, binder %d-%d"
    % (os.path.basename(TARGET), len(ORDER), len(HOTSPOT.split(",")), LEN_MIN, LEN_MAX))

# Kyte-Doolittle hydropathy. The competition records the mean of this per design, and a
# greasy sequence reads as an aggregation risk, so it is a selection variable not a note.
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5,
          L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)
def hydropathy(seq):
    vals = [KD[c] for c in seq if c in KD]
    return sum(vals) / len(vals) if vals else 0.0

def dedupe(chain):
    """AfDesign output carries each residue twice; keep the first copy of each number."""
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
    # Identify the target by its LENGTH, not by "whichever chain is smaller". The
    # extracellular target is 76 residues and a binder may be 80, so smaller-is-binder
    # would silently swap them and every number below would describe the wrong molecule.
    tid = None
    for cid, rs in kept.items():
        if len(rs) == len(ORDER):
            s = "".join(T31.get(r.get_resname(), "X") for r in rs)
            if s == TSEQ:
                tid = cid; break
    if tid is None:
        raise ValueError("no chain matches the target sequence; chains %s"
                         % {k: len(v) for k, v in kept.items()})
    bid = [c for c in kept if c != tid]
    if len(bid) != 1:
        raise ValueError("expected exactly one binder chain, got %s" % bid)
    bres = kept[bid[0]]
    if len(bres) != blen:
        raise ValueError("binder chain %d residues, expected %d" % (len(bres), blen))

    batoms = np.array([a.get_coord() for r in bres for a in r])
    hits = []
    for k, r in enumerate(kept[tid]):
        ra = np.array([a.get_coord() for a in r])
        if np.min(np.linalg.norm(ra[:, None, :] - batoms[None, :, :], axis=-1)) <= 5.0:
            hits.append(ORDER[k]["uniprot"])
    rec = dict(epitope=sorted(hits), on_site=sorted(set(hits) & WANT))

    # Helicity without DSSP: in an alpha helix CA(i) to CA(i+4) is about 6.2 A.
    cas = [r["CA"].get_coord() for r in bres if "CA" in r]
    if len(cas) > 4:
        d = [float(np.linalg.norm(cas[i] - cas[i + 4])) for i in range(len(cas) - 4)]
        rec["helicity"] = round(sum(1 for x in d if 5.5 <= x <= 7.0) / len(d), 3)

    if MODE == "thread":
        proj = batoms @ AXIS - AXIS_OFF
        nterm = np.array([a.get_coord() for a in bres[0]]) @ AXIS - AXIS_OFF
        cterm = np.array([a.get_coord() for a in bres[-1]]) @ AXIS - AXIS_OFF
        rec["deepest"] = round(float(np.max(proj)), 1)
        rec["nterm_depth"] = round(float(np.max(nterm)), 1)
        rec["cterm_depth"] = round(float(np.max(cterm)), 1)
        # the polarity screen: is the N-terminal end actually the deep end
        rec["threaded"] = bool(rec["nterm_depth"] >= rec["deepest"] - 5.0
                               and rec["nterm_depth"] > rec["cterm_depth"])
        tm = ecd = 0
        for k, r in enumerate(kept[tid]):
            ra = np.array([a.get_coord() for a in r])
            if np.min(np.linalg.norm(ra[:, None, :] - batoms[None, :, :], axis=-1)) <= 4.5:
                if ORDER[k].get("region") == "tm":
                    tm += 1
                else:
                    ecd += 1
        rec["tm_contacts"], rec["ecd_contacts"] = tm, ecd
    return rec

outdir = "/kaggle/working/glp1r_%s" % MODE
os.makedirs(outdir, exist_ok=True)
rec_path = "%s/designs.jsonl" % outdir
rng = random.Random(9000 if MODE == "ecd" else 9500)

made = 0
i = 0
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
        design_iptm = float(model.aux["log"].get("i_ptm", 0.0))

        # The score that counts is an INDEPENDENT re-prediction of the discrete sequence,
        # never the optimiser's own number, which is optimistic by construction.
        try:
            model.predict(seq=seq, num_recycles=3, verbose=False)
        except TypeError:
            model.predict(seq=seq, verbose=False)
        lg = model.aux["log"]
        r = dict(name=name, mode=MODE, length=L, sequence=seq,
                 design_iptm=round(design_iptm, 3),
                 true_iptm=round(float(lg.get("i_ptm", lg.get("ptm", 0.0))), 3),
                 plddt=round(float(lg.get("plddt", 0.0)), 3),
                 hydropathy=round(hydropathy(seq), 3))

        out = "%s/%s.pdb" % (outdir, name)
        model.save_pdb(out)
        # Geometry is measured in its own guarded step: a failure here must not discard a
        # design whose confidence numbers are already in hand.
        try:
            r.update(measure(out, L))
        except Exception as e:
            r["measure_error"] = "%s %s" % (type(e).__name__, str(e)[:160])

        with open(rec_path, "a") as fh:
            fh.write(json.dumps(r) + "\n")
        made += 1
        extra = ""
        if MODE == "thread" and "deepest" in r:
            extra = " deep=%.1f nterm=%.1f %s %dTM" % (
                r["deepest"], r["nterm_depth"],
                "THREADED" if r["threaded"] else "backwards", r["tm_contacts"])
        log("%s len=%d i_ptm=%.3f pLDDT=%.3f kd=%+.2f hlx=%.2f on_site=%d%s"
            % (name, L, r["true_iptm"], r["plddt"], r["hydropathy"],
               r.get("helicity", 0.0), len(r.get("on_site", [])), extra))
    except Exception as e:
        log("%s FAILED %s %s" % (name, type(e).__name__, str(e)[:200]))

log("done: %d designs saved from %d attempts" % (made, i))
if made == 0:
    raise SystemExit("worker %s produced nothing - read the failures above" % MODE)
'''

os.makedirs("/kaggle/working/glp1r_ecd", exist_ok=True)
os.makedirs("/kaggle/working/glp1r_thread", exist_ok=True)
with open("/kaggle/working/glp1r_worker.py", "w") as fh:
    fh.write(textwrap.dedent(WORKER))

procs = []
for gpu, mode in (("0", "ecd"), ("1", "thread")):
    p = subprocess.Popen([sys.executable, "-u", "/kaggle/working/glp1r_worker.py",
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
threads = [threading.Thread(target=pump, args=(p,), daemon=True) for p in procs]
[t.start() for t in threads]
for p in procs:
    p.wait()
[t.join(timeout=10) for t in threads]

import glob, json
recs = []
for f in sorted(glob.glob("/kaggle/working/glp1r_*/designs.jsonl")):
    for line in open(f):
        line = line.strip()
        if line:
            recs.append(json.loads(line))
json.dump(recs, open("/kaggle/working/glp1r_designs.json", "w"), indent=1)
if not recs:
    raise SystemExit("no designs at all from either job - this run produced nothing")

for mode in ("ecd", "thread"):
    sub = sorted([r for r in recs if r["mode"] == mode], key=lambda r: -r["true_iptm"])
    print("\n==== %s: %d designs ====" % (mode.upper(), len(sub)))
    hdr = "%-18s %4s %7s %7s %7s %6s %7s" % ("name", "len", "i_ptm", "pLDDT",
                                             "kd", "hlx", "on_site")
    print(hdr + ("  depth  nterm  threaded" if mode == "thread" else ""))
    for r in sub[:25]:
        line = "%-18s %4d %7.3f %7.3f %+7.2f %6.2f %7d" % (
            r["name"], r["length"], r["true_iptm"], r["plddt"], r["hydropathy"],
            r.get("helicity", 0.0), len(r.get("on_site", [])))
        if mode == "thread":
            line += "  %5.1f  %5.1f  %s" % (r.get("deepest", 0.0),
                                            r.get("nterm_depth", 0.0),
                                            r.get("threaded", False))
        print(line)
    if sub:
        good = [r for r in sub if r["true_iptm"] >= 0.50]
        onsite = [r for r in sub if len(r.get("on_site", [])) >= 5]
        cool = [r for r in sub if r["hydropathy"] <= 0.0]
        print("%d at i_ptm>=0.50, %d contact 5+ intended residues, %d not greasy (kd<=0)"
              % (len(good), len(onsite), len(cool)))
        if mode == "thread":
            thr = [r for r in sub if r.get("threaded")]
            print("%d of %d are threaded the right way round (N-terminus deepest)"
                  % (len(thr), len(sub)))

print("\n===== COPY FROM HERE =====")
for r in recs:
    print(json.dumps(r, separators=(",", ":")))
print("===== COPY TO HERE =====")
