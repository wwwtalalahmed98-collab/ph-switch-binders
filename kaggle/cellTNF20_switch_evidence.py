# === CHALLENGE 2 (TNF-alpha) : EVIDENCE FOR THE SWITCH, AND ERROR BARS ===
#
# Everything submitted so far rests on binding at pH 7.4 against two species. Objective 1
# - the switch - has never been computed. No calculation in this project touches
# protonation. This run produces the first evidence about it, plus the error bars the
# ranking has been missing, plus an answer to a question the Proteinbase portal raised.
#
# Four measurements, on the 20 designs actually submitted (re-derived here from the
# switch5 output by the same rule, so no sequence is ever written into the repo):
#
# 1. CHARGE SENSITIVITY. The switch claim is that a protonated histidine cannot be
#    buried against a target cation. AlphaFold2 knows nothing about pH, but it does
#    respond to charge. So at each histidine position, predict Gln (neutral, similar
#    size) and Lys (charged, similar length). The DIFFERENCE isolates charge from
#    sterics. If Lys collapses the interface and Gln does not, the position genuinely
#    punishes positive charge and the switch should fire. If Lys is tolerated, our
#    switch argument for that design is unsupported. AF2 is a poor dG predictor, so
#    only large effects are interpretable - but this is falsifiable, which is new.
#
# 2. MODEL-ENSEMBLE ERROR BARS. Every i_ptm reported so far is a single prediction from
#    one parameter set. Across-model spread is typically +/-0.05 to 0.10, the same size
#    as several margins the ranking depends on. Five models, mean and sd.
#
# 3. MONOMER FOLDABILITY. The portal reported pLDDT 55-63 where our complex predictions
#    gave 60-77. The likely explanation is that it predicted the binder ALONE. If these
#    are target-dependent folders they are poor bets to express. Predicted here directly
#    rather than inferred.
#
# 4. BURIAL ON BINDING. The proposal geometry used Cbeta-to-cation distance with a
#    3.0-8.8 A window. That is the wrong variable. The switch is a binding-coupled
#    protonation: in the free binder the histidine is exposed and protonates at pH 6;
#    binding has to bury that charge, and the burial is what costs affinity. So measure
#    histidine side-chain rel_SASA free vs bound, the real Nd1/Ne2-to-cation distance
#    after folding, and a Debye-screened Coulomb estimate at 150 mM. CPU only.
#
# Nothing here redesigns anything. No sequence is changed except the deliberate
# Gln/Lys probes, which are never submitted.
import glob, json, os, subprocess, sys, textwrap, time

TIME_BUDGET_MIN = 210
N_SWITCH = 12
TOTAL_SLOTS = 20
N_MODELS = 5

WORKER = r"""
import glob, json, os, sys, time
from collections import Counter
GPU, SHARD, BUDGET = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
N_SWITCH, TOTAL_SLOTS, N_MODELS = int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
os.environ["CUDA_VISIBLE_DEVICES"] = GPU
t0 = time.time()
def log(m): print("[w%d] %s [%.1f min]" % (SHARD, m, (time.time()-t0)/60), flush=True)

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import protein_letters_3to1

HT = "/kaggle/working/targets/TNF_site_crop.pdb"
MT = "/kaggle/working/targets/TNF_mouse_site_crop.pdb"
CMAP = json.load(open("/kaggle/working/construct_map.json"))["order"]
H_HOTSPOT = "A86,A87,A88,A90,A91,A92,B29,B32,B33,B34,B35,B36"
M_HOTSPOT = open("/kaggle/working/mouse_hotspot.txt").read().strip()
CAT_N = {"ARG": ("NE", "NH1", "NH2"), "LYS": ("NZ",)}
MAXS = {"A":129,"R":274,"N":195,"D":193,"C":167,"E":223,"Q":225,"G":104,"H":224,
        "I":197,"L":201,"K":236,"M":224,"F":240,"P":159,"S":155,"T":172,"W":285,
        "Y":263,"V":174}
# Debye-Huckel at 150 mM monovalent, 298 K: kappa^-1 = 7.9 A. eps_eff = 30 is the
# usual compromise for a partially desolvated protein-protein interface; the absolute
# number is soft, the RANKING across positions is what this is for.
EPS_EFF, DEBYE = 30.0, 7.9

# ---- rebuild the submitted set from the switch5 output, by the same rule -----------
src = glob.glob("/kaggle/input/**/tnf_switch5.json", recursive=True)
if not src:
    raise SystemExit("tnf_switch5.json not found - attach the tnf_switch5_repair_fix "
                     "output. Do NOT read Kaggle's file PREVIEW for this: it truncates "
                     "at 102400 bytes and closes the JSON brackets, so a short read "
                     "parses as valid and silently drops records.")
recs = json.load(open(src[0]))
if len(recs) < 26:
    raise SystemExit("only %d records in tnf_switch5.json - expected 29; the file looks "
                     "truncated, refusing to select from it" % len(recs))

def rank_key(v):
    return (v["n_his"], v["n_anchor"], min(v["human_i_ptm"], v["mouse_i_ptm"]))

chosen = []
for r in recs:
    ok = [v for v in r.get("variants", []) if v.get("passes")]
    if not ok:
        continue
    best = max(ok, key=rank_key)
    best = dict(best); best["parent"] = r["name"]; best["control"] = r["control"]
    chosen.append(best)
chosen.sort(key=rank_key, reverse=True)
picked = chosen[:N_SWITCH]
n_ctl = TOTAL_SLOTS - len(picked)
clean = [v for v in picked if "H" not in v["control"]["sequence"]]
with_ctl = {v["parent"] for v in clean[:n_ctl]}
if SHARD == 0:
    log("%d lineages usable; %d switches, %d controls" % (len(chosen), len(picked), n_ctl))

jobs = []
for i, v in enumerate(picked, 1):
    jobs.append(dict(slot="sw_%02d" % i, role="switch", parent=v["parent"],
                     seq=v["sequence"], his=v["his_substitutions"],
                     ctl_seq=v["control"]["sequence"], ref_h=v["human_i_ptm"],
                     ref_m=v["mouse_i_ptm"]))
    if v["parent"] in with_ctl:
        jobs.append(dict(slot="ct_%02d" % i, role="control", parent=v["parent"],
                         seq=v["control"]["sequence"], his="",
                         ctl_seq=v["control"]["sequence"],
                         ref_h=v["control"]["human_i_ptm"],
                         ref_m=v["control"]["mouse_i_ptm"]))
# Sorting by length so each model rebuild is amortised over every job of that length.
jobs.sort(key=lambda j: (len(j["seq"]), j["slot"]))
mine = [j for k, j in enumerate(jobs) if k % 2 == SHARD]
log("this worker takes %d of %d designs: %s"
    % (len(mine), len(jobs), ", ".join(j["slot"] for j in mine)))

from colabdesign import mk_afdesign_model, clear_mem
_c = {}

def _build(species, L):
    # ColabDesign loads only num_models parameter sets, default 1. Asking predict() for
    # models=[1] against a single loaded set is what killed the whole of run 1 with an
    # IndexError, so the count is requested up front AND verified below, never assumed.
    def mk(**extra):
        if species == "alone":
            return mk_afdesign_model(protocol="hallucination", length=L,
                                     data_dir="/kaggle/working", **extra)
        af = mk_afdesign_model(protocol="binder", data_dir="/kaggle/working", **extra)
        af.prep_inputs(pdb_filename=HT if species == "human" else MT, chain="A,B",
                       binder_len=L,
                       hotspot=H_HOTSPOT if species == "human" else M_HOTSPOT)
        return af
    try:
        return mk(num_models=N_MODELS)
    except TypeError:
        return mk()

def _model(species, L):
    key = (species, L)
    if _c.get("key") != key:
        clear_mem()
        af = _build(species, L)
        n = len(getattr(af, "_model_params", None) or [])
        _c["key"], _c["af"], _c["n_models"] = key, af, n
        log("built %s model for length %d with %d parameter set(s) loaded"
            % (species, L, n))
    return _c["af"]

def n_models_available():
    return _c.get("n_models", 0)

def predict(seq, species, out_pdb=None, model=None):
    af = _model(species, len(seq))
    if model is not None and model >= max(1, n_models_available()):
        return None                # this parameter set is not loaded; not an error
    kw = dict(seq=seq, num_recycles=3, verbose=False)
    if model is not None:
        kw["models"] = [model]
    try:
        af.predict(**kw)
    except TypeError:
        kw.pop("models", None)
        kw.pop("num_recycles", None)
        af.predict(**kw)
        if model is not None:
            return None            # per-model selection unavailable in this build
    except IndexError:
        return None                # model index rejected; caller keeps what it has
    lg = af.aux["log"]
    if out_pdb:
        af.save_pdb(out_pdb)
    return (round(float(lg.get("i_ptm", lg.get("ptm", 0.0))), 3),
            round(float(lg.get("plddt", 0.0)), 3))

# ---- burial on binding, from the complex PDB --------------------------------------
# A chain written twice overlaps itself exactly, so every one of its residues reports
# as fully buried and the whole measurement silently inverts. This project has already
# been bitten by duplicated chains in a saved PDB, so every chain is de-duplicated
# BEFORE any SASA is computed, not after.
def dedupe(chain):
    res, seen = [], set()
    for r in list(chain):
        if r.id[0] != " ":
            chain.detach_child(r.id); continue
        if r.id[1] in seen:
            chain.detach_child(r.id); continue
        seen.add(r.id[1]); res.append(r)
    return res

def burial(pdb_path, blen, his_idx):
    # rel_SASA of each histidine side chain in the complex vs in the binder alone,
    # the real Nd1/Ne2-to-cation distance after folding, and a screened Coulomb term.
    mdl = PDBParser(QUIET=True).get_structure("c", pdb_path)[0]
    kept = {ch.id: dedupe(ch) for ch in mdl}
    bnd_id = min(kept, key=lambda k: len(kept[k]))
    bres = kept[bnd_id]
    if len(bres) != blen:
        raise ValueError("binder chain %s has %d unique residues, expected %d"
                         % (bnd_id, len(bres), blen))
    # cationic nitrogens on the target, from the complex as folded
    cats = []
    for cid, rs in kept.items():
        if cid == bnd_id:
            continue
        for r in rs:
            for an in CAT_N.get(r.get_resname(), ()):
                if an in r:
                    cats.append((cid, r.id[1], an, r[an].get_coord()))
    if not cats:
        raise ValueError("no Arg/Lys nitrogens found on the target chains")
    # SASA twice: whole de-duplicated complex, then the binder chain on its own
    ShrakeRupley().compute(mdl, level="R")
    bound = {}
    for k, r in enumerate(bres):
        aa = protein_letters_3to1.get(r.get_resname())
        if aa:
            bound[k] = r.sasa / MAXS.get(aa, 200)
    import copy
    solo = copy.deepcopy(mdl[bnd_id])
    sres = [r for r in solo if r.id[0] == " "]
    ShrakeRupley().compute(solo, level="R")
    free = {}
    for k, r in enumerate(sres):
        aa = protein_letters_3to1.get(r.get_resname())
        if aa:
            free[k] = r.sasa / MAXS.get(aa, 200)
    out = {}
    for k in his_idx:
        r = bres[k]
        best = None
        for an in ("ND1", "NE2"):
            if an not in r:
                continue
            for cid, cnum, catn, xyz in cats:
                d = float(np.linalg.norm(r[an].get_coord() - xyz))
                if best is None or d < best[0]:
                    best = (d, "%s%d%s" % (cid, cnum, catn), an)
        ddg = None
        if best is not None:
            # 332 * q1*q2 / (eps * r) * exp(-r/lambda), kcal/mol, both charges +1
            ddg = round(332.0 / (EPS_EFF * best[0]) * np.exp(-best[0] / DEBYE), 3)
        out[k] = dict(rel_sasa_free=round(free.get(k, -1), 3),
                      rel_sasa_bound=round(bound.get(k, -1), 3),
                      buried_on_binding=round(free.get(k, 0) - bound.get(k, 0), 3),
                      nearest_cation=best[1] if best else None,
                      cation_dist=round(best[0], 2) if best else None,
                      which_n=best[2] if best else None,
                      screened_coulomb_kcal=ddg)
    return out

os.makedirs("/kaggle/working/ev_pdb", exist_ok=True)
rec_path = "/kaggle/working/tnf_evidence_w%d.jsonl" % SHARD

def emit(d):
    with open(rec_path, "a") as fh:
        fh.write(json.dumps(d) + chr(10))

for j in mine:
    if (time.time() - t0) / 60 > BUDGET:
        log("%s not started: time budget reached" % j["slot"])
        emit(dict(slot=j["slot"], role=j["role"], parent=j["parent"],
                  incomplete="time budget reached before this design"))
        continue
    seq, L = j["seq"], len(j["seq"])
    rec = dict(slot=j["slot"], role=j["role"], parent=j["parent"], length=L,
               his=j["his"], ref_human_i_ptm=j["ref_h"], ref_mouse_i_ptm=j["ref_m"])
    pdb0 = "/kaggle/working/ev_pdb/%s_human.pdb" % j["slot"]

    # Run 1 put the OPTIONAL measurement first and let it take down the irreplaceable
    # one: a model-index IndexError in the error-bar loop killed every design before the
    # charge test ran, and Kaggle still reported "Successfully ran". So the order is now
    # most-important-first, and every step is isolated - one failure costs that step
    # only, and whatever succeeded is still emitted.
    def step(name, fn):
        try:
            fn()
        except Exception as e:
            rec["%s_error" % name] = "%s %s" % (type(e).__name__, str(e)[:160])
            log("%s %s FAILED %s" % (j["slot"], name, rec["%s_error" % name]))

    # --- 0. baseline, single model, and the complex structure burial needs
    def _baseline():
        h = predict(seq, "human", out_pdb=pdb0)
        m = predict(seq, "mouse")
        rec["human_i_ptm"], rec["human_plddt"] = h
        rec["mouse_i_ptm"], rec["mouse_plddt"] = m
        log("%s baseline human %.3f mouse %.3f (switch5 said %.3f/%.3f)"
            % (j["slot"], h[0], m[0], j["ref_h"], j["ref_m"]))
    step("baseline", _baseline)

    hidx = []
    if j["role"] == "switch" and j["his"]:
        hidx = [int(s[1:-1]) - 1 for s in j["his"].split(",")]

    # --- 1. charge sensitivity: the only falsifiable test of objective 1 we have
    def _charge():
        probes = []
        for pos in hidx:
            if seq[pos] != "H":
                raise ValueError("expected H at %d, found %s" % (pos + 1, seq[pos]))
            for aa, label in (("Q", "neutral"), ("K", "charged")):
                mut = seq[:pos] + aa + seq[pos+1:]
                hp = predict(mut, "human")
                mp = predict(mut, "mouse")
                probes.append(dict(position=pos+1, residue=aa, kind=label,
                                   human_i_ptm=hp[0], mouse_i_ptm=mp[0]))
        rec["charge_probes"] = probes
        for pos in hidx:
            q = next((p for p in probes if p["position"] == pos+1 and p["residue"] == "Q"), None)
            kk = next((p for p in probes if p["position"] == pos+1 and p["residue"] == "K"), None)
            if q and kk:
                d = round(kk["human_i_ptm"] - q["human_i_ptm"], 3)
                rec.setdefault("charge_sensitivity", {})["H%d" % (pos+1)] = d
                log("%s H%d charge test: Gln %.3f, Lys %.3f -> %+.3f %s"
                    % (j["slot"], pos+1, q["human_i_ptm"], kk["human_i_ptm"], d,
                       "punishes positive charge" if d <= -0.10
                       else "tolerates positive charge - switch UNSUPPORTED"))
    if hidx:
        step("charge", _charge)

    # --- 2. burial on binding: the variable we should have filtered on
    def _burial():
        rec["burial"] = {("H%d" % (k+1)): v for k, v in burial(pdb0, L, hidx).items()}
        for k, v in rec["burial"].items():
            log("%s %s burial: free %.3f -> bound %.3f (buried %.3f), nearest cation "
                "%s at %s A, screened %s kcal/mol"
                % (j["slot"], k, v["rel_sasa_free"], v["rel_sasa_bound"],
                   v["buried_on_binding"], v["nearest_cation"], v["cation_dist"],
                   v["screened_coulomb_kcal"]))
    if hidx:
        step("burial", _burial)

    # --- 3. does it fold without the target
    def _alone():
        a = predict(seq, "alone")
        if a:
            rec["alone_plddt"] = a[1]
            if "human_plddt" in rec:
                rec["plddt_drop_alone"] = round(rec["human_plddt"] - a[1], 3)
            log("%s alone pLDDT %.3f vs %.3f bound (drop %s)"
                % (j["slot"], a[1], rec.get("human_plddt", 0),
                   rec.get("plddt_drop_alone")))
    step("alone", _alone)

    # --- 4. error bars, last because it is the one we can afford to lose
    def _ensemble():
        for species in ("human", "mouse"):
            vals = []
            for k in range(N_MODELS):
                r = predict(seq, species, model=k)
                if r is None:
                    break
                vals.append(r)
            if len(vals) > 1:
                ip = [v[0] for v in vals]
                rec["%s_i_ptm_mean" % species] = round(float(np.mean(ip)), 3)
                rec["%s_i_ptm_sd" % species] = round(float(np.std(ip)), 3)
                rec["%s_i_ptm_all" % species] = ip
            elif vals:
                rec["%s_models_note" % species] = (
                    "only %d parameter set loaded, no spread available"
                    % n_models_available())
        if "human_i_ptm_sd" in rec:
            log("%s spread over %d models: human %.3f+/-%.3f  mouse %.3f+/-%.3f"
                % (j["slot"], len(rec.get("human_i_ptm_all", [])),
                   rec.get("human_i_ptm_mean", 0), rec.get("human_i_ptm_sd", 0),
                   rec.get("mouse_i_ptm_mean", 0), rec.get("mouse_i_ptm_sd", 0)))
    step("ensemble", _ensemble)
    emit(rec)
log("worker done")
"""

with open("/kaggle/working/tnf_evidence_worker.py", "w") as fh:
    fh.write(textwrap.dedent(WORKER))

subprocess.run("pip install -q biopython", shell=True)
NEEDED = ["/kaggle/working/targets/TNF_site_crop.pdb",
          "/kaggle/working/targets/TNF_mouse_site_crop.pdb",
          "/kaggle/working/construct_map.json",
          "/kaggle/working/construct_map_mouse.json",
          "/kaggle/working/mouse_hotspot.txt"]
missing = [p for p in NEEDED if not os.path.exists(p)]
if missing:
    raise SystemExit("setup cells have not run - missing %s" % missing)

procs = []
for shard, gpu in enumerate(["0", "1"]):
    p = subprocess.Popen([sys.executable, "-u", "/kaggle/working/tnf_evidence_worker.py",
                          gpu, str(shard), str(TIME_BUDGET_MIN), str(N_SWITCH),
                          str(TOTAL_SLOTS), str(N_MODELS)],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         text=True, bufsize=1)
    procs.append(p)
    print("launched evidence worker %d on GPU %s" % (shard, gpu), flush=True)

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
for f in sorted(glob.glob("/kaggle/working/tnf_evidence_w*.jsonl")):
    for line in open(f):
        line = line.strip()
        if line:
            rows.append(json.loads(line))
rows.sort(key=lambda r: r.get("slot", ""))
json.dump(rows, open("/kaggle/working/tnf_evidence.json", "w"), indent=1)
if not rows:
    raise SystemExit("nothing was measured - read the worker log above, do not read an "
                     "empty table as 'the switch does not work'")

# Run 1 printed twenty rows of 0.000 under a green "Successfully ran". A cell that
# measured nothing must say so by failing, not by printing an empty table.
measured = [r for r in rows if "human_i_ptm" in r]
charged = [r for r in rows if r.get("charge_sensitivity")]
if not measured:
    raise SystemExit("no design was scored at all - every baseline failed. Read the "
                     "per-design *_error fields above. Do NOT read this as a result.")
print("\n==== BASELINE AND FOLD WITHOUT THE TARGET ====")
print("%-6s %-8s %8s %8s %8s %8s %8s %-14s"
      % ("slot", "role", "human", "mouse", "pLDDT", "alone", "drop", "spread h"))
for r in rows:
    sd = ("+/-%.3f" % r["human_i_ptm_sd"]) if "human_i_ptm_sd" in r else "-"
    print("%-6s %-8s %8s %8s %8s %8s %8s %-14s"
          % (r.get("slot"), r.get("role", ""),
             ("%.3f" % r["human_i_ptm"]) if "human_i_ptm" in r else "FAILED",
             ("%.3f" % r["mouse_i_ptm"]) if "mouse_i_ptm" in r else "-",
             ("%.3f" % r["human_plddt"]) if "human_plddt" in r else "-",
             ("%.3f" % r["alone_plddt"]) if "alone_plddt" in r else "-",
             ("%+.3f" % r["plddt_drop_alone"]) if "plddt_drop_alone" in r else "-",
             sd))
drops = [r["plddt_drop_alone"] for r in rows if "plddt_drop_alone" in r]
if drops:
    drops.sort()
    print("\npLDDT drop when the target is removed: median %+.3f over %d designs"
          % (drops[len(drops)//2], len(drops)))
    print("A large positive drop means the binder needs the target to hold its shape,")
    print("which is a poor sign for expressing it on its own.")

print("\n==== CHARGE SENSITIVITY: does the position punish a positive charge? ====")
print("A negative number means Lys costs more than Gln at that position, i.e. the")
print("position rejects positive charge and the switch has a reason to fire.")
print("%-6s %-14s %6s %8s %8s %8s  %s"
      % ("slot", "lineage", "his", "Gln", "Lys", "Lys-Gln", "reading"))
supported, unsupported = [], []
for r in rows:
    for key, d in (r.get("charge_sensitivity") or {}).items():
        pos = int(key[1:])
        pr = {p["residue"]: p for p in r.get("charge_probes", []) if p["position"] == pos}
        good = d <= -0.10
        (supported if good else unsupported).append((r["slot"], key, d))
        print("%-6s %-14s %6s %8.3f %8.3f %+8.3f  %s"
              % (r["slot"], r.get("parent", ""), key,
                 pr.get("Q", {}).get("human_i_ptm", 0),
                 pr.get("K", {}).get("human_i_ptm", 0), d,
                 "punishes +" if good else "tolerates + : UNSUPPORTED"))
print("\n%d histidine positions punish a positive charge; %d tolerate it"
      % (len(supported), len(unsupported)))
if unsupported:
    print("switch claim NOT supported at: %s"
          % ", ".join("%s %s" % (s, k) for s, k, _ in unsupported))
by_slot = {}
for s, k, d in supported + unsupported:
    by_slot.setdefault(s, []).append(d <= -0.10)
clean = sorted(s for s, v in by_slot.items() if all(v))
print("designs where EVERY histidine rejects positive charge: %s"
      % (", ".join(clean) if clean else "none"))
if not charged:
    print("\nWARNING: not one charge test completed. Objective 1 is still untested and "
          "nothing below should be read as evidence about the switch.")

print("\n==== BURIAL ON BINDING (the variable we should have filtered on) ====")
print("%-6s %6s %8s %8s %8s %12s %7s %9s"
      % ("slot", "his", "free", "bound", "buried", "cation", "dist", "kcal/mol"))
for r in rows:
    for key, v in (r.get("burial") or {}).items():
        print("%-6s %6s %8.3f %8.3f %8.3f %12s %7s %9s"
              % (r["slot"], key, v["rel_sasa_free"], v["rel_sasa_bound"],
                 v["buried_on_binding"], v["nearest_cation"] or "-",
                 v["cation_dist"] or "-", v["screened_coulomb_kcal"] or "-"))

print("\n===== COPY FROM HERE =====")
for r in rows:
    print(json.dumps(r, separators=(",", ":")))
print("===== COPY TO HERE =====")
