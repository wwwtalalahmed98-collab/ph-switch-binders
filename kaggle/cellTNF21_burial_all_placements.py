# === CHALLENGE 2 (TNF-alpha) : BURIAL AND CHARGE FOR EVERY VERIFIED PLACEMENT ===
#
# The evidence run measured the 12 placements we happened to submit and found that
# burial on binding predicts the charge effect at r = -0.58 while the Cbeta-to-cation
# distance we actually designed with predicts it at r = +0.24. Two submitted designs
# turned out to have ZERO burial, meaning their histidines sit in open solvent when
# bound and there is no mechanism for them to switch at all.
#
# This measures all 35 verified placements across the 14 usable lineages, not just the
# 12 chosen, so the set can be re-ranked on the variable that actually predicts the
# switch rather than the one that does not. No new design work: every sequence here was
# already generated and already passed the both-species test in switch5.
#
# Three defects from the evidence run are fixed here:
#   1. The monomer prediction failed on all 20 designs because protocol="hallucination"
#      needs its own prep_inputs call. Whether our binders fold WITHOUT the target is
#      still unknown, and it is the question the Proteinbase pLDDT numbers raised.
#   2. Cation identity labels read residue numbers straight out of the predicted PDB
#      instead of mapping through the construct map, so the names were in the wrong
#      numbering system. Now mapped, and asserted.
#   3. The run rebuilt its AF2 model 152 times for ~212 predictions, because it
#      alternated human/mouse/monomer per design and each switch forces a rebuild.
#      Work is now done in passes - every human prediction, then every mouse one - so
#      rebuilds drop to roughly one per distinct binder length per species.
#
# Error bars are dropped: the evidence run showed only 2 of 5 parameter sets ever load,
# and we already know the spread is about 0.02 on human and 0.03 on mouse.
import glob, json, os, subprocess, sys, textwrap, time

TIME_BUDGET_MIN = 200

WORKER = r"""
import glob, json, os, sys, time
GPU, SHARD, BUDGET = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
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
EPS_EFF, DEBYE = 30.0, 7.9
ANCHORS = {("A", 166), ("B", 108)}   # the conserved, exposed cations we designed against

src = glob.glob("/kaggle/input/**/tnf_switch5.json", recursive=True)
if not src:
    raise SystemExit("tnf_switch5.json not found - attach tnf_switch5_repair_fix")
recs = json.load(open(src[0]))
if len(recs) < 26:
    raise SystemExit("tnf_switch5.json has only %d records, expected 29; it looks "
                     "truncated (Kaggle's preview cuts at 102400 bytes and closes the "
                     "JSON brackets). Refusing to select from it." % len(recs))

jobs = []
for r in recs:
    for v in r.get("variants", []):
        if not v.get("passes"):
            continue
        jobs.append(dict(parent=r["name"], his=v["his_substitutions"], seq=v["sequence"],
                         kind=v["kind"], n_his=v["n_his"], n_anchor=v["n_anchor"],
                         switch_score=v["switch_score"], ref_h=v["human_i_ptm"],
                         ref_m=v["mouse_i_ptm"], ref_mo=len(v.get("mouse_on_site") or []),
                         ctl=r["control"]["sequence"]))
jobs.sort(key=lambda j: (len(j["seq"]), j["parent"], j["his"]))
mine = [j for k, j in enumerate(jobs) if k % 2 == SHARD]
if SHARD == 0:
    log("%d verified placements across %d lineages"
        % (len(jobs), len({j["parent"] for j in jobs})))
log("this worker takes %d placements" % len(mine))

from colabdesign import mk_afdesign_model, clear_mem
_c = {}

def _model(species, L):
    key = (species, L)
    if _c.get("key") != key:
        clear_mem()
        if species == "alone":
            # The evidence run omitted this prep_inputs call and every monomer
            # prediction raised AssertionError: inputs not set {'length': L}.
            af = mk_afdesign_model(protocol="hallucination", data_dir="/kaggle/working")
            af.prep_inputs(length=L)
        else:
            af = mk_afdesign_model(protocol="binder", data_dir="/kaggle/working")
            af.prep_inputs(pdb_filename=HT if species == "human" else MT, chain="A,B",
                           binder_len=L,
                           hotspot=H_HOTSPOT if species == "human" else M_HOTSPOT)
        _c["key"], _c["af"] = key, af
        log("built %s model for length %d" % (species, L))
    return _c["af"]

def predict(seq, species, out_pdb=None):
    af = _model(species, len(seq))
    try:
        af.predict(seq=seq, num_recycles=3, verbose=False)
    except TypeError:
        af.predict(seq=seq, verbose=False)
    lg = af.aux["log"]
    if out_pdb:
        af.save_pdb(out_pdb)
    return (round(float(lg.get("i_ptm", lg.get("ptm", 0.0))), 3),
            round(float(lg.get("plddt", 0.0)), 3))

# ---- burial, with cation names mapped through the construct map ------------------
def dedupe(chain):
    res, seen = [], set()
    for r in list(chain):
        if r.id[0] != " " or r.id[1] in seen:
            chain.detach_child(r.id); continue
        seen.add(r.id[1]); res.append(r)
    return res

def burial(pdb_path, blen, his_idx):
    mdl = PDBParser(QUIET=True).get_structure("c", pdb_path)[0]
    kept = {ch.id: dedupe(ch) for ch in mdl}
    bnd_id = min(kept, key=lambda k: len(kept[k]))
    bres = kept[bnd_id]
    if len(bres) != blen:
        raise ValueError("binder chain %s has %d residues, expected %d"
                         % (bnd_id, len(bres), blen))
    # target residues in file order; colabdesign writes the two target chains merged,
    # so the order must line up with construct_map or every label would be wrong.
    tres = []
    for cid in sorted(k for k in kept if k != bnd_id):
        tres.extend(kept[cid])
    if len(tres) != len(CMAP):
        raise ValueError("target has %d residues, construct map has %d"
                         % (len(tres), len(CMAP)))
    seen_aa = "".join(protein_letters_3to1.get(r.get_resname(), "X") for r in tres)
    want_aa = "".join(o["aa"] for o in CMAP)
    if seen_aa != want_aa:
        raise ValueError("target sequence does not match construct_map")
    cats = []
    for k, r in enumerate(tres):
        for an in CAT_N.get(r.get_resname(), ()):
            if an in r:
                cats.append((CMAP[k]["chain"], CMAP[k]["uniprot"], an,
                             r[an].get_coord()))
    if not cats:
        raise ValueError("no Arg/Lys nitrogens on the target")
    ShrakeRupley().compute(mdl, level="R")
    bound = {}
    for k, r in enumerate(bres):
        aa = protein_letters_3to1.get(r.get_resname())
        if aa: bound[k] = r.sasa / MAXS.get(aa, 200)
    import copy
    solo = copy.deepcopy(mdl[bnd_id])
    sres = [r for r in solo if r.id[0] == " "]
    ShrakeRupley().compute(solo, level="R")
    free = {}
    for k, r in enumerate(sres):
        aa = protein_letters_3to1.get(r.get_resname())
        if aa: free[k] = r.sasa / MAXS.get(aa, 200)
    out = {}
    for k in his_idx:
        r = bres[k]
        best = None
        for an in ("ND1", "NE2"):
            if an not in r: continue
            for cid, unum, catn, xyz in cats:
                d = float(np.linalg.norm(r[an].get_coord() - xyz))
                if best is None or d < best[0]:
                    best = (d, cid, unum, catn)
        ddg = None
        if best is not None:
            ddg = round(332.0 / (EPS_EFF * best[0]) * np.exp(-best[0] / DEBYE), 3)
        out[k] = dict(rel_sasa_free=round(free.get(k, -1), 3),
                      rel_sasa_bound=round(bound.get(k, -1), 3),
                      buried_on_binding=round(free.get(k, 0) - bound.get(k, 0), 3),
                      cation="%s%d%s" % (best[1], best[2], best[3]) if best else None,
                      cation_is_anchor=bool(best and (best[1], best[2]) in ANCHORS),
                      cation_dist=round(best[0], 2) if best else None,
                      clash=bool(best and best[0] < 2.5),
                      screened_coulomb_kcal=ddg)
    return out

os.makedirs("/kaggle/working/bur_pdb", exist_ok=True)

# Build the full prediction list up front, then run it in passes grouped by
# (species, length). The evidence run rebuilt the model 152 times for ~212
# predictions because it interleaved species per design; this does not.
for j in mine:
    j["hidx"] = [int(s[1:-1]) - 1 for s in j["his"].split(",")]
    j["key"] = "%s__%s" % (j["parent"], j["his"].replace(",", "_"))
    j["mut"] = {}
    for pos in j["hidx"]:
        if j["seq"][pos] != "H":
            log("%s expected H at %d, found %s - skipping its probes"
                % (j["key"], pos + 1, j["seq"][pos]))
            continue
        for aa in ("Q", "K"):
            j["mut"][(pos, aa)] = j["seq"][:pos] + aa + j["seq"][pos+1:]

def pass_species(species):
    todo = []
    for j in mine:
        todo.append((j, "base", None, j["seq"]))
        for (pos, aa), s in j["mut"].items():
            todo.append((j, "mut", (pos, aa), s))
    todo.sort(key=lambda x: len(x[3]))
    for j, what, tag, s in todo:
        if (time.time()-t0)/60 > BUDGET:
            log("time budget reached during the %s pass" % species); return
        try:
            out = ("/kaggle/working/bur_pdb/%s.pdb" % j["key"]) if (
                what == "base" and species == "human") else None
            r = predict(s, species, out_pdb=out)
            if what == "base":
                j.setdefault("base", {})[species] = r
            else:
                j.setdefault("probe", {})[(species,) + tag] = r[0]
        except Exception as e:
            j.setdefault("err", []).append("%s %s %s" % (species, what,
                                                         type(e).__name__))
            log("%s %s %s FAILED %s %s" % (j["key"], species, what,
                                           type(e).__name__, str(e)[:100]))

for sp in ("human", "mouse"):
    log("---- %s pass ----" % sp)
    pass_species(sp)

log("---- monomer pass ----")
for j in sorted(mine, key=lambda x: len(x["seq"])):
    if (time.time()-t0)/60 > BUDGET: break
    try:
        j["alone"] = predict(j["seq"], "alone")[1]
    except Exception as e:
        j.setdefault("err", []).append("alone " + type(e).__name__)
        log("%s alone FAILED %s %s" % (j["key"], type(e).__name__, str(e)[:120]))

log("---- burial pass (cpu) ----")
rec_path = "/kaggle/working/tnf_burial_w%d.jsonl" % SHARD
for j in mine:
    rec = dict(parent=j["parent"], his=j["his"], kind=j["kind"], n_his=j["n_his"],
               n_anchor=j["n_anchor"], switch_score=j["switch_score"],
               length=len(j["seq"]), ref_human=j["ref_h"], ref_mouse=j["ref_m"],
               ref_mouse_on_site=j["ref_mo"], errors=j.get("err", []))
    b = j.get("base", {})
    if "human" in b: rec["human_i_ptm"], rec["human_plddt"] = b["human"]
    if "mouse" in b: rec["mouse_i_ptm"], rec["mouse_plddt"] = b["mouse"]
    if "alone" in j:
        rec["alone_plddt"] = j["alone"]
        if "human_plddt" in rec:
            rec["plddt_drop_alone"] = round(rec["human_plddt"] - j["alone"], 3)
    pr = j.get("probe", {})
    chg = {}
    for pos in j["hidx"]:
        q, k = pr.get(("human", pos, "Q")), pr.get(("human", pos, "K"))
        if q is not None and k is not None:
            chg["H%d" % (pos+1)] = dict(gln=q, lys=k, lys_minus_gln=round(k-q, 3))
        qm, km = pr.get(("mouse", pos, "Q")), pr.get(("mouse", pos, "K"))
        if qm is not None and km is not None:
            chg.setdefault("H%d" % (pos+1), {})["mouse_lys_minus_gln"] = round(km-qm, 3)
    if chg: rec["charge"] = chg
    try:
        pdb0 = "/kaggle/working/bur_pdb/%s.pdb" % j["key"]
        rec["burial"] = {("H%d" % (k+1)): v
                         for k, v in burial(pdb0, len(j["seq"]), j["hidx"]).items()}
    except Exception as e:
        rec["burial_error"] = "%s %s" % (type(e).__name__, str(e)[:160])
    bs = [v["buried_on_binding"] for v in (rec.get("burial") or {}).values()]
    rec["min_burial"] = min(bs) if bs else None
    rec["sum_burial"] = round(sum(bs), 3) if bs else None
    cs = [v["lys_minus_gln"] for v in chg.values() if "lys_minus_gln" in v]
    rec["worst_charge"] = max(cs) if cs else None     # least negative = weakest link
    with open(rec_path, "a") as fh:
        fh.write(json.dumps(rec) + chr(10))
    log("%-28s %-12s burial min %-6s sum %-6s charge worst %-7s human %s mouse %s"
        % (j["parent"], j["his"], rec["min_burial"], rec["sum_burial"],
           rec["worst_charge"], rec.get("human_i_ptm"), rec.get("mouse_i_ptm")))
log("worker done")
"""

with open("/kaggle/working/tnf_burial_worker.py", "w") as fh:
    fh.write(textwrap.dedent(WORKER))

subprocess.run("pip install -q biopython", shell=True)
NEEDED = ["/kaggle/working/targets/TNF_site_crop.pdb",
          "/kaggle/working/targets/TNF_mouse_site_crop.pdb",
          "/kaggle/working/construct_map.json",
          "/kaggle/working/mouse_hotspot.txt"]
missing = [p for p in NEEDED if not os.path.exists(p)]
if missing:
    raise SystemExit("setup cells have not run - missing %s" % missing)

procs = []
for shard, gpu in enumerate(["0", "1"]):
    p = subprocess.Popen([sys.executable, "-u", "/kaggle/working/tnf_burial_worker.py",
                          gpu, str(shard), str(TIME_BUDGET_MIN)],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         text=True, bufsize=1)
    procs.append(p)
    print("launched burial worker %d on GPU %s" % (shard, gpu), flush=True)

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
for f in sorted(glob.glob("/kaggle/working/tnf_burial_w*.jsonl")):
    for line in open(f):
        line = line.strip()
        if line:
            rows.append(json.loads(line))
json.dump(rows, open("/kaggle/working/tnf_burial.json", "w"), indent=1)
if not rows:
    raise SystemExit("nothing was measured - read the worker log above. Do NOT read an "
                     "empty table as 'no placement is buried'.")
scored = [r for r in rows if r.get("min_burial") is not None]
if not scored:
    raise SystemExit("no placement produced a burial number - every burial call failed. "
                     "Check burial_error above; this is not a result.")

# rank on the variable that actually predicted the charge effect
scored.sort(key=lambda r: (-(r["min_burial"] or 0), r.get("worst_charge") or 0))
print("\n==== ALL VERIFIED PLACEMENTS, RANKED BY WORST-HISTIDINE BURIAL ====")
print("The evidence run found burial predicts the charge effect at r = -0.58 and the")
print("Cbeta-to-cation distance we designed with predicts it at r = +0.24.")
print("%-14s %-12s %3s %7s %7s %8s %7s %7s %5s"
      % ("lineage", "histidines", "nH", "burial", "sum", "charge", "human", "mouse", "site"))
for r in scored:
    print("%-14s %-12s %3d %7s %7s %8s %7s %7s %5d"
          % (r["parent"], r["his"], r["n_his"], r["min_burial"], r["sum_burial"],
             r["worst_charge"], r.get("human_i_ptm"), r.get("mouse_i_ptm"),
             r["ref_mouse_on_site"]))

dead = [r for r in scored if (r["min_burial"] or 0) < 0.02]
print("\n%d of %d placements have a histidine with essentially no burial (<0.02):"
      % (len(dead), len(scored)))
for r in dead:
    print("  %-14s %-12s  no mechanism for that histidine to switch"
          % (r["parent"], r["his"]))

best = {}
for r in scored:
    if r["parent"] not in best:
        best[r["parent"]] = r
print("\n==== BEST PLACEMENT PER LINEAGE, BY BURIAL ====")
print("%-14s %-12s %3s %7s %8s %7s %7s  %s"
      % ("lineage", "histidines", "nH", "burial", "charge", "human", "mouse", "vs submitted"))
for p, r in sorted(best.items(), key=lambda kv: -(kv[1]["min_burial"] or 0)):
    print("%-14s %-12s %3d %7s %8s %7s %7s"
          % (p, r["his"], r["n_his"], r["min_burial"], r["worst_charge"],
             r.get("human_i_ptm"), r.get("mouse_i_ptm")))

drops = [r["plddt_drop_alone"] for r in rows if "plddt_drop_alone" in r]
if drops:
    drops.sort()
    print("\n==== DOES THE BINDER FOLD WITHOUT THE TARGET? ====")
    print("pLDDT bound minus pLDDT alone, over %d placements: median %+.3f, "
          "range %+.3f to %+.3f" % (len(drops), drops[len(drops)//2], drops[0], drops[-1]))
    print("A large positive number means the binder needs the target to hold its shape.")
else:
    print("\nthe monomer pass produced nothing - whether these fold alone is STILL unknown")

clashes = [(r["parent"], r["his"], k, v["cation_dist"])
           for r in rows for k, v in (r.get("burial") or {}).items() if v.get("clash")]
if clashes:
    print("\n%d histidine(s) sit closer than 2.5 A to a cation, which is a clash rather "
          "than a contact - the predicted structure is wrong there:" % len(clashes))
    for c in clashes:
        print("  %-14s %-12s %s at %.2f A" % c)

anc = sum(1 for r in rows for v in (r.get("burial") or {}).values()
          if v.get("cation_is_anchor"))
tot = sum(1 for r in rows for v in (r.get("burial") or {}).values())
print("\n%d of %d histidines sit nearest to one of the two CONSERVED anchors "
      "(K166A / R108B); the rest face some other cation." % (anc, tot))

print("\n===== COPY FROM HERE =====")
for r in rows:
    print(json.dumps(r, separators=(",", ":")))
print("===== COPY TO HERE =====")
