# === CHALLENGE 2 (TNF-alpha) : CONFIRMATION OF THE SUBMITTED DESIGNS ===
#
# The switch stage scored each variant by prediction but did not save the predicted
# complex, so two things were asserted rather than measured:
#
#   1. that a switch variant still binds the SAME epitope as its parent. Retaining
#      i_ptm does not prove that - a histidine could shift the binding mode and keep a
#      similar score. On EGFR the designs drifted off the intended patch while looking
#      fine, and that was only found at the end.
#   2. that each designed histidine actually ends up facing the cation it was placed
#      against. Placement was computed on the PARENT backbone; whether the geometry
#      survives re-folding with three residues changed is a separate question, and it
#      is the one the whole acid-OFF hypothesis rests on.
#
# This re-predicts the 16 submitted sequences WITH structures saved and measures both,
# plus the control for each pair, so any epitope shift can be attributed to the
# histidines rather than to prediction variance between two runs.
import glob, json, os, subprocess, sys, textwrap, time

TIME_BUDGET_MIN = 60

WORKER = r"""
import glob, json, os, sys, time
from collections import Counter
GPU, SHARD, BUDGET = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
os.environ["CUDA_VISIBLE_DEVICES"] = GPU
t0 = time.time()
def log(m): print("[w%d] %s [%.1f min]" % (SHARD, m, (time.time()-t0)/60), flush=True)

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1

TARGET = "/kaggle/working/targets/TNF_site_crop.pdb"
CMAP = json.load(open("/kaggle/working/construct_map.json"))["order"]
N = len(CMAP)
SITE = {95,105,106,107,108,109,110,111,112,155,159,161,162,163,164,165,166,167,168,
        201,203,207,220,221,223}
CAT_N = {"ARG": ("NE","NH1","NH2"), "LYS": ("NZ",)}
# index into the construct for each partner named in the switch record, e.g. "R108B"
PIDX = {}
for i, o in enumerate(CMAP):
    PIDX["%s%d%s" % (o["aa"], o["uniprot"], o["chain"])] = i

def vcb(r):
    if "CB" in r: return r["CB"].get_coord()
    try: n,a,c = r["N"].get_coord(), r["CA"].get_coord(), r["C"].get_coord()
    except KeyError: return None
    b, cc = a-n, c-a
    return -0.58273431*np.cross(b,cc) + 0.56802827*b - 0.54067466*cc + a

def trim(chain, n, what):
    res = [r for r in chain if r.id[0] == " "]
    if len(res) > n and len(res) % n == 0:
        for r in res[n:]:
            chain.detach_child(r.id)
        res = res[:n]
    if len(res) != n:
        raise ValueError("%s chain has %d residues, expected %d" % (what, len(res), n))
    return res

def measure(pdb_path, blen):
    mdl = PDBParser(QUIET=True).get_structure("c", pdb_path)[0]
    sizes = Counter({ch.id: sum(1 for r in ch if r.id[0]==" ") for ch in mdl})
    tgt, bnd = max(sizes, key=lambda k: sizes[k]), min(sizes, key=lambda k: sizes[k])
    tres = trim(mdl[tgt], N, "target")
    bres = trim(mdl[bnd], blen, "binder")
    seen = "".join(protein_letters_3to1.get(r.get_resname(), "X") for r in tres)
    if seen != "".join(o["aa"] for o in CMAP):
        raise ValueError("target sequence does not match construct_map")
    batoms = np.array([a.get_coord() for r in bres for a in r])
    epi = []
    for i, r in enumerate(tres):
        ra = np.array([a.get_coord() for a in r])
        if np.min(np.linalg.norm(ra[:,None,:]-batoms[None,:,:], axis=-1)) <= 5.0:
            epi.append(CMAP[i]["uniprot"])
    return sorted(set(epi)), tres, bres

def his_contacts(tres, bres, his_pos, partners):
    # closest Cbeta->(cation N) distance for each designed histidine, against the
    # partners it was placed to face. Same measurement the placement used.
    out = []
    for p in his_pos:
        r = bres[p-1]
        cb = vcb(r)
        aa = protein_letters_3to1.get(r.get_resname(), "X")
        best, who = 99.0, "-"
        for name in partners:
            tr = tres[PIDX[name]]
            at = [tr[a].get_coord() for a in CAT_N.get(tr.get_resname(), ()) if a in tr]
            for x in at:
                d = float(np.linalg.norm(cb - x))
                if d < best:
                    best, who = d, name
        out.append((p, aa, who, round(best, 2)))
    return out

from colabdesign import mk_afdesign_model, clear_mem
HOTSPOT = "A86,A87,A88,A90,A91,A92,B29,B32,B33,B34,B35,B36"
_c = {}
def predict(seq, out_pdb):
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

rj = glob.glob("/kaggle/input/**/tnf_designs_final.json", recursive=True)
if not rj:
    raise SystemExit("tnf_designs_final.json not found - attach the tnf_switch output")
recs = {r["name"]: r for r in json.load(open(rj[0]))}
usable = sorted([r for r in recs.values()
                 if r["role"] == "switch" and r["i_ptm"] >= 0.45
                 and r["iptm_retained"] >= 0.70],
                key=lambda r: -r["i_ptm"])
# pairs stay on one worker: the control shares its switch's length, so the model is
# built once per pair instead of twice
mine = [r for i, r in enumerate(usable) if i % 2 == SHARD]
if SHARD == 0:
    log("%d records, %d usable switches to confirm" % (len(recs), len(usable)))
log("this worker takes %d: %s" % (len(mine), ", ".join(r["name"] for r in mine)))

os.makedirs("/kaggle/working/confirm", exist_ok=True)
rec_path = "/kaggle/working/confirm_w%d.jsonl" % SHARD
for sw in mine:
    if (time.time()-t0)/60 > BUDGET:
        log("time budget reached"); break
    ctl = recs.get(sw["source_backbone"] + "_ctl")
    parent_epi = set(sw["epitope_residues"])
    his_pos = [int(x[1:-1]) for x in sw["his_substitutions"].split(",")]
    row = dict(name=sw["name"], source=sw["source_backbone"], length=sw["length"],
               his_substitutions=sw["his_substitutions"],
               intended_partners=sw["his_partners"],
               stage3_i_ptm=sw["i_ptm"], parent_epitope=sorted(parent_epi))
    try:
        for tag, r in (("sw", sw), ("ctl", ctl)):
            if r is None:
                continue
            out = "/kaggle/working/confirm/%s_%s.pdb" % (sw["source_backbone"], tag)
            ip, pl = predict(r["sequence"], out)
            epi, tres, bres = measure(out, r["length"])
            keep = sorted(set(epi) & parent_epi)
            row[tag] = dict(i_ptm=ip, plddt=pl, epitope=epi,
                            on_site=sorted(set(epi) & SITE),
                            kept_parent=len(keep),
                            frac_parent=round(len(keep)/max(len(parent_epi), 1), 3))
            if tag == "sw":
                row["his_geometry"] = his_contacts(tres, bres, his_pos,
                                                   sw["his_partners"])
        g = row.get("his_geometry", [])
        near = sum(1 for _, _, _, d in g if d <= 8.8)
        still_his = sum(1 for _, aa, _, _ in g if aa == "H")
        log("%-14s i_ptm %.3f (stage3 %.3f)  epitope %d, %d of %d parent (%.0f%%), "
            "%d on site | His kept %d/%d, facing partner %d/%d, Cb %s"
            % (sw["name"], row["sw"]["i_ptm"], sw["i_ptm"], len(row["sw"]["epitope"]),
               row["sw"]["kept_parent"], len(parent_epi),
               100*row["sw"]["frac_parent"], len(row["sw"]["on_site"]),
               still_his, len(g), near, len(g),
               ", ".join("%s%d->%s %.1f" % (aa, p, who, d) for p, aa, who, d in g)))
        if "ctl" in row:
            log("   %-11s control i_ptm %.3f  epitope %d, %d of %d parent (%.0f%%)"
                % (sw["source_backbone"], row["ctl"]["i_ptm"],
                   len(row["ctl"]["epitope"]), row["ctl"]["kept_parent"],
                   len(parent_epi), 100*row["ctl"]["frac_parent"]))
    except Exception as e:
        row["error"] = "%s %s" % (type(e).__name__, str(e)[:200])
        log("%s FAILED %s" % (sw["name"], row["error"]))
    with open(rec_path, "a") as fh:
        fh.write(json.dumps(row) + chr(10))
log("worker done")
"""

with open("/kaggle/working/confirm_worker.py", "w") as fh:
    fh.write(textwrap.dedent(WORKER))

procs = []
for shard, gpu in enumerate(["0", "1"]):
    p = subprocess.Popen([sys.executable, "-u", "/kaggle/working/confirm_worker.py",
                          gpu, str(shard), str(TIME_BUDGET_MIN)],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         text=True, bufsize=1)
    procs.append(p)
    print("launched confirm worker %d on GPU %s" % (shard, gpu), flush=True)

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
for f in sorted(glob.glob("/kaggle/working/confirm_w*.jsonl")):
    for line in open(f):
        line = line.strip()
        if line:
            rows.append(json.loads(line))
rows.sort(key=lambda r: -r.get("sw", {}).get("i_ptm", 0))
json.dump(rows, open("/kaggle/working/tnf_confirm.json", "w"), indent=1)

print("\n==== CONFIRMATION OF THE 8 SUBMITTED SWITCHES ====")
print("%-14s %7s %7s %8s %8s %6s %6s %-5s"
      % ("design", "i_ptm", "stage3", "epi kept", "ctl kept", "site", "His", "faces"))
for r in rows:
    if "sw" not in r:
        print("%-14s  ERROR %s" % (r["name"], r.get("error", "")))
        continue
    g = r.get("his_geometry", [])
    print("%-14s %7.3f %7.3f %7.0f%% %7s %6d %5d/%d %5d"
          % (r["name"], r["sw"]["i_ptm"], r["stage3_i_ptm"],
             100*r["sw"]["frac_parent"],
             ("%.0f%%" % (100*r["ctl"]["frac_parent"])) if "ctl" in r else "-",
             len(r["sw"]["on_site"]),
             sum(1 for _, aa, _, _ in g if aa == "H"), len(g),
             sum(1 for _, _, _, d in g if d <= 8.8)))

ok = [r for r in rows if "sw" in r and r["sw"]["frac_parent"] >= 0.60
      and len(r["sw"]["on_site"]) >= 5]
geo = [r for r in rows if "sw" in r
       and all(d <= 8.8 for _, _, _, d in r.get("his_geometry", []))]
print("\n%d of %d hold >=60%% of the parent epitope AND 5+ on-site contacts"
      % (len(ok), len(rows)))
print("%d of %d have EVERY designed histidine still within reach of its cation"
      % (len(geo), len(rows)))
print("\n===== COPY FROM HERE =====")
for r in rows:
    print(json.dumps(r, separators=(",", ":")))
print("===== COPY TO HERE =====")
