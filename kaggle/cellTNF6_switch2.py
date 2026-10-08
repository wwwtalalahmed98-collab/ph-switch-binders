# === CHALLENGE 2 (TNF-alpha) : acid-OFF SWITCH PLACEMENT + LIABILITY REPAIR, ROUND 2 ===
#
# Takes the 27 round-1 designs, keeps the ones that both bind and bind the RIGHT
# place, puts histidines where protonation at pH 6.0 will cost the interface, repairs
# sequence liabilities around them, and re-predicts everything so that every step's
# cost is measured rather than assumed.
#
# Three things carried over from Challenge 1, each because it cost something there:
#
#   1. Nothing is trusted to a residue NUMBER. AfDesign concatenates the two target
#      chains into one and renumbers from 1, so position in construct_map.json IS the
#      output residue index. A three-residue numbering shift silently emptied the
#      switch-site list for two whole rounds on EGFR and returned no error. Here
#      every mapping is by order, with assertions.
#   2. The partner whitelist is NOT restricted to the two anchors chosen in Phase 1.
#      On EGFR the whitelist was set to four carboxylates when the same rule applied
#      across the construct found seven, and widening it was what produced working
#      switches. So: every exposed Arg/Lys in the crop is a candidate, and the two
#      conserved anchors are flagged rather than privileged.
#   3. Wholesale redesign is not used. ProteinMPNN destroyed 12 of 14 interfaces on
#      EGFR. Liabilities are repaired residue by residue, histidines protected, and
#      the repair is kept only if it holds its i_ptm.
import glob, json, os, subprocess, sys, textwrap, time

TIME_BUDGET_MIN = 60           # round 2 is ~half the designs and the stage is cheap
MIN_IPTM = 0.45                # round-1 re-scored i_ptm
MIN_ONSITE = 4                 # contacts on the intended 25-residue site
KEEP_FRAC = 0.70               # a switch must hold this much of the parent i_ptm
REPAIR_KEEP = 0.85

WORKER = r"""
import glob, json, os, sys, time
from collections import Counter
GPU, SHARD, BUDGET = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
os.environ["CUDA_VISIBLE_DEVICES"] = GPU          # must precede the jax import
MIN_IPTM, MIN_ONSITE = float(sys.argv[4]), int(sys.argv[5])
KEEP_FRAC, REPAIR_KEEP = float(sys.argv[6]), float(sys.argv[7])

t0 = time.time()
def log(m): print("[w%d] %s [%.1f min]" % (SHARD, m, (time.time()-t0)/60), flush=True)

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import protein_letters_3to1

TARGET = "/kaggle/working/targets/TNF_site_crop.pdb"
CMAP = json.load(open("/kaggle/working/construct_map.json"))["order"]
N = len(CMAP)
SITE = {95,105,106,107,108,109,110,111,112,155,159,161,162,163,164,165,166,167,168,
        201,203,207,220,221,223}
ANCHORS = {("A",166), ("B",108)}      # conserved, exposed, from PHASE1_REPORT

# Constants from tnf/scripts/acidoff_switch.py, calibrated on 17 structures (278
# His-carboxylate and 279 His-cation closest approaches) and validated against the
# FcRn/IgG acid-ON switch. The window is WIDER than Challenge 1's 3.5-7.0: that one
# was set by eye and would have missed both known FcRn switch residues.
CB_MIN, CB_MAX, COS_MIN = 3.0, 8.8, -0.30
CAT_N = {"ARG": ("NE","NH1","NH2"), "LYS": ("NZ",)}
NOMUT = {"PRO","GLY","CYS"}
# A cationic position is the WORST acid-OFF site, not the best: swapping Arg/Lys for
# His keeps a positive charge at pH 6.0 and loses it at 7.4, which is backwards.
PREF = {"S":1.3,"T":1.3,"N":1.2,"Q":1.2,"A":1.1,"V":0.9,"L":0.8,"I":0.8,"M":0.7,
        "Y":0.6,"F":0.6,"W":0.4,"D":0.5,"E":0.5,"K":0.0,"R":0.0,"H":0.0}
MAXS = {"A":129,"R":274,"N":195,"D":193,"C":167,"E":223,"Q":225,"G":104,"H":224,
        "I":197,"L":201,"K":236,"M":224,"F":240,"P":159,"S":155,"T":172,"W":285,
        "Y":263,"V":174}
MIN_PARTNER_EXPOSURE = 0.20

def vcb(r):
    if "CB" in r: return r["CB"].get_coord()
    try: n,a,c = r["N"].get_coord(), r["CA"].get_coord(), r["C"].get_coord()
    except KeyError: return None
    b, cc = a-n, c-a
    return -0.58273431*np.cross(b,cc) + 0.56802827*b - 0.54067466*cc + a

def rsasa_chain(entity):
    ShrakeRupley().compute(entity, level="R")
    out = {}
    for r in entity.get_residues():
        aa = protein_letters_3to1.get(r.get_resname())
        if r.id[0] == " " and aa:
            out[r.id[1]] = r.sasa / MAXS.get(aa, 200)
    return out

# Free-state exposure of every target residue, measured on the cropped target ITSELF
# rather than on a chain pulled out of the complex. The complex file duplicates its
# target residues, and computing SASA over the duplicate would report everything as
# buried. Keyed by (chain, pdb number), then joined to construct order through CMAP.
_free = {}
for _ch in PDBParser(QUIET=True).get_structure("t", TARGET)[0]:
    _one = PDBParser(QUIET=True).get_structure("x", TARGET)[0][_ch.id]
    for _k, _v in rsasa_chain(_one).items():
        _free[(_ch.id, _k)] = _v
TFREE = [_free[(o["chain"], o["pdb"])] for o in CMAP]

# Every exposed Arg/Lys in the crop is a candidate partner; anchors merely flagged.
CATIONS = []
for i, o in enumerate(CMAP):
    if o["aa"] in ("R","K") and TFREE[i] >= MIN_PARTNER_EXPOSURE:
        CATIONS.append(dict(idx=i, chain=o["chain"], uni=o["uniprot"], aa=o["aa"],
                            anchor=(o["chain"], o["uniprot"]) in ANCHORS))
if SHARD == 0:
    log("construct %d residues; %d exposed cations (%d conserved anchors): %s"
        % (N, len(CATIONS), sum(c["anchor"] for c in CATIONS),
           ", ".join("%s%d%s%s" % (c["aa"], c["uni"], c["chain"],
                                   "*" if c["anchor"] else "") for c in CATIONS)))
got = {(c["chain"], c["uni"]) for c in CATIONS if c["anchor"]}
assert got == ANCHORS, "a conserved anchor is missing or buried: %s" % (ANCHORS - got)

# Cut a chain back to n residues when save_pdb has written it more than once.
#
# The smoke test recorded chains {'A': 608, 'B': 140} for a 304-residue target and a
# 70-residue binder: EVERY chain is written twice. Generation was unaffected because
# it only measured distances, and a duplicate at identical coordinates does not change
# a minimum distance. Solvent accessibility is a different matter - a chain
# overlapping itself reads as fully buried, which would have silently found zero
# histidine sites and produced a clean-looking "no switch possible" result. So the
# duplication is removed here, for BOTH chains, and anything that is not an exact
# multiple is an error rather than a guess.
def trim(chain, n, what):
    res = [r for r in chain if r.id[0] == " "]
    if len(res) > n and len(res) % n == 0:
        for r in res[n:]:
            chain.detach_child(r.id)
        res = res[:n]
    if len(res) != n:
        raise ValueError("%s chain has %d residues, expected %d" % (what, len(res), n))
    return res

def load(pdb_path, blen):
    mdl = PDBParser(QUIET=True).get_structure("c", pdb_path)[0]
    sizes = Counter({ch.id: sum(1 for r in ch if r.id[0]==" ") for ch in mdl})
    tgt, bnd = max(sizes, key=lambda k: sizes[k]), min(sizes, key=lambda k: sizes[k])
    tres = trim(mdl[tgt], N, "target")
    bres = trim(mdl[bnd], blen, "binder")
    seen = "".join(protein_letters_3to1.get(r.get_resname(), "X") for r in tres)
    want = "".join(o["aa"] for o in CMAP)
    if seen != want:
        raise ValueError("target sequence does not match construct_map")
    return mdl, tgt, bnd, tres, bres

def analyse(pdb_path, seq):
    mdl, tgt, bnd, tres, bres = load(pdb_path, len(seq))
    # free-state binder exposure, on its own freshly de-duplicated copy
    _m2, _t2, _b2, _, _br2 = load(pdb_path, len(seq))
    _b = rsasa_chain(_m2[_b2])
    bfree = {i: _b.get(r.id[1], 0.0) for i, r in enumerate(_br2)}
    ShrakeRupley().compute(mdl, level="R")
    tcplx = [r.sasa / MAXS.get(protein_letters_3to1.get(r.get_resname(), "A"), 200)
             for r in tres]

    batoms = np.array([a.get_coord() for r in bres for a in r])
    epitope = []
    for i, r in enumerate(tres):
        ra = np.array([a.get_coord() for a in r])
        if np.min(np.linalg.norm(ra[:,None,:]-batoms[None,:,:], axis=-1)) <= 5.0:
            epitope.append(CMAP[i]["uniprot"])

    partners = []
    for c in CATIONS:
        r = tres[c["idx"]]
        at = [r[a].get_coord() for a in CAT_N[r.get_resname()] if a in r]
        if at:
            p = dict(c)
            p["at"] = at
            p["buried"] = max(TFREE[c["idx"]] - tcplx[c["idx"]], 0.0)
            partners.append(p)

    props = []
    for i, r in enumerate(bres):
        if "CA" not in r or r.get_resname() in NOMUT:
            continue
        aa = protein_letters_3to1.get(r.get_resname())
        if aa is None or PREF.get(aa, 0.9) <= 0 or bfree[i] < 0.15:
            continue
        cb = vcb(r)
        if cb is None:
            continue
        d = cb - r["CA"].get_coord()
        nn = float(np.linalg.norm(d))
        if nn < 1e-6:
            continue
        d = d / nn
        for p in partners:
            ds = [float(np.linalg.norm(cb - x)) for x in p["at"]]
            dcb = min(ds)
            if not (CB_MIN <= dcb <= CB_MAX):
                continue
            to = p["at"][int(np.argmin(ds))] - cb
            to = to / np.linalg.norm(to)
            if float(np.dot(d, to)) < COS_MIN:
                continue
            # Repulsion grows monotonically as the charges approach - unlike the
            # attractive acid-ON case there is no optimum distance - and burial
            # removes the solvent screening that would otherwise soften it.
            q = (CB_MAX - dcb) / (CB_MAX - CB_MIN) * (1 + 2.0 * p["buried"])
            props.append(dict(i=i, aa=aa, d=round(dcb,2), qgeom=round(q,3),
                              q=round(q * PREF.get(aa,0.9), 3), anchor=p["anchor"],
                              partner="%s%d%s" % (p["aa"], p["uni"], p["chain"]),
                              pkey=(p["chain"], p["uni"])))

    chosen, used_i, used_p = [], set(), set()
    for c in sorted(props, key=lambda x: (-x["anchor"], -x["q"])):
        if c["i"] in used_i:
            continue
        if c["pkey"] in used_p and len(chosen) >= 2:
            continue
        chosen.append(c); used_i.add(c["i"]); used_p.add(c["pkey"])
        if len(chosen) >= 4:
            break
    return chosen, epitope

# ---- liability repair (unchanged from Challenge 1, which unit-tested it) ----------
SUB = {"C":"S","L":"I","I":"L","V":"I","E":"Q","Q":"E","K":"R","R":"K","P":"A",
       "A":"S","S":"T","T":"S","G":"A","N":"Q","D":"N","M":"L","F":"Y","Y":"F",
       "W":"F","H":"H"}
RUN_MIN = 5

def repair(seq, protect=()):
    s, notes, protect = list(seq), [], set(protect)
    for i, c in enumerate(s):
        if c == "C" and i not in protect:
            s[i] = "S"; notes.append("C%dS" % (i+1))
    i = 0
    while i < len(s):
        j = i
        while j+1 < len(s) and s[j+1] == s[i]:
            j += 1
        if j-i+1 >= RUN_MIN:
            for k in range(i+2, j+1, 2):
                if k in protect:
                    continue
                alt = SUB.get(s[k], "A")
                if alt == s[k]:
                    alt = "A"
                notes.append("%s%d%s" % (s[k], k+1, alt)); s[k] = alt
        i = j+1
    return "".join(s), notes

def flags(s):
    runs, i = [], 0
    while i < len(s):
        j = i
        while j+1 < len(s) and s[j+1] == s[i]:
            j += 1
        if j-i+1 >= 4:
            runs.append("%s%d" % (s[i], j-i+1))
        i = j+1
    return s.count("C"), runs

# ---- scoring ---------------------------------------------------------------------
from colabdesign import mk_afdesign_model, clear_mem
HOTSPOT = "A86,A87,A88,A90,A91,A92,B29,B32,B33,B34,B35,B36"
_c = {}
def score(seq):
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
    return (round(float(lg.get("i_ptm", lg.get("ptm", 0.0))), 3),
            round(float(lg.get("plddt", 0.0)), 3))

rj = glob.glob("/kaggle/input/**/tnf_round2.json", recursive=True)
if not rj:
    raise SystemExit("tnf_round2.json not found - attach the tnf_generate output")
designs = json.load(open(rj[0]))
pdb_by = {}
for p in glob.glob("/kaggle/input/**/tnf_out2/*.pdb", recursive=True):
    pdb_by.setdefault(os.path.basename(p)[:-4], p)
kept = [d for d in designs
        if d["true_iptm"] >= MIN_IPTM and len(d["on_site"]) >= MIN_ONSITE
        and d["name"] in pdb_by]
# Sorting by LENGTH, not by score: every change of binder length forces a model
# rebuild and recompile, which is the dominant cost when only predicting.
kept.sort(key=lambda d: (d["length"], d["name"]))
mine = [d for i, d in enumerate(kept) if i % 2 == SHARD]
if SHARD == 0:
    log("%d designs, %d pass i_ptm>=%.2f and %d+ on-site contacts"
        % (len(designs), len(kept), MIN_IPTM, MIN_ONSITE))
log("this worker takes %d: %s" % (len(mine), ", ".join(d["name"] for d in mine)))

rec_path = "/kaggle/working/tnf_switch2_w%d.jsonl" % SHARD
for d in mine:
    if (time.time()-t0)/60 > BUDGET:
        log("time budget reached"); break
    name, seq, par = d["name"], d["sequence"], d["true_iptm"]
    try:
        chosen, epi = analyse(pdb_by[name], seq)
    except Exception as e:
        log("%s ANALYSIS FAILED %s %s" % (name, type(e).__name__, str(e)[:200]))
        continue
    base = dict(source_backbone=name, length=len(seq), parent_i_ptm=par,
                parent_plddt=d["plddt"], epitope_residues=epi,
                on_site=sorted(set(epi) & SITE),
                generator="AfDesign binder hallucination (ColabDesign)")

    rep, notes = repair(seq)
    if rep == seq:
        c_iptm, c_plddt = par, d["plddt"]
    else:
        c_iptm, c_plddt = score(rep)
    cys, runs = flags(rep)
    log("%s CONTROL repair %s -> i_ptm %.3f vs %.3f (%.0f%%) cys=%d runs=%s %s"
        % (name, ",".join(notes) or "none needed", c_iptm, par, 100*c_iptm/par,
           cys, runs or "none", "OK" if c_iptm >= REPAIR_KEEP*par else "REPAIR COSTLY"))
    out = [dict(name=name+"_ctl", sequence=rep, role="control", n_his=0, n_anchor=0,
                switch_score=0.0, his_substitutions="none", his_partners=[],
                repair=",".join(notes) or "none", i_ptm=c_iptm, plddt=c_plddt,
                iptm_retained=round(c_iptm/par, 3), free_cys=cys,
                homopolymer_runs=runs, **base)]

    if not chosen:
        log("   %s: no acid-OFF site found -> control only" % name)
    else:
        ss, subs, hidx = list(seq), [], []
        for c in chosen:
            if ss[c["i"]] != "H":
                subs.append("%s%dH" % (ss[c["i"]], c["i"]+1))
                ss[c["i"]] = "H"; hidx.append(c["i"])
        sw, sw_notes = repair("".join(ss), protect=hidx)
        s_iptm, s_plddt = score(sw)
        cys, runs = flags(sw)
        ret = round(s_iptm/par, 3) if par else 0.0
        nanc = sum(1 for c in chosen if c["anchor"])
        log("   %s SWITCH %s (%d anchor) vs %s + repair %s -> i_ptm %.3f (%.0f%% of "
            "parent) cys=%d runs=%s %s"
            % (name, ",".join(subs), nanc,
               "/".join(sorted({c["partner"] for c in chosen})),
               ",".join(sw_notes) or "none", s_iptm, 100*ret, cys, runs or "none",
               "KEEP" if ret >= KEEP_FRAC else "switch breaks interface"))
        out.append(dict(name=name+"_sw", sequence=sw, role="switch",
                        n_his=len(chosen), n_anchor=nanc,
                        switch_score=round(sum(c["q"] for c in chosen), 3),
                        his_substitutions=",".join(subs) or "none",
                        his_partners=sorted({c["partner"] for c in chosen}),
                        repair=",".join(sw_notes) or "none", i_ptm=s_iptm,
                        plddt=s_plddt, iptm_retained=ret, free_cys=cys,
                        homopolymer_runs=runs, **base))
    with open(rec_path, "a") as fh:
        for r in out:
            fh.write(json.dumps(r) + chr(10))
log("worker done")
"""

with open("/kaggle/working/tnf_switch2_worker.py", "w") as fh:
    fh.write(textwrap.dedent(WORKER))

args = [str(MIN_IPTM), str(MIN_ONSITE), str(KEEP_FRAC), str(REPAIR_KEEP)]
procs = []
for shard, gpu in enumerate(["0", "1"]):
    p = subprocess.Popen([sys.executable, "-u", "/kaggle/working/tnf_switch2_worker.py",
                          gpu, str(shard), str(TIME_BUDGET_MIN)] + args,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         text=True, bufsize=1)
    procs.append(p)
    print("launched switch worker %d on GPU %s" % (shard, gpu), flush=True)

import threading
def pump(p):
    for line in p.stdout:
        print(line.rstrip(), flush=True)
threads = [threading.Thread(target=pump, args=(p,), daemon=True) for p in procs]
[t.start() for t in threads]
for p in procs:
    p.wait()
[t.join(timeout=10) for t in threads]

records = []
for f in sorted(glob.glob("/kaggle/working/tnf_switch2_w*.jsonl")):
    for line in open(f):
        line = line.strip()
        if line:
            records.append(json.loads(line))
json.dump(records, open("/kaggle/working/tnf_designs_final2.json", "w"), indent=1)

sw = [r for r in records if r["role"] == "switch"]
good = [r for r in sw if r["iptm_retained"] >= KEEP_FRAC and r["i_ptm"] >= MIN_IPTM]
anch = [r for r in good if r["n_anchor"] >= 1]
clean = [r for r in records if r["free_cys"] == 0 and not r["homopolymer_runs"]]
print("\n==== TNF SWITCH STAGE ====")
print("%-18s %5s %4s %4s %7s %7s %7s %5s %-10s"
      % ("name", "role", "His", "anc", "i_ptm", "parent", "kept", "cys", "runs"))
for r in sorted(records, key=lambda x: (-(x["role"] == "switch"), -x["i_ptm"])):
    print("%-18s %5s %4d %4d %7.3f %7.3f %6.0f%% %5d %-10s"
          % (r["name"], r["role"][:5], r["n_his"], r["n_anchor"], r["i_ptm"],
             r["parent_i_ptm"], 100*r["iptm_retained"], r["free_cys"],
             ",".join(r["homopolymer_runs"]) or "-"))
print("\n%d records: %d switch (%d usable, %d engaging a conserved anchor), "
      "%d control; %d sequences fully clean"
      % (len(records), len(sw), len(good), len(anch), len(records)-len(sw), len(clean)))
print("\n===== COPY FROM HERE =====")
for r in sorted(records, key=lambda x: (-(x["role"] == "switch"), -x["i_ptm"])):
    print(json.dumps(r, separators=(",", ":")))
print("===== COPY TO HERE =====")
