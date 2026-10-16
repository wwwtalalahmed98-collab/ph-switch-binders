# === CHALLENGE 2 (TNF-alpha) : MOUSE-AWARE acid-OFF SWITCH PLACEMENT ===
#
# Why this run exists.
#
# The parent screen measured 31 parent binders against mouse TNF-alpha and found 20 of
# them bind it. Only 7 of the switches built from those parents do. Comparing the 17
# lineages that have both numbers gives the reason:
#
#     the same histidine substitutions keep a median 90% of HUMAN binding
#     and a median 49% of MOUSE binding
#
# That asymmetry is not a property of histidine. It is the selection criterion showing
# through. Every switch stage so far scored candidate positions on the human structure
# and kept or rejected them on human i_ptm retention alone, so a position that is
# tolerable in human and destructive in mouse passed the filter by construction.
# tnf2_w0_007 is the clearest case: parent mouse i_ptm 0.719 with 14 mouse-site
# contacts, switch 0.148 with none, while keeping 83% of its human binding.
#
# So this run changes one thing and only one thing: a histidine placement is kept only
# if it holds BOTH species. Everything else - the geometry, the scoring constants, the
# repair, the paired control, the honest re-prediction at 3 recycles - is carried over
# unchanged from cellTNF15 so that the comparison is clean.
#
# Two further changes follow from that:
#
#   1. Parents come from the measured mouse screen, not from human score alone. Only
#      lineages that already bind mouse are offered a histidine. Four of them
#      (tnf3_w0_004, tnf3_w0_005, tnf3_w1_009, tnf3_w1_019) have never been through a
#      switch stage at all because round 3's ran before any mouse number existed.
#   2. Each parent is offered SEVERAL histidine placements rather than one. Previously
#      the greedy best-scoring set was taken and that was the design. Here the top
#      single positions and the greedy pair are all built and measured, and the variant
#      that holds both species best is the one kept. The geometry proposes; the
#      measurement disposes.
#
# Cost control: the two targets are scored in two separate passes, each sorted by
# binder length, so only one AF model is resident at a time and a rebuild happens once
# per length per species rather than on every alternation.
import glob, json, os, subprocess, sys, textwrap, time

TIME_BUDGET_MIN = 230          # inside the GPU quota left this week
HUMAN_KEEP = 0.70              # a switch must hold this much of the parent's HUMAN i_ptm
MOUSE_KEEP = 0.70              # and this much of its MOUSE i_ptm - the new criterion
MIN_MOUSE_ONSITE = 5
REPAIR_KEEP = 0.85
PARENT_HUMAN_MIN = 0.45   # the screen cut at 0.50 and left 9 lineages unexamined
MOUSE_FLOOR = 0.45        # absolute, so retaining 70% of a bad mouse number cannot pass

WORKER = r"""
import glob, json, os, sys, time
from collections import Counter
GPU, SHARD, BUDGET = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
os.environ["CUDA_VISIBLE_DEVICES"] = GPU          # must precede the jax import
HUMAN_KEEP, MOUSE_KEEP = float(sys.argv[4]), float(sys.argv[5])
MIN_MOUSE_ONSITE, REPAIR_KEEP = int(sys.argv[6]), float(sys.argv[7])
PARENT_HUMAN_MIN, MOUSE_FLOOR = float(sys.argv[8]), float(sys.argv[9])

t0 = time.time()
def log(m): print("[w%d] %s [%.1f min]" % (SHARD, m, (time.time()-t0)/60), flush=True)

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import protein_letters_3to1

# ---- the two targets -------------------------------------------------------------
HT = "/kaggle/working/targets/TNF_site_crop.pdb"
MT = "/kaggle/working/targets/TNF_mouse_site_crop.pdb"
CMAP = json.load(open("/kaggle/working/construct_map.json"))["order"]
N = len(CMAP)
MCM = json.load(open("/kaggle/working/construct_map_mouse.json"))
MMAP = MCM["order"]
MN = len(MMAP)
MSITE = set(MCM["mouse_site"])
M2H = {int(k): v for k, v in MCM["mouse_to_human_site"].items()}
H_HOTSPOT = "A86,A87,A88,A90,A91,A92,B29,B32,B33,B34,B35,B36"
M_HOTSPOT = open("/kaggle/working/mouse_hotspot.txt").read().strip()

SITE = {95,105,106,107,108,109,110,111,112,155,159,161,162,163,164,165,166,167,168,
        201,203,207,220,221,223}
ANCHORS = {("A",166), ("B",108)}      # conserved, exposed, from PHASE1_REPORT

# Constants unchanged from cellTNF15, calibrated on 17 structures (278 His-carboxylate
# and 279 His-cation closest approaches) and validated against the FcRn/IgG acid-ON
# switch. Nothing here is retuned - the point of this run is to change the ACCEPTANCE
# test, not the proposal geometry, so that any difference is attributable.
MAX_HIS = 2
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
N_VARIANTS = 4          # up to 3 single-histidine options plus the greedy pair

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

# Free-state exposure of every human target residue, measured on the cropped target
# ITSELF rather than on a chain pulled out of the complex: the complex file duplicates
# its target residues and SASA over the duplicate reports everything as buried.
_free = {}
for _ch in PDBParser(QUIET=True).get_structure("t", HT)[0]:
    _one = PDBParser(QUIET=True).get_structure("x", HT)[0][_ch.id]
    for _k, _v in rsasa_chain(_one).items():
        _free[(_ch.id, _k)] = _v
TFREE = [_free[(o["chain"], o["pdb"])] for o in CMAP]

CATIONS = []
for i, o in enumerate(CMAP):
    if o["aa"] in ("R","K") and TFREE[i] >= MIN_PARTNER_EXPOSURE:
        CATIONS.append(dict(idx=i, chain=o["chain"], uni=o["uniprot"], aa=o["aa"],
                            anchor=(o["chain"], o["uniprot"]) in ANCHORS))
if SHARD == 0:
    log("human construct %d residues; %d exposed cations (%d conserved anchors)"
        % (N, len(CATIONS), sum(c["anchor"] for c in CATIONS)))
got = {(c["chain"], c["uni"]) for c in CATIONS if c["anchor"]}
assert got == ANCHORS, "a conserved anchor is missing or buried: %s" % (ANCHORS - got)

# Cut a chain back to n residues when save_pdb has written it more than once. A chain
# overlapping itself reads as fully buried, which would silently find zero histidine
# sites and produce a clean-looking "no switch possible" result.
def trim(chain, n, what):
    res = [r for r in chain if r.id[0] == " "]
    if len(res) > n and len(res) % n == 0:
        for r in res[n:]:
            chain.detach_child(r.id)
        res = res[:n]
    if len(res) != n:
        raise ValueError("%s chain has %d residues, expected %d" % (what, len(res), n))
    return res

def load(pdb_path, blen, cmap):
    n = len(cmap)
    mdl = PDBParser(QUIET=True).get_structure("c", pdb_path)[0]
    sizes = Counter({ch.id: sum(1 for r in ch if r.id[0]==" ") for ch in mdl})
    tgt, bnd = max(sizes, key=lambda k: sizes[k]), min(sizes, key=lambda k: sizes[k])
    tres = trim(mdl[tgt], n, "target")
    bres = trim(mdl[bnd], blen, "binder")
    seen = "".join(protein_letters_3to1.get(r.get_resname(), "X") for r in tres)
    if seen != "".join(o["aa"] for o in cmap):
        raise ValueError("target sequence does not match its construct map")
    return mdl, tgt, bnd, tres, bres

def contacts(pdb_path, blen, cmap):
    _m, _t, _b, tres, bres = load(pdb_path, blen, cmap)
    batoms = np.array([a.get_coord() for r in bres for a in r])
    hits = []
    for i, r in enumerate(tres):
        ra = np.array([a.get_coord() for a in r])
        if np.min(np.linalg.norm(ra[:,None,:]-batoms[None,:,:], axis=-1)) <= 5.0:
            hits.append(cmap[i]["uniprot"])
    return sorted(set(hits))

def propose(pdb_path, seq):
    mdl, tgt, bnd, tres, bres = load(pdb_path, len(seq), CMAP)
    _m2, _t2, _b2, _, _br2 = load(pdb_path, len(seq), CMAP)
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
    return props, epitope

def variants(props):
    # The greedy pair is what every previous run would have built; it is kept as one
    # option among several rather than as the answer.
    ranked = sorted(props, key=lambda x: (-x["anchor"], -x["q"]))
    pair, used_i, used_p = [], set(), set()
    for c in ranked:
        if c["i"] in used_i:
            continue
        if c["pkey"] in used_p and len(pair) >= 2:
            continue
        pair.append(c); used_i.add(c["i"]); used_p.add(c["pkey"])
        if len(pair) >= MAX_HIS:
            break
    out = []
    if len(pair) >= 2:
        out.append(("pair", pair))
    seen_i = set()
    for c in ranked:
        if c["i"] in seen_i:
            continue
        seen_i.add(c["i"])
        out.append(("single", [c]))
        if len(out) >= N_VARIANTS:
            break
    return out[:N_VARIANTS]

# ---- liability repair, in two stages that are paid for separately ----------------
#
# The previous version charged both liabilities at once and never asked what they cost.
# On tnf3_w1_019 the homopolymer rule turned four glutamates into glutamines, the human
# i_ptm went 0.714 -> 0.472 and the repaired control's mouse i_ptm fell to 28% of the
# parent, so the switch stage never got a fair attempt at a design that binds both
# species. Three other lineages lost 17-35% the same way.
#
# Two changes. First, the substitution table no longer throws away charge: E->Q removed
# a negative charge at every one of those four positions. E and D now swap for each
# other, K and R for each other, Q and N for each other.
#
# Second, the two liabilities are separated. A free cysteine is a real liability - it
# scrambles disulfides - and is always removed. Breaking a homopolymer run is a
# synthesis convenience, so it is applied, measured, and REVERTED if it costs more than
# REPAIR_KEEP of the parent's i_ptm. The run is then recorded as an unrepaired liability
# rather than silently paid for.
SUB = {"C":"S",
       "E":"D","D":"E","K":"R","R":"K","Q":"N","N":"Q",
       "L":"I","I":"L","V":"I","M":"L","F":"Y","Y":"F","W":"F",
       "A":"S","S":"T","T":"S","G":"A","P":"A","H":"H"}
RUN_MIN = 5

def _cys_plan(seq):
    return [(i, "S") for i, c in enumerate(seq) if c == "C"]

def _run_plan(seq):
    s, plan, i = list(seq), [], 0
    while i < len(s):
        j = i
        while j+1 < len(s) and s[j+1] == s[i]:
            j += 1
        if j-i+1 >= RUN_MIN:
            for k in range(i+2, j+1, 2):
                alt = SUB.get(s[k], "A")
                if alt == s[k]:
                    alt = "A"
                plan.append((k, alt))
        i = j+1
    return plan

# A plan is a list of (index, new residue). Applying the SAME plan to the control and to
# each switch - skipping only the positions that carry a histidine - is what guarantees
# the pair differs at the histidines and nowhere else. The old code re-derived the repair
# from each sequence separately, so introducing a histidine could change which run was
# broken and quietly desynchronise the twins.
def apply_plan(seq, plan, skip=()):
    s, skip = list(seq), set(skip)
    for i, a in plan:
        if i not in skip:
            s[i] = a
    return "".join(s)

def plan_notes(seq, plan, skip=()):
    skip = set(skip)
    return ["%s%d%s" % (seq[i], i+1, a) for i, a in plan if i not in skip and seq[i] != a]

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

# ---- scoring: one model resident at a time ---------------------------------------
from colabdesign import mk_afdesign_model, clear_mem
_c = {}
def score(seq, species, out_pdb=None):
    key = (species, len(seq))
    if _c.get("key") != key:
        clear_mem()
        af = mk_afdesign_model(protocol="binder", data_dir="/kaggle/working")
        af.prep_inputs(pdb_filename=HT if species == "human" else MT, chain="A,B",
                       binder_len=len(seq),
                       hotspot=H_HOTSPOT if species == "human" else M_HOTSPOT)
        _c["key"], _c["af"] = key, af
    af = _c["af"]
    try:
        af.predict(seq=seq, num_recycles=3, verbose=False)
    except TypeError:
        af.predict(seq=seq, verbose=False)
    lg = af.aux["log"]
    if out_pdb:
        af.save_pdb(out_pdb)
    return (round(float(lg.get("i_ptm", lg.get("ptm", 0.0))), 3),
            round(float(lg.get("plddt", 0.0)), 3))

# ---- the parents: measured mouse binders only ------------------------------------
mp = glob.glob("/kaggle/input/**/tnf_mouse_parents.json", recursive=True)
if not mp:
    raise SystemExit("tnf_mouse_parents.json not found - attach the tnf_mouse_parents2 "
                     "output; this run selects parents on MEASURED mouse binding")
screen = {r["name"]: r for r in json.load(open(mp[0]))}
seqs = {}
for f in sorted(glob.glob("/kaggle/input/**/tnf_round*.json", recursive=True)):
    for r in json.load(open(f)):
        seqs.setdefault(r["name"], r)
pdb_by = {}
for p in glob.glob("/kaggle/input/**/tnf_out*/*.pdb", recursive=True):
    pdb_by.setdefault(os.path.basename(p)[:-4], p)

kept = [dict(screen[n], parent_mouse_known=True) for n in screen
        if screen[n]["mouse_i_ptm"] >= 0.45
        and len(screen[n]["mouse_on_site"]) >= MIN_MOUSE_ONSITE
        and n in seqs and n in pdb_by]

# The screen filtered parents at human i_ptm >= 0.50, so every lineage between 0.45 and
# 0.50 was never measured against mouse at all and never offered a histidine here. One
# of them, tnf2_w0_001 at human 0.481, is in the current submission and its SWITCH
# measures mouse 0.654 with 12 mouse-site contacts - better than most of the set. These
# have no mouse parent number, so the worker measures one inline before placing.
have = {r["name"] for r in kept}
band = [dict(name=r["name"], length=r["length"], round="?",
             human_i_ptm=r["true_iptm"], human_on_site=len(r.get("on_site", [])),
             mouse_i_ptm=None, mouse_on_site=[], parent_mouse_known=False)
        for r in seqs.values()
        if r["name"] not in have and r["name"] in pdb_by
        and PARENT_HUMAN_MIN <= r.get("true_iptm", 0) < 0.50
        and len(r.get("on_site", [])) >= MIN_MOUSE_ONSITE]
kept += band
# Sorting by LENGTH, not by score: every change of binder length forces a model
# rebuild and recompile, which is the dominant cost when only predicting.
kept.sort(key=lambda r: (r["length"], r["name"]))
mine = [r for i, r in enumerate(kept) if i % 2 == SHARD]
if SHARD == 0:
    log("%d measured mouse binders + %d unscreened parents at human %.2f-0.50 = %d"
        % (len(kept)-len(band), len(band), PARENT_HUMAN_MIN, len(kept)))
    log("unscreened band: %s" % ", ".join(sorted(r["name"] for r in band)))
missing = [n for n in screen if screen[n]["mouse_i_ptm"] >= 0.45
           and len(screen[n]["mouse_on_site"]) >= MIN_MOUSE_ONSITE
           and (n not in seqs or n not in pdb_by)]
if SHARD == 0 and missing:
    log("WARNING mouse-binding parents with no backbone/sequence attached: %s" % missing)
log("this worker takes %d: %s" % (len(mine), ", ".join(r["name"] for r in mine)))

# ---- PASS 1: propose histidines, score every variant against HUMAN ---------------
work = []
for d in mine:
    if (time.time()-t0)/60 > BUDGET * 0.55:
        log("human pass budget reached"); break
    name = d["name"]
    seq, par_h = seqs[name]["sequence"], seqs[name]["true_iptm"]
    par_m = d["mouse_i_ptm"]
    try:
        props, epi = propose(pdb_by[name], seq)
    except Exception as e:
        log("%s PROPOSE FAILED %s %s" % (name, type(e).__name__, str(e)[:200]))
        continue
    vs = variants(props)
    if not vs:
        log("%s: no acid-OFF site passes the geometry -> control only" % name)

    # Decide the repair LEVEL once, on the control, and apply the identical plan to
    # every switch built from this parent. Deciding per sequence is what let the twins
    # drift apart, and paying for the run-breaking unconditionally is what cost
    # tnf3_w1_019 a third of its binding before any histidine was placed.
    cys_plan = _cys_plan(seq)
    cys_seq = apply_plan(seq, cys_plan)
    c_cys = (par_h, seqs[name]["plddt"]) if cys_seq == seq else score(cys_seq, "human")
    run_plan = _run_plan(cys_seq)
    if not run_plan:
        level, plan, (c_h, c_hp) = "cys-only", cys_plan, c_cys
        why = "no run to break"
    else:
        full_seq = apply_plan(cys_seq, run_plan)
        c_full = score(full_seq, "human")
        if c_full[0] >= REPAIR_KEEP * par_h or c_full[0] >= c_cys[0]:
            level, plan, (c_h, c_hp) = "full", cys_plan + run_plan, c_full
            why = "run-breaking held (%.3f)" % c_full[0]
        else:
            level, plan, (c_h, c_hp) = "cys-only", cys_plan, c_cys
            why = ("run-breaking REVERTED: %.3f vs %.3f cys-only, below %.0f%% of parent"
                   % (c_full[0], c_cys[0], 100*REPAIR_KEEP))
    ctl = apply_plan(seq, plan)
    ctl_notes = plan_notes(seq, plan)
    cys, runs = flags(ctl)
    log("%s CONTROL repair[%s] %s -> human %.3f vs %.3f (%.0f%%) | %s"
        % (name, level, ",".join(ctl_notes) or "none needed", c_h, par_h,
           100*c_h/par_h, why))

    rec = dict(name=name, length=len(seq), parent_human=par_h, parent_mouse=par_m,
               parent_mouse_known=d["parent_mouse_known"],
               parent_mouse_on_site=len(d["mouse_on_site"]),
               repair_level=level, repair_note=why,
               human_epitope=epi, human_on_site=sorted(set(epi) & SITE),
               control=dict(sequence=ctl, repair=",".join(ctl_notes) or "none",
                            human_i_ptm=c_h, human_plddt=c_hp,
                            free_cys=cys, homopolymer_runs=runs),
               variants=[])
    for kind, chosen in vs:
        ss, subs, hidx = list(seq), [], []
        for c in chosen:
            if ss[c["i"]] != "H":
                subs.append("%s%dH" % (ss[c["i"]], c["i"]+1))
                ss[c["i"]] = "H"; hidx.append(c["i"])
        if not subs:
            continue
        sw = apply_plan("".join(ss), plan, skip=hidx)
        sw_notes = plan_notes("".join(ss), plan, skip=hidx)
        # the twins must differ at the histidines and nowhere else
        diff = [k for k in range(len(sw)) if sw[k] != ctl[k]]
        if sorted(diff) != sorted(hidx):
            log("   %s %s PAIRING BROKEN: differs at %s, expected %s"
                % (name, ",".join(subs), diff, sorted(hidx)))
            continue
        s_h, s_hp = score(sw, "human")
        hret = round(s_h/par_h, 3) if par_h else 0.0
        cys, runs = flags(sw)
        log("   %s %-6s %s -> human %.3f (%3.0f%% of parent) %s"
            % (name, kind, ",".join(subs), s_h, 100*hret,
               "ok" if hret >= HUMAN_KEEP else "breaks human interface"))
        rec["variants"].append(dict(
            kind=kind, sequence=sw, n_his=len(chosen),
            n_anchor=sum(1 for c in chosen if c["anchor"]),
            switch_score=round(sum(c["q"] for c in chosen), 3),
            his_substitutions=",".join(subs),
            his_partners=sorted({c["partner"] for c in chosen}),
            his_geometry=[[c["i"]+1, c["partner"], c["d"]] for c in chosen],
            repair=",".join(sw_notes) or "none",
            human_i_ptm=s_h, human_plddt=s_hp, human_retained=hret,
            free_cys=cys, homopolymer_runs=runs))
    work.append(rec)

# ---- PASS 2: score the human survivors, and every control, against MOUSE ---------
# Every variant that held human is measured against mouse. Nothing is pre-filtered on
# a human-derived guess about which will survive - that guess is exactly what this run
# exists to stop making.
os.makedirs("/kaggle/working/switch5_mouse_pdb", exist_ok=True)
work.sort(key=lambda r: (r["length"], r["name"]))

def emit(rec):
    with open("/kaggle/working/tnf_switch5_w%d.jsonl" % SHARD, "a") as fh:
        fh.write(json.dumps(rec) + chr(10))

for pos, rec in enumerate(work):
    if (time.time()-t0)/60 > BUDGET:
        # Write the rest out with their human numbers rather than dropping them. An
        # earlier version broke here and silently discarded every record the mouse pass
        # had not reached, including work that had already been paid for.
        for leftover in work[pos:]:
            leftover["best"] = None
            leftover["incomplete"] = "mouse pass not reached before the time budget"
            emit(leftover)
        log("mouse pass budget reached; %d records written with human data only"
            % len(work[pos:]))
        break
    name = rec["name"]
    # Parents from the unscreened 0.45-0.50 band have no mouse number yet. Measure the
    # PARENT, not the control, so the denominator means the same thing it does for the
    # screened lineages. The model is already on this species and length, so it costs
    # one prediction and no rebuild.
    if not rec["parent_mouse_known"]:
        try:
            out = "/kaggle/working/switch5_mouse_pdb/%s_parent.pdb" % name
            pm, _pp = score(seqs[name]["sequence"], "mouse", out)
            pepi = contacts(out, rec["length"], MMAP)
            pon = sorted(set(pepi) & MSITE)
            rec["parent_mouse"] = pm
            rec["parent_mouse_on_site"] = len(pon)
            log("%s PARENT measured against mouse: %.3f, %d on mouse site%s"
                % (name, pm, len(pon),
                   "" if (pm >= 0.45 and len(pon) >= MIN_MOUSE_ONSITE)
                   else "  (does not bind mouse - variants will be recorded, not kept)"))
        except Exception as e:
            log("%s parent mouse FAILED %s %s" % (name, type(e).__name__, str(e)[:200]))
    par_m = rec["parent_mouse"]
    try:
        out = "/kaggle/working/switch5_mouse_pdb/%s_ctl.pdb" % name
        m, mp_ = score(rec["control"]["sequence"], "mouse", out)
        epi = contacts(out, rec["length"], MMAP)
        on = sorted(set(epi) & MSITE)
        rec["control"].update(mouse_i_ptm=m, mouse_plddt=mp_, mouse_on_site=on,
                              mouse_retained=round(m/par_m, 3) if par_m else None)
        log("%s CONTROL mouse %.3f vs parent %.3f (%3.0f%%), %d on mouse site"
            % (name, m, par_m, 100*m/par_m if par_m else 0, len(on)))
    except Exception as e:
        log("%s control mouse FAILED %s %s" % (name, type(e).__name__, str(e)[:200]))
    for v in rec["variants"]:
        if v["human_retained"] < HUMAN_KEEP:
            v["mouse_i_ptm"] = None
            v["skipped"] = "human retention below %.2f" % HUMAN_KEEP
            continue
        try:
            out = ("/kaggle/working/switch5_mouse_pdb/%s_%s.pdb"
                   % (name, v["his_substitutions"].replace(",", "_")))
            m, mp_ = score(v["sequence"], "mouse", out)
            epi = contacts(out, rec["length"], MMAP)
            on = sorted(set(epi) & MSITE)
            mret = round(m/par_m, 3) if par_m else 0.0
            v.update(mouse_i_ptm=m, mouse_plddt=mp_, mouse_on_site=on,
                     on_site_human_equiv=sorted(M2H[x] for x in on if x in M2H),
                     mouse_retained=mret,
                     passes=bool(v["human_retained"] >= HUMAN_KEEP
                                 and mret >= MOUSE_KEEP
                                 and m >= MOUSE_FLOOR
                                 and len(on) >= MIN_MOUSE_ONSITE))
            log("   %s %-6s %s -> mouse %.3f (%3.0f%% of parent), %2d on mouse site  %s"
                % (name, v["kind"], v["his_substitutions"], m, 100*mret, len(on),
                   "KEEP" if v["passes"] else "loses mouse"))
        except Exception as e:
            log("   %s %s mouse FAILED %s %s"
                % (name, v["his_substitutions"], type(e).__name__, str(e)[:200]))

    # The variant kept is the one whose WORSE species is best. Maximising the sum would
    # let a design that is excellent in human and broken in mouse win, which is the
    # failure this whole run is about.
    ok = [v for v in rec["variants"] if v.get("passes")]
    rec["best"] = (max(ok, key=lambda v: min(v["human_retained"], v["mouse_retained"]))
                   ["his_substitutions"] if ok else None)
    emit(rec)
log("worker done")
"""

with open("/kaggle/working/tnf_switch_worker5.py", "w") as fh:
    fh.write(textwrap.dedent(WORKER))

# The parent screen ran its workers before the setup cell had installed biopython and
# reported "Successfully ran" over an empty table. Refuse to launch without the files.
subprocess.run("pip install -q biopython", shell=True)
NEEDED = ["/kaggle/working/targets/TNF_site_crop.pdb",
          "/kaggle/working/targets/TNF_mouse_site_crop.pdb",
          "/kaggle/working/construct_map.json",
          "/kaggle/working/construct_map_mouse.json",
          "/kaggle/working/mouse_hotspot.txt"]
missing = [p for p in NEEDED if not os.path.exists(p)]
if missing:
    raise SystemExit("a setup cell has not run - missing %s" % missing)

args = [str(HUMAN_KEEP), str(MOUSE_KEEP), str(MIN_MOUSE_ONSITE), str(REPAIR_KEEP),
        str(PARENT_HUMAN_MIN), str(MOUSE_FLOOR)]
procs = []
for shard, gpu in enumerate(["0", "1"]):
    p = subprocess.Popen([sys.executable, "-u", "/kaggle/working/tnf_switch_worker5.py",
                          gpu, str(shard), str(TIME_BUDGET_MIN)] + args,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         text=True, bufsize=1)
    procs.append(p)
    print("launched mouse-aware switch worker %d on GPU %s" % (shard, gpu), flush=True)

import threading
def pump(p):
    for line in p.stdout:
        print(line.rstrip(), flush=True)
threads = [threading.Thread(target=pump, args=(p,), daemon=True) for p in procs]
[t.start() for t in threads]
for p in procs:
    p.wait()
[t.join(timeout=10) for t in threads]

recs = []
for f in sorted(glob.glob("/kaggle/working/tnf_switch5_w*.jsonl")):
    for line in open(f):
        line = line.strip()
        if line:
            recs.append(json.loads(line))
json.dump(recs, open("/kaggle/working/tnf_switch5.json", "w"), indent=1)
if not recs:
    raise SystemExit("no parent produced a record - read the worker log above rather "
                     "than treating this as 'no switch is possible'")

print("\n==== MOUSE-AWARE SWITCH PLACEMENT ====")
print("%-14s %-12s %3s %7s %7s %7s %7s %5s %s"
      % ("parent", "his", "anc", "humanP", "human", "mouseP", "mouse", "site", "verdict"))
n_pass = 0
for r in sorted(recs, key=lambda x: -x["parent_mouse"]):
    for v in r["variants"]:
        m = v.get("mouse_i_ptm")
        print("%-14s %-12s %3d %7.3f %7.3f %7.3f %7s %5s %s"
              % (r["name"], v["his_substitutions"], v["n_anchor"],
                 r["parent_human"], v["human_i_ptm"], r["parent_mouse"],
                 ("%.3f" % m) if m is not None else "-",
                 str(len(v.get("mouse_on_site", []))) if m is not None else "-",
                 ("KEEP" if v.get("passes") else
                  v.get("skipped", "loses mouse")) +
                 (" <- best" if r.get("best") == v["his_substitutions"] else "")))
        n_pass += bool(v.get("passes"))
winners = [r for r in recs if r.get("best")]
print("\n%d parents offered histidines; %d variants hold BOTH species; "
      "%d parents yield at least one usable switch"
      % (len(recs), n_pass, len(winners)))
print("\nbest variant per parent (human retained / mouse retained):")
for r in sorted(winners, key=lambda x: -x["parent_mouse"]):
    v = next(v for v in r["variants"] if v["his_substitutions"] == r["best"])
    print("  %-14s %-12s human %.3f (%3.0f%%)  mouse %.3f (%3.0f%%)  %2d mouse-site"
          % (r["name"], v["his_substitutions"], v["human_i_ptm"],
             100*v["human_retained"], v["mouse_i_ptm"], 100*v["mouse_retained"],
             len(v["mouse_on_site"])))
print("\n===== COPY FROM HERE =====")
for r in recs:
    print(json.dumps(r, separators=(",", ":")))
print("===== COPY TO HERE =====")
