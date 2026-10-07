# === CHALLENGE 1 : ACID-ON SWITCH PLACEMENT + VARIANT RE-SCORING ===
#
# Round 1 emitted His variants without ever checking them. That is the mistake this
# cell exists to avoid: substituting 2-4 interface residues to histidine can destroy
# the interface it is meant to make pH-dependent, so every variant is re-predicted
# and the cost of the switch is recorded next to it.
#
# Geometry is the calibration from results/PHASE3_METHOD.md: measured on 23 genuine
# His-carboxylate salt bridges, validated by recovering barnase His102-Asp39 and by
# finding nothing in the cetuximab/EGFR interface.

import glob, json, os, time
from collections import Counter
import numpy as np

MIN_IPTM = 0.45          # a backbone must bind before a pH switch on it means anything
KEEP_FRAC = 0.70         # a variant must retain this much of the parent i_ptm
TIME_BUDGET_MIN = 300
t0 = time.time()

CB_MIN, CB_MAX, CB_IDEAL, COS_MIN = 3.5, 7.0, 5.8, -0.30
CARBOX = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}
# Every conserved, solvent-exposed carboxylate in the construct, not just the four on
# the Phase 1 patch. The selection rule was always "conserved + exposed + at the
# interface"; restricting it to the ligand face was an assumption about where the
# binders would land, and the epitope report showed they land elsewhere. Derived in
# scripts/phase1_conservation.py: Asp/Glu, identical in mouse, rel_SASA >= 0.20.
ACIDS = [21, 22, 51, 73, 90, 110, 155]   # UniProt 45,46,75,97,114,134,179
# human/mouse identity per construct residue, PDB numbering (1 = identical)
CONSERVED = json.loads(r'''{"4":1,"5":1,"6":1,"7":1,"8":1,"9":1,"10":1,"11":1,"12":1,"13":0,"14":1,"15":1,"16":1,"17":1,"18":1,"19":1,"20":1,"21":1,"22":1,"23":1,"24":1,"25":1,"26":1,"27":1,"28":1,"29":1,"30":1,"31":0,"32":1,"33":1,"34":1,"35":1,"36":1,"37":1,"38":1,"39":1,"40":1,"41":1,"42":1,"43":1,"44":1,"45":1,"46":1,"47":1,"48":1,"49":1,"50":1,"51":1,"52":1,"53":1,"54":1,"55":1,"56":1,"57":1,"58":1,"59":1,"60":1,"61":1,"62":1,"63":1,"64":1,"65":1,"66":1,"67":1,"68":1,"69":1,"70":1,"71":1,"72":1,"73":1,"74":1,"75":1,"76":1,"77":1,"78":1,"79":1,"80":1,"81":1,"82":1,"83":1,"84":1,"85":1,"86":1,"87":0,"88":0,"89":1,"90":1,"91":1,"92":0,"93":1,"94":1,"95":1,"96":1,"97":0,"98":1,"99":1,"100":1,"101":1,"102":0,"103":0,"104":1,"105":0,"106":1,"107":1,"108":1,"109":0,"110":1,"111":1,"112":1,"113":1,"114":1,"115":1,"116":1,"117":1,"118":1,"119":1,"120":1,"121":0,"122":1,"123":1,"124":1,"125":1,"126":1,"127":1,"128":1,"129":1,"130":1,"131":0,"132":1,"133":1,"134":1,"135":0,"136":0,"137":0,"138":1,"139":1,"140":1,"141":1,"142":1,"143":1,"144":1,"145":0,"146":0,"147":0,"148":1,"149":0,"150":1,"151":1,"152":1,"153":1,"154":1,"155":1,"156":0,"157":1,"158":0,"159":1,"160":0,"161":0,"162":1,"163":1,"164":0}''')
NOMUT = {"PRO", "GLY", "CYS"}
PREF = {"K": 1.5, "R": 1.5, "Q": 1.2, "N": 1.2, "S": 1.15, "T": 1.15, "A": 1.0,
        "D": 1.1, "E": 1.1, "H": 0.0, "M": 0.7, "L": 0.6, "I": 0.6, "V": 0.6,
        "F": 0.4, "Y": 0.4, "W": 0.3}
MAXS = {"A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "E": 223, "Q": 225,
        "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240,
        "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174}

# ---- inputs: round 2's designs and their predicted complexes ----
rj = glob.glob("/kaggle/input/**/round2.json", recursive=True)
if not rj:
    raise SystemExit("round2.json not found - attach the round2 notebook's output.")
designs = json.load(open(rj[0]))
pdb_by = {}
for p in glob.glob("/kaggle/input/**/out2/*.pdb", recursive=True):
    pdb_by.setdefault(os.path.basename(p).replace(".pdb", ""), p)
kept = [d for d in designs if d["true_iptm"] >= MIN_IPTM and d["name"] in pdb_by]
kept.sort(key=lambda d: -d["true_iptm"])
print("%d designs, %d with i_ptm >= %.2f and a structure"
      % (len(designs), len(kept), MIN_IPTM), flush=True)

from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import protein_letters_3to1


def vcb(r):
    if "CB" in r:
        return r["CB"].get_coord()
    try:
        n, a, c = r["N"].get_coord(), r["CA"].get_coord(), r["C"].get_coord()
    except KeyError:
        return None
    b, cc = a - n, c - a
    return -0.58273431 * np.cross(b, cc) + 0.56802827 * b - 0.54067466 * cc + a


def rsasa(entity):
    ShrakeRupley().compute(entity, level="R")
    out = {}
    for r in entity.get_residues():
        aa = protein_letters_3to1.get(r.get_resname())
        if r.id[0] == " " and aa:
            out[r.id[1]] = r.sasa / MAXS.get(aa, 200)
    return out


# ---- map our carboxylate numbering onto whatever the saved complex uses ----
# The target construct is chain A of 6ARU residues 4-164. AfDesign renumbers each
# chain of the complex it writes, so residue 21 of the saved file is NOT E45 - it is
# three residues off. Round 1 and the first run of this stage both hard-coded the
# source numbers against the saved file, found no ASP/GLU there, and reported
# "0 carboxylates" on designs that bind perfectly well. Map by position in the chain
# instead of trusting the numbers, and assert the result really is ASP/GLU.
SRC_TARGET = "/kaggle/working/targets/EGFR_A_domainI.pdb"
_src_chain = list(PDBParser(QUIET=True).get_structure("s", SRC_TARGET)[0])[0]
_src_res = [r for r in _src_chain if r.id[0] == " "]
_src_num = [r.id[1] for r in _src_res]
ACID_IDX = [_src_num.index(n) for n in ACIDS]
print("source target: %d residues, %d-%d; carboxylates %s at chain positions %s (%s)"
      % (len(_src_res), _src_num[0], _src_num[-1], ACIDS, ACID_IDX,
         [_src_res[i].get_resname() for i in ACID_IDX]), flush=True)
assert all(_src_res[i].get_resname() in CARBOX for i in ACID_IDX)
UNIPROT_OFFSET = 24          # PDB numbering + 24 = UniProt P00533 numbering


def place_switch(pdb_path, seq):
    """Pick His positions facing the conserved carboxylates. Returns
    (chosen, n distinct carboxylates, score, n candidate contacts)."""
    pr = PDBParser(QUIET=True)
    mdl = pr.get_structure("c", pdb_path)[0]
    sizes = Counter({ch.id: sum(1 for r in ch if r.id[0] == " ") for ch in mdl})
    tgt_ch = max(sizes, key=lambda k: sizes[k])
    bnd_ch = min(sizes, key=lambda k: sizes[k])
    bfree = rsasa(pr.get_structure("b", pdb_path)[0][bnd_ch])
    tfree = rsasa(pr.get_structure("t", pdb_path)[0][tgt_ch])
    ShrakeRupley().compute(mdl, level="R")
    tcplx = {r.id[1]: r.sasa / MAXS.get(protein_letters_3to1.get(r.get_resname(), "A"), 200)
             for r in mdl[tgt_ch] if r.id[0] == " "
             and r.get_resname() in protein_letters_3to1}

    tgt_res = [r for r in mdl[tgt_ch] if r.id[0] == " "]
    if len(tgt_res) != len(_src_res):
        raise ValueError("target chain is %d residues, source is %d - cannot map"
                         % (len(tgt_res), len(_src_res)))
    acids = []
    for src_n, i in zip(ACIDS, ACID_IDX):
        r = tgt_res[i]
        if r.get_resname() not in CARBOX:
            raise ValueError("position %d maps to %s, not a carboxylate" % (i, r.get_resname()))
        ox = [r[a].get_coord() for a in CARBOX[r.get_resname()] if a in r]
        if ox:
            acids.append(dict(num=src_n, ox=ox,
                              des=max(tfree.get(r.id[1], 0) - tcplx.get(r.id[1], 0), 0.0)))

    # Which epitope did the binder actually land on? A high i_ptm says the interface
    # is confident, not that it is the interface we asked for, and the whole
    # cross-reactivity and pH-switch rationale rests on it being the domain I ligand
    # face. Reported in UniProt numbering so it can be checked against Phase 1.
    bnd_atoms = np.array([a.get_coord() for r in mdl[bnd_ch] if r.id[0] == " "
                          for a in r])
    epitope = []
    for i, r in enumerate(tgt_res):
        ra = np.array([a.get_coord() for a in r])
        if np.min(np.linalg.norm(ra[:, None, :] - bnd_atoms[None, :, :], axis=-1)) <= 5.0:
            epitope.append(_src_num[i] + UNIPROT_OFFSET)

    props = []
    for r in mdl[bnd_ch]:
        if r.id[0] != " " or "CA" not in r or r.get_resname() in NOMUT:
            continue
        if bfree.get(r.id[1], 0) < 0.15:
            continue
        cb = vcb(r)
        if cb is None:
            continue
        d = cb - r["CA"].get_coord()
        nn = np.linalg.norm(d)
        if nn < 1e-6:
            continue
        d = d / nn
        for ac in acids:
            ds = [np.linalg.norm(cb - o) for o in ac["ox"]]
            dcb = float(min(ds))
            if not (CB_MIN <= dcb <= CB_MAX):
                continue
            to = ac["ox"][int(np.argmin(ds))] - cb
            to = to / np.linalg.norm(to)
            if float(np.dot(d, to)) < COS_MIN:
                continue
            q = max(0.0, 1 - abs(dcb - CB_IDEAL) / max(CB_MAX - CB_IDEAL, CB_IDEAL - CB_MIN))
            props.append(dict(pos=r.id[1], tgt=ac["num"], q=q * (1 + ac["des"])))

    def rank(x):
        j = x["pos"] - 1
        return -(x["q"] * (PREF.get(seq[j], 0.9) if 0 <= j < len(seq) else 0))

    chosen, used_pos, used_acid = [], set(), set()
    for cand in sorted(props, key=rank):
        j = cand["pos"] - 1
        if cand["pos"] in used_pos:
            continue
        if not (0 <= j < len(seq)) or PREF.get(seq[j], 0.9) <= 0:
            continue
        if cand["tgt"] in used_acid and len(chosen) >= 2:
            continue
        chosen.append(cand)
        used_pos.add(cand["pos"])
        used_acid.add(cand["tgt"])
        if len(chosen) >= 4:
            break
    nac = len({c["tgt"] for c in chosen})
    return chosen, nac, round(sum(c["q"] for c in chosen), 3), len(props), epitope


# ---- AF2 for re-scoring the variants ----
from colabdesign import mk_afdesign_model, clear_mem

TARGET = "/kaggle/working/targets/EGFR_A_domainI.pdb"
HOTSPOT = "20,21,22,25,28,50,51"
_cache = {}


def score(seq):
    L = len(seq)
    if _cache.get("L") != L:
        clear_mem()
        af = mk_afdesign_model(protocol="binder", data_dir="/kaggle/working")
        af.prep_inputs(pdb_filename=TARGET, chain="A", binder_len=L, hotspot=HOTSPOT)
        _cache["L"], _cache["af"] = L, af
    af = _cache["af"]
    try:
        af.predict(seq=seq, num_recycles=3, verbose=False)
    except TypeError:
        af.predict(seq=seq, verbose=False)
    lg = af.aux["log"]
    return (round(float(lg.get("i_ptm", lg.get("ptm", 0.0))), 3),
            round(float(lg.get("plddt", 0.0)), 3))


records = []
for d in kept:
    if (time.time() - t0) / 60 > TIME_BUDGET_MIN:
        print("time budget reached", flush=True)
        break
    name, seq = d["name"], d["sequence"]
    try:
        chosen, nac, swsc, nprop, epi = place_switch(pdb_by[name], seq)
    except Exception as e:
        print("placement failed", name, type(e).__name__, str(e)[:200], flush=True)
        continue

    # Phase 1 selected F44,E45,D46,L49,Q52,Y74,D75 as the epitope. How much of the
    # observed interface actually falls there?
    on_target = sorted(set(epi) & {44, 45, 46, 49, 52, 74, 75})
    # Requirement 2 is mouse cross-reactivity, and what decides it is whether THIS
    # design's interface is conserved - not whether it sits on the patch we picked.
    ident = [CONSERVED.get(str(u - UNIPROT_OFFSET), 0) for u in epi]
    pct = 100.0 * sum(ident) / len(ident) if ident else 0.0
    nonid = [u for u, k in zip(epi, ident) if not k]
    print("%s i_ptm=%.3f epitope %s (%d res, %.0f%% conserved in mouse, "
          "substitutions at %s, %d on the Phase 1 patch: %s)"
          % (name, d["true_iptm"], "%d-%d" % (min(epi), max(epi)) if epi else "none",
             len(epi), pct, nonid or "none", len(on_target), on_target or "-"),
          flush=True)

    base = dict(epitope_residues=epi, epitope_on_target=on_target,
                epitope_conservation_pct=round(pct, 1),
                epitope_mouse_substitutions=nonid,
                generator="AfDesign", source_backbone=name,
                parent_i_ptm=d["true_iptm"], length=d["length"])
    records.append(dict(name=name + "_ctl", sequence=seq, role="control",
                        n_his_contacts=0, switch_score=0.0, his_substitutions="none",
                        i_ptm=d["true_iptm"], plddt=d["plddt"], iptm_retained=1.0, **base))

    if nac < 2:
        print("%s: %d candidate contacts, %d carboxylates -> control only"
              % (name, nprop, nac), flush=True)
        continue

    ss, subs = list(seq), []
    for c in chosen:
        j = c["pos"] - 1
        if 0 <= j < len(ss) and ss[j] != "H":
            subs.append("%s%dH" % (ss[j], c["pos"]))
            ss[j] = "H"
    var = "".join(ss)
    v_iptm, v_plddt = score(var)
    ret = round(v_iptm / d["true_iptm"], 3) if d["true_iptm"] else 0.0
    verdict = "KEEP" if ret >= KEEP_FRAC else "switch breaks the interface"
    print("%s: %d carboxylates, subs %s -> i_ptm %.3f vs parent %.3f (%.0f%%) %s [%.1f min]"
          % (name, nac, ",".join(subs), v_iptm, d["true_iptm"], 100 * ret, verdict,
             (time.time() - t0) / 60), flush=True)
    records.append(dict(name=name + "_sw", sequence=var, role="switch",
                        n_his_contacts=nac, switch_score=swsc,
                        his_substitutions=",".join(subs) or "none",
                        i_ptm=v_iptm, plddt=v_plddt, iptm_retained=ret, **base))

json.dump(records, open("/kaggle/working/designs.json", "w"), indent=1)
sw = [r for r in records if r["role"] == "switch"]
good = [r for r in sw if r["iptm_retained"] >= KEEP_FRAC and r["i_ptm"] >= 0.45]
print("\n%d records: %d switch (%d keeping >=%.0f%% of parent i_ptm), %d control"
      % (len(records), len(sw), len(good), 100 * KEEP_FRAC, len(records) - len(sw)))

print("\n===== COPY FROM HERE =====")
print("name\trole\tnHis\tswScore\tsubs\tplddt\ti_ptm\tparent\tretained\tsequence")
for r in sorted(records, key=lambda x: (-(x["role"] == "switch"), -x["i_ptm"])):
    print("\t".join(str(x) for x in [r["name"], r["role"], r["n_his_contacts"],
                                     r["switch_score"], r["his_substitutions"],
                                     r["plddt"], r["i_ptm"], r["parent_i_ptm"],
                                     r["iptm_retained"], r["sequence"]]))
print("===== COPY TO HERE =====")
