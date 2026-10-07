# === CHALLENGE 1 : SEQUENCE-LIABILITY REPAIR + SWITCH + FINAL SCORING ===
#
# Stage 5 settled the AfDesign-vs-ProteinMPNN question and the answer was not the one
# the literature led me to expect. MPNN produced clean, well-folded sequences - pLDDT
# rose to 0.82-0.87 - but it destroyed the interface on 12 of 14 backbones, taking
# i_ptm from 0.50-0.79 down to ~0.10. Only r2_w0_007 (0.718 -> 0.803) and r2_w0_009
# (0.695 -> 0.766) survived redesign, and placing histidines on r2_w0_007's MPNN
# sequence then cost it everything (0.803 -> 0.091).
#
# So wholesale redesign is the wrong instrument. The AfDesign sequences bind and take
# the pH switch; their problem is a short list of specific, local liabilities - free
# cysteines and homopolymer runs. This cell repairs those residues in place with
# conservative substitutions, protects the switch histidines, and re-predicts, so the
# cost of every repair is measured rather than assumed.

import glob, json, os, time
from collections import Counter
import numpy as np

MIN_IPTM = 0.45
KEEP_FRAC = 0.70
REPAIR_KEEP = 0.85       # a repair must hold this much of the unrepaired i_ptm
TIME_BUDGET_MIN = 260
t0 = time.time()

CB_MIN, CB_MAX, CB_IDEAL, COS_MIN = 3.5, 7.0, 5.8, -0.30
CARBOX = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}
ACIDS = [21, 22, 51, 73, 90, 110, 155]      # UniProt 45,46,75,97,114,134,179
NOMUT = {"PRO", "GLY", "CYS"}
PREF = {"K": 1.5, "R": 1.5, "Q": 1.2, "N": 1.2, "S": 1.15, "T": 1.15, "A": 1.0,
        "D": 1.1, "E": 1.1, "H": 0.0, "M": 0.7, "L": 0.6, "I": 0.6, "V": 0.6,
        "F": 0.4, "Y": 0.4, "W": 0.3}
MAXS = {"A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "E": 223, "Q": 225,
        "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240,
        "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174}
UNIPROT_OFFSET = 24

rj = glob.glob("/kaggle/input/**/round2.json", recursive=True)
if not rj:
    raise SystemExit("round2.json not found - attach the round2 notebook output.")
designs = json.load(open(rj[0]))
pdb_by = {}
for p in glob.glob("/kaggle/input/**/out2/*.pdb", recursive=True):
    pdb_by.setdefault(os.path.basename(p).replace(".pdb", ""), p)
kept = sorted([d for d in designs if d["true_iptm"] >= MIN_IPTM and d["name"] in pdb_by],
              key=lambda d: -d["true_iptm"])
print("%d designs, %d with i_ptm >= %.2f and a structure"
      % (len(designs), len(kept), MIN_IPTM), flush=True)

from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import protein_letters_3to1

SRC_TARGET = "/kaggle/working/targets/EGFR_A_domainI.pdb"
_src_chain = list(PDBParser(QUIET=True).get_structure("s", SRC_TARGET)[0])[0]
_src_res = [r for r in _src_chain if r.id[0] == " "]
_src_num = [r.id[1] for r in _src_res]
ACID_IDX = [_src_num.index(n) for n in ACIDS]
assert all(_src_res[i].get_resname() in CARBOX for i in ACID_IDX)
print("carboxylates %s at chain positions %s (%s)"
      % (ACIDS, ACID_IDX, [_src_res[i].get_resname() for i in ACID_IDX]), flush=True)


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


def analyse(pdb_path, seq):
    pr = PDBParser(QUIET=True)
    mdl = pr.get_structure("c", pdb_path)[0]
    sizes = Counter({ch.id: sum(1 for r in ch if r.id[0] == " ") for ch in mdl})
    tgt_ch = max(sizes, key=lambda k: sizes[k])
    bnd_ch = min(sizes, key=lambda k: sizes[k])
    tgt_res = [r for r in mdl[tgt_ch] if r.id[0] == " "]
    if len(tgt_res) != len(_src_res):
        raise ValueError("target chain %d residues vs source %d"
                         % (len(tgt_res), len(_src_res)))
    bfree = rsasa(pr.get_structure("b", pdb_path)[0][bnd_ch])
    tfree = rsasa(pr.get_structure("t", pdb_path)[0][tgt_ch])
    ShrakeRupley().compute(mdl, level="R")
    tcplx = {r.id[1]: r.sasa / MAXS.get(protein_letters_3to1.get(r.get_resname(), "A"), 200)
             for r in mdl[tgt_ch] if r.id[0] == " "
             and r.get_resname() in protein_letters_3to1}

    bnd_atoms = np.array([a.get_coord() for r in mdl[bnd_ch] if r.id[0] == " " for a in r])
    epitope = []
    for i, r in enumerate(tgt_res):
        ra = np.array([a.get_coord() for a in r])
        if np.min(np.linalg.norm(ra[:, None, :] - bnd_atoms[None, :, :], axis=-1)) <= 5.0:
            epitope.append(_src_num[i] + UNIPROT_OFFSET)

    acids = []
    for src_n, i in zip(ACIDS, ACID_IDX):
        r = tgt_res[i]
        ox = [r[a].get_coord() for a in CARBOX[r.get_resname()] if a in r]
        if ox:
            acids.append(dict(num=src_n, ox=ox,
                              des=max(tfree.get(r.id[1], 0) - tcplx.get(r.id[1], 0), 0.0)))

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
    return chosen, nac, round(sum(c["q"] for c in chosen), 3), epitope


# ---- liability repair -------------------------------------------------------
# Conservative, single-residue substitutions: same broad class, so the fold has the
# best chance of tolerating them. Measured afterwards either way.
SUB = {"C": "S", "L": "I", "I": "L", "V": "I", "E": "Q", "Q": "E", "K": "R",
       "R": "K", "P": "A", "A": "S", "S": "T", "T": "S", "G": "A", "N": "Q",
       "D": "N", "M": "L", "F": "Y", "Y": "F", "W": "F", "H": "H"}
RUN_MIN = 5


def repair(seq, protect=()):
    """Replace free cysteines and break homopolymer runs, leaving `protect` alone."""
    s, notes = list(seq), []
    protect = set(protect)
    for i, c in enumerate(s):
        if c == "C" and i not in protect:
            s[i] = "S"
            notes.append("C%dS" % (i + 1))
    i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]:
            j += 1
        if j - i + 1 >= RUN_MIN:
            for k in range(i + 2, j + 1, 2):
                if k in protect:
                    continue
                alt = SUB.get(s[k], "A")
                if alt == s[k]:
                    alt = "A"
                notes.append("%s%d%s" % (s[k], k + 1, alt))
                s[k] = alt
        i = j + 1
    return "".join(s), notes


def flags(s):
    runs, i = [], 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]:
            j += 1
        if j - i + 1 >= 4:
            runs.append("%s%d" % (s[i], j - i + 1))
        i = j + 1
    return s.count("C"), runs


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
        chosen, nac, swsc, epi = analyse(pdb_by[name], seq)
    except Exception as e:
        print("analysis failed", name, type(e).__name__, str(e)[:200], flush=True)
        continue

    base = dict(source_backbone=name, epitope_residues=epi,
                afdesign_i_ptm=d["true_iptm"], length=len(seq),
                generator="AfDesign binder hallucination (ColabDesign)")

    # --- control: repaired, no switch ---
    rep, notes = repair(seq)
    if rep == seq:
        c_iptm, c_plddt, c_notes = d["true_iptm"], d["plddt"], []
    else:
        c_iptm, c_plddt = score(rep)
        c_notes = notes
    cys, runs = flags(rep)
    keep_c = c_iptm >= REPAIR_KEEP * d["true_iptm"]
    print("%s CONTROL repair %s -> i_ptm %.3f vs %.3f (%.0f%%) cys=%d runs=%s %s"
          % (name, ",".join(c_notes) or "none needed", c_iptm, d["true_iptm"],
             100 * c_iptm / d["true_iptm"], cys, runs or "none",
             "OK" if keep_c else "REPAIR COSTLY"), flush=True)
    records.append(dict(name=name + "_ctl", sequence=rep, role="control",
                        n_his_contacts=0, switch_score=0.0, his_substitutions="none",
                        repair=",".join(c_notes) or "none", i_ptm=c_iptm, plddt=c_plddt,
                        iptm_retained=round(c_iptm / d["true_iptm"], 3),
                        free_cys=cys, homopolymer_runs=runs, **base))

    if nac < 2:
        print("   %s: %d carboxylates -> no switch variant" % (name, nac), flush=True)
        continue

    # --- switch: histidines first, then repair everything except them ---
    ss, subs, his_idx = list(seq), [], []
    for c in chosen:
        j = c["pos"] - 1
        if 0 <= j < len(ss) and ss[j] != "H":
            subs.append("%s%dH" % (ss[j], c["pos"]))
            ss[j] = "H"
            his_idx.append(j)
    sw = "".join(ss)
    sw_rep, sw_notes = repair(sw, protect=his_idx)
    s_iptm, s_plddt = score(sw_rep)
    cys, runs = flags(sw_rep)
    ret = round(s_iptm / d["true_iptm"], 3) if d["true_iptm"] else 0.0
    print("   %s SWITCH %s + repair %s -> i_ptm %.3f vs parent %.3f (%.0f%%) "
          "cys=%d runs=%s %s [%.1f min]"
          % (name, ",".join(subs), ",".join(sw_notes) or "none", s_iptm,
             d["true_iptm"], 100 * ret, cys, runs or "none",
             "KEEP" if ret >= KEEP_FRAC else "switch breaks interface",
             (time.time() - t0) / 60), flush=True)
    records.append(dict(name=name + "_sw", sequence=sw_rep, role="switch",
                        n_his_contacts=nac, switch_score=swsc,
                        his_substitutions=",".join(subs) or "none",
                        repair=",".join(sw_notes) or "none", i_ptm=s_iptm,
                        plddt=s_plddt, iptm_retained=ret,
                        free_cys=cys, homopolymer_runs=runs, **base))

json.dump(records, open("/kaggle/working/designs_final.json", "w"), indent=1)
sw = [r for r in records if r["role"] == "switch"]
good = [r for r in sw if r["iptm_retained"] >= KEEP_FRAC and r["i_ptm"] >= 0.45]
clean = [r for r in records if r["free_cys"] == 0 and not r["homopolymer_runs"]]
print("\n%d records: %d switch (%d usable), %d control; %d fully clean sequences"
      % (len(records), len(sw), len(good), len(records) - len(sw), len(clean)))

print("\n===== COPY FROM HERE =====")
for r in sorted(records, key=lambda x: (-(x["role"] == "switch"), -x["i_ptm"])):
    print(json.dumps(r, separators=(",", ":")))
print("===== COPY TO HERE =====")
