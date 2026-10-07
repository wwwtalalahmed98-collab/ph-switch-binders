# Add as a NEW cell at the end of the Kaggle notebook.
#
# Runs the Phase 3 acid-ON switch analysis in place, so the heavy PDBs never have to
# leave Kaggle. Writes a small designs.json (a few KB) and prints a compact table that
# can simply be copy-pasted if the download is awkward.
#
# Mechanism: histidines on the binder are placed against CONSERVED target carboxylates.
# At pH 6.5 the protonated imidazolium pairs with the carboxylate (binding); at pH 7.4
# the neutral imidazole loses that contact, leaving a desolvated unpaired carboxylate.
#
# Geometry thresholds were calibrated on 23 genuine His-carboxylate salt bridges across
# 12 PDB structures: CB->carboxylate-O 3.55-6.77 A (mean 5.79), direction cosine
# -0.33..0.81 (p5 -0.26), N-O contact 2.52-3.27 A (median 2.84).

import glob
import json
import os
from collections import Counter

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import protein_letters_3to1

# ---- calibrated constants -------------------------------------------------------
CB_O_MIN, CB_O_MAX, CB_O_IDEAL = 3.5, 7.0, 5.8
COS_MIN = -0.30
MAX_HIS, MIN_CONTACTS = 4, 2

# Epitope A carboxylates, PDB numbering (UniProt - 24): E45, D46, D75, E97
TARGET_CARBOXYLATES = [21, 22, 51, 73]
EPITOPE_NAME = "A_domainI_ligand_face"

MAX_SASA = {"A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "E": 223, "Q": 225,
            "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240,
            "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174}
CARBOXYL = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}
NO_MUTATE = {"PRO", "GLY", "CYS"}
# Lys/Arg: converts a permanent salt bridge into a pH-dependent one -> preferred.
# Aromatics/hydrophobics: substituting costs packing and may be load-bearing.
PREF = {"K": 1.5, "R": 1.5, "Q": 1.2, "N": 1.2, "S": 1.15, "T": 1.15, "A": 1.0,
        "D": 1.1, "E": 1.1, "H": 0.0, "M": 0.7, "L": 0.6, "I": 0.6, "V": 0.6,
        "F": 0.4, "Y": 0.4, "W": 0.3}


def virtual_cb(r):
    if "CB" in r:
        return r["CB"].get_coord()
    try:
        n, ca, c = r["N"].get_coord(), r["CA"].get_coord(), r["C"].get_coord()
    except KeyError:
        return None
    b, cc = ca - n, c - ca
    a = np.cross(b, cc)
    return -0.58273431 * a + 0.56802827 * b - 0.54067466 * cc + ca


def rel_sasa(entity):
    ShrakeRupley().compute(entity, level="R")
    out = {}
    for r in entity.get_residues():
        if r.id[0] != " ":
            continue
        aa = protein_letters_3to1.get(r.get_resname())
        if aa:
            out[r.id[1]] = r.sasa / MAX_SASA.get(aa, 200)
    return out


def quality(d):
    span = max(CB_O_MAX - CB_O_IDEAL, CB_O_IDEAL - CB_O_MIN)
    return max(0.0, 1.0 - abs(d - CB_O_IDEAL) / span)


def chains_of(pdb):
    c = Counter()
    for line in open(pdb):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            c[line[21]] += 1
    if len(c) < 2:
        return None, None
    o = sorted(c.items(), key=lambda kv: -kv[1])
    return o[-1][0], o[0][0]          # binder (smallest), target (largest)


def analyse(pdb, binder_ch, target_ch):
    p = PDBParser(QUIET=True)
    model = p.get_structure("c", pdb)[0]
    b_free = rel_sasa(p.get_structure("b", pdb)[0][binder_ch])
    t_free = rel_sasa(p.get_structure("t", pdb)[0][target_ch])
    ShrakeRupley().compute(model, level="R")
    t_cplx = {r.id[1]: r.sasa / MAX_SASA.get(
        protein_letters_3to1.get(r.get_resname(), "A"), 200)
        for r in model[target_ch] if r.id[0] == " "
        and r.get_resname() in protein_letters_3to1}

    acids = []
    for r in model[target_ch]:
        if r.id[0] != " " or r.get_resname() not in CARBOXYL:
            continue
        if r.id[1] not in TARGET_CARBOXYLATES:
            continue
        oxys = [r[a].get_coord() for a in CARBOXYL[r.get_resname()] if a in r]
        if oxys:
            acids.append(dict(num=r.id[1],
                              aa=protein_letters_3to1[r.get_resname()], oxys=oxys,
                              desolv=max(t_free.get(r.id[1], 0)
                                         - t_cplx.get(r.id[1], 0), 0.0)))

    props = []
    for r in model[binder_ch]:
        if r.id[0] != " " or "CA" not in r or r.get_resname() in NO_MUTATE:
            continue
        if b_free.get(r.id[1], 0) < 0.15:
            continue
        cb, ca = virtual_cb(r), r["CA"].get_coord()
        if cb is None:
            continue
        d = cb - ca
        nrm = np.linalg.norm(d)
        if nrm < 1e-6:
            continue
        d = d / nrm
        for ac in acids:
            dists = [np.linalg.norm(cb - o) for o in ac["oxys"]]
            dcb = float(min(dists))
            if not (CB_O_MIN <= dcb <= CB_O_MAX):
                continue
            to = ac["oxys"][int(np.argmin(dists))] - cb
            to = to / np.linalg.norm(to)
            cos = float(np.dot(d, to))
            if cos < COS_MIN:
                continue
            q = quality(dcb)
            props.append(dict(pos=r.id[1], target=ac["num"], target_aa=ac["aa"],
                              dcb=round(dcb, 2), cos=round(cos, 2),
                              score=round(q * (1 + ac["desolv"]), 3)))
    return props


def pick(props, seq, offset=1):
    def adj(p):
        i = p["pos"] - offset
        return p["score"] * (PREF.get(seq[i], 0.9) if 0 <= i < len(seq) else 0.9)
    chosen, used_p, used_a = [], set(), set()
    for p in sorted(props, key=lambda x: -adj(x)):
        if adj(p) <= 0 or p["pos"] in used_p:
            continue
        if p["target"] in used_a and len(chosen) >= 2:
            continue
        chosen.append(p)
        used_p.add(p["pos"])
        used_a.add(p["target"])
        if len(chosen) >= MAX_HIS:
            break
    return chosen


# ---- load quality metrics -------------------------------------------------------
sc = {}
if os.path.exists("/kaggle/working/selfconsistency.json"):
    for r in json.load(open("/kaggle/working/selfconsistency.json")):
        sc[r["name"]] = r
    print(f"loaded {len(sc)} self-consistency records "
          f"({sum(1 for r in sc.values() if r.get('designable'))} designable)")

# ---- run ------------------------------------------------------------------------
records = []
for fa in sorted(glob.glob("/kaggle/working/seqs/seqs/*.fa")):
    stem = os.path.basename(fa).replace(".fa", "")
    pdb = f"/kaggle/working/out/{stem}.pdb"
    if not os.path.exists(pdb):
        continue
    bch, tch = chains_of(pdb)
    if not bch:
        print("skip (chains)", stem)
        continue

    entries, hdr = [], None
    for line in open(fa):
        line = line.strip()
        if line.startswith(">"):
            hdr = line
        elif line:
            entries.append((hdr, line))
    seqs = entries[1:]                      # entry 0 is the input backbone
    if not seqs:
        continue

    props = analyse(pdb, bch, tch)
    chosen = pick(props, seqs[0][1])
    n_acids = len({c["target"] for c in chosen})
    sw_score = round(sum(c["score"] for c in chosen), 3)
    positions = [c["pos"] for c in chosen]
    print(f"{stem}: {len(props)} contacts, {n_acids} distinct carboxylates "
          f"-> {'SWITCH' if n_acids >= MIN_CONTACTS else 'control-only'}")

    for i, (h, s) in enumerate(seqs, start=1):
        name = f"{stem}_s{i}"
        q = sc.get(name, {})
        base = dict(epitope=EPITOPE_NAME, epitope_conservation=1.0,
                    generator="RFdiffusion", source_backbone=stem,
                    plddt=q.get("plddt"), sc_rmsd=q.get("sc_rmsd"),
                    mpnn_score=q.get("mpnn_score"), designable=q.get("designable"))
        if n_acids >= MIN_CONTACTS:
            ss, applied = list(s), []
            for pos in positions:
                j = pos - 1
                if 0 <= j < len(ss) and ss[j] != "H":
                    applied.append(f"{ss[j]}{pos}H")
                    ss[j] = "H"
            records.append(dict(name=f"{name}_sw", sequence="".join(ss),
                                role="switch", n_his_contacts=n_acids,
                                switch_score=sw_score,
                                his_substitutions=",".join(applied) or "none",
                                **base))
        if i == 1:
            records.append(dict(name=f"{name}_ctl", sequence=s, role="control",
                                n_his_contacts=0, switch_score=0.0,
                                his_substitutions="none", **base))

out = "/kaggle/working/designs.json"
json.dump(records, open(out, "w"), indent=1)
nsw = sum(1 for r in records if r["role"] == "switch")
print(f"\n{len(records)} records ({nsw} switch, {len(records)-nsw} control) -> {out}")
print(f"size: {os.path.getsize(out)/1024:.1f} KB")

# ---- compact paste-able table ---------------------------------------------------
print("\n===== COPY FROM HERE =====")
print("name\trole\tnHis\tswScore\tsubs\tplddt\tscRMSD\tsequence")
for r in sorted(records, key=lambda x: (-x["n_his_contacts"],
                                        x.get("sc_rmsd") or 99)):
    print(f"{r['name']}\t{r['role']}\t{r['n_his_contacts']}\t{r['switch_score']}\t"
          f"{r['his_substitutions']}\t{r.get('plddt')}\t{r.get('sc_rmsd')}\t"
          f"{r['sequence']}")
print("===== COPY TO HERE =====")
