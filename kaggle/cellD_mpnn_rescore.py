# === CHALLENGE 1 : MPNN REDESIGN + AF2 RE-SCORE + SWITCH PLACEMENT ===
# AfDesign's hallucinated sequences lose interface confidence when discretised
# (i_ptm collapses at the hard stage). The established remedy, which BindCraft uses,
# is to resample sequences on the same backbone with ProteinMPNN and re-predict.
# A single AF2 prediction costs ~1 min vs ~18 min for a full design, so this is cheap.
#
# Produces /kaggle/working/designs.json ready for the submission pipeline.

NSEQ = 8                 # MPNN sequences per backbone
TIME_BUDGET_MIN = 300    # hard stop so the run always completes and saves
MIN_IPTM_BACKBONE = 0.25 # only redesign backbones worth the compute

import glob, json, os, subprocess, time
from collections import Counter
import numpy as np

t0 = time.time()

# Backbones live in the previous version's output, attached as a Kaggle input.
# Search both locations so the cell works whether re-run in place or from an input.
cands = (glob.glob("/kaggle/working/out/egfr_a_*.pdb")
         + glob.glob("/kaggle/input/**/out/egfr_a_*.pdb", recursive=True)
         + glob.glob("/kaggle/input/**/egfr_a_*.pdb", recursive=True))
seen, pdbs = set(), []
for c in sorted(cands):
    b = os.path.basename(c)
    if b not in seen and "traj" not in c:
        seen.add(b); pdbs.append(c)
print(len(pdbs), "backbone PDBs found", flush=True)
if not pdbs:
    raise SystemExit("No backbones. Attach the previous version's output as an input.")

# Round 1 i_ptm, from the generation log. Only redesign backbones worth the compute:
# a backbone whose interface AF2 cannot see at all is unlikely to be rescued by
# resampling its sequence, and each prediction costs ~1 min.
ROUND1 = {"egfr_a_15":0.725,"egfr_a_18":0.604,"egfr_a_9":0.601,"egfr_a_12":0.531,
          "egfr_a_20":0.477,"egfr_a_16":0.472,"egfr_a_8":0.462,"egfr_a_5":0.461,
          "egfr_a_7":0.442,"egfr_a_0":0.422,"egfr_a_4":0.416,"egfr_a_17":0.305,
          "egfr_a_11":0.289,"egfr_a_21":0.177,"egfr_a_3":0.153,"egfr_a_2":0.126,
          "egfr_a_13":0.124,"egfr_a_14":0.117,"egfr_a_23":0.117,"egfr_a_1":0.109,
          "egfr_a_6":0.108,"egfr_a_22":0.105,"egfr_a_10":0.102,"egfr_a_19":0.098}
pdbs = [p for p in pdbs
        if ROUND1.get(os.path.basename(p).replace(".pdb",""), 0) >= MIN_IPTM_BACKBONE]
pdbs.sort(key=lambda p: -ROUND1.get(os.path.basename(p).replace(".pdb",""), 0))
print("redesigning %d backbones with round-1 i_ptm >= %.2f (best first)"
      % (len(pdbs), MIN_IPTM_BACKBONE), flush=True)

# ---- parent i_ptm from the generation log, if the .seq files carry it ----
def chain_ids(p):
    c = Counter()
    for line in open(p):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            c[line[21]] += 1
    o = sorted(c.items(), key=lambda kv: -kv[1])
    return o[-1][0], o[0][0], dict(c)

BINDER_CH, TARGET_CH, counts = chain_ids(pdbs[0])
print("binder chain", BINDER_CH, "target chain", TARGET_CH, counts, flush=True)

# ---- 1. ProteinMPNN on each backbone ----
os.makedirs("/kaggle/working/seqs", exist_ok=True)
for p in pdbs:
    r = subprocess.run(["python", "/kaggle/working/ProteinMPNN/protein_mpnn_run.py",
                        "--pdb_path", p, "--pdb_path_chains", BINDER_CH,
                        "--out_folder", "/kaggle/working/seqs",
                        "--num_seq_per_target", str(NSEQ), "--sampling_temp", "0.1",
                        "--omit_AAs", "CX", "--use_soluble_model", "--seed", "37"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        print("MPNN failed on", os.path.basename(p), r.stderr[-400:], flush=True)
print("MPNN done [%.1f min]" % ((time.time() - t0) / 60), flush=True)

# ---- 2. AF2 re-score each sequence (single prediction, no design) ----
from colabdesign import mk_afdesign_model, clear_mem
TARGET = "/kaggle/working/targets/EGFR_A_domainI.pdb"
HOTSPOT = "20,21,22,25,28,50,51"

scores = {}
for fa in sorted(glob.glob("/kaggle/working/seqs/seqs/*.fa")):
    if (time.time() - t0) / 60 > TIME_BUDGET_MIN:
        print("time budget reached during scoring", flush=True); break
    stem = os.path.basename(fa).replace(".fa", "")
    if stem not in {os.path.basename(p).replace(".pdb","") for p in pdbs}:
        continue
    ent, hdr = [], None
    for line in open(fa):
        line = line.strip()
        if line.startswith(">"): hdr = line
        elif line: ent.append((hdr, line))
    seqs = ent[1:]
    if not seqs:
        continue
    L = len(seqs[0][1])
    try:
        clear_mem()
        af = mk_afdesign_model(protocol="binder", data_dir="/kaggle/working")
        af.prep_inputs(pdb_filename=TARGET, chain="A", binder_len=L, hotspot=HOTSPOT)
    except Exception as e:
        print("prep failed", stem, str(e)[:200], flush=True); continue
    for i, (h, s) in enumerate(seqs, start=1):
        name = stem + "_s" + str(i)
        try:
            af.predict(seq=s, verbose=False)
            log = af.aux["log"]
            scores[name] = dict(sequence=s, length=len(s),
                                plddt=round(float(log.get("plddt", 0)), 3),
                                i_ptm=round(float(log.get("i_ptm", log.get("ptm", 0))), 3),
                                backbone=stem)
        except Exception as e:
            print("predict failed", name, str(e)[:150], flush=True)
    done = [v for k, v in scores.items() if v["backbone"] == stem]
    if done:
        print("%s: best i_ptm %.3f of %d seqs [%.1f min]"
              % (stem, max(d["i_ptm"] for d in done), len(done),
                 (time.time() - t0) / 60), flush=True)

json.dump(scores, open("/kaggle/working/rescored.json", "w"), indent=1)
print("\nre-scored %d sequences [%.1f min]" % (len(scores), (time.time() - t0) / 60),
      flush=True)

# ---- 3. acid-ON switch placement (CPU, calibrated geometry) ----
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import protein_letters_3to1

CB_MIN, CB_MAX, CB_IDEAL, COS_MIN = 3.5, 7.0, 5.8, -0.30
CARBOX = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}
ACIDS = [21, 22, 51, 73]            # E45, D46, D75, E97 (UniProt - 24)
NOMUT = {"PRO", "GLY", "CYS"}
PREF = {"K": 1.5, "R": 1.5, "Q": 1.2, "N": 1.2, "S": 1.15, "T": 1.15, "A": 1.0,
        "D": 1.1, "E": 1.1, "H": 0.0, "M": 0.7, "L": 0.6, "I": 0.6, "V": 0.6,
        "F": 0.4, "Y": 0.4, "W": 0.3}
MAXS = {"A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "E": 223, "Q": 225,
        "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240,
        "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174}

def vcb(r):
    if "CB" in r: return r["CB"].get_coord()
    try: n, a, c = r["N"].get_coord(), r["CA"].get_coord(), r["C"].get_coord()
    except KeyError: return None
    b, cc = a - n, c - a
    return -0.58273431 * np.cross(b, cc) + 0.56802827 * b - 0.54067466 * cc + a

def rsasa(e):
    ShrakeRupley().compute(e, level="R")
    d = {}
    for r in e.get_residues():
        aa = protein_letters_3to1.get(r.get_resname())
        if r.id[0] == " " and aa: d[r.id[1]] = r.sasa / MAXS.get(aa, 200)
    return d

records = []
for p in pdbs:
    stem = os.path.basename(p).replace(".pdb", "")
    mine = {k: v for k, v in scores.items() if v["backbone"] == stem}
    if not mine:
        continue
    pr = PDBParser(QUIET=True)
    mdl = pr.get_structure("c", p)[0]
    bfree = rsasa(pr.get_structure("b", p)[0][BINDER_CH])
    tfree = rsasa(pr.get_structure("t", p)[0][TARGET_CH])
    ShrakeRupley().compute(mdl, level="R")
    tcplx = {r.id[1]: r.sasa / MAXS.get(protein_letters_3to1.get(r.get_resname(), "A"), 200)
             for r in mdl[TARGET_CH] if r.id[0] == " "
             and r.get_resname() in protein_letters_3to1}
    acids = []
    for r in mdl[TARGET_CH]:
        if r.id[0] != " " or r.get_resname() not in CARBOX or r.id[1] not in ACIDS:
            continue
        ox = [r[a].get_coord() for a in CARBOX[r.get_resname()] if a in r]
        if ox:
            acids.append(dict(num=r.id[1], ox=ox,
                              des=max(tfree.get(r.id[1], 0) - tcplx.get(r.id[1], 0), 0.0)))
    props = []
    for r in mdl[BINDER_CH]:
        if r.id[0] != " " or "CA" not in r or r.get_resname() in NOMUT: continue
        if bfree.get(r.id[1], 0) < 0.15: continue
        cb = vcb(r)
        if cb is None: continue
        d = cb - r["CA"].get_coord(); nn = np.linalg.norm(d)
        if nn < 1e-6: continue
        d = d / nn
        for ac in acids:
            ds = [np.linalg.norm(cb - o) for o in ac["ox"]]
            dcb = float(min(ds))
            if not (CB_MIN <= dcb <= CB_MAX): continue
            to = ac["ox"][int(np.argmin(ds))] - cb; to = to / np.linalg.norm(to)
            cs = float(np.dot(d, to))
            if cs < COS_MIN: continue
            q = max(0.0, 1 - abs(dcb - CB_IDEAL) / max(CB_MAX - CB_IDEAL, CB_IDEAL - CB_MIN))
            props.append(dict(pos=r.id[1], tgt=ac["num"], q=q * (1 + ac["des"])))

    best = sorted(mine.items(), key=lambda kv: -kv[1]["i_ptm"])
    rep = best[0][1]["sequence"]
    chosen, up, ua = [], set(), set()
    for pr2 in sorted(props, key=lambda x: -(x["q"] * PREF.get(rep[x["pos"] - 1], 0.9)
                                             if 0 <= x["pos"] - 1 < len(rep) else 0)):
        if pr2["pos"] in up: continue
        if 0 <= pr2["pos"] - 1 < len(rep) and PREF.get(rep[pr2["pos"] - 1], 0.9) <= 0: continue
        if pr2["tgt"] in ua and len(chosen) >= 2: continue
        chosen.append(pr2); up.add(pr2["pos"]); ua.add(pr2["tgt"])
        if len(chosen) >= 4: break
    nac = len({c["tgt"] for c in chosen})
    swsc = round(sum(c["q"] for c in chosen), 3)
    print("%s: %d contacts, %d carboxylates -> %s"
          % (stem, len(props), nac, "SWITCH" if nac >= 2 else "control-only"), flush=True)

    for rank, (name, v) in enumerate(best):
        base = dict(epitope="A_domainI_ligand_face", epitope_conservation=1.0,
                    generator="AfDesign+ProteinMPNN", source_backbone=stem,
                    plddt=v["plddt"], i_ptm=v["i_ptm"])
        if nac >= 2:
            ss, app = list(v["sequence"]), []
            for c in chosen:
                j = c["pos"] - 1
                if 0 <= j < len(ss) and ss[j] != "H":
                    app.append("%s%dH" % (ss[j], c["pos"])); ss[j] = "H"
            records.append(dict(name=name + "_sw", sequence="".join(ss), role="switch",
                                n_his_contacts=nac, switch_score=swsc,
                                his_substitutions=",".join(app) or "none", **base))
        if rank == 0:
            records.append(dict(name=name + "_ctl", sequence=v["sequence"],
                                role="control", n_his_contacts=0, switch_score=0.0,
                                his_substitutions="none", **base))

json.dump(records, open("/kaggle/working/designs.json", "w"), indent=1)
nsw = sum(1 for r in records if r["role"] == "switch")
print("\n%d records (%d switch, %d control) -> designs.json [%.1f min]"
      % (len(records), nsw, len(records) - nsw, (time.time() - t0) / 60))

print("\n===== COPY FROM HERE =====")
print("name\trole\tnHis\tswScore\tsubs\tplddt\ti_ptm\tsequence")
for r in sorted(records, key=lambda x: (-x["n_his_contacts"], -(x.get("i_ptm") or 0)))[:60]:
    print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (r["name"], r["role"], r["n_his_contacts"],
          r["switch_score"], r["his_substitutions"], r.get("plddt"), r.get("i_ptm"),
          r["sequence"]))
print("===== COPY TO HERE =====")
