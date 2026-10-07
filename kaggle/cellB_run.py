# === CHALLENGE 1 : GENERATE + SCORE + SWITCH ANALYSIS ===
# Set NUM = 2 for a smoke test, then NUM = 60 for the real run.
NUM = 2
PREFIX = "/kaggle/working/out/egfr_a"

import glob, json, os, subprocess, time
from collections import Counter
import numpy as np

os.environ["DGLBACKEND"] = "pytorch"
t0 = time.time()

# ---------------- 1. RFdiffusion ----------------
os.makedirs("/kaggle/working/out", exist_ok=True)
HOTSPOTS = "[A20,A21,A22,A25,A28,A50,A51]"   # F44,E45,D46,L49,Q52,Y74,D75 (UniProt-24)
CONTIGS = "[A1-164/0 55-95]"
cmd = ("cd /kaggle/working/RFdiffusion && python scripts/run_inference.py"
       " inference.output_prefix=" + PREFIX +
       " inference.input_pdb=/kaggle/working/targets/EGFR_A_domainI.pdb"
       " 'contigmap.contigs=" + CONTIGS + "'"
       " 'ppi.hotspot_res=" + HOTSPOTS + "'"
       " inference.num_designs=" + str(NUM) +
       " denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0")
print(cmd, flush=True)
os.system(cmd)
pdbs = sorted(p for p in glob.glob("/kaggle/working/out/egfr_a*.pdb") if "traj" not in p)
print("\nbackbones:", len(pdbs), "| %.1f min" % ((time.time() - t0) / 60), flush=True)
assert pdbs, "RFdiffusion produced no backbones"

# ---------------- 2. ProteinMPNN ----------------
def chain_counts(p):
    c = Counter()
    for line in open(p):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            c[line[21]] += 1
    return c

cc = chain_counts(pdbs[0])
order = sorted(cc.items(), key=lambda kv: -kv[1])
TARGET_CH, BINDER_CH = order[0][0], order[-1][0]
print("chains:", dict(cc), "-> target", TARGET_CH, "binder", BINDER_CH, flush=True)

os.makedirs("/kaggle/working/seqs", exist_ok=True)
for p in pdbs:
    r = subprocess.run(["python", "/kaggle/working/ProteinMPNN/protein_mpnn_run.py",
                        "--pdb_path", p, "--pdb_path_chains", BINDER_CH,
                        "--out_folder", "/kaggle/working/seqs",
                        "--num_seq_per_target", "8", "--sampling_temp", "0.1",
                        "--omit_AAs", "CX", "--use_soluble_model", "--seed", "37"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        print("MPNN FAILED on", os.path.basename(p)); print(r.stderr[-800:]); break
print("sequences done | %.1f min" % ((time.time() - t0) / 60), flush=True)

# ---------------- 3. ESMFold self-consistency ----------------
import torch
from transformers import AutoTokenizer, EsmForProteinFolding

def ca(path, chain=None):
    out = []
    for line in open(path):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            if chain and line[21] != chain:
                continue
            out.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(out)

def rmsd(P, Q):
    if len(P) != len(Q) or len(P) == 0:
        return float("nan")
    P = P - P.mean(0); Q = Q - Q.mean(0)
    V, S, W = np.linalg.svd(P.T @ Q)
    if np.linalg.det(V) * np.linalg.det(W) < 0:
        V[:, -1] = -V[:, -1]
    P = P @ (V @ W)
    return float(np.sqrt(((P - Q) ** 2).sum() / len(P)))

print("loading ESMFold ...", flush=True)
tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
esm = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1",
                                           low_cpu_mem_usage=True).cuda().eval()
esm.esm = esm.esm.half()
os.makedirs("/kaggle/working/folded", exist_ok=True)

quality = {}
for fa in sorted(glob.glob("/kaggle/working/seqs/seqs/*.fa")):
    stem = os.path.basename(fa).replace(".fa", "")
    parent = "/kaggle/working/out/" + stem + ".pdb"
    if not os.path.exists(parent):
        continue
    ref = ca(parent, BINDER_CH)
    ent, hdr = [], None
    for line in open(fa):
        line = line.strip()
        if line.startswith(">"): hdr = line
        elif line: ent.append((hdr, line))
    for i, (h, s) in enumerate(ent[1:], start=1):
        name = stem + "_s" + str(i)
        with torch.no_grad():
            tt = tok([s], return_tensors="pt", add_special_tokens=False)
            out = esm(**{k: v.cuda() for k, v in tt.items()})
        open("/kaggle/working/folded/%s.pdb" % name, "w").write(esm.output_to_pdb(out)[0])
        p = float(out["plddt"][0, :, 1].mean())
        p = p * 100 if p <= 1.0 else (p / 100 if p > 100 else p)   # normalise scale
        rr = rmsd(ca("/kaggle/working/folded/%s.pdb" % name), ref)
        quality[name] = dict(plddt=round(p, 1), sc_rmsd=round(rr, 2),
                             designable=bool(p >= 80 and rr <= 2.0))
print("folded %d sequences | %.1f min" % (len(quality), (time.time() - t0) / 60), flush=True)
ok = sum(1 for v in quality.values() if v["designable"])
print("designable (pLDDT>=80, scRMSD<=2):", ok, "/", len(quality), flush=True)

# ---------------- 4. acid-ON switch placement ----------------
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import protein_letters_3to1

CB_MIN, CB_MAX, CB_IDEAL, COS_MIN = 3.5, 7.0, 5.8, -0.30
CARBOX = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}
ACIDS = [21, 22, 51, 73]          # E45, D46, D75, E97
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
    b, cc2 = a - n, c - a
    return -0.58273431 * np.cross(b, cc2) + 0.56802827 * b - 0.54067466 * cc2 + a

def rsasa(e):
    ShrakeRupley().compute(e, level="R")
    d = {}
    for r in e.get_residues():
        aa = protein_letters_3to1.get(r.get_resname())
        if r.id[0] == " " and aa: d[r.id[1]] = r.sasa / MAXS.get(aa, 200)
    return d

records = []
for fa in sorted(glob.glob("/kaggle/working/seqs/seqs/*.fa")):
    stem = os.path.basename(fa).replace(".fa", "")
    pdb = "/kaggle/working/out/" + stem + ".pdb"
    if not os.path.exists(pdb): continue
    ent, hdr = [], None
    for line in open(fa):
        line = line.strip()
        if line.startswith(">"): hdr = line
        elif line: ent.append((hdr, line))
    seqs = ent[1:]
    if not seqs: continue

    pr = PDBParser(QUIET=True)
    mdl = pr.get_structure("c", pdb)[0]
    bfree = rsasa(pr.get_structure("b", pdb)[0][BINDER_CH])
    tfree = rsasa(pr.get_structure("t", pdb)[0][TARGET_CH])
    ShrakeRupley().compute(mdl, level="R")
    tcplx = {r.id[1]: r.sasa / MAXS.get(protein_letters_3to1.get(r.get_resname(), "A"), 200)
             for r in mdl[TARGET_CH] if r.id[0] == " " and r.get_resname() in protein_letters_3to1}

    acids = []
    for r in mdl[TARGET_CH]:
        if r.id[0] != " " or r.get_resname() not in CARBOX or r.id[1] not in ACIDS: continue
        ox = [r[a].get_coord() for a in CARBOX[r.get_resname()] if a in r]
        if ox: acids.append(dict(num=r.id[1], aa=protein_letters_3to1[r.get_resname()], ox=ox,
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

    rep = seqs[0][1]
    chosen, up, ua = [], set(), set()
    for p in sorted(props, key=lambda x: -(x["q"] * PREF.get(rep[x["pos"] - 1], 0.9)
                                           if 0 <= x["pos"] - 1 < len(rep) else 0)):
        if p["pos"] in up: continue
        if PREF.get(rep[p["pos"] - 1], 0.9) <= 0: continue
        if p["tgt"] in ua and len(chosen) >= 2: continue
        chosen.append(p); up.add(p["pos"]); ua.add(p["tgt"])
        if len(chosen) >= 4: break
    nac = len({c["tgt"] for c in chosen})
    sw_score = round(sum(c["q"] for c in chosen), 3)
    print("%s: %d contacts, %d carboxylates -> %s"
          % (stem, len(props), nac, "SWITCH" if nac >= 2 else "control-only"), flush=True)

    for i, (h, s) in enumerate(seqs, start=1):
        nm = stem + "_s" + str(i)
        qd = quality.get(nm, {})
        base = dict(epitope="A_domainI_ligand_face", epitope_conservation=1.0,
                    generator="RFdiffusion", source_backbone=stem,
                    plddt=qd.get("plddt"), sc_rmsd=qd.get("sc_rmsd"),
                    designable=qd.get("designable"))
        if nac >= 2:
            ss, app = list(s), []
            for c in chosen:
                j = c["pos"] - 1
                if 0 <= j < len(ss) and ss[j] != "H":
                    app.append("%s%dH" % (ss[j], c["pos"])); ss[j] = "H"
            records.append(dict(name=nm + "_sw", sequence="".join(ss), role="switch",
                                n_his_contacts=nac, switch_score=sw_score,
                                his_substitutions=",".join(app) or "none", **base))
        if i == 1:
            records.append(dict(name=nm + "_ctl", sequence=s, role="control",
                                n_his_contacts=0, switch_score=0.0,
                                his_substitutions="none", **base))

json.dump(records, open("/kaggle/working/designs.json", "w"), indent=1)
nsw = sum(1 for r in records if r["role"] == "switch")
print("\n%d records (%d switch, %d control) | total %.1f min"
      % (len(records), nsw, len(records) - nsw, (time.time() - t0) / 60))

print("\n===== COPY FROM HERE =====")
print("name\trole\tnHis\tswScore\tsubs\tplddt\tscRMSD\tsequence")
for r in sorted(records, key=lambda x: (-x["n_his_contacts"], x.get("sc_rmsd") or 99)):
    print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (r["name"], r["role"], r["n_his_contacts"],
          r["switch_score"], r["his_substitutions"], r.get("plddt"), r.get("sc_rmsd"),
          r["sequence"]))
print("===== COPY TO HERE =====")
