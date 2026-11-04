"""Kaggle - Challenge 5 phase 2. Generate, with the two hotspot sets run head to head.

WHY THIS IS NOT A SEPARATE PROBE. At roughly 19 s per optimisation step on this target,
one design costs about 18 minutes. A forty-minute probe buys two designs, one per arm,
which is noise. So the arm comparison runs inside the real generation run instead: both
arms, one per GPU, same budget, same seeds policy. Whichever arm wins, its designs are
keepable output rather than diagnostics, so answering the question costs nothing extra.

    arm A   hotspots C5,C6,C7   num_pos 3   GPU 0   the phase 1 choice
    arm B   hotspots C4..C9     num_pos 6   GPU 1   the obvious set

THE QUESTION. Asp6 is 11% exposed; Gly7 beside it is 73% exposed and identical in the
wild type. Phase 1 argued arm B would be satisfied by Gly7 and never read the mutation.
This measures it.

THE METRIC IS on_mutation, NOT on_peptide. A design with twenty peptide contacts and
none on residue 6 reads nothing the wild type lacks and is worth zero. pmhc/test_measure.py
proves the measurement can tell those two cases apart, against the exact chain layout
save_pdb produces.

FOUR FAILURES THIS CARRIES GUARDS FOR, ALL OF WHICH ALREADY HAPPENED
  1. i_con's default cutoff 21.69 accepts contacts out to 21.4 A. TfR1 round 1 aimed at
     nothing. Set to 8.0, and the worker REFUSES TO RUN if the setting did not take.
  2. No AlphaFold params -> "no model params defined", 47 s in. Fetched and counted here.
  3. save_pdb MERGES the two target chains into one. Probe v2 died on it after 32
     minutes. measure() handles the merged layout and verifies by sequence.
  4. save_pdb renumbers from 1, so nothing is identified by residue number anywhere.

Records are appended per design, so a crash loses one design and not the run.
"""
import json
import os
import subprocess
import sys
import time

W = "/kaggle/working"
BUDGET_MIN = 320          # per worker; both run concurrently
LEN_MIN, LEN_MAX = 65, 95
SOFT, TEMP, HARD = 35, 15, 5      # 55 steps, about 18 min per design

PEP_G12D = "VVVGADGVGK"
DESIGN_RADIUS = 20.0
MIN_SEGMENT = 3
EXPECTED_A = ("SHSMRYFYTSVSRPGGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDQETRNVKAQSQ"
              "TDRVDLGTLRGYYNQSEDGSHTIQIMYGCDVGPDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQI"
              "TKRKWEAAHAAEQQRAYLEGRCVEWLRRYLENGKETLQR")

ARMS = {"A_tight": dict(hotspot=[5, 6, 7], num_pos=3, gpu=0, seed=7200),
        "B_broad": dict(hotspot=[4, 5, 6, 7, 8, 9], num_pos=6, gpu=1, seed=7400)}


def sh(c):
    subprocess.run(c, shell=True, check=True)


# ---------------------------------------------------------------- setup (parent)
try:
    import colabdesign  # noqa: F401
    print("colabdesign already present")
except ImportError:
    sh("pip -q install git+https://github.com/sokrypton/ColabDesign.git@v1.1.1 biopython")

import glob  # noqa: E402

import numpy as np  # noqa: E402
from Bio import PDB  # noqa: E402

if not os.path.exists("%s/params/params_model_1_ptm.npz" % W):
    os.makedirs("%s/params" % W, exist_ok=True)
    sh("wget -q -O %s/af2.tar https://storage.googleapis.com/alphafold/"
       "alphafold_params_2022-12-06.tar" % W)
    sh("tar -xf %s/af2.tar -C %s/params" % (W, W))
    os.remove("%s/af2.tar" % W)
n_params = len(glob.glob("%s/params/params_model_*.npz" % W))
print("params: %d files" % n_params)
if n_params == 0:
    raise SystemExit("no AlphaFold parameters - refusing to start")

T3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
      "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
      "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
      "TYR": "Y", "VAL": "V"}


def residues(chain):
    out, seen = [], set()
    for r in chain:
        if r.id[0] == " " and r.get_resname() in T3 and r.id[1] not in seen:
            seen.add(r.id[1])
            out.append(r)
    return out


def seq(rs):
    return "".join(T3[r.get_resname()] for r in rs)


def build_target(path):
    import urllib.request
    for pid in ("7OW6", "7OW5"):
        f = "%s/%s.cif" % (W, pid)
        if not os.path.exists(f):
            urllib.request.urlretrieve(
                "https://files.rcsb.org/download/%s.cif" % pid, f)
    p = PDB.MMCIFParser(QUIET=True)
    m6 = p.get_structure("6", "%s/7OW6.cif" % W)[0]
    m5 = p.get_structure("5", "%s/7OW5.cif" % W)[0]
    hv, pep = residues(m6["A"]), residues(m6["C"])
    if seq(pep) != PEP_G12D:
        raise SystemExit("7OW6 chain C is %s" % seq(pep))
    allowed = {r.id[1] for r in hv} & {r.id[1] for r in residues(m5["A"])}
    pep_atoms = np.concatenate([np.array([a.coord for a in r if a.element != "H"])
                                for r in pep])
    near = [r.id[1] for r in hv if r.id[1] in allowed
            and np.linalg.norm(
                np.array([a.coord for a in r if a.element != "H"])[:, None, :]
                - pep_atoms[None, :, :], axis=-1).min() <= DESIGN_RADIUS]
    runs, cur = [], []
    for n in sorted(near):
        if cur and n == cur[-1] + 1:
            cur.append(n)
        else:
            if len(cur) >= MIN_SEGMENT:
                runs += cur
            cur = [n]
    if len(cur) >= MIN_SEGMENT:
        runs += cur
    keep = {"A": set(runs), "C": {r.id[1] for r in pep}}

    class Sel(PDB.Select):
        def accept_model(self, m):
            return m.id == 0

        def accept_chain(self, c):
            return c.id in keep

        def accept_residue(self, r):
            return r.id[0] == " " and r.id[1] in keep[r.get_parent().id]

        def accept_atom(self, a):
            return (not a.is_disordered()) or a.get_altloc() in (" ", "A")

    io = PDB.PDBIO()
    io.set_structure(m6)
    io.save(path, Sel())
    chk = PDB.PDBParser(QUIET=True).get_structure("t", path)[0]
    if seq(residues(chk["A"])) != EXPECTED_A or seq(residues(chk["C"])) != PEP_G12D:
        raise SystemExit("rebuilt target does not match phase 0 - stopping")
    print("target rebuilt and verified: chain A %d aa, peptide %s"
          % (len(EXPECTED_A), PEP_G12D))
    return path


TARGET = build_target("%s/target_design.pdb" % W)
json.dump(dict(arms=ARMS, expected_a=EXPECTED_A, peptide=PEP_G12D,
               soft=SOFT, temp=TEMP, hard=HARD, len_min=LEN_MIN, len_max=LEN_MAX),
          open("%s/pmhc_cfg.json" % W, "w"), indent=0)

WORKER = r'''
import json, os, sys, time
ARM = sys.argv[1]
os.environ["CUDA_VISIBLE_DEVICES"] = sys.argv[2]
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ["XLA_PYTHON_CLIENT_ALLOCATOR"] = "platform"
BUDGET = float(sys.argv[3])

import numpy as np
from Bio import PDB
from colabdesign import clear_mem, mk_afdesign_model

W = "/kaggle/working"
CFG = json.load(open("%s/pmhc_cfg.json" % W))
CF = CFG["arms"][ARM]
EXPECTED_A, PEP = CFG["expected_a"], CFG["peptide"]
MUT_POS, CONTACT = 6, 4.5
T3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
      "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
      "THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
t0 = time.time()
def log(m): print("[%s] %s [%.1f min]" % (ARM, m, (time.time()-t0)/60), flush=True)

def residues(ch):
    out, seen = [], set()
    for r in ch:
        if r.id[0] == " " and r.get_resname() in T3 and r.id[1] not in seen:
            seen.add(r.id[1]); out.append(r)
    return out
def seq(rs): return "".join(T3[r.get_resname()] for r in rs)
def heavy(r): return np.array([a.coord for a in r if a.element != "H"])

def measure(path):
    """save_pdb merges the target chains; find the merged one by sequence."""
    m = PDB.PDBParser(QUIET=True).get_structure("d", path)[0]
    chains = {c.id: residues(c) for c in m}
    tgt = pep_rs = None
    for cid, rs in chains.items():
        s = seq(rs)
        if s == EXPECTED_A + PEP:
            tgt, pep_rs = cid, rs[len(EXPECTED_A):]; break
        if s == PEP + EXPECTED_A:
            tgt, pep_rs = cid, rs[:len(PEP)]; break
    if tgt is None:
        raise RuntimeError("no chain reproduces the crop; chains are %s"
                           % {c: len(rs) for c, rs in chains.items()})
    if seq(pep_rs) != PEP:
        raise RuntimeError("located peptide reads %s" % seq(pep_rs))
    groove = [r for r in chains[tgt] if r not in pep_rs]
    other = [c for c in chains if c != tgt]
    if len(other) != 1:
        raise RuntimeError("expected one binder chain, got %s" % other)
    b = np.concatenate([heavy(r) for r in chains[other[0]]])
    def touch(rs):
        return sum(1 for r in rs
                   if np.linalg.norm(heavy(r)[:,None,:]-b[None,:,:],axis=-1).min() <= CONTACT)
    mut = pep_rs[MUT_POS-1]
    atoms = int(np.sum(np.linalg.norm(heavy(mut)[:,None,:]-b[None,:,:],axis=-1) <= CONTACT))
    return dict(on_mutation=touch([mut]), on_peptide=touch(pep_rs),
                on_groove=touch(groove), mutation_atoms=atoms)

rec_path = "%s/pmhc_%s.jsonl" % (W, ARM)
hot = ",".join("C%d" % p for p in CF["hotspot"])
i, made, failed = 0, 0, 0
while (time.time()-t0)/60.0 < BUDGET:
    seed = CF["seed"] + i
    rng = np.random.RandomState(seed)
    L = int(rng.randint(CFG["len_min"], CFG["len_max"]+1))
    name = "pmhc_%s_%03d" % (ARM, i); i += 1
    try:
        clear_mem()
        model = mk_afdesign_model(protocol="binder", data_dir=W)
        model.prep_inputs(pdb_filename="%s/target_design.pdb" % W, chain="A,C",
                          binder_len=L, hotspot=hot)
        before = dict(model.opt["i_con"])
        if "cutoff" not in before:
            raise SystemExit("no i_con cutoff in this colabdesign: %s" % before)
        model.set_opt(i_con=dict(cutoff=8.0, num=1, num_pos=CF["num_pos"]),
                      set_defaults=True)
        now = dict(model.opt["i_con"])
        if abs(float(now["cutoff"]) - 8.0) > 0.01:
            raise SystemExit("i_con cutoff still %s - refusing to repeat TfR1 round 1"
                             % now["cutoff"])
        if i == 1:
            log("i_con %s -> %s" % (before, now))
        model.restart(seed=seed)
        model.design_3stage(soft_iters=CFG["soft"], temp_iters=CFG["temp"],
                            hard_iters=CFG["hard"])
        out = "%s/%s.pdb" % (W, name)
        model.save_pdb(out)
        pae = model.aux.get("pae")
        if pae is not None:
            np.save("%s/%s_pae.npy" % (W, name), np.asarray(pae))
        met = measure(out)
        lg = model.aux.get("log", {})
        rec = dict(name=name, arm=ARM, seed=seed, length=L, hotspot=hot,
                   num_pos=CF["num_pos"], sequence=model.get_seq()[0],
                   pae_saved=pae is not None,
                   i_ptm=round(float(lg.get("i_ptm", 0)), 3),
                   plddt=round(float(lg.get("plddt", 0)), 3), **met)
        with open(rec_path, "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        made += 1
        log("%s len %d | ON-MUTATION %d (%d atoms) | peptide %d groove %d | "
            "i_ptm %.2f plddt %.2f | pae %s"
            % (name, L, rec["on_mutation"], rec["mutation_atoms"], rec["on_peptide"],
               rec["on_groove"], rec["i_ptm"], rec["plddt"],
               "yes" if rec["pae_saved"] else "NO"))
    except SystemExit:
        raise
    except Exception as e:
        failed += 1
        log("%s FAILED %s %s" % (name, type(e).__name__, str(e)[:160]))
        if failed >= 3 and made == 0:
            raise SystemExit("%s: three failures and nothing made - stopping" % ARM)
log("done: %d designs, %d failures" % (made, failed))
'''

open("%s/pmhc_worker.py" % W, "w").write(WORKER)

print("\nGPU state at launch (index, used, free):")
print(subprocess.run("nvidia-smi --query-gpu=index,memory.used,memory.free "
                     "--format=csv,noheader", shell=True, capture_output=True,
                     text=True).stdout.strip())
print("\nbudget %d min per worker, %d steps per design (~18 min each)\n"
      % (BUDGET_MIN, SOFT + TEMP + HARD))

procs = []
for arm, cfg in ARMS.items():
    p = subprocess.Popen([sys.executable, "%s/pmhc_worker.py" % W, arm,
                          str(cfg["gpu"]), str(BUDGET_MIN)],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         text=True, bufsize=1)
    procs.append((arm, p))
    print("launched arm %s on GPU %d (pid %d)" % (arm, cfg["gpu"], p.pid), flush=True)

for arm, p in procs:
    for line in p.stdout:
        print(line.rstrip(), flush=True)
    p.wait()

# ---------------------------------------------------------------- aggregate
records = []
for arm in ARMS:
    f = "%s/pmhc_%s.jsonl" % (W, arm)
    if os.path.exists(f):
        records += [json.loads(l) for l in open(f) if l.strip()]
json.dump(records, open("%s/pmhc_gen.json" % W, "w"), indent=0)

print("\n==================== ARM COMPARISON ====================")
print("THE QUESTION: does the tight set (5,6,7) reach residue 6 more often than the")
print("obvious set (4-9)? Asp6 is 11%% exposed; Gly7 beside it is 73%% and is")
print("identical in the wild type, so reaching the peptide is NOT the same as")
print("reaching the mutation.\n")
print("%-10s %5s %13s %11s %10s %10s" % ("arm", "n", "reached mut", "mean on-mut",
                                         "mean pept", "mut share"))
for arm in ARMS:
    rs = [r for r in records if r["arm"] == arm]
    if not rs:
        print("%-10s %5d  (nothing completed)" % (arm, 0))
        continue
    reached = sum(1 for r in rs if r["on_mutation"] > 0)
    share = [r["on_mutation"] / r["on_peptide"] for r in rs if r["on_peptide"] > 0]
    print("%-10s %5d %8d/%-4d %11.2f %10.2f %9s"
          % (arm, len(rs), reached, len(rs),
             np.mean([r["on_mutation"] for r in rs]),
             np.mean([r["on_peptide"] for r in rs]),
             "%.0f%%" % (100 * np.mean(share)) if share else "n/a"))
print("\nPAE saved for %d of %d designs." % (sum(1 for r in records if r["pae_saved"]),
                                             len(records)))
print("\n===== COPY FROM HERE =====")
for r in records:
    print(json.dumps(r))
print("===== COPY TO HERE =====")
