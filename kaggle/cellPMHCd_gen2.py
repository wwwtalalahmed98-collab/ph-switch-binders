"""Kaggle - Challenge 5 generation, run 2. Rebalanced so folding is not crowded out.

WHAT RUN 1 SHOWED. 28 designs, 4 of them touching Asp6 with a single residue each, and
pLDDT sitting around 0.40 - not folded proteins. The per-step lines said why:

    loss 4.59  i_con 4.53

The contact term was about 99% of the total loss, starting near 9.2 and plateauing near
4.5 without ever reaching satisfaction. Asp6 is 11% exposed at the bottom of a groove, so
that contact is geometrically hard; the optimiser spent the whole budget reaching for it
with nothing left to make the design a protein.

TWO THINGS CHANGE, AND THEY ARE TESTED SEPARATELY

  arm W   rebalanced weights, contact cutoff 8 A    GPU 0
  arm R   rebalanced weights, contact cutoff 10 A   GPU 1

Both arms get fold pressure back. They differ only in whether the contact being demanded
is geometrically reachable. If arm W still plateaus and arm R converges, the 8 A demand
was the problem; if both plateau, the method is wrong for this target and no amount of
tuning fixes it. That is a real possibility and the run is designed to show it.

THE WEIGHTS ARE PRINTED BEFORE THEY ARE OVERRIDDEN. I do not know ColabDesign's binder
defaults from memory and will not guess them into a 5-hour run. The worker prints what it
found, prints what it set, and refuses to continue if the override did not take - the
same guard that caught the i_con cutoff on TfR1.

THE SUCCESS CRITERION IS STATED IN ADVANCE. i_con must come down and keep coming down,
rather than stalling near 4.5. The worker records i_con at the first and last step of
every design so this is answerable from the records, not from reading the log by eye.

THE PIPE BUG FROM RUN 1 IS FIXED. The parent drained arm A's stream to completion before
touching arm B's, and ColabDesign logs a line per step - about 60 kB over 11 designs,
which fills Linux's 64 kB pipe. Arm B blocked on write for two hours and produced 11
designs against arm A's 17. Each worker now writes to its own file and the parent waits
on both, so neither can block the other.
"""
import glob
import json
import os
import subprocess
import sys
import time

W = "/kaggle/working"
BUDGET_MIN = 300
LEN_MIN, LEN_MAX = 70, 95
SOFT, TEMP, HARD = 35, 15, 5

PEP_G12D = "VVVGADGVGK"
DESIGN_RADIUS = 20.0
MIN_SEGMENT = 3
EXPECTED_A = ("SHSMRYFYTSVSRPGGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDQETRNVKAQSQ"
              "TDRVDLGTLRGYYNQSEDGSHTIQIMYGCDVGPDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQI"
              "TKRKWEAAHAAEQQRAYLEGRCVEWLRRYLENGKETLQR")

# Rebalanced. Run 1 left the contact term as essentially the whole loss; these put real
# weight back on "is this a protein" (plddt) and "is the interface confident" (i_pae),
# and on binder self-contact (con) so it stays compact rather than unspooling toward the
# groove. Printed against the defaults at run time so the change is visible, not implied.
WEIGHTS = dict(i_con=0.5, con=0.5, plddt=0.5, i_pae=0.4, pae=0.1)

ARMS = {"W_cut8": dict(cutoff=8.0, gpu=0, seed=8100),
        "R_cut10": dict(cutoff=10.0, gpu=1, seed=8300)}
HOTSPOT = [5, 6, 7]          # phase 1's set, adopted; num_pos = 3
NUM_POS = 3


def sh(c):
    subprocess.run(c, shell=True, check=True)


try:
    import colabdesign  # noqa: F401
    print("colabdesign already present")
except ImportError:
    sh("pip -q install git+https://github.com/sokrypton/ColabDesign.git@v1.1.1 biopython")

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


def residues(ch):
    out, seen = [], set()
    for r in ch:
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


build_target("%s/target_design.pdb" % W)
json.dump(dict(arms=ARMS, weights=WEIGHTS, hotspot=HOTSPOT, num_pos=NUM_POS,
               expected_a=EXPECTED_A, peptide=PEP_G12D, soft=SOFT, temp=TEMP,
               hard=HARD, len_min=LEN_MIN, len_max=LEN_MAX),
          open("%s/pmhc_cfg2.json" % W, "w"), indent=0)

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
CFG = json.load(open("%s/pmhc_cfg2.json" % W))
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

rec_path = "%s/pmhc2_%s.jsonl" % (W, ARM)
hot = ",".join("C%d" % p for p in CFG["hotspot"])
i, made, failed = 0, 0, 0
while (time.time()-t0)/60.0 < BUDGET:
    seed = CF["seed"] + i
    rng = np.random.RandomState(seed)
    L = int(rng.randint(CFG["len_min"], CFG["len_max"]+1))
    name = "pmhc2_%s_%03d" % (ARM, i); i += 1
    try:
        clear_mem()
        model = mk_afdesign_model(protocol="binder", data_dir=W)
        try:
            model.prep_inputs(pdb_filename="%s/target_design.pdb" % W, chain="A,C",
                              binder_len=L, hotspot=hot, rm_aa="C")
            rm_cys = True
        except TypeError:
            model.prep_inputs(pdb_filename="%s/target_design.pdb" % W, chain="A,C",
                              binder_len=L, hotspot=hot)
            rm_cys = False

        before_con = dict(model.opt["i_con"])
        if "cutoff" not in before_con:
            raise SystemExit("no i_con cutoff in this colabdesign: %s" % before_con)
        model.set_opt(i_con=dict(cutoff=CF["cutoff"], num=1, num_pos=CFG["num_pos"]),
                      set_defaults=True)
        now_con = dict(model.opt["i_con"])
        if abs(float(now_con["cutoff"]) - CF["cutoff"]) > 0.01:
            raise SystemExit("i_con cutoff is %s, wanted %s - refusing to run"
                             % (now_con["cutoff"], CF["cutoff"]))

        before_w = {k: float(v) for k, v in model.opt["weights"].items()}
        model.set_weights(**CFG["weights"])
        after_w = {k: float(v) for k, v in model.opt["weights"].items()}
        for k, v in CFG["weights"].items():
            if k not in after_w or abs(after_w[k] - float(v)) > 1e-6:
                raise SystemExit("weight %s did not take: wanted %s, got %s"
                                 % (k, v, after_w.get(k)))
        if i == 1:
            log("cysteines removed: %s" % rm_cys)
            log("i_con %s -> %s" % (before_con, now_con))
            log("weights BEFORE %s" % before_w)
            log("weights AFTER  %s" % after_w)

        model.restart(seed=seed)
        model.design_3stage(soft_iters=CFG["soft"], temp_iters=CFG["temp"],
                            hard_iters=CFG["hard"])

        # did the contact term actually converge, or stall as it did in run 1?
        hist = [h.get("i_con") for h in model._tmp.get("log", [])
                if isinstance(h, dict) and h.get("i_con") is not None]
        i_con_first = round(float(hist[0]), 3) if hist else None
        i_con_last = round(float(hist[-1]), 3) if hist else None
        i_con_min = round(float(min(hist)), 3) if hist else None

        out = "%s/%s.pdb" % (W, name)
        model.save_pdb(out)
        pae = model.aux.get("pae")
        if pae is not None:
            np.save("%s/%s_pae.npy" % (W, name), np.asarray(pae))
        met = measure(out)
        lg = model.aux.get("log", {})
        # _tmp["log"] is not guaranteed to exist; aux["log"]["i_con"] is. Record both,
        # because the convergence question is the whole point of this run and must not
        # depend on an attribute that may quietly be absent.
        i_con_final = (round(float(lg["i_con"]), 3) if lg.get("i_con") is not None
                       else None)
        rec = dict(name=name, arm=ARM, seed=seed, length=L, hotspot=hot,
                   i_con_final=i_con_final,
                   cutoff=CF["cutoff"], num_pos=CFG["num_pos"], rm_cys=rm_cys,
                   sequence=model.get_seq()[0], pae_saved=pae is not None,
                   i_con_first=i_con_first, i_con_last=i_con_last, i_con_min=i_con_min,
                   i_ptm=round(float(lg.get("i_ptm", 0)), 3),
                   plddt=round(float(lg.get("plddt", 0)), 3), **met)
        with open(rec_path, "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        made += 1
        log("%s len %d | ON-MUT %d (%d atoms) | pept %d groove %d | i_ptm %.2f "
            "plddt %.2f | i_con %s->%s final %s"
            % (name, L, rec["on_mutation"], rec["mutation_atoms"], rec["on_peptide"],
               rec["on_groove"], rec["i_ptm"], rec["plddt"], i_con_first, i_con_last,
               i_con_final))
    except SystemExit:
        raise
    except Exception as e:
        failed += 1
        log("%s FAILED %s %s" % (name, type(e).__name__, str(e)[:160]))
        if failed >= 3 and made == 0:
            raise SystemExit("%s: three failures and nothing made - stopping" % ARM)
log("done: %d designs, %d failures" % (made, failed))
'''

open("%s/pmhc_worker2.py" % W, "w").write(WORKER)

print("\nGPU state at launch (index, used, free):")
print(subprocess.run("nvidia-smi --query-gpu=index,memory.used,memory.free "
                     "--format=csv,noheader", shell=True, capture_output=True,
                     text=True).stdout.strip())
print("\nrebalanced weights to apply: %s" % WEIGHTS)
print("hotspots %s, num_pos %d, %d steps per design\n"
      % (HOTSPOT, NUM_POS, SOFT + TEMP + HARD))

# Each worker writes to its OWN FILE. In run 1 the parent drained one pipe to completion
# before touching the other, and ColabDesign's per-step logging filled the 64 kB buffer,
# so the second worker blocked on write for two hours. Files cannot block.
procs = []
for arm, cfg in ARMS.items():
    fh = open("%s/log_%s.txt" % (W, arm), "w")
    p = subprocess.Popen([sys.executable, "%s/pmhc_worker2.py" % W, arm,
                          str(cfg["gpu"]), str(BUDGET_MIN)],
                         stdout=fh, stderr=subprocess.STDOUT, text=True)
    procs.append((arm, p, fh))
    print("launched arm %s on GPU %d (pid %d) -> log_%s.txt"
          % (arm, cfg["gpu"], p.pid, arm), flush=True)

while any(p.poll() is None for _, p, _ in procs):
    time.sleep(120)
    for arm, _, _ in procs:
        f = "%s/pmhc2_%s.jsonl" % (W, arm)
        n = sum(1 for _ in open(f)) if os.path.exists(f) else 0
        print("  %s: %d designs so far" % (arm, n), flush=True)
    sys.stdout.flush()
for _, p, fh in procs:
    p.wait()
    fh.close()

for arm in ARMS:
    print("\n---------- log for arm %s ----------" % arm)
    f = "%s/log_%s.txt" % (W, arm)
    if os.path.exists(f):
        for line in open(f):
            if line.startswith("[") or "FAILED" in line or "done:" in line:
                print(line.rstrip())

records = []
for arm in ARMS:
    f = "%s/pmhc2_%s.jsonl" % (W, arm)
    if os.path.exists(f):
        records += [json.loads(l) for l in open(f) if l.strip()]
json.dump(records, open("%s/pmhc_gen2.json" % W, "w"), indent=0)

print("\n==================== REBALANCE RESULT ====================")
print("RUN 1 BASELINE: 4 of 28 touched the mutation, mean pLDDT about 0.40, and the")
print("contact term stalled near 4.5 instead of converging.\n")
if not records:
    raise SystemExit("nothing was produced")
print("%-9s %4s %12s %9s %9s %11s %11s"
      % ("arm", "n", "reached mut", "mean mut", "mean pept", "mean plddt", "i_con last"))
for arm in ARMS:
    rs = [r for r in records if r["arm"] == arm]
    if not rs:
        print("%-9s %4d  (nothing completed)" % (arm, 0))
        continue
    lasts = [r.get("i_con_last") if r.get("i_con_last") is not None
             else r.get("i_con_final") for r in rs]
    lasts = [x for x in lasts if x is not None]
    print("%-9s %4d %7d/%-4d %9.2f %9.2f %11.2f %11s"
          % (arm, len(rs), sum(1 for r in rs if r["on_mutation"] > 0), len(rs),
             np.mean([r["on_mutation"] for r in rs]),
             np.mean([r["on_peptide"] for r in rs]),
             np.mean([r["plddt"] for r in rs]),
             "%.2f" % np.mean(lasts) if lasts else "n/a"))
print("\nTHE QUESTION THE WEIGHTS WERE CHANGED TO ANSWER: is mean pLDDT above run 1's")
print("0.40? If not, fold pressure was not the binding constraint.")
print("THE QUESTION THE CUTOFF ARM ANSWERS: does i_con still stall near 4.5 at 10 A?")
print("If it stalls at both cutoffs, the method is wrong for this target.\n")
print("PAE saved for %d of %d." % (sum(1 for r in records if r["pae_saved"]),
                                   len(records)))
print("\n===== COPY FROM HERE =====")
for r in records:
    print(json.dumps(r))
print("===== COPY TO HERE =====")
