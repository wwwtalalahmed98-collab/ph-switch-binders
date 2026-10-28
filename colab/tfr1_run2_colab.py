"""Colab - TfR1 Challenge 4, generation run 2 on a second platform.

WHY THIS EXISTS. Kaggle gives 30 GPU-hours a week and we are nearly out with three days
left. Colab has its own separate free allowance, so running the SAME generation here adds
designs without touching Kaggle's quota. The measurement code is byte-identical to
kaggle/cellTFR1f_run2.py so results from the two platforms are directly comparable; only the
paths, the worker count and the random seeds differ.

WHAT IS DIFFERENT FROM THE KAGGLE VERSION, and it is only this:
  - writes to /content instead of /kaggle/working
  - ONE worker, because free Colab gives one GPU, not two
  - SITE picks which epitope to work on; Kaggle is running both, so this adds depth to one
  - seeds 6100/6300 instead of 5100/5300, so designs do not duplicate Kaggle's
  - a shorter budget, because free Colab sessions are dropped without warning

The correction this run exists to test is unchanged:
    set_opt(i_con=dict(cutoff=8.0, num=1, num_pos=6), set_defaults=True)
ColabDesign's interface contact term defaults to a 21.6875 A cutoff, which passes 63 of
AlphaFold's 64 distance bins and counts anything under 21.4 A as a satisfied contact. The
cut face of our cropped target sits 6.5 A from the epitope, so binders gripping the wrong
surface already satisfied the objective. See tfr1/results/PROBE_RESULT.md.

DISCONNECTS. Every design is appended to a .jsonl AND printed as it finishes, so a dropped
session loses only the designs still in flight - whatever has printed is recoverable from
the cell output. Keep the tab open and the machine awake.

UNTESTED ON COLAB. Colab ships its own jax, and ColabDesign v1.1.1 may disagree with it.
If the import fails on a jax version, that is the first thing to look at - not the design
logic, which has run on Kaggle.
"""
import os

os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ["XLA_PYTHON_CLIENT_ALLOCATOR"] = "platform"

import glob
import json
import subprocess
import sys
import textwrap
import urllib.request

W = "/content"
BUDGET_MIN = 150                # free Colab sessions get dropped; a shorter budget means
                                # the summary actually prints. Raise it if the session holds.
SITE = "e559"                   # which epitope this platform works on. Kaggle runs both.
LEN_MIN, LEN_MAX = 55, 80
I_CON = dict(cutoff=8.0, num=1, num_pos=6)   # the correction; see the header

# Run 1's crop, unchanged. One worker per epitope, as in run 1.
ARMS = {
    "e559": dict(radius=20.0, min_seg=3),
    "l537": dict(radius=20.0, min_seg=3),
}
EPITOPES = {
    "e559": [287, 289, 486, 491, 492, 493, 495, 517, 559, 560, 561, 562, 579],
    "l537": [261, 262, 266, 269, 303, 533, 534, 536, 537, 538, 695, 696, 704, 705],
}


def sh(c):
    print("$", c, flush=True)
    subprocess.run(c, shell=True, check=True)


try:
    import colabdesign  # noqa: F401
    print("colabdesign already present")
except ImportError:
    sh("pip -q install git+https://github.com/sokrypton/ColabDesign.git@v1.1.1")

# Colab does not ship biopython. The Kaggle image did only because ColabDesign pulled it in
# there; the CPU remeasure notebook hit exactly this and failed on "No module named 'Bio'".
try:
    import Bio  # noqa: F401
except ImportError:
    sh("pip -q install biopython")

if not os.path.exists("%s/params/params_model_1_ptm.npz" % W):
    os.makedirs("%s/params" % W, exist_ok=True)
    src = glob.glob("/content/drive/**/params/params_model_*.npz", recursive=True)
    if src:
        import shutil
        for f in src:
            shutil.copy(f, "%s/params/%s" % (W, os.path.basename(f)))
    else:
        sh("wget -q -O %s/af2.tar https://storage.googleapis.com/alphafold/"
           "alphafold_params_2022-12-06.tar" % W)
        sh("tar -xf %s/af2.tar -C %s/params" % (W, W))
        os.remove("%s/af2.tar" % W)
print("params:", len(glob.glob("%s/params/params_model_*.npz" % W)), "files")

CIF = "%s/11QC.cif" % W
FA = "%s/P02786.fasta" % W
if not os.path.exists(CIF):
    urllib.request.urlretrieve("https://files.rcsb.org/download/11QC.cif", CIF)
if not os.path.exists(FA):
    urllib.request.urlretrieve("https://rest.uniprot.org/uniprotkb/P02786.fasta", FA)
HUMAN = "".join(l.strip() for l in open(FA) if not l.startswith(">"))

import numpy as np                                    # noqa: E402
from Bio import PDB                                   # noqa: E402
from Bio.PDB.SASA import ShrakeRupley                 # noqa: E402

T31 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
       "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
       "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}

st = PDB.MMCIFParser(QUIET=True).get_structure("x", CIF)[0]
chains = [c for c in st if sum(1 for r in c if r.id[0] == " ") > 100]
REF = chains[0].id
bad = sum(1 for c in chains for r in c
          if r.id[0] == " " and r.get_resname() in T31
          and 1 <= r.id[1] <= len(HUMAN) and T31[r.get_resname()] != HUMAN[r.id[1] - 1])
print("numbering check against P02786: %d mismatches" % bad)
if bad:
    raise SystemExit("PDB numbering disagrees with UniProt - stop and fix before designing")

ShrakeRupley().compute(st, level="R")
FULL = {r.id[1]: r.sasa for r in st[REF] if r.id[0] == " " and r.get_resname() in T31}
POS = {r.id[1]: r["CA"].get_coord() for r in st[REF]
       if r.id[0] == " " and "CA" in r and r.get_resname() in T31}
HYDRO = set("AVLIMFWYC")
# the epitope's own exposed area in the intact protein - the yardstick the cut face is
# compared against, because "3427 A^2 of new surface" means nothing without it
EP_AREA = {}   # per epitope; the cut face is only meaningful against its own target


class Sel(PDB.Select):
    def __init__(self, keep):
        self.keep = set(keep)

    def accept_model(self, m):
        return m.id == 0

    def accept_chain(self, c):
        return c.id == REF

    def accept_residue(self, r):
        return r.id[0] == " " and r.id[1] in self.keep

    def accept_atom(self, a):
        return (not a.is_disordered()) or a.get_altloc() in (" ", "A")


MAPS = {}
for arm, cfg in sorted(ARMS.items()):
    EPITOPE = EPITOPES[arm]
    EP_AREA[arm] = sum(FULL[n] for n in EPITOPE if n in FULL)
    centre = np.mean([POS[n] for n in EPITOPE], axis=0)
    keep = sorted(n for n, x in POS.items() if np.linalg.norm(x - centre) <= cfg["radius"])
    runs, cur = [], [keep[0]]
    for a, b in zip(keep, keep[1:]):
        if b - a <= 2:
            cur.append(b)
        else:
            runs.append(cur)
            cur = [b]          # a NEW list; append-then-clear stores a reference and
    runs.append(cur)           # then empties what was stored
    keep = sorted(n for r_ in runs if len(r_) >= cfg["min_seg"] for n in r_)
    lost = [n for n in EPITOPE if n not in keep]
    if lost:
        raise SystemExit("%s: crop lost epitope residues %s" % (arm, lost))

    out = "%s/TfR1_%s.pdb" % (W, arm)
    io_ = PDB.PDBIO()
    io_.set_structure(st)
    io_.save(out, Sel(keep))

    cr = PDB.PDBParser(QUIET=True).get_structure("c", out)
    ShrakeRupley().compute(cr, level="R")
    new = ep_new = phobic = 0.0
    for r in cr[0][REF]:
        d = r.sasa - FULL.get(r.id[1], r.sasa)
        if d > 1.0:
            new += d
            if T31.get(r.get_resname()) in HYDRO:
                phobic += d
            if r.id[1] in EPITOPE:
                ep_new += d
    print("%-6s radius %.0f A: %d residues, cut created %.0f A^2 (%.0f A^2 hydrophobic "
          "= %.1fx the epitope's own %.0f A^2), %.0f A^2 of it on the epitope"
          % (arm, cfg["radius"], len(keep), new, phobic, phobic / EP_AREA[arm], EP_AREA[arm], ep_new))
    if ep_new > 150:
        raise SystemExit("%s: the crop changed the epitope itself" % arm)
    MAPS[arm] = dict(arm=arm, target=out, chain=REF, epitope=EPITOPE, kept=keep,
                     hotspot=",".join("%s%d" % (REF, n) for n in EPITOPE),
                     i_con=I_CON, new_surface_A2=round(new, 1),
                     new_hydrophobic_A2=round(phobic, 1))

json.dump(MAPS, open("%s/tfr1_run2_targets.json" % W, "w"), indent=0)
print("\nBaseline to beat, from run 1: 3 of 47 reach 4+ epitope contacts, mean 0.96")
print("i_con correction in effect: %s\n" % I_CON)

WORKER = r'''
import json, os, random, sys, time

ARM = sys.argv[1]; GPU = sys.argv[2]; BUDGET = float(sys.argv[3])
LEN_MIN = int(sys.argv[4]); LEN_MAX = int(sys.argv[5])
os.environ["CUDA_VISIBLE_DEVICES"] = GPU            # must precede the jax import
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ["XLA_PYTHON_CLIENT_ALLOCATOR"] = "platform"
t0 = time.time()
def log(m): print("[%s] %s [%.1f min]" % (ARM, m, (time.time()-t0)/60), flush=True)

import numpy as np
from Bio.PDB import PDBParser
from colabdesign import mk_afdesign_model, clear_mem

W = "/content"
M = json.load(open("%s/tfr1_run2_targets.json" % W))[ARM]
TARGET, CHAIN, KEPT = M["target"], M["chain"], M["kept"]
EPITOPE = set(M["epitope"])

def dedupe(chain):
    out, seen = [], set()
    for r in chain:
        if r.id[0] == " " and r.id[1] not in seen:
            seen.add(r.id[1]); out.append(r)
    return out

# The crop's residue types in order. save_pdb renumbers from 1, so the target chain is found
# by reproducing this list - not by looking for epitope residue numbers, which is exactly
# what left all 47 of run 1's designs unmeasured.
CROP_SEQ = [r.get_resname() for r in dedupe(list(PDBParser(QUIET=True).get_structure("t", TARGET)[0])[0])]
if len(CROP_SEQ) != len(KEPT):
    raise SystemExit("crop PDB %d residues, kept list %d" % (len(CROP_SEQ), len(KEPT)))

def longest_run(s):
    best = cur = 1
    for i in range(1, len(s)):
        cur = cur + 1 if s[i] == s[i-1] else 1
        best = max(best, cur)
    return best

KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5,
          L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)

def measure(pdb_path, blen):
    mdl = PDBParser(QUIET=True).get_structure("c", pdb_path)[0]
    kept = {ch.id: dedupe(ch) for ch in mdl}
    tid = None
    for cid, rs in kept.items():
        if [r.get_resname() for r in rs] == CROP_SEQ:
            tid = cid; break
    if tid is None:
        raise ValueError("no chain reproduces the crop; %s" % {k: len(v) for k, v in kept.items()})
    bid = [c for c in kept if c != tid]
    if len(bid) != 1:
        raise ValueError("expected one binder chain, got %s" % bid)
    bres = kept[bid[0]]
    ba = np.array([a.get_coord() for r in bres for a in r])
    on = off = 0; touched = []
    for i, r in enumerate(kept[tid]):
        ra = np.array([a.get_coord() for a in r])
        if np.min(np.linalg.norm(ra[:, None, :] - ba[None, :, :], axis=-1)) <= 4.5:
            if KEPT[i] in EPITOPE:
                on += 1; touched.append(KEPT[i])
            else:
                off += 1
    return dict(on_site=on, off_site=off, touched=sorted(touched),
                on_fraction=round(on / float(max(on + off, 1)), 3))

rec_path = "%s/run2_%s.jsonl" % (W, ARM)
rng = random.Random(6100 if ARM == "e559" else 6300)   # different seeds per arm
made = i = 0

while (time.time() - t0) / 60 < BUDGET:
    i += 1
    name = "tfr1c_%s_%03d" % (ARM, i)
    L = rng.randint(LEN_MIN, LEN_MAX)
    try:
        clear_mem()
        model = mk_afdesign_model(protocol="binder", data_dir=W)
        kw = dict(pdb_filename=TARGET, chain=CHAIN, binder_len=L, hotspot=M["hotspot"])
        try:
            model.prep_inputs(rm_aa="C", **kw)
        except TypeError:
            model.prep_inputs(**kw)

        before = dict(model.opt["i_con"])
        if "cutoff" not in before:
            raise SystemExit("this colabdesign has no i_con cutoff; i_con is %s - the"
                             " correction cannot be applied" % before)
        model.set_opt(i_con=dict(M["i_con"]), set_defaults=True)
        now = dict(model.opt["i_con"])
        if abs(float(now["cutoff"]) - 8.0) > 0.01:
            raise SystemExit("i_con cutoff is still %s - refusing to spend the quota"
                             " repeating run 1" % now["cutoff"])
        if i == 1:
            log("i_con now %s (was %s)" % (now, before))
            log("weights %s" % dict(model.opt["weights"]))

        model.design_3stage(soft_iters=60, temp_iters=30, hard_iters=10)
        seq = model.get_seqs()[0]
        try:
            model.predict(seq=seq, num_recycles=3, verbose=False)
        except TypeError:
            model.predict(seq=seq, verbose=False)
        lg = model.aux["log"]
        out = "%s/%s.pdb" % (W, name)
        model.save_pdb(out)
        r = dict(name=name, arm=ARM, length=L, sequence=seq,
                 true_iptm=round(float(lg.get("i_ptm", lg.get("ptm", 0.0))), 3),
                 plddt=round(float(lg.get("plddt", 0.0)), 3),
                 hydropathy=round(sum(KD[c] for c in seq) / len(seq), 3),
                 n_cys=seq.count("C"), longest_run=longest_run(seq))
        try:
            r.update(measure(out, L))
        except Exception as e:
            r["measure_error"] = "%s %s" % (type(e).__name__, str(e)[:160])
        with open(rec_path, "a") as fh:
            fh.write(json.dumps(r) + "\n")
        made += 1
        log("%s len=%d | ON-SITE %d/%d, off-site %d | i_ptm %.3f pLDDT %.3f run %d"
            % (name, L, r.get("on_site", -1), len(EPITOPE), r.get("off_site", -1),
               r["true_iptm"], r["plddt"], r["longest_run"]))
    except SystemExit:
        raise
    except Exception as e:
        log("%s FAILED %s %s" % (name, type(e).__name__, str(e)[:200]))

log("done: %d designs from %d attempts" % (made, i))
if made == 0:
    raise SystemExit("arm %s produced nothing - read the failures above" % ARM)
'''

with open("%s/tfr1_run2_worker.py" % W, "w") as fh:
    fh.write(textwrap.dedent(WORKER))

print("GPU state at launch (index, used, free):")
print(subprocess.run("nvidia-smi --query-gpu=index,memory.used,memory.free "
                     "--format=csv,noheader", shell=True,
                     capture_output=True, text=True).stdout.strip(), flush=True)

if SITE not in ARMS:
    raise SystemExit("SITE is %r but the targets built are %s" % (SITE, sorted(ARMS)))

procs = []
p = subprocess.Popen([sys.executable, "%s/tfr1_run2_worker.py" % W, SITE, "0",
                      str(BUDGET_MIN), str(LEN_MIN), str(LEN_MAX)],
                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                     bufsize=1)
procs.append(p)
print("launched site %s on the single GPU (pid %d)" % (SITE, p.pid), flush=True)

import threading


def pump(p):
    for line in p.stdout:
        sys.stdout.write(line)
        sys.stdout.flush()


threads = [threading.Thread(target=pump, args=(p,), daemon=True) for p in procs]
for t in threads:
    t.start()
for p in procs:
    p.wait()
for t in threads:
    t.join(timeout=10)

recs = []
for f in sorted(glob.glob("%s/run2_*.jsonl" % W)):
    for line in open(f):
        line = line.strip()
        if line:
            recs.append(json.loads(line))
json.dump(recs, open("%s/tfr1_designs2_colab.json" % W, "w"), indent=0)

print("\n============== RUN 2 RESULT (COLAB) ==============")
print("run 1 baseline: 3 of 47 reach 4+ epitope contacts (6%), mean 0.96\n")
print("%-6s %5s %7s %7s %7s %8s" % ("arm", "n", "mean on", ">=2", ">=4", "mean off"))
for arm in sorted(set(r["arm"] for r in recs)):
    rs = [r for r in recs if r["arm"] == arm and "on_site" in r]
    if not rs:
        print("%-6s %5d  no measured designs" % (arm, 0))
        continue
    on = [r["on_site"] for r in rs]
    print("%-6s %5d %7.2f %7d %7d %8.1f"
          % (arm, len(rs), sum(on) / float(len(on)), sum(1 for x in on if x >= 2),
             sum(1 for x in on if x >= 4),
             sum(r["off_site"] for r in rs) / float(len(rs))))
print("\nJudge on mean on-site and the >=4 count, NOT on i_ptm (r = +0.103 on this target).")

print("\n===== COPY FROM HERE =====")
for r in recs:
    print(json.dumps(r))
print("===== COPY TO HERE =====")
