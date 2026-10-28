"""Kaggle - TfR1 Challenge 4. A one-hour probe: why did the designs bind the cut face?

RUN THIS ALONE. It builds its own targets and needs no attached input. One hour of GPU time
out of the 6 h 50 m left this week; the remaining ~5 h goes to whichever arm wins.

THE PROBLEM IT IS TESTING
-------------------------
Generation run 1 produced 47 designs, of which 32 made ZERO contacts with the chosen
epitope, including all five of the most confident ones (i_ptm 0.66 to 0.78, 13 to 31
contacts elsewhere). The crop exposes 3,427 A^2 of interior surface, 1,722 A^2 of it
hydrophobic, and the designs went there instead. The hotspot argument biases the optimiser;
it does not forbid binding elsewhere, and on this target the bias lost.

There are two candidate fixes and they are not the same fix:

  ARM A, "pull harder" - same 20 A crop, but the interface-contact term weighted up. When a
    hotspot is given, AfDesign scores i_con against those target positions only, so raising
    its weight should pay the optimiser specifically for epitope contacts. That is a claim
    about the library, so the weights actually in effect are PRINTED, and the run refuses if
    the key is not there rather than quietly training with defaults.

  ARM B, "offer less cut face" - default weights, but a TIGHTER crop (18 A). The geometry
    was swept on the local structure before choosing this, and the sweep overturned the
    obvious guess. The epitope's entire exposed surface in the intact protein is 1012 A^2;
    run 1's 20 A crop presented 1722 A^2 of oily cut face, a competitor 1.7 times the size
    of the target itself. Widening the crop makes that WORSE, not better - 22 A gives 2.4x,
    26 A gives 2.6x - because a bigger crop has a bigger rim. Tightening it to 18 A gives
    1.3x, puts only 30 A^2 of new surface on the epitope itself, and costs 0.68 of run 1's
    time per design, so this arm also returns MORE designs, not fewer.

    What tightening gives up is distance: the epitope sits 3.8 A from the nearest cut edge
    instead of 6.5 A. So the two arms together also separate the two mechanisms - if the cut
    face wins by sheer area, the tight arm improves; if it wins by being next door, the
    tight arm gets worse. Either result is informative.

Each arm changes exactly ONE thing relative to run 1, so the result is attributable.

THE READOUT is on-site contacts, measured here in-process with the corrected chain matching
(the target is found by reproducing the crop residue for residue, not by its numbering -
save_pdb renumbers from 1, which is what left run 1 unmeasured). Run 1's baseline, to beat:

    at least 4 epitope contacts:  3 of 47  (6%)
    at least 2 epitope contacts: 15 of 47  (32%)
    mean epitope contacts:       0.96

NOT i_ptm. Confidence versus epitope contacts measured r = +0.103 on this target, so it
carries no information about where a binder sits and is recorded only for the record.

NO GATE CELL. Run 1 already proved AlphaFold runs on these targets, and the parent process
here never touches the GPU, so the preallocation hazard that needed the gate does not exist.
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

W = "/kaggle/working"
PROBE_MIN = 55                  # wall clock per worker; the two run in parallel
LEN_MIN, LEN_MAX = 55, 75
I_CON_WEIGHT = 3.0              # arm A only
SITE = "e559"                   # one site, so the two arms are comparable

# Run 1 (radius 20, segments >=3, default weights) is the shared control. Each arm changes
# exactly one thing relative to it, so a difference is attributable.
ARMS = {
    "pull":  dict(radius=20.0, min_seg=3, i_con=I_CON_WEIGHT),
    "tight": dict(radius=18.0, min_seg=3, i_con=None),
}
EPITOPE = [287, 289, 486, 491, 492, 493, 495, 517, 559, 560, 561, 562, 579]


def sh(c):
    print("$", c, flush=True)
    subprocess.run(c, shell=True, check=True)


try:
    import colabdesign  # noqa: F401
    print("colabdesign already present")
except ImportError:
    sh("pip -q install git+https://github.com/sokrypton/ColabDesign.git@v1.1.1")

if not os.path.exists("%s/params/params_model_1_ptm.npz" % W):
    os.makedirs("%s/params" % W, exist_ok=True)
    src = glob.glob("/kaggle/input/**/params/params_model_*.npz", recursive=True)
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
EP_AREA = sum(FULL[n] for n in EPITOPE if n in FULL)
print("epitope exposed area in the intact protein: %.0f A^2" % EP_AREA)


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
centre = np.mean([POS[n] for n in EPITOPE], axis=0)
for arm, cfg in sorted(ARMS.items()):
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
          % (arm, cfg["radius"], len(keep), new, phobic, phobic / EP_AREA, EP_AREA, ep_new))
    if ep_new > 150:
        raise SystemExit("%s: the crop changed the epitope itself" % arm)
    MAPS[arm] = dict(arm=arm, target=out, chain=REF, epitope=EPITOPE, kept=keep,
                     hotspot=",".join("%s%d" % (REF, n) for n in EPITOPE),
                     i_con=cfg["i_con"], new_surface_A2=round(new, 1),
                     new_hydrophobic_A2=round(phobic, 1))

json.dump(MAPS, open("%s/tfr1_probe_targets.json" % W, "w"), indent=0)
print("\nBaseline to beat, from run 1: 3 of 47 reach 4+ epitope contacts, mean 0.96\n")

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

W = "/kaggle/working"
M = json.load(open("%s/tfr1_probe_targets.json" % W))[ARM]
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

rec_path = "%s/probe_%s.jsonl" % (W, ARM)
rng = random.Random(4100 if ARM == "pull" else 4300)   # different seeds per arm
made = i = 0

while (time.time() - t0) / 60 < BUDGET:
    i += 1
    name = "probe_%s_%03d" % (ARM, i)
    L = rng.randint(LEN_MIN, LEN_MAX)
    try:
        clear_mem()
        model = mk_afdesign_model(protocol="binder", data_dir=W)
        kw = dict(pdb_filename=TARGET, chain=CHAIN, binder_len=L, hotspot=M["hotspot"])
        try:
            model.prep_inputs(rm_aa="C", **kw)
        except TypeError:
            model.prep_inputs(**kw)

        if M["i_con"] is not None:
            before = dict(model.opt["weights"])
            if "i_con" not in before:
                raise SystemExit("this colabdesign has no i_con weight; keys are %s - the"
                                 " arm cannot be run as described" % sorted(before))
            model.set_weights(i_con=M["i_con"])
            if i == 1:
                log("weights now in effect: %s (were %s)"
                    % (dict(model.opt["weights"]), before))
        elif i == 1:
            log("default weights: %s" % dict(model.opt["weights"]))

        # half of run 1's iteration counts, to buy designs instead of polish - the readout is
        # WHERE the binder sits, which is decided early, not how refined it is
        model.design_3stage(soft_iters=30, temp_iters=15, hard_iters=5)
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

with open("%s/tfr1_probe_worker.py" % W, "w") as fh:
    fh.write(textwrap.dedent(WORKER))

print("GPU state at launch (index, used, free):")
print(subprocess.run("nvidia-smi --query-gpu=index,memory.used,memory.free "
                     "--format=csv,noheader", shell=True,
                     capture_output=True, text=True).stdout.strip(), flush=True)

procs = []
for gpu, arm in enumerate(sorted(ARMS)):
    p = subprocess.Popen([sys.executable, "%s/tfr1_probe_worker.py" % W, arm, str(gpu),
                          str(PROBE_MIN), str(LEN_MIN), str(LEN_MAX)],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                         bufsize=1)
    procs.append(p)
    print("launched arm %s on GPU %d (pid %d)" % (arm, gpu, p.pid), flush=True)

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
for f in sorted(glob.glob("%s/probe_*.jsonl" % W)):
    for line in open(f):
        line = line.strip()
        if line:
            recs.append(json.loads(line))
json.dump(recs, open("%s/tfr1_probe.json" % W, "w"), indent=0)

print("\n==================== PROBE RESULT ====================")
print("run 1 baseline: 3 of 47 reach 4+ epitope contacts (6%), mean 0.96\n")
print("%-6s %5s %7s %7s %7s %8s" % ("arm", "n", "mean on", ">=2", ">=4", "mean off"))
for arm in sorted(ARMS):
    rs = [r for r in recs if r["arm"] == arm and "on_site" in r]
    if not rs:
        print("%-6s %5d  no measured designs" % (arm, 0))
        continue
    on = [r["on_site"] for r in rs]
    print("%-6s %5d %7.2f %7d %7d %8.1f"
          % (arm, len(rs), sum(on) / float(len(on)), sum(1 for x in on if x >= 2),
             sum(1 for x in on if x >= 4),
             sum(r["off_site"] for r in rs) / float(len(rs))))
print("\nn is small; read a flat zero as 'no signal', not as 'slightly worse'.")
print("Judge on mean on-site and the >=4 count, NOT on i_ptm (r = +0.103 here).")

print("\n===== COPY FROM HERE =====")
for r in recs:
    print(json.dumps(r))
print("===== COPY TO HERE =====")
