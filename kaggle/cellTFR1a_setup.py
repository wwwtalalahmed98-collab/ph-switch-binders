"""Kaggle cell 1 of 2 - TfR1 Challenge 4. Install, build both targets, prove they are right.

Nothing is designed here. This cell exists so that the eight-hour run in cell 2 starts from
a target that has been checked rather than assumed. Three separate failures in this project
came from a target that was wrong in a way nothing printed: a chain renumbered by ordinal,
two numbering systems that silently disagreed, and a crop that invented a binding cavity.

What it proves before cell 2 is allowed to start:
  * the PDB's residue numbering matches UniProt P02786, residue by residue
  * every epitope residue survived the crop
  * the epitope looks the same cropped as it does in the intact dimer
  * AlphaFold actually runs, on a two-minute design, not on faith

GPU MEMORY. The gate below runs AlphaFold in THIS process, and this process does not exit -
it goes on to launch the workers in cell 2. Without the two environment variables set before
jax is imported, JAX preallocates most of GPU 0 and holds it, and every worker on GPU 0 then
dies out of memory. That cost this project twenty-odd failed runs on the previous challenge
before the cause was found, so the allocator settings and the nvidia-smi print below are not
optional decoration.
"""
import os

os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ["XLA_PYTHON_CLIENT_ALLOCATOR"] = "platform"

import glob
import json
import re
import shutil
import subprocess
import sys
import urllib.request

W = "/kaggle/working"


def sh(c):
    print("$", c, flush=True)
    subprocess.run(c, shell=True, check=True)


# ---------------------------------------------------------------- dependencies
try:
    import colabdesign  # noqa: F401
    print("colabdesign already present")
except ImportError:
    sh("pip -q install git+https://github.com/sokrypton/ColabDesign.git@v1.1.1")

if not os.path.exists("%s/params/params_model_1_ptm.npz" % W):
    src = glob.glob("/kaggle/input/**/params/params_model_*.npz", recursive=True)
    os.makedirs("%s/params" % W, exist_ok=True)
    if src:
        print("using AF2 params from an attached dataset (%d files)" % len(src))
        for f in src:
            shutil.copy(f, "%s/params/%s" % (W, os.path.basename(f)))
    else:
        sh("wget -q -O %s/af2.tar https://storage.googleapis.com/alphafold/"
           "alphafold_params_2022-12-06.tar" % W)
        sh("tar -xf %s/af2.tar -C %s/params" % (W, W))
        os.remove("%s/af2.tar" % W)
print("params:", len(glob.glob("%s/params/params_model_*.npz" % W)), "files")

# ---------------------------------------------------------------- inputs
os.makedirs("%s/targets" % W, exist_ok=True)
CIF = "%s/11QC.cif" % W
FA = "%s/P02786.fasta" % W
if not os.path.exists(CIF):
    urllib.request.urlretrieve("https://files.rcsb.org/download/11QC.cif", CIF)
if not os.path.exists(FA):
    urllib.request.urlretrieve("https://rest.uniprot.org/uniprotkb/P02786.fasta", FA)

HUMAN = "".join(l.strip() for l in open(FA) if not l.startswith(">"))
print("human TfR1 %d aa" % len(HUMAN))

# Epitopes chosen in tfr1/results/EPITOPE_CHOICE.md after conservation, glycan and
# transferrin-site filtering. UniProt numbering.
SITES = {
    "e559": [287, 289, 486, 491, 492, 493, 495, 517, 559, 560, 561, 562, 579],
    "l537": [261, 262, 266, 269, 303, 533, 534, 536, 537, 538, 695, 696, 704, 705],
}
RADIUS = 20.0

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
for tag, epitope in SITES.items():
    centre = np.mean([POS[n] for n in epitope], axis=0)
    keep = sorted(n for n, x in POS.items() if np.linalg.norm(x - centre) <= RADIUS)
    runs, cur = [], [keep[0]]
    for a, b in zip(keep, keep[1:]):
        if b - a <= 2:
            cur.append(b)
        else:
            runs.append(cur)
            cur = [b]          # a NEW list: append-then-clear stores a reference and
    runs.append(cur)           # then empties what was stored
    keep = sorted(n for r_ in runs if len(r_) >= 3 for n in r_)
    lost = [n for n in epitope if n not in keep]
    if lost:
        raise SystemExit("%s crop lost epitope residues %s" % (tag, lost))

    out = "%s/targets/TfR1_%s.pdb" % (W, tag)
    io_ = PDB.PDBIO()
    io_.set_structure(st)
    io_.save(out, Sel(keep))

    cr = PDB.PDBParser(QUIET=True).get_structure("c", out)
    ShrakeRupley().compute(cr, level="R")
    new = ep_new = 0.0
    for r in cr[0][REF]:
        d = r.sasa - FULL.get(r.id[1], r.sasa)
        if d > 1.0:
            new += d
            if r.id[1] in epitope:
                ep_new += d
    print("%s: %d residues kept, %d epitope, cut created %.0f A^2 of surface,"
          " %.0f A^2 of it on the epitope" % (tag, len(keep), len(epitope), new, ep_new))
    if ep_new > 150:
        raise SystemExit("%s: the crop changed the epitope itself - it would be designed"
                         " against a surface that does not exist" % tag)
    MAPS[tag] = dict(tag=tag, target=out, chain=REF, epitope=epitope, kept=keep,
                     hotspot=",".join("%s%d" % (REF, n) for n in epitope),
                     new_surface_A2=round(new, 1))

json.dump(MAPS, open("%s/tfr1_targets.json" % W, "w"), indent=0)
print("wrote %s/tfr1_targets.json" % W)

# ---------------------------------------------------------------- the gate
print("\nGATE: a short design, to prove AlphaFold runs before committing eight hours.")
from colabdesign import mk_afdesign_model, clear_mem                      # noqa: E402

tag = "e559"
clear_mem()
m = mk_afdesign_model(protocol="binder", data_dir=W)
try:
    m.prep_inputs(pdb_filename=MAPS[tag]["target"], chain=REF, binder_len=60,
                  hotspot=MAPS[tag]["hotspot"], rm_aa="C")
    RM_AA_OK = True
except TypeError:
    m.prep_inputs(pdb_filename=MAPS[tag]["target"], chain=REF, binder_len=60,
                  hotspot=MAPS[tag]["hotspot"])
    RM_AA_OK = False
print("cysteine can be excluded at prep_inputs:", RM_AA_OK)
m.design_3stage(soft_iters=8, temp_iters=4, hard_iters=2)
seq = m.get_seqs()[0]
print("gate design:", seq)
print("gate i_ptm %.3f  contains cysteine: %s"
      % (float(m.aux["log"].get("i_ptm", 0.0)), "C" in seq))
json.dump(dict(rm_aa_ok=RM_AA_OK), open("%s/tfr1_gate.json" % W, "w"))

del m
clear_mem()
import gc                                                                # noqa: E402
gc.collect()
print("\nGPU state handed to the workers (index, used, free):")
print(subprocess.run("nvidia-smi --query-gpu=index,memory.used,memory.free "
                     "--format=csv,noheader", shell=True,
                     capture_output=True, text=True).stdout.strip(), flush=True)
print("\nSetup OK. Run cell 2.")
