"""Generate a self-contained Kaggle notebook that produces Challenge 1 binder designs."""
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def code(src):
    return {"cell_type": "code", "execution_count": None, "metadata": {},
            "outputs": [], "source": src.rstrip("\n").split("\n")}


def md(src):
    return {"cell_type": "markdown", "metadata": {},
            "source": src.rstrip("\n").split("\n")}


cells = []

cells.append(md("""# Challenge 1 - EGFR binder generation (Kaggle)

Anthropic x Adaptyv Protein Design Competition. Epitope and hotspots come from the
Phase 1 cross-species conservation analysis.

**Before running:** Settings -> Accelerator -> **GPU P100**, and Settings -> **Internet: On**
(internet requires phone verification on your Kaggle account).

Pipeline: RFdiffusion (backbones) -> ProteinMPNN soluble (sequences) -> FASTA.
Run cells in order. Cell 1 fails fast if the environment is wrong."""))

cells.append(md("## 1. Environment check (run first - fails fast)"))
cells.append(code(
    'import subprocess, urllib.request\n'
    'print(subprocess.run(["nvidia-smi","--query-gpu=name,memory.total","--format=csv"],\n'
    '                     capture_output=True, text=True).stdout)\n'
    'try:\n'
    '    urllib.request.urlopen("https://pypi.org", timeout=15)\n'
    '    print("internet: OK")\n'
    'except Exception as e:\n'
    '    raise SystemExit("INTERNET IS OFF - enable it in Settings. " + str(e))\n'
    'import torch\n'
    'print("torch", torch.__version__, "cuda:", torch.cuda.is_available())\n'
    'if not torch.cuda.is_available():\n'
    '    raise SystemExit("No GPU - set Accelerator to GPU P100 in Settings.")'))

cells.append(md("""## 2. Build the target construct

Downloads PDB 6ARU and trims to EGFR domain I (chain A, PDB residues 1-164) with the
cetuximab Fab removed. Nothing to upload.

Epitope A = F44, E45, D46, L49, Q52, Y74, D75 (UniProt) = PDB 20, 21, 22, 25, 28, 50, 51.
100% identical human/mouse, overlaps the EGF ligand site, and verified to have the same
solvent exposure in this trimmed construct as in the intact receptor."""))
cells.append(code(
    'import urllib.request, os\n'
    'from Bio.PDB import PDBParser, PDBIO, Select\n'
    '\n'
    'os.makedirs("/kaggle/working/targets", exist_ok=True)\n'
    'raw = "/kaggle/working/6ARU.pdb"\n'
    'if not os.path.exists(raw):\n'
    '    data = urllib.request.urlopen(\n'
    '        "https://files.rcsb.org/download/6ARU.pdb", timeout=120).read()\n'
    '    open(raw, "wb").write(data)\n'
    '\n'
    'class DomainI(Select):\n'
    '    def accept_chain(self, c): return c.id == "A"\n'
    '    def accept_residue(self, r): return r.id[0] == " " and 1 <= r.id[1] <= 164\n'
    '\n'
    'st = PDBParser(QUIET=True).get_structure("x", raw)\n'
    'io = PDBIO(); io.set_structure(st)\n'
    'TARGET = "/kaggle/working/targets/EGFR_A_domainI.pdb"\n'
    'io.save(TARGET, DomainI())\n'
    'n = sum(1 for l in open(TARGET)\n'
    '        if l.startswith("ATOM") and l[12:16].strip() == "CA")\n'
    'print("wrote", TARGET, "with", n, "residues")\n'
    'assert 150 < n < 170, "unexpected residue count"'))

cells.append(md("## 3. Install RFdiffusion (~10 min)"))
cells.append(code(
    '%%bash\n'
    'set -e\n'
    'cd /kaggle/working\n'
    '[ -d RFdiffusion ] || git clone -q https://github.com/RosettaCommons/RFdiffusion.git\n'
    'cd RFdiffusion && mkdir -p models && cd models\n'
    'for f in 6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt \\\n'
    '         e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt \\\n'
    '         f572d396fae9206628714fb2ce00f72e/Complex_beta_ckpt.pt ; do\n'
    '  name=$(basename $f)\n'
    '  [ -f "$name" ] || wget -q "http://files.ipd.uw.edu/pub/RFdiffusion/$f"\n'
    'done\n'
    'ls -la'))

cells.append(code(
    '%%bash\n'
    'set -e\n'
    'pip install -q hydra-core omegaconf icecream pyrsistent e3nn opt_einsum_fx 2>&1 | tail -2\n'
    'pip install -q --no-deps dgl -f https://data.dgl.ai/wheels/torch-2.4/cu124/repo.html 2>&1 | tail -2\n'
    'cd /kaggle/working/RFdiffusion/env/SE3Transformer && pip install -q --no-deps . 2>&1 | tail -2\n'
    'cd /kaggle/working/RFdiffusion && pip install -q --no-deps -e . 2>&1 | tail -2\n'
    'echo INSTALL_DONE'))

cells.append(md("""## 4. Run RFdiffusion

`contigs=[A1-164/0 55-95]` keeps the target and diffuses a 55-95 aa binder, which keeps
every design inside the competition's minibinder category (40-100 aa). `noise_scale=0`
is the standard setting for binder design.

Start with `NUM = 20` to confirm the pipeline runs end to end, then raise it. Kaggle
sessions cap at 12 hours."""))
cells.append(code(
    'import os\n'
    'os.environ["DGLBACKEND"] = "pytorch"\n'
    '\n'
    'HOTSPOTS = "[A20,A21,A22,A25,A28,A50,A51]"\n'
    'CONTIGS  = "[A1-164/0 55-95]"\n'
    'NUM      = 20      # raise once a first batch succeeds\n'
    'PREFIX   = "/kaggle/working/out/egfr_a"\n'
    'CKPT     = ""      # beta batch: see the next cell\n'
    '\n'
    'os.makedirs("/kaggle/working/out", exist_ok=True)\n'
    'cmd = ("cd /kaggle/working/RFdiffusion && python scripts/run_inference.py"\n'
    '       " inference.output_prefix=" + PREFIX +\n'
    '       " inference.input_pdb=/kaggle/working/targets/EGFR_A_domainI.pdb"\n'
    '       " \'contigmap.contigs=" + CONTIGS + "\'"\n'
    '       " \'ppi.hotspot_res=" + HOTSPOTS + "\'"\n'
    '       " inference.num_designs=" + str(NUM) +\n'
    '       " denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0 " + CKPT)\n'
    'print(cmd)\n'
    'get_ipython().system(cmd)'))

cells.append(md("""### Optional second batch: the beta checkpoint

Only ~15 of 354 binders in Anthropic's campaign had design models with >=20% beta strand,
because the default methods overwhelmingly produce all-alpha folds. Track 3 selection
weighs design novelty, so a minority of beta-containing backbones is cheap
differentiation. Rerun cell 4 with these two values changed."""))
cells.append(code(
    '# PREFIX = "/kaggle/working/out/egfr_a_beta"\n'
    '# CKPT = ("inference.ckpt_override_path="\n'
    '#         "/kaggle/working/RFdiffusion/models/Complex_beta_ckpt.pt")\n'
    'print("edit cell 4 with the two lines above, then rerun it")'))

cells.append(md("""## 5. ProteinMPNN - sequence design

Designs the binder chain (B) only, with target chain A held fixed. Soluble weights, no
cysteines.

Histidines are deliberately NOT placed here. Switch positions are chosen in Phase 3
against specific conserved carboxylates (E45, D46, D75, E97), because the published
pH-switch work found that random histidine incorporation rarely produces a working
switch."""))
cells.append(code(
    '%%bash\n'
    'set -e\n'
    'cd /kaggle/working\n'
    '[ -d ProteinMPNN ] || git clone -q https://github.com/dauparas/ProteinMPNN.git\n'
    'echo PROTEINMPNN_READY'))
cells.append(code(
    'import glob, os, subprocess\n'
    '\n'
    'pdbs = sorted(glob.glob("/kaggle/working/out/egfr_a*.pdb"))\n'
    'print(len(pdbs), "backbones from RFdiffusion")\n'
    'os.makedirs("/kaggle/working/seqs", exist_ok=True)\n'
    'for p in pdbs:\n'
    '    cmd = ["python", "/kaggle/working/ProteinMPNN/protein_mpnn_run.py",\n'
    '           "--pdb_path", p, "--pdb_path_chains", "B",\n'
    '           "--out_folder", "/kaggle/working/seqs",\n'
    '           "--num_seq_per_target", "8", "--sampling_temp", "0.1",\n'
    '           "--omit_AAs", "CX", "--use_soluble_model", "--seed", "37"]\n'
    '    r = subprocess.run(cmd, capture_output=True, text=True)\n'
    '    if r.returncode != 0:\n'
    '        print(os.path.basename(p), "FAILED:", r.stderr[-400:])\n'
    '        break\n'
    'print("sequence design done")'))

cells.append(md("## 6. Collect into one FASTA and zip for download"))
cells.append(code(
    'import glob, os, zipfile\n'
    '\n'
    'records = []\n'
    'for f in sorted(glob.glob("/kaggle/working/seqs/seqs/*.fa")):\n'
    '    name = os.path.basename(f).replace(".fa", "")\n'
    '    entries, hdr = [], None\n'
    '    for line in open(f):\n'
    '        line = line.strip()\n'
    '        if line.startswith(">"):\n'
    '            hdr = line\n'
    '        elif line:\n'
    '            entries.append((hdr, line))\n'
    '    # entry 0 is the input backbone sequence; the designs follow\n'
    '    for i, (h, s) in enumerate(entries[1:], start=1):\n'
    '        records.append((name + "_s" + str(i), s))\n'
    '\n'
    'out = "/kaggle/working/designs.fasta"\n'
    'with open(out, "w") as fh:\n'
    '    for n, s in records:\n'
    '        fh.write(">" + n + "\\n" + s + "\\n")\n'
    'print(len(records), "sequences ->", out)\n'
    'if records:\n'
    '    print("example:", records[0][0], records[0][1][:60], "...")\n'
    '\n'
    'with zipfile.ZipFile("/kaggle/working/challenge1_outputs.zip", "w",\n'
    '                     zipfile.ZIP_DEFLATED) as z:\n'
    '    for p in glob.glob("/kaggle/working/out/*.pdb"):\n'
    '        z.write(p, "pdb/" + os.path.basename(p))\n'
    '    if os.path.exists(out):\n'
    '        z.write(out, "designs.fasta")\n'
    'print("zipped -> challenge1_outputs.zip (download from the Output panel)")'))

nb = {"cells": cells,
      "metadata": {
          "kernelspec": {"display_name": "Python 3", "language": "python",
                         "name": "python3"},
          "language_info": {"name": "python"},
          "accelerator": "GPU"},
      "nbformat": 4, "nbformat_minor": 5}

os.makedirs(os.path.join(ROOT, "kaggle"), exist_ok=True)
path = os.path.join(ROOT, "kaggle", "challenge1_egfr_kaggle.ipynb")
json.dump(nb, open(path, "w"), indent=1)
print("wrote", path, "with", len(cells), "cells")
