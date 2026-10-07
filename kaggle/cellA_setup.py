# === CHALLENGE 1 : SETUP (run once, ~15 min) ===
import os, subprocess, time, urllib.request
t0 = time.time()

def sh(cmd, check=True):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        print("FAILED:", cmd)
        print(r.stdout[-1200:])
        print(r.stderr[-1200:])
        if check:
            raise SystemExit("setup failed")
    return r

print(subprocess.run("nvidia-smi --query-gpu=name,memory.total --format=csv",
                     shell=True, capture_output=True, text=True).stdout)
import torch
TORCH = torch.__version__
print("torch", TORCH, "| cuda available:", torch.cuda.is_available())
assert torch.cuda.is_available(), "GPU not enabled"
urllib.request.urlopen("https://pypi.org", timeout=20)
print("internet: OK", flush=True)

# Biopython is not in this Kaggle image.
sh("pip install -q biopython", check=False)

# ---- target construct: EGFR domain I from 6ARU, Fab removed ----
from Bio.PDB import PDBParser, PDBIO, Select
os.makedirs("/kaggle/working/targets", exist_ok=True)
raw = "/kaggle/working/6ARU.pdb"
if not os.path.exists(raw):
    open(raw, "wb").write(urllib.request.urlopen(
        "https://files.rcsb.org/download/6ARU.pdb", timeout=180).read())

class DomainI(Select):
    def accept_chain(self, c):
        return c.id == "A"
    def accept_residue(self, r):
        return r.id[0] == " " and 1 <= r.id[1] <= 164

st = PDBParser(QUIET=True).get_structure("x", raw)
io = PDBIO(); io.set_structure(st)
TARGET = "/kaggle/working/targets/EGFR_A_domainI.pdb"
io.save(TARGET, DomainI())
n = sum(1 for l in open(TARGET) if l.startswith("ATOM") and l[12:16].strip() == "CA")
print("target construct residues:", n, flush=True)
assert 150 < n < 170, "unexpected residue count"

# ---- RFdiffusion ----
if not os.path.isdir("/kaggle/working/RFdiffusion"):
    sh("cd /kaggle/working && git clone -q https://github.com/RosettaCommons/RFdiffusion.git")
os.makedirs("/kaggle/working/RFdiffusion/models", exist_ok=True)
for f in ["6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt",
          "e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt",
          "f572d396fae9206628714fb2ce00f72e/Complex_beta_ckpt.pt"]:
    name = f.split("/")[-1]
    p = "/kaggle/working/RFdiffusion/models/" + name
    if not os.path.exists(p):
        print("downloading", name, flush=True)
        sh("wget -q -O %s http://files.ipd.uw.edu/pub/RFdiffusion/%s" % (p, f))
print("model weights present", flush=True)

sh("pip install -q hydra-core omegaconf icecream pyrsistent e3nn opt_einsum_fx", check=False)

# dgl is the fragile dependency: its wheels are built per torch version, and this
# image ships a much newer torch than any published dgl wheel targets. Try the
# generic index first, then torch-specific repos, and report what actually imports.
def dgl_ok():
    r = subprocess.run("python -c \"import dgl; print(dgl.__version__)\"",
                       shell=True, capture_output=True, text=True)
    return (r.returncode == 0, (r.stdout + r.stderr).strip()[-200:])

ok, msg = dgl_ok()
if not ok:
    for src in ["dgl",
                "dgl -f https://data.dgl.ai/wheels/torch-2.4/cu124/repo.html",
                "dgl -f https://data.dgl.ai/wheels/torch-2.3/cu121/repo.html"]:
        print("trying: pip install", src, flush=True)
        sh("pip install -q --no-deps " + src, check=False)
        ok, msg = dgl_ok()
        if ok:
            break
print("dgl import:", "OK " + msg if ok else "FAILED -> " + msg, flush=True)

sh("cd /kaggle/working/RFdiffusion/env/SE3Transformer && pip install -q --no-deps .", check=False)
sh("cd /kaggle/working/RFdiffusion && pip install -q --no-deps -e .", check=False)

if not os.path.isdir("/kaggle/working/ProteinMPNN"):
    sh("cd /kaggle/working && git clone -q https://github.com/dauparas/ProteinMPNN.git")

# ---- verify RFdiffusion can actually import before we rely on it ----
chk = subprocess.run(
    "cd /kaggle/working/RFdiffusion && python -c \"import dgl, se3_transformer; "
    "from rfdiffusion.inference import model_runners; print('RFDIFFUSION IMPORT OK')\"",
    shell=True, capture_output=True, text=True)
print(chk.stdout[-800:])
if chk.returncode != 0:
    print("RFDIFFUSION IMPORT FAILED:")
    print(chk.stderr[-1500:])

print("SETUP DONE in %.1f min" % ((time.time() - t0) / 60))
