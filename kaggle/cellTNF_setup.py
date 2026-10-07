# === CHALLENGE 2 (TNF-alpha) : SETUP ===
# Reuses the AlphaFold2 parameters from an attached notebook output so we skip the
# 3.5 GB download, and rebuilds the target construct from the PDB rather than
# depending on a file upload.
import glob, os, shutil, subprocess, sys, time
t0 = time.time()


def sh(cmd, check=False):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        print("FAILED:", cmd)
        print((r.stdout + r.stderr)[-1000:])
        if check:
            raise SystemExit("setup failed")
    return r


print("python", sys.version.split()[0], flush=True)
print(subprocess.run("nvidia-smi --query-gpu=name --format=csv,noheader",
                     shell=True, capture_output=True, text=True).stdout.strip(), flush=True)

os.makedirs("/kaggle/working/params", exist_ok=True)
src = glob.glob("/kaggle/input/**/params/params_model_*.npz", recursive=True)
if src:
    print("reusing %d AF2 param files from the attached input" % len(src), flush=True)
    for f in src:
        dst = "/kaggle/working/params/" + os.path.basename(f)
        if not os.path.exists(dst):
            shutil.copy(f, dst)
else:
    print("downloading AF2 params (~3.5 GB) ...", flush=True)
    sh("wget -q -O /kaggle/working/af2.tar "
       "https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar")
    sh("tar -xf /kaggle/working/af2.tar -C /kaggle/working/params")
    sh("rm -f /kaggle/working/af2.tar")
print("AF2 params present:", len(glob.glob("/kaggle/working/params/*.npz")), "files", flush=True)

# ---- target: 1TNF chains A and B (two protomers of the trimer) ------------------
# The receptor site sits BETWEEN protomers, so a single chain cannot present it. All
# three chains would cost ~5x the AlphaFold compute of Challenge 1's target and leave
# about nine designs in budget; two chains cost ~2.7x and allow roughly eighteen.
# Measured residue by residue, dropping chain C creates no artificially exposed
# surface at this particular site (tnf/results/PHASE1_REPORT.md).
sh("pip install -q biopython")
import urllib.request
from Bio.PDB import PDBParser, PDBIO, Select

os.makedirs("/kaggle/working/targets", exist_ok=True)
raw = "/kaggle/working/1TNF.pdb"
if not os.path.exists(raw):
    open(raw, "wb").write(urllib.request.urlopen(
        "https://files.rcsb.org/download/1TNF.pdb", timeout=180).read())


class AB(Select):
    def accept_chain(self, c): return c.id in ("A", "B")
    def accept_residue(self, r): return r.id[0] == " "
    def accept_atom(self, a): return a.element != "H"


st = PDBParser(QUIET=True).get_structure("t", raw)
io = PDBIO(); io.set_structure(st)
TARGET = "/kaggle/working/targets/TNF_AB.pdb"
io.save(TARGET, AB())
counts = {}
for line in open(TARGET):
    if line.startswith("ATOM") and line[12:16].strip() == "CA":
        counts[line[21]] = counts.get(line[21], 0) + 1
print("target construct chains:", counts, "total", sum(counts.values()), flush=True)
assert sorted(counts) == ["A", "B"] and sum(counts.values()) > 290

# ---- JAX pinned for ColabDesign -------------------------------------------------
PROBE = ("python -c \"import jax, jax.lib; assert hasattr(jax.lib,'xla_bridge'); "
         "print(jax.__version__)\"")
r = subprocess.run(PROBE, shell=True, capture_output=True, text=True)
if r.returncode != 0 or not r.stdout.strip().startswith("0.4.38"):
    print("installing jax 0.4.38 ...", flush=True)
    sh("pip install -q 'jax[cuda12]==0.4.38' 'jaxlib==0.4.38'")
sh("pip install -q git+https://github.com/sokrypton/ColabDesign.git")

# ---- the thing actually being smoke-tested: a TWO-CHAIN target with chain-prefixed
# hotspots. Challenge 1 used a single chain, so neither has been exercised here.
HOTSPOT = "A86,A87,A88,A90,A91,A92,B29,B32,B33,B34,B35,B36"
chk = subprocess.run(
    "cd /kaggle/working && python -c \""
    "import jax; from colabdesign import mk_afdesign_model; "
    "m = mk_afdesign_model(protocol='binder', data_dir='/kaggle/working'); "
    "m.prep_inputs(pdb_filename='/kaggle/working/targets/TNF_AB.pdb', "
    "chain='A,B', binder_len=70, hotspot='%s'); "
    "print('READY', jax.__version__, jax.devices()); "
    "print('target length seen by the model:', m._target_len); "
    "print('binder length:', m._binder_len)\"" % HOTSPOT,
    shell=True, capture_output=True, text=True)
print(chk.stdout[-700:])
if chk.returncode != 0:
    print("TWO-CHAIN PREP FAILED:")
    print(chk.stderr[-2000:])
print("TNF SETUP DONE in %.1f min" % ((time.time() - t0) / 60))
