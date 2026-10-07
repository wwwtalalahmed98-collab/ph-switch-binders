# === CHALLENGE 1 : REDESIGN NOTEBOOK — SETUP ===
# Reuses the generation run's output (attached as an input): backbones AND the
# AlphaFold2 parameters, so we skip the 3.5 GB download.
import glob, os, shutil, subprocess, sys, time
t0 = time.time()

def sh(cmd, check=False):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        print("FAILED:", cmd); print((r.stdout + r.stderr)[-1000:])
        if check: raise SystemExit("setup failed")
    return r

print("python", sys.version.split()[0], flush=True)
print(subprocess.run("nvidia-smi --query-gpu=name --format=csv,noheader",
                     shell=True, capture_output=True, text=True).stdout.strip(), flush=True)

# ---- what did we actually get attached? ----
roots = sorted(glob.glob("/kaggle/input/*"))
print("inputs:", roots, flush=True)
bb = glob.glob("/kaggle/input/**/egfr_a_*.pdb", recursive=True)
print("backbone PDBs visible:", len(bb), flush=True)
if not bb:
    raise SystemExit("No backbones found — attach the generation notebook's output.")

# ---- AF2 params: reuse from the input if present, else download ----
os.makedirs("/kaggle/working/params", exist_ok=True)
src = glob.glob("/kaggle/input/**/params/params_model_*.npz", recursive=True)
if src:
    print("reusing %d AF2 param files from input" % len(src), flush=True)
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
print("AF2 params present:",
      len(glob.glob("/kaggle/working/params/*.npz")), "files", flush=True)

# ---- target construct (small, rebuild rather than depend on the input layout) ----
sh("pip install -q biopython")
import urllib.request
from Bio.PDB import PDBParser, PDBIO, Select
os.makedirs("/kaggle/working/targets", exist_ok=True)
raw = "/kaggle/working/6ARU.pdb"
if not os.path.exists(raw):
    open(raw, "wb").write(urllib.request.urlopen(
        "https://files.rcsb.org/download/6ARU.pdb", timeout=180).read())
class DomainI(Select):
    def accept_chain(self, c): return c.id == "A"
    def accept_residue(self, r): return r.id[0] == " " and 1 <= r.id[1] <= 164
st = PDBParser(QUIET=True).get_structure("x", raw)
io = PDBIO(); io.set_structure(st)
io.save("/kaggle/working/targets/EGFR_A_domainI.pdb", DomainI())
n = sum(1 for l in open("/kaggle/working/targets/EGFR_A_domainI.pdb")
        if l.startswith("ATOM") and l[12:16].strip() == "CA")
print("target construct residues:", n, flush=True)
assert 150 < n < 170

# ---- JAX pinned for ColabDesign, then ColabDesign and ProteinMPNN ----
PROBE_V = ("python -c \"import jax, jax.lib; assert hasattr(jax.lib,'xla_bridge'); "
           "print(jax.__version__)\"")
r = subprocess.run(PROBE_V, shell=True, capture_output=True, text=True)
if r.returncode != 0 or not r.stdout.strip().startswith("0.4.38"):
    print("installing jax 0.4.38 ...", flush=True)
    sh("pip install -q 'jax[cuda12]==0.4.38' 'jaxlib==0.4.38'")
sh("pip install -q git+https://github.com/sokrypton/ColabDesign.git")
if not os.path.isdir("/kaggle/working/ProteinMPNN"):
    sh("cd /kaggle/working && git clone -q https://github.com/dauparas/ProteinMPNN.git")

chk = subprocess.run(
    "cd /kaggle/working && python -c \""
    "import jax; from colabdesign import mk_afdesign_model; "
    "m = mk_afdesign_model(protocol='binder', data_dir='/kaggle/working'); "
    "m.prep_inputs(pdb_filename='/kaggle/working/targets/EGFR_A_domainI.pdb', "
    "chain='A', binder_len=60, hotspot='20,21,22,25,28,50,51'); "
    "print('READY', jax.__version__, jax.devices())\"",
    shell=True, capture_output=True, text=True)
print(chk.stdout[-400:])
if chk.returncode != 0:
    print("SETUP CHECK FAILED:"); print(chk.stderr[-1500:])
print("REDESIGN SETUP DONE in %.1f min" % ((time.time() - t0) / 60))
