# === CHALLENGE 2 (TNF-alpha) : PRODUCTION SETUP ===
# Builds the CROPPED target. The smoke test measured 43.3 s per iteration on the full
# two-chain construct (304 + 70 residues) - 7.2x Challenge 1's rate for a 1.7x length
# increase, so the cost scales nearer L^3.7 than L^2. At that rate the whole GPU budget
# buys about six designs, which is not worth running.
#
# Cropping to a 22 A sphere around the switch site keeps all 25 site residues, both
# switch anchors (K166 of chain A, R108 of chain B) and - measured against the intact
# trimer - introduces ZERO artificially exposed surface at the site, fewer than the
# uncropped two-chain construct, because cropping also removes residues that chain C
# had been covering.
import glob, json, os, shutil, subprocess, sys, time
t0 = time.time()
OFF = 76
RADIUS = 22.0


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

sh("pip install -q biopython")
import urllib.request
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import protein_letters_3to1

os.makedirs("/kaggle/working/targets", exist_ok=True)
raw = "/kaggle/working/1TNF.pdb"
if not os.path.exists(raw):
    open(raw, "wb").write(urllib.request.urlopen(
        "https://files.rcsb.org/download/1TNF.pdb", timeout=180).read())

P = PDBParser(QUIET=True)
mdl = P.get_structure("m", raw)[0]
ca = {}
for ch in mdl:
    for r in ch:
        if r.id[0] == " " and "CA" in r:
            ca[(ch.id, r.id[1])] = r["CA"].get_coord()
centre = (ca[("A", 166 - OFF)] + ca[("B", 108 - OFF)]) / 2
keep = {k for k, v in ca.items() if k[0] in ("A", "B")
        and float(np.linalg.norm(v - centre)) <= RADIUS}


class Crop(Select):
    def accept_chain(self, c): return c.id in ("A", "B")
    def accept_residue(self, r): return r.id[0] == " " and (r.parent.id, r.id[1]) in keep
    def accept_atom(self, a): return a.element != "H"


TARGET = "/kaggle/working/targets/TNF_site_crop.pdb"
st = P.get_structure("t", raw)
io = PDBIO(); io.set_structure(st); io.save(TARGET, Crop())

# The map from construct order to the original chain and UniProt number. AfDesign
# CONCATENATES the target chains into one chain on output and renumbers from 1 - the
# smoke test showed chains {'A': 608, 'B': 140} where three chains were expected - so
# position in this list IS the output residue number minus one. Everything downstream
# maps through this rather than assuming a numbering.
order = []
for ch in P.get_structure("c", TARGET)[0]:
    for r in ch:
        if r.id[0] == " ":
            order.append(dict(chain=ch.id, pdb=r.id[1], uniprot=r.id[1] + OFF,
                              aa=protein_letters_3to1.get(r.get_resname(), "X")))
json.dump(dict(radius=RADIUS, n=len(order), order=order),
          open("/kaggle/working/construct_map.json", "w"))

SITE = [95, 105, 106, 107, 108, 109, 110, 111, 112, 155, 159, 161, 162, 163, 164, 165,
        166, 167, 168, 201, 203, 207, 220, 221, 223]
have = {(o["chain"], o["uniprot"]) for o in order}
print("cropped target: %d residues (A %d, B %d)"
      % (len(order), sum(1 for o in order if o["chain"] == "A"),
         sum(1 for o in order if o["chain"] == "B")), flush=True)
print("site residues retained: %d of %d"
      % (len({o["uniprot"] for o in order} & set(SITE)), len(SITE)), flush=True)
assert ("A", 166) in have and ("B", 108) in have, "a switch anchor was cropped away"
assert len({o["uniprot"] for o in order} & set(SITE)) == len(SITE), "site residue lost"
print("switch anchors present: K166(A), R108(B)", flush=True)

HOTSPOT = "A86,A87,A88,A90,A91,A92,B29,B32,B33,B34,B35,B36"

PROBE = ("python -c \"import jax, jax.lib; assert hasattr(jax.lib,'xla_bridge'); "
         "print(jax.__version__)\"")
r = subprocess.run(PROBE, shell=True, capture_output=True, text=True)
if r.returncode != 0 or not r.stdout.strip().startswith("0.4.38"):
    print("installing jax 0.4.38 ...", flush=True)
    sh("pip install -q 'jax[cuda12]==0.4.38' 'jaxlib==0.4.38'")
sh("pip install -q git+https://github.com/sokrypton/ColabDesign.git")

chk = subprocess.run(
    "cd /kaggle/working && python -c \""
    "import jax; from colabdesign import mk_afdesign_model; "
    "m = mk_afdesign_model(protocol='binder', data_dir='/kaggle/working'); "
    "m.prep_inputs(pdb_filename='%s', chain='A,B', binder_len=65, hotspot='%s'); "
    "print('READY', jax.__version__, jax.devices()); "
    "print('target length seen by the model:', m._target_len)\"" % (TARGET, HOTSPOT),
    shell=True, capture_output=True, text=True)
print(chk.stdout[-500:])
if chk.returncode != 0:
    print("PREP FAILED:"); print(chk.stderr[-1500:])
print("TNF PRODUCTION SETUP DONE in %.1f min" % ((time.time() - t0) / 60))
