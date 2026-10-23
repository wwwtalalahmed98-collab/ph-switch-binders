# === CHALLENGE 3 (GLP1R) : TARGET SETUP ===
#
# Builds the cropped active-state receptor from PDB 6X18, and saves the geometry needed
# to ask the only question that matters here: how DEEP does a candidate reach.
#
# GLP1R is a class B1 GPCR and works by a two-domain mechanism. The hormone's C-terminal
# half grips the extracellular domain; its first eight residues thread ~40 A down into
# the seven-helix bundle, and that is what activates the receptor. Removing those eight
# residues gives exendin 9-39, a real drug that binds well and ANTAGONISES. So a design
# that only binds the extracellular groove - which is what the challenge page recommends,
# and what a folded minibinder would naturally do - cannot win, because the winner has to
# activate.
#
# Measured from 6X18 itself: GLP-1 spans 39.6 A of depth, makes 22 extracellular contacts
# and 30 transmembrane ones, and 30 residues of alpha helix rises about 45 A. The agonist
# is essentially one long helix threading a channel, not a globular domain sitting on a
# surface.
import glob, json, os, shutil, subprocess, sys, time
t0 = time.time()

def sh(cmd, check=False):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        print("FAILED:", cmd)
        print((r.stdout + r.stderr)[-1000:])
        if check:
            raise SystemExit("setup failed")
    return r

print(sys.version.split()[0], flush=True)
sh("nvidia-smi --query-gpu=name --format=csv,noheader")

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
       "https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar", check=True)
    sh("tar -xf /kaggle/working/af2.tar -C /kaggle/working/params")
print("AF2 params present:", len(glob.glob("/kaggle/working/params/*.npz")), "files", flush=True)

sh("pip install -q biopython")
import urllib.request
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import protein_letters_3to1 as T31

os.makedirs("/kaggle/working/targets", exist_ok=True)
raw = "/kaggle/working/6X18.pdb"
if not os.path.exists(raw):
    open(raw, "wb").write(urllib.request.urlopen(
        "https://files.rcsb.org/download/6X18.pdb", timeout=180).read())

RADIUS, MAXGAP, ECD_END = 10.0, 4, 145
mdl = PDBParser(QUIET=True).get_structure("x", raw)[0]
rec = [r for r in mdl["R"] if r.id[0] == " "]
pep = [r for r in mdl["P"] if r.id[0] == " "]
native_seq = "".join(T31.get(r.get_resname(), "X") for r in pep)
print("6X18 chain R: %d residues; chain P (GLP-1): %d-%d %s"
      % (len(rec), pep[0].id[1], pep[-1].id[1], native_seq), flush=True)
if native_seq[0] != "H" or pep[0].id[1] != 7:
    raise SystemExit("chain P does not start at His7 - the depth reference would be wrong")

patoms = np.array([a.get_coord() for r in pep for a in r])
mind, ratoms = {}, {}
for r in rec:
    ra = np.array([a.get_coord() for a in r])
    ratoms[r.id[1]] = ra
    mind[r.id[1]] = float(np.min(np.linalg.norm(ra[:, None, :] - patoms[None, :, :], axis=-1)))

keep = sorted(k for k, v in mind.items() if v <= RADIUS)
present = {r.id[1] for r in rec}
filled = set(keep)
for a, b in zip(keep, keep[1:]):
    if 1 < b - a <= MAXGAP + 1 and all(x in present for x in range(a + 1, b)):
        filled.update(range(a + 1, b))
kept = sorted(filled)

contacts = set()
for r in pep:
    ra = np.array([a.get_coord() for a in r])
    for rid, rra in ratoms.items():
        if np.min(np.linalg.norm(ra[:, None, :] - rra[None, :, :], axis=-1)) <= 4.5:
            contacts.add(rid)
lost = sorted(contacts - filled)
if lost:
    raise SystemExit("crop drops GLP-1 contacts %s - widen RADIUS" % lost)
print("cropped target: %d residues, all %d GLP-1 contacts retained"
      % (len(kept), len(contacts)), flush=True)

class Keep(Select):
    def accept_chain(self, c): return c.id == "R"
    def accept_residue(self, r): return r.id[0] == " " and r.id[1] in filled
    def accept_atom(self, a): return a.element != "H"

TARGET = "/kaggle/working/targets/GLP1R_crop.pdb"
io = PDBIO(); io.set_structure(mdl); io.save(TARGET, Keep())

# Depth axis: extracellular-domain centroid -> transmembrane centroid, so one number per
# atom says how far into the bundle it sits. Saved so every later run uses the same axis.
def ca(r): return r["CA"].get_coord() if "CA" in r else None
kept_res = [r for r in rec if r.id[1] in filled]
ecd_c = np.mean([ca(r) for r in kept_res if r.id[1] <= ECD_END and ca(r) is not None], axis=0)
tm_c = np.mean([ca(r) for r in kept_res if r.id[1] > ECD_END and ca(r) is not None], axis=0)
axis = (tm_c - ecd_c) / np.linalg.norm(tm_c - ecd_c)

native_depth = {}
for r in pep:
    c = ca(r)
    if c is not None:
        native_depth[r.id[1]] = round(float(np.dot(c - ecd_c, axis)), 2)
ref_deep = max(native_depth.values())
print("GLP-1 native depth: H7 at %.1f A, C-terminus at %.1f A, span %.1f A"
      % (native_depth[7], native_depth[pep[-1].id[1]], ref_deep - min(native_depth.values())),
      flush=True)

cmap = [dict(uniprot=r.id[1], aa=T31.get(r.get_resname(), "X"),
             depth=round(float(np.dot(ca(r) - ecd_c, axis)), 2) if ca(r) is not None else None,
             region="ecd" if r.id[1] <= ECD_END else "tm",
             contact=r.id[1] in contacts)
        for r in kept_res]
json.dump(dict(order=cmap, ecd_end=ECD_END, radius=RADIUS,
               ecd_centroid=[float(x) for x in ecd_c], axis=[float(x) for x in axis],
               native_glp1=dict(seq=native_seq, start=pep[0].id[1], depth=native_depth),
               tm_contacts=sorted(x for x in contacts if x > ECD_END),
               ecd_contacts=sorted(x for x in contacts if x <= ECD_END)),
          open("/kaggle/working/glp1r_map.json", "w"), indent=1)

TM_HOT = ",".join("R%d" % x for x in sorted(x for x in contacts if x > ECD_END))
ECD_HOT = ",".join("R%d" % x for x in sorted(x for x in contacts if x <= ECD_END))
open("/kaggle/working/glp1r_tm_hotspot.txt", "w").write(TM_HOT)
open("/kaggle/working/glp1r_ecd_hotspot.txt", "w").write(ECD_HOT)
print("activation hotspots (transmembrane): %d" % len(TM_HOT.split(",")), flush=True)
print("binding hotspots (extracellular)   : %d" % len(ECD_HOT.split(",")), flush=True)

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
    "m.prep_inputs(pdb_filename='%s', chain='R', binder_len=30, hotspot='%s'); "
    "print('READY', jax.__version__, jax.devices()); "
    "print('target length seen by the model:', m._target_len)\"" % (TARGET, TM_HOT),
    shell=True, capture_output=True, text=True)
print(chk.stdout[-600:])
if chk.returncode != 0:
    print("PREP FAILED:"); print(chk.stderr[-1500:])
print("GLP1R SETUP DONE in %.1f min" % ((time.time() - t0) / 60))
