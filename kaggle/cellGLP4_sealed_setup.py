# === CHALLENGE 3 (GLP1R) : SEALED TARGET, AND A GATE BEFORE WE SPEND ON IT ===
#
# Why this exists. Ten of eleven threading designs on the previous target reached into
# space that only existed because we cropped the bottom off the receptor. Below 45 A our
# crop held 4 atoms where the intact receptor holds 211. Confidence correlated -0.963 with
# landing on the hormone's own site: the optimiser found the flaw in our preparation
# faster and more consistently than it found the binding mode.
#
# The fix is to keep the bundle core below the hormone's reach, so the floor is protein
# instead of vacuum. Measured on 6X18: adding every residue whose CA lies between 28 and
# 62 A deep and within 18 A of the receptor axis takes the target from 192 to 283 residues
# and lifts depth-band coverage from 91/55/2/0% to 100/100/96/74%.
#
# THE SECOND HALF OF THIS CELL IS A GATE. Sealing the floor could disturb the groove
# docking it is meant to protect, and we would not find out for six hours. So the native
# hormone and the blocker-shaped truncation are predicted against the new target FIRST,
# and if the hormone no longer docks the way the crystal structure says it should, this
# cell raises and the generation cell never runs. Cost of being wrong: 25 minutes instead
# of a night.
import glob, json, os, shutil, subprocess, sys, time
t0 = time.time()

# THIS is why the two-domain worker died 20+ times across three attempts, and it was
# never the binder length. The gate below runs AlphaFold IN THIS PROCESS, and this
# process does not exit - it goes on to launch the workers and wait on them. With jax's
# default allocator that gate permanently reserves most of GPU 0. The worker handed
# GPU 0 then has ~6 GB of headroom and asks for 6-7 GB; the worker on the untouched
# GPU 1 has the whole card and runs fine. Identical target, overlapping binder lengths,
# one fails always and one never does, and every failure names "device ordinal 0".
# Shrinking the binder from 45-70 to 40-50 to 30-38 could not fix a fixed deficit.
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ["XLA_PYTHON_CLIENT_ALLOCATOR"] = "platform"


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
    sh("wget -q -O /kaggle/working/af2.tar "
       "https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar", check=True)
    sh("tar -xf /kaggle/working/af2.tar -C /kaggle/working/params")
print("AF2 params present:", len(glob.glob("/kaggle/working/params/*.npz")), flush=True)

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
SEAL_R, SEAL_LO, SEAL_HI = 18.0, 28.0, 62.0

mdl = PDBParser(QUIET=True).get_structure("x", raw)[0]
rec = [r for r in mdl["R"] if r.id[0] == " "]
pep = [r for r in mdl["P"] if r.id[0] == " "]
native_seq = "".join(T31.get(r.get_resname(), "X") for r in pep)
if native_seq[0] != "H" or pep[0].id[1] != 7:
    raise SystemExit("chain P does not start at His7 - the depth reference would be wrong")

patoms = np.array([a.get_coord() for r in pep for a in r])
present = {r.id[1] for r in rec}
mind, ratoms = {}, {}
for r in rec:
    ra = np.array([a.get_coord() for a in r])
    ratoms[r.id[1]] = ra
    mind[r.id[1]] = float(np.min(np.linalg.norm(ra[:, None, :] - patoms[None, :, :], axis=-1)))

base = sorted(k for k, v in mind.items() if v <= RADIUS)
filled = set(base)
for a, b in zip(base, base[1:]):
    if 1 < b - a <= MAXGAP + 1 and all(x in present for x in range(a + 1, b)):
        filled.update(range(a + 1, b))


def ca(r):
    return r["CA"].get_coord() if "CA" in r else None


prelim = [r for r in rec if r.id[1] in filled]
ecd_c = np.mean([ca(r) for r in prelim if r.id[1] <= ECD_END and ca(r) is not None], axis=0)
tm_c = np.mean([ca(r) for r in prelim if r.id[1] > ECD_END and ca(r) is not None], axis=0)
AXIS = (tm_c - ecd_c) / np.linalg.norm(tm_c - ecd_c)
OFF = float(ecd_c @ AXIS)


def depth(x):
    return float(np.dot(x - ecd_c, AXIS))


def radial(x):
    v = x - ecd_c
    return float(np.linalg.norm(v - np.dot(v, AXIS) * AXIS))


# seal the floor
for r in rec:
    c = ca(r)
    if c is not None and SEAL_LO <= depth(c) <= SEAL_HI and radial(c) <= SEAL_R:
        filled.add(r.id[1])
ks = sorted(filled)
for a, b in zip(ks, ks[1:]):
    if 1 < b - a <= MAXGAP + 1 and all(x in present for x in range(a + 1, b)):
        filled.update(range(a + 1, b))

contacts = set()
for r in pep:
    ra = np.array([a.get_coord() for a in r])
    for rid, rra in ratoms.items():
        if np.min(np.linalg.norm(ra[:, None, :] - rra[None, :, :], axis=-1)) <= 4.5:
            contacts.add(rid)
lost = sorted(contacts - filled)
if lost:
    raise SystemExit("sealed crop drops GLP-1 contacts %s" % lost)


class Keep(Select):
    def accept_chain(self, c): return c.id == "R"
    def accept_residue(self, r): return r.id[0] == " " and r.id[1] in filled
    def accept_atom(self, a): return a.element != "H"


TARGET = "/kaggle/working/targets/GLP1R_sealed.pdb"
io = PDBIO(); io.set_structure(mdl); io.save(TARGET, Keep())
kept_res = [r for r in rec if r.id[1] in filled]
print("sealed target: %d residues, all %d GLP-1 contacts retained"
      % (len(kept_res), len(contacts)), flush=True)

# Prove the floor is actually sealed, against the intact receptor, in the output itself.
def band(residues, lo, hi):
    xyz = np.array([a.get_coord() for r in residues for a in r])
    d = xyz @ AXIS - OFF
    return int(((d >= lo) & (d < hi)).sum())


print("\ndepth band    sealed target    intact receptor    kept")
seal_ok = True
for lo, hi in ((33.2, 40), (40, 45), (45, 50), (50, 60)):
    nc, nf = band(kept_res, lo, hi), band(rec, lo, hi)
    pct = 100.0 * nc / nf if nf else 0.0
    print("%5.1f-%4.1f A %13d %18d %7.0f%%" % (lo, hi, nc, nf, pct))
    if hi <= 50 and pct < 70:
        seal_ok = False
if not seal_ok:
    raise SystemExit("the floor is still open below 50 A - widen SEAL_R")

native_depth = {r.id[1]: round(depth(ca(r)), 2) for r in pep if ca(r) is not None}
cmap = [dict(uniprot=r.id[1], aa=T31.get(r.get_resname(), "X"),
             depth=round(depth(ca(r)), 2) if ca(r) is not None else None,
             region="ecd" if r.id[1] <= ECD_END else "tm",
             contact=r.id[1] in contacts)
        for r in kept_res]
TM_C = sorted(x for x in contacts if x > ECD_END)
ECD_C_ = sorted(x for x in contacts if x <= ECD_END)
json.dump(dict(order=cmap, ecd_end=ECD_END, radius=RADIUS, seal_r=SEAL_R,
               ecd_centroid=[float(x) for x in ecd_c], axis=[float(x) for x in AXIS],
               native_glp1=dict(seq=native_seq, start=pep[0].id[1], depth=native_depth),
               tm_contacts=TM_C, ecd_contacts=ECD_C_),
          open("/kaggle/working/glp1r_sealed_map.json", "w"), indent=1)

TM_HOT = ",".join("R%d" % x for x in TM_C)
ECD_HOT = ",".join("R%d" % x for x in ECD_C_)
BOTH_HOT = ",".join("R%d" % x for x in sorted(contacts))
open("/kaggle/working/hot_tm.txt", "w").write(TM_HOT)
open("/kaggle/working/hot_ecd.txt", "w").write(ECD_HOT)
open("/kaggle/working/hot_both.txt", "w").write(BOTH_HOT)
print("\nhotspots: %d activation, %d binding, %d combined"
      % (len(TM_C), len(ECD_C_), len(contacts)), flush=True)

PROBE = ("python -c \"import jax, jax.lib; assert hasattr(jax.lib,'xla_bridge'); "
         "print(jax.__version__)\"")
r = subprocess.run(PROBE, shell=True, capture_output=True, text=True)
if r.returncode != 0 or not r.stdout.strip().startswith("0.4.38"):
    print("installing jax 0.4.38 ...", flush=True)
    sh("pip install -q 'jax[cuda12]==0.4.38' 'jaxlib==0.4.38'")
sh("pip install -q git+https://github.com/sokrypton/ColabDesign.git")

# ---------------------------------------------------------------- THE GATE
print("\n==== VALIDATING THE SEALED TARGET BEFORE SPENDING A NIGHT ON IT ====",
      flush=True)
from colabdesign import mk_afdesign_model, clear_mem

ORDER = cmap
WANT_TM = set(TM_C)


def dedupe(chain):
    res, seen = [], set()
    for r in list(chain):
        if r.id[0] != " " or r.id[1] in seen:
            chain.detach_child(r.id); continue
        seen.add(r.id[1]); res.append(r)
    return res


def check(seq, tag):
    clear_mem()
    af = mk_afdesign_model(protocol="binder", data_dir="/kaggle/working")
    af.prep_inputs(pdb_filename=TARGET, chain="R", binder_len=len(seq), hotspot=TM_HOT)
    try:
        af.predict(seq=seq, num_recycles=3, verbose=False)
    except TypeError:
        af.predict(seq=seq, verbose=False)
    lg = af.aux["log"]
    out = "/kaggle/working/gate_%s.pdb" % tag
    af.save_pdb(out)
    m = PDBParser(QUIET=True).get_structure("c", out)[0]
    kept = {ch.id: dedupe(ch) for ch in m}
    rid = next(c for c, v in kept.items() if len(v) == len(ORDER))
    bid = [c for c in kept if c != rid][0]
    ba = np.array([a.get_coord() for r in kept[bid] for a in r])
    proj = ba @ AXIS - OFF
    nt = np.array([a.get_coord() for a in kept[bid][0]]) @ AXIS - OFF
    on = 0
    for k, rr in enumerate(kept[rid]):
        ra = np.array([a.get_coord() for a in rr])
        if np.min(np.linalg.norm(ra[:, None, :] - ba[None, :, :], axis=-1)) <= 4.5:
            if ORDER[k]["uniprot"] in WANT_TM:
                on += 1
    return dict(tag=tag, i_ptm=round(float(lg.get("i_ptm", 0.0)), 3),
                plddt=round(float(lg.get("plddt", 0.0)), 3),
                deepest=round(float(np.max(proj)), 1),
                nterm=round(float(np.max(nt)), 1), on_site=on)


g = [check(native_seq, "native"), check(native_seq[2:], "trunc")]
for r in g:
    print("%-8s i_ptm %.3f pLDDT %.3f deepest %5.1f A nterm %5.1f A on_site %d"
          % (r["tag"], r["i_ptm"], r["plddt"], r["deepest"], r["nterm"], r["on_site"]),
          flush=True)
json.dump(g, open("/kaggle/working/glp1r_gate.json", "w"), indent=1)

nat = g[0]
REF = max(native_depth.values())
fail = []
if nat["i_ptm"] < 0.45:
    fail.append("native hormone i_ptm %.3f < 0.45 - the sealed target does not dock it"
                % nat["i_ptm"])
if not (25.0 <= nat["deepest"] <= 45.0):
    fail.append("native hormone reaches %.1f A, expected 25-45 (crystal says %.1f)"
                % (nat["deepest"], REF))
if nat["nterm"] < nat["deepest"] - 5.0:
    fail.append("native hormone went in backwards on the sealed target "
                "(N-term %.1f vs deepest %.1f)" % (nat["nterm"], nat["deepest"]))
if nat["on_site"] < 8:
    fail.append("native hormone touches only %d of its own %d bundle contacts"
                % (nat["on_site"], len(WANT_TM)))

if fail:
    for f in fail:
        print("GATE FAILED: %s" % f, flush=True)
    raise SystemExit("Sealing the target broke the docking it was meant to protect. "
                     "The generation cell will not run. Submit the landing-pad binders.")

print("\nGATE PASSED. The sealed target still docks the natural hormone correctly, and "
      "the floor below it is solid. Generation may proceed.", flush=True)

# Hand the cards back before the workers start, and SAY what is actually free on each.
# Three attempts were spent inferring a memory problem from design failures; one line of
# nvidia-smi here would have shown it in the first minute.
clear_mem()
import gc
gc.collect()
print("\nGPU state handed to the workers (index, used, free):")
print(subprocess.run("nvidia-smi --query-gpu=index,memory.used,memory.free "
                     "--format=csv,noheader", shell=True,
                     capture_output=True, text=True).stdout.strip(), flush=True)
print("SEALED SETUP DONE in %.1f min" % ((time.time() - t0) / 60))
