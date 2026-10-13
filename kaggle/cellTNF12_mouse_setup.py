# === CHALLENGE 2 (TNF-alpha) : MOUSE TARGET SETUP ===
#
# Objective 2 is mouse cross-reactivity, and it is ranked above affinity. Until now it
# has only been inferred from sequence conservation across each design's human
# epitope. This builds the mouse equivalent of the human target so binding can be
# predicted directly instead.
#
# PDB 2TNF is mouse TNF-alpha at 1.4 A (UniProt P06804, chains A/B/C, a trimer) -
# higher resolution than the human 1TNF used for design.
#
# The numbering is the part that has to be right. Human 1TNF is numbered UniProt-76;
# mouse 2TNF is numbered UniProt-79. Worse, the human->mouse alignment is NOT a
# constant shift: across the 25 site residues the offset is +2 at some positions and
# +3 at others, so there is an indel between the orthologs. Every mouse residue here
# therefore comes from the pairwise alignment, never from arithmetic on a human
# number. On EGFR a three-residue numbering shift silently emptied a result list for
# two rounds, so each step below asserts what it believes.
#
# CORRECTED after the first attempt failed on its own assertion with
# "mouse 168 is GLU, expected LYS". 2TNF's author numbering is NOT a constant offset
# from UniProt: the chain runs 9-157, which is 149 positions, while its DBREF claims
# 1-156, and the observed offset takes BOTH the values 78 and 79 within a single
# chain. Mouse K168 is residue 90, not 89. Residue numbers are therefore derived by
# aligning each observed chain sequence to P06804 and are never computed by adding a
# constant - the same rule already applied to the human construct, which this cell
# originally failed to follow.
import glob, json, os, shutil, subprocess, sys, time
t0 = time.time()
RADIUS = 22.0                  # same sphere as the human crop

# UniProt P06804, mouse TNF-alpha precursor, 235 aa. Embedded so the mapping does not
# depend on a second network fetch at run time.
P06804 = (
    "MSTESMIRDVELAEEALPQKMGGFQNSRRCLCLSLFSFLLVAGATTLFCLLNFGVIGPQR"
    "DEKFPNGLPLISSMAQTLTLRSSSQNSSDKPVAHVVANHQVEEQLEWLSQRANALLANGM"
    "DLKDNQLVVPADGLYLVYSQVLFKGQGCPDYVLLTHTVSRFAISYQEKVNLLSAVKSPCP"
    "KDTPEGAELKPWYEPIYLGGVFQLEKGDQLSAEVNLPKYLDFAESGQVYFGVIAL")

# From the pairwise alignment of P01375 and P06804 (tnf/results/human_to_mouse_index.json).
# Mouse UniProt numbers for the 25 human site residues, in human order.
MOUSE_SITE = [98, 108, 109, 110, 111, 112, 113, 114, 115, 157, 161, 163, 164, 165,
              166, 167, 168, 169, 170, 203, 205, 209, 222, 223, 225]
HUMAN_SITE = [95, 105, 106, 107, 108, 109, 110, 111, 112, 155, 159, 161, 162, 163,
              164, 165, 166, 167, 168, 201, 203, 207, 220, 221, 223]
ANCHOR_K, ANCHOR_R = 168, 111  # mouse equivalents of human K166 and R108


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
raw = "/kaggle/working/2TNF.pdb"
if not os.path.exists(raw):
    open(raw, "wb").write(urllib.request.urlopen(
        "https://files.rcsb.org/download/2TNF.pdb", timeout=180).read())

from Bio import Align
_al = Align.PairwiseAligner()
_al.open_gap_score, _al.extend_gap_score = -11, -1
_al.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")

P = PDBParser(QUIET=True)
mdl = P.get_structure("m", raw)[0]

# residue number -> UniProt index, per chain, by ALIGNMENT not arithmetic
R2U = {}
for cid in ("A", "B"):
    res = [r for r in mdl[cid] if r.id[0] == " "]
    seq = "".join(protein_letters_3to1.get(r.get_resname(), "X") for r in res)
    aln = _al.align(P06804, seq)[0]
    mp = {}
    for (us, ue), (ss, se) in zip(aln.aligned[0], aln.aligned[1]):
        for k in range(int(ue - us)):
            mp[res[int(ss + k)].id[1]] = int(us + k) + 1
    bad = [(n, u) for n, u in mp.items()
           if P06804[u - 1] != protein_letters_3to1.get(mdl[cid][n].get_resname(), "X")]
    assert not bad, "chain %s: %d residues disagree with P06804, e.g. %s" % (cid, len(bad), bad[:3])
    R2U[cid] = mp
    print("chain %s: %d residues mapped to P06804, 0 mismatches, offsets present %s"
          % (cid, len(mp), sorted({u - n for n, u in mp.items()})), flush=True)

U2R = {cid: {u: n for n, u in mp.items()} for cid, mp in R2U.items()}
ca, name3 = {}, {}
for ch in mdl:
    for r in ch:
        if r.id[0] == " " and "CA" in r:
            ca[(ch.id, r.id[1])] = r["CA"].get_coord()
            name3[(ch.id, r.id[1])] = r.get_resname()

assert ANCHOR_K in U2R["A"], "mouse K%d is absent from 2TNF chain A" % ANCHOR_K
assert ANCHOR_R in U2R["B"], "mouse R%d is absent from 2TNF chain B" % ANCHOR_R
kA, rB = ("A", U2R["A"][ANCHOR_K]), ("B", U2R["B"][ANCHOR_R])
assert name3[kA] == "LYS", "mouse %d (resi A%d) is %s, expected LYS" % (ANCHOR_K, kA[1], name3[kA])
assert name3[rB] == "ARG", "mouse %d (resi B%d) is %s, expected ARG" % (ANCHOR_R, rB[1], name3[rB])
print("anchors confirmed in 2TNF: K%d = resi A%d, R%d = resi B%d"
      % (ANCHOR_K, kA[1], ANCHOR_R, rB[1]), flush=True)

centre = (ca[kA] + ca[rB]) / 2
keep = {k for k, v in ca.items() if k[0] in ("A", "B")
        and float(np.linalg.norm(v - centre)) <= RADIUS}


class Crop(Select):
    def accept_chain(self, c): return c.id in ("A", "B")
    def accept_residue(self, r): return r.id[0] == " " and (r.parent.id, r.id[1]) in keep
    def accept_atom(self, a): return a.element != "H"


TARGET = "/kaggle/working/targets/TNF_mouse_site_crop.pdb"
st = P.get_structure("t", raw)
io = PDBIO(); io.set_structure(st); io.save(TARGET, Crop())

order = []
for ch in P.get_structure("c", TARGET)[0]:
    for r in ch:
        if r.id[0] == " ":
            order.append(dict(chain=ch.id, pdb=r.id[1], uniprot=R2U[ch.id][r.id[1]],
                              aa=protein_letters_3to1.get(r.get_resname(), "X")))
H2M = dict(zip(HUMAN_SITE, MOUSE_SITE))
M2H = {m: h for h, m in H2M.items()}
json.dump(dict(species="mouse", pdb="2TNF", radius=RADIUS, n=len(order),
               numbering="per-residue, from alignment to P06804; offset is not constant",
               order=order, mouse_site=MOUSE_SITE, human_site=HUMAN_SITE,
               mouse_to_human_site=M2H),
          open("/kaggle/working/construct_map_mouse.json", "w"))

have = {o["uniprot"] for o in order}
print("cropped mouse target: %d residues (A %d, B %d)"
      % (len(order), sum(1 for o in order if o["chain"] == "A"),
         sum(1 for o in order if o["chain"] == "B")), flush=True)
print("mouse site residues retained: %d of %d" % (len(have & set(MOUSE_SITE)),
                                                  len(MOUSE_SITE)), flush=True)
missing = sorted(set(MOUSE_SITE) - have)
assert not missing, "mouse site residues cropped away: %s" % missing
assert ("A", ANCHOR_K) in {(o["chain"], o["uniprot"]) for o in order}
assert ("B", ANCHOR_R) in {(o["chain"], o["uniprot"]) for o in order}
print("switch anchors present: K%d(A), R%d(B)" % (ANCHOR_K, ANCHOR_R), flush=True)

# Sanity check the alignment itself: every site residue identity must match what the
# human/mouse alignment predicted, or the mapping is wrong somewhere.
seqmap = {o["uniprot"]: o["aa"] for o in order}
HUMAN_AA = dict(zip(HUMAN_SITE, "NLNRRANALTIVSYQTKVNQERFAS"))
bad = [(m, seqmap[m]) for m in MOUSE_SITE if m in seqmap and seqmap[m] not in "ACDEFGHIKLMNPQRSTVWY"]
assert not bad, "unparsable residues %s" % bad
print("mouse site sequence: %s" % "".join(seqmap[m] for m in MOUSE_SITE), flush=True)
print("human site sequence: %s" % "".join(HUMAN_AA[h] for h in HUMAN_SITE), flush=True)

HOTSPOT = ",".join(
    sorted({"%s%d" % (o["chain"], o["pdb"]) for o in order
            if o["uniprot"] in set(MOUSE_SITE)}))
open("/kaggle/working/mouse_hotspot.txt", "w").write(HOTSPOT)
print("mouse hotspot string: %s" % HOTSPOT, flush=True)

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
    "m.prep_inputs(pdb_filename='%s', chain='A,B', binder_len=65); "
    "print('READY', jax.__version__, jax.devices()); "
    "print('mouse target length seen by the model:', m._target_len)\"" % TARGET,
    shell=True, capture_output=True, text=True)
print(chk.stdout[-500:])
if chk.returncode != 0:
    print("PREP FAILED:"); print(chk.stderr[-1500:])
print("MOUSE TARGET SETUP DONE in %.1f min" % ((time.time() - t0) / 60))
