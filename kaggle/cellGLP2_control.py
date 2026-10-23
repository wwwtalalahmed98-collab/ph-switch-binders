# === CHALLENGE 3 (GLP1R) : CAN ALPHAFOLD SEE THE AGONIST MECHANISM AT ALL? ===
#
# Before designing anything, test whether the method can distinguish an agonist from an
# antagonist. Nature has already run the experiment: GLP-1(7-36) activates the receptor,
# and the same peptide missing its first two residues is the shape of exendin 9-39, a
# real drug that binds the extracellular groove and BLOCKS activation.
#
# So predict both against the cropped receptor and measure how deep each one reaches.
#   - If 7-36 threads deep into the bundle and 9-36 sits shallower, the depth proxy has
#     support and designing agonists is worth the window.
#   - If both land in the same place, AlphaFold is blind to the thing the whole challenge
#     turns on, and we should fall back to extracellular binders and say so plainly.
#
# This is deliberately done FIRST. In the TNF challenge we selected on a geometric proxy
# nobody had validated and found out afterwards that it correlated +0.24 with what we
# actually wanted. The cost of finding that out early is one short run.
import json, os, random
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1 as T31

M = json.load(open("/kaggle/working/glp1r_map.json"))
TARGET = "/kaggle/working/targets/GLP1R_crop.pdb"
ECD_END = M["ecd_end"]
ECD_C = np.array(M["ecd_centroid"]); AXIS = np.array(M["axis"])
TM_HOT = open("/kaggle/working/glp1r_tm_hotspot.txt").read().strip()
NATIVE = M["native_glp1"]["seq"]          # GLP-1(7-36), 30 residues
NAT_DEPTH = {int(k): v for k, v in M["native_glp1"]["depth"].items()}
REF_DEEP = max(NAT_DEPTH.values())
print("native GLP-1(7-36): %s" % NATIVE)
print("its deepest residue sits at %.1f A; that is the bar an agonist has to clear\n"
      % REF_DEEP, flush=True)

random.seed(0)
scram = list(NATIVE); random.shuffle(scram); scram = "".join(scram)
CASES = [
    ("GLP-1(7-36) native agonist",      NATIVE,      "should reach deep"),
    ("GLP-1(9-36) antagonist-shaped",   NATIVE[2:],  "should sit shallower"),
    ("scrambled GLP-1, same residues",  scram,       "should bind poorly"),
    ("poly-alanine, same length",       "A" * 30,    "floor control"),
]

from colabdesign import mk_afdesign_model, clear_mem
_c = {}
def predict(seq, out_pdb):
    if _c.get("L") != len(seq):
        clear_mem()
        af = mk_afdesign_model(protocol="binder", data_dir="/kaggle/working")
        af.prep_inputs(pdb_filename=TARGET, chain="R", binder_len=len(seq), hotspot=TM_HOT)
        _c["L"], _c["af"] = len(seq), af
    af = _c["af"]
    try:
        af.predict(seq=seq, num_recycles=3, verbose=False)
    except TypeError:
        af.predict(seq=seq, verbose=False)
    lg = af.aux["log"]
    af.save_pdb(out_pdb)
    return (round(float(lg.get("i_ptm", lg.get("ptm", 0.0))), 3),
            round(float(lg.get("plddt", 0.0)), 3))

def dedupe(chain):
    res, seen = [], set()
    for r in list(chain):
        if r.id[0] != " " or r.id[1] in seen:
            chain.detach_child(r.id); continue
        seen.add(r.id[1]); res.append(r)
    return res

def measure(pdb_path, blen):
    """Depth of the binder's deepest atom, and how many receptor contacts it makes in
    the transmembrane bundle versus the extracellular domain."""
    mdl = PDBParser(QUIET=True).get_structure("c", pdb_path)[0]
    kept = {ch.id: dedupe(ch) for ch in mdl}
    bid = min(kept, key=lambda k: len(kept[k]))
    bres = kept[bid]
    if len(bres) != blen:
        raise ValueError("binder chain has %d residues, expected %d" % (len(bres), blen))
    tres = [r for cid, rs in kept.items() if cid != bid for r in rs]
    order = M["order"]
    if len(tres) != len(order):
        raise ValueError("target has %d residues, map has %d" % (len(tres), len(order)))
    batoms = np.array([a.get_coord() for r in bres for a in r])
    deep = float(np.max(batoms @ AXIS - float(ECD_C @ AXIS)))
    nterm = np.array([a.get_coord() for a in bres[0]])
    nterm_depth = float(np.max(nterm @ AXIS - float(ECD_C @ AXIS)))
    tm = ecd = 0
    for k, r in enumerate(tres):
        ra = np.array([a.get_coord() for a in r])
        if np.min(np.linalg.norm(ra[:, None, :] - batoms[None, :, :], axis=-1)) <= 4.5:
            if order[k]["region"] == "tm":
                tm += 1
            else:
                ecd += 1
    return deep, nterm_depth, tm, ecd

os.makedirs("/kaggle/working/ctrl_pdb", exist_ok=True)
rows = []
for label, seq, expect in CASES:
    out = "/kaggle/working/ctrl_pdb/%s.pdb" % label.split()[0].replace("(", "_").replace(")", "")
    try:
        ip, pl = predict(seq, out)
        deep, nd, tm, ecd = measure(out, len(seq))
        rows.append(dict(label=label, seq=seq, i_ptm=ip, plddt=pl, deepest=round(deep, 1),
                         nterm_depth=round(nd, 1), tm_contacts=tm, ecd_contacts=ecd))
        print("%-32s i_ptm %.3f pLDDT %.3f | deepest %5.1f A, N-term %5.1f A | "
              "%2d TM / %2d ECD contacts   (%s)"
              % (label, ip, pl, deep, nd, tm, ecd, expect), flush=True)
    except Exception as e:
        rows.append(dict(label=label, error="%s %s" % (type(e).__name__, str(e)[:160])))
        print("%-32s FAILED %s %s" % (label, type(e).__name__, str(e)[:160]), flush=True)

json.dump(rows, open("/kaggle/working/glp1r_control.json", "w"), indent=1)
ok = [r for r in rows if "deepest" in r]
if not ok:
    raise SystemExit("every control failed - this is not a verdict, read the errors above")

print("\n==== VERDICT ====")
ag = next((r for r in ok if r["label"].startswith("GLP-1(7")), None)
an = next((r for r in ok if r["label"].startswith("GLP-1(9")), None)
if ag and an:
    gap = ag["deepest"] - an["deepest"]
    tmgap = ag["tm_contacts"] - an["tm_contacts"]
    print("native agonist reaches %.1f A with %d transmembrane contacts"
          % (ag["deepest"], ag["tm_contacts"]))
    print("antagonist-shaped      %.1f A with %d"
          % (an["deepest"], an["tm_contacts"]))
    print("difference: %.1f A of depth, %d transmembrane contacts" % (gap, tmgap))
    print("native reference from the crystal structure: %.1f A" % REF_DEEP)
    recovered = ag["deepest"] >= REF_DEEP - 8
    separates = gap >= 4 or tmgap >= 4
    if recovered and separates:
        print("\nGO. AlphaFold recovers the deep binding mode AND separates the agonist")
        print("from the antagonist-shaped truncation. Depth is worth designing on.")
    elif recovered:
        print("\nPARTIAL. The deep mode is recovered but the truncation is not clearly")
        print("separated, so depth ranks candidates only weakly. Treat as a weak proxy.")
    else:
        print("\nNO GO. AlphaFold does not even reproduce the native deep binding mode")
        print("(%.1f A against a reference of %.1f A). Designing agonists on this signal"
              % (ag["deepest"], REF_DEEP))
        print("would be selecting on noise. Fall back to extracellular binders and say so.")
print("\n===== COPY FROM HERE =====")
for r in rows:
    print(json.dumps(r, separators=(",", ":")))
print("===== COPY TO HERE =====")
