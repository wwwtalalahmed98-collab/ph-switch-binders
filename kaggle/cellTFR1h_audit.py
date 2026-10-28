"""Kaggle - TfR1 Challenge 4. Two checks we have never run, both free.

RUNS ON CPU. No GPU quota. Attach tfr1_generate, tfr1_remeasure, and the MPNN notebook's
output.

CHECK 1 - WHERE THE OFF-SITE CONTACTS ACTUALLY ARE
We have been counting off-site contacts as a single number, which hides a distinction that
matters. "Off-site" covers two very different things:

  the cut face   surface the crop invented. Harmless in itself - it does not exist in the
                 real protein - but a binder gripping it is binding nothing.
  the TRANSFERRIN / HFE site   the receptor's actual job. A binder that blocks it is not a
                 missed shot, it is a functional problem, and the whole point of three
                 revisions of epitope filtering was to stay clear of it.

We excluded that site from the epitope and then never checked whether any design touches it.
The footprint is the 45 residues contacted in 1SUV (transferrin) and 1DE4 (HFE), computed
earlier by tfr1/find_patches.py and embedded below so this needs no network.

CHECK 2 - COMPOSITION, BECAUSE THE BUILDABILITY GATE HAS A HOLE
The MPNN pass improved all 39 sequences by MPNN's own score and removed every homopolymer
run. Reading the output, some are not proteins:

    ITVKLPKPDKKKGPKKPKPKKKKKEKEKKEKKKKEEEKKKEEERKKKEEEERRRKKEEEE
    83% charged, 5% hydrophobic, a 41-residue unbroken run of charged residues

The existing gate looks for runs of the SAME residue, so K-E-K-K-E-R never trips it. A
charged run is just as unbuildable and does not fold. Natural globular proteins sit near 25%
charged and 30-35% hydrophobic.

This also says something about where the problem is: the MPNN rewrites that look like real
sequences come from the parents with the highest confidence, and the polyelectrolytes come
from the junk parents. MPNN cannot repair a backbone that is not a protein. That points at
the backbones, not the sequences, as the thing still to fix.

NOTHING HERE IS A PREDICTION. Both checks are arithmetic on structures and sequences we
already have.
"""
import glob
import json
import os
import urllib.request

import numpy as np
from Bio import PDB
from Bio.PDB.SASA import ShrakeRupley

W = "/kaggle/working"
ROOTS = os.environ.get("TFR1_ROOT", "/kaggle/input|/kaggle/working").split("|")
CONTACT = 4.5
RIM_SASA = 1.0           # A^2 of surface the crop created at a residue before it counts as rim

# Residues of human TfR1 contacted by transferrin (1SUV) or HFE (1DE4), from
# tfr1/find_patches.py. Embedded rather than recomputed so this needs no network and cannot
# drift from the list the epitope choice was actually made against.
TRANSFERRIN_SITE = {
    123, 124, 125, 126, 318, 527, 528, 529, 604, 615, 616, 619, 620, 622, 623, 624, 626,
    627, 629, 630, 638, 640, 641, 643, 644, 646, 648, 650, 651, 653, 654, 655, 657, 658,
    659, 661, 662, 663, 664, 665, 667, 669, 671, 674, 759,
}

CHARGED = set("DEKR")
HYDROPHOBIC = set("AVLIMFWYC")
# Graded rather than pass/fail, deliberately. A single strict cut rejected all 39 MPNN
# sequences including the most protein-like one, which is a verdict rather than a ranking
# and leaves nothing to choose between. Anchors: natural globular proteins sit near 25%
# charged and 30-35% hydrophobic; de novo helical bundles run more charged than natural for
# solubility, commonly 30-40%. "good" is near natural, "marginal" is designed-protein
# territory, "fail" is past anything that folds.
TIERS = [
    ("good",     dict(charge_run=8,  charged_frac=0.40, hydrophobic_frac=0.25)),
    ("marginal", dict(charge_run=14, charged_frac=0.55, hydrophobic_frac=0.18)),
]

T31 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
       "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
       "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def find(name):
    hits = [h for root in ROOTS for h in glob.glob("%s/**/%s" % (root, name), recursive=True)]
    return hits[0] if hits else None


def find_all(*patterns):
    for pat in patterns:
        hits = [h for root in ROOTS
                for h in glob.glob("%s/**/%s" % (root, pat), recursive=True)]
        if hits:
            return hits
    return []


def residues(chain):
    out, seen = [], set()
    for r in chain:
        if r.id[0] == " " and r.get_resname() in T31 and r.id[1] not in seen:
            seen.add(r.id[1])
            out.append(r)
    return out


def class_run(seq, group):
    best = cur = 0
    for c in seq:
        cur = cur + 1 if c in group else 0
        best = max(best, cur)
    return best


def composition(seq):
    n = float(len(seq))
    return dict(charged_frac=round(sum(c in CHARGED for c in seq) / n, 3),
                hydrophobic_frac=round(sum(c in HYDROPHOBIC for c in seq) / n, 3),
                charge_run=class_run(seq, CHARGED))


def composition_tier(c):
    """Returns (tier, why) where tier is good / marginal / fail."""
    for tier, lim in TIERS:
        if (c["charge_run"] <= lim["charge_run"]
                and c["charged_frac"] <= lim["charged_frac"]
                and c["hydrophobic_frac"] >= lim["hydrophobic_frac"]):
            return tier, ""
    loose = TIERS[-1][1]
    why = []
    if c["charge_run"] > loose["charge_run"]:
        why.append("%d charged in a row" % c["charge_run"])
    if c["charged_frac"] > loose["charged_frac"]:
        why.append("%.0f%% charged" % (100 * c["charged_frac"]))
    if c["hydrophobic_frac"] < loose["hydrophobic_frac"]:
        why.append("only %.0f%% hydrophobic" % (100 * c["hydrophobic_frac"]))
    return "fail", "; ".join(why)


# ---------------------------------------------------------------- inputs
measured = find("tfr1_designs_measured.json")
if measured is None:
    raise SystemExit("no tfr1_designs_measured.json - attach the tfr1_remeasure output")
DESIGNS = json.load(open(measured))
TARGETS = json.load(open(find("tfr1_targets.json")))

CIF = "%s/11QC.cif" % W
if not os.path.exists(CIF):
    hit = find("11QC.cif")
    if hit:
        CIF = hit
    else:
        urllib.request.urlretrieve("https://files.rcsb.org/download/11QC.cif", CIF)

full = PDB.MMCIFParser(QUIET=True).get_structure("x", CIF)[0]
REF = [c for c in full if sum(1 for r in c if r.id[0] == " ") > 100][0].id
ShrakeRupley().compute(full, level="R")
FULL_SASA = {r.id[1]: r.sasa for r in full[REF] if r.id[0] == " " and r.get_resname() in T31}

# which kept residues gained surface when the crop was cut: that is the rim
RIM, CROP = {}, {}
for tag, M in sorted(TARGETS.items()):
    path = find(os.path.basename(M["target"]))
    st = PDB.PDBParser(QUIET=True).get_structure("t", path)
    ShrakeRupley().compute(st[0], level="R")
    rs = residues(list(st[0])[0])
    CROP[tag] = [r.get_resname() for r in rs]
    RIM[tag] = {r.id[1] for r in rs
                if r.sasa - FULL_SASA.get(r.id[1], r.sasa) > RIM_SASA}
    ep = set(M["epitope"])
    print("%s: crop %d residues, %d on the rim, %d of the rim also epitope, "
          "%d of the crop in the transferrin site"
          % (tag, len(rs), len(RIM[tag]), len(RIM[tag] & ep),
             len(set(M["kept"]) & TRANSFERRIN_SITE)))

# ---------------------------------------------------------------- check 1
print("\n==== WHERE THE CONTACTS LAND ====")
print("%-18s %6s %6s %6s %6s  %s" % ("design", "epitope", "transf", "rim", "other", "flag"))
rows = []
for d in sorted(DESIGNS, key=lambda d: -d.get("on_site", 0)):
    hits = find_all("tfr1_%s/%s.pdb" % (d["site"], d["name"]), "%s.pdb" % d["name"])
    if not hits:
        continue
    M = TARGETS[d["site"]]
    kept, ep = M["kept"], set(M["epitope"])
    model = PDB.PDBParser(QUIET=True).get_structure("c", hits[0])[0]
    chains = {c.id: residues(c) for c in model}
    tid = next((cid for cid, rs in chains.items()
                if [r.get_resname() for r in rs] == CROP[d["site"]]), None)
    if tid is None:
        continue
    bid = [c for c in chains if c != tid]
    if len(bid) != 1:
        continue
    ba = np.array([a.get_coord() for r in chains[bid[0]] for a in r if a.element != "H"])
    n_ep = n_tf = n_rim = n_other = 0
    tf_hit = []
    for i, r in enumerate(chains[tid]):
        uni = kept[i]
        ra = np.array([a.get_coord() for a in r if a.element != "H"])
        if np.min(np.linalg.norm(ra[:, None, :] - ba[None, :, :], axis=-1)) > CONTACT:
            continue
        if uni in ep:
            n_ep += 1
        elif uni in TRANSFERRIN_SITE:
            n_tf += 1
            tf_hit.append(uni)
        elif uni in RIM[d["site"]]:
            n_rim += 1
        else:
            n_other += 1
    flag = "BLOCKS TRANSFERRIN SITE: %s" % tf_hit if tf_hit else ""
    rows.append(dict(name=d["name"], site=d["site"], epitope=n_ep, transferrin=n_tf,
                     rim=n_rim, other=n_other, transferrin_residues=tf_hit))
    print("%-18s %6d %6d %6d %6d  %s" % (d["name"], n_ep, n_tf, n_rim, n_other, flag))

tf_bad = [r for r in rows if r["transferrin"] > 0]
print("\n%d of %d designs touch the transferrin/HFE site" % (len(tf_bad), len(rows)))
if not tf_bad:
    print("  (good: the epitope filtering held - no design blocks the receptor's own job)")
rim_bound = [r for r in rows if r["rim"] > r["epitope"]]
print("%d of %d designs make more contacts on the cut face than on the epitope"
      % (len(rim_bound), len(rows)))

# ---------------------------------------------------------------- check 2
print("\n==== COMPOSITION ====")
print("natural globular proteins: about 25%% charged, 30-35%% hydrophobic\n")
mp = find("tfr1_mpnn.json")
pool = []
for d in DESIGNS:
    pool.append(("run1", d["name"], d["sequence"], d.get("on_site", 0)))
if mp:
    for r in json.load(open(mp)):
        pool.append(("mpnn", r["name"], r["sequence"], r.get("on_site_parent", 0)))
else:
    print("(no tfr1_mpnn.json attached - checking the original designs only)\n")

print("%-24s %6s %7s %8s %8s %9s  %s"
      % ("sequence", "source", "charged", "hydroph", "chg-run", "tier", "why"))
graded = []
for src, name, seq, on in sorted(pool, key=lambda x: -x[3]):
    c = composition(seq)
    tier, why = composition_tier(c)
    graded.append(dict(source=src, name=name, sequence=seq, on_site=on, tier=tier, **c))
    if src == "mpnn" or on >= 2:
        print("%-24s %6s %6.0f%% %7.0f%% %8d %9s  %s"
              % (name[:24], src, 100 * c["charged_frac"], 100 * c["hydrophobic_frac"],
                 c["charge_run"], tier, why))

print("\ncomposition tiers, by source:")
for src in ("run1", "mpnn"):
    sub = [g for g in graded if g["source"] == src]
    if sub:
        print("  %-5s %d sequences: %d good, %d marginal, %d fail"
              % (src, len(sub), sum(g["tier"] == "good" for g in sub),
                 sum(g["tier"] == "marginal" for g in sub),
                 sum(g["tier"] == "fail" for g in sub)))

clean = [g for g in graded if g["tier"] != "fail"]
print("\n%d of %d sequences are good or marginal" % (len(clean), len(graded)))
if not clean:
    print("  Nothing survives. These thresholds are GENEROUS against natural proteins at")
    print("  25% charged, so this says the backbones are the problem, not the sequences.")

print("\nleast-bad by composition, regardless of tier:")
for g in sorted(graded, key=lambda g: (g["charge_run"], g["charged_frac"]))[:8]:
    print("  %-24s %6s  %3.0f%% charged, %3.0f%% hydrophobic, run %2d  [%s]"
          % (g["name"][:24], g["source"], 100 * g["charged_frac"],
             100 * g["hydrophobic_frac"], g["charge_run"], g["tier"]))

out = dict(contacts=rows, clean=clean)
outpath = "%s/tfr1_audit.json" % (W if os.path.isdir(W) else ".")
json.dump(out, open(outpath, "w"), indent=0)
print("\nwrote %s" % outpath)

print("\n===== COPY FROM HERE =====")
for c in clean:
    print(json.dumps(c))
print("===== COPY TO HERE =====")
