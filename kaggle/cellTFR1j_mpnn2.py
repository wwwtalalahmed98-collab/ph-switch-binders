"""Kaggle - TfR1 Challenge 4. Rewrite the sequences of run 2's WELL-AIMED backbones.

RUNS ON CPU. No GPU quota. Attach the run 2 notebook's output (notebook2e4cbf5468).

WHY THESE BACKBONES AND NOT YESTERDAY'S. Yesterday's rewrite ran on run 1's shapes, which
were well folded and badly aimed. Run 2 inverted that: at the e559 site it reaches 4+ epitope
contacts in 41% of designs against 11% before, confirmed on two platforms. Those shapes have
never been rewritten. Their sequences are the worst we have - runs of up to 22 identical
residues - which is exactly the defect a sequence rewrite repairs.

THE PROBLEM THIS IS TRYING TO BREAK. Across 102 designs from three runs, every sequence with
plausible composition makes ZERO epitope contacts, and every design that reaches the epitope
fails composition. The two lists have never overlapped. Aim and foldability have been trading
against each other, and this is the step that is supposed to get both at once: take a shape
that is already aimed, and put a real sequence on it.

WHAT IS DIFFERENT FROM THE FIRST MPNN PASS
  - run 2's designs carry on_site already (measured in-process), so no remeasure is needed
  - more sequences per backbone, at two temperatures, because this is the last attempt and
    diversity is cheap on a processor
  - selection is on COMPOSITION, two-sided, not on MPNN's own score. MPNN's score rewards
    what MPNN likes, which on a poor backbone is a charged repeat: the first pass improved
    all 39 sequences by that score and produced sequences that were 83% charged.

SCORE DIRECTION. MPNN's score is an average negative log-likelihood, so LOWER IS BETTER. It
is reported, but it does not drive selection.

WHAT THIS DOES NOT DO. It does not check that the rewritten sequence still folds or still
binds. That needs an AlphaFold pass per candidate and is the next step, not this one.
"""
import glob
import json
import os

import numpy as np
from Bio import PDB

W = "/kaggle/working"
ROOTS = os.environ.get("TFR1_ROOT", "/kaggle/input|/kaggle/working").split("|")
MIN_ON_SITE = 3          # 3+ epitope contacts before a backbone is worth rewriting
N_SEQ = 24               # per backbone per temperature
TEMPERATURES = (0.1, 0.2)
KEEP_PER_BACKBONE = 6

CHARGED = set("DEKR")
HYDROPHOBIC = set("AVLIMFWYC")
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5,
          L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)
# two-sided, same limits as kaggle/cellTFR1h_audit.py. The one-sided version of this scored
# a 64%-hydrophobic aggregate with AAAAAAAA as "good".
TIERS = [
    ("good", dict(same_run=5, charge_run=8, charged=(0.18, 0.40),
                  hydrophobic=(0.25, 0.42), kd=(-1.5, 0.0))),
    ("marginal", dict(same_run=6, charge_run=14, charged=(0.10, 0.55),
                      hydrophobic=(0.20, 0.50), kd=(-2.2, 0.3))),
]
T31 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
       "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
       "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def find(name):
    hits = [h for root in ROOTS for h in glob.glob("%s/**/%s" % (root, name), recursive=True)]
    return hits[0] if hits else None


def residues(chain):
    out, seen = [], set()
    for r in chain:
        if r.id[0] == " " and r.get_resname() in T31 and r.id[1] not in seen:
            seen.add(r.id[1])
            out.append(r)
    return out


def same_run(s):
    best = cur = 1
    for i in range(1, len(s)):
        cur = cur + 1 if s[i] == s[i - 1] else 1
        best = max(best, cur)
    return best


def class_run(s, group):
    best = cur = 0
    for c in s:
        cur = cur + 1 if c in group else 0
        best = max(best, cur)
    return best


def composition(seq):
    n = float(len(seq))
    return dict(charged_frac=round(sum(c in CHARGED for c in seq) / n, 3),
                hydrophobic_frac=round(sum(c in HYDROPHOBIC for c in seq) / n, 3),
                charge_run=class_run(seq, CHARGED), same_run=same_run(seq),
                kd=round(sum(KD[c] for c in seq) / n, 3))


def _fails(c, lim):
    why = []
    if c["same_run"] > lim["same_run"]:
        why.append("%d identical" % c["same_run"])
    if c["charge_run"] > lim["charge_run"]:
        why.append("%d charged in a row" % c["charge_run"])
    lo, hi = lim["charged"]
    if c["charged_frac"] > hi:
        why.append("%.0f%% charged" % (100 * c["charged_frac"]))
    elif c["charged_frac"] < lo:
        why.append("only %.0f%% charged" % (100 * c["charged_frac"]))
    lo, hi = lim["hydrophobic"]
    if c["hydrophobic_frac"] > hi:
        why.append("%.0f%% hydrophobic" % (100 * c["hydrophobic_frac"]))
    elif c["hydrophobic_frac"] < lo:
        why.append("only %.0f%% hydrophobic" % (100 * c["hydrophobic_frac"]))
    lo, hi = lim["kd"]
    if c["kd"] > hi:
        why.append("kd %+.2f net hydrophobic" % c["kd"])
    elif c["kd"] < lo:
        why.append("kd %+.2f extremely polar" % c["kd"])
    return why


def tier_of(c):
    for tier, lim in TIERS:
        if not _fails(c, lim):
            return tier, ""
    return "fail", "; ".join(_fails(c, TIERS[-1][1]))


def distance(c):
    """Distance from natural composition: 25% charged, 32% hydrophobic, kd about -0.6."""
    return round(abs(c["charged_frac"] - 0.25) / 0.25
                 + abs(c["hydrophobic_frac"] - 0.32) / 0.32
                 + abs(c["kd"] + 0.6) / 1.0
                 + max(0, c["same_run"] - 5) / 5.0
                 + max(0, c["charge_run"] - 8) / 8.0, 3)


src = find("tfr1_designs2.json")
if src is None:
    raise SystemExit("no tfr1_designs2.json - attach the run 2 notebook output")
DESIGNS = json.load(open(src))
TARGETS = json.load(open(find("tfr1_run2_targets.json")))
print("%d designs from run 2, sites %s" % (len(DESIGNS), sorted(TARGETS)))

CROP = {}
for tag, M in sorted(TARGETS.items()):
    path = find(os.path.basename(M["target"]))
    rs = residues(list(PDB.PDBParser(QUIET=True).get_structure("t", path)[0])[0])
    CROP[tag] = [r.get_resname() for r in rs]
    if len(CROP[tag]) != len(M["kept"]):
        raise SystemExit("%s: crop %d vs kept %d" % (tag, len(CROP[tag]), len(M["kept"])))

cands = sorted([d for d in DESIGNS if d.get("on_site", 0) >= MIN_ON_SITE],
               key=lambda d: (d["on_site"], -d.get("off_site", 0)), reverse=True)
print("%d of %d backbones reach %d+ epitope contacts\n"
      % (len(cands), len(DESIGNS), MIN_ON_SITE))
if not cands:
    raise SystemExit("no backbone is on the epitope - nothing worth rewriting")

from colabdesign.mpnn import mk_mpnn_model                      # noqa: E402

mpnn = mk_mpnn_model(model_name="v_48_020")
print("ProteinMPNN v_48_020 loaded\n")

records = []
for d in cands:
    pdb_path = find("%s.pdb" % d["name"])
    if pdb_path is None:
        print("%s: no structure found, skipped" % d["name"])
        continue
    model = PDB.PDBParser(QUIET=True).get_structure("c", pdb_path)[0]
    chains = {c.id: residues(c) for c in model}
    tid = next((cid for cid, rs in chains.items()
                if [r.get_resname() for r in rs] == CROP[d["site"]]), None)
    if tid is None:
        print("%s: no chain reproduces the crop, skipped" % d["name"])
        continue
    bid = [c for c in chains if c != tid]
    if len(bid) != 1:
        print("%s: expected one binder chain, got %s" % (d["name"], bid))
        continue
    bid = bid[0]

    try:
        mpnn.prep_inputs(pdb_filename=pdb_path, chain="%s,%s" % (tid, bid),
                         fix_pos=tid, rm_aa="C")
    except Exception as e:
        print("%s: prep failed, %s %s" % (d["name"], type(e).__name__, str(e)[:120]))
        continue

    orig = d["sequence"]
    oc = composition(orig)
    rows = []
    for temp in TEMPERATURES:
        out = mpnn.sample(num=N_SEQ, batch=1, temperature=temp)
        seqs = [s.split("/")[-1] if "/" in s else s
                for s in np.asarray(out["seq"]).reshape(-1)]
        scores = np.asarray(out["score"]).reshape(-1)
        for s, sc in zip(seqs, scores):
            if len(s) != len(orig):
                continue
            c = composition(s)
            tier, why = tier_of(c)
            rows.append(dict(sequence=s, mpnn_score=round(float(sc), 4), temperature=temp,
                             tier=tier, why=why, distance=distance(c), **c))

    # selection is on composition, not on MPNN's own score
    rows.sort(key=lambda r: (r["tier"] == "fail", r["distance"], r["mpnn_score"]))
    keep = [r for r in rows if r["tier"] != "fail"][:KEEP_PER_BACKBONE]
    tiers = {t: sum(1 for r in rows if r["tier"] == t) for t in ("good", "marginal", "fail")}

    print("%-20s on-site %d/%d off %-3d | original d=%.2f [%s] | %d sampled: "
          "%d good %d marginal %d fail | best d=%s"
          % (d["name"], d["on_site"], d.get("n_epitope", 0), d.get("off_site", 0),
             distance(oc), tier_of(oc)[0], len(rows), tiers["good"], tiers["marginal"],
             tiers["fail"], ("%.2f" % keep[0]["distance"]) if keep else "none"))

    for i, r in enumerate(keep, 1):
        r.update(name="%s_r%02d" % (d["name"], i), parent=d["name"], site=d["site"],
                 length=len(r["sequence"]), on_site_parent=d["on_site"],
                 off_site_parent=d.get("off_site"), parent_iptm=d.get("true_iptm"),
                 parent_plddt=d.get("plddt"), parent_sequence=orig,
                 parent_distance=distance(oc), structure=os.path.basename(pdb_path))
        records.append(r)

outpath = "%s/tfr1_mpnn2.json" % (W if os.path.isdir(W) else ".")
json.dump(records, open(outpath, "w"), indent=0)
print("\nwrote %s" % outpath)

print("\n==================== REWRITE RESULT ====================")
print("Selection is on composition (two-sided), not on MPNN's own score.")
print("Reference: 25%% charged, 32%% hydrophobic, Kyte-Doolittle about -0.6.\n")
if records:
    good = [r for r in records if r["tier"] == "good"]
    print("%d candidates from %d backbones: %d good, %d marginal"
          % (len(records), len(set(r["parent"] for r in records)), len(good),
             len(records) - len(good)))
    print("\nTHE QUESTION THIS RUN EXISTS TO ANSWER: is any design now BOTH on the epitope")
    print("and plausible as a protein? Previous answer across 102 designs: none.\n")
    print("%-26s %5s %6s %6s %6s %6s %5s %9s %6s"
          % ("candidate", "on", "chg", "phob", "kd", "dist", "mpnn", "tier", "T"))
    for r in sorted(records, key=lambda r: (r["tier"] == "marginal", r["distance"])):
        print("%-26s %5d %5.0f%% %5.0f%% %+6.2f %6.2f %5.2f %9s %6.1f"
              % (r["name"], r["on_site_parent"], 100 * r["charged_frac"],
                 100 * r["hydrophobic_frac"], r["kd"], r["distance"], r["mpnn_score"],
                 r["tier"], r["temperature"]))
    if good:
        print("\n%d sequences are BOTH on a well-aimed backbone and composition-good." % len(good))
        print("That is the first time in this project. Still unverified: none has been")
        print("folded, so none is yet shown to keep the binding its backbone had.")
    else:
        print("\nNone reaches 'good'. The marginal ones are the whole pool.")
else:
    print("nothing survived composition on any backbone - the shapes are likely not proteins")

print("\n===== COPY FROM HERE =====")
for r in records:
    print(json.dumps(r))
print("===== COPY TO HERE =====")
