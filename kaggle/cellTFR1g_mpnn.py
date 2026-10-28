"""Kaggle - TfR1 Challenge 4. Rewrite the sequences of the backbones that are ON TARGET.

RUNS ON CPU. No GPU quota. Attach BOTH notebook outputs as inputs:
    tfr1_generate   - the design structures
    tfr1_remeasure  - tfr1_designs_measured.json, which says which ones hit the epitope

WHY THIS EXISTS, and why it should have existed three challenges ago. Run 1's measurement
split the designs two ways and nothing landed in the middle:

    well-formed, wrong place   32 of 47 made zero epitope contacts, including all five
                               of the most confident (i_ptm 0.66-0.78)
    on target, malformed       the 3 that reached 4+ contacts score i_ptm 0.15-0.46 and
                               carry runs like SPPSPPPPPPPP and NNNNNNNNNSYNNTNN

The second group is not a lost cause: a homopolymer run is a SEQUENCE defect, and those
designs' BACKBONES are sitting on the epitope. ProteinMPNN reads a fixed backbone and writes
a sequence for it, which is exactly the repair needed. Hallucinate-then-MPNN is the standard
shape of a binder pipeline; this project has been doing only the first half.

WHAT IS FIXED AND WHAT IS DESIGNED. The target chain is fixed (fix_pos), so only the binder
is rewritten. The target chain is identified by reproducing the crop residue for residue, not
by its numbering - save_pdb renumbers from 1, which is what left run 1 unmeasured.

SCORE DIRECTION. ProteinMPNN's score is an average negative log-likelihood of the sequence
given the backbone, so LOWER IS BETTER. The original hallucinated sequence is scored on the
same backbone for comparison, which measures how bad it actually was rather than relying on
"that looks like junk to me".

WHAT THIS DOES NOT DO. It does not change where the binder sits, and it does not verify that
the rewritten sequence still folds and still binds. Every sequence out of here is a CANDIDATE
that needs one AlphaFold pass to confirm, which is GPU work and is not done here. MPNN
agreeing with a backbone is not evidence the complex is real.
"""
import glob
import json
import os

import numpy as np
from Bio import PDB

W = "/kaggle/working"
ROOTS = os.environ.get("TFR1_ROOT", "/kaggle/input|/kaggle/working").split("|")
MIN_ON_SITE = 2          # epitope contacts a backbone needs before it is worth rewriting
N_SEQ = 16               # sequences sampled per backbone
TEMPERATURE = 0.1        # low: take MPNN's confident choices rather than explore
KEEP_PER_BACKBONE = 4    # best-scoring sequences carried forward per backbone
MAX_RUN = 6              # homopolymer run we refuse to carry forward, as in challenge 3
T31 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
       "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
       "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5,
          L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)


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


def longest_run(s):
    best = cur = 1
    for i in range(1, len(s)):
        cur = cur + 1 if s[i] == s[i - 1] else 1
        best = max(best, cur)
    return best


def residues(chain):
    out, seen = [], set()
    for r in chain:
        if r.id[0] == " " and r.get_resname() in T31 and r.id[1] not in seen:
            seen.add(r.id[1])
            out.append(r)
    return out


measured = find("tfr1_designs_measured.json")
if measured is None:
    raise SystemExit("no tfr1_designs_measured.json - attach the tfr1_remeasure notebook's "
                     "output as well as tfr1_generate's. Without it there is no way to tell "
                     "which backbones are on the epitope, and rewriting the sequence of a "
                     "backbone stuck to the cut face just produces a tidier wrong answer.")
DESIGNS = json.load(open(measured))
TARGETS = json.load(open(find("tfr1_targets.json")))
print("%d designs measured, sites %s" % (len(DESIGNS), sorted(TARGETS)))

CROP = {}
for tag, M in sorted(TARGETS.items()):
    path = find(os.path.basename(M["target"]))
    rs = residues(list(PDB.PDBParser(QUIET=True).get_structure("t", path)[0])[0])
    CROP[tag] = [r.get_resname() for r in rs]
    if len(CROP[tag]) != len(M["kept"]):
        raise SystemExit("%s: crop PDB %d residues, kept list %d"
                         % (tag, len(CROP[tag]), len(M["kept"])))
print("crops loaded: %s" % {k: len(v) for k, v in sorted(CROP.items())})

cands = sorted([d for d in DESIGNS if d.get("on_site", 0) >= MIN_ON_SITE],
               key=lambda d: (d["on_site"], -d.get("off_site", 0)), reverse=True)
print("\n%d of %d backbones make at least %d epitope contacts"
      % (len(cands), len(DESIGNS), MIN_ON_SITE))
if not cands:
    raise SystemExit("no backbone is on the epitope - there is nothing worth rewriting, and "
                     "the problem is targeting, not sequence")

from colabdesign.mpnn import mk_mpnn_model                      # noqa: E402

mpnn = mk_mpnn_model(model_name="v_48_020")
print("ProteinMPNN v_48_020 loaded (ships inside ColabDesign, no download)\n")

records = []
for d in cands:
    hits = find_all("tfr1_%s/%s.pdb" % (d["site"], d["name"]), "%s.pdb" % d["name"])
    if not hits:
        print("%s: no structure found, skipped" % d["name"])
        continue
    pdb_path = hits[0]

    model = PDB.PDBParser(QUIET=True).get_structure("c", pdb_path)[0]
    chains = {c.id: residues(c) for c in model}
    tid = None
    for cid, rs in chains.items():
        if [r.get_resname() for r in rs] == CROP[d["site"]]:
            tid = cid
            break
    if tid is None:
        print("%s: no chain reproduces the crop, skipped" % d["name"])
        continue
    bid = [c for c in chains if c != tid]
    if len(bid) != 1:
        print("%s: expected one binder chain, got %s, skipped" % (d["name"], bid))
        continue
    bid = bid[0]

    try:
        # fix the whole target chain; a bare chain letter selects all of it and asserts if
        # the chain is absent, so a wrong id fails loudly rather than designing the target
        mpnn.prep_inputs(pdb_filename=pdb_path, chain="%s,%s" % (tid, bid),
                         fix_pos=tid, rm_aa="C")
    except Exception as e:
        print("%s: prep failed, %s %s" % (d["name"], type(e).__name__, str(e)[:140]))
        continue

    orig = d["sequence"]
    try:
        base = float(np.asarray(mpnn.score(seq=orig)["score"]).reshape(-1)[0])
    except Exception as e:
        base = float("nan")
        print("  (scoring the original failed: %s %s)" % (type(e).__name__, str(e)[:90]))

    out = mpnn.sample(num=N_SEQ, batch=1, temperature=TEMPERATURE)
    seqs = [s.split("/")[-1] if "/" in s else s for s in np.asarray(out["seq"]).reshape(-1)]
    scores = np.asarray(out["score"]).reshape(-1)
    seqids = np.asarray(out["seqid"]).reshape(-1)

    rows = []
    for s, sc, si in zip(seqs, scores, seqids):
        if len(s) != len(orig):
            continue                      # the binder chain is what we asked for; skip odd ones
        rows.append(dict(sequence=s, mpnn_score=round(float(sc), 4),
                         seqid_to_original=round(float(si), 3),
                         hydropathy=round(sum(KD[c] for c in s) / len(s), 3),
                         n_cys=s.count("C"), longest_run=longest_run(s)))
    rows.sort(key=lambda r: r["mpnn_score"])          # lower is better
    good = [r for r in rows if r["longest_run"] <= MAX_RUN and r["n_cys"] == 0]
    keep = good[:KEEP_PER_BACKBONE]

    print("%-18s on-site %d/%d  original score %.3f run %d  ->  %d of %d rewrites clean, "
          "best score %.3f run %d"
          % (d["name"], d["on_site"], d.get("n_epitope", 0), base, longest_run(orig),
             len(good), len(rows), keep[0]["mpnn_score"] if keep else float("nan"),
             keep[0]["longest_run"] if keep else -1))

    for i, r in enumerate(keep, 1):
        r.update(name="%s_mpnn%02d" % (d["name"], i), parent=d["name"], site=d["site"],
                 length=len(r["sequence"]), on_site_parent=d["on_site"],
                 off_site_parent=d.get("off_site"), parent_iptm=d.get("true_iptm"),
                 parent_sequence=orig, parent_mpnn_score=round(base, 4),
                 parent_longest_run=longest_run(orig), structure=os.path.basename(pdb_path))
        records.append(r)

outpath = "%s/tfr1_mpnn.json" % (W if os.path.isdir(W) else ".")
json.dump(records, open(outpath, "w"), indent=0)
print("\nwrote %s" % outpath)

print("\n==================== MPNN RESULT ====================")
print("ProteinMPNN score is average negative log-likelihood: LOWER IS BETTER.\n")
if records:
    improved = [r for r in records
                if r["parent_mpnn_score"] == r["parent_mpnn_score"]
                and r["mpnn_score"] < r["parent_mpnn_score"]]
    print("%d candidate sequences from %d backbones; %d score better than the hallucinated "
          "original on the same backbone"
          % (len(records), len(set(r["parent"] for r in records)), len(improved)))
    print("\n%-26s %6s %8s %8s %6s %5s %4s" % ("candidate", "on-site", "score", "was",
                                               "seqid", "kd", "run"))
    for r in sorted(records, key=lambda r: (-r["on_site_parent"], r["mpnn_score"])):
        print("%-26s %6d %8.3f %8.3f %6.2f %5.2f %4d"
              % (r["name"], r["on_site_parent"], r["mpnn_score"], r["parent_mpnn_score"],
                 r["seqid_to_original"], r["hydropathy"], r["longest_run"]))
    print("\nNONE OF THESE IS VERIFIED. A rewritten sequence still has to fold and still has")
    print("to bind; that needs one AlphaFold pass per candidate and has not been run.")
else:
    print("no clean rewrite survived the buildability filter (max run %d, no cysteine)"
          % MAX_RUN)

print("\n===== COPY FROM HERE =====")
for r in records:
    print(json.dumps(r))
print("===== COPY TO HERE =====")
