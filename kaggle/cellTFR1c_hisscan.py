"""Kaggle - TfR1 Challenge 4, pH-switch pass. Fold histidine mutants and score the switch.

RUN AFTER the generation notebook. Attach its output as an input dataset; this script finds
it under /kaggle/input and will say so if it cannot.

WHY A SEPARATE PASS. Designing a binder and designing a pH switch at once gave neither on
TNF-alpha. Generating a good binder first and then scanning histidine placements onto it
gave 35 measured placements across 14 lineages, which is where the useful numbers came from.
This repeats that, with one change: the switch is now scored by a thermodynamic calculation
that was validated against complexes whose pH behaviour is known, instead of a geometric
proxy. The first proxy this project used had the wrong sign.

THE PIPELINE
  1. rank the generated binders by how much of the chosen epitope they actually touch
  2. pick histidine positions on CPU (tfr1/pick_his_sites.py) - buried on binding, facing
     neutral or basic surface, not already salt-bridged, and spaced far enough apart to
     titrate independently
  3. fold each mutant complex with AlphaFold (the expensive step, hence steps 1 and 2)
  4. score the switch with tfr1/pka_switch.py: pKa in the complex versus pKa in the free
     binder, combined through the linkage equation

SIGN. A POSITIVE score means weaker binding at pH 6.0, which is what this challenge wants.
Negative means tighter in acid - the wrong way, and the easy mistake.

BUDGET. GPU hours are the constraint, not ideas. The loop stops on a wall clock, writes
whatever it has, and prints it. A partial result is useful; an overrun that loses the week's
quota is not.
"""
import os

os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ["XLA_PYTHON_CLIENT_ALLOCATOR"] = "platform"

import glob
import json
import subprocess
import sys
import time
import urllib.request

W = "/kaggle/working"
SHA = "main"                      # pin to a commit when running for real
BASE = ("https://raw.githubusercontent.com/wwwtalalahmed98-collab/"
        "ph-switch-binders/%s/tfr1/" % SHA)
TIME_BUDGET_MIN = 200             # sized against the GPU quota left, not against ambition
TOP_BINDERS = 10                  # how many generated binders to scan
MIN_ON_SITE = 4                   # epitope contacts a binder needs before it is worth it
t0 = time.time()


def log(m):
    print("[hisscan] %s [%.1f min]" % (m, (time.time() - t0) / 60), flush=True)


for f in ("pick_his_sites.py", "pka_switch.py", "propka_compat.py"):
    urllib.request.urlretrieve(BASE + f, "%s/%s" % (W, f))
sys.path.insert(0, W)
try:
    import propka  # noqa: F401
except ImportError:
    subprocess.run([sys.executable, "-m", "pip", "-q", "install", "propka"], check=True)

from Bio import PDB                                            # noqa: E402

import pick_his_sites                                         # noqa: E402
import pka_switch                                             # noqa: E402


def score_switch(pdb_path, label):
    """Binder and target chains are worked out from the file, then the switch is scored."""
    model = PDB.PDBParser(QUIET=True).get_structure("x", pdb_path)[0]
    binder, target = pick_his_sites._resolve(model, None, None)
    return pka_switch.analyse(pdb_path, binder, target, label=label)

# ---------------------------------------------------------------- find the previous run
src = (glob.glob("/kaggle/input/**/tfr1_designs.json", recursive=True)
       + glob.glob("%s/tfr1_designs.json" % W))
if not src:
    raise SystemExit("no tfr1_designs.json found. Attach the tfr1_generate notebook's "
                     "output as an input dataset and run again.")
DESIGNS = json.load(open(src[0]))
ROOT = os.path.dirname(src[0])
log("loaded %d designs from %s" % (len(DESIGNS), src[0]))

# Does predicted confidence track epitope accuracy on THIS task? On GLP-1R it anti-correlated
# (-0.96, -0.67), but that was a buried channel. This is a surface epitope and the finding may
# simply not transfer, so it is measured rather than assumed in either direction.
import numpy as np                                            # noqa: E402
ok = [d for d in DESIGNS if "on_site" in d and "true_iptm" in d]
if len(ok) > 4:
    r = float(np.corrcoef([d["true_iptm"] for d in ok], [d["on_site"] for d in ok])[0, 1])
    log("confidence vs on-site contacts across %d designs: r = %+.3f" % (len(ok), r))
    log("  (on GLP-1R this was negative; if it is positive here, the GLP-1R finding was"
        " about buried channels and does not generalise to surface epitopes)")

cands = sorted([d for d in DESIGNS if d.get("on_site", 0) >= MIN_ON_SITE],
               key=lambda d: (d["on_site"], d["true_iptm"]), reverse=True)[:TOP_BINDERS]
log("%d binders clear %d on-site contacts; scanning the top %d"
    % (len(cands), MIN_ON_SITE, len(cands)))
if not cands:
    raise SystemExit("no generated binder touches the epitope - nothing to put a switch on")

TARGETS = json.load(open(glob.glob("/kaggle/input/**/tfr1_targets.json", recursive=True)[0]
                         if glob.glob("/kaggle/input/**/tfr1_targets.json", recursive=True)
                         else "%s/tfr1_targets.json" % W))

from colabdesign import mk_afdesign_model, clear_mem          # noqa: E402

outdir = "%s/tfr1_hisscan" % W
os.makedirs(outdir, exist_ok=True)
records = []

for d in cands:
    if (time.time() - t0) / 60 > TIME_BUDGET_MIN:
        log("budget reached, stopping with %d records" % len(records))
        break
    wt_pdb = os.path.join(ROOT, "tfr1_%s" % d["site"], "%s.pdb" % d["name"])
    if not os.path.exists(wt_pdb):
        hits = glob.glob("/kaggle/input/**/%s.pdb" % d["name"], recursive=True)
        if not hits:
            log("%s: no structure found, skipped" % d["name"])
            continue
        wt_pdb = hits[0]

    try:
        plan = pick_his_sites.analyse(wt_pdb, quiet=True)
    except Exception as e:
        log("%s: site picking failed, %s %s" % (d["name"], type(e).__name__, str(e)[:120]))
        continue
    singles = plan["singles"]
    pairs = plan["pairs"][:3]
    if not singles:
        log("%s: no usable histidine position" % d["name"])
        continue

    seq = d["sequence"]
    first = min(r["pos"] for r in plan["interface"]) if plan["interface"] else 1
    muts = []
    for s in singles:
        i = s["pos"] - first
        if 0 <= i < len(seq) and seq[i] == s["wt"]:
            muts.append(([s["pos"]], seq[:i] + "H" + seq[i + 1:]))
    for p in pairs:
        ia, ib = p["a"] - first, p["b"] - first
        if 0 <= ia < len(seq) and 0 <= ib < len(seq):
            m = list(seq)
            m[ia] = m[ib] = "H"
            muts.append(([p["a"], p["b"]], "".join(m)))
    log("%s: %d mutants (%d single, %d pair)"
        % (d["name"], len(muts), len(singles), len(pairs)))

    M = TARGETS[d["site"]]
    clear_mem()
    model = mk_afdesign_model(protocol="binder", data_dir=W)
    model.prep_inputs(pdb_filename=M["target"], chain=M["chain"],
                      binder_len=d["length"], hotspot=M["hotspot"])

    for positions, mseq in [([], seq)] + muts:
        if (time.time() - t0) / 60 > TIME_BUDGET_MIN:
            break
        tag = "%s_%s" % (d["name"], "wt" if not positions
                         else "H" + "_".join(str(p) for p in positions))
        try:
            try:
                model.predict(seq=mseq, num_recycles=3, verbose=False)
            except TypeError:
                model.predict(seq=mseq, verbose=False)
            lg = model.aux["log"]
            out = "%s/%s.pdb" % (outdir, tag)
            model.save_pdb(out)
            rec = dict(name=tag, parent=d["name"], site=d["site"], positions=positions,
                       sequence=mseq, length=len(mseq),
                       true_iptm=round(float(lg.get("i_ptm", lg.get("ptm", 0.0))), 3),
                       plddt=round(float(lg.get("plddt", 0.0)), 3))
            try:
                sw = score_switch(out, tag)
                rec["log10_ratio_his"] = sw["log10_ratio_his"]
                rec["log10_ratio_all"] = sw["log10_ratio_all"]
            except Exception as e:
                rec["pka_error"] = "%s %s" % (type(e).__name__, str(e)[:140])
            records.append(rec)
            with open("%s/hisscan.jsonl" % outdir, "a") as fh:
                fh.write(json.dumps(rec) + "\n")
            log("  %-34s i_ptm %.3f  switch %+.2f"
                % (tag, rec["true_iptm"], rec.get("log10_ratio_his", 0.0)))
        except Exception as e:
            log("  %s FAILED %s %s" % (tag, type(e).__name__, str(e)[:160]))

json.dump(records, open("%s/tfr1_hisscan.json" % W, "w"), indent=0)
log("done: %d structures scored" % len(records))

if records:
    switches = [r for r in records if r.get("log10_ratio_his", 0) > 0 and r["positions"]]
    switches.sort(key=lambda r: r["log10_ratio_his"], reverse=True)
    print("\nbest predicted acid-release switches (positive = releases in acid):")
    print("%-34s %8s %8s %8s" % ("design", "switch", "i_ptm", "pLDDT"))
    for r in switches[:15]:
        print("%-34s %+8.2f %8.3f %8.3f"
              % (r["name"], r["log10_ratio_his"], r["true_iptm"], r["plddt"]))
    print("\n%d of %d placements switch the right way" % (len(switches), len(records)))
    print("\n===== COPY FROM HERE =====")
    for r in records:
        print(json.dumps(r))
    print("===== COPY TO HERE =====")
