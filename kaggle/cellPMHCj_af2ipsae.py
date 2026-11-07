"""Three-state AlphaFold2 screen of the designs, with ipSAE computed at last.

WHY AF2 AND NOT BOLTZ-2. The plan was Boltz-2, because Adaptyv select designs for the
wet lab on Boltz-2 ipSAE. It does not install on Kaggle: the image is Python 3.13, and
boltz 2.0.3 pins dm-tree==0.1.8, which has no cp313 wheel and fails to compile. The
same pins force source builds of numpy 1.26 and scipy 1.13, which is what made the
first attempt sit silently for twenty-two minutes before failing.

Rather than spend the remaining time on dependency surgery, this uses AlphaFold2
through ColabDesign, which this project already runs reliably on Kaggle, and computes
ipSAE from its output. That is a real downgrade and it is stated as one: the
meta-analysis of 3,766 experimentally characterised binders ranked AF3 ipSAE_min best,
with Boltz-1 close behind, and AF2 below them. But AF2 ipSAE is the mode ipsae.py was
originally written for, and it is enormously better than what we have now, which is no
structure prediction at all.

WHAT IT DOES. For each candidate, three states, two AF2 models each:

    g12d   the target as crystallised, peptide VVVGADGVGK   <- should bind
    wt     the same background, Asp6 truncated to Gly       <- should bind worse
    apo    the same background, no peptide                  <- should bind worse

Two models per state is not redundancy, it is the noise floor. A margin between states
only counts if it is bigger than the disagreement between two models of the SAME state.
Without that, we would read model-to-model scatter as specificity.

All three states are built from 7OW6 so only the peptide differs between them; 7OW5 is
used solely to confirm the wild-type sequence. The cost of that choice is that wild type
and G12D share one backbone, which makes the comparison cleaner than reality, so the
margins are an upper bound.

WHAT IT SAVES, AND WHY IT MATTERS. The PAE matrix and the PER-RESIDUE pLDDT, for every
prediction. Four challenges in a row reported "cannot compute ipSAE" for the single
reason that the PAE matrix was thrown away. ipsae.py substitutes zeros for a missing
pLDDT and runs anyway, producing a number that looks fine and means less than it
appears to - so the per-residue array is saved and its presence is asserted, not hoped
for.
"""
import glob
import json
import os
import subprocess
import sys
import time

import numpy as np

W = "/kaggle/working"
REPO = os.path.join(W, "ph-switch-binders")
OUR_REPO = "https://github.com/wwwtalalahmed98-collab/ph-switch-binders.git"
IPSAE_URL = "https://raw.githubusercontent.com/DunbrackLab/IPSAE/main/ipsae.py"
CAND_GLOB = "/kaggle/input/*/*candidates*.csv"

MODELS = (0, 1)
MAX_CANDIDATES = 24
BUDGET_MIN = 400
PAE_CUTOFF, DIST_CUTOFF = 10, 10

PEP_G12D = "VVVGADGVGK"
PEP_WT = "VVVGAGGVGK"
T3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
      "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
      "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
      "TYR": "Y", "VAL": "V"}


def sh(cmd, check=True, stream=False):
    print("$ %s" % cmd, flush=True)
    if stream:
        r = subprocess.run(cmd, shell=True)
        out = ""
    else:
        r = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, text=True)
        out = r.stdout or ""
        if out.strip():
            print(out[-2000:], flush=True)
    sys.stdout.flush()
    if check and r.returncode != 0:
        raise SystemExit("command failed (%d): %s" % (r.returncode, cmd))
    return out


def setup():
    print("=" * 72)
    print("SETUP")
    print("=" * 72)
    try:
        import colabdesign  # noqa: F401
        print("colabdesign already present")
    except ImportError:
        sh("pip install git+https://github.com/sokrypton/ColabDesign.git@v1.1.1 "
           "biopython", stream=True)

    sh("wget -q -O %s/ipsae.py %s" % (W, IPSAE_URL))
    if not os.path.exists("%s/ipsae.py" % W):
        raise SystemExit("could not download ipsae.py")
    print("ipsae.py: %d bytes" % os.path.getsize("%s/ipsae.py" % W))

    if not os.path.isdir(REPO):
        sh("git clone --depth 1 %s %s" % (OUR_REPO, REPO))
    d = os.path.join(REPO, "data", "pmhc")
    os.makedirs(d, exist_ok=True)
    import urllib.request
    for pdb_id in ("7OW6", "7OW5"):
        p = os.path.join(d, "%s.cif" % pdb_id)
        if not os.path.exists(p):
            urllib.request.urlretrieve(
                "https://files.rcsb.org/download/%s.cif" % pdb_id, p)

    if not os.path.exists("%s/params/params_model_1_ptm.npz" % W):
        os.makedirs("%s/params" % W, exist_ok=True)
        sh("wget -q -O %s/af2.tar https://storage.googleapis.com/alphafold/"
           "alphafold_params_2022-12-06.tar" % W)
        sh("tar -xf %s/af2.tar -C %s/params" % (W, W))
        os.remove("%s/af2.tar" % W)
    n = len(glob.glob("%s/params/params_model_*.npz" % W))
    print("AlphaFold parameter files: %d" % n)
    if n == 0:
        raise SystemExit("no AlphaFold parameters - refusing to start")

    # cd into the repo: build_states.py resolves its data paths relative to its own
    # location, but the script is invoked by relative path, so the working directory
    # has to be the repo root or nothing is found.
    sh("cd %s && %s -W ignore pmhc/build_states.py" % (REPO, sys.executable))
    return d


def build_states():
    """Three states from 7OW6, each verified before use."""
    from Bio import PDB
    sh("cd %s && %s -W ignore pmhc/build_counter_states.py"
       % (REPO, sys.executable))
    src = os.path.join(REPO, "data", "pmhc", "states")
    states = {}
    for tag, chains in (("g12d", "A,B,C"), ("wt", "A,B,C"), ("apo", "A,B")):
        path = os.path.join(src, "target_state_%s.pdb" % tag)
        if not os.path.exists(path):
            raise SystemExit("missing %s" % path)
        m = PDB.PDBParser(QUIET=True).get_structure("s", path)[0]
        ch = {c.id: [r for r in c if r.id[0] == " " and r.get_resname() in T3]
              for c in m}
        pseq = ("".join(T3[r.get_resname()] for r in ch["C"])
                if "C" in ch and ch["C"] else None)
        states[tag] = dict(path=path, chains=chains, peptide=pseq,
                           n=sum(len(v) for v in ch.values()))
    if states["g12d"]["peptide"] != PEP_G12D:
        raise SystemExit("g12d peptide is %r" % states["g12d"]["peptide"])
    if states["wt"]["peptide"] != PEP_WT:
        raise SystemExit("wt peptide is %r - the mutation was not applied"
                         % states["wt"]["peptide"])
    if states["apo"]["peptide"]:
        raise SystemExit("apo state still carries a peptide")
    for tag, s in states.items():
        print("state %-5s %3d residues  peptide %s" % (tag, s["n"], s["peptide"]))
    return states


def load_candidates():
    import csv
    for p in glob.glob(CAND_GLOB):
        rows = list(csv.DictReader(open(p)))
        print("candidates from %s: %d" % (p, len(rows)))
        AA = set("ACDEFGHIKLMNPQRSTVWY")
        bad = [r["name"] for r in rows
               if set(r["seq"]) - AA or not (10 <= len(r["seq"]) <= 250)]
        if bad:
            raise SystemExit("invalid sequences: %s" % bad[:6])
        if len(rows) > MAX_CANDIDATES:
            print("taking the leading %d (file is ordered best-first)"
                  % MAX_CANDIDATES)
            rows = rows[:MAX_CANDIDATES]
        return rows
    raise SystemExit("no candidate csv found at %s" % CAND_GLOB)


def run_ipsae(pdb, pae_npy, plddt_npy):
    """Convert AF2 output to what ipsae.py reads, then run it."""
    sys.path.insert(0, os.path.join(REPO, "pmhc"))
    import make_ipsae_inputs as mk
    try:
        js, n, has_plddt = mk.convert(pdb, pae_npy, plddt_npy)
    except SystemExit as e:
        return None, "convert refused: %s" % e
    if not has_plddt:
        return None, "pLDDT missing - ipSAE would run on zeros"
    out = sh("cd %s && python ipsae.py %s %s %d %d"
             % (W, js, pdb, PAE_CUTOFF, DIST_CUTOFF), check=False)
    rows = []
    for f in glob.glob(os.path.splitext(js)[0] + "*.txt") + \
            glob.glob(os.path.join(os.path.dirname(pdb), "*ipsae*.txt")):
        for line in open(f):
            parts = line.split()
            if len(parts) > 6 and parts[0] not in ("Chn1", "#"):
                rows.append(parts)
    if not rows:
        return None, "ipsae.py produced no rows: %s" % out[-300:]
    return rows, None


def main():
    setup()
    states = build_states()
    cands = load_candidates()

    sys.path.insert(0, os.path.join(REPO, "pmhc"))
    from colabdesign import clear_mem, mk_afdesign_model

    print()
    print("=" * 72)
    print("PREDICTION - %d candidates x 3 states x %d models = %d runs"
          % (len(cands), len(MODELS), len(cands) * 3 * len(MODELS)))
    print("=" * 72)

    results = []
    t0 = time.time()
    ipsae_fail = 0
    for i, c in enumerate(cands, 1):
        if (time.time() - t0) / 60 > BUDGET_MIN:
            print("budget reached after %d candidates" % (i - 1))
            break
        rec = {"name": c["name"], "seq": c["seq"], "length": len(c["seq"]),
               "states": {}}
        try:
            for tag, st in states.items():
                per_model = []
                for mi in MODELS:
                    clear_mem()
                    model = mk_afdesign_model(protocol="binder", data_dir=W)
                    model.prep_inputs(pdb_filename=st["path"], chain=st["chains"],
                                      binder_len=len(c["seq"]))
                    model.predict(seq=c["seq"], models=[mi], verbose=False)
                    lg = model.aux.get("log", {})
                    stem = "%s/%s_%s_m%d" % (W, c["name"], tag, mi)
                    model.save_pdb(stem + ".pdb")
                    pae = model.aux.get("pae")
                    pr = model.aux.get("plddt")
                    if pae is None or pr is None:
                        raise RuntimeError("AF2 returned no pae/plddt")
                    np.save(stem + "_pae.npy", np.asarray(pae))
                    np.save(stem + "_plddt.npy", np.asarray(pr).reshape(-1))
                    entry = dict(model=mi,
                                 i_ptm=round(float(lg.get("i_ptm", 0)), 4),
                                 ptm=round(float(lg.get("ptm", 0)), 4),
                                 plddt=round(float(lg.get("plddt", 0)), 4))
                    rows, err = run_ipsae(stem + ".pdb", stem + "_pae.npy",
                                          stem + "_plddt.npy")
                    if rows is None:
                        ipsae_fail += 1
                        entry["ipsae_error"] = err
                    else:
                        entry["ipsae_rows"] = rows
                    per_model.append(entry)
                iptms = [x["i_ptm"] for x in per_model]
                rec["states"][tag] = dict(
                    per_model=per_model,
                    i_ptm=round(float(np.mean(iptms)), 4),
                    i_ptm_spread=round(float(np.ptp(iptms)), 4))
        except Exception as e:                       # noqa: BLE001
            print("%s FAILED %s %s" % (c["name"], type(e).__name__, str(e)[:200]))
            continue

        g, w, a = rec["states"]["g12d"], rec["states"]["wt"], rec["states"]["apo"]
        noise = max(g["i_ptm_spread"], w["i_ptm_spread"], a["i_ptm_spread"])
        rec["margin_wt"] = round(g["i_ptm"] - w["i_ptm"], 4)
        rec["margin_apo"] = round(g["i_ptm"] - a["i_ptm"], 4)
        rec["noise"] = noise
        # a margin only counts if it clears the disagreement between two models of
        # the SAME state; otherwise it is scatter being read as specificity
        rec["measurable"] = bool(min(rec["margin_wt"], rec["margin_apo"]) > noise)
        results.append(rec)
        print("  %2d/%d %-24s i_ptm %.3f  margin_wt %+.3f  apo %+.3f  noise %.3f %s"
              % (i, len(cands), c["name"][:24], g["i_ptm"], rec["margin_wt"],
                 rec["margin_apo"], noise,
                 "MEASURABLE" if rec["measurable"] else ""), flush=True)
        json.dump(results, open("%s/af2_ipsae_results.json" % W, "w"), indent=1)

    json.dump(results, open("%s/af2_ipsae_results.json" % W, "w"), indent=1)
    print()
    print("=" * 72)
    print("RESULTS (%d designs, ipSAE failures %d)" % (len(results), ipsae_fail))
    print("=" * 72)
    for r in results:
        print(json.dumps({k: v for k, v in r.items() if k != "seq"}))
    print()
    print("NOT DONE HERE: the chain-pair selection and the final ranking are done")
    print("locally by pmhc/parse_boltz.py, so the arithmetic that decides the")
    print("submission is reviewable rather than buried in a notebook log.")
    print()
    print("AF2, not Boltz-2. Adaptyv select on Boltz-2 ipSAE; boltz will not install")
    print("on this image (Python 3.13 vs a dm-tree==0.1.8 pin). These numbers are")
    print("AF2 ipSAE and should not be compared against Boltz-2 thresholds.")


main()
