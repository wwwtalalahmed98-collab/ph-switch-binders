"""Challenge 5, optional extra: predict the final designs against all three states.

THIS IS NOT REQUIRED FOR THE SUBMISSION. Track 3 needs name, sequence and
molecule_class, ranked. Structure models, ipTM and ipSAE are "encouraged", not
mandatory, and our ranking criterion - discrimination between the G12D peptide and
wild type - is already measured on CPU by pmhc/specificity_screen.py. This script adds
independent, structure-based evidence on top. If it cannot be run, the submission still
stands; it is just thinner.

WHAT IT DOES. For each of the 20 submitted sequences, AlphaFold2 predicts the complex
three times:

    g12d   the target as crystallised, peptide VVVGADGVGK      <- what we want bound
    wt     the same background with Asp6 truncated to Gly      <- must bind worse
    apo    the same background with no peptide at all          <- must bind worse

and records pLDDT, ipTM, pTM and the full PAE matrix for each. The ratios the
competition ranks on are KD(wt)/KD(g12d) and KD(apo)/KD(g12d); AF2 cannot give a KD, so
what we report is the confidence margin between states, clearly labelled as such.

WHY ALL THREE STATES SHARE ONE BACKGROUND. If the wild-type state came from 7OW5 and
the G12D state from 7OW6, any difference could be a difference between two crystals
rather than a difference our binder created. All three states are built from 7OW6, so
the only thing that varies is the peptide. 7OW5 is used solely to confirm the wild-type
sequence is what we think it is. This is the same decision recorded for the earlier
counter-screen, and it has a cost: the comparison is cleaner than reality, so the
measured margins are an upper bound.

THE APO STATE IS THE EASY ONE. An empty groove is the most different of the three, so
a large apo margin is weaker evidence than an equally large wild-type margin. Both are
reported; only the wild-type margin should carry much weight.

LESSONS BAKED IN FROM EARLIER RUNS:
  - /kaggle/working does not survive a committed run, so every number is PRINTED.
  - Child processes write to log FILES, never pipes: a 64 kB pipe buffer with an
    undrained parent stalled a worker on run 2 and cost about two GPU hours.
  - The script refuses to start if the AlphaFold parameters are missing, rather than
    discovering it half an hour in.
  - PAE and per-residue pLDDT are saved, because four challenges in a row reported
    "cannot compute ipSAE" for the single reason that we never saved the PAE matrix.

HOW TO USE. Paste the 20 sequences into SEQUENCES below (the submission detail table
has them), then run on a Kaggle GPU session.
"""
import json
import os
import subprocess
import sys

W = "/kaggle/working"
REPO = os.path.join(W, "ph-switch-binders")
OUR_REPO = "https://github.com/wwwtalalahmed98-collab/ph-switch-binders.git"

# ---------------------------------------------------------------------------
# PASTE THE SUBMISSION HERE: list of (name, sequence) from
# submission/challenge5_kras_g12d.csv
SEQUENCES = [
    # ("krasG12D_01", "..."),
]
# ---------------------------------------------------------------------------

MODELS = (0, 1)          # two AF2 models per state; enough to see disagreement
RECYCLES = 3
BUDGET_MIN = 330         # stop cleanly before the Kaggle session limit


def sh(cmd, cwd=None, check=True):
    print("$ %s" % cmd, flush=True)
    r = subprocess.run(cmd, shell=True, cwd=cwd, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, text=True)
    if (r.stdout or "").strip():
        print(r.stdout[-3000:], flush=True)
    if check and r.returncode != 0:
        raise SystemExit("command failed (%d): %s" % (r.returncode, cmd))


def setup():
    print("=" * 72)
    print("SETUP")
    print("=" * 72)
    if not SEQUENCES:
        raise SystemExit("SEQUENCES is empty - paste the submission set in first")
    print("sequences to predict: %d" % len(SEQUENCES))
    bad = [n for n, s in SEQUENCES if not s or set(s) - set("ACDEFGHIKLMNPQRSTVWY")]
    if bad:
        raise SystemExit("not plain protein sequences: %s" % bad[:5])

    try:
        import colabdesign  # noqa: F401
    except ImportError:
        sh("pip -q install git+https://github.com/sokrypton/ColabDesign.git@v1.1.1 "
           "biopython")

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
        print("  %s.cif %d bytes" % (pdb_id, os.path.getsize(p)))

    if not os.path.exists("%s/params/params_model_1_ptm.npz" % W):
        os.makedirs("%s/params" % W, exist_ok=True)
        sh("wget -q -O %s/af2.tar https://storage.googleapis.com/alphafold/"
           "alphafold_params_2022-12-06.tar" % W)
        sh("tar -xf %s/af2.tar -C %s/params" % (W, W))
        os.remove("%s/af2.tar" % W)
    import glob
    n = len(glob.glob("%s/params/params_model_*.npz" % W))
    print("AlphaFold parameter files: %d" % n)
    if n == 0:
        raise SystemExit("no AlphaFold parameters - refusing to start")

    # Build all three states from 7OW6, verified by the build script itself.
    sh("%s -W ignore pmhc/build_states.py" % sys.executable, cwd=REPO)
    sh("%s -W ignore pmhc/build_counter_states.py" % sys.executable, cwd=REPO)
    states = os.path.join(REPO, "data", "pmhc", "states")
    need = ["target_state_g12d.pdb", "target_state_wt.pdb", "target_state_apo.pdb"]
    for f in need:
        p = os.path.join(states, f)
        if not os.path.exists(p):
            raise SystemExit("missing state file %s" % p)
        print("  %s %d bytes" % (f, os.path.getsize(p)))
    return states


def main():
    states = setup()
    print()
    print("=" * 72)
    print("THREE-STATE PREDICTION")
    print("=" * 72)
    print("This script is a harness: the per-design prediction loop is the same")
    print("three-state screen already unit-tested in pmhc/test_screen.py, which")
    print("covers the merged-chain layout, the apo state reporting no phantom")
    print("peptide contacts, and a deliberately mismatched state being rejected.")
    print()
    print("states directory: %s" % states)
    print("models: %s, recycles: %d, budget: %d min"
          % (list(MODELS), RECYCLES, BUDGET_MIN))
    print()
    print("NOT YET WIRED: the prediction loop itself still has to be connected to")
    print("kaggle/cellPMHCc_screen.py's scoring functions. This file currently does")
    print("the setup and verification only, and stops here rather than pretending")
    print("to have produced numbers.")
    results_path = os.path.join(W, "pmhc_threestate.json")
    json.dump({"status": "setup_only", "n_sequences": len(SEQUENCES)},
              open(results_path, "w"), indent=1)
    print("wrote %s" % results_path)


main()
