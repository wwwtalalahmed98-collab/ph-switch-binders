"""Score the designs with Boltz-2 and ipSAE - the metric Adaptyv actually selects on.

WHY THIS, AND WHY NOW.

Two findings changed what this submission should optimise:

1. Adaptyv, who run the wet lab, say in their own Nipah competition write-up that they
   chose which designs to test using **Boltz-2 ipSAE**: "we opted to use the Boltz-2
   ipSAE score for filtering and selected the 600 most promising designs based on this
   metric." Earlier rounds used AF2 iPAE/ipTM; they moved to Boltz-2 ipSAE.

2. A meta-analysis of 3,766 experimentally characterised binders (bioRxiv
   2025.08.14.670059) found **ipSAE_min is the single best predictor of in vitro
   binding**, beating ipTM and iPAE across 15 targets. Its maximum-F1 threshold for
   AF3 was 0.61.

This project has written "cannot compute ipSAE" in four consecutive challenges, for the
single reason that we never saved a PAE matrix. Everything we have ranked on so far has
been a proxy for a proxy. This computes the real thing.

WHAT IT PREDICTS. For every candidate, two states, each as three chains:

    A  MHC-I platform, 179 aa (HLA-A*11:01, from 7OW6)
    B  peptide - VVVGADGVGK (G12D) or VVVGAGGVGK (wild type)
    C  the designed binder

From each prediction we take ipSAE for the binder-peptide pair and the binder-MHC pair
separately. That separation is the whole point of the design strategy: we want a high
binder-peptide score and a comparatively low binder-MHC score, because a binder that
grips the MHC binds every peptide the MHC presents.

The specificity readout is then ipSAE(binder,peptide) on G12D minus the same on wild
type - a confidence margin between two states, which is as close to KD(WT)/KD(G12D) as
a structure predictor can get. It is NOT a KD ratio and is not reported as one.

SINGLE-SEQUENCE MODE, AND THE CHECK THAT JUSTIFIES IT. Every chain is given an empty
MSA. For the designed binder that is mandatory - a de novo sequence has no homologs, and
feeding it one would be meaningless. For the MHC it is a choice: it makes all
predictions identical in setup, fast, and independent of an external MSA server that
could rate-limit us halfway through. MHC-I is among the most structurally represented
folds in the PDB, so Boltz should manage it unaided - but "should" is not evidence, so
the script FIRST predicts the target alone and superposes it on the 7OW6 crystal
structure. If that RMSD is bad, the groove is wrong and every number afterwards is
meaningless, so the script says so loudly and stops.

Run on a Kaggle GPU session with internet enabled.
"""
import glob
import json
import os
import shutil
import subprocess
import sys
import time

W = "/kaggle/working"
REPO = os.path.join(W, "ph-switch-binders")
OUR_REPO = "https://github.com/wwwtalalahmed98-collab/ph-switch-binders.git"
IPSAE_URL = "https://raw.githubusercontent.com/DunbrackLab/IPSAE/main/ipsae.py"

# Candidates are uploaded as a Kaggle dataset because the sequences are withheld from
# GitHub until the submission window closes.
# Candidate sequences reach this notebook in one of two ways, never through GitHub:
# a JSON file staged by a preceding cell (the usual route - the sequences are withheld
# while the submission window is open), or a private Kaggle dataset.
CAND_JSON = "/kaggle/working/cands.json"
CAND_GLOB = "/kaggle/input/*/*candidates*.csv"

PEP_G12D = "VVVGADGVGK"
PEP_WT = "VVVGAGGVGK"
PAE_CUTOFF = 10
DIST_CUTOFF = 10
RMSD_LIMIT = 2.5          # target-alone prediction vs 7OW6, in angstroms
BUDGET_MIN = 400          # stop cleanly before the session limit
# 56 candidates x 2 states is ~112 predictions, which risks running past a Kaggle
# session. The CSV is ordered best-first by the shuffle-controlled score, so taking
# the leading slice keeps the candidates most likely to be submitted and still leaves
# room to promote ones currently ranked below 20.
MAX_CANDIDATES = 32


def sh(cmd, check=True, quiet=False, stream=False):
    """Run a shell command.

    stream=True lets the child write straight to our stdout instead of being
    captured. That matters for anything slow: the first run of this script sat on
    `pip install boltz` for over thirty minutes with no output at all, because the
    output was captured and only printed on completion, so there was no way to tell
    a working download from a wedged dependency resolver. Long commands stream.
    """
    if not quiet:
        print("$ %s" % cmd, flush=True)
    if stream:
        r = subprocess.run(cmd, shell=True)
        out = ""
    else:
        r = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, text=True)
        out = r.stdout or ""
        if out.strip() and not quiet:
            print(out[-2500:], flush=True)
    sys.stdout.flush()
    if check and r.returncode != 0:
        raise SystemExit("command failed (%d): %s" % (r.returncode, cmd))
    return out


def setup():
    print("=" * 72)
    print("SETUP")
    print("=" * 72)
    # streamed: this is the step that previously ran blind for 30+ min
    sh("pip install boltz biopython", stream=True)
    sh("wget -q -O %s/ipsae.py %s" % (W, IPSAE_URL))
    if not os.path.exists(os.path.join(W, "ipsae.py")):
        raise SystemExit("could not download ipsae.py")
    print("ipsae.py: %d bytes" % os.path.getsize(os.path.join(W, "ipsae.py")))

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
    sh("%s -W ignore pmhc/build_states.py" % sys.executable, check=True)
    return d


def target_sequence():
    """The 179-residue MHC platform we designed against, read from the built file."""
    from Bio import PDB
    T3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
          "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
          "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
          "TYR": "Y", "VAL": "V"}
    p = os.path.join(REPO, "data", "pmhc", "states", "target_design.pdb")
    m = PDB.PDBParser(QUIET=True).get_structure("t", p)[0]
    heavy = "".join(T3[r.get_resname()] for r in m["A"]
                    if r.id[0] == " " and r.get_resname() in T3)
    pep = "".join(T3[r.get_resname()] for r in m["C"]
                  if r.id[0] == " " and r.get_resname() in T3)
    if pep != PEP_G12D:
        raise SystemExit("target peptide is %r, expected %r" % (pep, PEP_G12D))
    if len(heavy) != 179:
        raise SystemExit("MHC platform is %d aa, expected 179" % len(heavy))
    print("MHC platform %d aa, peptide %s" % (len(heavy), pep))
    return heavy


def write_yaml(path, chains):
    """chains: list of (id, sequence). Empty MSA on every chain - see the module
    docstring for why, and for the check that justifies it."""
    lines = ["version: 1", "sequences:"]
    for cid, seq in chains:
        lines += ["  - protein:",
                  "      id: %s" % cid,
                  "      sequence: %s" % seq,
                  "      msa: empty"]
    open(path, "w").write("\n".join(lines) + "\n")
    return path


def run_boltz(yaml_path, out_dir, verbose=False):
    # mmCIF, not PDB: ipsae.py's Boltz mode is documented as
    #   python ipsae.py <boltz_pae_npz> <boltz_cif> <pae_cutoff> <dist_cutoff>
    # and it infers the format from the extension. Handing it a .pdb alongside a
    # Boltz .npz is not a combination it advertises, and a scorer that silently
    # misreads chains would give us confident nonsense - the exact failure this
    # project has had three times.
    sh("boltz predict %s --out_dir %s --output_format mmcif --diffusion_samples 1 "
       "--override" % (yaml_path, out_dir), check=False, quiet=not verbose,
       stream=verbose)
    return out_dir


def find_outputs(out_dir, stem):
    pdb = glob.glob(os.path.join(out_dir, "boltz_results_%s" % stem,
                                 "predictions", stem, "*.pdb"))
    cif = glob.glob(os.path.join(out_dir, "boltz_results_%s" % stem,
                                 "predictions", stem, "*.cif"))
    pae = glob.glob(os.path.join(out_dir, "boltz_results_%s" % stem,
                                 "predictions", stem, "pae_*.npz"))
    conf = glob.glob(os.path.join(out_dir, "boltz_results_%s" % stem,
                                  "predictions", stem, "confidence_*.json"))
    # cif first: that is what ipsae.py's Boltz mode expects.
    return (cif or pdb), pae, conf


def sanity_check_target(heavy, work):
    """Predict the target alone and compare it to the crystal structure.

    If single-sequence Boltz cannot reproduce the MHC groove, every binder score
    computed afterwards is scoring the wrong pocket. This is cheap and decisive, so it
    runs before anything else.
    """
    print()
    print("=" * 72)
    print("SANITY CHECK - can single-sequence Boltz reproduce the MHC groove?")
    print("=" * 72)
    stem = "targetonly"
    y = write_yaml(os.path.join(work, "%s.yaml" % stem),
                   [("A", heavy), ("B", PEP_G12D)])
    t0 = time.time()
    run_boltz(y, work, verbose=True)
    print("prediction took %.0f s" % (time.time() - t0))
    struct, pae, conf = find_outputs(work, stem)
    if not struct:
        raise SystemExit("Boltz produced no structure for the target-only check - "
                         "fix that before scoring anything")
    print("structure: %s" % os.path.basename(struct[0]))
    if conf:
        c = json.load(open(conf[0]))
        print("confidence: %s" % {k: round(v, 3) for k, v in c.items()
                                  if isinstance(v, (int, float))})

    from Bio import PDB
    import numpy as np
    ref = PDB.PDBParser(QUIET=True).get_structure(
        "r", os.path.join(REPO, "data", "pmhc", "states", "target_design.pdb"))[0]
    parser = (PDB.MMCIFParser(QUIET=True) if struct[0].endswith(".cif")
              else PDB.PDBParser(QUIET=True))
    pred = parser.get_structure("p", struct[0])[0]
    ref_ca = [r["CA"] for r in ref["A"] if r.id[0] == " " and "CA" in r]
    pred_chain = list(pred)[0]
    pred_ca = [r["CA"] for r in pred_chain if r.id[0] == " " and "CA" in r]
    n = min(len(ref_ca), len(pred_ca))
    sup = PDB.Superimposer()
    sup.set_atoms(ref_ca[:n], pred_ca[:n])
    print()
    print("MHC platform vs 7OW6 crystal: RMSD %.2f A over %d CA atoms"
          % (sup.rms, n))
    if sup.rms > RMSD_LIMIT:
        raise SystemExit(
            "RMSD %.2f A exceeds the %.1f A limit. Single-sequence Boltz is NOT "
            "reproducing the groove, so binder scores against it would be "
            "meaningless. Re-run with --use_msa_server before trusting anything."
            % (sup.rms, RMSD_LIMIT))
    print("PASS - the groove is reproduced, so binder scores against it mean "
          "something.")
    return True


def score_ipsae(struct, pae):
    """Run the official ipsae.py and return its per-chain-pair rows."""
    out = sh("cd %s && python ipsae.py %s %s %d %d"
             % (W, pae, struct, PAE_CUTOFF, DIST_CUTOFF), check=False, quiet=True)
    base = os.path.splitext(pae)[0]
    rows = []
    for f in glob.glob(base + "*.txt") + glob.glob(
            os.path.join(os.path.dirname(struct), "*ipsae*.txt")):
        for line in open(f):
            parts = line.split()
            if len(parts) > 6 and parts[0] not in ("Chn1", "#"):
                rows.append(parts)
    if not rows:
        print("    ipsae.py produced no rows; its output was:")
        print("    %s" % out[-600:].replace("\n", "\n    "))
    return rows


def main():
    setup()
    heavy = target_sequence()
    work = os.path.join(W, "boltz")
    os.makedirs(work, exist_ok=True)

    sanity_check_target(heavy, work)

    cands = []
    if os.path.exists(CAND_JSON):
        raw = json.load(open(CAND_JSON))
        cands = [{"name": n, "seq": s} for n, s in raw]
        print()
        print("candidates staged by a preceding cell: %d" % len(cands))
    else:
        import csv
        for p in glob.glob(CAND_GLOB):
            cands = list(csv.DictReader(open(p)))
            print()
            print("candidates loaded from %s: %d" % (p, len(cands)))
            break
    if not cands:
        print()
        print("No candidates found (neither %s nor %s)." % (CAND_JSON, CAND_GLOB))
        print("The sanity check above has already answered the important question:")
        print("whether single-sequence Boltz reproduces this target's groove.")
        return

    # Refuse to spend GPU hours on anything that is not a plain protein sequence.
    AA = set("ACDEFGHIKLMNPQRSTVWY")
    bad = [c["name"] for c in cands
           if not c.get("seq") or set(c["seq"]) - AA or not (10 <= len(c["seq"]) <= 250)]
    if bad:
        raise SystemExit("not valid designs: %s" % bad[:8])
    print("all %d sequences validated: standard residues, 10-250 aa" % len(cands))
    if len(cands) > MAX_CANDIDATES:
        print("taking the leading %d of %d (the file is ordered best-first)"
              % (MAX_CANDIDATES, len(cands)))
        cands = cands[:MAX_CANDIDATES]

    print()
    print("=" * 72)
    print("SCORING - two states per candidate (G12D and wild type)")
    print("=" * 72)
    results = []
    t0 = time.time()
    for i, c in enumerate(cands, 1):
        if (time.time() - t0) / 60 > BUDGET_MIN:
            print("budget reached after %d candidates" % (i - 1))
            break
        row = {"name": c["name"], "seq": c["seq"]}
        for state, pep in (("g12d", PEP_G12D), ("wt", PEP_WT)):
            stem = "%s_%s" % (c["name"], state)
            y = write_yaml(os.path.join(work, "%s.yaml" % stem),
                           [("A", heavy), ("B", pep), ("C", c["seq"])])
            run_boltz(y, work)
            struct, pae, conf = find_outputs(work, stem)
            if not struct or not pae:
                row["%s_ok" % state] = False
                continue
            row["%s_ok" % state] = True
            if conf:
                cj = json.load(open(conf[0]))
                for k in ("confidence_score", "ptm", "iptm", "complex_plddt"):
                    if k in cj:
                        row["%s_%s" % (state, k)] = round(float(cj[k]), 4)
            row["%s_ipsae_rows" % state] = score_ipsae(struct[0], pae[0])
        results.append(row)
        if i % 5 == 0 or i == len(cands):
            print("  %d/%d candidates  (%.0f min)"
                  % (i, len(cands), (time.time() - t0) / 60), flush=True)
            json.dump(results, open(os.path.join(W, "boltz_results.json"), "w"),
                      indent=1)

    json.dump(results, open(os.path.join(W, "boltz_results.json"), "w"), indent=1)
    print()
    print("=" * 72)
    print("RESULTS - printed in full, because /kaggle/working does not persist")
    print("=" * 72)
    for r in results:
        print(json.dumps({k: v for k, v in r.items() if k != "seq"}))
    print()
    print("NOT DONE: ipSAE rows above are raw ipsae.py output. The chain-pair")
    print("selection (binder-peptide vs binder-MHC) and the G12D minus wild-type")
    print("margin are computed locally by pmhc/parse_boltz.py, so that the")
    print("arithmetic is reviewable rather than buried in a notebook log.")


main()
