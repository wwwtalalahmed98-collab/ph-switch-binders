"""Challenge 5, step 3: redesign the docked scaffold binders with ProteinMPNN.

WHAT THIS IS. Step 3 of the Liu et al. pMHC-I binder pipeline (Science 389:386-391,
2025), run on their published scaffolds after our dock_scaffolds.py has placed them on
the KRAS-G12D / HLA-A*11:01 target. Every amino acid in the binder is replaced; the
backbone is kept.

WHY IT RUNS ENTIRELY FROM PUBLIC SOURCES. We still have no Kaggle API token, so nothing
can be uploaded to this notebook. Instead it clones our own public repository and the
Baker lab's, downloads the two crystal structures from RCSB, and rebuilds the target and
the docked poses from scratch. That is slower than uploading 63 files but it also means
this notebook is a complete, reproducible record: anyone can run it and get our inputs.

WHY VANILLA ProteinMPNN AND NOT ColabDesign'S. Two reasons. The published pipeline uses
the vanilla release with the v_48_020 weights, which are committed in the Baker repo, so
this matches the method we are citing. And on the previous run ColabDesign's wrapper
reported `cysteines removed: True` while emitting sequences containing cysteine - the
flag was silently ignored and I reported it as working without checking the output. The
vanilla script's --omit_AAs is checked here against the actual sequences, not trusted.

LESSONS FROM EARLIER RUNS THAT ARE BAKED IN:
  - /kaggle/working does not survive a committed run, so every result is PRINTED as
    well as saved. The printed log is the artifact.
  - Child processes get their own log FILES, never pipes. A 64 kB pipe buffer with an
    undrained parent cost ~2 GPU hours on run 2, a risk I had explicitly dismissed.
  - The script refuses to start if any input is missing, rather than discovering it
    30 minutes in.
  - The screen is imported from pmhc/screen.py, not copy-pasted. Four divergent copies
    are how the proline blind spot survived two days.

THIS IS A CPU JOB. It uses no GPU quota. Of the 30 h weekly Kaggle GPU allowance,
~14 h remain and all of it is reserved for the three-state prediction pass.
"""
import glob
import json
import os
import shutil
import subprocess
import sys
import urllib.request

W = "/kaggle/working"
REPO = os.path.join(W, "ph-switch-binders")
BAKER = os.path.join(W, "pMHCI_binder_design")

OUR_REPO = "https://github.com/wwwtalalahmed98-collab/ph-switch-binders.git"
BAKER_REPO = "https://github.com/77forest/pMHCI_binder_design.git"

# ProteinMPNN settings. Temperatures and counts follow the ranges used in the paper's
# pipeline; KEEP_PER_POSE limits how many sequences from one backbone can reach the GPU
# step, so that 20 final designs are not 20 variants of one pose.
N_SEQ = 16              # sequences per pose per temperature
TEMPERATURES = (0.1, 0.2, 0.3)
KEEP_PER_POSE = 6
OMIT_AA = "CX"          # cysteine: unpaired cysteines misbehave when expressed

PEP_G12D = "VVVGADGVGK"
MUT_POS = 6
PEP_LEN = len(PEP_G12D)


def sh(cmd, cwd=None, check=True):
    print("$ %s" % cmd, flush=True)
    r = subprocess.run(cmd, shell=True, cwd=cwd, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, text=True)
    out = r.stdout or ""
    if out.strip():
        print(out[-4000:], flush=True)
    if check and r.returncode != 0:
        raise SystemExit("command failed (%d): %s" % (r.returncode, cmd))
    return out


def setup():
    print("=" * 72)
    print("SETUP - everything from public sources, nothing uploaded")
    print("=" * 72)
    sh("pip -q install biopython")

    for url, dest in ((OUR_REPO, REPO), (BAKER_REPO, BAKER)):
        if not os.path.isdir(dest):
            sh("git clone --depth 1 %s %s" % (url, dest))

    # The two crystal structures. build_states.py reads them from data/pmhc/.
    d = os.path.join(REPO, "data", "pmhc")
    os.makedirs(d, exist_ok=True)
    for pdb_id in ("7OW6", "7OW5"):
        p = os.path.join(d, "%s.cif" % pdb_id)
        if not os.path.exists(p):
            url = "https://files.rcsb.org/download/%s.cif" % pdb_id
            print("downloading %s" % url, flush=True)
            urllib.request.urlretrieve(url, p)
        print("  %s.cif  %d bytes" % (pdb_id, os.path.getsize(p)), flush=True)

    # Guard: refuse to continue if the pieces we depend on are absent.
    need = [
        os.path.join(REPO, "pmhc", "build_states.py"),
        os.path.join(REPO, "pmhc", "dock_scaffolds.py"),
        os.path.join(REPO, "pmhc", "shortlist_poses.py"),
        os.path.join(REPO, "pmhc", "screen.py"),
        os.path.join(BAKER, "software", "mpnn", "protein_mpnn_run.py"),
        os.path.join(BAKER, "software", "mpnn", "vanilla_model_weights", "v_48_020.pt"),
    ]
    missing = [p for p in need if not os.path.exists(p)]
    if missing:
        for p in missing:
            print("MISSING: %s" % p)
        raise SystemExit("refusing to start with missing inputs")
    n_scaf = len(glob.glob(os.path.join(BAKER, "scaffolds", "scaf*.pdb")))
    print("scaffolds available: %d" % n_scaf, flush=True)
    if n_scaf < 100:
        raise SystemExit("expected ~152 scaffolds, found %d" % n_scaf)
    print("all inputs present", flush=True)


def build_and_dock():
    print()
    print("=" * 72)
    print("BUILD TARGET, DOCK SCAFFOLDS, SHORTLIST POSES")
    print("=" * 72)
    py = sys.executable
    sh("%s -W ignore pmhc/build_states.py" % py, cwd=REPO)
    sh("%s -W ignore pmhc/dock_scaffolds.py %s"
       % (py, os.path.join(BAKER, "scaffolds")), cwd=REPO)
    sh("%s -W ignore pmhc/shortlist_poses.py" % py, cwd=REPO)

    sl = os.path.join(REPO, "pmhc", "results", "shortlist_poses.json")
    data = json.load(open(sl))
    poses = data["poses"]
    print()
    print("shortlisted poses: %d" % len(poses), flush=True)
    if not poses:
        raise SystemExit("no poses to redesign")
    return poses


def run_mpnn(poses):
    print()
    print("=" * 72)
    print("ProteinMPNN - redesigning 100% of each binder, chain A held fixed")
    print("=" * 72)
    pose_dir = os.path.join(REPO, "data", "pmhc", "docked")
    keep = {"%s_docked" % p["name"] for p in poses}
    work = os.path.join(W, "mpnn_in")
    if os.path.isdir(work):
        shutil.rmtree(work)
    os.makedirs(work)
    for f in glob.glob(os.path.join(pose_dir, "*_docked.pdb")):
        if os.path.basename(f)[:-4] in keep:
            shutil.copy(f, work)
    n_in = len(glob.glob(os.path.join(work, "*.pdb")))
    print("pose files staged for redesign: %d" % n_in, flush=True)
    if n_in != len(poses):
        raise SystemExit("staged %d files for %d poses - mismatch"
                         % (n_in, len(poses)))

    mpnn = os.path.join(BAKER, "software", "mpnn")
    weights = os.path.join(mpnn, "vanilla_model_weights")
    out_root = os.path.join(W, "mpnn_out")
    os.makedirs(out_root, exist_ok=True)
    log = os.path.join(W, "mpnn_log.txt")

    # --sampling_temp takes a space-separated list and generates num_seq_per_target
    # sequences at EACH temperature, so one invocation per pose covers all three.
    temps = " ".join(str(t) for t in TEMPERATURES)
    expected_per_pose = N_SEQ * len(TEMPERATURES)

    print()
    print("settings: --pdb_path_chains B, so chain B is designed and chain A is held")
    print("          fixed (ProteinMPNN derives the fixed list as 'all chains minus")
    print("          designed', verified in protein_mpnn_run.py lines 168-172)")
    print("          temperatures %s, %d sequences each = %d per pose"
          % (temps, N_SEQ, expected_per_pose))
    print("          omit_AAs=%r, weights v_48_020 (the published pipeline's)"
          % OMIT_AA)
    print("          expected total: %d sequences" % (expected_per_pose * n_in))
    print()

    cmd = (
        "%s %s/protein_mpnn_run.py"
        " --pdb_path_chains B"
        " --path_to_model_weights %s --model_name v_48_020"
        " --num_seq_per_target %d --sampling_temp \"%s\" --batch_size 1"
        " --omit_AAs %s --seed 37"
        " --out_folder %s --pdb_path %s"
    )
    # Child output goes to a FILE, never a pipe. On run 2 a 64 kB pipe buffer with an
    # undrained parent stalled a worker and cost ~2 GPU hours; I had explicitly
    # dismissed that risk at launch by counting only my own log lines.
    failed = []
    files = sorted(glob.glob(os.path.join(work, "*.pdb")))
    with open(log, "w") as lf:
        for i, f in enumerate(files, 1):
            c = cmd % (sys.executable, mpnn, weights, N_SEQ, temps, OMIT_AA,
                       out_root, f)
            if i == 1:
                print("first command, in full:")
                print("  %s" % c)
                print()
            r = subprocess.run(c, shell=True, stdout=lf,
                               stderr=subprocess.STDOUT, text=True)
            if r.returncode != 0:
                failed.append(os.path.basename(f))
            if i % 10 == 0 or i == len(files):
                print("  %d/%d poses redesigned" % (i, len(files)), flush=True)

    if failed:
        print()
        print("FAILED on %d pose(s): %s" % (len(failed), failed[:10]))
        print("last 2000 characters of the ProteinMPNN log:")
        txt = open(log).read()
        print(txt[-2000:])
        if len(failed) == len(files):
            raise SystemExit("ProteinMPNN failed on every pose")
    return out_root, [log]


def parse_fasta(path):
    """ProteinMPNN writes one FASTA per input PDB: first record is the native
    sequence, the rest are samples with score/temperature in the header."""
    recs = []
    name = None
    buf = []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                recs.append((name, "".join(buf)))
            name, buf = line[1:], []
        elif line:
            buf.append(line)
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def collect(out_root, poses):
    print()
    print("=" * 72)
    print("SCREEN - composition, proline/glycine gate, and parent identity")
    print("=" * 72)
    sys.path.insert(0, os.path.join(REPO, "pmhc"))
    import screen

    by_name = {"%s_docked" % p["name"]: p for p in poses}
    rows = []
    n_seen = 0
    n_cys = 0

    for fa in sorted(glob.glob(os.path.join(out_root, "seqs", "*.fa"))):
        stem = os.path.basename(fa)[:-3]
        pose = by_name.get(stem)
        if pose is None:
            continue
        recs = parse_fasta(fa)
        if not recs:
            continue
        native = recs[0][1]
        for header, seq in recs[1:]:
            n_seen += 1
            # ProteinMPNN emits the full designed chain; chain B only, since that is
            # the only chain we asked it to design.
            seq = seq.split("/")[-1] if "/" in seq else seq
            if "C" in seq:
                n_cys += 1
                continue
            r = screen.report(seq)
            r["name"] = stem
            r["parent"] = pose["name"]
            r["seq"] = seq
            r["header"] = header
            r["parent_identity"] = screen.identity(seq, pose["scaffold_seq"])
            r["native_identity"] = screen.identity(seq, native)
            for k in ("n_mut_contacts", "n_pep_contacts", "n_mhc_contacts",
                      "pep_fraction", "binder_len", "align_rmsd"):
                r[k] = pose[k]
            # score from the FASTA header, e.g. "...,score=1.234,..."
            r["mpnn_score"] = None
            for part in header.split(","):
                if part.strip().startswith("score="):
                    try:
                        r["mpnn_score"] = float(part.split("=")[1])
                    except ValueError:
                        pass
            rows.append(r)

    print("sequences generated        : %d" % n_seen)
    print("dropped for cysteine       : %d  (omit_AAs=%r %s)"
          % (n_cys, OMIT_AA,
             "WORKED" if n_cys == 0 else "WAS IGNORED - checked, not trusted"))
    print("sequences screened         : %d" % len(rows))
    print()
    if not rows:
        raise SystemExit("no sequences survived collection")

    from collections import Counter
    tiers = Counter(r["tier"] for r in rows)
    print("--- tier counts (denominator %d) ---" % len(rows))
    for t in ("good", "marginal", "fail", "fold-breaker"):
        print("    %-14s %5d" % (t, tiers.get(t, 0)))
    print()
    print("THE COMPARISON THAT MATTERS: the previous approach produced 720 sequences")
    print("across 15 hallucinated backbones with ZERO good and ZERO marginal.")
    print()

    ident = sorted(r["parent_identity"] for r in rows)
    print("--- identity to the parent scaffold (the de novo declaration) ---")
    print("    min %.2f  median %.2f  max %.2f"
          % (ident[0], ident[len(ident) // 2], ident[-1]))
    print("    Every design's parent and identity is recorded and will be stated in")
    print("    the submission methods text. Track 3 says 'de novo only; modifying an")
    print("    existing binder is not allowed'. These backbones are de novo designs")
    print("    built for DIFFERENT pMHC targets, and 100% of the sequence is")
    print("    rewritten here. We declare it rather than hide it.")
    print()
    return rows


def pick(rows):
    """Keep the best few per pose, preferring good over marginal, then MPNN score."""
    print("--- selection rule, stated before the numbers are read ---")
    print("    1. tier 'good' before 'marginal'; 'fail' and 'fold-breaker' discarded")
    print("    2. within a tier, lower ProteinMPNN score (better sequence recovery)")
    print("    3. at most %d per pose, so the final set is not one backbone repeated"
          % KEEP_PER_POSE)
    print()
    rank = {"good": 0, "marginal": 1}
    usable = [r for r in rows if r["tier"] in rank]
    if not usable:
        raise SystemExit("nothing reached 'good' or 'marginal' - this approach has "
                         "failed the same test the last one did")

    by_pose = {}
    for r in usable:
        by_pose.setdefault(r["parent"], []).append(r)
    out = []
    for parent, rs in by_pose.items():
        rs.sort(key=lambda r: (rank[r["tier"]],
                               r["mpnn_score"] if r["mpnn_score"] is not None else 9e9,
                               r["distance"]))
        out.extend(rs[:KEEP_PER_POSE])

    out.sort(key=lambda r: (rank[r["tier"]], -r["n_mut_contacts"],
                            -r["pep_fraction"],
                            r["mpnn_score"] if r["mpnn_score"] is not None else 9e9))
    print("poses contributing: %d" % len(by_pose))
    print("candidates carried forward: %d" % len(out))
    return out


def main():
    setup()
    poses = build_and_dock()
    out_root, logs = run_mpnn(poses)
    rows = collect(out_root, poses)
    picked = pick(rows)

    print()
    print("=" * 72)
    print("CANDIDATES - printed in full because /kaggle/working does not persist")
    print("=" * 72)
    print("%-22s %-8s %5s %4s %6s %6s %6s %6s"
          % ("name", "tier", "len", "mut", "pepfr", "ident", "score", "dist"))
    for r in picked:
        print("%-22s %-8s %5d %4d %6.2f %6.2f %6s %6.2f"
              % (r["name"][:22], r["tier"], r["length"], r["n_mut_contacts"],
                 r["pep_fraction"], r["parent_identity"],
                 ("%.3f" % r["mpnn_score"]) if r["mpnn_score"] is not None else "-",
                 r["distance"]))
    print()
    print("--- sequences, TSV, copy this block ---")
    print("name\tparent\ttier\tlen\tmut\tpepfr\tident\tmpnn\tseq")
    for i, r in enumerate(picked, 1):
        print("%s_%03d\t%s\t%s\t%d\t%d\t%.3f\t%.3f\t%s\t%s"
              % (r["parent"], i, r["parent"], r["tier"], r["length"],
                 r["n_mut_contacts"], r["pep_fraction"], r["parent_identity"],
                 ("%.4f" % r["mpnn_score"]) if r["mpnn_score"] is not None else "NA",
                 r["seq"]))

    save = os.path.join(W, "pmhc_scafmpnn_candidates.json")
    json.dump(picked, open(save, "w"), indent=1)
    print()
    print("saved %s (will not survive a committed run - the log above is the record)"
          % save)
    print()
    print("NOT DONE AND NOT VERIFIED:")
    print("  - No structure prediction has been run on any of these sequences.")
    print("  - No specificity measurement yet: whether a design prefers Asp6 over")
    print("    Gly6 is the next step and is the thing we are actually ranked on.")
    print("  - Composition tier says 'could fold and could be expressed'. It says")
    print("    nothing about binding.")


main()
