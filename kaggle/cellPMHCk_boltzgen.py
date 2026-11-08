"""BoltzGen generation for Challenge 5: ask for peptide contact instead of filtering for it.

WHY THIS EXISTS. Every previous attempt generated binder backbones against the whole
pMHC surface and then threw away the ones that ignored the peptide. Liu et al. say why
that is expensive: "only a small fraction of trajectories yield backbones that interact
extensively with the peptide while largely avoiding the MHC." BoltzGen's design
specification language states that requirement up front:

    binding:      peptide residues 4-8   (the G12D aspartate at 6 and its neighbours)
    not_binding:  heavy-chain 182-276    (the alpha-3 domain, away from the groove)

so the diffusion trajectory is steered toward the peptide rather than filtered after
the fact.

WHY GENERATION ONLY. The protein-anything protocol scores with ESMFold2, which needs a
CUDA 13 driver and ~33 GB of weights, and BoltzGen's own docs say "missing, non-finite,
or stale ESMFold2 scores raise errors; there is no implicit fallback to Boltz2 iPTM".
Kaggle's T4 nodes are not on a CUDA 13 driver. So this runs --steps design
inverse_folding and stops. Scoring is done afterwards by the three-state AF2 + ipSAE
screen we already have working, which has the side benefit that BoltzGen designs and
the scaffold-docked designs are ranked on ONE scale instead of two incomparable ones.

WHAT IT DOES NOT DO. It does not refold, score, rank or filter. It emits sequences.
Everything downstream is the existing pipeline.

PROBE, DO NOT ASSERT. Three Kaggle runs in a row died on an assumed path or an assumed
dependency. This one prints the driver version, the free disk, and `boltzgen run --help`
BEFORE committing to any flag, so a wrong assumption shows up in the log as a fact
rather than as a crash forty minutes later.
"""
import glob
import json
import os
import re
import shutil
import subprocess
import sys
import time

W = "/kaggle/working"
REPO = os.path.join(W, "ph-switch-binders")
OUR_REPO = "https://github.com/wwwtalalahmed98-collab/ph-switch-binders.git"

# Weights are large and must not land in /kaggle/working, which is capped and is
# persisted as notebook output.
SCRATCH = "/tmp/bg"
CACHE = os.path.join(SCRATCH, "boltz_cache")
SPEC_DIR = os.path.join(SCRATCH, "spec")
OUT_DIR = os.path.join(SCRATCH, "run")

NUM_DESIGNS = 120          # halved from 240: adding B2M makes each design slower
                           # and the deadline is 9 Nov 11:59 UTC
BINDER_LEN = "60..110"     # our docked scaffolds ran 54-118
BUDGET_MIN = 400

PEPTIDE_CHAIN = "C"
MHC_CHAIN = "A"
# Run 4 left beta-2 microglobulin out, and all 172 of its designs came back AS
# beta-2 microglobulin: the model filled the obvious hole instead of gripping the
# peptide. B2M is now part of the target, so there is no hole left to fill.
B2M_CHAIN = "B"
BIND_PEPTIDE_AUTHOR = (4, 8)        # G12D sits at peptide position 6
NO_BIND_MHC_AUTHOR = (182, 276)     # alpha-3 domain, remote from the groove

T3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
      "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
      "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
      "TYR": "Y", "VAL": "V"}

T0 = time.time()


def mins():
    return (time.time() - T0) / 60.0


def sh(cmd, check=True, stream=False, tail=4000):
    print("$ %s" % cmd, flush=True)
    if stream:
        r = subprocess.run(cmd, shell=True)
        out = ""
    else:
        r = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, text=True)
        out = r.stdout or ""
        if out.strip():
            print(out[-tail:], flush=True)
    sys.stdout.flush()
    if check and r.returncode != 0:
        raise SystemExit("command failed (%d): %s" % (r.returncode, cmd))
    return out


def banner(s):
    print()
    print("=" * 72)
    print("%s   [t+%.1f min]" % (s, mins()))
    print("=" * 72, flush=True)


# --------------------------------------------------------------------------
def probe():
    banner("PROBE - what this machine actually is")
    sh("nvidia-smi --query-gpu=name,driver_version,memory.total "
       "--format=csv,noheader", check=False)
    print("python: %s" % sys.version.split()[0])
    sh("df -h /tmp /kaggle/working | sed -n '1,8p'", check=False)
    # CUDA 13 needs driver >= 580. Report it rather than guess later.
    drv = sh("nvidia-smi --query-gpu=driver_version --format=csv,noheader",
             check=False).strip().splitlines()
    if drv:
        try:
            major = int(drv[0].split(".")[0])
            print()
            print("driver major %d -> ESMFold2 (needs CUDA 13, driver >=580) is %s"
                  % (major, "possible" if major >= 580 else "NOT available"))
            print("this script does not use ESMFold2 either way")
        except ValueError:
            pass


def install():
    banner("INSTALL boltzgen")
    try:
        import boltzgen  # noqa: F401
        print("boltzgen already present")
    except ImportError:
        sh("pip install -q boltzgen", stream=True)
    sh("python -c \"import boltzgen, sys; print('boltzgen', "
       "getattr(boltzgen,'__version__','?'))\"", check=False)
    print()
    print("--- boltzgen run --help (read the real flags, do not assume them) ---")
    sh("boltzgen run --help", check=False, tail=12000)


def get_target():
    banner("TARGET")
    if not os.path.isdir(REPO):
        sh("git clone --depth 1 %s %s" % (OUR_REPO, REPO))
    # data/ is gitignored, so the crystal structures are NOT in the clone. The first
    # attempt assumed they were and died on a missing 7OW6.cif after installing
    # boltzgen for five minutes. Fetch them the same way the AF2 script does.
    import urllib.request
    d = os.path.join(REPO, "data", "pmhc")
    os.makedirs(d, exist_ok=True)
    for pdb_id in ("7OW6", "7OW5"):
        p = os.path.join(d, "%s.cif" % pdb_id)
        if not os.path.exists(p):
            urllib.request.urlretrieve(
                "https://files.rcsb.org/download/%s.cif" % pdb_id, p)
        print("%s.cif %d bytes" % (pdb_id, os.path.getsize(p)))
    sh("cd %s && python -W ignore pmhc/build_states.py" % REPO, check=False)
    sh("cd %s && python -W ignore pmhc/build_counter_states.py" % REPO, check=False)
    src = os.path.join(REPO, "data", "pmhc", "states", "target_state_g12d.pdb")
    if not os.path.exists(src):
        raise SystemExit("no target state at %s" % src)
    os.makedirs(SPEC_DIR, exist_ok=True)
    dst = os.path.join(SPEC_DIR, "target.pdb")
    shutil.copy(src, dst)

    # Residue numbers in order of appearance, per chain. BoltzGen indexes residues by
    # ORDINAL POSITION in the chain, not by the author numbering in the file: the first
    # attempt wrote the alpha-3 domain as its author range 182..276 and got
    # "Specified end 182..276 is higher than the length of the chain", because chain A
    # is numbered 2..276 but holds 275 residues. So build the author->ordinal map here
    # and derive every index from it instead of writing numbers by hand.
    order, seqs = {}, {}
    for line in open(dst):
        if line.startswith(("ATOM", "HETATM")):
            c, num = line[21], int(line[22:26])
            lst = order.setdefault(c, [])
            if not lst or lst[-1] != num:
                if num not in lst:
                    lst.append(num)
                    seqs.setdefault(c, []).append(
                        T3.get(line[17:20].strip().upper(), "X"))
    print("chains in target.pdb:")
    for c in sorted(order):
        print("    %s  %d residues  author %d..%d  -> ordinals 1..%d"
              % (c, len(order[c]), order[c][0], order[c][-1], len(order[c])))
    for need in (MHC_CHAIN, PEPTIDE_CHAIN, B2M_CHAIN):
        if need not in order:
            raise SystemExit("chain %s missing from the target" % need)
    if len(order[PEPTIDE_CHAIN]) != 10:
        raise SystemExit("peptide chain %s has %d residues, expected 10"
                         % (PEPTIDE_CHAIN, len(order[PEPTIDE_CHAIN])))
    # Only the chains we actually feed BoltzGen count as target sequences. The
    # harvester subtracts these to find the designed chain, so beta-2
    # microglobulin, which is in the file but excluded from the spec, must NOT
    # be in this set.
    target_seqs = ["".join(seqs[c])
                   for c in (MHC_CHAIN, PEPTIDE_CHAIN, B2M_CHAIN)]
    for c, t in zip((MHC_CHAIN, PEPTIDE_CHAIN, B2M_CHAIN), target_seqs):
        print("target chain %s: %d aa" % (c, len(t)))
    if target_seqs[1] != "VVVGADGVGK":
        raise SystemExit("peptide reads %s, expected VVVGADGVGK"
                         % target_seqs[1])
    return dst, order, target_seqs


def ordinal_range(order, chain, author_lo, author_hi):
    """Author residue range -> 1-based ordinal range, clipped to what exists."""
    nums = order[chain]
    idx = [i + 1 for i, n in enumerate(nums) if author_lo <= n <= author_hi]
    if not idx:
        raise SystemExit("chain %s has no residues in author range %d..%d"
                         % (chain, author_lo, author_hi))
    return idx[0], idx[-1]


def write_spec(target_pdb, order):
    banner("DESIGN SPECIFICATION")
    bind_lo, bind_hi = ordinal_range(order, PEPTIDE_CHAIN, *BIND_PEPTIDE_AUTHOR)
    nb_lo, nb_hi = ordinal_range(order, MHC_CHAIN, *NO_BIND_MHC_AUTHOR)
    no_bind_b2m = "1..%d" % len(order[B2M_CHAIN])
    bind = "%d..%d" % (bind_lo, bind_hi)
    no_bind = "%d..%d" % (nb_lo, nb_hi)
    print("peptide bind  : author %d..%d -> ordinals %s"
          % (BIND_PEPTIDE_AUTHOR[0], BIND_PEPTIDE_AUTHOR[1], bind))
    print("MHC not_bind  : author %d..%d -> ordinals %s"
          % (NO_BIND_MHC_AUTHOR[0], NO_BIND_MHC_AUTHOR[1], no_bind))
    print("B2M not_bind : whole chain  -> ordinals %s" % no_bind_b2m)
    # beta-2 microglobulin is now IN the target. Leaving it out was the whole failure
    # of run 4: an HLA heavy chain with a B2M-shaped hole in it has one overwhelmingly
    # likely completion, and the model produced that completion 172 times out of 172.
    # It is marked not_binding over its full length, as is the alpha-3 domain, so the
    # underside is covered rather than offered as a landing site.
    spec = """entities:
  - protein:
      id: D
      sequence: %s

  - file:
      path: target.pdb
      include:
        - chain:
            id: %s
        - chain:
            id: %s
        - chain:
            id: %s
      binding_types:
        - chain:
            id: %s
            binding: "%s"
        - chain:
            id: %s
            not_binding: "%s"
        - chain:
            id: %s
            not_binding: "%s"
""" % (BINDER_LEN, MHC_CHAIN, PEPTIDE_CHAIN, B2M_CHAIN,
       PEPTIDE_CHAIN, bind, MHC_CHAIN, no_bind, B2M_CHAIN, no_bind_b2m)
    path = os.path.join(SPEC_DIR, "pmhc.yaml")
    open(path, "w").write(spec)
    print(spec)
    print("--- boltzgen check ---")
    sh("cd %s && boltzgen check pmhc.yaml" % SPEC_DIR, check=False)
    return path


def _run(steps, label):
    os.makedirs(CACHE, exist_ok=True)
    env = "HF_HOME=%s/hf UV_CACHE_DIR=%s/uv " % (SCRATCH, SCRATCH)
    cmd = ("cd %s && %s boltzgen run pmhc.yaml --output %s "
           "--protocol protein-anything --steps %s "
           "--num_designs %d --cache %s --reuse"
           % (SPEC_DIR, env, OUT_DIR, steps, NUM_DESIGNS, CACHE))
    sh(cmd, check=False, stream=True)
    print()
    print("[t+%.1f min] %s returned" % (mins(), label))
    sh("find %s -maxdepth 2 -type d | head -40" % OUT_DIR, check=False)
    sh("df -h /tmp | tail -1", check=False)


def run_design(spec):
    # STAGE 1 banks the sequences. Everything after it is a bonus, so a failure in
    # scoring cannot cost us the GPU time already spent generating.
    banner("STAGE 1 - design + inverse_folding (banks the sequences)")
    _run("design inverse_folding", "stage 1")


def run_scoring():
    # The driver on Kaggle's T4 nodes turned out to be 580, so ESMFold2 IS available
    # - the opposite of what I assumed when writing this. Its ipSAE ranking is
    # BoltzGen's own quality filter and is worth having. But it wants ~33 GB of
    # weights, so it runs second, with --reuse, after stage 1 is safely on disk.
    banner("STAGE 2 - folding + ESMFold2 scoring + analysis + filtering (optional)")
    _run("folding esmfold2_scoring analysis filtering", "stage 2")
    for p in sorted(glob.glob(os.path.join(OUT_DIR, "final_ranked_designs", "*.csv"))):
        print()
        print("--- %s ---" % os.path.basename(p))
        try:
            with open(p) as fh:
                for n, line in enumerate(fh):
                    if n > 25:
                        print("    ...")
                        break
                    print("    " + line.rstrip()[:200])
            shutil.copy(p, os.path.join(W, os.path.basename(p)))
        except OSError as e:
            print("could not read %s: %s" % (p, e))


# --------------------------------------------------------------------------
def seqs_from_cif(path):
    """Pull one sequence per chain out of an mmCIF, by CA atoms, in file order.

    Written by hand rather than with gemmi because gemmi is an optional dependency
    here and a missing import at harvest time would throw away the whole GPU run.
    """
    hdr, rows, inloop = [], [], False
    for line in open(path, errors="replace"):
        s = line.strip()
        if s.startswith("_atom_site."):
            hdr.append(s)
            inloop = True
            continue
        if inloop:
            if s.startswith(("ATOM", "HETATM")):
                rows.append(s.split())
            elif s.startswith("#") and rows:
                break
    if not hdr or not rows:
        return {}
    idx = {k: i for i, k in enumerate(hdr)}

    def col(*names):
        for n in names:
            k = "_atom_site." + n
            if k in idx:
                return idx[k]
        return None

    ca = col("label_atom_id", "auth_atom_id")
    comp = col("label_comp_id", "auth_comp_id")
    asym = col("auth_asym_id", "label_asym_id")
    seqid = col("auth_seq_id", "label_seq_id")
    if None in (ca, comp, asym, seqid):
        return {}
    out, seen = {}, set()
    for f in rows:
        if len(f) <= max(ca, comp, asym, seqid):
            continue
        if f[ca].strip('"') != "CA":
            continue
        key = (f[asym], f[seqid])
        if key in seen:
            continue
        seen.add(key)
        out.setdefault(f[asym], []).append(T3.get(f[comp].upper(), "X"))
    return {k: "".join(v) for k, v in out.items()}


def pick_design(ch, targets):
    """Which chain is the design? The one that is not a target chain.

    Identified by sequence, never by chain letter: BoltzGen does not carry the
    ids from the input yaml through to its output, and keying on letters is
    exactly what made the first completed run harvest 0 of 48 files.
    """
    cands = []
    for k, v in ch.items():
        if v in targets:
            continue
        if any(v and (v in t or t in v) for t in targets):
            continue          # a crop or extension of a target chain
        cands.append((k, v))
    if len(cands) == 1:
        return cands[0]
    # fall back on the length band we asked the model for
    inband = [(k, v) for k, v in cands if 40 <= len(v) <= 160]
    if len(inband) == 1:
        return inband[0]
    return None, None


def identity_to_targets(seq, targets):
    """Highest ungapped identity to any target chain, over the best offset.

    Run 4 harvested 172 designs and every one of them was beta-2 microglobulin; I
    only found that out by reading the sequences by hand afterwards. This makes the
    machine say it. Cruder than an alignment, and deliberately conservative.
    """
    best = 0.0
    for t in targets:
        a, b = (seq, t) if len(seq) <= len(t) else (t, seq)
        if not a:
            continue
        for off in range(0, len(b) - len(a) + 1):
            same = sum(1 for i, ch in enumerate(a) if ch == b[off + i])
            best = max(best, same / float(len(a)))
    return best


def harvest(target_seqs=None):
    banner("HARVEST")
    pats = [os.path.join(OUT_DIR, "intermediate_designs_inverse_folded", "*.cif"),
            os.path.join(OUT_DIR, "**", "*inverse_folded*", "*.cif")]
    files = []
    for p in pats:
        files = sorted(glob.glob(p, recursive=True))
        if files:
            break
    if not files:
        print("no inverse-folded CIFs found. What is actually there:")
        for root, dirs, fs in os.walk(OUT_DIR):
            print("  %s -> %d files %s" % (root, len(fs), fs[:4]))
        return []
    print("inverse-folded CIFs: %d" % len(files))

    # UNCONDITIONAL DIAGNOSTIC. The previous run produced all 48 CIFs and then
    # harvested 0 of them, because the designed chain was identified by chain letter
    # and BoltzGen does not keep the letters from the input yaml. Nothing printed, so
    # 70 minutes of GPU bought one line saying "skipped 48". Whatever happens below,
    # the chain map of the first files is now in the log.
    print()
    print("--- chain map of the first %d file(s), as written ---"
          % min(2, len(files)))
    for f in files[:2]:
        ch = seqs_from_cif(f)
        print("  %s" % os.path.basename(f))
        if not ch:
            print("    PARSED NOTHING - first 3 lines of the file:")
            with open(f, errors="replace") as fh:
                for n, line in enumerate(fh):
                    if n >= 3:
                        break
                    print("      %s" % line.rstrip()[:100])
        for k, v in ch.items():
            print("    chain %-3s %4d aa  %s%s"
                  % (k, len(v), v[:40], "..." if len(v) > 40 else ""))

    # Identify the designed chain by SEQUENCE, not by letter: it is the chain that is
    # not one of the target chains we fed in. target_seqs comes from the same file the
    # yaml points at, so this does not depend on how BoltzGen labels anything.
    targets = set(target_seqs or [])
    if targets:
        print()
        print("target chains to exclude: %s"
              % ", ".join("%d aa" % len(t) for t in sorted(targets, key=len)))

    designs, bad, reasons = [], 0, {}
    for f in files:
        ch = seqs_from_cif(f)
        if not ch:
            bad += 1
            reasons["unparsable"] = reasons.get("unparsable", 0) + 1
            continue
        cid, seq = pick_design(ch, targets)
        if cid is None:
            bad += 1
            reasons["chain ambiguous"] = reasons.get("chain ambiguous", 0) + 1
            continue
        if "X" in seq:
            bad += 1
            reasons["contains X"] = reasons.get("contains X", 0) + 1
            continue
        if not (30 <= len(seq) <= 250):
            bad += 1
            reasons["length %d" % len(seq)] = reasons.get(
                "length %d" % len(seq), 0) + 1
            continue
        designs.append({"name": os.path.splitext(os.path.basename(f))[0],
                        "chain": cid, "seq": seq, "length": len(seq)})

    print()
    print("parsed %d designs, skipped %d" % (len(designs), bad))
    for k, v in sorted(reasons.items(), key=lambda kv: -kv[1]):
        print("    skipped, %-18s %d" % (k, v))
    if not designs:
        return []

    uniq = {}
    for d in designs:
        uniq.setdefault(d["seq"], d)
    designs = list(uniq.values())
    for d in designs:
        d["identity_to_target"] = round(identity_to_targets(d["seq"], targets), 3)
    lens = sorted(d["length"] for d in designs)
    ncys = sum(1 for d in designs if "C" in d["seq"])
    print("unique sequences : %d" % len(designs))
    print("length           : min %d  median %d  max %d"
          % (lens[0], lens[len(lens) // 2], lens[-1]))
    print("contain cysteine : %d of %d" % (ncys, len(designs)))
    # THE CHECK RUN 4 DID NOT HAVE. If the model has simply re-drawn one of the
    # chains we handed it, say so here rather than three hours later by eye.
    ids = sorted(d["identity_to_target"] for d in designs)
    echo = [d for d in designs if d["identity_to_target"] > 0.6]
    print("identity to a target chain: min %.2f  median %.2f  max %.2f"
          % (ids[0], ids[len(ids) // 2], ids[-1]))
    print("more than 60%% identical to a target chain: %d of %d"
          % (len(echo), len(designs)))
    if len(echo) > len(designs) // 2:
        print()
        print("    WARNING: most of this run is a copy of the target, not a binder.")
        print("    That is what happened in run 4. Do not screen these.")

    out = os.path.join(W, "boltzgen_designs.json")
    json.dump({"n": len(designs), "num_designs_requested": NUM_DESIGNS,
               "binder_len_spec": BINDER_LEN,
               "binding_author": {"chain": PEPTIDE_CHAIN,
                                  "res": BIND_PEPTIDE_AUTHOR},
               "not_binding_author": {"chain": MHC_CHAIN,
                                      "res": NO_BIND_MHC_AUTHOR},
               "designs": designs}, open(out, "w"), indent=1)
    print()
    print("wrote %s" % out)

    csvp = os.path.join(W, "boltzgen_designs.csv")
    with open(csvp, "w", newline="") as fh:
        fh.write("name,sequence,length,identity_to_target\n")
        for d in designs:
            fh.write("%s,%s,%d,%.3f\n" % (d["name"], d["seq"],
                                           d["length"],
                                           d["identity_to_target"]))
    print("wrote %s" % csvp)
    return designs


def main():
    probe()
    install()
    tgt, order, target_seqs = get_target()
    spec = write_spec(tgt, order)
    if mins() > BUDGET_MIN * 0.6:
        print("too little budget left to start design; stopping after setup")
        return 0
    run_design(spec)
    got = harvest(target_seqs)
    if got and mins() < BUDGET_MIN * 0.55:
        try:
            run_scoring()
        except Exception as e:                               # noqa: BLE001
            print("stage 2 failed: %s %s" % (type(e).__name__, str(e)[:300]))
            print("the sequences from stage 1 are already written and stand.")
    else:
        print("skipping stage 2: %s"
              % ("no designs to score" if not got else "not enough budget left"))
    banner("DONE")
    print("elapsed %.1f min" % mins())
    return 0


if __name__ == "__main__":
    sys.exit(main())
