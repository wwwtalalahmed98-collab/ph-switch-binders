"""Run Adaptyv's own pre-submission screen on our designs, instead of assuming we pass.

WHAT THIS IS. The competition links its submission gate as source-available code at
github.com/adaptyvbio/binder-prescreen. It answers two questions from sequence alone:

  1. Is this already public?  MMseqs2 against PDB, SwissProt, PLAbDab (+ nano),
     Thera-SAbDab, THPdb, published Proteinbase designs, and 11.2M USPTO patent
     sequences.
  2. Is it a known binder of the target?  Whole-sequence and binding-region comparison
     against a curated set for that target.

For `--target kras-pmhc` that curated set is 144 references and they are ALL T-cell
receptors: "No TCR-mimic antibody or designed-binder sequences were available for this
pMHC." Our designs are miniproteins, so the target arm should be quiet. The arm that
can actually bite us is prior art - whether a ProteinMPNN redesign at 25-35% identity
to a published scaffold still reads as that scaffold.

WHAT IT DOES NOT DECIDE. Scaffold reuse. The README is explicit: "This pipeline does
not judge whether the scaffold is de novo. Reusing a published framework is allowed in
these competitions", and framework-reuse observations "never disqualify". The separate
novelty LEVEL, computed by ProteinTyper on structure as well as sequence, is not
checked here and remains the open risk in this submission.

WHY ON KAGGLE. It needs MMseqs2 and uv, which are painful on Windows, it is CPU-only
so it costs no GPU quota, and it runs third-party code - which belongs in a disposable
VM rather than on the project machine.

TWO PASSES, CHEAP FIRST. Pass 1 is --skip-prior-art: no downloads, ~25 s, exercises the
target arm and the verbatim-fragment check. Pass 2 builds the reference databases minus
USPTO (the patent arm alone is 42 minutes and 1.7 GB) and runs the full screen. If
pass 2 fails for any reason, pass 1's verdicts have already been printed.
"""
import glob
import os
import subprocess
import sys
import time

W = "/kaggle/working"
WORK = "/tmp/prescreen"
REPO = os.path.join(WORK, "binder-prescreen")
DBROOT = os.path.join(WORK, "dbs")
TARGET = "kras-pmhc"

# every arm except uspto, which is 1.7 GB and ~42 min on its own
ARMS = "thpdb therasabdab plabdab_nano proteinbase_public plabdab pdb swissprot"

T0 = time.time()


def mins():
    return (time.time() - T0) / 60.0


def sh(cmd, check=False, stream=False, tail=6000, cwd=None):
    print("$ %s" % cmd, flush=True)
    if stream:
        r = subprocess.run(cmd, shell=True, cwd=cwd)
        out = ""
    else:
        r = subprocess.run(cmd, shell=True, cwd=cwd, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, text=True)
        out = r.stdout or ""
        if out.strip():
            print(out[-tail:], flush=True)
    sys.stdout.flush()
    if check and r.returncode != 0:
        raise SystemExit("failed (%d): %s" % (r.returncode, cmd))
    return out


def banner(s):
    print()
    print("=" * 72)
    print("%s   [t+%.1f min]" % (s, mins()))
    print("=" * 72, flush=True)


def find_candidates():
    banner("INPUT")
    found = sorted(glob.glob("/kaggle/input/**/*.csv", recursive=True))
    print("csv files under /kaggle/input: %d" % len(found))
    for f in found:
        print("    %s" % f)
    if not found:
        for root, dirs, files in os.walk("/kaggle/input"):
            print("      %s -> dirs %s files %s" % (root, dirs[:6], files[:6]))
        raise SystemExit("no candidate csv attached")
    src = found[0]

    # The prescreen wants name,sequence[,molecule_class]. Ours may carry extra
    # columns under different names, so normalise rather than hope.
    import csv as _csv
    rows = list(_csv.DictReader(open(src)))
    if not rows:
        raise SystemExit("%s is empty" % src)
    cols = list(rows[0].keys())
    print()
    print("columns: %s" % cols)
    scol = next((c for c in cols if c.lower() in ("sequence", "seq")), None)
    ncol = next((c for c in cols if c.lower() in ("name", "id", "design",
                                                  "submit_name")), None)
    if scol is None:
        raise SystemExit("no sequence column among %s" % cols)
    out = os.path.join(WORK, "submissions.csv")
    os.makedirs(WORK, exist_ok=True)
    seen = set()
    n = 0
    with open(out, "w", newline="") as fh:
        w = _csv.writer(fh)
        w.writerow(["name", "sequence", "molecule_class"])
        for i, r in enumerate(rows, 1):
            s = (r[scol] or "").strip().upper()
            if not s or s in seen:
                continue
            seen.add(s)
            nm = (r[ncol].strip() if ncol and r.get(ncol) else "design_%03d" % i)
            w.writerow([nm, s, "miniprotein"])
            n += 1
    print("wrote %s with %d unique sequences" % (out, n))
    return out


def setup():
    banner("SETUP")
    os.makedirs(WORK, exist_ok=True)
    sh("pip install -q uv", stream=True)
    if not os.path.isdir(REPO):
        sh("git clone --depth 1 https://github.com/adaptyvbio/binder-prescreen.git %s"
           % REPO, check=True)
    # MMseqs2 static binary: no conda, no compile.
    mm = os.path.join(WORK, "mmseqs", "bin")
    if not os.path.exists(os.path.join(mm, "mmseqs")):
        sh("cd %s && curl -sSL https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz "
           "| tar xz" % WORK, stream=True)
    os.environ["PATH"] = mm + ":" + os.environ.get("PATH", "")
    sh("mmseqs version", check=False)
    sh("uv venv && uv pip install -e .", cwd=REPO, stream=True)
    print()
    sh("uv run prescreen --list-targets", cwd=REPO)


def pass1(subs):
    banner("PASS 1 - target arm only, no databases")
    out = os.path.join(WORK, "report_fast")
    sh("uv run prescreen %s --skip-prior-art --target %s -o %s"
       % (subs, TARGET, out), cwd=REPO, stream=True)
    show(out, "report_fast")


def pass2(subs):
    banner("PASS 2 - full screen, all arms except USPTO patents")
    os.makedirs(DBROOT, exist_ok=True)
    sh("bash scripts/build_reference_dbs.sh %s %s" % (DBROOT, ARMS),
       cwd=REPO, stream=True)
    os.environ["PRESCREEN_DB_ROOT"] = DBROOT
    sh("uv run python -c \"from prescreen import refdb; print(refdb.available())\"",
       cwd=REPO)
    out = os.path.join(WORK, "report_full")
    sh("PRESCREEN_DB_ROOT=%s uv run prescreen %s --target %s --no-patent -o %s"
       % (DBROOT, subs, TARGET, out), cwd=REPO, stream=True)
    show(out, "report_full")


def show(outdir, label):
    import csv as _csv
    import shutil
    p = os.path.join(outdir, "report.csv")
    if not os.path.exists(p):
        print("no report.csv in %s; contents: %s"
              % (outdir, os.listdir(outdir) if os.path.isdir(outdir) else "missing"))
        return
    rows = list(_csv.DictReader(open(p)))
    print()
    print("--- %s: %d rows ---" % (label, len(rows)))
    vc = {}
    for r in rows:
        v = r.get("verdict", "?")
        vc[v] = vc.get(v, 0) + 1
    for v, n in sorted(vc.items(), key=lambda kv: -kv[1]):
        flag = "" if v == "pass" else "   <-- not a pass"
        print("    %-32s %3d%s" % (v, n, flag))
    bad = [r for r in rows if r.get("verdict") != "pass"]
    if bad:
        print()
        print("    designs that did NOT pass:")
        for r in bad:
            print("      %-22s %-28s %s"
                  % (r.get("name", "?")[:22], r.get("verdict", "?"),
                     (r.get("prior_art_hit") or r.get("target_hit") or "")[:40]))
    for name in ("report.csv", "report.json"):
        s = os.path.join(outdir, name)
        if os.path.exists(s):
            shutil.copy(s, os.path.join(W, "%s_%s" % (label, name)))
            print("    saved %s_%s" % (label, name))


def main():
    subs = find_candidates()
    setup()
    try:
        pass1(subs)
    except SystemExit:
        raise
    except Exception as e:                                   # noqa: BLE001
        print("pass 1 failed: %s %s" % (type(e).__name__, str(e)[:300]))
    try:
        pass2(subs)
    except Exception as e:                                   # noqa: BLE001
        print("pass 2 failed: %s %s" % (type(e).__name__, str(e)[:300]))
        print("pass 1's verdicts above still stand.")
    banner("DONE")
    print("elapsed %.1f min" % mins())
    return 0


if __name__ == "__main__":
    sys.exit(main())
