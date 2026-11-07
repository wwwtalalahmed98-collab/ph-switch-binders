"""Check that threading a designed sequence onto a pose writes what we think it does.

specificity_screen.py takes a docked pose and rewrites chain B's residue names to carry
a ProteinMPNN-designed sequence, keeping only backbone atoms. Every number the
specificity screen produces depends on that file being correct. If chain B came out
shifted by one, or chain A were altered, or the peptide moved, the screen would still
run and still print confident numbers - which is precisely the failure mode that cost
this project three rounds of measurement on earlier challenges.

So: write a file, read it back, and check it residue by residue.

Run:  python pmhc/test_threading.py
"""
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import specificity_screen as ss  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCKED = os.path.join(ROOT, "data", "pmhc", "docked")
TMP = os.path.join(ROOT, "data", "pmhc", "spec", "_test")

T3 = {v: k for k, v in ss.T1.items()}
FAILURES = []


def read_chains(path):
    """(chain -> sequence, chain -> atom names seen)."""
    seqs, atoms, seen = {}, {}, {}
    for line in open(path):
        if not line.startswith("ATOM"):
            continue
        ch = line[21]
        num = int(line[22:26])
        name = line[12:16].strip()
        atoms.setdefault(ch, set()).add(name)
        key = (ch, num)
        if key not in seen:
            seen[key] = T3[line[17:20].strip()]
            seqs.setdefault(ch, []).append((num, seen[key]))
    return ({c: "".join(x[1] for x in sorted(v)) for c, v in seqs.items()}, atoms)


def fail(msg):
    print("    FAIL: %s" % msg)
    FAILURES.append(msg)


def main():
    poses = sorted(glob.glob(os.path.join(DOCKED, "*_docked.pdb")))
    if not poses:
        raise SystemExit("run pmhc/dock_scaffolds.py first")
    pose = poses[0]
    os.makedirs(TMP, exist_ok=True)
    print("pose under test: %s" % os.path.basename(pose))

    orig, _ = read_chains(pose)
    n_binder = len(orig["B"])
    print("    original: chain A %d res, chain B %d res" % (len(orig["A"]), n_binder))
    print()

    # A designed sequence that is deliberately nothing like the original, and that
    # uses a distinct residue at each end so an off-by-one would be visible.
    designed = ("W" + "AEKLAEKLAEKL" * 20)[:n_binder]
    designed = designed[:-1] + "H"
    print("--- threading a %d-residue designed sequence onto chain B ---"
          % len(designed))
    out = ss.thread_sequence(pose, designed, os.path.join(TMP, "threaded.pdb"))
    got, atoms = read_chains(out)

    if got.get("B") == designed:
        print("    chain B reads back exactly as designed (including both termini)")
    else:
        fail("chain B is %r, expected %r" % (got.get("B", "")[:40], designed[:40]))

    if got.get("A") == orig["A"]:
        print("    chain A is unchanged (%d residues)" % len(got["A"]))
    else:
        fail("chain A changed during threading")

    if got.get("A", "").endswith(ss.PEP):
        print("    chain A still ends in the peptide %s" % ss.PEP)
    else:
        fail("chain A no longer ends in %s" % ss.PEP)

    # The index the whole screen hangs on.
    a = got.get("A", "")
    if len(a) > ss.MUT_IDX and a[ss.MUT_IDX] == "D":
        print("    index %d is 'D', the G12D aspartate" % ss.MUT_IDX)
    else:
        fail("index %d is %r, not the aspartate"
             % (ss.MUT_IDX, a[ss.MUT_IDX] if len(a) > ss.MUT_IDX else None))

    extra = set()
    for ch, names in atoms.items():
        extra |= names - set(ss.BACKBONE)
    if not extra:
        print("    only backbone atoms written (%s), as ProteinMPNN requires"
              % ", ".join(ss.BACKBONE))
    else:
        fail("non-backbone atoms present: %s" % sorted(extra)[:8])
    print()

    print("--- the no-binder control ---")
    bare = ss.thread_sequence(pose, "", os.path.join(TMP, "bare.pdb"),
                              drop_binder=True)
    bg, _ = read_chains(bare)
    if "B" in bg:
        fail("the bare control still contains chain B")
    else:
        print("    chain B absent, as intended")
    if bg.get("A") == orig["A"]:
        print("    chain A identical to the pose (%d residues)" % len(bg["A"]))
    else:
        fail("chain A differs in the bare control")
    print()

    print("--- a too-short sequence must be refused, not silently truncated ---")
    try:
        ss.thread_sequence(pose, "AAA", os.path.join(TMP, "short.pdb"))
        fail("a 3-residue sequence was accepted for a %d-residue backbone" % n_binder)
    except SystemExit:
        print("    refused, with an explanation")
    print()

    if FAILURES:
        print("TEST FAILED (%d): %s" % (len(FAILURES), "; ".join(FAILURES)))
        return 1
    print("PASS - the threaded file carries the designed sequence on chain B, leaves")
    print("the target untouched, keeps the mutation where the screen looks for it,")
    print("and refuses a mismatched sequence")
    return 0


if __name__ == "__main__":
    sys.exit(main())
