"""Check the extracted screen module still behaves like the copies it replaces.

pmhc/screen.py was lifted out of four Kaggle scripts that each held their own copy. An
extraction that quietly changes a threshold would silently alter every result from here
on, so the cases below pin the behaviour that mattered in previous challenges.

Run:  python pmhc/test_screen_module.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import screen  # noqa: E402

# The one TfR1 candidate that survived both the composition screen and the proline gate
# (pmhc/../tfr1/results/mpnn2_candidates.tsv, row tfr1b_e559_003_r02). Recorded in
# checkpoint 24 as tier "good".
TFR1_SURVIVOR = ("SPTKTYVSIEKTKSKGGKGPKPKPKGRKKLVEKLAERLAEREAARAAAEAAARAAALAAIPA")

CASES = [
    # (label, sequence, expected tier or None for "just must not crash")
    ("polyproline rod", "P" * 40, "fold-breaker"),
    ("polyglycine", "G" * 40, "fold-breaker"),
    ("pro/gly mix over the joint gate",
     "PGPGPGPGPGAAAAAAAAAALLLLLLLLLLEEEEEEEEEE", "fold-breaker"),
    ("all-hydrophobic aggregate", "LLLLLVVVVVIIIIIFFFFFAAAAAWWWWWMMMMMYYYYY", None),
    ("all-charged", "DEKRDEKRDEKRDEKRDEKRDEKRDEKRDEKRDEKRDEKR", None),
]

FAILURES = []


def check(label, got, want):
    ok = got == want
    print("    %-40s %-14s %s" % (label, got, "ok" if ok else "EXPECTED %s" % want))
    if not ok:
        FAILURES.append(label)


def main():
    print("--- the proline/glycine gate must fire first ---")
    for label, seq, want in CASES:
        r = screen.report(seq)
        if want is None:
            print("    %-40s %-14s (no expectation, must not crash)"
                  % (label, r["tier"]))
            continue
        check(label, r["tier"], want)
    print()

    print("--- a hydrophobic aggregate must NOT pass as good ---")
    r = screen.report("LLLLLVVVVVIIIIIFFFFFAAAAAWWWWWMMMMMYYYYY")
    print("    hydrophobic fraction %.2f, tier %s" % (r["hydrophobic_frac"], r["tier"]))
    if r["tier"] == "good":
        print("    FAIL: the two-sided hydrophobicity limit is not being applied")
        FAILURES.append("hydrophobic aggregate passed as good")
    print()

    print("--- the one real TfR1 survivor should still pass ---")
    r = screen.report(TFR1_SURVIVOR)
    print("    length %d, %.0f%% charged, %.0f%% hydrophobic, %.0f%% Pro, %.0f%% Gly,"
          " kd %+.2f" % (r["length"], 100 * r["charged_frac"],
                         100 * r["hydrophobic_frac"], 100 * r["pro_frac"],
                         100 * r["gly_frac"], r["kd"]))
    print("    tier: %s %s" % (r["tier"], r["why"]))
    if r["tier"] == "fold-breaker":
        print("    FAIL: this sequence passed the gate when the gate was added in")
        print("    checkpoint 24; the thresholds have drifted during extraction")
        FAILURES.append("TfR1 survivor now fails")
    print()

    print("--- identity, used for the de novo declaration ---")
    cases = [("identical", "ACDEFG", "ACDEFG", 1.0),
             ("completely different", "AAAAAA", "LLLLLL", 0.0),
             ("half", "AAALLL", "AAAEEE", 0.5)]
    for label, a, b, want in cases:
        got = screen.identity(a, b)
        check("%s (%s vs %s)" % (label, a, b), got, want)
    print()

    if FAILURES:
        print("TEST FAILED: %s" % "; ".join(FAILURES))
        return 1
    print("PASS - the gate fires first, the two-sided limits hold, and the known")
    print("survivor still survives")
    return 0


if __name__ == "__main__":
    sys.exit(main())
