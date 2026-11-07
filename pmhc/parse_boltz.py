"""Turn raw ipsae.py output into the ranking, with the arithmetic in the open.

The Kaggle notebook prints ipsae.py's rows verbatim rather than reducing them there,
so the step that actually decides the submission order happens here, in a file that can
be read and argued with, instead of being buried in a notebook log.

WHAT WE TAKE FROM EACH PREDICTION. ipsae.py emits one row per ordered chain pair. Our
complexes have three chains:

    A  MHC-I platform (179 aa)
    B  peptide, VVVGADGVGK (G12D) or VVVGAGGVGK (wild type)
    C  the designed binder

so the rows we care about are C-B (binder against peptide) and C-A (binder against
MHC). We read the ipSAE column, and because ipsae.py reports each pair in both
directions with different values, we take the MINIMUM of the two. That is not an
arbitrary choice: the meta-analysis of 3,766 experimentally characterised binders
(bioRxiv 2025.08.14.670059) found "aggregating asymmetric predictions between chains
using the minimum, rather than the average or maximum, yielded the most discriminative
results", and ipSAE_min was the best single predictor of in vitro binding they tested.

THE THREE NUMBERS WE RANK ON:

    ipsae_pep        ipSAE_min for binder-peptide on the G12D state. Higher is better.
                     This is the closest thing we have to "will this bind in a tube".
    ipsae_margin     ipsae_pep on G12D minus the same on wild type. Higher is better.
                     A confidence margin between two states - NOT a KD ratio.
    pep_centricity   ipsae_pep / (ipsae_pep + ipsae_mhc). Higher is better. A binder
                     that grips the MHC binds every peptide the MHC presents, which is
                     the opposite of what this challenge scores.

The meta-analysis put the maximum-F1 threshold for AF3 ipSAE_min at 0.61. That number
came from AF3 on a different target distribution, so it is used here only as a
reference line to report against, never as a pass/fail.

Run:  python pmhc/parse_boltz.py <boltz_results.json>
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "pmhc", "results")

BINDER, PEPTIDE, MHC = "C", "B", "A"
REFERENCE_THRESHOLD = 0.61      # AF3 max-F1 from the meta-analysis; reference only


def parse_rows(rows):
    """ipsae.py rows -> {(chain1, chain2): ipSAE}.

    The header is: Chn1 Chn2 PAE Dist Type ipSAE ipSAE_d0chn ipSAE_d0dom ... and the
    'max' summary lines are the per-chain-pair values we want. Rows arrive here as
    pre-split token lists from the notebook.
    """
    out = {}
    for r in rows:
        if len(r) < 6:
            continue
        c1, c2 = r[0], r[1]
        if c1 not in (BINDER, PEPTIDE, MHC) or c2 not in (BINDER, PEPTIDE, MHC):
            continue
        # find the first parseable float after the 'Type' column
        val = None
        for tok in r[5:]:
            try:
                val = float(tok)
                break
            except ValueError:
                continue
        if val is None:
            continue
        key = (c1, c2)
        # keep the best row seen for a pair (ipsae.py prints asym and max rows)
        if key not in out or val > out[key]:
            out[key] = val
    return out


def ipsae_min(pairs, x, y):
    """Minimum of the two directed ipSAE values for a chain pair, or None."""
    a, b = pairs.get((x, y)), pairs.get((y, x))
    vals = [v for v in (a, b) if v is not None]
    return min(vals) if vals else None


def main(argv):
    path = argv[1] if len(argv) > 1 else os.path.join(RESULTS, "boltz_results.json")
    if not os.path.exists(path):
        raise SystemExit(
            "no Boltz results at %s.\n"
            "Run kaggle/cellPMHCi_boltz.py on a Kaggle GPU session first, then save "
            "its boltz_results.json here." % path)
    data = json.load(open(path))
    print("predictions loaded: %d" % len(data))

    rows = []
    skipped = []
    for d in data:
        g = parse_rows(d.get("g12d_ipsae_rows", []))
        w = parse_rows(d.get("wt_ipsae_rows", []))
        pep_g = ipsae_min(g, BINDER, PEPTIDE)
        mhc_g = ipsae_min(g, BINDER, MHC)
        pep_w = ipsae_min(w, BINDER, PEPTIDE)
        if pep_g is None:
            skipped.append(d["name"])
            continue
        r = {"name": d["name"], "seq": d.get("seq", "")}
        r["ipsae_pep"] = round(pep_g, 4)
        r["ipsae_mhc"] = round(mhc_g, 4) if mhc_g is not None else None
        r["ipsae_pep_wt"] = round(pep_w, 4) if pep_w is not None else None
        r["ipsae_margin"] = (round(pep_g - pep_w, 4) if pep_w is not None else None)
        if mhc_g is not None and (pep_g + mhc_g) > 0:
            r["pep_centricity"] = round(pep_g / (pep_g + mhc_g), 4)
        else:
            r["pep_centricity"] = None
        for k in ("g12d_iptm", "g12d_ptm", "g12d_complex_plddt",
                  "g12d_confidence_score"):
            if k in d:
                r[k] = d[k]
        rows.append(r)

    if skipped:
        print("no binder-peptide ipSAE for %d design(s): %s"
              % (len(skipped), skipped[:8]))
    if not rows:
        raise SystemExit("nothing parsed - check the ipsae.py output format")

    vals = sorted(r["ipsae_pep"] for r in rows)
    print()
    print("--- ipSAE_min, binder against peptide, G12D state (n=%d) ---" % len(rows))
    print("    min %.3f  median %.3f  max %.3f"
          % (vals[0], vals[len(vals) // 2], vals[-1]))
    above = sum(1 for v in vals if v >= REFERENCE_THRESHOLD)
    print("    at or above %.2f (the AF3 max-F1 reference line): %d of %d"
          % (REFERENCE_THRESHOLD, above, len(rows)))
    print("    That line came from AF3 on other targets, so it is a reference point,")
    print("    not a pass mark, and Boltz-2 values are not interchangeable with AF3's.")

    marg = [r["ipsae_margin"] for r in rows if r["ipsae_margin"] is not None]
    if marg:
        m = sorted(marg)
        print()
        print("--- specificity margin: G12D minus wild type (n=%d) ---" % len(m))
        print("    min %+.3f  median %+.3f  max %+.3f"
              % (m[0], m[len(m) // 2], m[-1]))
        print("    positive (prefers the mutant): %d of %d"
              % (sum(1 for v in m if v > 0), len(m)))
        # Only meaningful with enough designs to have a spread at all; with one or
        # two it fires on arithmetic, not on evidence.
        if len(m) >= 5 and m[-1] - m[0] < 0.05:
            print()
            print("    WARNING: the margins are flat (spread %.3f over %d designs)."
                  % (m[-1] - m[0], len(m)))
            print("    Boltz is then not separating the two states, and this must")
            print("    not be used to rank them. Report it as a null result.")

    rows.sort(key=lambda r: (-(r["ipsae_margin"] or -9), -r["ipsae_pep"]))
    print()
    print("--- ranked by specificity margin, then binder-peptide ipSAE ---")
    print("%-22s %9s %9s %9s %8s %8s"
          % ("design", "ipSAE_pep", "ipSAE_wt", "margin", "pepcent", "ipTM"))
    for r in rows:
        print("%-22s %9.3f %9s %9s %8s %8s"
              % (r["name"], r["ipsae_pep"],
                 "%.3f" % r["ipsae_pep_wt"] if r["ipsae_pep_wt"] is not None else "-",
                 "%+.3f" % r["ipsae_margin"] if r["ipsae_margin"] is not None else "-",
                 "%.2f" % r["pep_centricity"] if r["pep_centricity"] is not None
                 else "-",
                 "%.2f" % r["g12d_iptm"] if "g12d_iptm" in r else "-"))

    out = os.path.join(RESULTS, "boltz_parsed.json")
    json.dump(rows, open(out, "w"), indent=1)
    print()
    print("wrote %s" % out)
    print()
    print("NOT A KD. ipSAE is a confidence score over a predicted interface. The")
    print("margin between two states is a difference of confidences, not a ratio of")
    print("dissociation constants, and must not be presented as one.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
