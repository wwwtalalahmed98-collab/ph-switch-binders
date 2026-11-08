"""Run the organisers' own published expression checks, which we had never run.

WHERE THESE COME FROM. Adaptyv publish a protein-QC guide with explicit thresholds
for designed binders. This project has been checking charge, hydrophobic fraction and
mean Kyte-Doolittle since challenge 1, but four of their sequence-level criteria were
never checked at all: the instability index, the odd-cysteine rule, the deamidation
motifs, and the runs (K/R stretches and hydrophobic stretches). They are arithmetic on
the sequence, they cost nothing, and they are the organisers' own.

    instability index  < 40 standard, < 30 stringent   (Biopython)
    GRAVY              < 0.4 standard, < 0.2 stringent (Biopython; mean Kyte-Doolittle)
    odd cysteine count                    -> redesign
    NG / NS / NT                          -> flag, deamidation risk
    three or more consecutive K or R      -> flag, proteolysis risk
    hydrophobic run of six or more        -> redesign, aggregation risk

WHAT THIS IS NOT. None of it predicts binding, and Adaptyv say so themselves: these
work "as filters (eliminating failures)" rather than as predictors of success. A design
that passes every line here is not a better binder, only a likelier one to survive
being made.

Run:  python pmhc/qc_adaptyv.py [file.csv ...]
      default: the final twenty, then both pools it was drawn from
"""
import csv
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SUB = os.path.join(ROOT, "submission")
RESULTS = os.path.join(ROOT, "pmhc", "results")

HYDROPHOBIC = set("AVILMFWY")
DEAMIDATION = ("NG", "NS", "NT")

INSTABILITY_STD, INSTABILITY_STRICT = 40.0, 30.0
GRAVY_STD, GRAVY_STRICT = 0.4, 0.2


def longest_run(seq, members):
    best = run = 0
    for ch in seq:
        run = run + 1 if ch in members else 0
        best = max(best, run)
    return best


def liabilities(seq):
    out = {}
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    pa = ProteinAnalysis(seq)
    out["instability"] = round(pa.instability_index(), 2)
    out["gravy"] = round(pa.gravy(), 3)
    out["pi"] = round(pa.isoelectric_point(), 2)
    out["n_cys"] = seq.count("C")
    out["odd_cys"] = (out["n_cys"] % 2) == 1
    out["deamidation"] = sum(seq.count(m) for m in DEAMIDATION)
    out["kr_run"] = longest_run(seq, set("KR"))
    out["hydrophobic_run"] = longest_run(seq, HYDROPHOBIC)
    redesign, flag = [], []
    if out["odd_cys"]:
        redesign.append("odd cysteine (%d)" % out["n_cys"])
    if out["hydrophobic_run"] >= 6:
        redesign.append("hydrophobic run %d" % out["hydrophobic_run"])
    if out["instability"] >= INSTABILITY_STD:
        redesign.append("instability %.1f" % out["instability"])
    if out["gravy"] >= GRAVY_STD:
        redesign.append("GRAVY %+.2f" % out["gravy"])
    if out["deamidation"]:
        flag.append("%d deamidation motif(s)" % out["deamidation"])
    if out["kr_run"] >= 3:
        flag.append("K/R run %d" % out["kr_run"])
    if out["instability"] >= INSTABILITY_STRICT > out["instability"] - 1e-9 or \
            INSTABILITY_STRICT <= out["instability"] < INSTABILITY_STD:
        flag.append("instability above the stringent 30")
    if GRAVY_STRICT <= out["gravy"] < GRAVY_STD:
        flag.append("GRAVY above the stringent 0.2")
    out["redesign"] = redesign
    out["flag"] = flag
    out["verdict"] = "redesign" if redesign else ("flag" if flag else "clean")
    return out


def load(path):
    rows = list(csv.DictReader(open(path)))
    if not rows:
        return []
    cols = list(rows[0].keys())
    scol = next((c for c in cols if c.lower() in ("sequence", "seq")), None)
    ncol = next((c for c in cols if c.lower() in ("name", "id")), None)
    if scol is None:
        return []
    out = []
    for i, r in enumerate(rows, 1):
        s = (r[scol] or "").strip().upper()
        if s:
            out.append(((r[ncol].strip() if ncol and r.get(ncol)
                         else "row_%03d" % i), s))
    return out


def run(path):
    rows = load(path)
    print()
    print("=" * 78)
    print("%s   (%d sequences)" % (os.path.basename(path), len(rows)))
    print("=" * 78)
    if not rows:
        print("  nothing readable here")
        return {}
    hdr = "%-16s %4s %6s %7s %5s %4s %4s %4s %-9s"
    print(hdr % ("name", "len", "instab", "GRAVY", "pI", "Cys", "deam", "hrun",
                 "verdict"))
    tally, worst = {}, []
    for name, seq in rows:
        q = liabilities(seq)
        tally[q["verdict"]] = tally.get(q["verdict"], 0) + 1
        print(hdr % (name[:16], len(seq), "%.1f" % q["instability"],
                     "%+.3f" % q["gravy"], "%.1f" % q["pi"], q["n_cys"],
                     q["deamidation"], q["hydrophobic_run"], q["verdict"]))
        if q["redesign"]:
            worst.append((name, "; ".join(q["redesign"])))
    n = len(rows)
    print()
    for k in ("clean", "flag", "redesign"):
        v = tally.get(k, 0)
        print("    %-9s %3d of %d  (%4.1f%%)" % (k, v, n, 100.0 * v / n))
    if worst:
        print()
        print("    failing a redesign-level rule:")
        for name, why in worst:
            print("      %-20s %s" % (name[:20], why))
    return tally


def main(argv):
    paths = argv[1:] or [
        os.path.join(SUB, "challenge5_final20.csv"),
        os.path.join(SUB, "challenge5_candidates.csv"),
        os.path.join(RESULTS, "boltzgen_selected.csv"),
    ]
    first = None
    for p in paths:
        if not os.path.exists(p):
            print("skipping missing %s" % p)
            continue
        t = run(p)
        if first is None:
            first = t
    print()
    print("These are expression and manufacturability filters from the organisers'")
    print("own QC guide. They eliminate likely failures. They do not predict")
    print("binding, and Adaptyv say so explicitly.")
    if first and first.get("redesign"):
        print()
        print("THE SUBMITTED SET CONTAINS DESIGNS FAILING A REDESIGN-LEVEL RULE.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
