# Select the 20 GLP1R designs for Track 3.
#
# THE RULE IS WRITTEN HERE BEFORE THE DESIGNS EXIST. That is the whole point. Yesterday's
# control run was trustworthy because its verdict was defined before the numbers arrived;
# the TNF switch rule was not, and we selected on a proxy that turned out to correlate
# +0.24 with what we wanted. A rule written after seeing 42 designs is a rule tuned to
# pick the ones we already liked.
#
# TWO POPULATIONS, SELECTED SEPARATELY. They answer different questions and must not be
# ranked against each other.
#
#   LANDING-PAD BINDERS - 14 slots. Bind the extracellular domain. This is the job our
#   pipeline does well and the lab can measure reliably. By mechanism these are
#   ANTAGONISTS: exendin(9-39) is this exact molecule class, sold as a drug. We are not
#   pretending otherwise.
#     pass: i_ptm >= 0.50, on_site >= 8 of 22, hydropathy <= 0.0
#     rank: i_ptm, then on_site
#
#   CHANNEL THREADERS - 6 slots. Reach into the helical bundle, which is what actually
#   activates the receptor. Declared long shots. Their justification is the challenge's
#   own objective 2, which is explicit: "The binder has to be an agonist... A design that
#   binds the receptor but does not activate it does not meet this objective", and the
#   winner must meet both objectives. It is NOT, as previously claimed here, that the
#   selection prompt rewards long shots - the published prompt says no such thing.
#     pass: threaded (N-terminus is the deep end), i_ptm >= 0.40, tm_contacts >= 5,
#           hydropathy <= 0.5
#     rank: tm_contacts, then i_ptm
#
# WHY THOSE THRESHOLDS.
#   i_ptm is from an INDEPENDENT re-prediction at 3 recycles, never the optimiser's own
#     number. 0.50 is where our submitted TNF set sat (0.576-0.711 measured). Threaders
#     get 0.40 because threading a channel is a harder prediction problem and holding
#     them to the easier job's bar would simply empty that population.
#   on_site counts intended contacts actually made. Below 8 of 22 the design has drifted
#     off the groove and is binding somewhere we did not choose.
#   hydropathy is mean Kyte-Doolittle. CORRECTION, 24 Oct: this was previously justified
#     here as "one of the three screens the competition records". That is not true. The
#     published selection prompt was checked and contains no mention of pLDDT, ipSAE or
#     hydropathy anywhere; it weighs method novelty, design diversity and an unspecified
#     "small set of confidence metrics". The filter stays because it is defensible on its
#     own terms - a greasy de novo binder is an aggregation and solubility risk, and the
#     assay immobilises the design and flows detergent-solubilised receptor over it - but
#     it is OUR criterion, not the organisers'. <= 0.0 means not greasy. Threaders are
#     allowed to 0.5 because a chain that threads a partly buried channel cannot be as
#     polar as a surface binder; a stated concession, not an oversight.
#   threaded is the ONLY mechanism screen we have, and it is a weak one. Depth magnitude
#     was tested on 23 Oct and failed: the real hormone reached 34.1 A, the blocker-shaped
#     truncation 33.3 A, and a scrambled nonsense sequence 37.0 A. What did separate was
#     polarity - both real peptides put their N-terminus at the deep end to within 0.1 A
#     while the scrambled control sprawled in backwards. So this filter says "threaded at
#     all", NOT "activates". ACTIVATION IS UNSCREENED. The methods file must say so.
#
# IF FEWER THAN THE SLOTS PASS, THE SET IS SHORT. We do not pad with designs that failed
# our own filters just to reach 20. A short honest set beats a full dishonest one.
import csv, difflib, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.abspath(os.path.join(HERE, "..", "submission"))
N_ECD, N_THREAD = 14, 6
MAX_IDENT = 0.80          # near-duplicates waste a slot; keep the better one

ECD = dict(i_ptm=0.50, on_site=8, kd=0.0)
THR = dict(i_ptm=0.40, tm=5, kd=0.5, on_site=5, max_deep=40.0)

# AMENDMENT, 24 Oct, made after seeing 2 of ~12 threaders and recorded as an amendment
# rather than quietly folded into the rule above.
#
# glp1r_thread_001 came back with i_ptm 0.800 - the highest number this project has ever
# produced - reaching 47.4 A with its N-terminus at 46.6 A. The native hormone stops at
# 34.2 A. Measuring our own target explained it:
#
#     depth band      atoms in our crop      atoms in the intact receptor
#     33.2-40.0 A           259                        284   (we kept 91%)
#     40.0-45.0 A           107                        195   (we kept 55%)
#     45.0-50.0 A             4                        211   (we kept  2%)
#
# Our crop ends at 46.8 A. The design is threading into a hole that exists only because
# we cut the bottom off the receptor; in the intact protein that space is solid helix and
# the design would clash into it. It is an artifact of target preparation, not a find.
#
# Hence a depth ceiling of 40 A, where the crop still represents at least half of the
# real receptor, and a requirement to contact at least 5 residues the hormone itself uses
# - thread_001 contacted 2. This amendment REMOVES the best-scoring design we have, which
# is the opposite direction from tuning a rule to flatter the results, but it is still an
# amendment made with data in hand and should be read as one.


def ident(a, b):
    return difflib.SequenceMatcher(None, a, b).ratio()


def dedupe(cands):
    """Keep the highest-ranked of any near-identical group. Two designs that differ by a
    couple of residues measure as one result in the lab and waste a slot."""
    kept, dropped = [], []
    for c in cands:
        twin = next((k for k in kept if ident(k["sequence"], c["sequence"]) >= MAX_IDENT),
                    None)
        if twin:
            dropped.append((c, twin["name"]))
        else:
            kept.append(c)
    return kept, dropped


def screen(recs, mode):
    passed, failed = [], []
    for r in recs:
        why = []
        kd = r.get("hydropathy")
        ip = r.get("true_iptm")
        if mode == "ecd":
            if ip is None or ip < ECD["i_ptm"]:
                why.append("i_ptm %.3f < %.2f" % (ip or 0, ECD["i_ptm"]))
            if len(r.get("on_site", [])) < ECD["on_site"]:
                why.append("on_site %d < %d" % (len(r.get("on_site", [])), ECD["on_site"]))
            if kd is None or kd > ECD["kd"]:
                why.append("hydropathy %+.2f > %+.2f" % (kd or 0, ECD["kd"]))
        else:
            if not r.get("threaded"):
                why.append("not threaded (N-terminus is not the deep end)")
            if ip is None or ip < THR["i_ptm"]:
                why.append("i_ptm %.3f < %.2f" % (ip or 0, THR["i_ptm"]))
            if r.get("tm_contacts", 0) < THR["tm"]:
                why.append("tm_contacts %d < %d" % (r.get("tm_contacts", 0), THR["tm"]))
            if kd is None or kd > THR["kd"]:
                why.append("hydropathy %+.2f > %+.2f" % (kd or 0, THR["kd"]))
            if len(r.get("on_site", [])) < THR["on_site"]:
                why.append("touches only %d residues the hormone itself uses, needs %d"
                           % (len(r.get("on_site", [])), THR["on_site"]))
            if (r.get("deepest") or 0) > THR["max_deep"]:
                why.append("deepest %.1f A > %.1f - past where our crop still represents "
                           "the real receptor" % (r.get("deepest") or 0, THR["max_deep"]))
        (failed if why else passed).append((r, why))
    return [p[0] for p in passed], failed


def main(path):
    recs = json.load(open(path))
    for r in recs:
        r.setdefault("on_site", [])
    ecd = [r for r in recs if r.get("mode") == "ecd"]
    thr = [r for r in recs if r.get("mode") == "thread"]
    print("read %d designs: %d landing-pad, %d threaders\n" % (len(recs), len(ecd), len(thr)))

    chosen, report = [], []
    for mode, pool, want, rank in (
            ("ecd", ecd, N_ECD,
             lambda r: (r["true_iptm"], len(r["on_site"]))),
            ("thread", thr, N_THREAD,
             lambda r: (r.get("tm_contacts", 0), r["true_iptm"]))):
        ok, bad = screen(pool, mode)
        ok.sort(key=rank, reverse=True)
        kept, dropped = dedupe(ok)
        take = kept[:want]
        label = "landing-pad binders" if mode == "ecd" else "channel threaders"
        print("=== %s: %d of %d passed the filters, %d near-duplicates removed, "
              "taking %d of %d slots ===" % (label, len(ok), len(pool), len(dropped),
                                             len(take), want))
        if len(take) < want:
            print("  SHORT by %d. Not padding with designs that failed our own filters."
                  % (want - len(take)))
        for r in take:
            extra = ("%2d TM, deep %5.1f, nterm %5.1f"
                     % (r.get("tm_contacts", 0), r.get("deepest", 0),
                        r.get("nterm_depth", 0))) if mode == "thread" else \
                    ("%2d of 22 on site" % len(r["on_site"]))
            print("  %-18s len %3d  i_ptm %.3f  kd %+5.2f  hlx %.2f  %s"
                  % (r["name"], r["length"], r["true_iptm"], r["hydropathy"],
                     r.get("helicity", 0), extra))
        for r, tw in dropped:
            print("  dup  %-18s ~ %s" % (r["name"], tw))
        for r, why in sorted(bad, key=lambda x: -(x[0].get("true_iptm") or 0))[:8]:
            print("  out  %-18s %s" % (r["name"], "; ".join(why)))
        print()
        chosen.extend((mode, r) for r in take)
        report.append((mode, len(pool), len(ok), len(dropped), len(take), want))

    if not chosen:
        raise SystemExit("nothing passed. Do not submit; read the failures above.")

    os.makedirs(OUT, exist_ok=True)
    prop = os.path.join(OUT, "challenge3_proposed.csv")
    mets = os.path.join(OUT, "challenge3_metrics.csv")
    with open(prop, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "sequence", "molecule_class"])
        for i, (mode, r) in enumerate(chosen, 1):
            tag = "glp1r_ecd_%02d" if mode == "ecd" else "glp1r_thr_%02d"
            n = sum(1 for m, _ in chosen[:i] if m == mode)
            r["sub_name"] = tag % n
            # CORRECTION, 24 Oct: threaders were being written as molecule_class
            # "peptide". The challenge page allows exactly four formats - single chain
            # protein, nanobody, scFv, Fab - and "peptide" is not one of them, so that
            # would have been an invalid submission for six of the twenty slots. A
            # 30-45 residue chain is a single chain protein; the stated length range is
            # 10-250 amino acids, so it qualifies on length too.
            w.writerow([r["sub_name"], r["sequence"], "single_chain"])
    with open(mets, "w", newline="", encoding="utf-8") as fh:
        cols = ["sub_name", "name", "mode", "length", "true_iptm", "design_iptm", "plddt",
                "hydropathy", "helicity", "n_on_site", "tm_contacts", "ecd_contacts",
                "deepest", "nterm_depth", "cterm_depth", "threaded"]
        w = csv.writer(fh)
        w.writerow(cols)
        for mode, r in chosen:
            w.writerow([r.get(c, "") if c != "n_on_site" else len(r["on_site"])
                        for c in cols])

    print("wrote %s and %s  (%d designs)" % (prop, mets, len(chosen)))
    print("\nBoth files are embargoed: they carry sequences and must stay out of the "
          "public repo until the window closes. Verify with git check-ignore.")
    print("\nSTILL UNSCREENED: whether any threader activates the receptor. Nothing in "
          "this script tests that, and nothing we have can.")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else
         os.path.abspath(os.path.join(HERE, "..", "data", "glp1r", "designs.json")))
