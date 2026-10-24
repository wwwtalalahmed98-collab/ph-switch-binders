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
#   activates the receptor. Declared long shots.
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
#   hydropathy is mean Kyte-Doolittle, one of the three screens the competition records.
#     <= 0.0 means not greasy. Threaders are allowed to 0.5 because a peptide that
#     threads a partly buried channel cannot be as polar as a surface binder; this is a
#     stated concession, not an oversight.
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
THR = dict(i_ptm=0.40, tm=5, kd=0.5)


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
            # a 30-45 residue threader is a peptide, not a folded single chain; the
            # competition records nine molecule classes and the wrong one misdescribes
            # what the lab is being asked to make
            w.writerow([r["sub_name"], r["sequence"],
                        "single_chain" if mode == "ecd" else "peptide"])
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
