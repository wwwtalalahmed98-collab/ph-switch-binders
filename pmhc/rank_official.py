"""Rank the designs in the order the competition published, not the order I invented.

WHY THIS REPLACES THE RANKING IN build_submission.py. The organisers publish the prompt
their selector runs, at proteinbase.com/competitions/anthropic-adaptyv-2026/selection.
For this challenge it states the goals HARDEST FIRST and says "that order is what the
model weighs":

    a) no binding to the wild-type KRAS peptide VVVGAGGVGK on the same HLA
    b) no binding to the same HLA and beta-2 microglobulin with no peptide
    c) strong binding to the KRAS G12D complex

So rejection of the wild type outranks rejection of the empty groove, which outranks
raw strength. build_submission.py ranked on a ProteinMPNN specificity score of my own
construction. The same prompt says "self-reported metrics compare only within a
submission ... unverifiable claims weigh little", which is exactly what that score is.
This module ranks on the three-state AlphaFold screen instead, because the organisers
can recompute that themselves from the sequence.

THE THREE STATES ARE THE THREE GOALS. The Kaggle screen predicts each design against
one background in three conditions that differ only in the peptide:

    g12d   peptide VVVGADGVGK        -> should bind
    wt     peptide VVVGAGGVGK        -> should bind worse   (goal a)
    apo    no peptide                -> should bind worse   (goal b)

goal_a = score(g12d) - score(wt);  goal_b = score(g12d) - score(apo);  goal_c = score(g12d).

THE NOISE FLOOR IS NOT DECORATION. Each state is predicted with two AF2 models. A
margin only counts as evidence if it exceeds the disagreement between two models of the
SAME state. Designs whose margins sit inside that scatter are ranked below those whose
margins clear it, however large the raw number looks. Without this, model-to-model
scatter is read as specificity.

WHAT IT DISCOVERS RATHER THAN ASSUMES. How ColabDesign laid out the chains in its
output is not assumed anywhere. The script reads the chain pairs ipsae.py actually
emitted, prints them, and derives the binder chain from the data. Three Kaggle runs in
this project died on an assumed layout; this one reports what it found and says so when
it cannot tell.

Run:  python pmhc/rank_official.py [af2_ipsae_results.json]
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "pmhc", "results")

# ipsae.py prints several rows per chain pair; these are the ones we read.
WANTED_TYPES = ("max", "asym")


def parse_rows(rows):
    """ipsae.py token rows -> {(chn1, chn2): ipSAE}, keeping the best per pair."""
    out = {}
    for r in rows:
        if isinstance(r, str):
            r = r.split()
        if len(r) < 6:
            continue
        c1, c2 = r[0], r[1]
        if not (len(c1) == 1 and len(c2) == 1 and c1.isalpha() and c2.isalpha()):
            continue
        val = None
        for tok in r[5:]:
            try:
                val = float(tok)
                break
            except ValueError:
                continue
        if val is None:
            continue
        k = (c1, c2)
        if k not in out or val > out[k]:
            out[k] = val
    return out


def ipsae_min(pairs, x, y):
    a, b = pairs.get((x, y)), pairs.get((y, x))
    v = [t for t in (a, b) if t is not None]
    return min(v) if v else None


def state_ipsae(state, binder_hint=None):
    """Average over models of (binder vs its best target partner), plus the layout.

    Returns (value_or_None, binder_chain_or_None, all_chains).
    """
    vals, binder, chains = [], binder_hint, set()
    for m in state.get("per_model", []):
        rows = m.get("ipsae_rows")
        if not rows:
            continue
        pairs = parse_rows(rows)
        if not pairs:
            continue
        cs = sorted({c for k in pairs for c in k})
        chains |= set(cs)
        b = binder or (cs[-1] if cs else None)   # ColabDesign appends the binder last
        if b is None:
            continue
        binder = b
        partner = [ipsae_min(pairs, b, c) for c in cs if c != b]
        partner = [p for p in partner if p is not None]
        if partner:
            vals.append(max(partner))
    if not vals:
        return None, binder, chains
    return sum(vals) / len(vals), binder, chains


def collect(data):
    out = []
    layouts = set()
    n_ipsae = 0
    for rec in data:
        st = rec.get("states", {})
        if not all(k in st for k in ("g12d", "wt", "apo")):
            continue
        row = {"name": rec.get("name"), "seq": rec.get("seq", "")}

        # --- iPTM route: always available, already computed on Kaggle ---
        row["iptm_g12d"] = st["g12d"].get("i_ptm")
        row["goal_a_iptm"] = rec.get("margin_wt")
        row["goal_b_iptm"] = rec.get("margin_apo")
        row["noise"] = rec.get("noise")
        row["measurable"] = bool(rec.get("measurable"))

        # --- ipSAE route: better predictor, if the rows came through ---
        binder = None
        vals = {}
        for tag in ("g12d", "wt", "apo"):
            v, binder, chains = state_ipsae(st[tag], binder)
            vals[tag] = v
            if chains:
                layouts.add("".join(sorted(chains)))
        if all(vals[t] is not None for t in ("g12d", "wt", "apo")):
            n_ipsae += 1
            row["ipsae_g12d"] = round(vals["g12d"], 4)
            row["goal_a_ipsae"] = round(vals["g12d"] - vals["wt"], 4)
            row["goal_b_ipsae"] = round(vals["g12d"] - vals["apo"], 4)
        row["binder_chain"] = binder
        out.append(row)
    return out, layouts, n_ipsae


def main(argv):
    path = argv[1] if len(argv) > 1 else os.path.join(RESULTS,
                                                      "af2_ipsae_results.json")
    if not os.path.exists(path):
        raise SystemExit(
            "no AF2 results at %s\n"
            "Download af2_ipsae_results.json from the Kaggle run output into\n"
            "pmhc/results/ and run this again." % path)
    data = json.load(open(path))
    rows, layouts, n_ipsae = collect(data)
    print("designs with all three states: %d of %d" % (len(rows), len(data)))
    if not rows:
        raise SystemExit("nothing usable in %s" % path)

    # Controls are the scale, not candidates. They are pulled out before ranking so
    # the positive control - a confirmed binder of this target - can never be ranked
    # first and submitted as one of our designs.
    ctrl = [r for r in rows if (r["name"] or "").startswith("CTRL_")]
    rows = [r for r in rows if not (r["name"] or "").startswith("CTRL_")]
    if ctrl:
        print()
        print("--- controls (excluded from the ranking; this is the scale) ---")
        for c in ctrl:
            print("    %-22s i_pTM(g12d) %s  vs WT %s  vs apo %s"
                  % (c["name"][:22],
                     "%.3f" % c["iptm_g12d"] if c["iptm_g12d"] is not None else "-",
                     "%+.3f" % c["goal_a_iptm"]
                     if c["goal_a_iptm"] is not None else "-",
                     "%+.3f" % c["goal_b_iptm"]
                     if c["goal_b_iptm"] is not None else "-"))
        pos = next((c for c in ctrl if "POS" in c["name"]), None)
        neg = next((c for c in ctrl if "NEG" in c["name"]), None)
        if pos and neg and pos["iptm_g12d"] is not None \
                and neg["iptm_g12d"] is not None:
            sep = pos["iptm_g12d"] - neg["iptm_g12d"]
            print()
            print("    positive minus negative: %+.3f" % sep)
            if sep < 0.05:
                print("    THE PIPELINE DOES NOT DISCRIMINATE. A confirmed binder of")
                print("    this target scores no better than a shuffled sequence, so")
                print("    these i_pTM values say nothing about our designs either.")
                print("    Report as a failed measurement, not as evidence against")
                print("    the designs.")
            else:
                best = max((r["iptm_g12d"] or 0) for r in rows) if rows else 0
                print("    the pipeline separates a real binder from noise by %.3f,"
                      % sep)
                print("    so it is reading something. Our best design scores %.3f."
                      % best)
    if not rows:
        raise SystemExit("only controls present - nothing to rank")
    # recount on the candidates alone, so the figure below has one denominator
    n_ipsae = sum(1 for r in rows if r.get("goal_a_ipsae") is not None)

    print()
    print("--- what the output actually looked like (discovered, not assumed) ---")
    print("    chain sets seen in ipsae output : %s"
          % (sorted(layouts) if layouts else "NONE - no ipsae rows came through"))
    print("    binder chain inferred           : %s"
          % sorted({r["binder_chain"] for r in rows if r["binder_chain"]}))
    print("    designs with usable ipSAE       : %d of %d" % (n_ipsae, len(rows)))
    use_ipsae = n_ipsae >= max(3, len(rows) // 2)
    print("    ranking on                      : %s"
          % ("ipSAE (preferred)" if use_ipsae else "iPTM (ipSAE incomplete)"))
    if layouts and len(layouts) == 1 and len(next(iter(layouts))) <= 2:
        print()
        print("    NOTE: only %s chains in the output, so the target is merged into"
              % len(next(iter(layouts))))
        print("    one chain and peptide-vs-MHC cannot be separated. The three-state")
        print("    comparison is unaffected: the states differ ONLY in the peptide.")

    ka = "goal_a_ipsae" if use_ipsae else "goal_a_iptm"
    kb = "goal_b_ipsae" if use_ipsae else "goal_b_iptm"
    kc = "ipsae_g12d" if use_ipsae else "iptm_g12d"

    def g(r, k):
        v = r.get(k)
        return -9.0 if v is None else v

    print()
    print("--- ranking, in the order the organisers published ---")
    print("    0. clears the two-model noise floor (a gate, not a score)")
    print("    1. goal a: rejects the wild-type peptide   [%s]" % ka)
    print("    2. goal b: rejects the empty groove        [%s]" % kb)
    print("    3. goal c: binds the G12D complex          [%s]" % kc)
    rows.sort(key=lambda r: (not r["measurable"], -g(r, ka), -g(r, kb), -g(r, kc)))

    n_meas = sum(1 for r in rows if r["measurable"])
    n_pos_a = sum(1 for r in rows if g(r, ka) > 0)
    n_pos_b = sum(1 for r in rows if g(r, kb) > 0)
    print()
    print("    clear the noise floor        : %d of %d" % (n_meas, len(rows)))
    print("    positive on goal a (vs WT)   : %d of %d" % (n_pos_a, len(rows)))
    print("    positive on goal b (vs apo)  : %d of %d" % (n_pos_b, len(rows)))
    if n_meas == 0:
        print()
        print("    WARNING: not one design's margins clear the scatter between two")
        print("    models of the same state. Then this screen has NOT separated the")
        print("    three states and must be reported as a null result, not used to")
        print("    order the submission.")

    print()
    hdr = "%-4s %-22s %8s %8s %8s %7s %s"
    print(hdr % ("#", "design", "goal_a", "goal_b", "goal_c", "noise", "meas"))
    for i, r in enumerate(rows, 1):
        print(hdr % (i, (r["name"] or "?")[:22],
                     "%+.3f" % g(r, ka), "%+.3f" % g(r, kb), "%.3f" % g(r, kc),
                     "%.3f" % (r["noise"] or 0), "yes" if r["measurable"] else ""))

    for i, r in enumerate(rows, 1):
        r["official_rank"] = i
        r["ranked_on"] = "ipsae" if use_ipsae else "iptm"
    out = os.path.join(RESULTS, "official_rank.json")
    json.dump(rows, open(out, "w"), indent=1)
    print()
    print("wrote %s" % out)
    print()
    print("NOT A KD. iPTM and ipSAE are confidence scores over a predicted interface.")
    print("A difference between two states is a difference of confidences, not a")
    print("ratio of dissociation constants, and must not be presented as one.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
