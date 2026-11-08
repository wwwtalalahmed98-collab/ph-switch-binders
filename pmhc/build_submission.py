"""Assemble the final 20 designs for Challenge 5, ranked, with the methods text.

RANKING. The competition ranks Track 3 entries by, in order:
  1. KD(wild-type peptide) / KD(G12D peptide)   - discrimination against wild type
  2. KD(peptide-free HLA)  / KD(G12D peptide)   - discrimination against empty groove
  3. affinity, as a tie-break only

So the ordering here puts discrimination first and everything else behind it. Our proxy
for criterion 1 is the ProteinMPNN specificity margin from pmhc/specificity_screen.py:
how much more the pocket prefers Asp (G12D) over Gly (wild type) at peptide position 6,
measured against the bare-groove baseline. It is a sequence model's opinion, not a KD,
and the methods text says so.

DIVERSITY. At most MAX_PER_PARENT designs may come from the same parent backbone.
Twenty variations on one pose is one bet, not twenty, and the organisers select for
diversity as well as confidence.

THE DE NOVO DECLARATION. Track 3 says "de novo only; modifying an existing binder is
not allowed". These backbones come from the scaffold set published with Liu et al.
(Science 389:386-391, 2025), which are themselves de novo designs built for DIFFERENT
pMHC targets. Every residue of the binder was rewritten by ProteinMPNN against our
target. Each design's identity to its parent scaffold is computed, printed, and stated
in the methods text, and anything at or above MAX_PARENT_IDENTITY is dropped. We declare
the provenance rather than hide it; the organisers can then judge it themselves.

Run:  python pmhc/build_submission.py
"""
import csv
import json
import os
import sys
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "pmhc", "results")
SUB = os.path.join(ROOT, "submission")

N_FINAL = 20
MAX_PER_PARENT = 3
MAX_PARENT_IDENTITY = 0.45     # above this the sequence is too close to its parent
MIN_LEN, MAX_LEN = 10, 250     # Track 3 limits

PEP = "VVVGADGVGK"


def load():
    spec = os.path.join(RESULTS, "specificity.json")
    if os.path.exists(spec):
        d = json.load(open(spec))
        return d["designs"], d.get("baseline_margin"), True
    cand = os.path.join(RESULTS, "mpnn_candidates.json")
    if not os.path.exists(cand):
        raise SystemExit("run pmhc/run_mpnn.py (and ideally "
                         "pmhc/specificity_screen.py) first")
    print("WARNING: specificity.json not found. Ranking will fall back to interface")
    print("geometry and buildability only, which is NOT the criterion we are scored")
    print("on. Run pmhc/specificity_screen.py before submitting if at all possible.")
    print()
    return json.load(open(cand)), None, False


def main():
    rows, baseline, have_spec = load()
    print("candidates available: %d" % len(rows))

    # --- hard rules first ---
    # Controls can never be submitted. The positive control is a confirmed binder of
    # this exact target taken from the PDB, so submitting it would be both a rule
    # breach and a lie. Refuse loudly rather than filter quietly.
    ctrl = [r for r in rows
            if str(r.get("design", r.get("name", ""))).startswith("CTRL_")]
    if ctrl:
        raise SystemExit(
            "refusing to build a submission from a set containing controls: %s\n"
            "These are measurement scales, not designs. The positive control is a "
            "known binder of this target." % [r.get("design") for r in ctrl][:4])
    n0 = len(rows)
    rows = [r for r in rows if MIN_LEN <= r["length"] <= MAX_LEN]
    n_len = len(rows)
    rows = [r for r in rows if "C" not in r["seq"]]
    n_cys = len(rows)
    rows = [r for r in rows if r["parent_identity"] < MAX_PARENT_IDENTITY]
    n_id = len(rows)
    print()
    print("--- hard rules ---")
    print("    start                                  : %d" % n0)
    print("    length within %d-%d aa                  : %d" % (MIN_LEN, MAX_LEN, n_len))
    print("    no cysteine                            : %d" % n_cys)
    print("    parent identity below %.0f%%             : %d"
          % (100 * MAX_PARENT_IDENTITY, n_id))
    if not rows:
        raise SystemExit("nothing survived the hard rules")

    # --- ranking ---
    print()
    print("--- ranking, stated before the numbers are read ---")
    if have_spec:
        print("    1. shuffle-controlled specificity score (proxy for KD(WT)/KD(G12D))")
        print("    2. number of heavy-atom contacts to the Asp6 mutation")
        print("    3. peptide share of the interface (peptide-centric = specific)")
        print("    4. ProteinMPNN score, as a tie-break only")
        print()
        print("    The binder-fit measure is NOT used for ranking. On the real data")
        print("    its whole range was 0.16 log units against 1.15 for the measure")
        print("    above, and it correlated -0.03 with it. A measure with almost no")
        print("    dynamic range that agrees with nothing is not evidence; it is")
        print("    reported in the detail table and left out of the ordering.")
        rows.sort(key=lambda r: (-r.get("delta_vs_shuffled", 0.0),
                                 -r["n_mut_contacts"], -r["pep_fraction"],
                                 r["mpnn_score"]))
    else:
        print("    NO SPECIFICITY DATA. Falling back to:")
        print("    1. contacts to the Asp6 mutation")
        print("    2. peptide share of the interface")
        print("    3. ProteinMPNN score")
        rows.sort(key=lambda r: (-r["n_mut_contacts"], -r["pep_fraction"],
                                 r["mpnn_score"]))
    print("    then at most %d per parent backbone, taking the best of each"
          % MAX_PER_PARENT)

    seen = Counter()
    final = []
    for r in rows:
        if seen[r["parent"]] >= MAX_PER_PARENT:
            continue
        seen[r["parent"]] += 1
        final.append(r)
        if len(final) == N_FINAL:
            break

    print()
    print("    distinct parent backbones in the final set: %d"
          % len({r["parent"] for r in final}))
    print("    designs selected                          : %d" % len(final))
    if len(final) < N_FINAL:
        print("    NOTE: fewer than %d available; submitting %d."
              % (N_FINAL, len(final)))

    os.makedirs(SUB, exist_ok=True)

    # --- the submission file itself ---
    path = os.path.join(SUB, "challenge5_kras_g12d.csv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["name", "sequence", "molecule_class"])
        for i, r in enumerate(final, 1):
            r["rank"] = i
            r["submit_name"] = "krasG12D_%02d" % i
            w.writerow([r["submit_name"], r["seq"], "protein"])
    print()
    print("wrote %s" % path)

    # --- the full table, for us and for the methods text ---
    detail = os.path.join(SUB, "challenge5_kras_g12d_detail.tsv")
    cols = ["rank", "submit_name", "parent", "length", "tier", "parent_identity",
            "n_mut_contacts", "n_pep_contacts", "n_mhc_contacts", "pep_fraction",
            "mpnn_score", "logp_asp6", "logp_gly6", "margin", "delta_vs_baseline", "margin_shuffled", "delta_vs_shuffled",
            "binder_score_g12d", "binder_score_wt", "binder_fit_delta",
            "charged_frac", "hydrophobic_frac", "pro_frac", "gly_frac", "kd", "seq"]
    with open(detail, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for r in final:
            w.writerow([r.get(c, "") for c in cols])
    print("wrote %s" % detail)

    print()
    print("--- the 20 ---")
    hdr = "%-14s %-10s %4s %-8s %5s %4s %8s %8s"
    print(hdr % ("name", "parent", "len", "tier", "ident", "mut",
                 "vs shuf", "vs bare"))
    for r in final:
        print(hdr % (r["submit_name"], r["parent"], r["length"], r["tier"],
                     "%.2f" % r["parent_identity"], r["n_mut_contacts"],
                     ("%+.3f" % r["delta_vs_shuffled"]) if have_spec else "-",
                     ("%+.3f" % r["delta_vs_baseline"]) if have_spec else "-"))

    ident = [r["parent_identity"] for r in final]
    lens = [r["length"] for r in final]
    write_methods(final, baseline, have_spec, ident, lens)
    return 0


def write_methods(final, baseline, have_spec, ident, lens):
    p = os.path.join(SUB, "challenge5_methods.md")
    spec_line = (
        "Specificity was estimated with the ProteinMPNN peptide-recovery test used by "
        "Liu et al.: with the binder, the MHC and the other nine peptide residues held "
        "fixed, we read the conditional probability of each amino acid at peptide "
        "position 6 and took log P(Asp) - log P(Gly), i.e. G12D over wild type. The "
        "same quantity measured on the bare pMHC with no binder present (%+.3f) is "
        "subtracted, so the reported margin is the binder's own contribution rather "
        "than the groove's intrinsic preference.\n\n"
        "That subtraction alone is not enough, and we checked: a nonsense poly-AEKL "
        "helix threaded onto one of these backbones scored +0.52 against the "
        "bare-groove baseline, so a large part of it rewards mere occupancy of the "
        "site. The score designs are ranked on is therefore measured against a "
        "SHUFFLE of each design's own sequence on its own backbone, holding pose, "
        "length and amino-acid composition fixed and varying only the order. What "
        "survives that subtraction is attributable to the design.\n\n"
        "A second, independent measure was attempted: the ProteinMPNN score of the "
        "binder sequence against the G12D background and against a wild-type "
        "background built by renaming the peptide's Asp6 to Gly on the same backbone. "
        "It did not work. Across the candidate set its entire range was 0.16 log units "
        "(-0.105 to +0.059) against 1.15 for the measure above, and it correlated "
        "-0.03 with it. We report it in the detail table as `binder_fit_delta` and do "
        "NOT use it for ranking: a measure with almost no dynamic range that agrees "
        "with nothing is not corroboration."
        % baseline if (have_spec and baseline is not None) else
        "Specificity was NOT measured for this submission."
    )
    txt = """# Challenge 5 - KRAS G12D / HLA-A*11:01 - methods

## Summary

%d de novo mini-binders designed to read the G12D aspartate of the KRAS peptide
VVVGADGVGK presented by HLA-A*11:01 (PDB 7OW6), while making limited contact with the
MHC itself. Lengths %d-%d aa. No cysteines in any design.

## Target preparation

Chains A (heavy), B (beta-2-microglobulin) and C (peptide) of 7OW6; the TCR chains D
and E were discarded. A crop sweep over eight candidate trimmings compared buried
hydrophobic surface newly exposed by each cut; the chosen 20 A crop around the peptide
(179 heavy-chain residues + the 10-mer peptide = 189) exposed 250 A^2, 1.04x the intact
complex, the least of any option tested. Solvent exposure of the peptide was measured
rather than assumed: Asp6 is 11%% exposed and Gly7 73%%, with the anchor residues 1, 2,
3 and 10 at 1-3%%.

## Backbone generation

Backbones were taken from the 152-scaffold set published with Liu et al., *Design of
high-specificity binders for peptide-MHC-I complexes*, Science 389:386-391 (2025)
(Zenodo 10.5281/zenodo.15169815), each of which is a de novo designed protein shown to
arc over an MHC-I groove, built against pMHC targets **other than ours**. Each scaffold
was docked onto our target by superposing its own MHC platform onto ours with a local
sequence alignment driving the residue correspondence; superposition RMSD was 0.47-1.50
A over matched CA atoms (median 0.63).

Poses were then filtered on criteria fixed before any sequence existed: no backbone
clash with the target, at least 4 heavy-atom contacts to the Asp6 side chain, and at
least 20%% of the interface on the peptide rather than the MHC. 151 of 152 scaffolds
docked; 111 were free of backbone clashes; 75 contacted Asp6; 63 were peptide-centric
and went forward. All 63 contact peptide positions 6, 7 and 8, and only 2 contact
position 1 - the peptide-centric footprint the source paper describes.

## Sequence design

Every residue of each binder was redesigned with ProteinMPNN (v_48_020), with the entire
target chain held fixed, at temperatures 0.1, 0.2 and 0.3, 16 sequences each. Cysteine
was excluded and the output was checked for it rather than the flag being trusted.

## Filtering

Sequences passed a composition screen developed over the four previous challenges:
two-sided limits on charged and hydrophobic fraction and on mean Kyte-Doolittle
hydropathy, limits on identical and same-charge runs, and a proline/glycine gate
(P <= 12%%, G <= 12%%, P+G <= 20%%) applied first.

%s

## Ranking

Ordered by specificity margin first, then contacts to the Asp6 mutation, then the
peptide share of the interface, with the ProteinMPNN score as a tie-break. At most 3
designs share a parent backbone, so the set covers %d distinct backbones.

## De novo declaration

Track 3 requires de novo designs and forbids modifying an existing binder. The
backbones here are de novo designs from the published scaffold set cited above, built
against different pMHC targets, and 100%% of each binder sequence was rewritten for this
target. Sequence identity of each design to the scaffold it came from is %.0f-%.0f%%
(median %.0f%%); designs at or above %.0f%% identity were excluded. Per-design identities
and parent scaffold names are in the accompanying detail table. We state this
explicitly so the organisers can judge eligibility directly.

## Reproducibility

Target preparation, scaffold docking, pose filtering and all screens are
deterministic and reproducible from the public inputs (PDB 7OW6/7OW5 and Zenodo
10.5281/zenodo.15169815). ProteinMPNN sampling is stochastic, and the run that produced
these sequences was made before a fixed seed was added to our script, so the exact
sampled set cannot be regenerated bit-for-bit. All sequences that run produced are
recorded alongside the code, so the set is auditable even though it is not
reproducible; later runs are seeded.

## What we did not do

- No experimental data of any kind.
- Affinity was not predicted. No computational metric used here predicts binding
  strength, and the source paper relies on experiment for affinity.
- The specificity margin is a sequence model's conditional probability on a fixed
  backbone. It is not a KD ratio and should not be read as one.
- The wild-type comparison is implicit in the margin above; we did not build and
  predict a separate wild-type complex for these designs.
- Side chains were never packed or minimised; all interface counts are from the docked
  backbone with ProteinMPNN-assigned identities.
""" % (len(final), min(lens), max(lens), spec_line,
       len({r["parent"] for r in final}),
       100 * min(ident), 100 * max(ident),
       100 * sorted(ident)[len(ident) // 2], 100 * MAX_PARENT_IDENTITY)
    open(p, "w").write(txt)
    print()
    print("wrote %s" % p)


if __name__ == "__main__":
    sys.exit(main())
