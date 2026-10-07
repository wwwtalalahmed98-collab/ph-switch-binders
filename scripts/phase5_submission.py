"""Phase 5: assemble the competition submission.

Produces the ranked CSV Proteinbase expects (name, sequence, molecule_class) plus a
metrics table for the supporting metadata package.

Ranking follows the challenge's own stated priority, not ours:
    1. pH-selectivity   2. mouse cross-reactivity   3. affinity
The FAQ is explicit that a weak but clearly pH-sensitive binder may be considered more
impactful than a high-affinity binder that is not pH-sensitive, so switch evidence
outranks predicted affinity here by design.

Usage:
    python scripts/phase5_submission.py designs.json --out submission/
"""
import argparse
import csv
import json
import os

MAX_DESIGNS = 20          # Track 3 limit
MIN_LEN, MAX_LEN = 10, 250
NOVELTY_LIMIT = 30.0      # Proteinbase de novo identity threshold


def molecule_class(seq):
    """Returns (molecule_class for the CSV, reporting category from FAQ 2).

    The submission form and the FAQ disagree: FAQ 5 says the value is "protein",
    the upload template and the form itself say "single_chain". The template is
    what the uploader parses, so it wins. The FAQ 2 category is carried alongside
    for the metrics table, since the organisers stratify results by it.
    """
    n = len(seq)
    if n < 40:
        return "single_chain", "microbinder (<40 aa)"
    if n <= 100:
        return "single_chain", "minibinder (40-100 aa)"
    return "single_chain", "large protein binder (>100 aa)"


def validate(d):
    """Hard rules that would otherwise cost us the design at screening."""
    problems = []
    seq = d["sequence"].strip().upper()
    if not (MIN_LEN <= len(seq) <= MAX_LEN):
        problems.append(f"length {len(seq)} outside {MIN_LEN}-{MAX_LEN}")
    if set(seq) - set("ACDEFGHIKLMNPQRSTVWY"):
        problems.append("non-standard amino acid")
    ident = d.get("max_identity")
    if ident is not None and ident > NOVELTY_LIMIT:
        problems.append(f"novelty {ident:.1f}% > {NOVELTY_LIMIT}% (de novo rule)")
    if d.get("n_his_contacts") is not None and d.get("role") == "switch":
        if d["n_his_contacts"] < 2:
            problems.append("switch design with <2 His-carboxylate contacts")
    # A histidine substitution that guts the interface is worse than no switch at all:
    # the design then scores zero on pH-selectivity AND on affinity. Its own control
    # is the better use of the slot.
    ret = d.get("iptm_retained")
    if ret is not None and ret < 0.50:
        problems.append(f"His substitutions cost {100 * (1 - ret):.0f}% of parent i_ptm")
    iptm = d.get("i_ptm")
    if iptm is not None and iptm < 0.35:
        problems.append(f"i_ptm {iptm:.2f} below the 0.35 floor for a plausible binder")
    return problems


def rank_key(d):
    """Sort best-first on the challenge's stated priority order.

    Epitope conservation is identical across our designs (all target the same fully
    conserved patch), so it does not discriminate here; cross-species performance was
    designed in at epitope selection rather than ranked for. The discriminating signals
    are therefore switch evidence first, then structural quality.

    i_ptm, not pLDDT, carries the binding signal: across round 1's designs the two
    correlated at r = 0.048, so pLDDT cannot stand in for interface quality. It is
    kept as the last tiebreak because a binder AlphaFold cannot fold confidently is
    an expression risk in the wet lab even when its interface scores well.
    """
    return (
        -(d.get("n_his_contacts") or 0),          # pH switch evidence first
        -(d.get("switch_score") or 0.0),
        -(d.get("i_ptm") or 0.0),                 # predicted binding
        -(d.get("plddt") or 0.0),                 # fold confidence / expression risk
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("designs", help="JSON list of design dicts")
    ap.add_argument("--out", default="submission")
    ap.add_argument("--max-per-backbone", type=int, default=2,
                    help="cap designs drawn from any one backbone (diversity)")
    ap.add_argument("--allow-invalid", action="store_true",
                    help="include designs that fail validation (flagged in the report)")
    a = ap.parse_args()

    designs = json.load(open(a.designs))
    os.makedirs(a.out, exist_ok=True)

    for d in designs:
        d["sequence"] = d["sequence"].strip().upper()
        d["_problems"] = validate(d)

    ok = [d for d in designs if not d["_problems"]]
    bad = [d for d in designs if d["_problems"]]
    pool = designs if a.allow_invalid else ok

    # Backbone diversity: ProteinMPNN produces 8 near-identical sequences per backbone,
    # so an unconstrained top-20 would be a handful of folds repeated. Spreading the 20
    # slots across distinct backbones hedges against any single fold simply not working,
    # which is the dominant failure mode in published binder campaigns.
    selected, per_backbone = [], {}
    for d in sorted(pool, key=rank_key):
        bb = d.get("source_backbone", d["name"])
        if per_backbone.get(bb, 0) >= a.max_per_backbone:
            continue
        selected.append(d)
        per_backbone[bb] = per_backbone.get(bb, 0) + 1
        if len(selected) >= MAX_DESIGNS:
            break

    # If diversity capping left slots unfilled, top up with the best remaining.
    if len(selected) < MAX_DESIGNS:
        chosen = {id(d) for d in selected}
        for d in sorted(pool, key=rank_key):
            if id(d) in chosen:
                continue
            selected.append(d)
            if len(selected) >= MAX_DESIGNS:
                break
    pool = selected

    csv_path = os.path.join(a.out, "submission.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "sequence", "molecule_class"])
        for d in pool:
            cls, _ = molecule_class(d["sequence"])
            w.writerow([d["name"], d["sequence"], cls])

    metrics_path = os.path.join(a.out, "design_metrics.csv")
    fields = ["rank", "name", "length", "category", "role", "epitope",
              "n_his_contacts", "switch_score", "his_substitutions",
              "epitope_conservation", "plddt", "sc_rmsd", "mpnn_score",
              "max_identity", "generator", "source_backbone"]
    with open(metrics_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for i, d in enumerate(pool, 1):
            _, cat = molecule_class(d["sequence"])
            row = dict(d)
            row.update(rank=i, length=len(d["sequence"]), category=cat)
            w.writerow(row)

    print(f"{len(pool)} designs written to {csv_path}")
    print(f"metrics -> {metrics_path}")

    n_switch = sum(1 for d in pool if d.get("role") == "switch")
    n_ctrl = sum(1 for d in pool if d.get("role") == "control")
    cats = {}
    for d in pool:
        cats[molecule_class(d["sequence"])[1]] = \
            cats.get(molecule_class(d["sequence"])[1], 0) + 1
    print(f"\ncomposition: {n_switch} switch designs, {n_ctrl} pH-insensitive controls")
    for c, n in sorted(cats.items()):
        print(f"   {c}: {n}")
    gens = {}
    for d in pool:
        gens[d.get("generator", "unknown")] = gens.get(d.get("generator", "unknown"), 0) + 1
    print("   generators:", ", ".join(f"{k}={v}" for k, v in sorted(gens.items())))
    top = max(gens.values()) / len(pool) if pool else 0
    if top > 0.5:
        print("   NOTE: one generator contributes >50% of designs. Anthropic's own "
              "protocol capped this to avoid method-specific failure modes.")

    if bad:
        print(f"\n{len(bad)} design(s) failed validation and were EXCLUDED:")
        for d in bad:
            print(f"   {d['name']}: {'; '.join(d['_problems'])}")


if __name__ == "__main__":
    main()
