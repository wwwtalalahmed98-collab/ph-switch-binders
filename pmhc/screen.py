"""The buildability screen, in one place.

This is the composition screen built over challenges 1-4 and corrected twice. It has
been copy-pasted into four Kaggle scripts, which is how the proline blind spot survived
as long as it did - a fix in one copy was not a fix in the others. Both the local
scripts and the Kaggle scripts now import this module instead.

WHAT IT CHECKS, and why each one is there:

  Proline and glycine, checked FIRST. On the TfR1 challenge, 4 of 12 sequences that
  passed every other check were 25-32% proline. A polyproline rod has unremarkable
  charge and hydrophobicity, so it sails through composition limits while being
  physically incapable of folding into a binder. The gate is P <= 12%, G <= 12%,
  P+G <= 20%.

  Two-sided limits on charge and hydrophobicity. An earlier one-sided version passed a
  64%-hydrophobic sequence, which is an aggregate, not a protein. Being far below
  natural composition is as diagnostic as being far above it.

  Runs. Long identical or same-charge stretches are what a sequence design model
  produces when the backbone underneath it is not a real fold.

  Kyte-Doolittle mean hydropathy. A single number that catches sequences which are
  individually in range on every fraction but globally wrong.

WHAT IT DOES NOT CHECK. Nothing here says anything about whether a sequence binds
anything. It is a filter for "could this be expressed and could it fold", applied before
spending GPU time. Passing it is necessary, not sufficient, and a "good" tier is not a
prediction of success.
"""

CHARGED = set("DEKR")
HYDROPHOBIC = set("AVLIMFWYC")

# Kyte-Doolittle hydropathy.
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5,
          L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)

TIERS = [
    ("good", dict(same_run=5, charge_run=8, charged=(0.18, 0.40),
                  hydrophobic=(0.25, 0.42), kd=(-1.5, 0.0))),
    ("marginal", dict(same_run=6, charge_run=14, charged=(0.10, 0.55),
                      hydrophobic=(0.20, 0.50), kd=(-2.2, 0.3))),
]
P_MAX, G_MAX, PG_MAX = 0.12, 0.12, 0.20

# Natural-composition reference point, used only for ranking within a tier.
NAT_CHARGED, NAT_HYDRO, NAT_KD = 0.25, 0.32, -0.6


def same_run(s):
    best = cur = 1
    for i in range(1, len(s)):
        cur = cur + 1 if s[i] == s[i - 1] else 1
        best = max(best, cur)
    return best


def class_run(s, group):
    best = cur = 0
    for c in s:
        cur = cur + 1 if c in group else 0
        best = max(best, cur)
    return best


def composition(seq):
    n = float(len(seq))
    return dict(charged_frac=round(sum(c in CHARGED for c in seq) / n, 3),
                hydrophobic_frac=round(sum(c in HYDROPHOBIC for c in seq) / n, 3),
                charge_run=class_run(seq, CHARGED), same_run=same_run(seq),
                pro_frac=round(seq.count("P") / n, 3),
                gly_frac=round(seq.count("G") / n, 3),
                kd=round(sum(KD[c] for c in seq) / n, 3))


def _fails(c, lim):
    why = []
    if c["same_run"] > lim["same_run"]:
        why.append("%d identical" % c["same_run"])
    if c["charge_run"] > lim["charge_run"]:
        why.append("%d charged in a row" % c["charge_run"])
    lo, hi = lim["charged"]
    if c["charged_frac"] > hi:
        why.append("%.0f%% charged" % (100 * c["charged_frac"]))
    elif c["charged_frac"] < lo:
        why.append("only %.0f%% charged" % (100 * c["charged_frac"]))
    lo, hi = lim["hydrophobic"]
    if c["hydrophobic_frac"] > hi:
        why.append("%.0f%% hydrophobic" % (100 * c["hydrophobic_frac"]))
    elif c["hydrophobic_frac"] < lo:
        why.append("only %.0f%% hydrophobic" % (100 * c["hydrophobic_frac"]))
    lo, hi = lim["kd"]
    if c["kd"] > hi:
        why.append("kd %+.2f net hydrophobic" % c["kd"])
    elif c["kd"] < lo:
        why.append("kd %+.2f extremely polar" % c["kd"])
    return why


def tier_of(c):
    """Proline/glycine first: a polyproline rod passes every other check."""
    if c["pro_frac"] > P_MAX or c["gly_frac"] > G_MAX or \
            (c["pro_frac"] + c["gly_frac"]) > PG_MAX:
        return "fold-breaker", ("%.0f%% Pro, %.0f%% Gly"
                                % (100 * c["pro_frac"], 100 * c["gly_frac"]))
    for tier, lim in TIERS:
        if not _fails(c, lim):
            return tier, ""
    return "fail", "; ".join(_fails(c, TIERS[-1][1]))


def distance(c):
    """Distance from natural composition. Ranking only, never a pass/fail."""
    return round(abs(c["charged_frac"] - NAT_CHARGED) / NAT_CHARGED
                 + abs(c["hydrophobic_frac"] - NAT_HYDRO) / NAT_HYDRO
                 + abs(c["kd"] - NAT_KD) / 1.0
                 + max(0, c["same_run"] - 5) / 5.0
                 + max(0, c["charge_run"] - 8) / 8.0, 3)


def identity(a, b):
    """Sequence identity over the aligned prefix, as a fraction.

    Used for the Track 3 "de novo only" declaration: every design's identity to the
    scaffold its backbone came from is measured and reported. Designs are the same
    length as their parent here (ProteinMPNN keeps the backbone), so a positional
    comparison is the right one and no alignment is needed.
    """
    if not a or not b:
        return 0.0
    n = min(len(a), len(b))
    return round(sum(1 for i in range(n) if a[i] == b[i]) / float(n), 4)


def report(seq):
    """Everything the screen has to say about one sequence."""
    c = composition(seq)
    tier, why = tier_of(c)
    c["tier"] = tier
    c["why"] = why
    c["distance"] = distance(c)
    c["length"] = len(seq)
    return c
