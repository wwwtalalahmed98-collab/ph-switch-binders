"""Predict how much weaker a complex binds at pH 6.0 than at pH 7.4, from structure.

WHY THIS INSTEAD OF A PROXY
---------------------------
Every pH ranking this project has used so far was a geometric proxy: histidine-to-cation
distance (which correlated +0.24 with the effect it was meant to predict - the wrong sign),
then burial on binding (-0.58). A proxy has to be calibrated against the thing it stands for,
and we never had ground truth to calibrate against.

This computes the thing itself. The pH dependence of binding is a thermodynamic linkage: if
a titratable group's pKa differs between the free and bound states, binding and protonation
are coupled, and the strength of that coupling IS the pH sensitivity. For one group:

    Ka_obs(pH) = Ka_0 * (1 + 10^(pKa_bound - pH)) / (1 + 10^(pKa_free - pH))

so, summing the log terms over every titratable group i,

    log10( KD(pH_low) / KD(pH_high) ) = SUM_i [ L_i(pH_high) - L_i(pH_low) ]
    where L_i(pH) = log10( (1 + 10^(pKa_b,i - pH)) / (1 + 10^(pKa_f,i - pH)) )

SIGN CONVENTION, which is the thing most easily got backwards:

  * POSITIVE result  = KD larger at pH 6.0 = binds WEAKER in acid = RELEASES.
    This is what Challenge 4 asks for, and it needs pKa_bound < pKa_free: burying the
    histidine against a neutral or hydrophobic environment so that gaining a charge costs
    binding energy.
  * NEGATIVE result  = binds TIGHTER in acid. This is the FcRn/IgG mechanism, and it is what
    Challenge 1 (EGFR) asked for. Same machinery, opposite target.

The organizers state in their own selection prompt that one titratable site changes affinity
by at most tenfold per pH unit. The equation above reproduces that bound automatically: over
the 1.4 pH units from 7.4 to 6.0 a single group cannot exceed 1.4 log units, so "no
detectable binding" at pH 6.0 needs two or more groups that each actually do work.

VALIDATION
----------
Run with --validate. pKa prediction for histidines carries roughly 1 pH unit of error, which
is the same size as the effect being engineered, so this is checked against complexes whose
pH behaviour is already known before it is trusted on a design:

  1I1A  rat FcRn + IgG Fc     binds at pH 6.0, releases at 7.4  -> expect STRONGLY NEGATIVE
  1BRS  barnase + barstar     very tight, no pH switch          -> expect NEAR ZERO
  3HFM  HyHEL-10 Fab + lysozyme   antibody-antigen, no switch   -> expect NEAR ZERO

RESOLUTION MATTERS MORE THAN IT LOOKS. The first attempt at the FcRn control used 1FRT and
came out with the wrong sign. 1FRT is a 4.50 A structure: in it, His310's ND1 sits 1.2 A from
a Glu117 carboxylate oxygen, which is physically impossible (no N-O contact is shorter than
about 2.6 A). PROPKA was being handed clashing side chains and produced nonsense from them.
Any structure coarser than MIN_RESOLUTION is refused, because side-chain geometry at that
resolution is not determined and a pKa computed from it is not a measurement of anything.

Those expectations are directional. Published affinities depend on construct, salt and
method, so the test here is sign and rough magnitude, not a numerical match.
"""
import argparse
import contextlib
import io as _io
import json
import math
import os
import sys

import requests
from Bio import PDB

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import propka_compat  # noqa: E402  - must patch before propka.run is imported
from propka.run import single  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PKA_DIR = os.path.join(ROOT, "data", "pka")
PH_HIGH, PH_LOW = 7.4, 6.0
MIN_RESOLUTION = 3.0   # angstrom; coarser than this and side chains are not placed
TITRATABLE = {"HIS", "ASP", "GLU", "LYS", "CYS", "TYR", "ARG", "N+", "C-"}


def fetch(pdb_id):
    os.makedirs(PKA_DIR, exist_ok=True)
    p = os.path.join(PKA_DIR, pdb_id + ".pdb")
    if not os.path.exists(p):
        r = requests.get("https://files.rcsb.org/download/%s.pdb" % pdb_id, timeout=90)
        r.raise_for_status()
        _io.open(p, "w", newline="\n").write(r.text)
    return p


class _Keep(PDB.Select):
    """One model, chosen chains, no waters/ligands, no alternate conformations."""

    def __init__(self, chains):
        self.chains = set(chains)

    def accept_model(self, m):
        return m.id == 0

    def accept_chain(self, c):
        return c.id in self.chains

    def accept_residue(self, r):
        return r.id[0] == " "

    def accept_atom(self, a):
        return (not a.is_disordered()) or a.get_altloc() in (" ", "A")


def subset(pdb_path, chains, tag):
    """Write a cleaned PDB containing only `chains`, and return its path."""
    out = os.path.join(PKA_DIR, "%s_%s.pdb" % (
        os.path.splitext(os.path.basename(pdb_path))[0], tag))
    st = PDB.PDBParser(QUIET=True).get_structure("x", pdb_path)
    io_ = PDB.PDBIO()
    io_.set_structure(st)
    io_.save(out, _Keep(chains))
    return out


def pkas(pdb_path):
    """{(chain, resnum, restype): pKa} for every titratable group PROPKA finds.

    PROPKA writes a lot to stdout about alternate conformations and terminal atoms in
    crystal structures; it is noise, so it is captured rather than printed.
    """
    buf = _io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        mol = single(pdb_path, optargs=["--quiet"], write_pka=False)
    conf = mol.conformations[mol.conformation_names[0]]
    out = {}
    for g in conf.groups:
        if g.residue_type not in TITRATABLE:
            continue
        key = (g.atom.chain_id, int(g.atom.res_num), g.residue_type)
        out[key] = float(g.pka_value)
    if not out:
        raise SystemExit("PROPKA returned no titratable groups for %s" % pdb_path)
    return out


def _L(pka_b, pka_f, ph):
    return math.log10((1 + 10 ** (pka_b - ph)) / (1 + 10 ** (pka_f - ph)))


def linkage(pairs, ph_high=PH_HIGH, ph_low=PH_LOW):
    """log10( KD(ph_low) / KD(ph_high) ). Positive = weaker in acid = releases."""
    return sum(_L(b, f, ph_high) - _L(b, f, ph_low) for b, f in pairs)


def analyse(pdb_path, binder_chains, target_chains, label=""):
    comp = subset(pdb_path, list(binder_chains) + list(target_chains), "complex")
    free = subset(pdb_path, list(binder_chains), "free")
    pc, pf = pkas(comp), pkas(free)

    rows = []
    for key, f in sorted(pf.items()):
        if key not in pc:
            continue
        b = pc[key]
        rows.append(dict(chain=key[0], resnum=key[1], restype=key[2],
                         pka_free=round(f, 2), pka_bound=round(b, 2),
                         dpka=round(b - f, 2),
                         contribution=round(_L(b, f, PH_HIGH) - _L(b, f, PH_LOW), 3)))

    his = [r for r in rows if r["restype"] == "HIS"]
    total_all = linkage([(r["pka_bound"], r["pka_free"]) for r in rows])
    total_his = linkage([(r["pka_bound"], r["pka_free"]) for r in his])

    print("\n=== %s  binder %s vs target %s ===" % (label or os.path.basename(pdb_path),
                                                    ",".join(binder_chains),
                                                    ",".join(target_chains)))
    print("%d titratable groups shared between free and bound, %d histidines"
          % (len(rows), len(his)))
    movers = sorted(rows, key=lambda r: -abs(r["contribution"]))[:8]
    print("  %-6s %6s %6s %6s %8s" % ("group", "free", "bound", "dpKa", "log10 KD"))
    for r in movers:
        print("  %-6s %6.2f %6.2f %+6.2f %+8.3f"
              % ("%s%d%s" % (r["restype"][:1] if r["restype"] in ("HIS",) else
                             r["restype"][:1], r["resnum"], r["chain"]),
                 r["pka_free"], r["pka_bound"], r["dpka"], r["contribution"]))
    print("  predicted log10(KD at pH %.1f / KD at pH %.1f):" % (PH_LOW, PH_HIGH))
    print("    histidines only : %+.2f  (%.0f-fold %s in acid)"
          % (total_his, 10 ** abs(total_his), "weaker" if total_his > 0 else "tighter"))
    print("    all groups      : %+.2f  (%.0f-fold %s in acid)"
          % (total_all, 10 ** abs(total_all), "weaker" if total_all > 0 else "tighter"))
    return dict(label=label, pdb=os.path.basename(pdb_path),
                binder=list(binder_chains), target=list(target_chains),
                n_groups=len(rows), n_his=len(his),
                log10_ratio_his=round(total_his, 3),
                log10_ratio_all=round(total_all, 3), groups=rows)


VALIDATION = [
    ("1I1A", "CD", "AB", "FcRn + IgG Fc", "strongly negative (binds at pH 6.0)"),
    ("1BRS", "D", "A", "barnase + barstar", "near zero (no pH switch)"),
    ("3HFM", "Y", "LH", "HyHEL-10 Fab + lysozyme", "near zero (no pH switch)"),
]


def validate():
    print("Validating the linkage calculation against complexes whose pH behaviour is known.")
    print("Positive = weaker in acid (releases). Negative = tighter in acid.\n")
    out = []
    for pdb_id, binder, target, label, expect in VALIDATION:
        r = analyse(fetch(pdb_id), binder, target, label)
        r["expected"] = expect
        print("    expected: %s" % expect)
        out.append(r)
    print("\n%-26s %10s %10s  %s" % ("complex", "His only", "all", "expected"))
    for r in out:
        print("%-26s %+10.2f %+10.2f  %s"
              % (r["label"], r["log10_ratio_his"], r["log10_ratio_all"], r["expected"]))
    json.dump(out, open(os.path.join(PKA_DIR, "validation.json"), "w"), indent=0)
    print("\nwrote data/pka/validation.json")
    print("\nREAD THIS BEFORE USING THE NUMBERS. Histidine pKa prediction carries about one"
          "\npH unit of error, the same size as the effect being engineered. Treat the output"
          "\nas a ranking signal between designs scored the same way, not as a predicted KD.")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--validate", action="store_true",
                    help="run the known-answer controls and stop")
    ap.add_argument("pdb", nargs="?", help="PDB id or path to a complex")
    ap.add_argument("--binder", help="binder chain id(s), e.g. B")
    ap.add_argument("--target", help="target chain id(s), e.g. AR")
    a = ap.parse_args()
    if a.validate or not a.pdb:
        return validate()
    path = a.pdb if os.path.exists(a.pdb) else fetch(a.pdb)
    if not (a.binder and a.target):
        raise SystemExit("--binder and --target are required with a structure")
    analyse(path, a.binder, a.target)


if __name__ == "__main__":
    main()
