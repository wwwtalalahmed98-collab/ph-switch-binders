"""Kaggle - Challenge 5. Rewrite the shortlisted backbones' sequences with ProteinMPNN.

WHY. Hallucination produces a backbone and a sequence at the same time, and the sequence
is the part it does badly. The raw run-2 output looks like this:

    ...GGGPGGGGGGGGCCGHHHHHHHHHHHRRRRRRRRRRRRRRSGSG
    CCGGRGGGGGGGGGFDDWDDWFDDYYDFFWPYNDDDDPGGG...RRRRRRRRRRC
    MKRRRRAAAAAAAAMLAFWAAWWAGYAGY...          <- the pLDDT 0.63 one

Eleven histidines in a row, fourteen arginines, eight alanines. Those are optimiser
artefacts, not proteins, and feeding them to the counter-screen would be scoring noise
carefully. ProteinMPNN holds the backbone fixed and writes a real sequence for it. It is
bundled inside ColabDesign, needs no GPU, and takes minutes.

This is the same gap that went unnoticed for three challenges on TfR1: doing the shape
half of the pipeline and skipping the sequence half.

ALSO NOTE: run 2 passed rm_aa="C" to prep_inputs, the call was accepted without error,
and the sequences still contain cysteines. The argument was silently ignored. MPNN's own
rm_aa is applied here AND the result is checked - if a cysteine survives, the design is
dropped rather than reported as cysteine-free.

THE COMPOSITION SCREEN IS THE ONE FROM TfR1, INCLUDING THE PROLINE GATE. It has been
wrong twice before and each fix is in it:
  - two-sided limits, because a one-sided version scored a 64%-hydrophobic aggregate
    with AAAAAAAA as "good";
  - a proline/glycine gate, because without it a 30%-proline rod scored "good" and 11 of
    12 "good" candidates on TfR1 did not survive adding it.
It is a list of ways to be obviously wrong, not a prediction that anything folds.

SELECTION IS ON COMPOSITION, NOT ON MPNN's OWN SCORE. MPNN's score rewards what MPNN
likes, which on a poor backbone is a charged repeat.
"""
import glob
import json
import os
import subprocess

W = "/kaggle/working"
MIN_ON_PEPTIDE = 1      # a backbone that touches nothing is not worth rewriting
N_SEQ = 24              # per backbone per temperature
TEMPERATURES = (0.1, 0.2)
KEEP_PER_BACKBONE = 4

PEP_G12D = "VVVGADGVGK"
EXPECTED_A = ("SHSMRYFYTSVSRPGGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDQETRNVKAQSQ"
              "TDRVDLGTLRGYYNQSEDGSHTIQIMYGCDVGPDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQI"
              "TKRKWEAAHAAEQQRAYLEGRCVEWLRRYLENGKETLQR")

CHARGED = set("DEKR")
HYDROPHOBIC = set("AVLIMFWYC")
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5,
          L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)
TIERS = [
    ("good", dict(same_run=5, charge_run=8, charged=(0.18, 0.40),
                  hydrophobic=(0.25, 0.42), kd=(-1.5, 0.0))),
    ("marginal", dict(same_run=6, charge_run=14, charged=(0.10, 0.55),
                      hydrophobic=(0.20, 0.50), kd=(-2.2, 0.3))),
]
P_MAX, G_MAX, PG_MAX = 0.12, 0.12, 0.20

T3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
      "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
      "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
      "TYR": "Y", "VAL": "V"}


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
    """Distance from natural composition: 25% charged, 32% hydrophobic, kd about -0.6."""
    return round(abs(c["charged_frac"] - 0.25) / 0.25
                 + abs(c["hydrophobic_frac"] - 0.32) / 0.32
                 + abs(c["kd"] + 0.6) / 1.0
                 + max(0, c["same_run"] - 5) / 5.0
                 + max(0, c["charge_run"] - 8) / 8.0, 3)


def main():
    try:
        import colabdesign  # noqa: F401
    except ImportError:
        subprocess.run("pip -q install git+https://github.com/sokrypton/ColabDesign.git"
                       "@v1.1.1 biopython", shell=True, check=True)

    import numpy as np
    from Bio import PDB

    def residues(ch):
        out, seen = [], set()
        for r in ch:
            if r.id[0] == " " and r.get_resname() in T3 and r.id[1] not in seen:
                seen.add(r.id[1])
                out.append(r)
        return out

    def seq(rs):
        return "".join(T3[r.get_resname()] for r in rs)

    # ---- find the generation records and the backbones, and say what was found ----
    recs, srcs = [], []
    for pat in ("pmhc_gen.json", "pmhc_gen2.json"):
        for f in sorted(glob.glob("/kaggle/input/**/%s" % pat, recursive=True)) + \
                 (["%s/%s" % (W, pat)] if os.path.exists("%s/%s" % (W, pat)) else []):
            got = [d for d in json.load(open(f)) if d.get("sequence")]
            recs += got
            srcs.append("%s (%d)" % (f, len(got)))
            break
    pdbs = {os.path.basename(p)[:-4]: p
            for p in glob.glob("/kaggle/input/**/*.pdb", recursive=True)
            if os.path.basename(p).startswith("pmhc")}
    print("records: %s" % ("; ".join(srcs) if srcs else "NONE FOUND"))
    print("backbone structures found: %d" % len(pdbs))
    if not recs:
        raise SystemExit("no pmhc_gen*.json under /kaggle/input - attach the generation "
                         "notebook's output as an input")
    if not pdbs:
        raise SystemExit("no pmhc*.pdb under /kaggle/input - the records are there but "
                         "the structures are not, and MPNN needs the structures")

    seen, uniq = set(), []
    for d in recs:
        if d["name"] not in seen:
            seen.add(d["name"])
            uniq.append(d)
    recs = uniq

    cands = [d for d in recs if d.get("on_peptide", 0) >= MIN_ON_PEPTIDE
             and d["name"] in pdbs]
    cands.sort(key=lambda d: (-d.get("on_mutation", 0), -d.get("on_peptide", 0),
                              -d.get("plddt", 0)))
    print("\n%d of %d designs touch the peptide and have a structure; rewriting those."
          % (len(cands), len(recs)))
    missing = [d["name"] for d in recs
               if d.get("on_peptide", 0) >= MIN_ON_PEPTIDE and d["name"] not in pdbs]
    if missing:
        print("  (%d qualify but have no structure attached: %s)"
              % (len(missing), ", ".join(missing[:6])))
    if not cands:
        raise SystemExit("nothing to rewrite")

    from colabdesign.mpnn import mk_mpnn_model
    mpnn = mk_mpnn_model(model_name="v_48_020")
    print("ProteinMPNN v_48_020 loaded\n")

    out = []
    for d in cands:
        path = pdbs[d["name"]]
        model = PDB.PDBParser(QUIET=True).get_structure("c", path)[0]
        chains = {c.id: residues(c) for c in model}
        tgt = next((c for c, rs in chains.items()
                    if seq(rs) in (EXPECTED_A + PEP_G12D, PEP_G12D + EXPECTED_A)), None)
        if tgt is None:
            print("%s: no chain reproduces the target, skipped (chains %s)"
                  % (d["name"], {c: len(r) for c, r in chains.items()}))
            continue
        others = [c for c in chains if c != tgt]
        if len(others) != 1:
            print("%s: expected one binder chain, got %s" % (d["name"], others))
            continue
        binder = others[0]
        try:
            mpnn.prep_inputs(pdb_filename=path, chain="%s,%s" % (tgt, binder),
                             fix_pos=tgt, rm_aa="C")
        except Exception as e:
            print("%s: prep failed, %s %s" % (d["name"], type(e).__name__, str(e)[:110]))
            continue

        orig = d["sequence"]
        rows = []
        for temp in TEMPERATURES:
            s = mpnn.sample(num=N_SEQ, batch=1, temperature=temp)
            seqs = [x.split("/")[-1] if "/" in x else x
                    for x in np.asarray(s["seq"]).reshape(-1)]
            scores = np.asarray(s["score"]).reshape(-1)
            for x, sc in zip(seqs, scores):
                if len(x) != len(orig) or "C" in x:
                    continue          # rm_aa was ignored once already; verify, do not trust
                c = composition(x)
                tier, why = tier_of(c)
                rows.append(dict(sequence=x, mpnn_score=round(float(sc), 4),
                                 temperature=temp, tier=tier, why=why,
                                 distance=distance(c), **c))
        rows.sort(key=lambda r: (r["tier"] in ("fail", "fold-breaker"), r["distance"],
                                 r["mpnn_score"]))
        keep = [r for r in rows if r["tier"] in ("good", "marginal")][:KEEP_PER_BACKBONE]
        oc = composition(orig)
        counts = {t: sum(1 for r in rows if r["tier"] == t)
                  for t in ("good", "marginal", "fail", "fold-breaker")}
        print("%-22s on-mut %d pept %2d | original [%s] d=%.2f | %d sampled: "
              "%d good %d marginal %d fail %d fold-breaker | best d=%s"
              % (d["name"], d.get("on_mutation", 0), d.get("on_peptide", 0),
                 tier_of(oc)[0], distance(oc), len(rows), counts["good"],
                 counts["marginal"], counts["fail"], counts["fold-breaker"],
                 ("%.2f" % keep[0]["distance"]) if keep else "none"))
        for i, r in enumerate(keep, 1):
            r.update(name="%s_m%02d" % (d["name"], i), parent=d["name"],
                     length=len(r["sequence"]), parent_on_mutation=d.get("on_mutation"),
                     parent_on_peptide=d.get("on_peptide"),
                     parent_plddt=d.get("plddt"), parent_i_ptm=d.get("i_ptm"),
                     parent_sequence=orig, structure=os.path.basename(path))
            out.append(r)

    json.dump(out, open("%s/pmhc_mpnn.json" % W, "w"), indent=0)
    print("\nwrote %s/pmhc_mpnn.json" % W)

    print("\n==================== REWRITE RESULT ====================")
    if not out:
        raise SystemExit("nothing survived composition on any backbone")
    good = [r for r in out if r["tier"] == "good"]
    print("%d candidates from %d backbones: %d good, %d marginal"
          % (len(out), len(set(r["parent"] for r in out)), len(good),
             len(out) - len(good)))
    print("\n%-26s %5s %6s %6s %6s %6s %6s %6s %10s"
          % ("candidate", "p-mut", "chg", "phob", "Pro", "Gly", "kd", "dist", "tier"))
    for r in sorted(out, key=lambda r: (r["tier"] != "good", r["distance"]))[:30]:
        print("%-26s %5d %5.0f%% %5.0f%% %5.0f%% %5.0f%% %+6.2f %6.2f %10s"
              % (r["name"], r["parent_on_mutation"] or 0, 100 * r["charged_frac"],
                 100 * r["hydrophobic_frac"], 100 * r["pro_frac"], 100 * r["gly_frac"],
                 r["kd"], r["distance"], r["tier"]))
    print("\nNone of these has been folded. The contacts quoted belong to the PARENT")
    print("backbone; whether a rewritten sequence still binds there is what the")
    print("counter-screen has to decide.")
    print("\n===== COPY FROM HERE =====")
    for r in out:
        print(json.dumps(r))
    print("===== COPY TO HERE =====")


main()
