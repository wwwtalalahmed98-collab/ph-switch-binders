"""Two things the organisers' own guide says matter, which we had not computed.

(2) ESM2 PSEUDO-LOG-LIKELIHOOD. Adaptyv's QC guide reports ROC AUCs for predicting
whether a design binds: ipTM 0.64, PAE 0.65, ESM2 PLL 0.72 - "best single metric" -
and a composite around 0.75. We have never computed the best one. Their threshold is
written as "PLL > 0.0", which cannot be a raw log-likelihood, since those are always
negative; they must normalise it somehow and they do not say how. So this reports BOTH
the raw mean per-residue log-likelihood AND the excess over the same sequence shuffled,
which is a quantity that is genuinely positive when a sequence is more plausible than
its own composition would give by chance. We use our own definition and say so rather
than guess theirs.

(3) PAE_INTERACTION, RESCUED FROM A RUN WE ALREADY PAID FOR. The three-state AlphaFold
screen stored a full PAE matrix for every one of its 156 predictions and we then ranked
on i-pTM, which is the weaker signal: Bennett et al. (Nat Commun 2023) filter on
pae_interaction < 10, and Adaptyv list < 12 standard, < 10 stringent. No GPU is needed
to go back and read those matrices.

THE POINT OF THIS IS THE CONTROLS, NOT THE DESIGNS. That run carried a TCR known to
bind this exact pMHC and a shuffled sequence. On i-pTM they differed by +0.024 against
0.022 of noise, which is why that screen was declared blind. If pae_interaction
separates them, the five GPU-hours come back and we have a ruler. If it does not, we
have confirmed the blindness on a second metric and we stop pretending.

Caveat stated up front: these predictions were made WITHOUT an AF2 initial guess, which
is the protocol the pae_interaction < 10 threshold was calibrated on. A null result
here is evidence about our pipeline, not about pae_interaction.

Inputs to attach: the AF2 screen's notebook output (for *_pae.npy), plus any CSVs of
sequences to score.

Run:  python cellPMHCm_qc.py
"""
import csv
import glob
import json
import os
import random
import re
import subprocess
import sys
import time

W = "/kaggle/working"
T0 = time.time()
ESM_MODEL = "esm2_t33_650M_UR50D"
N_SHUFFLES = 3


def mins():
    return (time.time() - T0) / 60.0


def banner(s):
    print()
    print("=" * 72)
    print("%s   [t+%.1f min]" % (s, mins()))
    print("=" * 72, flush=True)


def sh(cmd, check=False):
    print("$ %s" % cmd, flush=True)
    r = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, text=True)
    if r.stdout.strip():
        print(r.stdout[-3000:], flush=True)
    if check and r.returncode:
        raise SystemExit("failed: %s" % cmd)
    return r.stdout


# --------------------------------------------------------------------------
def load_sequences():
    banner("SEQUENCES")
    found = sorted(glob.glob("/kaggle/input/**/*.csv", recursive=True))
    print("csv files under /kaggle/input: %d" % len(found))
    seqs = {}
    for p in found:
        try:
            rows = list(csv.DictReader(open(p)))
        except OSError:
            continue
        if not rows:
            continue
        cols = list(rows[0].keys())
        scol = next((c for c in cols if c.lower() in ("sequence", "seq")), None)
        ncol = next((c for c in cols if c.lower() in ("name", "id")), None)
        if scol is None:
            continue
        n0 = len(seqs)
        for i, r in enumerate(rows, 1):
            s = (r[scol] or "").strip().upper()
            if not s or not set(s) <= set("ACDEFGHIKLMNPQRSTVWY"):
                continue
            nm = (r[ncol].strip() if ncol and r.get(ncol) else "row_%03d" % i)
            seqs.setdefault(nm, s)
        print("    %-70s +%d" % (p[-70:], len(seqs) - n0))
    print("unique named sequences: %d" % len(seqs))
    return seqs


def esm_pll(seqs):
    banner("ESM2 PSEUDO-LOG-LIKELIHOOD")
    try:
        import torch
        import esm
    except ImportError:
        sh("pip install -q fair-esm")
        import torch
        import esm
    print("torch %s  cuda %s" % (torch.__version__, torch.cuda.is_available()))
    model, alphabet = getattr(esm.pretrained, ESM_MODEL)()
    model.eval()
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.to(dev)
    bc = alphabet.get_batch_converter()

    def score(batch):
        """Mean per-residue log-likelihood, one forward pass per sequence."""
        _, _, toks = bc([("x%d" % i, s) for i, s in enumerate(batch)])
        toks = toks.to(dev)
        with torch.no_grad():
            lg = model(toks)["logits"].log_softmax(-1)
        out = []
        for i, s in enumerate(batch):
            # token 0 is BOS; residues occupy 1..len(s)
            idx = toks[i, 1:len(s) + 1]
            lp = lg[i, 1:len(s) + 1].gather(-1, idx.unsqueeze(-1)).squeeze(-1)
            out.append(float(lp.mean()))
        return out

    names = sorted(seqs)
    rnd = random.Random(20261108)
    res = {}
    B = 4
    for i in range(0, len(names), B):
        chunk = names[i:i + B]
        real = score([seqs[n] for n in chunk])
        nulls = [[] for _ in chunk]
        for _ in range(N_SHUFFLES):
            sh_seqs = []
            for n in chunk:
                t = list(seqs[n])
                rnd.shuffle(t)
                sh_seqs.append("".join(t))
            for j, v in enumerate(score(sh_seqs)):
                nulls[j].append(v)
        for j, n in enumerate(chunk):
            null = sum(nulls[j]) / len(nulls[j])
            res[n] = {"pll": round(real[j], 4),
                      "pll_shuffled": round(null, 4),
                      "pll_excess": round(real[j] - null, 4)}
        if (i // B) % 5 == 0:
            print("    %d/%d  [t+%.1f min]" % (i + len(chunk), len(names), mins()),
                  flush=True)
    vals = sorted(r["pll_excess"] for r in res.values())
    print()
    print("excess over own shuffle: min %+.3f  median %+.3f  max %+.3f"
          % (vals[0], vals[len(vals) // 2], vals[-1]))
    print("above zero: %d of %d" % (sum(1 for v in vals if v > 0), len(vals)))
    print()
    print("NOTE: 'excess over own shuffle' is our definition, not Adaptyv's. Their")
    print("guide writes the threshold as PLL > 0 without saying how they normalise,")
    print("and a raw log-likelihood is always negative. Do not compare the two.")
    return res


# --------------------------------------------------------------------------
def pae_interaction(seqs):
    banner("PAE_INTERACTION, FROM THE MATRICES WE ALREADY SAVED")
    import numpy as np
    files = sorted(glob.glob("/kaggle/input/**/*_pae.npy", recursive=True))
    print("PAE matrices found: %d" % len(files))
    if not files:
        print("none attached - attach the AF2 screen's notebook output")
        return {}
    pat = re.compile(r"^(?P<name>.+)_(?P<state>g12d|wt|apo)_m(?P<model>\d+)_pae$")
    out = {}
    skipped = {}
    for f in files:
        stem = os.path.basename(f)[:-4]
        m = pat.match(stem)
        if not m:
            skipped["unparsable name"] = skipped.get("unparsable name", 0) + 1
            continue
        name = m.group("name")
        seq = seqs.get(name)
        if seq is None:
            skipped["no sequence for %s" % name] = 1
            continue
        a = np.load(f)
        if a.ndim != 2 or a.shape[0] != a.shape[1]:
            skipped["odd shape"] = skipped.get("odd shape", 0) + 1
            continue
        L = len(seq)
        n = a.shape[0]
        if n <= L + 20:
            skipped["matrix too small"] = skipped.get("matrix too small", 0) + 1
            continue
        # ColabDesign's binder protocol appends the binder after the target.
        ab = float(a[:n - L, n - L:].mean())
        ba = float(a[n - L:, :n - L].mean())
        rec = out.setdefault(name, {})
        rec.setdefault(m.group("state"), []).append(min(ab, ba))
    print("designs with PAE read: %d" % len(out))
    for k, v in sorted(skipped.items(), key=lambda kv: -kv[1])[:6]:
        print("    skipped, %-28s %d" % (k, v))

    res = {}
    for name, st in out.items():
        if not all(k in st for k in ("g12d", "wt", "apo")):
            continue
        g = sum(st["g12d"]) / len(st["g12d"])
        w = sum(st["wt"]) / len(st["wt"])
        p = sum(st["apo"]) / len(st["apo"])
        spread = (max(st["g12d"]) - min(st["g12d"])) if len(st["g12d"]) > 1 else 0.0
        res[name] = {"pae_g12d": round(g, 3), "pae_wt": round(w, 3),
                     "pae_apo": round(p, 3),
                     "pae_margin_wt": round(w - g, 3),    # positive = prefers G12D
                     "pae_margin_apo": round(p - g, 3),
                     "pae_noise": round(spread, 3)}
    print("designs with all three states: %d" % len(res))
    return res


def controls(res):
    pos = next((k for k in res if "CTRL_POS" in k), None)
    neg = next((k for k in res if "CTRL_NEG" in k), None)
    print()
    print("--- the controls, on pae_interaction this time ---")
    if not (pos and neg):
        print("    controls not present in this set - nothing can be concluded")
        return
    for k in (pos, neg):
        r = res[k]
        print("    %-24s pae(g12d) %6.2f   vs WT %+6.2f   vs apo %+6.2f   "
              "noise %5.2f"
              % (k[:24], r["pae_g12d"], r["pae_margin_wt"], r["pae_margin_apo"],
                 r["pae_noise"]))
    sep = res[neg]["pae_g12d"] - res[pos]["pae_g12d"]   # lower PAE is better
    noise = max(res[pos]["pae_noise"], res[neg]["pae_noise"])
    print()
    print("    shuffle minus real binder: %+.2f   (lower PAE is better, so a"
          " positive number here means the real binder scored better)" % sep)
    print("    worst within-state noise : %.2f" % noise)
    if sep <= noise:
        print()
        print("    STILL BLIND. pae_interaction separates a known binder of this")
        print("    target from a shuffled sequence by no more than the noise. Two")
        print("    metrics now agree that this screen cannot see binding, which is")
        print("    a statement about our pipeline - no AF2 initial guess - and not")
        print("    about the designs. Do not rank on it.")
    else:
        print()
        print("    IT SEPARATES. A known binder beats a shuffle by more than the")
        print("    noise, so pae_interaction is reading something that i-pTM was")
        print("    not. This can be used to rank, with the caveat that the")
        print("    threshold of 10 was calibrated with an AF2 initial guess and")
        print("    these predictions had none.")
        ok = [k for k, r in res.items()
              if not k.startswith("CTRL_") and r["pae_g12d"] < 12]
        print("    designs under the standard threshold of 12: %d of %d"
              % (len(ok), sum(1 for k in res if not k.startswith("CTRL_"))))


def main():
    seqs = load_sequences()
    if not seqs:
        raise SystemExit("no sequences found under /kaggle/input")
    pae = pae_interaction(seqs)
    if pae:
        controls(pae)
    pll = esm_pll(seqs)

    merged = {}
    for n in seqs:
        rec = {"name": n, "length": len(seqs[n])}
        rec.update(pll.get(n, {}))
        rec.update(pae.get(n, {}))
        merged[n] = rec
    out = os.path.join(W, "qc_metrics.json")
    json.dump(merged, open(out, "w"), indent=1)
    with open(os.path.join(W, "qc_metrics.csv"), "w", newline="") as fh:
        keys = ["name", "length", "pll", "pll_shuffled", "pll_excess",
                "pae_g12d", "pae_wt", "pae_apo", "pae_margin_wt",
                "pae_margin_apo", "pae_noise"]
        w = csv.writer(fh)
        w.writerow(keys)
        for n in sorted(merged):
            w.writerow([merged[n].get(k, "") for k in keys])
    banner("DONE")
    print("wrote qc_metrics.json and qc_metrics.csv")
    print("Neither metric is an affinity. Adaptyv's own guide says these work as")
    print("filters that eliminate failures, not as predictors of binding.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
