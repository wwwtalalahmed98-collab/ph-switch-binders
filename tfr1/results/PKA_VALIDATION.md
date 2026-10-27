# The acid-release calculation, and what it can and cannot do

Built and validated 27 Oct 2026. Reproduce with:

    python tfr1/pka_switch.py --validate

## What it computes

Challenge 4 asks for a binder that holds human TfR1 at pH 7.4 and shows **no detectable
binding** at pH 6.0. Every pH ranking this project has used before was a geometric proxy, and
the first one we tried had the wrong sign (+0.24, see the TNF-alpha work). This replaces the
proxy with the quantity itself.

If a titratable group's pKa differs between the free and the bound state, binding and
protonation are thermodynamically coupled, and the size of that coupling *is* the pH
sensitivity. For each group i:

    L_i(pH) = log10( (1 + 10^(pKa_bound,i - pH)) / (1 + 10^(pKa_free,i - pH)) )
    log10( KD(6.0) / KD(7.4) ) = SUM_i [ L_i(7.4) - L_i(6.0) ]

pKa values come from PROPKA 3.5.1, run twice: once on the complex and once on the binder
alone.

**Sign convention.** Positive means weaker binding in acid, which is what this challenge
wants, and it requires pKa_bound < pKa_free - the histidine must be buried somewhere that
makes gaining a charge costly. Negative means tighter in acid, which is the FcRn mechanism and
what Challenge 1 (EGFR) asked for. Same machinery, opposite sign.

## Why validation was not optional

Published error on histidine pKa prediction is about one pH unit, and the whole effect being
engineered spans 1.4 pH units. A method with error that size can easily be useless and still
produce confident-looking numbers. So it was tested on complexes whose pH behaviour is known
before being pointed at any design.

| complex | His only | all groups | expected | verdict |
|---|---|---|---|---|
| **1I1A** FcRn + IgG Fc | **−0.94** | **−1.02** | strongly negative — binds at pH 6.0 | correct sign |
| **1BRS** barnase + barstar | +0.04 | −0.07 | near zero — no switch | correct |
| **3HFM** HyHEL-10 Fab + lysozyme | +0.20 | +0.21 | near zero — no switch | correct |

A separation of about one log unit between the known switch and the two non-switches. That is
enough to rank designs against each other.

## The first attempt failed, and the reason is worth keeping

The FcRn control was first run on **1FRT** and came out **+0.86** — the wrong sign, on the
textbook case. Rather than accept or discard the method, the structure was checked:

- His310 and His435, the two residues that actually drive the FcRn switch, were present and
  correctly numbered.
- PROPKA reported their pKa *dropping* on binding (6.81 → 2.20 and 5.11 → 3.47), the opposite
  of the known mechanism, in which they protonate at pH 6.0 and salt-bridge to FcRn
  carboxylates.
- Measuring the actual geometry: **His310's ND1 sits 1.2 angstrom from a Glu117 carboxylate
  oxygen.** No nitrogen-oxygen contact is shorter than about 2.6 angstrom. The side chains are
  clashing.
- **1FRT is a 4.50 angstrom structure.** Side-chain positions are not determined at that
  resolution.

So the failure was the input, not the equation. Re-run on **1I1A at 2.80 angstrom**, the same
calculation gives the right answer.

`pka_switch.py` now refuses any structure coarser than **3.0 angstrom** and says why, because
this failure mode is silent: a bad structure does not error, it returns a confident number
with the wrong sign.

## What it gets wrong

- **It underestimates magnitude by roughly one to two log units.** FcRn-IgG is reported in the
  literature as a hundredfold to thousandfold difference between pH 6.0 and 7.4 — two to three
  log units. This calculation gives one. So treat the output as **a ranking signal between
  designs scored the same way, not a predicted KD**.
- **One validation point for the positive case.** Two of the three controls are non-switches.
  A single known switch getting the right sign is encouraging, not proof.
- **PROPKA is empirical.** It ignores salt effects and, on the evidence of the 1FRT run, is
  sensitive to side-chain placement in a way that makes input quality decisive.
- **Designs are predicted structures.** Everything above used crystal structures. Our own
  designs will be scored on AlphaFold2 models, which have their own side-chain error, and
  nothing here measures how much that degrades the result. This is the largest untested gap.
- **The free state is approximated by deleting the target.** The binder is not relaxed after
  deletion, so any conformational change on binding is ignored.

## How it will be used

As a filter with a direction and a floor, not as a predicted affinity:

- require a **positive** predicted log ratio — the design must release in acid, not grip
  harder;
- prefer designs whose signal comes from **two or more separated histidines**, since the
  organizers' own stated physics caps one titratable site at tenfold per pH unit, and 1.4 pH
  units therefore allows at most about 25-fold from a single site, while "no detectable
  binding" needs around a hundredfold;
- report the per-residue contribution for every design so a reader can see whether the number
  rests on one histidine or several.
