# First evidence about the switch itself

Run: `tnf_switch_evidence` version 2, Kaggle, GPU T4 x2, 78m 42s, completed.
Script: [kaggle/cellTNF20_switch_evidence.py](../../kaggle/cellTNF20_switch_evidence.py).
Version 1 of this notebook measured nothing — see "The run that measured nothing" below.

This is the first time anything in this project has produced evidence about **objective 1**,
the acid-OFF switch. Everything before it measured binding at pH 7.4 only.

## Headline

> **Of 18 histidine positions in the submitted set, 2 punish a positive charge and 16
> tolerate it.**

And, more usefully:

> **Burial on binding predicts the charge effect (r = -0.58, n = 17). The
> Cbeta-to-cation distance we actually designed with does not (r = +0.24).**

Excluding one physically impossible contact (0.86 A, see defects) the correlations are
**-0.53 and +0.10** over 16 points. The conclusion survives; the margin is wide either way.

## What was measured

Four things on the 20 submitted designs. Three worked.

**1. Charge sensitivity.** AlphaFold2 knows nothing about pH, but it responds to charge.
At each histidine, the position was predicted as Gln (neutral, similar size) and as Lys
(charged, similar length). The difference isolates charge from sterics. A strongly
negative `Lys - Gln` means the position rejects a positive charge, which is the only
reason a protonated histidine would push the binder off at pH 6.0.

| | Lys - Gln (human i_ptm) |
|---|---|
| best | **-0.198** (sw_04 H37, `tnf2_w0_001` Y37H) |
| second | **-0.105** (sw_06 H18, `tnf3_w0_017` Y18H) |
| median | -0.051 |
| weakest | -0.003 |
| pass at <= -0.10 | **2 of 18** |

Every value is negative, so a positive charge always costs *something*. Almost none of
it is large.

**2. Burial on binding.** Histidine side-chain relative SASA in the free binder minus the
same in the complex. Median **0.143**, range 0 to 0.355. **Three positions are buried by
0.00** — `sw_01` H57, `sw_10` H55, `sw_12` H45 are as exposed when bound as when free, so
protonating them costs nothing and there is no mechanism for them to switch at all.

**3. Error bars.** Only **2 of the 5** requested parameter sets loaded, so the spread is
over two models, not five. Human i_ptm sd: median 0.021, max 0.074. Mouse: median 0.023,
max **0.305**.

**4. Fold without the target.** Failed on all 20 — see defects. **Still unknown.**

## Reproducibility, measured by accident

Every design was re-predicted from the same sequence that `switch5` had scored. The two
runs should agree exactly and do not:

| | median absolute difference | worst |
|---|---|---|
| human i_ptm | 0.021 | 0.074 |
| mouse i_ptm | 0.026 | **0.305** |

A mouse number that moves by 0.305 on re-prediction is not a measurement we should be
ranking on. Several margins in the submitted set are smaller than this.

## What this means

**The switch design, as built, has little in-silico support.** We placed histidines by a
geometric rule and 16 of 18 of those positions do not noticeably object to a positive
charge sitting there.

**But this is absence of evidence, not evidence of absence.** AlphaFold2 is a poor
predictor of free-energy changes and is known to respond weakly to point mutations. A null
result here cannot cleanly separate "the switch will not fire" from "AF2 cannot see it".
Lysine is also an imperfect stand-in for imidazolium — longer, flexible, differently
shaped. The honest statement is that **we looked for support for objective 1 and did not
find it at 16 of 18 positions.**

**The part that is actionable is the correlation.** This is the same shape of finding as
the mouse one: a cheap geometric proxy carrying almost no information, beaten by a direct
measurement.

| variable | correlation with charge effect | was it used to design? |
|---|---|---|
| burial on binding | **-0.58** | no |
| screened Coulomb | -0.28 | no |
| Cbeta-to-cation distance | +0.24 | **yes, this was the rule** |

The 3.0-8.8 A distance window we filtered on is roughly as uninformative for the switch
as epitope conservation was for mouse binding. Burial on binding is the better variable
and it is cheap — it needs no GPU, only the complex structure we already save.

**Two designs have no switch mechanism at all.** `sw_10` and `sw_12` have zero burial and
among the weakest charge effects (-0.028, -0.021). Their histidines sit in open solvent
when bound. `sw_01` is half-dead: H54 is buried 0.204 and contributes, H57 is buried 0.000
and does not.

## Defects found in this run

- **The monomer prediction failed on all 20 designs.** `protocol="hallucination"` needs its
  own `prep_inputs` call, which the cell does not make, so every attempt raised
  `AssertionError: the following inputs were not set: {'length': ...}`. **Whether our
  binders fold without the target is still unknown** — the question the Proteinbase pLDDT
  numbers raised remains open. Step isolation worked: this cost only this measurement.
- **Only 2 of 5 parameter sets loaded** despite requesting `num_models=5`, across all 152
  model builds. The error bars are therefore over two models and understate the true
  spread.
- **One physically impossible contact.** `sw_06` H18 reports a 0.86 A nitrogen-to-nitrogen
  distance, which is a clash, not a contact — the predicted structure is wrong there. It
  produced the largest screened-Coulomb value in the set (11.5 kcal/mol), which should be
  discarded. All correlations above are reported with and without this point.
- **The cation identity labels may be in the wrong numbering system.** The burial code
  reads residue numbers straight from the predicted PDB rather than mapping through the
  construct map, so a label like `A199NH1` may not be UniProt 199. The distances and
  burial fractions are geometric and unaffected; only the names are suspect. Not yet
  fixed.

## The run that measured nothing

Version 1 reported **"Successfully ran"** and produced 20 rows of `0.000`. Every design
died with `IndexError` before any science happened.

The cause was a wrong assumption about the library: ColabDesign lazily loads only
`num_models` parameter sets, default 1, so `predict(models=[0])` works and `models=[1]`
raises. Each design failed about 9 s into prediction — exactly one prediction's worth.

The cause that mattered was structural and mine. The **optional** measurement, error bars,
ran first inside a single try block covering everything, so when it threw it took down the
charge test, which was the entire point of the run. The fix was to order the measurements
most-important-first and isolate each one. That is why version 2 still returned the charge
test and burial data despite the monomer step failing on all 20 designs.

The cell now refuses to exit quietly when nothing was measured.

Cost of version 1: 11 minutes of the 30-hour weekly quota.

## What this does not change

The submitted set is unaffected and remains banked and designated. Nothing here says the
designs do not bind; binding at pH 7.4 is the one thing we have consistently measured.
What it says is that our **reason for believing they will release at pH 6.0 is weak**, and
that we now have a better variable for choosing placements than the one we used.

## What it suggests doing next

Measure burial and charge sensitivity for **all 35 verified placements**, not only the 12
that were chosen, then re-rank on burial. Several lineages have alternative placements we
never tested, and the two designs with zero burial should be replaced by placements that
are actually buried on binding. This needs no new design work and costs roughly 2 hours of
GPU.
