# Burial and charge for every verified placement

Run: `tnf_burial_all` version 1, Kaggle, GPU T4 x2, **40m 28s**, completed — 35 placements
across 14 lineages. Script:
[kaggle/cellTNF21_burial_all_placements.py](../../kaggle/cellTNF21_burial_all_placements.py).

Half the wall-clock of the 20-design evidence run for nearly twice the work, because the
model is no longer rebuilt between every prediction.

## Three findings, in order of how much they matter

### 1. Only one placement in 35 shows a charge effect above noise

| | Lys - Gln, worst histidine |
|---|---|
| strongest | **-0.264** (`tnf2_w0_001` Y37H) |
| median | -0.019 |
| weakest | **+0.083** (a positive charge *helps*) |
| below -0.10 | **1 of 35** |
| actually positive | **8 of 35** |

`tnf2_w0_001` Y37H is the only position clearly above the noise floor, and it is the only
one that **reproduced**: the evidence run measured -0.198 for the same position in a
different context. Everything else sits within the ±0.02–0.03 spread we measured between
identical re-predictions.

This is the clearest statement we can make about objective 1: **across every histidine
placement this project produced, one has in-silico support and 34 do not.**

The caveat from the evidence run stands and is important — AlphaFold2 responds weakly to
point mutations and lysine is an imperfect stand-in for imidazolium, so this cannot
separate "the switch will not fire" from "AF2 cannot see it". It is absence of evidence.
But we now have that absence across 35 placements rather than 18, and the one positive is
consistent across two runs.

### 2. I was wrong about target-dependent folding

The monomer prediction worked this time. pLDDT **bound minus alone**, over 35 placements:

> median **-0.076**, range -0.144 to +0.074

The sign is negative, which means the binders score **higher on their own than in
complex**. They are not target-dependent folders. The hypothesis I drew from Proteinbase's
pLDDT of 55–63 against our 60–77 was wrong, and the difference must come from their using
a different predictor or settings rather than from predicting the monomer.

Recorded as a reversal. The reason it is a reversal and not a mistake we acted on is that
the inference was flagged as unverified when it was made and measured rather than trusted.

### 3. Every two-histidine placement has a dead partner

Ranking all 35 by the burial of their **worst** histidine, singles occupy the whole top of
the table and pairs the bottom:

| | best worst-histidine burial |
|---|---|
| best single | **0.404** (`tnf_w0_012` A28H) |
| best pair | **0.100** (`tnf3_w1_009` I22H,L18H) |
| submitted sw_01 pair (`tnf3_w0_012` E54H,W57H) | **0.000** |

**7 of 35 placements have a histidine buried by less than 0.02** — it sits in open solvent
when bound, so protonating it costs nothing and it cannot contribute to release.

This matters because the whole argument for two histidines is Henderson-Hasselbalch: one
histidine leaves 50% of molecules uncharged and still binding at pH 6.0, two leaves 25%.
**That argument assumes both histidines are doing work.** If one of the pair has zero
burial, the effective count is one — and you have paid the pH 7.4 binding cost of two to
get the pH 6.0 push of roughly one. On this evidence the pairs we submitted are worse than
the best single from the same lineage, not better.

## Best placement per lineage, by burial

| lineage | histidines | burial | charge | human | mouse | currently submitted? |
|---|---|---|---|---|---|---|
| tnf_w0_012 | A28H | 0.404 | -0.002 | 0.710 | 0.638 | no — we sent the pair |
| tnf2_w1_001 | L13H | 0.383 | -0.034 | 0.576 | 0.576 | no — we sent the pair |
| tnf3_w0_012 | L39H | 0.383 | -0.035 | 0.694 | 0.608 | no — we sent E54H,W57H |
| tnf3_w0_017 | Y18H | 0.382 | -0.036 | 0.606 | 0.517 | no — we sent the pair |
| tnf_w0_011 | Q19H | 0.327 | +0.024 | 0.651 | 0.638 | no — we sent L18H |
| tnf2_w0_001 | **Y37H** | 0.271 | **-0.264** | 0.498 | 0.630 | no — we sent the pair |
| tnf2_w0_007 | L30H | 0.266 | -0.065 | 0.585 | 0.656 | **yes** |
| tnf3_w1_009 | I22H | 0.263 | -0.045 | 0.559 | 0.587 | no — we sent the pair |
| tnf3_w0_016 | L9H | 0.184 | -0.032 | 0.476 | 0.452 | no — not submitted |
| tnf3_w0_007 | E56H | 0.180 | +0.018 | 0.491 | 0.516 | no — we sent E55H |
| tnf2_w0_008 | E2H | 0.154 | -0.015 | 0.598 | 0.628 | **yes** |
| tnf_w0_007 | Q46H | 0.059 | -0.019 | 0.486 | 0.502 | yes, but weak |
| tnf3_w1_010 | A30H | 0.005 | +0.023 | 0.453 | 0.519 | no — effectively dead |
| tnf_w1_004 | Q45H | 0.005 | -0.022 | 0.663 | 0.477 | yes — effectively dead |

Only **2 of 12** submitted placements are the best their lineage has to offer.

## Two smaller results

**30 of 41 histidines face one of the two conserved anchors** (K166 chain A, R108 chain B).
The remaining 11 face some other cation, which is not guaranteed to exist in mouse — so
for those the switch, if it works at all, may be human-only. This is the first time the
cation identities have been labelled in the right numbering system; the evidence run read
residue numbers straight out of the predicted PDB and got the names wrong.

**One histidine still sits 0.86 A from a cation**, which is a clash, not a contact. The
predicted structure for `tnf3_w0_017` Y18H,E19H is wrong there, and any energy computed
from that geometry should be discarded.

## What this does not say

- It does not say any design binds or fails to bind. Binding at pH 7.4 is the one thing
  we have measured consistently, and it is unchanged.
- It does not say the switch will not work. It says we looked for in-silico support at 35
  placements and found it at one.
- Burial is a *correlate* of the charge effect (r = -0.58, n = 17), not a measurement of
  switching. The chain from burial to a real pH response has two unvalidated links in it.
- Nothing here has been physically tested.
