# Challenge 2 (TNF-alpha) - confirmation of the submitted designs

Kaggle `tnf_confirm`, GPU T4 x2, **ran 421.8 s (7.0 min), completed normally**. The 8
usable switches and their 8 matched controls re-predicted with structures saved, so
that two things previously asserted could be measured.

## 1. i_ptm is reproducible

| | |
|---|---|
| max deviation from the switch stage | **0.026** |
| median deviation | **0.001** |
| within 0.005 | **7 of 8** |

Independent re-prediction of the same sequence returns the same number. This was not
guaranteed and it matters retroactively: every i_ptm quoted in this project is a
stable measurement rather than a draw from a noisy distribution. It also means the two
switches that scored *above* their parents (102%, 115%) are genuinely small real
differences in the model's opinion, not run-to-run scatter - though still far too
small to claim the histidines improved binding.

## 2. Seven of eight switches still bind the intended epitope

Measured on each switch's **own** predicted structure rather than its parent's.

| design | i_ptm | stage 3 | epitope kept | control kept | on site | His in reach |
|---|---|---|---|---|---|---|
| tnf_w1_012_sw | 0.725 | 0.725 | 93% | 97% | 14 | 3/3 |
| tnf_w0_011_sw | 0.651 | 0.651 | 92% | 100% | 5 | 2/2 |
| tnf_w0_012_sw | 0.622 | 0.621 | 88% | 96% | 12 | 2/2 |
| tnf_w0_003_sw | 0.600 | 0.604 | 96% | 96% | 10 | 3/3 |
| tnf_w0_010_sw | 0.561 | 0.561 | **68%** | 92% | **4** | 2/2 |
| tnf_w0_008_sw | 0.544 | 0.570 | 82% | 93% | 5 | **2/3** |
| tnf_w1_004_sw | 0.489 | 0.493 | 90% | 100% | 10 | 2/2 |
| tnf_w0_004_sw | 0.487 | 0.487 | 88% | 92% | 14 | 2/2 |

Four designs gained on-site contacts when measured on their own structure rather than
the parent's (`tnf_w1_012` 11 -> 14, `tnf_w1_004` 7 -> 10, `tnf_w0_004` 13 -> 14,
`tnf_w0_012` 11 -> 12), so the earlier parent-based figures were if anything
conservative.

**`tnf_w0_010_sw` fails and has been dropped.** It keeps only 68% of the parent
epitope and makes 4 on-site contacts, where its own control keeps 92% and makes 9. The
control is the reason this is attributable: the histidines moved this design off the
site, and nothing else did. It was already the weakest on independent grounds - lowest
pLDDT of the set at 0.533, and one of the two whose i_ptm rose above parent. This is
exactly the failure the confirmation run existed to catch, and without it the design
would have gone into the submission looking fine on i_ptm alone.

## 3. The designed histidines mostly survive re-folding

This is the measurement the acid-OFF hypothesis actually rests on: placement was
computed on the parent backbone, and re-folding with 2-3 residues changed could move
the histidines anywhere.

**7 of 8 designs have every designed histidine still within the 8.8 A window** of a
target cation. The single exception is `tnf_w0_008_sw`, whose H54 sits 10.2 A from
R108B - out of reach, so that histidine contributes nothing. Its other two
(H49 6.5 A, H46 5.1 A) are fine, so the design is kept with two working histidines
rather than three.

Many contacts are tight: H33 3.11 A, H27 4.01 A, H32 4.35 A, H28 4.41 A, H26 4.74 A,
H17 4.92 A. Against the calibration - natural His-cation closest approaches have a
median of 7.74 A and only 9.7% fall within 4.0 A - these sit in the region nature
avoids, which is precisely where protonation at pH 6.0 should cost energy.

## A caveat the run surfaced on its own

Of the 7 kept designs, **5 place a histidine on a conserved anchor** (K166A or R108B).
Three designs - `tnf_w0_011`, `tnf_w0_003`, `tnf_w0_010` - ended up with every
histidine facing **R107B instead**, which the widened whitelist allowed but which is
not one of the conserved anchors from Phase 1. Note that the switch stage reported
`n_anchor` from the parent-backbone placement; re-folding shifted which partner each
histidine is actually closest to, so those stage-3 anchor counts were optimistic.

This costs nothing for the competition, which assays human TNF-alpha. It only matters
for the mouse cross-reactivity argument used to pick anchors in Phase 1 - so it is a
narrowing of a robustness claim, not a defect in the designs.

## Submission updated

`submission/challenge2_tnf_designs.csv`: **14 designs - 7 confirmed switches and their
7 matched controls**, format-validated. Six of the 20 slots remain free.

Every switch in it now has: an i_ptm reproduced across two independent runs, a
confirmed epitope on its own predicted structure, 5+ contacts on the intended site,
and at least two histidines within reach of a target cation.
