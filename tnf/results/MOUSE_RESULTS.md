# Mouse cross-reactivity, measured

Source: `tnf_mouse2` on Kaggle, GPU T4 x2, 23.7 min. 68 designs - every switch and
control from the cap-2 stages of all three rounds - predicted against a mouse TNF-alpha
target built from PDB 2TNF, same model, same 3 recycles, same scoring path as the human
numbers. Only the target changes.

Objective 2 is mouse cross-reactivity and it is ranked above affinity. Until this run it
was being inferred from sequence conservation across each design's human epitope.

## The numbering fix held

    chain A: 148 residues mapped to P06804, 0 mismatches, offsets present [78, 79]
    chain B: 148 residues mapped to P06804, 0 mismatches, offsets present [78, 79]
    anchors confirmed in 2TNF: K168 = resi A90, R111 = resi B32
    mouse site residues retained: 25 of 25

`offsets present [78, 79]` is the thing that killed the first attempt: the offset from
UniProt numbering to 2TNF author numbering genuinely takes two different values inside a
single chain, so no constant could have been right. The site sequences differ at 6 of 25
positions:

    mouse  NLSQRANALTFISYQEKVNQEQFAS
    human  NLNRRANALTIVSYQTKVNQERFAS

## The conservation proxy did not work

Across the 11 designs that had both a proxy score and a measured value, the correlation
between epitope conservation and measured mouse i_ptm is **r = +0.272**. That is close
to no information. The design the proxy ranked first, `tnf2_w0_007_sw` at 82.6%
conserved, has the worst measured mouse binding in the whole set: i_ptm 0.148, 30% of
its human value, and **zero** contacts on the mouse site.

| design | proxy cons% | measured mouse | human | retained |
|---|---|---|---|---|
| tnf2_w0_007_sw | 82.6 | 0.148 | 0.487 | 30% |
| tnf2_w0_008_sw | 82.4 | 0.628 | 0.598 | 105% |
| tnf_w0_012_sw | 82.1 | 0.592 | 0.620 | 96% |
| tnf2_w0_001_sw | 80.8 | 0.654 | 0.470 | 139% |
| tnf_w1_012_sw | 80.6 | 0.326 | 0.744 | 44% |
| tnf_w0_004_sw | 79.2 | 0.189 | 0.487 | 39% |
| tnf_w0_008_sw | 79.2 | 0.284 | 0.710 | 40% |
| tnf2_w1_001_sw | 78.6 | 0.566 | 0.553 | 102% |
| tnf_w0_011_sw | 75.0 | 0.644 | 0.651 | 99% |
| tnf_w0_003_sw | 73.1 | 0.351 | 0.650 | 54% |
| tnf_w1_004_sw | 61.5 | 0.222 | 0.493 | 45% |

Conservation counts residues that are shared. It says nothing about whether the six that
differ are the ones the binder leans on, and that is what decides the outcome.

## Six of the ten submitted switches fail objective 2

| submitted | internal | human | mouse | retained | mouse on-site | verdict |
|---|---|---|---|---|---|---|
| tnf_sw_01 | tnf_w0_012_sw  | 0.620 | 0.592 |  96% | 12 | holds |
| tnf_sw_07 | tnf2_w1_001_sw | 0.553 | 0.566 | 102% |  9 | holds |
| tnf_sw_09 | tnf2_w0_008_sw | 0.598 | 0.628 | 105% | 10 | holds |
| tnf_sw_10 | tnf2_w0_001_sw | 0.470 | 0.654 | 139% | 12 | holds |
| tnf_sw_02 | tnf_w1_012_sw  | 0.744 | 0.326 |  44% |  **0** | fails |
| tnf_sw_04 | tnf_w0_003_sw  | 0.650 | 0.351 |  54% |  7 | fails |
| tnf_sw_06 | tnf_w0_008_sw  | 0.710 | 0.284 |  40% |  9 | fails |
| tnf_sw_08 | tnf_w1_004_sw  | 0.493 | 0.222 |  45% |  7 | fails |
| tnf_sw_03 | tnf_w0_004_sw  | 0.487 | 0.189 |  39% |  6 | fails |
| tnf_sw_05 | tnf2_w0_007_sw | 0.487 | 0.148 |  30% |  **0** | fails |

`tnf_sw_02` was the highest human affinity in the submission at 0.744. It does not touch
the mouse site at all.

## What clears both bars

Switches with human i_ptm >= 0.50, mouse i_ptm >= 0.45 and 5+ mouse on-site contacts:

| design | human | mouse | retained | mouse on-site | submitted? |
|---|---|---|---|---|---|
| tnf3_w0_012_sw | 0.722 | 0.669 |  93% | 16 | no |
| tnf2_w0_008_sw | 0.598 | 0.628 | 105% | 10 | yes |
| tnf_w0_012_sw  | 0.620 | 0.592 |  96% | 12 | yes |
| tnf2_w1_001_sw | 0.553 | 0.566 | 102% |  9 | yes |
| tnf3_w0_017_sw | 0.561 | 0.499 |  89% |  7 | no |
| tnf_w0_011_sw  | 0.651 | 0.644 |  99% |  4 | no (4 on-site, just under) |
| tnf2_w0_001_sw | 0.470 | 0.654 | 139% | 12 | yes (human just under 0.50) |

`tnf3_w0_012_sw` is the best design the project has produced on every axis at once:
highest confirmed human i_ptm among switches (0.722), 93% retained against mouse, 16
contacts on the mouse site, both histidines in the geometric window facing an arginine
conserved between the species.

## Two things this does not establish

- **i_ptm is not affinity.** A design that scores well against both targets has been
  shown to be predicted to fold against both, nothing more. No pH response has been
  measured anywhere in this project; the switch is a geometric hypothesis.
- **`tnf3_w0_017_sw` beat my reasoning, and I do not know why.** Its histidines point at
  R107, which is glutamine in mouse, and it has only 3 on-site contacts against human.
  It was expected to fail objective 2 and it measured 0.499 with 7 mouse-site contacts.
  Being right by accident is not the same as being right.
