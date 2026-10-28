# Generation run 1: the designs bind the wrong surface

Measured 2026-10-28 by `kaggle/cellTFR1d_remeasure.py` (notebook `tfr1_remeasure`, version 2,
CPU, 43 s). Both crops verified first: e559 108 residues / 13 epitope, l537 139 residues /
14 epitope, kept list agrees with the crop PDB in both cases. 47 of 47 designs measured, no
refusals.

## The result

| | out of 47 |
|---|---|
| zero contacts with the chosen epitope | **32** |
| at least 1 epitope contact | 15 |
| at least 2 | 15 |
| at least 4 | 3 |
| at least 6 | 1 |

## The part that decides the next step

The five most confident designs all make **zero** epitope contacts:

| design | i_ptm | on-site | off-site |
|---|---|---|---|
| tfr1_l537_014 | 0.782 | 0 | 31 |
| tfr1_e559_001 | 0.740 | 0 | 18 |
| tfr1_l537_002 | 0.735 | 0 | 13 |
| tfr1_e559_006 | 0.697 | 0 | 21 |
| tfr1_l537_016 | 0.664 | 0 | 15 |

They are confident, well-formed, and stuck to the face the crop exposed. This is precisely
the failure the on-site/off-site measurement was written to catch, and it is the majority
outcome, not an edge case. `cellTFR1a_setup.py` checked that the crop did not distort the
epitope (it did not: +5 A^2) but nothing stopped the 3,427 A^2 of new surface - 1,722 A^2 of
it hydrophobic - from being a more attractive place to bind than the epitope itself.

The hotspot argument to `prep_inputs` biases the optimiser; it does not constrain it. On this
target the bias lost.

## The designs that did hit the epitope are not usable either

| design | on/epitope | off | i_ptm | pLDDT | longest run | note |
|---|---|---|---|---|---|---|
| tfr1_e559_011 | 7/13 | 14 | 0.456 | 0.742 | 8 | `SPPSPPPPPPPPAAAAAAT...` |
| tfr1_e559_010 | 5/13 | 12 | 0.246 | 0.678 | 9 | `...NNNNNNNNNSYNNTNN...` |
| tfr1_e559_026 | 4/13 | 7 | 0.149 | 0.615 | 7 | `SAAAIILLLLLGGGGGGGPG...` |

The designs that reach the epitope are the ones the optimiser gave up on - low confidence,
homopolymer runs, hydrophobic. The designs that look like proteins bind the wrong place.
There is no design in this set that is both.

## Confidence does not predict epitope accuracy here

Over the 20 designs with at least one epitope contact, i_ptm versus on-site contacts gives

    r = +0.103

Effectively zero. On GLP-1R the same comparison gave -0.96 (n=11) and -0.667 (n=15), which is
why it was measured again rather than assumed.

**The GLP-1R anti-correlation does not transfer.** That finding was about a buried channel,
where high confidence meant the binder had found a pocket that was not the real one. On this
surface epitope, confidence carries no information about where the binder sits, in either
direction. It cannot be used to screen, and neither can its inverse.

## What this costs

GPU quota at the time of writing: **23:10 used of 30 this week, 6 h 50 m left.** Generation
run 1 spent 5 h 48 m of it and produced nothing selectable. Three things still want GPU time:
a second generation attempt, the histidine scan, and the mouse cross-reactivity check.

## What was learned, as distinct from what was spent

- The measurement pipeline works end to end and is verified against a synthetic case.
- The crop's cut face is a real competitor for the binding site, quantitatively, not in
  principle. Any second attempt has to deal with it rather than hope.
- Confidence screening is unavailable on this target, in both directions.
