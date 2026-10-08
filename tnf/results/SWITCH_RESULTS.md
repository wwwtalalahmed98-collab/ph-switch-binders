# Challenge 2 (TNF-alpha) - acid-OFF switch placement and liability repair

Kaggle `tnf_switch`, GPU T4 x2, **ran 550 s (9.2 min), completed normally**. 13 designs
in, 26 records out (one switch and one matched control each), zero failures.

It was far quicker than the 40-55 min estimated: setup took 2.1 min rather than 7
because the AF2 parameters came from the attached input, and a 3-recycle prediction on
the 122-residue cropped target costs about 10 s, not the 60 s assumed from the
uncropped smoke test. Total cost 0.15 of the 4 GPU-hours remaining.

## The widened partner list mattered

```
construct 122 residues; 6 exposed cations (2 conserved anchors):
R120A, K166A*, K204A, R207A, R107B, R108B*
```

Phase 1 chose two conserved anchors. Applying the same exposure rule across the whole
crop finds six, four of them inside the intended site. **R107B - not an anchor - is
used by 7 of the 13 switches**, and `tnf_w0_011`, the second-best switch, places one of
its two histidines there. A whitelist restricted to the two anchors would have given
those designs fewer sites or none. This is the EGFR carboxylate lesson repeating
exactly: the narrow list was the limiting factor there too.

By contrast R120A and K204A, which sit outside the site, were never selected - the
distance and orientation filters excluded them on their own. Widening the list did not
mean accepting junk.

## Repair is free; the histidines are what cost

This is what the control arm was built to separate, and it answers cleanly.

**Repair alone is neutral.** Across 13 controls, i_ptm retention ranges 88-111% with a
median at 100%; five went *up*. Removing seven free cysteines from `tnf_w0_006` moved
it 0.497 -> 0.550. So the cysteine substitutions and run-breaking cost nothing
measurable, and any loss in a switch variant is attributable to the histidines rather
than to the repair. On EGFR this could not be said, because no control was run.

**The histidines cost a real but mostly acceptable amount.** Of 13 switches, 8 hold
>= 70% of parent i_ptm *and* stay above 0.45; 3 clearly break the interface
(`tnf_w1_008` 42%, `tnf_w0_009` 54%, `tnf_w0_007` 56%).

| design | His | anchors | partners | i_ptm | parent | kept | on site |
|---|---|---|---|---|---|---|---|
| tnf_w1_012_sw | 3 | 2 | K166A/R107B/R108B | **0.725** | 0.758 | 96% | 11 |
| tnf_w0_011_sw | 2 | 1 | R107B/R108B | **0.651** | 0.638 | 102% | 6 |
| tnf_w0_012_sw | 2 | 2 | R108B | **0.621** | 0.691 | 90% | 11 |
| tnf_w0_003_sw | 3 | 2 | R107B/R108B | **0.604** | 0.717 | 84% | 9 |
| tnf_w0_008_sw | 3 | 2 | R107B/R108B | **0.570** | 0.682 | 84% | 10 |
| tnf_w0_010_sw | 2 | 1 | R107B/R108B | **0.561** | 0.489 | 115% | 8 |
| tnf_w1_004_sw | 2 | 2 | R108B | **0.493** | 0.632 | 78% | 7 |
| tnf_w0_004_sw | 2 | 2 | R108B | **0.487** | 0.597 | 82% | 13 |
| tnf_w0_005_sw | 3 | 2 | R107B/R108B | 0.374 | 0.486 | 77% | 11 |
| tnf_w0_006_sw | 2 | 2 | R108B | 0.352 | 0.497 | 71% | 13 |
| tnf_w0_009_sw | 2 | 1 | R107B/R108B | 0.347 | 0.647 | 54% | 4 |
| tnf_w0_007_sw | 2 | 2 | R108B | 0.282 | 0.504 | 56% | 12 |
| tnf_w1_008_sw | 3 | 2 | R107B/R108B | 0.199 | 0.477 | 42% | 7 |

Two switches scored *above* their parent (`tnf_w0_011` 102%, `tnf_w0_010` 115%). That
is within prediction noise and should not be read as the histidines improving binding.

**Every one of the 26 sequences is free of cysteines**, and all residual homopolymer
runs are length 4, below the liability threshold.

## Compared with EGFR

On EGFR, ProteinMPNN redesign destroyed 12 of 14 interfaces, and placing histidines on
the one surviving MPNN sequence cost it everything (0.803 -> 0.091). Here, targeted
repair plus histidine placement retains a usable interface on 8 of 13. The difference
is not the switch engine - it is not doing wholesale redesign.

## Limitation worth stating

`on_site` and the epitope lists are measured on the **parent** structure, because the
switch variants were scored by prediction without saving their complexes. So we know
each switch retains its i_ptm, but have not confirmed it still binds the same epitope.
Nothing in the result rules out a histidine shifting the binding mode while keeping a
similar score. Re-predicting the 8 usable switches with structures saved would settle
it and costs about 5 minutes of GPU - worth doing before the deadline.

The pH response itself is a structural hypothesis, not a measurement: AlphaFold has no
concept of pH or protonation state, so i_ptm cannot report what happens at pH 6.0. The
claim rests on the calibration (His-cation close packing is 2.3x rarer in nature than
His-carboxylate, over 557 measured contacts) and on the FcRn/IgG control, where the
same engine recovered a known acid-ON switch in the right direction and did not fire
acid-OFF on it.

## Submission staged

`submission/challenge2_tnf_designs.csv`: **16 designs - the 8 usable switches and
their 8 matched controls**, format-validated. The pairs are the point: a switch that
binds at pH 7.4 and not at 6.0 is only interpretable if its histidine-free parent,
differing by 2-3 residues, does not show the same drop. Four slots of the 20 are left
free for a possible second generation round.
