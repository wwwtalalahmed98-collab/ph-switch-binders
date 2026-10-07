# Round 2 generation - complete

Notebook `round2` Version 2, GPU T4 x2, ran 21989 s (6.1 h) and finished cleanly.
24 trajectories across two workers (13 from worker 0, 11 from worker 1), no failures. Every i_ptm below is an independent
re-prediction of the discrete sequence at 3 recycles, not the optimiser's own number.

**13 designs at i_ptm >= 0.50, 18 at >= 0.35, best 0.786.**
Round 1's best, scored honestly, was 0.273 - it would rank 15th here.

| design | len | i_ptm | (design said) | pLDDT |
|---|---|---|---|---|
| r2_w0_002 | 55 | 0.786 | 0.635 | 0.705 |
| r2_w0_013 | 63 | 0.748 | 0.722 | 0.634 |
| r2_w0_007 | 57 | 0.718 | 0.622 | 0.624 |
| r2_w0_006 | 66 | 0.701 | 0.604 | 0.545 |
| r2_w0_009 | 65 | 0.695 | 0.583 | 0.547 |
| r2_w0_010 | 56 | 0.684 | 0.595 | 0.622 |
| r2_w0_003 | 64 | 0.678 | 0.590 | 0.717 |
| r2_w1_004 | 69 | 0.592 | 0.709 | 0.628 |
| r2_w1_001 | 53 | 0.580 | 0.516 | 0.642 |
| r2_w0_005 | 54 | 0.560 | 0.610 | 0.430 |
| r2_w1_009 | 70 | 0.560 | 0.620 | 0.636 |
| r2_w1_006 | 57 | 0.557 | 0.357 | 0.630 |
| r2_w0_001 | 65 | 0.501 | 0.443 | 0.477 |
| r2_w0_004 | 63 | 0.498 | 0.100 | 0.483 |
| r2_w0_012 | 59 | 0.407 | 0.505 | 0.552 |
| r2_w1_005 | 64 | 0.407 | 0.446 | 0.596 |
| r2_w1_002 | 58 | 0.372 | 0.186 | 0.537 |
| r2_w1_007 | 65 | 0.350 | 0.491 | 0.675 |
| r2_w1_008 | 63 | 0.343 | 0.524 | 0.454 |
| r2_w0_008 | 69 | 0.292 | 0.432 | 0.590 |
| r2_w1_003 | 54 | 0.250 | 0.600 | 0.508 |
| r2_w0_011 | 59 | 0.243 | 0.541 | 0.482 |
| r2_w1_011 | 64 | 0.152 | 0.173 | 0.753 |
| r2_w1_010 | 72 | 0.126 | 0.100 | 0.810 |

## What the two columns say about the round-1 failure

The design-time and re-scored values disagree in both directions - r2_w0_004 went
0.100 -> 0.498, r2_w1_003 went 0.600 -> 0.250. The optimiser's own number is simply
not a measurement, which is the single most important thing round 1 got wrong.

The improvement over round 1 cannot be attributed to iteration depth alone. Round 1's
re-score ran at 0 recycles and this one at 3, and recycling raises i_ptm by itself
(r2_w0_004's jump is a recycle effect, not a depth effect). The two changes are
confounded. What is defensible is the conclusion, not the mechanism: these designs
pass a standard filter and round 1's did not.

## Weakness to carry into the submission

pLDDT and i_ptm diverge badly at both ends. r2_w1_010 has the highest pLDDT in the
set (0.810) and the lowest i_ptm (0.126); r2_w0_005 has the second-lowest pLDDT
(0.430) with a solid 0.560. ColabDesign's binder defaults are `i_con 1.0, plddt 0.1`
- the objective barely penalises a poorly folded binder. A design AlphaFold cannot
fold confidently is an expression risk in the wet lab however good its interface
looks, so pLDDT stays in the ranking as a tiebreak and the low-pLDDT designs are a
real caveat rather than something to present as a clean result.
