# Challenge 2 (TNF-alpha) - Round 1 generation

Kaggle `tnf_generate` v2, GPU T4 x2, **ran 21742 s (6.04 h), completed normally**:
7.2 min setup + 5.5 h of design on two workers. **27 trajectories**, none failed.

## The construct was correct

All three setup assertions passed in the live run, which is the thing that cost a
round on EGFR:

```
cropped target: 122 residues (A 68, B 54)
site residues retained: 25 of 25
switch anchors present: K166(A), R108(B)
```

Cropping to a 22 A sphere was the decision that made this affordable. The smoke
test measured 43.3 s/iteration on the uncropped 304+70 construct; the cropped
122-residue target ran at **~26 min per complete design** (150 iterations plus a
3-recycle re-prediction), giving 27 designs where the uncropped target would have
given about six.

## Results

Every design is scored by an **independent re-prediction of the discrete sequence
at 3 recycles** (`true_iptm`), not by the optimiser's own number (`design_iptm`).
`on_site` counts how many of the 25 intended site residues the design actually
contacts within 5 A, mapped through `construct_map.json`.

- **9 of 27 at i_ptm >= 0.50**, best **0.758**
- **17 of 27 at i_ptm >= 0.35**
- **8 of 27 clear both i_ptm >= 0.50 and 5+ on-site contacts**

Full table: `round1_summary.csv`. The eight that clear both bars:

| name | len | design_iptm | **true_iptm** | pLDDT | epitope | on site |
|---|---|---|---|---|---|---|
| tnf_w1_012 | 73 | 0.587 | **0.758** | 0.830 | 29 | 11 |
| tnf_w0_003 | 79 | 0.506 | **0.717** | 0.852 | 24 | 9 |
| tnf_w0_012 | 73 | 0.556 | **0.691** | 0.718 | 26 | 11 |
| tnf_w0_008 | 63 | 0.368 | **0.682** | 0.751 | 27 | 10 |
| tnf_w0_011 | 60 | 0.582 | **0.638** | 0.775 | 24 | 6 |
| tnf_w1_004 | 80 | 0.391 | **0.632** | 0.746 | 10 | 7 |
| tnf_w0_004 | 70 | 0.408 | **0.597** | 0.745 | 24 | 13 |
| tnf_w0_007 | 79 | 0.453 | **0.504** | 0.668 | 16 | 12 |

## Three things this round establishes

**1. The epitope held.** This is the finding that matters most. On EGFR the designs
drifted off the intended patch and we only discovered it at the end. Here the top
designs put 9-13 of their contacts on the 25-residue site, and 14 of 27 designs
reach 5+ on-site contacts. Cropping the target to a sphere around the site is the
most plausible reason: there is much less off-target surface to drift onto.

**2. Re-scoring moved i_ptm UP, not down.** On EGFR round 1 the honest re-score
destroyed the result (0.725 -> 0.100). Here it improves almost every good design
(0.368 -> 0.682, 0.416 -> 0.647, 0.293 -> 0.489). The discipline is the same and
still necessary - it is what makes the number trustworthy in either direction -
but the direction is not a property of the method. The likely cause is the recycle
count: the optimiser runs at 0 recycles, the re-prediction at 3, and this target
benefits from recycling. The two numbers correlate only loosely (design_iptm 0.587
-> 0.758, but also 0.424 -> 0.384), so the optimiser's score remains unusable for
ranking.

**3. The two workers diverged sharply.** Worker 0 produced 8 of its 13 designs
above 0.48; worker 1 produced 11 of 14 below 0.20 - yet also the single best design
(tnf_w1_012, 0.758). Same settings, different RNG seed (2000+shard). This is
trajectory luck, not a settings difference, and it argues for more trajectories
rather than longer ones.

## Sequence liabilities

Much milder than EGFR's. The top four designs carry 0-1 cysteines and runs of at
most 6; only tnf_w0_007 (5 Cys) and tnf_w0_006 (7 Cys) are badly liable, and both
rank 9th-10th. The targeted repair from `kaggle/cellJ_repair.py` (cysteine -> Ser,
break runs >= 5, re-score and keep only what retains its i_ptm) should cost little
here.

## Not yet done

- acid-OFF switch placement (`tnf/scripts/acidoff_switch.py`) on the surviving designs
- liability repair + re-scoring
- novelty screen
