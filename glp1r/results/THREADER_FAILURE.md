# The threading job failed, and confidence pointed the wrong way

Measured 24 Oct 2026 from the `glp1r_generate` run, 11 threading designs complete.

## Result

**One of eleven threaders is usable.** The other ten reach into space that exists only
because we cropped the bottom off the receptor to make the target affordable.

| # | i_ptm | mean KD | on-site | deepest | N-term | verdict |
|---|---|---|---|---|---|---|
| 001 | 0.800 | +0.25 | 2 | 47.4 Å | 46.6 | past crop, off-site |
| **002** | **0.442** | **−0.28** | **27** | **34.7 Å** | **33.0** | **usable** |
| 003 | 0.720 | +1.00 | 0 | 43.6 Å | 3.5 | backwards, greasy |
| 004 | 0.757 | −0.23 | 0 | 44.2 Å | 40.9 | past crop, off-site |
| 005 | 0.740 | +0.79 | 0 | 47.7 Å | 45.5 | past crop, greasy |
| 006 | 0.481 | −0.04 | 25 | 56.8 Å | 56.8 | 10 Å beyond the target entirely |
| 007 | 0.740 | +0.38 | 0 | 46.1 Å | 41.0 | backwards |
| 008 | 0.441 | −0.14 | 25 | 65.2 Å | 7.0 | backwards, far past target |
| 009 | 0.664 | −0.21 | 2 | 46.3 Å | 44.8 | past crop, off-site |
| 010 | 0.672 | +0.20 | 3 | 48.8 Å | 45.6 | past crop, off-site |
| 011 | 0.757 | −0.10 | 1 | 44.3 Å | 44.3 | past crop, off-site |

"On-site" counts how many of the 30 receptor residues the natural hormone contacts inside
the bundle are also contacted by the design. The native hormone's deepest point is 34.2 Å;
our crop ends at 46.8 Å.

## The finding worth keeping

**Confidence anti-correlates with binding in the right place: r = −0.963 over 11 designs.**

Every design scoring above 0.66 touches at most 3 of the hormone's own contact residues.
The only two designs that land on the hormone's site (002 with 27, 006 with 25) score
0.442 and 0.481 — the two lowest confidences in the set.

This is not AlphaFold being unreliable in general. It is doing exactly what it was asked:
find a sequence that packs confidently against the target it was given. We gave it a
target with a hole in the bottom, and a hole is an easy place to pack confidently. The
optimiser found the flaw in our preparation faster and more consistently than it found the
binding mode we wanted.

**So for this job, selecting on confidence would have selected almost perfectly wrongly.**
Ranking the 11 by i_ptm and taking the top 6 would have chosen six designs with 0, 0, 0,
1, 2 and 2 on-site contacts, and rejected the only one that works.

## What caused it

Depth bands, comparing our cropped target with the intact receptor from 6X18:

| depth | atoms in our crop | atoms in the intact receptor | kept |
|---|---|---|---|
| 33.2–40.0 Å | 259 | 284 | 91% |
| 40.0–45.0 Å | 107 | 195 | 55% |
| 45.0–50.0 Å | 4 | 211 | 2% |
| beyond 50 Å | 0 | 825 | 0% |

Below about 45 Å our target is essentially absent. In the intact receptor that region is
solid helix, so all ten rejected designs would clash into protein that our target did not
contain.

## What we are doing about it

Filtering, and saying so. The threader rule requires the deepest atom within 40 Å, where
the crop still represents at least half the real receptor, and at least 5 contacts with
the hormone's own site.

Re-cropping deeper and re-running was considered and rejected on time: the window closes
25 Oct and the landing-pad half of the run — the half that can compete on binding
affinity — is working well.

## What should be done differently next time

Seal the target before designing into it. A crop that opens a cavity is not a neutral
simplification; it is an invitation the optimiser will accept. Either keep the full bundle,
or cap the opening with a few residues of real structure, or add a repulsive term below the
depth the natural ligand reaches. Validating the *target* deserves the same effort we spent
validating the selection metric.
