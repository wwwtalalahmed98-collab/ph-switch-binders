# Phase 2 - generation, and the arm comparison

Kaggle `notebook650f87c002` version 3, GPU T4 x2, 5 h 28 m. 28 designs, 0 failures.

## Result

| arm | hotspots | n | reached the mutation | mean on-mutation | mean on-peptide |
|---|---|---|---|---|---|
| A_tight | C5,C6,C7, num_pos 3 | 17 | **4 (24%)** | 0.24 | 2.06 |
| B_broad | C4..C9, num_pos 6 | 11 | 1 (9%) | 0.09 | 1.64 |

**Arm A is 2.6x better and the difference is not significant.** Fisher exact two-sided
**p = 0.62**. Reaching p < 0.05 at these rates needs about 60 designs per arm, roughly
19 hours of GPU - most of what remains.

**Decision: adopt arm A, stop testing the question.** It is directionally better and it
is the one phase 1 argued for on the exposure numbers. That is a judgement call on a
non-significant difference, and it is recorded as one rather than as a finding.

## The comparison is unbalanced by my own bug

Arm B logged its last design at its own 204-minute mark, then sat idle until 326. The
parent drains arm A's output stream to completion before touching arm B's, and
ColabDesign prints a line per optimisation step - roughly 60 kB over 11 designs, which
fills Linux's 64 kB pipe buffer. Arm B blocked on write for about two hours.

I had explicitly dismissed this risk when launching, on the grounds that arm B "produces
far too little output to fill a pipe buffer" - counting only my own log lines and
forgetting ColabDesign's. So arm B's 11 designs against arm A's 17 is my error, not a
property of the arms, and about two GPU-hours were lost.

Fix for the next run: drain both streams concurrently, or send each worker's output to
its own file and tail it.

## What the designs actually look like

Only 4 of 28 touch Asp6 at all, each with a single residue in contact. Most have pLDDT
near 0.40, which is not a folded protein. The two worth keeping:

```
A_tight_010  len 93 | ON-MUTATION 1 (7 atoms) | peptide 6 groove 9 | i_ptm 0.56 plddt 0.64
A_tight_013  len 81 | ON-MUTATION 1 (6 atoms) | peptide 6 groove 7 | i_ptm 0.32 plddt 0.64
```

One design, `B_broad_002`, makes **42 groove contacts and 4 peptide contacts** - the
exact failure the challenge is built to punish, a binder that reads the HLA and not the
peptide. It would bind all three states.

## The diagnosis from the loss, unchanged

Per-step lines read `loss 4.59 i_con 4.53`: the contact term is about 99% of the total
loss. It starts near 9.2 and plateaus around 4.5 without ever reaching satisfaction,
because the peptide is recessed and ringed by groove walls. So the optimiser spends the
whole budget chasing an unreachable contact with no pressure left to fold, which is why
pLDDT sits at 0.40.

The next run has to rebalance this. The weights will be printed before being overridden
rather than assumed.

## First thing in four challenges that worked first time

**PAE saved for 28 of 28 designs.** ipSAE has been impossible in every previous
challenge for one reason - we never saved the matrix it needs. It is computable now.

## Not done

- No design is a credible binder. 4 of 28 touch the mutation, with one residue each.
- Nothing has been through the three-state counter-screen.
- No ipSAE computed yet, though the inputs now exist.
- Nothing submitted. Deadline 8 Nov, final 7 Nov.
