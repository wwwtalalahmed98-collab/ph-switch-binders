# Phase 1 - the hotspots, and how small the target really is

Run 2026-11-04, CPU only, no quota spent. Script: [phase1_hotspots.py](../phase1_hotspots.py).

## The number that defines this challenge

Peptide exposure in the complex, with the TCR removed - what a binder actually sees:

| position | residue | exposed A^2 | of maximum |
|---|---|---|---|
| 1 | Val | 1.2 | 1% |
| 2 | Val | 2.1 | 1% |
| 3 | Val | 5.9 | 3% |
| 4 | Gly | 37.0 | 36% |
| 5 | Ala | 26.9 | 21% |
| **6** | **Asp** | **21.8** | **11%** |
| 7 | Gly | 75.5 | 73% |
| 8 | Val | 46.6 | 27% |
| 9 | Gly | 15.7 | 15% |
| 10 | Lys | 6.7 | 3% |

**Asp6 - the only thing that distinguishes the target from the counter-target - is 11%
exposed.** Positions 1, 2, 3 and 10 are the anchors holding the peptide in the groove and
are effectively buried, as expected.

That single number is the difficulty of this challenge stated plainly. There is about
22 A^2 of mutation to read, sitting beside 75 A^2 of Gly7 that is identical in the wild
type. Everything downstream has to be built around that imbalance.

## It is readable, though

7OW6 ships with a real T-cell receptor bound to this exact complex, doing the job we are
being asked to do. Its footprint at 4.5 A:

- peptide: **G4, A5, D6, G7, V8, G9** - it does touch the mutation
- heavy chain: E58, Q62, R65, N66, K68, A69, Q72, T73, R75, V76, A150, H151, E154, Q155, A158, R163

So an 11%-exposed residue is enough for a natural binder. The premise holds.

**Only positions were taken from the receptor, never sequence.** The challenge requires
de novo designs and bans modifying an existing binder. Contact geometry is public
structural information about the target, the same as the coordinates; the receptor's
sequence appears in no output of this script and goes nowhere near the design stage.

## The hotspots, and the trap avoided

The obvious set is the six peptide positions the receptor reads, 4 through 9. **That set
is a trap.** The contact term is averaged over its positions, and Gly7 at 73% exposure is
six times easier to reach than Asp6 at 11%. A design would satisfy the term by piling
onto Gly7 and never touching the mutation - and a binder that reads only Gly7 reads
nothing the wild type lacks. It would bind both and score zero on the first two ranking
steps.

**Chosen: peptide positions 5, 6, 7 with `num_pos=3`, `i_con` cutoff 8.0 A.**

Every member is adjacent to the mutation, so a design satisfying all three sits over it
rather than beside it. `num_pos=3` of 3 means none can be quietly dropped. This is the
TfR1 round-1 lesson - a contact term spread over 13 positions with a 21.7 A cutoff
pointed nowhere - applied before the run instead of diagnosed after it.

Groove residues Q155, T73, A69, N66 (3.8-6.9 A from the mutation, all in the receptor
footprint) are recorded for register but **are not hotspots**. Adding them would reopen
the same escape route one level up.

## The aim metric for phase 2

On TfR1 two whole rounds were measured with a broken ruler. Defining it now, before
generating anything:

- `on_mutation` - binder heavy atoms within 4.5 A of **peptide residue 6**
- `on_peptide` - within 4.5 A of any peptide residue
- `on_groove` - within 4.5 A of the heavy chain
- `mutation_share` = `on_mutation` / `on_peptide`

**`on_mutation` is the number that matters, not `on_peptide`.** A design with 20 peptide
contacts and zero on the mutation is a failure that will look like a success in any
summary that reports peptide contacts alone.

The target chain is identified by reproducing its sequence residue for residue, never by
residue numbering, because `save_pdb` renumbers from 1.

## Honest limits

- A hotspot term can push a design toward Asp6. It **cannot guarantee a contact that
  confers specificity.** With an 8 A cutoff a design can satisfy the term while making
  its real interface on Gly7. Phase 4's three-state screen is what actually tests this,
  and it is a prediction, not a measurement.
- 11% exposure may simply be too little for a de novo binder. Natural receptors reach it
  with a long loop evolved for the job; a hallucinated 60-80 residue protein may not.
- Nothing has been generated. No GPU has been used on this challenge.
