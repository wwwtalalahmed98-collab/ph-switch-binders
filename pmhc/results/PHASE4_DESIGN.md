# Phase 4 - the three-state counter-screen. Built and tested, not yet run.

Built 2026-11-04 while the generation run occupied the GPUs. CPU parts tested; the GPU
part has not executed yet.

## What it scores

The competition ranks on two ratios before it looks at strength, so the quantity to
predict is a *margin* between near-identical complexes, not an affinity.

| state | contents | required |
|---|---|---|
| g12d | HLA + b2m + VVVGADGVGK | must bind |
| wt | HLA + b2m + VVVGAGGVGK | must not bind |
| apo | HLA + b2m | must not bind |

All three share one background - the 275 heavy-chain and 99 b2m residues that both 7OW6
and 7OW5 model. Verified by [build_counter_states.py](../build_counter_states.py):

```
target_state_g12d   275  99  10  VVVGADGVGK
target_state_apo    275  99   0  -
target_state_wt     275  99  10  VVVGAGGVGK
```

## The control that makes the margin mean anything

AlphaFold's confidence differs between its own parameter sets on identical input. If
i_ptm moves by 0.04 between model 1 and model 2 on the *same* state, a 0.03 difference
between g12d and wt is noise with a story attached.

So every state is predicted with two parameter sets and **the model-to-model spread is
reported beside the margin**. A design is marked `MEASURABLE` only when both its margins
exceed its own spread. Everything else is reported as within noise - not as a small
effect.

No previous challenge in this project had a noise floor. It is the difference between a
ranking and a ranking that can be defended.

## Tested on a processor, before any GPU time

[test_screen.py](../test_screen.py) builds the merged-chain layout `save_pdb` really
produces - a 384-residue target chain with the peptide as its last ten residues - and
checks five cases:

```
g12d, binder on Asp6   -> on_mutation 1 on_peptide 5 on_groove 5
wt,   binder on pos 6  -> on_mutation 1 on_peptide 5 on_groove 5
apo,  binder in groove -> on_mutation 0 on_peptide 0 on_groove 5
g12d, binder 60 A away -> all zero
g12d structure offered as wild type -> rejected
```

The last two matter most: the apo state must not invent peptide contacts it has no
peptide for, and a state mix-up must raise rather than quietly produce a number.

[test_measure.py](../test_measure.py) covers the same class of failure for the
generation run, where it has already caught one real bug.

## Limits, for the methods text

- **The wild-type state is modelled, not crystallographic.** It is 7OW6 with Asp6 cut
  back to Gly, so the substitution is the only difference between the two states. The
  cost is that it carries the mutant's backbone. 7OW5 exists and is the real wild type;
  using it would confound the margin with a whole second crystal. This trade buys an
  interpretable margin at the price of a small structural inaccuracy.
- **The apo state is the weakest of the three.** Deleting a peptide leaves a groove
  still shaped as though one were in it. Real peptide-free HLA is unstable and partly
  collapsed, so our apo test is easier than the assay's, and passing it does not show a
  design ignores empty HLA.
- **AlphaFold is not known to predict specificity margins between near-identical
  complexes.** A predicted margin is a hypothesis. The noise floor bounds how small a
  margin we are entitled to report; it does not make the prediction correct.

PAE is saved for every prediction, so ipSAE can be computed afterwards on a processor.
