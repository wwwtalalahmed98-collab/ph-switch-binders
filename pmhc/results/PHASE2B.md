# Generation run 2 - the rebalance, and the answer it was built to give

Kaggle `notebook650f87c002` version 4, GPU T4 x2, 5 h 09 m. 34 designs, 0 failures.

## The defaults, now measured rather than inferred

```
weights BEFORE {'con': 0.0, 'exp_res': 0.0, 'helix': 0.0, 'i_con': 1.0,
                'i_pae': 0.0, 'pae': 0.0, 'plddt': 0.1, 'seq_ent': 0.0}
weights AFTER  {'con': 0.5, 'exp_res': 0.0, 'helix': 0.0, 'i_con': 0.5,
                'i_pae': 0.4, 'pae': 0.1, 'plddt': 0.5, 'seq_ent': 0.0}
```

ColabDesign's binder default is the contact term at 1.0 against pLDDT at 0.1, with every
other term at zero. The run-1 diagnosis - one instruction crowding out all others - was
correct, and is now a measurement rather than a reading of the loss.

## Result, against criteria fixed before the run

| arm | n | reached the mutation | mean on-mut | mean on-peptide | mean pLDDT | i_con last |
|---|---|---|---|---|---|---|
| W_cut8 | 18 | 4/18 | 0.22 | 1.44 | 0.45 | 4.62 |
| R_cut10 | 16 | 2/16 | 0.12 | 1.12 | 0.45 | 4.30 |
| *run 1, for comparison* | 28 | 4/28 | 0.14 | 1.86 | ~0.40 | ~4.5 |

**Question 1 - was fold pressure the binding constraint?** Mean pLDDT moved from 0.40 to
0.45. Marginal. **No.**

**Question 2 - was the 8 A demand geometrically impossible?** i_con stalls at 4.62 with an
8 A cutoff and at 4.30 with 10 A. Relaxing the demand by 2 A changed essentially nothing.
**No.**

**Reaching the mutation did not improve**: 4/18 and 2/16 against run 1's 4/28.

## The conclusion, taken as written in advance

The run's own header said: *"If it stalls at both cutoffs, the method is wrong for this
target and no amount of tuning fixes it."* It stalled at both.

**AfDesign hallucination cannot reliably place a binder on Asp6 in this groove.** Two
rounds, 62 designs, three objective configurations, and the rate of touching the mutation
sits between 11% and 22% with one residue of contact when it happens. Further tuning of
this method is not indicated, and continuing to tune it after setting that criterion
would be moving the goalposts.

The likely reason is structural rather than numerical: Asp6 is 11% exposed at the bottom
of a groove whose walls are identical in the wild type. Reaching it needs a protruding
loop, which is what natural receptors use and what purpose-built TCR-mimic pipelines
grow deliberately. Hallucinating a compact blob and hoping it reaches in is the wrong
shape of method.

## What the two rounds did produce

62 designs, all with PAE saved. A handful are materially better than anything earlier:

```
W_cut8_015   ON-MUT 1 | pept 5 groove 16 | i_ptm 0.47 plddt 0.73
W_cut8_003   ON-MUT 1 | pept 6 groove 11 | i_ptm 0.50 plddt 0.63
A_tight_010  ON-MUT 1 | pept 6 groove  9 | i_ptm 0.56 plddt 0.64   (run 1)
A_tight_013  ON-MUT 1 | pept 6 groove  7 | i_ptm 0.32 plddt 0.64   (run 1)
W_cut8_011   ON-MUT 1 | pept 5 groove 10 | i_ptm 0.29 plddt 0.44
```

`W_cut8_015` at pLDDT 0.73 is the first design in this challenge that is plausibly a
folded protein *and* touches the mutation. One design out of 62.

One design, `R_cut10_013`, makes 35 groove contacts and 4 peptide contacts - the failure
the challenge exists to punish.

## What this does not say

It does not say the designs are useless. It says this generator cannot be *driven* onto
the mutation reliably. Whether the handful that did land there survive the three-state
counter-screen is a separate question and is unmeasured.

## Next

Counter-screen the shortlist rather than all 62 - three states and two parameter sets per
design is six predictions each, and 62 designs would not fit the remaining budget.
