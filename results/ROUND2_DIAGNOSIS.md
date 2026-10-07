# Redesign run (notebook9e2d948101) — result and diagnosis

Ran 1227 s on T4 x2. 13 backbones x 8 ProteinMPNN sequences = 104 AF2 re-scores.
All steps completed; nothing crashed.

## Re-scored i_ptm vs round-1 i_ptm

| backbone  | round-1 i_ptm (logged by AfDesign) | honest re-score (best of 8 MPNN seqs) |
|-----------|------|-------|
| egfr_a_15 | 0.725 | 0.100 |
| egfr_a_18 | 0.604 | 0.170 |
| egfr_a_9  | 0.601 | 0.143 |
| egfr_a_12 | 0.531 | 0.273 |
| egfr_a_20 | 0.477 | 0.126 |
| egfr_a_16 | 0.472 | 0.129 |
| egfr_a_8  | 0.462 | 0.257 |
| egfr_a_5  | 0.461 | 0.143 |
| egfr_a_7  | 0.442 | 0.099 |
| egfr_a_0  | 0.422 | 0.110 |
| egfr_a_4  | 0.416 | 0.179 |
| egfr_a_17 | 0.305 | 0.100 |
| egfr_a_11 | 0.289 | 0.126 |

Switch placement: 13/13 backbones returned `0 contacts, 0 carboxylates` ->
0 switch designs, 13 controls.

## What this means

1. The round-1 numbers were NOT independent measurements. They came from
   `model.aux["log"]` at the end of `design_3stage`, i.e. the value the optimiser
   was maximising, on a sequence that had only 5 hard-stage iterations and so was
   still partly continuous. Rank correlation with the honest re-score is near zero.

2. The honest test (`af.predict` on the discrete sequence, binder folded from
   sequence alone against the target template) is the same test BindCraft and the
   RFdiffusion papers use for filtering. Typical pass threshold is i_ptm >= 0.5.
   Our best is 0.273. By AF2's own assessment none of these designs bind.

3. Zero carboxylate contacts is a downstream consequence, not a separate failure:
   the binders are not docking on the intended epitope, so no binder position is
   within 3.5-7.0 A of E45/D46/D75/E97. The pH-switch mechanism cannot be placed
   on a binder that does not bind.

## Ruled out as explanations

- Corrupted backbones. The log line `{'A': 322, 'B': 118}` is exactly 2x the
  expected 161 + 59, i.e. AfDesign wrote 2 MODELs per PDB and the naive CA-line
  counter double-counted. ProteinMPNN handled it correctly: the output sequences
  are 95/81/84 aa, single-copy, inside the 55-95 design range. Chain assignment
  (binder = B, target = A) was still correct.
- Re-score harness error. 104/104 predictions completed, lengths and chains
  consistent, control/switch bookkeeping intact.

## Root cause

Under-optimisation. Round 1 ran `design_3stage(soft=50, temp=25, hard=5)` - 80
iterations, roughly a third of the ColabDesign default, chosen to fit 24 designs
into the 12 h Kaggle batch limit. That saved wall-clock at the cost of the thing
being optimised. The hard stage in particular (5 iterations) never converged the
sequence to a discrete one, which is precisely why the logged i_ptm was inflated.

Secondary: scoring used `num_recycles=0` (ColabDesign default). BindCraft designs
at 1 recycle and filters at 3. This depresses i_ptm somewhat but cannot account
for 0.725 -> 0.100.

## Round 2 — what is now running

Notebook `notebookb5547d81b2`, Version 1, GPU T4 x2. Code: `kaggle/cellF_round2_generate.py`.
Kaggle reported 20 h of weekly GPU quota remaining; quota is billed on session
wall-clock, not per GPU, so the 7 h budget leaves ~13 h in reserve.

| | round 1 | round 2 |
|---|---|---|
| GPUs used | 1 | 2 (one worker process each) |
| iterations | 50 soft / 25 temp / 5 hard | 90 / 40 / 20 |
| binder length | 55-95 | 52-72 |
| early abort | none | drop trajectory if soft-stage i_ptm < 0.25 |
| recorded i_ptm | optimiser's own value | independent re-prediction of the discrete sequence at 3 recycles |
| failure visibility | at the end | every design appended to JSONL as it finishes |

The hard stage going 5 -> 20 is the substantive change; the honest in-line scoring
is what makes the result trustworthy either way.

## Decision point

Results stream to `/kaggle/working/out2/worker*.jsonl`. After roughly an hour there
should be 3-5 designs per worker. If `TRUE_iptm` is still clustered near 0.1 then
depth was not the binding constraint and continuing the run is not worth the quota -
stop it and fall back to submitting the best available 20 with the metrics reported
honestly. If the numbers have moved into 0.3-0.5 the run should go to completion.

## Round 2, first attempt — aborted after 56 min

Two faults, found by checking the log an hour in rather than at the end:

1. `AttributeError: 'mk_af_model' object has no attribute 'design_temp'`. The guard
   tested `hasattr(design_soft) and hasattr(design_hard)`, both of which exist, but
   `design_temp` does not exist in this ColabDesign build. The one trajectory that
   passed the quality gate crashed immediately afterwards - the worst possible case.
2. The early-abort gate was mis-calibrated. Soft-stage i_ptm came in at 0.085-0.137
   on every trajectory, all below the 0.25 threshold, so 5 of 6 were discarded. A
   soft-stage value is not comparable to a final one; the gate was measuring the
   wrong thing and would have produced ~0 designs in 6 hours.

Fixes: use `design_3stage`, the only staged entry point that exists; remove the gate
entirely; leave the loss weights at ColabDesign's binder defaults so that depth is
the single variable under test. Confirmed from the log that the weight keys are
`con, exp_res, helix, i_con, i_pae, pae, plddt, seq_ent`.

## Round 2, running

Notebook `round2` (imported from `kaggle/round2.ipynb`), Version 2, GPU T4 x2,
6 h budget, 19 h of weekly quota left. Two workers, `design_3stage(90, 40, 20)`,
binder length 52-72, every design scored by independent re-prediction at 3 recycles
and appended to `/kaggle/working/out2/worker*.jsonl` as it completes.

Note on tooling: Kaggle's web editor stopped rendering notebook cells partway
through this session, and clipboard writes from the page were blocked. Importing a
locally built `.ipynb` through Kaggle's import dialog turned out to be the reliable
way in, and is how the notebook should be updated from here on.
