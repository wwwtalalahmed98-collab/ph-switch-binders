# The probe: neither arm worked, and reading the source said why

Notebook `tfr1_probe` version 1, 2026-10-28. GPU T4 x2, **1 h 00 m 39 s, successful** - on
budget, no overrun. 19 designs, 0 failures.

## Result

| arm | n | mean on-site | >=2 | >=4 | mean off-site |
|---|---|---|---|---|---|
| run 1 baseline | 47 | 0.96 | 15 (32%) | 3 (6%) | - |
| pull (i_con weight 1 -> 3) | 8 | 1.25 | 3 (38%) | 1 (12%) | 9.6 |
| tight (radius 18 instead of 20) | 11 | **0.00** | 0 | 0 | 3.1 |

**Neither arm is a win.** `pull` is inside noise of the baseline at n=8 - 1 design of 8
reaching 4+ contacts against 3 of 47 is not a signal. `tight` is a clear loss: zero epitope
contacts across 11 designs, and only 3.1 off-site contacts, meaning those binders barely
touched the target at all.

Both hypotheses - "the pull is too weak" and "the cut face is too big" - are wrong.

## The actual cause, from the library source

The weights the run printed are the first clue:

    {'con': 0.0, 'exp_res': 0.0, 'helix': 0.0, 'i_con': 1.0,
     'i_pae': 0.0, 'pae': 0.0, 'plddt': 0.1, 'seq_ent': 0.0}

`i_con` is already the only term that matters. Scaling it 1 -> 3 scales a loss that is
nearly the whole objective, which changes almost nothing about the direction of descent.
That explains why `pull` did nothing, but not why run 1 failed in the first place.

Reading ColabDesign v1.1.1 settles that. The hotspot IS honoured - `af/loss.py:42` restricts
the interface-contact term to the hotspot positions when one is given, as assumed. The
problem is what that term counts as a contact. From `af/model.py:51`:

    "i_con": {"num":1, "cutoff":21.6875, "binary":False, "num_pos":float("inf")}

`_get_con_loss` keeps the distogram bins below `cutoff`, and the bins are
`append(0, linspace(2.3125, 21.6875, 63))`. A cutoff of 21.6875 therefore keeps **63 of 64
bins** and counts any distance **up to 21.4 A** as a satisfied contact. It is the loosest
setting the distogram can express; only "further apart than 21.7 A" is penalised at all.

For scale, the same options set the intra-binder `con` cutoff to 14.0 A, a real contact is
about 8 A between CB atoms, and we measure success at 4.5 A between heavy atoms.

**The cut face sits 6.5 A from the epitope.** A binder gripping the cut face is therefore
within 21.4 A of every epitope residue and fully satisfies the objective. There was never a
gradient pulling it onto the epitope.

So run 1's designs did not ignore the instruction. They satisfied it. We asked for "be in
the same neighbourhood as the epitope" and then measured "touch the epitope", and those are
different questions. Both probe arms also follow: `tight` made the neighbourhood smaller, so
the condition got easier still, while leaving less surface worth gripping.

This also means the earlier framing in GENERATION_RESULT.md - "the hotspot biases but does
not constrain, and the bias lost" - is wrong and should be read with this file. The bias did
not lose. It was satisfied at 21 A.

## The fix

    model.set_opt(i_con=dict(cutoff=8.0, num=1, num_pos=6), set_defaults=True)

- **cutoff 8.0** - a contact, by the usual CB-CB convention, instead of 21.4 A.
- **num_pos 6** - `min_k` keeps the k best-satisfied hotspot positions, so this asks for a
  genuine six-residue footprint on the epitope. The default `inf` averages over all 13, which
  a 65-residue binder cannot reach at once, and an unreachable target gives a flat gradient
  in the same way an already-satisfied one does.
- **set_defaults=True** - `restart()` inside the design stages resets options otherwise.

Not verified: that this produces on-target binders. It removes a specific, identified reason
they could not. The run below prints on-site contacts per design so it can be killed early
if the first few still read zero.

## Cost

Probe: 1 h of the 6 h that remained. Quota was 23:10 of 30 at launch.
