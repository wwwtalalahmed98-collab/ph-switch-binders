# The rewrite step, and what it proved about the backbones

Kaggle `notebook2717643405`, CPU only, no quota. 11 minutes.
Script: [kaggle/cellPMHCe_mpnn.py](../../kaggle/cellPMHCe_mpnn.py).

## Why this ran

The raw hallucinated sequences are not submittable:

```
...GGGPGGGGGGGGCCGHHHHHHHHHHHRRRRRRRRRRRRRRSGSG
CCGGRGGGGGGGGGFDDWDDWFDDYYDFFWPYNDDDD...RRRRRRRRRRC
MKRRRRAAAAAAAAMLAFWAAWWAGYAGY...            <- the pLDDT 0.63 one
```

Eleven histidines in a row, fourteen arginines, eight alanines. The composition screen
scores the best of them `fail (8 identical; 59% hydrophobic)`. ProteinMPNN holds the
backbone fixed and writes a real sequence for it, costs no GPU, and takes minutes. Not
running it was the same gap that went unnoticed for three challenges on TfR1.

## Result

15 backbones (every design touching the peptide, with a structure attached), 48
sequences each at two temperatures. **720 sequences. Zero good, zero marginal.**

```
pmhc2_W_cut8_003   48 sampled: 0 good 0 marginal 48 fail  0 fold-breaker
pmhc2_W_cut8_015   48 sampled: 0 good 0 marginal 48 fail  0 fold-breaker
pmhc2_W_cut8_011   48 sampled: 0 good 0 marginal  0 fail 48 fold-breaker
pmhc2_R_cut10_004  48 sampled: 0 good 0 marginal  9 fail 39 fold-breaker
pmhc2_R_cut10_013  48 sampled: 0 good 0 marginal  0 fail 48 fold-breaker
pmhc2_W_cut8_006   48 sampled: 0 good 0 marginal  5 fail 43 fold-breaker
pmhc2_W_cut8_012   48 sampled: 0 good 0 marginal 21 fail 27 fold-breaker
...
SystemExit: nothing survived composition on any backbone
```

## What the failure mode means

Most failures are **fold-breakers** - overwhelmingly proline and glycine.

MPNN places glycine and proline where the backbone's dihedral angles demand them. Give
it a real fold and it returns normal composition; give it a backbone that is not a fold
and it returns Gly/Pro soup because nothing else fits the geometry.

So this is **a readout of backbone quality, from a different program with a different
objective**. The earlier conclusion - "hallucination cannot aim at Asp6" - was true but
understated. The shapes are not proteins either.

That reframes the whole challenge's results: it was never only the sequence half.

## What it does not mean

It does not mean the composition screen is too strict. The screen's limits are the ones
that passed a real TfR1 candidate (`SPTKTYVSIEK...`, scored `good`) in the same test run
that rejected these. It is a list of ways to be obviously wrong, and these sequences hit
them.

## The one setting this points at

The weights printed in run 2 included `'helix': 0.0` - no bias toward secondary
structure at all, so nothing stops the optimiser producing backbones with neither
helices nor sheets. Run 3 sets `helix=0.3`.

Stated before that run: **it addresses foldability, not aim.** The 11-22% mutation-touch
rate has no reason to improve. The expected outcome is buildable designs that mostly
miss, against the current position of zero buildable designs of any kind.

Success criterion fixed in advance: MPNN on the new backbones must return something
other than Gly/Pro soup. If it is again zero good and zero marginal, this approach is
finished and the entry becomes the negative result plus whatever is least bad.
