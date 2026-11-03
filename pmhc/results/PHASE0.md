# Phase 0 - the three states, and which crop we design against

Run 2026-11-03, CPU only, no quota spent.
Scripts: [build_states.py](../build_states.py), [crop_sweep.py](../crop_sweep.py).
Outputs land in `data/pmhc/states/`, which is gitignored because the scripts rebuild it.

## Verified, not assumed

| check | result |
|---|---|
| 7OW6 chain C peptide | `VVVGADGVGK` - matches the spec |
| 7OW5 chain C peptide | `VVVGAGGVGK` - matches the spec |
| the two differ | at position 6 only, Asp vs Gly |
| chain roles | A heavy chain (276), B beta-2-microglobulin (99), C peptide (10) |
| TCR | chains D and E, discarded |

The script refuses to write anything if any of these fail.

## One background, not two

7OW6 models heavy-chain residue 1 (a `GS` cloning remnant) and 7OW5 does not. Using
each crystal for its own state would mean the mutant and wild-type backgrounds differed
by more than the mutation, and any specificity we measured could come from that instead.

So all three states share **one** background: the 275 heavy-chain residues and 99 b2m
residues that **both** structures model. 7OW5 is used only to confirm the wild-type
sequence and to tell us which residues are common. The wild-type state is the same
background with Asp6 changed to Gly.

## The crop decision, and a measurement error that reversed it

The full complex is 384 residues. Every design has to be folded three times, so the
full complex would buy very few designs in the days left. Cropping is not optional.

The cost of cropping is hydrophobic surface that is buried in the real protein and
exposed in the crop. On challenge 3 the optimiser found exactly such a face and dived
into it.

| crop | residues | new hydrophobic | vs peptide's own surface |
|---|---|---|---|
| full complex | 384 | 0 | 0.00x |
| **20 A around peptide** | **189** | **250 A^2** | **1.04x** |
| platform 1-180 | 189 | 255 A^2 | 1.07x |
| 18 A around peptide | 183 | 278 A^2 | 1.16x |
| platform + b2m | 288 | 309 A^2 | 1.29x |
| 16 A around peptide | 176 | 326 A^2 | 1.36x |
| 24 A around peptide | 202 | 478 A^2 | 1.99x |
| 14 A around peptide | 157 | 550 A^2 | 2.30x |

**The first version of this table was wrong and pointed the opposite way.** It counted
new exposure on the heavy chain only, and reported `platform + b2m` at 14 A^2 - an
apparently clean crop, which is what I was about to pick. Keeping b2m while dropping the
a3 domain exposes **b2m's own underside**, and counting only chain A hid it. Corrected to
count every chain, that crop is 309 A^2: both larger and dirtier than the alternatives,
so it is dominated rather than best.

**Chosen: the 20 A crop around the peptide**, 179 heavy-chain residues (2-181) plus the
10-residue peptide, 189 residues total. It is a coin-flip with the platform crop on the
numbers; it wins on being peptide-centred rather than domain-defined.

**No crop is clean.** All of them open roughly the size of the peptide's own footprint.
We accept that here for a reason that would not hold on an affinity-scored challenge:
the cut face is identical in all three states, so a design that grips it binds all three
and the counter-screen rejects it. It costs design budget, not correctness.

The three-state counter-screen in phase 4 runs against the **full** complex, not the crop.

## Hotspots measured

Peptide Asp6 plus the 18 heavy-chain residues within 8 A of it. Closest first:

```
R114 2.95   Q70  3.02   Q155 3.80   T73  4.03   Q156 4.22   W147 4.36
D77  4.76   I97  4.79   D116 5.01   D74  5.38   A69  6.50   A152 6.52
W133 6.67   N66  6.90   Y9   7.33   Y99  7.66   I95  7.92   S71  7.99
```

Choosing which of these to aim at is phase 1, not this phase. The shape of that choice
is already visible though: the specificity story is "we touch Asp6", so the hotspot set
must be weighted to the peptide rather than spread over 18 groove residues. On TfR1,
averaging the contact term over all 13 hotspots is part of why round 1 aimed at nothing.

## Not done

- Nothing has been generated. No GPU has been used on this challenge.
- The hotspot set is measured but not chosen.
- The apo (peptide-free) state is defined as a sequence set only. Real peptide-free
  HLA is conformationally unstable, and AlphaFold will happily predict a tidy empty
  groove. Our apo counter-screen therefore tests something easier than the assay does.
  That belongs in the methods text.
