# Challenge 3 - GLP-1R. Methods, and what they do not cover

20 designs, all single-chain, 32-80 residues. Generated de novo by AlphaFold2-based
hallucination (ColabDesign / AfDesign, three-stage logits -> softmax -> hard). No design
is derived from, threaded onto, or initialised from any natural sequence.

## ACTIVATION IS UNSCREENED

Read this before reading anything else.

The challenge asks for molecules that **activate** GLP-1R. Nothing in this pipeline tests
activation. We have no functional assay, no cAMP readout, and no model that predicts
receptor activation from structure. What we screened is **geometry**: whether a design
occupies the extracellular groove, whether it reaches into the transmembrane channel to
the depth the native hormone reaches, and which end of the chain points down the channel.
That geometry is a **necessary** condition for the known activation mechanism. It is not
sufficient, and we are not claiming it is.

Of the 20 designs:

- **4** have the shape the activation mechanism requires (groove contact + channel reach
  + N-terminus leading). These are the agonist candidates.
- **3** contact both sites extensively but enter C-terminus first. By mechanism these are
  **blockers, not agonists**. They are submitted as binders.
- **13** bind the extracellular groove only. By mechanism these are **antagonists** - the
  molecule class of exendin(9-39), an existing drug. They are submitted as binders.

So at most 4 of 20 are agonist candidates on mechanistic grounds, and none of the 4 has
been shown to activate anything.

## Target

Full-length human GLP-1R from PDB **6X18** (active state, peptide-bound), chain R,
283 residues retained, numbered 29-406.

An earlier target was cropped to a 12 angstrom shell around the binding sites. That crop
removed the bottom of the transmembrane bundle, and channel-targeted designs dived into
the resulting hole: at the depth they reached, the full receptor has 211 atoms and the
crop had 4. The target was rebuilt with the channel floor sealed. Before spending compute
on it, the native hormone was re-docked into the sealed target and reproduced its
published contacts - the only check that the target itself is not the thing being
optimised against. Build script: `glp1r/build_sealed_target.py`.

### Epitope

Two contact sets, both read off 6X18 rather than chosen by the optimiser:

- **groove (ECD), 22 residues** - the extracellular domain surface the hormone's
  C-terminal half grips.
- **channel (TM), 30 residues** - the residues lining the path the hormone's first ~8
  residues take, ~33 angstrom down into the 7TM bundle.

Cross-checked against 3IOL (the GLP-1R ECD structure distributed by Adaptyv's own design
tool): 16 of their 17 epitope residues are in our set, Jaccard 0.70, and **zero
residue-identity mismatches**, which also serves as a numbering check. Our set adds the
ECD C-terminal stalk (P137, E138, L141, L144, Y145), which does not exist in 3IOL and so
cannot be targeted by anyone working from that structure. Write-up:
`glp1r/results/EPITOPE_CROSSCHECK.md`.

## Scoring

Every number reported comes from an **independent AlphaFold2 re-prediction at 3 recycles**
of the design-target complex, run fresh after design. The optimiser's own score is
recorded separately as `design_iptm` and was never used to select.

Geometric measures, computed by the same code for every design in the project and for the
native hormone as a reference:

| measure | what it is |
|---|---|
| groove contacts | design heavy atoms within 5 A of a groove epitope residue |
| channel contacts | same, against the 30 channel residues |
| deepest | furthest penetration along the channel axis, angstrom |
| threaded | true when the N-terminus is the deep end (native polarity) |
| hydropathy | Kyte-Doolittle mean, as an aggregation proxy |
| helicity | fraction of residues with CA(i)-CA(i+4) in 5.5-7.0 A |

Native GLP-1 in this same target: reach 33.2 A, 26 on-site contacts, i_ptm 0.761.

## Selection, and the one place we went against the obvious

Selection code: `glp1r/select_v2.py`, which prints every pass and every failure.

**Predicted confidence was not used to select the channel-reaching designs.** Across the
15 sealed-target designs, predicted interface confidence (i_ptm) correlates **-0.651**
with touching the hormone's channel contacts, and **-0.804** among the two-domain designs
alone. An earlier run against a different target gave **-0.963** across 11 designs. Two
targets, two design tasks, the same inverted sign. The highest-confidence design in the
batch touches **1** of 30 channel residues; the lowest-confidence touches **19**. An i_ptm
floor would have deleted every agonist candidate and kept the ones that miss.

The groove population keeps its confidence floor, because that rule was written before any
design existed and is the one pre-registered rule here.

| population | slots | rule | written |
|---|---|---|---|
| landing-pad (groove) | 13 | i_ptm >= 0.50, groove >= 8 of 22, hydropathy <= 0.0 | before any design existed |
| two-domain | 4 (of 5) | groove >= 4, channel >= 5, N-terminus leading, hydropathy <= 0.5 | **after** seeing the designs |
| dual-site reversed | 3 | groove >= 5, channel >= 15, C-terminus leading | **after** seeing the designs |

Near-duplicates above 80% sequence identity were collapsed; none were found.

The two-domain pool filled 4 of its 5 slots. The spare slot was given to the groove
population rather than filled with a design that failed the two-domain filter. That costs
one agonist-shaped slot and buys one more design screened by the pre-registered rule.

## Threader-only designs get zero slots

An earlier class of design targeted the channel alone. It failed twice. The second failure
was a controlled experiment whose terms were fixed in advance: same job, same method, one
variable changed (the hole in the target closed). It failed again - 4 of 7 channel-only
designs touched **zero** channel contacts. The conclusion is that the method was at fault,
not only the broken crop. Write-up: `glp1r/results/THREADER_FAILURE.md`.

One design originally generated in that class (`glp1r_thread2_008`) nonetheless came out
with native polarity, 7 groove and 20 channel contacts and 33.9 A reach, and is included
in the two-domain population on its measurements. Its class got no reserved slots; it
earned a slot on geometry.

## Novelty

Checked locally against 9 natural class-B peptides (GLP-1, GLP-2, exendin-4,
exendin(9-39), glucagon, GIP, oxyntomodulin, secretin, liraglutide): the highest local
identity of any design to any of them is **0.45**, over a best-matching window.

This is **not** a database search. It is a comparison against 9 sequences chosen by us as
the plausible accidental-copy targets for this receptor. No BLAST against a sequence
database was run.

## What is not done, not verified, or known to be weak

- **Activation is unscreened.** Stated again because it is the central limitation.
- **No experimental result exists** for any design in this submission. Every number is a
  prediction or a geometric measurement on a predicted structure.
- **Foldability is not filtered.** The design matching the native hormone's depth most
  closely (33.2 A) has helicity 0.14 and pLDDT 0.458 and may simply not fold. It was kept
  rather than cut, because cutting on it would remove the best geometric match. Per-design
  values are in `submission/challenge3_metrics.csv`.
- **Two of the three selection rules were written after seeing the data**, which makes
  them weaker than the groove rule. Marked as such in the table above and in the code.
- **`design_iptm` is blank for the 13 groove designs** in the metrics CSV. The value
  exists in the run log; it was not carried across when the records were transcribed, and
  it was not used for selection.
- The agonist candidates are all **low-confidence** designs. That is deliberate and
  justified by the anti-correlation above, but it means they would be rejected by any
  conventional confidence filter, including, possibly, the organisers'.

## Files

- `submission/challenge3_proposed.csv` - name, sequence, molecule_class. The submission.
- `submission/challenge3_metrics.csv` - every measurement, per design, with population.
- `glp1r/select_v2.py` - selection rule, with its own limitations in the header.
- `glp1r/build_sealed_target.py` - rebuilds the target and prints proof of identity.
- `glp1r/crosscheck_3iol.py` - the 3IOL epitope cross-check.
