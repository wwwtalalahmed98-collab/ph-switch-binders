# Challenge 2 (TNF-alpha) — methods

**Track 3. Twenty single-chain de novo designs, 56–79 aa: twelve acid-OFF pH switches and
eight histidine-free controls, each control the exact twin of one switch.** Ranked
switch-first, per the challenge's stated order.

## Summary

The binders are AfDesign hallucinations against the receptor-binding site at the protomer
interface of human TNF-alpha. The switch is a histidine placed where **binding buries it**
next to a target arginine or lysine: in the free binder the histidine is solvent-exposed
and protonates normally at pH 6.0, and binding then has to bury that charge, which costs
affinity. Eight switches are submitted with a twin identical except at the histidine
positions, verified programmatically, so a measured pH response is attributable.

This submission contains **two negative results** that we think are more useful than the
designs. Both have the same shape: a cheap geometric proxy we had been selecting on turned
out to carry almost no information, and a direct measurement replaced it.

## Negative result 1: epitope conservation does not predict mouse binding

We first ranked mouse cross-reactivity using the fraction of contacted target residues
identical between human and mouse, reasoning that a conserved footprint should transfer.

Across the 11 designs with both a conservation score and a measured mouse prediction,
**Pearson r = +0.272**. The design the proxy ranked first had the worst measured mouse
binding, making **zero** contacts with the mouse target.

Measuring the parents separately located the cause. The parents were fine — 20 of 31 bind
mouse. The loss came from our own switch step: across 17 lineages the same histidine
substitutions retained a median of **90% of human binding and 49% of mouse**. We had
scored candidate positions on the human structure and accepted them on human results
alone, so positions human tolerates and mouse does not passed the filter by construction.
Rebuilt to test every placement against both species, one lineage moved from mouse i_ptm
0.148 with 0 mouse-site contacts to 0.655 with 15, by moving a single histidine.

## Negative result 2: distance to a cation does not predict the switch

Having fixed the species problem, we tested the switch rule itself. AlphaFold2 does not
model pH, but it does respond to charge, so each histidine position was predicted as
glutamine (neutral, similar size) and as lysine (charged, similar length). The difference
isolates charge from sterics: a position that rejects a positive charge is one where a
protonated histidine could drive release.

Across **35 verified placements**, exactly **one** shows an effect clearly above our noise
floor — `tnf2_w0_001` Y37H at Lys-Gln = **-0.264**, which reproduced independently at
-0.198 in a separate run. The median is -0.019, and in **8 of 35 a positive charge
helps**. Our measured reproducibility between identical re-predictions is 0.021 on human
and 0.026 on mouse, so most of these effects are indistinguishable from noise.

Correlating that charge effect against the variables available:

| variable | correlation with charge effect | used in the original design? |
|---|---|---|
| **burial on binding** | **-0.58** (n=17; -0.53 excluding one bad geometry) | no |
| screened Coulomb term | -0.28 | no |
| Cbeta-to-cation distance | **+0.24** | **yes — this was the rule** |

The 3.0–8.8 A distance window we had designed with is about as uninformative for the
switch as epitope conservation was for mouse binding. Burial on binding is the better
variable, it is geometry rather than a model's opinion, and it costs no GPU time.

Applying it exposed a concrete failure: **7 of 35 placements have a histidine buried by
less than 0.02**, meaning it sits in open solvent when bound, protonates at no cost, and
cannot contribute to release. Four such placements were in our previous submission.

## Selection rule for this set

1. **Conserved anchor first.** Keep only placements whose histidine faces one of the two
   cations present in both species (K166 chain A, R108 chain B). A histidine facing a
   human-only cation cannot push off the mouse target, and mouse ranks above affinity.
2. **Burial on binding**, worst histidine first.
3. Rank lineages by that burial and take the top twelve.

All 12 switches face a conserved anchor, **none has a histidine buried below 0.059**, and
every switch binds mouse at i_ptm >= 0.452 with at least 6 on-site mouse contacts.

## The cost of this rule, stated plainly

**This set contains no two-histidine designs**, where our previous one had six. That is a
real loss and we want it on the record.

The argument for two histidines is Henderson-Hasselbalch: at histidine pKa ~6.0, one
histidine leaves 50% of molecules uncharged and still binding at pH 6.0, two leaves 25%.
The challenge asks for *no detectable* binding. We have given that up.

We did so because that argument assumes both histidines are doing work, and in our data
none of the six pairs satisfies it. Ranked on the weaker of their two histidines, the best
pair reaches 0.100 where the best single reaches 0.404, and the pair in our previous
top-ranked design scored **0.000** — one of its two histidines was inert. With a dead
partner the effective count is one and the pH 7.4 binding cost is two. On this evidence
the best single from each lineage is better than the pair from the same lineage, and it is
also better on *measured* binding in six of the fourteen lineages.

If the burial reasoning is wrong, this trade was a mistake and the all-singles set will
leak at pH 6.0. That is the main risk in this submission.

## Pipeline

1. **Epitope and construct.** Human P01375 and mouse P06804 aligned; the mouse target
   built from PDB `2TNF` by alignment-based numbering, never arithmetic offsets. Target is
   the recommended receptor-binding site between two protomers, cropped to 122 residues
   for human and 123 for mouse with all 25 site residues retained in both.
2. **Generation.** AfDesign binder hallucination (ColabDesign), 3-stage, 2xT4, three
   rounds, 89 candidates at 56–80 aa.
3. **Scoring.** Every design re-predicted from sequence alone at 3 recycles, never the
   optimiser's own score. 31 of 89 reached human i_ptm >= 0.45 with 5+ on-site contacts.
4. **Mouse, measured.** All 31 parents predicted against the mouse construct with
   identical settings. 20 of 31 bind mouse.
5. **Switch placement.** Up to four placements per parent, each re-predicted against
   **both** species, kept only if it retains >= 70% of parent i_ptm on both with >= 5
   on-site mouse contacts. 35 placements across 14 lineages survived.
6. **Burial and charge.** Every surviving placement measured for histidine burial on
   binding, real Nd1/Ne2-to-cation distance after folding, and charge sensitivity.
7. **Liability repair.** Free cysteines always replaced; repetitive stretches broken with
   charge-preserving substitutions only when that costs less than 15% of parent i_ptm.
   The repair plan is decided once per lineage on the control and applied at identical
   positions to every switch, so each pair differs only at the histidines. Checked
   programmatically for all 20 submitted designs.

## Stated limitations

- **Nothing has been physically tested.** Every number is an AlphaFold2 prediction against
  a static structure. i_ptm is a confidence score, not an affinity, and is not independent
  validation: it comes from the same model that generated the backbones.
- **The switch is still not computed.** No calculation here models pH or protonation. The
  charge test is an AF2 proxy for chemistry AF2 does not simulate, and burial is a
  *correlate* of that proxy. Two unvalidated links. We looked for support at 35 placements
  and found it at one; AF2's known weak response to point mutations means this cannot
  separate "will not fire" from "cannot be seen".
- **Predicted fold confidence is low by this field's standards.** Median pLDDT 0.693 with
  none above 0.80, where published minibinder campaigns filter at 0.85 before ordering.
  The binders do fold without the target — pLDDT bound minus alone has median -0.076, so
  they score slightly *higher* alone — but the absolute confidence is modest.
- **No ProteinMPNN resequencing.** The standard step that would raise fold confidence and
  expression. We applied it incorrectly earlier in this competition, without the target in
  context, saw it destroy interfaces, and set it aside; re-running it properly would
  regenerate every sequence and invalidate all 35 placement measurements. This is the
  largest known gap.
- **Every design carries a composition liability** — poly-Ile/Leu faces, lysine runs,
  alternating low-complexity, hydrophobic stretches up to 14 residues. One such repeat was
  introduced by our own repair step, which broke a poly-glutamine run into an alternating
  QNQNQN pattern — reported as a defect in our method, not a fix.
- **Four of the twelve switches have no paired control**, because their parent sequence
  already contains a histidine and the twin could not serve as a pH-inert baseline.
- **11 of 41 histidines across all placements face a cation that is not one of the two
  conserved anchors.** Those are excluded from this set, but it means the conserved-anchor
  constraint is doing real work and was not satisfied by default.
- **Trimer avidity works against objective 1.** The assay presents three equivalent
  epitopes, so partial affinity loss at pH 6.0 can be rescued into "still detectable".
  With an all-singles set this risk is higher than it was for our previous submission.
- **One predicted structure is wrong**: a histidine 0.86 A from a cation is a clash, not a
  contact. That placement is not in this set.
