# Challenge 2 (TNF-α) — methods

**Track 3. Twenty single-chain de novo designs, 56–78 aa: ten acid-OFF pH switches, each
with its own histidine-free twin as a paired control.** Ranked switch-first, per the
challenge's stated order.

## Summary

The binders are AfDesign hallucinations against the receptor-binding site at the
protomer interface of human TNF-α. The pH switch is a histidine placed within reach of a
target arginine or lysine: neutral imidazole is tolerated at pH 7.4, imidazolium clashes
with the cation at pH 6.0 and pushes the binder off. Each switch is submitted alongside
a twin that is identical except for the histidine positions, so that a pH response
measured in the lab can be attributed to the histidines rather than to the scaffold.

The substantive methodological point of this submission is negative and is described in
full below: **we used sequence conservation as a proxy for mouse cross-reactivity, and it
did not work.** Replacing it with direct measurement changed half of our design set.

## Pipeline

1. **Epitope and construct.** Human P01375 and mouse P06804 aligned; the mouse target
   was built from PDB `2TNF` by alignment-based numbering, never arithmetic offsets.
   Target is the receptor-binding site between two adjacent protomers, as recommended,
   cropped to a two-chain construct: 122 residues for human (A 68, B 54), 123 for mouse
   (A 68, B 55), with all 25 site residues retained in both. The two conserved cationic
   anchors the switch relies on are K166/R108 in the human construct and K168/R111 in
   the mouse, confirmed present in both.
2. **Backbone and sequence generation.** AfDesign binder hallucination (ColabDesign),
   3-stage design, 2×T4 GPUs, three generation rounds, 89 candidate binders at lengths
   56–80.
3. **Scoring.** Every design re-predicted from sequence alone at 3 recycles — never the
   optimiser's own score. 31 of 89 reached human i_ptm ≥ 0.45 with at least 5 contacts on
   the intended site.
4. **Mouse cross-reactivity, measured.** All 31 parents predicted against the mouse
   construct using the identical model and settings, so human and mouse numbers are
   comparable. 20 of 31 bind mouse at i_ptm ≥ 0.45 with 5+ on-site contacts.
5. **Switch placement.** For each parent, interfacial positions whose virtual Cβ faces a
   solvent-exposed target cation at 3.0–8.8 Å with orientation cos ≥ −0.30 were mutated
   to histidine. Up to four placements were generated per parent — up to three single
   histidines plus the best-scoring pair — and **each was re-predicted against both
   species**. A placement is kept only if it retains ≥ 70% of the parent's i_ptm on human
   *and* on mouse, with ≥ 5 on-site mouse contacts. Across 29 parents, 35 placements held
   both species and 14 parents yielded at least one usable switch.
6. **Liability repair.** Free cysteines are always replaced (C→S). Repetitive stretches of
   5+ identical residues are broken with charge-preserving substitutions (E↔D, K↔R), but
   only if that costs less than 15% of the parent's human i_ptm; otherwise it is reverted.
   The repair plan is decided once per lineage on the control and applied at identical
   positions to every switch in it, so each switch/control pair differs *only* at the
   histidine positions. This is checked programmatically for all 20 submitted designs.

## Selection rule for the ten slots

Applied in the order the challenge ranks the objectives:

1. **Two histidines preferred over one.** At histidine pKa ≈ 6.0, a single histidine
   leaves 50% of molecules uncharged and still fully binding at pH 6.0; two leave 25%.
   The requirement is *no detectable* binding at pH 6.0.
2. **More histidines facing a conserved anchor preferred.** A histidine facing a cation
   present only in human cannot push off the mouse target, so the switch would be
   human-only and objective 2 would be met at pH 7.4 but not mechanistically shared.
3. **Higher worse-species i_ptm**, i.e. `min(human, mouse)`.

Six of the ten switches carry two histidines. Every switch binds mouse at i_ptm ≥ 0.480
with at least 6 on-site mouse contacts.

## The finding we think is worth reporting

Our earlier ranking used **epitope sequence conservation** — the fraction of contacted
target residues identical between human and mouse — as a stand-in for mouse binding,
on the reasoning that a conserved footprint should transfer between species.

Across the 11 designs for which we have both a conservation score and a measured mouse
prediction, the two barely track each other: **Pearson r = +0.272**. The design our
conservation score ranked first had the worst measured mouse binding in the set, making
**zero** contacts with the mouse target.

Measuring the parents separately located the cause. The parent binders were not the
problem — 20 of 31 bind mouse. The loss was introduced by our own switch step: across
the 17 lineages with both numbers, the same histidine substitutions retained a median of
**90% of human binding and 49% of mouse binding**. That asymmetry is not chemistry. We
had scored candidate histidine positions on the human structure and accepted or rejected
them on human results alone, so a position human tolerates and mouse does not passed the
filter by construction.

One lineage shows the size of the effect with everything else held constant — same
parent, same repair, same machinery, only the chosen histidine position differing:

| placement | human i_ptm | mouse i_ptm | mouse on-site contacts |
|---|---|---|---|
| L30H + D27H | 0.487 | 0.148 | 0 |
| L30H alone | 0.585 | 0.655 | 15 |
| E26H | 0.566 | 0.679 | 15 |

A design that binds mouse well was about to be submitted in a state where it barely
touches mouse at all, because of where one amino acid was placed. Testing placements
against both species, rather than against human and checking mouse afterwards, is the
one thing we would tell another team to do differently.

## Stated limitations

- **No design here has been physically tested.** Every number is an AlphaFold2
  prediction against a static structure. i_ptm is a confidence score, not an affinity,
  and is not independent validation: it comes from the same model that generated the
  backbones.
- **The pH switch has not been computed, only reasoned.** No calculation in this project
  models pH or protonation. The switch rests on measured interatomic geometry plus
  textbook histidine chemistry. It is the part of this submission most likely to fail.
- **Two of the ten controls carry a native histidine** (pairs sw_07 and sw_10); their
  parent sequences already contained one. For those two pairs a pH response cannot be
  attributed to our added histidine alone. The other eight controls are histidine-free.
- **The selection rule trades measured binding for an unvalidated switch argument.**
  Preferring two histidines cost real predicted binding in several lineages — sw_05 gives
  up 0.561→0.503 on human and 0.587→0.480 on mouse against its own single-histidine
  alternative. If the two-histidine reasoning is wrong, those were losses.
- **The set is optimised for binding while carrying histidines, not for switch strength
  directly.** Our geometric switch score was used only through the conserved-anchor
  count; the highest-scoring placements usually broke the interface, so several submitted
  designs score low on it.
- **sw_09 has only 6 on-site mouse contacts**, the weakest in the set. Its parent's
  best-binding placement was passed over because that histidine faces a non-conserved
  cation — a judgement about an unmeasured mechanism.
- **A free-cysteine repair can be expensive and is not optional.** Two otherwise strong
  lineages were lost to it: one fell from 0.714 to 0.463 human i_ptm on a single C→S, and
  one from 0.469 to 0.177 on six. We had assumed the repetitive-stretch substitutions
  were the costly part and built a revertible repair to address them; that assumption was
  wrong, and the fix did not recover either lineage.
- **Affinity is assayed against the soluble trimer**, which presents three copies of the
  epitope; designs were made against a two-chain crop of one interface. Avidity effects
  are not modelled.
- **Conservation-based ranking is reported above as a failure, and the designs it
  produced are not in this submission**, except where direct measurement independently
  supported them.
