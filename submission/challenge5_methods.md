# Challenge 5 - KRAS G12D / HLA-A*11:01 - methods

## Summary

20 de novo mini-binders designed to read the G12D aspartate of the KRAS peptide
VVVGADGVGK presented by HLA-A*11:01 (PDB 7OW6), while making limited contact with the
MHC itself. Lengths 54-103 aa. No cysteines in any design.

## Target preparation

Chains A (heavy), B (beta-2-microglobulin) and C (peptide) of 7OW6; the TCR chains D
and E were discarded. A crop sweep over eight candidate trimmings compared buried
hydrophobic surface newly exposed by each cut; the chosen 20 A crop around the peptide
(179 heavy-chain residues + the 10-mer peptide = 189) exposed 250 A^2, 1.04x the intact
complex, the least of any option tested. Solvent exposure of the peptide was measured
rather than assumed: Asp6 is 11% exposed and Gly7 73%, with the anchor residues 1, 2,
3 and 10 at 1-3%.

## Backbone generation

Backbones were taken from the 152-scaffold set published with Liu et al., *Design of
high-specificity binders for peptide-MHC-I complexes*, Science 389:386-391 (2025)
(Zenodo 10.5281/zenodo.15169815), each of which is a de novo designed protein shown to
arc over an MHC-I groove, built against pMHC targets **other than ours**. Each scaffold
was docked onto our target by superposing its own MHC platform onto ours with a local
sequence alignment driving the residue correspondence; superposition RMSD was 0.47-1.50
A over matched CA atoms (median 0.63).

Poses were then filtered on criteria fixed before any sequence existed: no backbone
clash with the target, at least 4 heavy-atom contacts to the Asp6 side chain, and at
least 20% of the interface on the peptide rather than the MHC. 151 of 152 scaffolds
docked; 111 were free of backbone clashes; 75 contacted Asp6; 63 were peptide-centric
and went forward. All 63 contact peptide positions 6, 7 and 8, and only 2 contact
position 1 - the peptide-centric footprint the source paper describes.

## Sequence design

Every residue of each binder was redesigned with ProteinMPNN (v_48_020), with the entire
target chain held fixed, at temperatures 0.1, 0.2 and 0.3, 16 sequences each. Cysteine
was excluded and the output was checked for it rather than the flag being trusted.

## Filtering

Sequences passed a composition screen developed over the four previous challenges:
two-sided limits on charged and hydrophobic fraction and on mean Kyte-Doolittle
hydropathy, limits on identical and same-charge runs, and a proline/glycine gate
(P <= 12%, G <= 12%, P+G <= 20%) applied first.

Specificity was estimated with the ProteinMPNN peptide-recovery test used by Liu et al.: with the binder, the MHC and the other nine peptide residues held fixed, we read the conditional probability of each amino acid at peptide position 6 and took log P(Asp) - log P(Gly), i.e. G12D over wild type. The same quantity measured on the bare pMHC with no binder present (+0.443) is subtracted, so the reported margin is the binder's own contribution rather than the groove's intrinsic preference.

That subtraction alone is not enough, and we checked: a nonsense poly-AEKL helix threaded onto one of these backbones scored +0.52 against the bare-groove baseline, so a large part of it rewards mere occupancy of the site. The score designs are ranked on is therefore measured against a SHUFFLE of each design's own sequence on its own backbone, holding pose, length and amino-acid composition fixed and varying only the order. What survives that subtraction is attributable to the design.

A second, independent measure was attempted: the ProteinMPNN score of the binder sequence against the G12D background and against a wild-type background built by renaming the peptide's Asp6 to Gly on the same backbone. It did not work. Across the candidate set its entire range was 0.16 log units (-0.105 to +0.059) against 1.15 for the measure above, and it correlated -0.03 with it. We report it in the detail table as `binder_fit_delta` and do NOT use it for ranking: a measure with almost no dynamic range that agrees with nothing is not corroboration.

## Ranking

Ordered by specificity margin first, then contacts to the Asp6 mutation, then the
peptide share of the interface, with the ProteinMPNN score as a tie-break. At most 3
designs share a parent backbone, so the set covers 12 distinct backbones.

## De novo declaration

Track 3 requires de novo designs and forbids modifying an existing binder. The
backbones here are de novo designs from the published scaffold set cited above, built
against different pMHC targets, and 100% of each binder sequence was rewritten for this
target. Sequence identity of each design to the scaffold it came from is 25-35%
(median 33%); designs at or above 45% identity were excluded. Per-design identities
and parent scaffold names are in the accompanying detail table. We state this
explicitly so the organisers can judge eligibility directly.

## Reproducibility

Target preparation, scaffold docking, pose filtering and all screens are
deterministic and reproducible from the public inputs (PDB 7OW6/7OW5 and Zenodo
10.5281/zenodo.15169815). ProteinMPNN sampling is stochastic, and the run that produced
these sequences was made before a fixed seed was added to our script, so the exact
sampled set cannot be regenerated bit-for-bit. All sequences that run produced are
recorded alongside the code, so the set is auditable even though it is not
reproducible; later runs are seeded.

## What we did not do

- No experimental data of any kind.
- Affinity was not predicted. No computational metric used here predicts binding
  strength, and the source paper relies on experiment for affinity.
- The specificity margin is a sequence model's conditional probability on a fixed
  backbone. It is not a KD ratio and should not be read as one.
- The wild-type comparison is implicit in the margin above; we did not build and
  predict a separate wild-type complex for these designs.
- Side chains were never packed or minimised; all interface counts are from the docked
  backbone with ProteinMPNN-assigned identities.
