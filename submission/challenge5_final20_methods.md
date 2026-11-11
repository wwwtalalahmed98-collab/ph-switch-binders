# Challenge 5 — KRAS G12D (VVVGADGVGK) on HLA-A\*11:01 — methods

## Summary

Twenty designs, from two independent routes, ten each.

| | route | count | lengths | fold origin |
|---|---|---|---|---|
| krasG12D_01–10 | scaffold-docked + ProteinMPNN resurfacing | 10 | 54–103 aa | de novo scaffolds published for *other* pMHC targets |
| krasG12D_11–20 | BoltzGen all-atom generation against this target | 10 | 61–109 aa | generated, no parent |

**The order of the twenty is not a confidence ranking, and we make no affinity claim
for any of them.** The reason is given under *Why there is no predicted affinity* below.
We report it rather than omit it, because an unverifiable self-reported metric is worth
less than an honest statement that we do not have one.

## Target preparation

Structures 7OW6 (G12D, `VVVGADGVGK`) and 7OW5 (wild type, `VVVGAGGVGK`) from the PDB.
Three states were built on one shared background that differ **only** in the peptide:

- `g12d` — heavy chain (275 aa) + β2-microglobulin (99 aa) + `VVVGADGVGK`
- `wt` — the same, with `VVVGAGGVGK`
- `apo` — the same, with no peptide

An automated check asserts the three states are identical outside the peptide and that
residue 6 carries the Asp side chain in `g12d` and not in `wt`. The epitope focus is
peptide position 6 — the G12D aspartate — plus the 18 heavy-chain residues within 8 Å
of it, which were computed rather than assumed.

## Route 1 — scaffold-docked (krasG12D_01–10)

Backbones are taken from the de novo scaffold library published with Liu et al.,
*Science* 389:386–391 (2025). Those scaffolds were designed for **different** pMHC
targets; none was a binder of KRAS G12D / HLA-A\*11:01. Each was docked over the
peptide groove and then **every residue of the binder was rewritten** by ProteinMPNN
against this target. Identity to the parent scaffold is computed per design, recorded
in the accompanying detail file, and anything above the cut was dropped.

We declare this provenance rather than obscure it. The sequence is ours; the fold is
reused, and the judges should weigh that themselves.

## Route 2 — BoltzGen (krasG12D_11–20)

BoltzGen 1 (`protein-anything` protocol, both released checkpoints, 120 designs,
design + inverse-folding steps) generating a 60–110 residue chain against the full
three-chain complex, with the specification stating the requirement up front rather
than filtering for it afterwards:

```
binding:      peptide chain, residues 4–8   (the G12D aspartate at 6 and its neighbours)
not_binding:  heavy chain, α3 domain
not_binding:  β2-microglobulin, whole chain
```

**An earlier run of this route failed completely and we report it.** β2-microglobulin
was omitted from the target, leaving a β2m-shaped gap in the complex, and all 172
harvested designs came back as β2-microglobulin — the model completed the obvious hole
instead of gripping the peptide. Including β2m fixed the gap but not the tendency: of
the 120 designs in the run used here, **59 are still substantially a copy of a chain we
supplied, and 18 are antibody or T-cell-receptor variable domains**. Both groups were
removed before any further consideration.

Echo detection does not use a flat identity threshold, because a short design slid
along a 275-residue chain scores ~0.30 identity by chance alone. Each design is scored
against its own shuffled self — same length, same composition, order destroyed — and
excluded on the excess over that null. Every one of the ten submitted sits at or below
+0.12 excess, most at or below chance.

TCR and antibody frameworks are excluded on principle, not on score: a TCR is the
natural binder of a pMHC, and submitting one would be a de novo rule problem regardless
of what it measured.

## Selection of the final twenty

Applied identically to both halves:

1. **Expression and manufacturability**, using Adaptyv's own published thresholds:
   instability index < 40, GRAVY < 0.4, no odd cysteine count, no hydrophobic run of
   six or more, plus flags for NG/NS/NT deamidation motifs and runs of three or more
   K/R. Anything failing a redesign-level rule was excluded outright. The submitted
   twenty run instability 8.3–38.8 and GRAVY −0.84 to +0.05, with no odd cysteine and
   no hydrophobic run of six. On top of that, our own charge-balance, hydrophobic
   fraction and proline/glycine fold-plausibility screen, used since challenge 1.
2. **Diversity** — no two designs in the final twenty exceed 0.60 identity, measured
   over the best ungapped register on the shorter sequence, which over-estimates
   similarity and therefore rejects more than an alignment would. Highest observed pair
   in the submitted set: **0.57**. The scaffold-docked half draws from ten distinct
   parent backbones.
3. **De novo standing** — `adaptyvbio/binder-prescreen`, full screen, all arms except
   the USPTO patent corpus.

## Pre-screen result

All twenty return `pass`. Prior-art similarity, design similarity and target similarity
are **0.0** for every design; ProteinTyper categorises all twenty as `miniprotein`,
agreeing with the declared class.

For context on how strict that is: run on all 120 raw BoltzGen designs, the same screen
returned 68 `pass`, 47 `near_known_sequence` and 5 `known_sequence`, and 11 of its
*passes* still carried prior-art similarity above 0.5. Our own filter is stricter than
the gate, and none of those 11 is in this submission.

## Ranking: ESM2 pseudo-log-likelihood

The twenty are ordered by ESM2 pseudo-log-likelihood, reported in Adaptyv's QC guide as
the best single predictor of binding at ROC AUC 0.72, ahead of ipTM at 0.64 and PAE at
0.65. Their guide writes the threshold as "PLL > 0", which cannot be a raw
log-likelihood, since those are always negative, and does not say how they normalise.
We therefore report our own normalisation and name it: the **excess over the same
sequence shuffled**, averaged over three shuffles, which is positive when a sequence is
more plausible than its own amino-acid composition would give by chance. All twenty
score positive, from +0.043 to +0.266.

An AUC of 0.72 is a weak signal. It is, however, the only criterion in this submission
with a published number attached, and it is a sequence metric the organisers can
recompute themselves.

## Why there is no predicted affinity

We built a three-state AlphaFold screen (ColabDesign binder protocol, two models per
state) and ran it twice, 300 predictions in total. Before trusting it we ran two
controls through the identical pipeline:

| | i-pTM (g12d) |
|---|---|
| JDI T-cell receptor from 7OW6 — a confirmed binder of this exact pMHC | 0.158 |
| A shuffled sequence — same length and composition, cannot bind | 0.134 |
| **Separation** | **+0.024** |
| Run-to-run disagreement between two models of the same state | 0.022 |

The separation is the size of the noise. The positive control's own preference for the
G12D peptide over wild type came out at **−0.001**: a protein that demonstrably
discriminates registered no discrimination. The second scoring method, ipSAE, failed on
all 300 predictions.

We then went back to the 156 stored PAE matrices and computed **pae_interaction**, the
metric Bennett et al. filter on at < 10 and Adaptyv list at < 12 standard. It gives the
same verdict, slightly worse:

| | pae_interaction (G12D state) |
|---|---|
| JDI TCR, confirmed binder | 27.62 |
| Shuffled sequence | 27.43 |
| **Shuffle minus real binder** | **−0.19** (the real binder is marginally *worse*) |
| Worst within-state noise | 2.16 |

Both sit far above the < 12 threshold, and the known binder does not beat the shuffle.
Two independent interface metrics therefore agree. The most likely cause is that these
predictions were run **without an AF2 initial guess** — re-predicting the complex with
the designed complex as the structural starting point — which is the protocol the
pae_interaction < 10 threshold was calibrated on. We did not have time to rebuild the
screen that way, and we would rather say so than publish a ranking from an instrument
we have shown cannot read.

We therefore treat this screen as a **failed measurement**, not as evidence about our
designs in either direction. An earlier internal conclusion that none of the designs
bound has been withdrawn as unmeasured. No number from it is used to order, filter or
describe this submission.

## What we did not do

- No affinity prediction of any kind is behind these twenty.
- No AF2 initial-guess screen was run. This is the single clearest gap in the work and
  the most likely reason our structure-based numbers are uninformative.
- No experimental data. Nothing here has been expressed or tested.
- ipSAE was never made to work; the Boltz-2 ipSAE pipeline used by the organisers was
  not run on our side, so our numbers should not be compared against Boltz-2 thresholds.
- Structural novelty (ProteinTyper level) was not measured by us. The scaffold-docked
  half reuses a published fold and may therefore score lower on structural novelty than
  the BoltzGen half; splitting the submission ten and ten is a deliberate hedge against
  that, not a claim that either half is better.
- The USPTO patent arm of the pre-screen was skipped (1.7 GB, ~42 min). Every other arm
  was run.

## Reproducibility

Code: `github.com/wwwtalalahmed98-collab/ph-switch-binders`. Design sequences were held
out of the public repository until the submission window closed on 9 November 2026
11:59 UTC, and were published with this methods file once it did. Every generation,
screening and selection step ran in a disposable cloud VM from a pinned commit; the
scripts that built the three states, generated the designs, triaged them, ran the
pre-screen and assembled this set are all in `pmhc/` and `kaggle/`.
