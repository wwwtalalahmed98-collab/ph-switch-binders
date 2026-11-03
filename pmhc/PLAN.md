# Challenge 5 - KRAS G12D pMHC. Plan.

Written 2026-11-03. Spec read from the official challenge page, not from memory.

## The target, exactly

| | |
|---|---|
| Peptide | `VVVGADGVGK` - KRAS residues 7-16 of P01116, Gly12 to Asp at **peptide position 6** |
| Presented by | HLA-A\*11:01 (allele A\*11:01:01:01, HLA00043) |
| Complex | HLA-A\*11:01 heavy chain + beta-2-microglobulin + peptide |
| Structure | PDB **7OW6**, chains A-C, TCR removed |
| Reference | 7OW5 (wild type, TCR removed); 7OW4 (G12D, no coordinates for peptide positions 5-6) |
| Recommended epitope | the peptide and the groove around it |

Counter-targets the design must **not** bind:

1. Wild-type complex - `VVVGAGGVGK` on the same HLA-A\*11:01
2. Peptide-free HLA-A\*11:01 + beta-2-microglobulin

Binding to the His/Avi tags or our own C-terminal Twin-Strep tag does not count.

## Submission rules

- Track 3: **up to 20 designs**, length 10-250 aa, single-chain protein / nanobody / scFv / Fab
- **De novo only.** Modifying an existing binder is not allowed.
- Fields: `name`, `sequence`, `molecule_class`; rows ranked, best first
- Encouraged: ipTM, ipSAE, self-consistency, structure models, methods text
- Deadline **Sunday 8 November 2026, 23:59 AoE**. One submission per 24 h, so
  **everything final Saturday 7 November.**

## What makes this different from the four challenges before it

Affinity is the *third* tie-break. The first two ranking steps are both ratios:

1. KD(wild type) / KD(G12D), largest first
2. KD(peptide-free) / KD(G12D), largest first

**The score is specificity, not strength.** A design that grips hard and binds all
three states ranks below a weak design that only binds the mutant.

The discrimination is a single Gly to Asp at peptide position 6, sitting in the
groove. Everything else about the three states is identical. A binder that reads
the HLA helices or beta-2-microglobulin instead of the peptide will bind all three
and score nothing.

On TfR1 we bolted the off-target check on at the end and discovered 3 of 47 designs
blocked the receptor's day job. Here the counter-screen cannot be an afterthought -
it is the objective.

## Plan

### Phase 0 - build the three states (today, no GPU, ~2 h)

Crop 7OW6 to heavy chain + b2m + peptide. Build the wild-type state from 7OW5 and
the peptide-free state by deleting the peptide chain. Three crops, same framing.

**Verification, from the TfR1 lesson:** every crop is identified downstream by
reproducing its sequence residue for residue, never by residue numbering and never
by "the bigger chain". `save_pdb` renumbers from 1 and that silently destroyed two
whole rounds of measurement on TfR1.

### Phase 1 - choose the hotspots (no GPU)

The entire specificity story is "we contact Asp6". Hotspots = the peptide Asp at
position 6 plus the groove residues flanking it. Objective uses the corrected
contact term from TfR1: `i_con = dict(cutoff=8.0, num=1, num_pos=<small>)`.

The default cutoff of 21.69 admits contacts out to 21.4 A and is why round 1 on
TfR1 aimed at nothing. That is fixed and the script refuses to run if the setting
did not take.

### Phase 2 - generate against the mutant state (GPU, ~6-8 h)

AfDesign hallucination against the G12D crop only. **Save the PAE matrix.** We have
failed to compute ipSAE in three consecutive challenges for exactly one reason: we
never saved the file it needs. Fixed at the start this time, not at the end.

### Phase 3 - rewrite the sequences (CPU, free)

ProteinMPNN with the backbone fixed, then the composition screen - now three-sided:
charge, hydrophobicity/Kyte-Doolittle, identical-residue runs, **and proline/glycine**.
The proline gate is new as of 29 Oct and it cut a reported 12 good candidates to 1.

### Phase 4 - the counter-screen (GPU, ~6 h)

Fold every surviving design against **all three states** and score the margin.
Rank on predicted margin, not on confidence. This triples the GPU cost per design,
so we will carry fewer designs further - that is the correct trade here, because
margin is what is scored.

### Phase 5 - novelty check and submit (CPU, Sat 7 Nov)

Run `adaptyvbio/binder-prescreen` on all 20 in its no-download mode, which checks
our sequences against a curated set of 144 known KRAS-G12D/HLA-A\*11:01 TCR
references. De novo is a hard requirement and we have never checked it in four
challenges. Then write the methods and submit.

## Honest assessment

- **We have never produced a verified binder.** Four challenges, nothing confirmed.
  TfR1's yield was 1 usable sequence out of 480 written.
- **This target is harder than TfR1, not easier.** Discriminating one buried residue
  is a TCR-mimic problem; the published work on it uses a dedicated pipeline.
- **The ranking rests on a quantity AlphaFold is not known to predict well.** A
  predicted specificity margin is not a measured one. This is the weakest joint in
  the plan and it should be stated that plainly in the methods text.
- **Five days, and the previous four rounds each lost most of their time to one
  silent bug.** The guards above exist because of specific bugs that happened; they
  do not rule out a new one.
- Realistic goal: 20 designs with an honest counter-screen and methods that state
  the limitations. Selection for wet lab weighs method novelty and diversity, which
  a disciplined three-state screen and a published failure record can actually speak
  to, whereas a claim of high affinity cannot.
