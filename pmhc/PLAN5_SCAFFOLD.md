# Challenge 5, revised plan: stop inventing backbones, reuse validated ones

Written 6 November 2026, after three generation runs produced zero usable designs.

## What went wrong, in one paragraph

We tried to invent binder backbones from scratch by back-propagating through
AlphaFold2 (AfDesign "hallucination"). Three runs, ~16 of 30 GPU hours, 62 candidate
structures, **zero** that passed our own buildability screen. The diagnosis was clean:
ProteinMPNN, asked to write sequences for the 15 best backbones, returned
glycine/proline soup — 720 sequences, zero survivors. ProteinMPNN places Gly and Pro
where local geometry demands them, so Gly/Pro soup is a readout meaning *this backbone
is not a foldable protein*. We were polishing noise.

## What we found instead

Liu *et al.*, **Science 389:386–391 (2025)**, "Design of high-specificity binders for
peptide–MHC-I complexes" — the Baker lab solved exactly this problem: de novo binders
to pMHC-I that discriminate a mutated peptide residue. Eleven targets, nanomolar
affinities, specificity confirmed by yeast display, CAR activation and a crystal
structure. They also state plainly why our route was so expensive:

> "only a small fraction of trajectories yield backbones that interact extensively
> with the peptide while largely avoiding the MHC"

Their answer: don't generate backbones from scratch. Reuse a set of scaffolds already
shown to arc over an MHC-I groove and dock them onto the new target. **They published
the scaffolds and the code** (Zenodo 10.5281/zenodo.15169815, GitHub
`77forest/pMHCI_binder_design`).

### What is in that repository

| Asset | What it gives us |
|---|---|
| `scaffolds/` — 152 PDBs | Binders 54–118 aa, **zero cysteines**, each paired with the pMHC it was designed against |
| `software/scaffold_alignment/align_chainB.py` | Docking by MHC superposition (needs PyMOL; we reimplemented in Biopython) |
| `contact_filter/` | Keep poses that read the peptide and avoid the MHC |
| `mpnn_spec_filter/` | The specificity score — see below |
| `pMHC_fold/` | Motmaen *et al.* fine-tuned AF2 for pMHC (PNAS 120:e2216697120) |
| `software/mpnn/` | Full ProteinMPNN with `v_48_020` weights included |

Their scaffold chain B is trimmed to 183–191 residues (MHC platform + peptide). Our
Phase 0 crop sweep independently arrived at **189** (179 platform + 10 peptide). The
same answer by a different route is reassuring about the earlier work.

## The one idea worth more than the scaffolds

The paper's specificity filter:

> "We then used ProteinMPNN to predict the sequence of the target peptide in the
> context of the designed binder and selected those designs for which the confidence
> of the on-target peptide sequence was higher than off-target sequences"

Turn it around for our target: hold the binder and the MHC fixed, delete the peptide
sequence, and ask ProteinMPNN what peptide belongs in that pocket. A design that truly
reads the G12D mutation will score **Asp** above **Gly** at peptide position 6. A design
that merely sits on the groove will not care.

This matters disproportionately because the competition ranks Track 3 on
**KD(WT)/KD(G12D)** first and affinity only as a tie-break. The screen is CPU-only,
free, runs in minutes, and measures the quantity we are actually scored on. We have had
no measure of specificity at all until now.

## Pipeline

| # | Step | Where | Cost | Status |
|---|---|---|---|---|
| 1 | Dock 152 scaffolds by MHC superposition | local CPU | minutes | **done** |
| 2 | Filter poses: peptide contact, Asp6 contact, no backbone clash | local CPU | seconds | **done** (same script) |
| 3 | ProteinMPNN redesign of the binder, many sequences per pose | CPU | ~1 h | to build |
| 4 | Composition + buildability screen (the TfR1 screen, Pro/Gly gate included) | local CPU | seconds | **reusable** |
| 5 | MPNN peptide-recovery specificity screen (D vs G at position 6) | CPU | minutes | to build |
| 6 | AF2 on all three states (G12D / WT / apo) in one pass | Kaggle GPU | 3–4 h | adapt built screen |
| 7 | Rank by specificity margin, then confidence; ipSAE; methods; submit | local | 1 h | converter **built** |

Steps 1, 2, 4 and the ipSAE converter already exist. Step 6 is the three-state
counter-screen we built and unit-tested but never ran.

### Merging steps 6 and 7 of the old plan

The previous plan validated designs with AF2 and *then* counter-screened them, two GPU
passes. One pass that predicts G12D, WT and apo for every candidate gives both the
confidence metrics and the specificity margin. That halves the GPU bill, which matters:
~14 of 30 hours remain this week.

## Compute and time

- **Deadline**: Sunday 8 November 23:59 AoE. One submission per 24 h, so the first
  submission must go in **Saturday 7 November** to leave a second attempt.
- **GPU remaining**: ~14 h of the 30 h weekly Kaggle quota.
- **GPU needed**: step 6 only, ~3–4 h on two GPUs. Everything else is CPU.

The critical path is deliberately short on GPU, because GPU is where all three previous
runs died.

## The rule risk, stated plainly

Track 3 says: **"De novo only. Modifying an existing binder is not allowed."**

These scaffolds are themselves de novo designed proteins, and they were designed against
*different* pMHC targets, not ours. Starting from a de novo scaffold library and
redesigning the sequence is standard practice in the field — it is what RFdiffusion's
fold-conditioning and the Cao *et al.* minibinder libraries do. But a strict reading of
"modifying an existing binder" could cover it, and I cannot resolve that from the rules
text.

Mitigation, which costs us nothing:

1. ProteinMPNN redesigns **100% of the binder sequence**, not just the interface.
2. We **measure and record sequence identity to the parent scaffold** for every design,
   and report it in the methods text.
3. We drop any design above a stated identity threshold to its parent.
4. The methods text names the scaffold set and the paper explicitly. If the organisers
   judge it ineligible, they can see exactly why from the submission itself.

Hiding the provenance would be both wrong and pointless, since the sequences are public.

## What this plan does not claim

- No design has been built yet. Step 1 produced docked *poses* carrying their original
  scaffold sequences; nothing has been redesigned or predicted.
- The specificity screen is the published method but we have not yet shown it
  discriminates anything on our target. It could come back flat.
- AF2 confidence metrics are not predictive of affinity. The BindCraft authors say so
  directly, and the Liu paper relies on experiment, not metrics, for affinity.
- Our WT counter-state is a **model**: we strip the Asp6 side chain and rename the
  residue to Gly. The real 7OW5 structure exists and we are not using it as the
  counter-state background, so the WT/G12D comparison shares one backbone by
  construction. That makes the comparison cleaner but less realistic.
- The apo state is the easiest of the three to predict well, so a good apo margin is
  weaker evidence than a good WT margin.
