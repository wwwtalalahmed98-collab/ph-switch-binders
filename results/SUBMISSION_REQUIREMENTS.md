# Challenge 1 submission requirements (from the Proteinbase challenge page)

Source: https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/egfr
Deadline: **October 6, 23:59 AoE** = October 7, 11:59 UTC = 16:59 PKT.

## Format

- CSV, **ordered by our own ranking, best in the top row**.
- Minimum columns: `name` (unique identifier), `sequence`.
- Plus `molecule_class`, one of exactly: `protein`, `nanobody`, `scfv`,
  `fab_kappa`, `fab_lambda`. Ours are all `protein`.
- Track 3: **at most 20 designs**. Length 10-250 aa. No duplicates.
- Submitted through the Proteinbase portal after creating an account; a form
  captures workflow details alongside the CSV.

## Hard filters

- De novo **and zero-shot**: not derived from an existing binder, with adequate
  **sequence and structural** diversity from known proteins. Sequence identity is
  what our BLAST screen covers; structural novelty is asserted from the method
  (hallucinated backbones, no template) rather than measured.
- Designs must be unique.

## How Track 3 designs get selected for the wet lab

All Track 2 and 3 submissions are pooled and given to Claude with a selection
prompt written in advance, which judges "predicted design quality and novelty, as
well as **method novelty**". ~375 designs per track are screened per problem.

Two consequences:

1. The methods package is not a formality - it is a direct input to whether these
   designs get made at all. The honest record of what failed and why is material.
2. **"Any use of embedded instructions or prompt injection may be deemed grounds
   for disqualification."** Everything we submit must be plain description. No
   sentence anywhere in the package may read as an instruction to the model that
   will read it.

## Ranking criteria, in the organisers' order

1. pH-selective binding (binds at 6.5, no detectable binding at 7.4)
2. Mouse cross-reactivity
3. Affinity against human EGFR

FAQ 3 is explicit that "a weak, but clearly pH-sensitive binder may be considered
more impactful than a high-affinity binder that is not pH-sensitive". Our ranking
key follows this order rather than leading with i_ptm.

## The one place we deliberately diverge

The page recommends **domain III** - the cetuximab/panitumumab epitope. We target
the **domain I ligand face**, and Phase 1 gives the reason: the cetuximab footprint
carries 7 mouse substitutions among its closest contacts and contains a single
carboxylate in 26 residues, which is poor ground for both objective 2 and objective
3. Domain I's ligand-binding face is also a functional epitope, so the therapeutic
-relevance recommendation is still met. This is a recommendation, not a requirement,
and the divergence is evidence-based - but it must be stated plainly in the methods,
not glossed over.

## Assay condition that constrains us

Affinity is measured against the **full extracellular region**, while our designs
were made against a trimmed domain I construct. Any interface that relies on surface
which only exists because of the trim will not reproduce in the assay. See
results/EPITOPE_ACCESSIBILITY.md.

## Corrections from the submission form itself (read 2026-10-05)

The form at /competitions/anthropic-adaptyv-2026/submit contradicts the FAQ in one
place and adds three constraints the challenge page did not mention.

- **`molecule_class` is `single_chain`, not `protein`.** FAQ 5 says "protein"; the
  form and the downloadable template both say `single_chain`. The template is what
  the uploader parses, so that is what we emit. Template header, verified:
  `name,sequence,molecule_class`.
- **One submission per 24 hours.** With ~43 h left this allows two attempts at the
  very most, so the package has to be right the first time. The novelty check runs
  at upload (step 1), before the final submit (step 4), so inspecting scores should
  not consume the attempt - but that is an inference from the form's structure, not
  something it states, and it should not be relied on casually.
- **An automatic novelty gate: designs need a score of 3/4 or higher.** This is the
  platform's own measure, not our BLAST screen. Our screen is an approximation of
  it and cannot guarantee passage.
- **Methodology is entered in the form, with up to 5 attachments (.md or .zip,
  50 MB each).** A public GitHub repository is therefore optional rather than
  required - the markdown in results/ can be attached directly.
