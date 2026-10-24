# De novo conditional protein binders

**Anthropic × Adaptyv Protein Design Competition 2026 — Track 3**

[![Licence: MIT](https://img.shields.io/badge/licence-MIT-1d9e75)](LICENSE)
[![Competition](https://img.shields.io/badge/competition-Anthropic%20%C3%97%20Adaptyv%202026-534ab7)](https://proteinbase.com/competitions/anthropic-adaptyv-2026)
[![Track](https://img.shields.io/badge/track-3%20(open)-888780)](https://proteinbase.com/competitions/anthropic-adaptyv-2026)
[![Designs](https://img.shields.io/badge/designs%20submitted-25-185fa5)](#results)

Designing small proteins that bind a chosen human target **under a condition** — a pH
window, a second species, or a mechanism — and submitting them for wet-lab testing.

[Results](#results) · [What is here](#what-is-here) · [Method](#method) ·
[Findings worth reusing](#findings-worth-reusing) · [Reproducing](#reproducing) ·
[Citation](#citation)

![Depth does not separate an activating peptide from a blocking one](assets/depth-validation.svg)

---

## Results

| | Target | The condition | Status |
|---|---|---|---|
| **Challenge 1** | EGFR ectodomain (6ARU) | bind at pH 6.5, not 7.4 | 5 designs submitted |
| **Challenge 2** | TNF-α trimer (1TNF) | bind at pH 7.4 not 6.0, **and** cross-react with mouse | 20 designs submitted |
| **Challenge 3** | GLP-1R (6X18) | bind **and activate** the receptor | in progress, closes 25 Oct |

No wet-lab results exist yet. Every number in this repository is a prediction or a
geometric measurement, and the distinction is kept explicit throughout.

## What is here

```
glp1r/        Challenge 3: target scoping, depth analysis, selection rule
tnf/          Challenge 2: epitope selection, switch placement, submission assembly
scripts/      conservation, epitope finding, novelty screen, pipeline stages
kaggle/       the notebook cells that ran on GPU, one per run, with their failure notes
viz/          builders for the interactive structure viewers
results/      reports and measurements, including the ones that went wrong
docs/         PROJECT-STORY.md — a plain-English record of every decision and reversal
assets/       figures
```

Design sequences are withheld until each submission window closes, by pattern in
`.gitignore`. Structures and reference data are downloaded by the scripts rather than
committed.

## Method

1. **Choose the epitope before designing.** Cross-species conservation and surface
   chemistry decide the site, not the optimiser.
2. **Generate** binders with AfDesign / ColabDesign hallucination against a cropped target.
3. **Score by independent re-prediction.** Every reported confidence comes from a fresh
   prediction of the discrete sequence at 3 recycles — never the optimiser's own number,
   which is optimistic by construction.
4. **Validate the selection rule before applying it**, against a case where the answer is
   already known.
5. **Select, declare, submit.** The rule is written down before the designs exist, and any
   later change is recorded as a dated amendment.

## Findings worth reusing

These are the results that cost the most to learn and are the most likely to save someone
else a run. Most of them are negative.

**A geometric proxy is worthless until it is tested against a known answer.** In Challenge
2 the switch designs were selected on histidine-to-cation distance. Measured afterwards,
that variable correlated **+0.24** with the effect it was supposed to predict — the wrong
sign. Burial on binding correlated **−0.58** and the set was rebuilt on it.

**Depth cannot distinguish a receptor agonist from an antagonist.** Nature supplies the
control: GLP-1(7-36) activates GLP-1R and the same peptide missing two residues is the
shape of a drug that blocks it. AlphaFold2 placed them 0.8 Å apart, and a scrambled
sequence reached *deeper than the real hormone* (figure above). What did separate sense
from nonsense was **polarity** — which end of the peptide points down the channel.

**A cropped target invents pockets, and the optimiser will find them.** Cutting the
receptor down to an affordable size left an opening where the helical bundle should be. At
that depth our crop held **4 atoms against the intact receptor's 211**. The three
highest-confidence threading designs all dived into it, including one scoring i_ptm 0.800
— the best number this project has produced, and an artifact.

**A model's confidence can anti-correlate with being right.** Across 11 designs aimed into
the GLP-1R helical bundle, predicted-interface confidence correlated **−0.963** with
touching the residues the natural hormone actually uses. Every design above 0.66 touched at
most 3 of them; the only two that landed on the hormone's site scored the two *lowest*
confidences. Selecting on confidence would have picked six wrong designs and rejected the
one that works — see [glp1r/results/THREADER_FAILURE.md](glp1r/results/THREADER_FAILURE.md).

**Check residue numbering against the structure, not against your own map.** A target
chain renumbered 1–322 with gaps is not indexable by ordinal. The same class of error
appeared three times here; the one that survived longest "passed" a check that compared a
map against itself.

**Two histidines are not twice one.** Henderson–Hasselbalch says one histidine leaves 50%
of molecules competent at pH 6.0 and two leaves 25% — but only if both are doing work.
Every pair measured had a partner buried by ≤0.10, so the effective count was one and the
cost at pH 7.4 was two.

## Reproducing

The GPU work ran on Kaggle (T4 ×2). Each file in `kaggle/` is one notebook cell, carrying
its parameters and a comment explaining what the previous version got wrong. The analysis
scripts run locally and need only `biopython` and `numpy`.

```bash
python glp1r/scope_6x18.py data/glp1r/6X18.pdb      # depth profile of the native hormone
python viz/build_tnf_viewer.py                      # interactive viewer for the designs
```

`docs/PROJECT-STORY.md` is the narrative record — what was chosen, what was rejected, and
which assumptions turned out to be wrong.

## Citation

```bibtex
@misc{ahmed2026conditionalbinders,
  author = {Ahmed, Talal},
  title  = {De novo conditional protein binders for the
            Anthropic x Adaptyv Protein Design Competition 2026},
  year   = {2026},
  url    = {https://github.com/wwwtalalahmed98-collab/ph-switch-binders}
}
```

## Licence

MIT — see [LICENSE](LICENSE).
