# Does the field's scoring actually work?

Eight read-only scripts that re-measure three public datasets of designed protein binders
against their **real laboratory outcomes**.

They exist because of a failure in this project's own work. The challenge 5 AlphaFold
screen was given two controls — a crystallographically confirmed binder and a shuffled
sequence — and could not tell them apart. That raised the question these scripts answer:
was it our pipeline, or are the metrics themselves weak?

## The data

Nothing here is committed; each script downloads what it needs. Combined, about 140 MB.

| dataset | designs | get it from |
|---|---|---|
| Proteinbase (Adaptyv Bio), ODC-BY | 4,280 | `https://storage.proteinbase.com/proteinbase_latest.csv` |
| Claude binder release (Anthropic), CC-BY-4.0 | 1,440 | `huggingface.co/datasets/Anthropic/claude-protein-binder-design` → `data/tables/design_summary.csv` |
| DTU meta-analysis, CC-BY-4.0 | 3,676 | `zenodo.org/records/15722219` → `final_dataset.csv` |

## The scripts

Run each with `python -I <script> <path-to-csv>`. Pure Python plus numpy, scipy and
scikit-learn; no GPU, minutes not hours.

| script | question it answers |
|---|---|
| `survey_pb.py` | what is actually in the Proteinbase export? |
| `auc_pb.py` | which computational metric predicts the measured binding outcome? |
| `conditioned.py` | how much of that is really predicting *expression* rather than binding? |
| `affinity.py` | does anything predict **how tightly** a binder binds? |
| `anth_survey.py` | ten structure predictors, identical designs, dual-CRO labels — which wins? |
| `dtu_check.py` | does the DTU conclusion survive its own target imbalance? |
| `pooled_vs_within.py` | are published pooled benchmarks inflated? |
| `feasibility.py` | can a trained model beat the best single score, held out by target? |
| `kinetics.py` | affinity is on-rate over off-rate — is either half predictable? |

## What they found

- **ipSAE scores 0.56, 0.74 and 0.81** on the three datasets. Same metric, same question,
  a 0.25 AUC spread. Thresholds do not transfer between pipelines.
- **Chai-1 (0.744) and OpenFold3 (0.741) beat AlphaFold3 (0.733)**; Boltz-2 is near the
  bottom at 0.692, on 1,320 identical designs with two labs agreeing.
- **Nothing predicts affinity.** Best single correlation −0.26; a model under
  leave-one-target-out does worse. Off-rate drives affinity (+0.70) but predictors track
  it at only −0.27, and they disagree about the *sign* on on-rate.
- **Expression is a hidden confound.** Boltz-2 complex pLDDT predicts whether a design can
  be made at all at 0.76, but only 0.62 for whether it binds once made. Filter on it and
  you are mostly selecting manufacturable proteins.
- **One hypothesis died.** Pooled benchmarks are *not* materially inflated relative to
  within-target ranking: the drop is 0.036 and 0.003.

Full write-up, with every claim marked as verified-here or relayed-from-a-paper:
[`../docs/RESEARCH-LANDSCAPE.md`](../docs/RESEARCH-LANDSCAPE.md).

## Caveats

Single-pass analyses. No confidence intervals, no correction for testing ~30 metrics,
replicates aggregated by majority vote rather than averaged, and target / design-method /
designer-skill are entangled throughout. Indicative, not publishable as they stand.

## Attribution

This work used **Proteinbase** by Adaptyv Bio under the ODC-BY licence; the **Claude
protein binder design** release by Anthropic under CC-BY-4.0; and the dataset for
*Predicting Experimental Success in De Novo Binder Design* (Zenodo 10.5281/zenodo.15722219)
under CC-BY-4.0.
