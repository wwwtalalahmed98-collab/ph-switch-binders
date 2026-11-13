# The de novo binder design landscape: data, tools, validation, gaps

Research dossier, compiled 10 November 2026. **No project is proposed here.** This is the
evidence base for choosing one later.

Everything below is marked **[verified]** when I downloaded it, ran it, or read the primary
source myself in this session, and **[reported]** when I am relaying a claim from a paper,
blog or vendor page without independent checking. That distinction matters: a lot of the
headline numbers in this field come from the people selling the method.


> **Correction, second pass (10 Nov 2026).** Section 5 of this dossier originally reported
> ipSAE at AUC 0.557 and implied the metric is near-useless. That was measured on
> Proteinbase's single-run Boltz-2 pipeline. On two other datasets, computed properly, the
> same metric reaches **0.69–0.74** and **0.79–0.81**. The metric is not the problem; **how
> it is computed, and how clean the labels are, moves it by 0.25 AUC**. See section 8.

---

## 1. The one-paragraph summary of where the field stands

In April 2026 David Baker's group published a review in *Nature* stating that "the
long-standing challenges of designing new protein structures, assemblies and protein
binders are close to being solved", and that "the main question in the protein binder
design area is no longer how to design binders, but **what to design binders to and
for**." **[reported]**

That is the strategic fact. Generation is commoditised. What is *not* solved, and what the
same literature says repeatedly, is everything downstream of generation: deciding which
designs to make, predicting how tightly they will bind, making them specific, and having
standard data to measure any of it against.

---

## 2. Part A — the data: who has produced validated binder measurements

### 2.1 Proteinbase (Adaptyv Bio) — the single-protocol corpus

**[verified]** Downloaded `https://storage.proteinbase.com/proteinbase_latest.csv`, 52 MB,
and counted it directly.

| | |
|---|---|
| proteins | **4,280** |
| with experimental measurements | 4,278 |
| with computational metrics | 4,251 |
| **with both** | **4,249** |
| experimental measurement rows | 64,026 |
| computational metric rows | 141,352 |
| distinct targets | 40 |
| distinct design methods | 140 |
| designs with a real K<sub>D</sub> | **1,062** (2,661 replicate measurements) |
| K<sub>D</sub> range | 57 pM to 32 µM |
| licence | ODC-BY (free use with attribution) |
| also at | `huggingface.co/datasets/adaptyvbio/proteinbase`, parquet, snapshot-pinnable |

Targets are concentrated: Nipah glycoprotein G, PD-L1, SARS-CoV-2 RBD, EGFR, VEGF110
carry most of the weight. Design methods include BoltzGen (519), BindCraft (174),
RFdiffusion (162).

**Why it matters:** one laboratory, one protocol, one set of definitions, across 40
targets. That is rare. Most datasets in this field are pooled from papers that each used
different assays and different binding thresholds.

**[verified]** HuggingFace metadata: 395 downloads, 3 likes, last updated 1 Nov 2026. It
is new and barely used yet.

### 2.2 The Anthropic / Claude binder release — the validation gold standard

**[verified]** Read the dataset documentation at
`huggingface.co/datasets/Anthropic/claude-protein-binder-design`.

| | |
|---|---|
| designs | **1,440** (50–120 residues), 16 targets |
| with reliable wet-lab measurements | 1,320 (mature GDF-8 excluded, antigen aggregated) |
| binders | **354** on 14 of 15 targets; hit rate 27% |
| validated by | **two independent CROs**: Adaptyv Bio and Twist Bioscience |
| licence | CC BY 4.0 (data), MIT (scripts) |
| files | ~100,000 — raw response curves, replicate reads, standardized re-fits, sensorgram images, ten co-fold predictors per design, PAE matrices in a companion release |
| **[verified]** uptake | 35,946 downloads, 240 likes |

**The number that makes this dataset unique:** of 1,235 designs measured by both CROs,
253 bound at both, 846 at neither, **69 at Adaptyv only and 67 at Twist only — 89.0%
agreement, Cohen's κ = 0.71**. **[reported, from the release docs]**

That 11% disagreement is a *direct, empirical measurement of label noise* in this field.
It sets a ceiling on how well any predictor can ever appear to perform against this kind
of label, and nobody else has published the equivalent. It also carries counter-screens:
mouse and cynomolgus orthologs, and an unrelated off-target control (CLEC12A).

### 2.3 The DTU meta-analysis — the current benchmark of record

**[verified]** Read the full text of *"Predicting Experimental Success in De Novo Binder
Design: A Meta-Analysis of 3,766 Experimentally Characterised Binders"*, bioRxiv
2025.08.14.670059.

- **3,766** binders across **15 targets** (including two pMHC complexes), pooled from
  multiple published sources.
- Only **436 (11.6%)** were binders — heavily imbalanced.
- Re-predicted every complex with AF2 initial guess, ColabFold, Boltz-1x and AF3; over
  **200 features** per design.
- **Finding: AF3 ipSAE_min was the best predictor across all targets**, a 1.4× increase in
  average precision over the ipAE score RFdiffusion filters on. ipTM performed well at low
  recall. Per-target precision at best-F1 ranged **0.1 to 1.0** — i.e. wildly target-dependent.
- Data: Zenodo `10.5281/zenodo.15722219`. Code: `github.com/DigBioLab/de_novo_binder_scoring`.

**Their own stated limitations, quoted because they define the gap precisely:**

> "For most targets only limited number of tested designs were available with only a
> handful of true binders... heterogeneous assay types and binding thresholds introduced
> label noise, and **affinity data was only available for a small subset and thus not
> incorporated**. More standardized, publicly available datasets linking binder structure
> to binding affinity are needed."

Proteinbase and the Anthropic release are, between them, exactly the datasets they are
asking for — and both postdate or sit outside their analysis.

### 2.4 The Adaptyv EGFR competition rounds

**[reported]** Rounds 1 and 2, data at `github.com/adaptyvbio/egfr_competition_1` and
`_2`. Over **1,800 designs submitted, 601 expressed and characterised by BLI, 60 bound**.
Hit rate rose from ~3% (round 1) to 13% (round 2). Community paper: Cotet et al., bioRxiv
2025.04.17.648362, which reported ipTM ROC AUC **0.64** and ESM2 log-likelihood correlating
with expression.

### 2.5 Larger but less accessible: the industrial-scale sets

- **[reported]** NVIDIA **Proteina-Complexa**: ~**1 million** binder candidates tested
  against **133 targets** by multiplexed phage screening, with partners including Manifold
  Bio, Novo Nordisk, Viva Biotech, Duke and Cambridge. Described as "one of the largest
  binder design benchmarks to date". *Public availability of the measurements is not
  established — I did not find an open release.*
- **[reported]** Baker lab, Cao et al. *Nature* 2022: experimental data on "nearly half a
  million computational designs and hundreds of thousands of point mutants".
- **[reported]** Chevalier et al. *Nature* 2017: **22,660** designed miniproteins plus
  **6,286** controls against influenza HA and botulinum neurotoxin B; 2,618 high-affinity
  binders.

These are the biggest numbers in the field, and they come from yeast-display or phage
selection, not from per-design kinetic measurement. Different data, different uses: huge
scale, much coarser labels.

### 2.6 BenchBB

**[reported]** Adaptyv's "Bench-tested Binder Benchmark" — a curated set of targets chosen
for ease of E. coli expression and availability of antibody controls (BHRF1, MBP and
others), intended as a standard benchmark for the community.

---

## 3. Part B — how trustworthy is this data?

This is the part that is usually skipped, and it determines what any model built on it can
possibly achieve.

**Assay format changes the answer.** Adaptyv expresses designs cell-free and immobilises
them as a monovalent ligand with the antigen as analyte. Twist expresses them as human
IgG1 Fc fusions in HEK293 cells, captured on an anti-Fc surface — a *bivalent* ligand.
Avidity effects differ. Many reported K<sub>D</sub> values on oligomeric targets are
"apparent", not true 1:1 affinities. **[reported, from the Anthropic release docs]**

**"Binding" is a threshold, not a fact.** The EGFR community paper notes that RFdiffusion
counted any design reaching 50% of the positive control signal as a hit, while BindCraft
used a different micromolar cut-off — so published hit rates are not comparable.
**[reported]**

**Inter-laboratory disagreement is ~11%** on the one dataset where it has been measured
(§2.2). Any classifier scoring above roughly that level of apparent accuracy should be
treated with suspicion.

**Expression is a hidden gate.** **[verified]** In Proteinbase, 3,918 of 4,278 designs
expressed (91.6%). A design that does not express is recorded as a non-binder, so any
metric that predicts manufacturability will masquerade as a binding predictor. I measured
this directly: Boltz-2 complex pLDDT scores **AUC 0.761 for expression** but only
**0.616 for binding among designs that expressed**.

**Self-reported metrics dominate the literature.** Adaptyv's own competition guidance says
plainly that unverifiable self-reported claims weigh little. Most published hit rates are
reported by the authors of the method being evaluated.

---

## 4. Part C — the tools the field actually runs

### 4.1 Generators

| Tool | Approach | Notes |
|---|---|---|
| **RFdiffusion** | backbone diffusion from target structure | the workhorse; Watson et al. 2023 reported ~19% experimental hit rate **[reported]** |
| **BindCraft** | AF2 hallucination, one-shot | reported experimental success 10–100% **[reported]**; 6 of 7 de novo leaderboard spots in EGFR round 2 **[reported]** |
| **AlphaProteo** (DeepMind) | generative, target-conditioned | reported 88% hit rate on BHRF1; **failed entirely on TNF-α** **[reported]** |
| **BoltzGen** | all-atom diffusion, any target type | 519 designs in Proteinbase **[verified]** |
| **Proteina-Complexa** (NVIDIA) | generative, trained on PDB + AFDB + PLINDER + Teddymer | ~1M candidates tested **[reported]** |
| **ProteinMPNN / LigandMPNN / SolubleMPNN** | inverse folding — sequence for a backbone | universal downstream step |

### 4.2 Scorers — what people filter on

| Metric | What it is | Measured performance |
|---|---|---|
| **ipAE / pAE_interaction** | AF2 interface error; the RFdiffusion filter (<10) | baseline; meta-analysis found AF3 ipSAE_min beats it 1.4× in AP **[reported]** |
| **ipTM** | interface predicted TM-score; BindCraft's filter | 0.64 AUC (EGFR paper) **[reported]**; **0.561** in Proteinbase **[verified]** |
| **ipSAE** | interface score from aligned errors | best in the meta-analysis **[reported]**; **0.557** in Proteinbase **[verified]** |
| **pDockQ / pDockQ2** | interface quality from AF2 | **0.454 — points backwards** in Proteinbase **[verified]** |
| **ESM2 pseudo-log-likelihood** | sequence plausibility | correlates with *expression* **[reported]** |
| **pLDDT** | per-residue fold confidence | mostly an expression signal (§3) **[verified]** |

### 4.3 The emerging specificity tools

**[reported]** **Odin-Multi** (bioRxiv 2026.09.08.749745, Jenkins group at DTU with
Johansen and Scapolo — the same people behind the *Science* 2025 pMHC binder work this
project built on in challenge 5) treats cross-reactivity and specificity as explicit
design objectives with counter-selection, rather than leaving them to downstream
screening. On a pMHC benchmark where target and off-target **differ by a single peptide
residue** — our challenge 5 problem exactly — counter-selection raised the fraction of
designs meeting both the target-confidence criterion and a 2.5× target-to-off-target ratio
from **6.0% to 14.2%**.

**It is an in-silico benchmark. No wet-lab validation is reported.**

---

## 5. Part D — what I measured myself on Proteinbase

Three read-only scripts, in [`benchmark/`](../benchmark/). All **[verified]**.

**Binary binding.** 4,271 labelled proteins, 1,082 binders (25.3%). ROC AUC against the
measured outcome:

| metric | AUC |
|---|---|
| ProteinMPNN sequence recovery | 0.727 |
| ESMFold pLDDT | 0.712 |
| Boltz-2 complex pLDDT | 0.639 |
| shape complementarity | 0.634 |
| Boltz-2 ipTM | 0.561 |
| Boltz-2 ipSAE | 0.557 |
| Boltz-2 pDockQ2 | 0.454 (backwards) |

**Conditioned on expression** (3,911 expressed and labelled, 1,082 binders, 27.7%): the
pLDDT family loses most of its advantage; ProteinMPNN sequence recovery holds at 0.713;
the interface metrics stay between 0.45 and 0.58.

**Per target, among expressed designs**, ipSAE ranges from **0.732 (TREM2, n=98)** down to
**0.593 (VEGF110, n=269)**. One global threshold cannot be right for all of them.

**Affinity — the sharpest result.** Among the **1,062** designs with a measured
K<sub>D</sub>, Spearman rank correlation against log K<sub>D</sub>:

| metric | rho |
|---|---|
| TM-score to AFDB50 | −0.350 |
| sequence identity to PDB100 | −0.348 |
| molecular weight | −0.219 |
| novelty | −0.154 |
| Boltz-2 complex pLDDT | +0.101 |
| **ipSAE, ipTM** | **absent from the top 22; ipTM reaches −0.14 on PD-L1 alone** |

**Nothing predicts affinity.** The strongest signals are *similarity to structures already
in the databases* — i.e. "resembles something known" beats every interface score at
predicting how tightly a design will bind. This is consistent with the EGFR community
paper, which found ipAE and ipTM had no significant correlation with K<sub>D</sub> and
"both trended opposite to the assumed direction". **[reported]**

**Caveats on my own numbers:** single pass, no confidence intervals, no multiple-comparison
correction across ~30 metrics, replicates aggregated by majority vote rather than averaged,
and target/design-method/designer-skill are entangled throughout. These are indicative, not
publishable as they stand.

---

## 6. Part E — the gaps, ranked by how open they are

**1. Affinity prediction. Wide open.** No published method and no metric in any dataset I
examined predicts K<sub>D</sub> usefully. The meta-analysis explicitly excluded affinity
for lack of data. Proteinbase has 1,062 affinity measurements under one protocol, and the
Anthropic release adds ~354 more with dual-CRO kinetics and raw curves. This gap now has
data behind it for the first time.

**2. Specificity and counter-selection. Barely started.** One tool (Odin-Multi), in silico
only, September 2026. The field's reviews name negative design and counter-screening as
required and under-served. Almost no public dataset pairs designs with *measured*
off-target behaviour — the Anthropic release's ortholog and CLEC12A counter-screens are a
rare exception.

**3. Target-conditional selection rules. Open.** Everyone applies one global threshold.
Both the meta-analysis (precision 0.1 to 1.0 across targets) and my own measurement
(ipSAE 0.59 to 0.73 across targets) show the right threshold is target-dependent. Nobody
has published a rule for choosing it.

**4. A standard benchmark. Named as missing.** A 2026 review states the field "lacks a
unified, ImageNet-like benchmark". Partially addressed by the DTU meta-analysis and
BenchBB, but neither is standardised across assay protocol.

**5. Expression / developability as a first-class objective. Under-served.** 8.4% of
Proteinbase designs never express, and expression confounds every binding benchmark that
does not control for it. Most design pipelines do not optimise for it at all.

**6. Everything after binding. Wide open and largely untouched by this community.**
Immunogenicity, pharmacokinetics, manufacturability. The Baker review names these as the
remaining obstacles; almost no open data exists.

---

## 7. What this research did NOT cover

Stated so the dossier is not read as more complete than it is.

- I did not download or verify the Zenodo deposit behind the DTU meta-analysis, the EGFR
  competition GitHub repos, or the Anthropic release's actual tables — only its
  documentation. Counts from those are **[reported]**.
- I did not establish whether the NVIDIA Proteina-Complexa measurements are public.
- I did not survey the classical affinity databases (SKEMPI, PDBbind, BindingDB) or the
  antibody-specific ones (SAbDab, FLAb). They cover natural and mutational data rather
  than de novo designs, and may or may not be relevant depending on the direction chosen.
- I did not read the Baker *Nature* 2026 review in full — only the abstract and the
  portions surfaced in search. Its quoted claims should be checked against the PDF before
  being relied on.
- No patent or freedom-to-operate research at all.
- Competitor analysis is shallow: I found three groups working the scoring-benchmark seam
  and did not systematically search for others.

---

## 8. Second research pass — the seven gaps from section 7, closed

Everything in this section is **[verified]**: I downloaded the files and ran the numbers.
Scripts in [`benchmark/`](../benchmark/).

### 8.1 The DTU Zenodo deposit — opened, and its conclusion stress-tested

Zenodo `10.5281/zenodo.15722219`, CC-BY-4.0, published 14 Aug 2025, 7,037 downloads.
Six files; `final_dataset.csv` is 82 MB with **3,676 rows and 312 columns** — every AF2,
ColabFold, Boltz-1 and AF3 confidence score, Rosetta energies, SASA, hydrogen bonds,
secondary structure, SAP scores.

Label column `binder`: **394 True, 3,275 False (10.7%)**.

**The problem the paper does not foreground: one target is 58% of the data.**

| target | rows |
|---|---|
| FGFR2 | **2,123** |
| EGFR | 434 |
| IL7Ra | 171 |
| TrkA | 128 |
| InsulinR | 117 |
| VirB8, SARS-CoV-2 RBD | 99 each |
| pMHC_SILSY1 | 96 |
| Mdm2 | 96 |
| Pdl1 | 95 |
| IL2Ra | 66 |
| sntx | 49 |
| pMHC_NY1 | 41 |
| LTK | 33 |
| IL10Ra | 22 |

So "15 diverse targets" is really one large target plus fourteen small ones, several with
fewer than ten binders.

**I re-ran their comparison three ways.** Their headline conclusion survives — and gets
stronger — when the dominant target is removed:

| metric | all rows | FGFR2 only | everything except FGFR2 |
|---|---|---|---|
| **AF3 ipSAE_min** | **0.787** | 0.753 (12th of 15) | **0.813** |
| Boltz-1 ipSAE_min | 0.768 | **0.781 (1st)** | 0.761 |
| AF3 ipAE | 0.762 | 0.759 | 0.746 |
| AF2 pae_interaction (the RFdiffusion filter) | 0.721 | 0.754 | 0.681 |
| AF2 pLDDT (binder) | 0.673 | 0.723 | 0.596 |

Note that on FGFR2 alone, **AF3 ipSAE is not the best metric** — Boltz-1 is. The published
ranking is a property of the pooled set, not a universal truth.

Per-target AUC for AF3 ipSAE_min ranges **0.573 (Mdm2, n=96) to 0.973 (sntx, n=49)**, with
EGFR 0.800, IL7Ra 0.802, PD-L1 0.876. Several of those rest on fewer than 15 binders.

### 8.2 The Anthropic tables — the single most useful file I found

`data/tables/design_summary.csv`, **1,440 rows and 65 columns**. It carries something no
other public dataset does: **ipSAE_min from ten different structure predictors on
identical designs**, scored against a dual-CRO integrated label.

Generators represented: PXDesign 387, RFdiffusion3 298, Genie3 213, BoltzGen 158,
FreeBindCraft 142, RFdiffusion 120, **Proteina-Complexa 104**, FoldCraft 14.

Labels: 1,320 usable, **354 binders (26.8%)**. Vendor agreement: both bind 253, neither
846, Twist only 67, Adaptyv only 69.

**Ten-predictor head-to-head, same designs, same labels:**

| predictor | ipSAE_min AUC | sc_dockQ AUC |
|---|---|---|
| **Chai-1** | **0.744** | 0.679 |
| **OpenFold3** | **0.741** | 0.679 |
| AlphaFold-Multimer 3 | 0.733 | **0.689** |
| ESMFold2 (full) | 0.721 | 0.632 |
| ESMFold2 (fast) | 0.720 | 0.656 |
| AF3/OF3 ensemble | 0.710 | 0.655 |
| Odde | 0.696 | 0.652 |
| Protenix v2 | 0.692 | 0.627 |
| **Boltz-2** | **0.692** | 0.682 |
| RF3 | 0.688 | 0.657 |

Chai-1 and OpenFold3 beat AlphaFold3 here; Boltz-2 sits near the bottom. **This
contradicts the common assumption that AF3 is the best scorer**, and it is the cleanest
head-to-head available because every predictor saw the identical design set.

**Affinity, again, fails.** Spearman rho of ipSAE_min against log K<sub>D</sub> among the
354 binders: best is AFM3 at **-0.264**, then OF3 -0.243; Chai-1 -0.059 and Boltz-2 -0.054
are effectively zero. **The predictor that best separates binders from non-binders is not
the one that best ranks affinity.**

**Counter-screen data — rarer than I realised.** The table has measured off-target
columns: mouse ortholog (157 binders / 562 non-binders tested), cynomolgus (170 / 429),
**CLEC12A as an unrelated control** (19 / 426), and Cas9-apo (11 / 79).

**Of the 354 confirmed binders, 0 also bound CLEC12A.** The designs were clean against an
unrelated protein. This is one of the only public datasets carrying measured specificity
on designed binders.

### 8.3 The EGFR competition repositories — exist, but note the licence

| repo | stars | licence | last push |
|---|---|---|---|
| `adaptyvbio/egfr_competition_1` | 37 | **NOASSERTION** | Dec 2024 |
| `adaptyvbio/egfr_competition_2` | 22 | **NOASSERTION** | Mar 2025 |
| `adaptyvbio/binder-prescreen` | 7 | NOASSERTION | 9 Nov 2026 |
| `DigBioLab/de_novo_binder_scoring` | 33 | **MIT** | Apr 2026 |

"NOASSERTION" means GitHub could not identify a standard licence. Anything built on the
EGFR competition data needs the licence read properly first — unlike Proteinbase (ODC-BY)
and the Anthropic release (CC-BY-4.0), which are unambiguous.

### 8.4 NVIDIA Proteina-Complexa — models open, measurements not

**The models and weights are public**: `NV-Proteina-Complexa-Protein-Target-160M-v1`
(3,766 downloads), ligand-target and AME variants, plus the Teddymer training dataset and
code on GitHub. Licence is `license:other` — NVIDIA's own terms, not a standard open
licence.

**The million-candidate screening data is not.** I searched NVIDIA's HuggingFace datasets
and found only `RFD4-Proteina`; no binder measurement release. Their reproducibility
statement promises "source code, model weights, and the novel Teddymer dataset" — it does
not promise the wet-lab measurements.

So the largest binder screen ever run (1M+ candidates, 133 targets, up to 63.5% hit rates
**[reported]**) contributes a model to the commons and keeps its data.

### 8.5 The classical affinity databases — smaller than their reputations

| resource | size | what it actually is | access |
|---|---|---|---|
| **SKEMPI 2.0** | 7,085 records, but only **345 crystal structures** | affinity *changes on mutation*, mostly mutants of a few complexes | free |
| **AB-Bind** | subset of SKEMPI | antibody mutations | free |
| **SAbDab** | >7,000 antibody–antigen structures | **14,148 records carry no affinity at all** | free |
| **PDBbind v2020** | 23,496 complexes, of which **2,789 protein–protein** | natural complexes; does not mark which chain is ligand | **v2025 is commercial**; free users get v2020 only |
| **Affinity Benchmark v5.5** | **207** protein–protein samples | curated docking benchmark | free |
| **BindingDB** | 1,342,167 measurements, 3,176 targets | overwhelmingly protein–**small molecule** | free |
| **PPB-Affinity** | integrates the above | largest public PP affinity set | free |

**The implication for any affinity-modelling idea:** the entire curated world of
*structurally resolved protein–protein* affinity is a few thousand complexes, mostly
natural, mostly mutational series around a handful of scaffolds. Proteinbase's **1,062**
de novo designs with K<sub>D</sub> under one protocol, plus the Anthropic release's
**354** with dual-CRO kinetics, are therefore not small by the standards of this field —
they are a meaningful fraction of what exists, and far more homogeneous.

### 8.6 Patents — one broad claim sits across this whole area

**[verified]** `WO2023288191A1`, "Multiple de novo designed protein binding proteins",
inventors Stewart, **Baker**, Cao, Coventry, Goreshnik, Huang, Bennett, Strauch,
Schlichthaerle (University of Washington).

Claim 1 covers a polypeptide **"at least 35% ... identical to the amino acid sequence
selected from the group consisting of SEQ ID NO:1-1559 and 1561-1570"**.

A **35% identity floor across about 1,570 designed miniprotein sequences** is a wide net.
Designed miniproteins are short and share a limited fold vocabulary, so incidental 35%
identity is not hard to hit. Anyone intending to commercialise a de novo miniprotein
should take advice before relying on "we designed it ourselves".

This also explains a detail of our own work: Adaptyv's pre-screen has a **USPTO patent
corpus arm** — the one arm we skipped, on grounds of download size and runtime. A
reasonable call for a competition deadline; not a reasonable call for anything commercial.

I did **not** do a systematic patent search. One landscape analysis and one representative
patent is not freedom-to-operate work, and nothing here is legal advice.

### 8.7 The competitor picture

Groups actively working the scoring and benchmark seam, as far as I found:

1. **DTU DigBioLab / Jenkins group** (Denmark) — the 3,766-binder meta-analysis, the
   `de_novo_binder_scoring` pipeline, **and** Odin-Multi for specificity. They also
   co-authored the *Science* 2025 pMHC work this project built on in challenge 5. The most
   directly overlapping group, and they are moving fast.
2. **Adaptyv Bio** — generating the data, publishing the community paper, running
   Proteinbase and the prescreen.
3. **Anthropic** — the 1,440-design dual-CRO release.
4. One-off efforts: a nanobody folding-model comparison, and "Revisiting Structure-Based
   Filtering for AI-Designed Antibodies".

The seam is occupied but not closed. What nobody has published, as far as I can tell: a
cross-dataset reconciliation of why the same metric scores 0.56, 0.74 and 0.81 on three
different datasets.

### 8.8 The finding that ties this pass together

The same metric, ipSAE, measured against real laboratory outcomes:

| dataset | how it was computed | labels | AUC |
|---|---|---|---|
| Proteinbase | single Boltz-2 pipeline run | one lab, one call | **0.557** |
| Anthropic release | best-of-five seeds, ipSAE_min, ten predictors | **two labs, integrated** | **0.69–0.74** |
| DTU meta-analysis | AF3, re-predicted, ipSAE_min | pooled from papers | **0.787–0.813** |

That is a **0.25 AUC spread on the same quantity**. Published filter thresholds are
therefore not transferable between pipelines, and a paper reporting "ipSAE works" or
"ipSAE does not work" is substantially reporting its own compute protocol and label
quality.

Confounds I have not separated: the three datasets also differ in target mix, binder
rate (10.7% / 26.8% / 25.3%) and design method mix, any of which could move AUC on its
own. The spread is real; the attribution to protocol and labels is a hypothesis, not a
result.

---

## 9. Third pass — feasibility tested, not speculated

Section 6 ranked the gaps. This section asks a harder question about the top one:
**is the affinity gap fillable with the data and features that exist?** Three tests, run
with leave-one-target-out cross-validation so no model can cheat by memorising a target.
All **[verified]**; scripts `benchmark/feasibility.py`, `pooled_vs_within.py`,
`kinetics.py`.

### 9.1 The confirmed blocker

**[verified, from the Boltz-2 documentation]** *"Boltz-2 does not support computation of
protein–protein binding affinity; only small molecules are supported for binding-affinity
predictions."* The team state publicly that they are working on extending it.

So the single most capable affinity model in biology today does not cover this problem,
and the gap is real rather than an artefact of my searching.

### 9.2 Test 1 — can a trained combination beat the best single score? No.

Anthropic release, 1,154 designs with complete features across 15 targets, 27.0% binders,
31 features (ten predictors' ipSAE_min, ten sc_dockQ, length, epitope size, nine sequence
composition features).

| approach | AUC |
|---|---|
| best single score, pooled (optimistic, no held-out penalty) | **0.746** (Chai-1 ipSAE) |
| logistic regression, leave-one-target-out | 0.731 |
| gradient boosting, leave-one-target-out | 0.663 |

**A trained model does not beat the single best number.** Gradient boosting is clearly
overfitting. The ceiling is in the features, not in the modelling.

### 9.3 Test 2 — are published pooled benchmarks inflated? Mostly no. *My hypothesis was wrong.*

I suspected that pooled AUCs overstate real usefulness, because a designer only ever ranks
candidates **within one target**, and pooling can score well purely from between-target
differences in mean score. I tested it on both datasets.

| dataset | mean pooled AUC | mean within-target AUC | difference |
|---|---|---|---|
| Anthropic release (20 metrics) | 0.687 | 0.652 | **0.036** |
| DTU meta-analysis (11 metrics) | 0.741 | 0.738 | **0.003** |

**The inflation is small and on the DTU set essentially zero.** The metrics really do rank
within a target. I record this because I had drafted the opposite claim from a partial
view of the data, and it did not survive being computed properly.

The honest within-target picture: Chai-1 ipSAE averages **0.698** across 12 Anthropic
targets (range 0.51 to 0.85, two targets below 0.60, none backwards); AF3 ipSAE_min
averages **0.771** across 11 DTU targets (range 0.57 to 0.97). Useful. Not decisive.

### 9.4 Test 3 — affinity, including splitting it into its two halves

Among confirmed binders with a measured K<sub>D</sub>:

| approach | Spearman rho vs log K<sub>D</sub> |
|---|---|
| best single feature, pooled (optimistic) | **−0.260** (OpenFold3 ipSAE) |
| ridge regression, leave-one-target-out | −0.121 |
| gradient boosting, leave-one-target-out | −0.176 |

**Again the trained model is worse than the single feature.** There is no combination of
these scores that predicts affinity.

K<sub>D</sub> is a ratio of two separately measurable things — how fast a binder finds its
target (k<sub>on</sub>) and how fast it falls off (k<sub>off</sub>). Almost nobody reports
these separately, so I checked whether one half hides a signal the ratio obscures. 327
binders with full Adaptyv kinetics, 272 with full Twist kinetics:

- **k<sub>off</sub> drives affinity in this data**: rho with K<sub>D</sub> **+0.70**
  (Adaptyv) and **+0.73** (Twist), against −0.62 and −0.61 for k<sub>on</sub>.
- But the predictors track k<sub>off</sub> only weakly: best is **−0.27**.
- On k<sub>on</sub>, **the predictors disagree about the sign**: AFM3 +0.25 and OpenFold3
  +0.19 against Boltz-2 −0.26 and Protenix −0.17, on identical designs. They are not
  measuring the same thing.

**Splitting affinity into its physical halves does not open the door either.**

### 9.5 What these three tests establish

1. **The binary selection problem is at its feature ceiling.** Within-target AUC of about
   0.70 to 0.77 is what current structure-prediction confidence gives you, and combining
   those confidences does not improve it. A better classifier is not the project.
2. **The affinity gap is real, confirmed, and not closable by modelling existing scores.**
   Every honest held-out attempt lands between −0.12 and −0.26. Closing it needs a
   different representation (a learned affinity head, as Boltz-2 did for small molecules)
   or different data — not better regression on confidence numbers.
3. **The people best placed to close it have said they are doing so.** The Boltz team is
   extending Boltz-2 to protein–protein affinity.

### 9.6 What was NOT tested

- No structure-derived features beyond what the two releases ship: no Rosetta energies, no
  buried surface area, no hydrogen-bond counts computed by me. The DTU deposit has 312
  such columns and I used it only for metric comparison, not for modelling.
- No learned representation of any kind — no fine-tuning, no embeddings from a protein
  language model.
- No use of the raw sensorgram curves, which both Adaptyv and Twist ship per replicate.
  Those contain information the fitted K<sub>D</sub> throws away.
- Only one dataset (Anthropic) was used for the modelling tests, because it is the only
  one carrying ten predictors, dual-CRO labels and full kinetics together.

Each of those is a live route rather than a closed door. The tests above rule out the
cheap version of the affinity project, not the expensive one.

---

## 10. Fourth pass — the recent literature, and a fourth dataset that changes the question

Sections 1–9 asked whether the field's scores predict **binding**. This pass started from
the November 2026 literature and ended somewhere else: at a dataset where binding is not
the outcome being measured, and where the scores stop working entirely.

### 10.1 What the recent literature says is blocking the field

Forty-plus papers from 2025–2026 were surveyed (search tools, abstracts read; full texts
read only where cited below). The bottlenecks the field names for itself, in its own
words:

| named bottleneck | where it is stated | status |
|---|---|---|
| Affinity is not predictable | our own tests, §9; no method found | **open** |
| Scores disagree depending on which predictor verifies them | ProtDBench, arXiv 2605.04118 **[verified, abstract]**: "substantial verifier-dependent bias and limited agreement under identical filtering protocols" | **named, now being worked** |
| The right filtering threshold is target-dependent and nobody knows how to pick it | DTU meta-analysis full text **[verified]**: precision "ranged from 0.1 to 1.0 across targets at the maximum F1 threshold"; "choosing the right threshold ... remains challenging" | **open, explicitly unsolved** |
| Specificity and counter-selection are left to downstream screening | Odin-Multi (Sep 2026) **[reported]**: deep-learning design methods "typically optimise one interaction at a time, leaving cross-reactivity and specificity to emerge during downstream screening" | **barely started** |
| Some targets yield zero binders from unconstrained methods | MIMOSA (Sep 2026) **[reported]**: 20–26% hit rates on KEAP1, uPA and IL-8, "on which unconstrained methods produce zero or near-zero binders" | **being worked** |
| Designs are made far faster than they can be tested | cell-free FIDA pipeline (Oct 2026) **[reported]**: existing validation "take 2 to 8 weeks" | wet-lab, not ours |
| Shortlisting from large pools is the practical bottleneck | LLM-shortlisting paper (Aug 2026) **[reported]**: Recall@10 of 0.589 against 0.571 for the best single score | **attempted, barely moved** |
| **Evaluation stops at binding and ignores higher-order function** | Bits to Binders (Mol Syst Biol, Mar 2026) **[verified, full data]**: "most evaluations emphasize binding affinity rather than higher-order biological function" | **open — and now has data** |

Two of these deserve flagging because they overlap our own work. ProtDBench's
verifier-dependent bias **is the same phenomenon as our 0.25 AUC spread** in §8.8, found
independently and published first. And the LLM-shortlisting result — a trained ranking
policy beating the best single score by 0.018 — **replicates our §9.2 finding** that
combining confidence scores does not help.

### 10.2 The fourth dataset: Bits to Binders

**[verified]** Downloaded and analysed here. `github.com/kosonocky/bits-to-binders`, MIT
licence; also Zenodo 20451404. Scripts: `benchmark/b2b_*.py`.

| | |
|---|---|
| designs | 12,000 CD20-targeting 80-mers |
| teams | 29, from 42 countries |
| features | **433 computed columns** — Boltz-2, Chai-1, ESMFold at four linker positions, Rosetta interface energies, ESM2 likelihoods, MPNN scores, solubility, composition |
| outcome | not binding — **CAR-T proliferation in a pooled screen** |
| labels | 707 Up, 5,941 Not Sig, 166 Down, 5,186 never recovered |
| plus | 101 designs with Adaptyv SPR (K<sub>D</sub>, k<sub>on</sub>, k<sub>off</sub>), and full cytotoxicity / cytokine / expansion assays on the top ten |

This is the **mirror image** of the three datasets in §2: those are many targets measured
by few methods; this is one target measured by 29 methods, scored against function rather
than affinity.

### 10.3 Test 4 — does measured binding predict measured function? No detectable relation.

> **Corrected in §10.8.** The 100 designs were sampled by outcome, so the percentages below
> are not population rates. The association test stands; the framing does not.

100 designs carry both an SPR call and a CAR-T outcome.

| SPR says | drove proliferation (Up) | did not | rate |
|---|---|---|---|
| binds | 5 | 9 | **36%** |
| does not bind | 19 | 50 | **28%** |
| did not express for SPR | 6 | 11 | 35% |

Fisher exact p = **0.53**. The designs that bind and the designs that do not bind succeed
as CARs at indistinguishable rates.

**Read this carefully, in both directions.** The sample is tiny — 14 binders — so this
does not *prove* the two are unrelated; it shows no evidence that they are related, with
nothing like the power to establish equivalence. The 101 were also selected rather than
random, which can flatten any real difference. And there is an innocent explanation that
is probably part of the truth: SPR measures a soluble protein in a cuvette, while the CAR
assay measures the same sequence anchored in a membrane at high density, where avidity and
geometry do work that affinity does not.

But that explanation is itself the point. The number the whole field optimises is measured
in the cuvette, and the thing anybody actually wants happens in the membrane.

### 10.4 Test 5 — do the 433 scores predict function? The structural ones do not.

> **This is a replication, not a finding — see §10.8.** The Bits to Binders paper reports
> the same result for the same metric list.

Among the 6,814 recovered designs, 707 Up:

| metric | AUC for function | the same family's AUC for *binding* elsewhere in this dossier |
|---|---|---|
| Boltz-2 ipSAE (binder vs whole target) | **0.477** | 0.56 / 0.69 / 0.79 |
| Chai-1 ipSAE | **0.507** | 0.744 (best in the Anthropic set) |
| Rosetta ΔG separated | **0.441** | — |
| Boltz SAP (aggregation) | 0.447 | — |
| binder pLDDT | 0.584 | 0.712 |
| ESM2 binder log-likelihood in CAR context | 0.585 | 0.72 reported by Adaptyv |

**Every interface-confidence score is at or below chance.** Not weak — *absent*. The same
metrics that reach 0.69–0.81 on binding sit at 0.44–0.51 here.

### 10.5 Test 6 — and the features that *do* score are mostly a team fingerprint

> **Largely superseded — see §10.8.** The source paper found this signal first, at 0.91 for
> recovery, and gave it a molecular explanation this section does not have.

The top of the function table is not structural at all. It is composition:
`ratio_KE` 0.726, `ratio_E` 0.708, poly-E stretches 0.673, poly-A DNA repeats 0.682.

That looked like a finding until it was tested against the obvious confound. Team hit
rates run from under 1% to **48.2%** (Nucleate UK London, 159 of 330). Any feature
describing a strong team's house style will score well when teams are pooled.

| feature | pooled AUC | mean **within-team** AUC (25 teams) |
|---|---|---|
| `ratio_KE` | 0.726 | **0.576** (range 0.43–0.75) |
| Boltz ipSAE | 0.477 | 0.498 (range 0.37–0.78) |

**Three quarters of the composition signal is team identity, not biology.** What survives
— 0.576 — is weak but not nothing, and it is the strongest honest signal in the table.

Note this cuts against §9.3, where pooling across *targets* inflated AUC by only 0.036.
Pooling across *designers* inflates it by 0.15. Designer skill is a bigger confound than
target identity, and no published benchmark I have seen controls for it.

### 10.6 What this pass establishes

1. **The field's scoring stack predicts binding and does not predict function.** That is
   not a claim about one bad metric; it is twelve metrics from four predictors, all at
   chance, on 6,814 designs with 707 positives.
2. **Measured binding itself did not predict function** in the only sample where both were
   measured, though that sample is too small to settle it.
3. **The confound that matters most in multi-designer benchmarks is who designed it**, and
   it is larger than the target confound everybody does control for.
4. The gaps the literature names as open and unclaimed are: **target-conditional
   thresholds** (the DTU authors say plainly they could not solve it), **specificity with
   measured off-target data**, and **anything downstream of binding**.

### 10.7 What was NOT done in this pass

- The Bits to Binders paper's own analysis has not been read in full — only its abstract
  and the data. Some of §10.4–10.6 may duplicate analysis already in their figures, and
  that needs checking before any of it is called new.
- "Beyond Structure and Affinity" (Apr 2026) already mined this same dataset with
  biology-informed sequence features and reports lifting hit rate 13.8% → 38.6%. Its
  relationship to §10.5 is unexamined; it may already cover the composition finding, and
  it is not clear whether it controlled for team.
- No statistical correction across the 433 features tested. The top AUCs are not corrected
  for multiple comparisons.
- The 5,186 never-recovered designs are excluded from every function number above. That is
  43% of the dataset, dropped on an outcome that is itself predictable.
- No reading of the `misc/` team-methods files, which would say what each team actually ran
  and would let the team confound be replaced with a *method* variable.

### 10.8 Correction — the two papers were read, and most of §10.3–10.6 is not new

Section 10.7 said the Bits to Binders figures were unread and that §10.4 might be a
replication. Both papers have now been read in full. **It is a replication, and so is most
of §10.5.** Recorded here rather than quietly edited, because the previous sections were
written as findings and were circulated as findings.

**On §10.4 (structure metrics at chance).** The Bits to Binders paper states it directly,
with a supplementary table behind it **[verified, full text]**:

> "global confidence statistics from structure prediction models, such as ipTM and pLDDT,
> did not predict experimental outcomes (Table SA1). Likewise, ipSAE, LIS score, SAP
> score, Rosetta InterfaceAnalyzer metrics, PDockQ..."

That is the same metric list and the same conclusion. §10.4 is an **independent
replication** with per-metric AUCs, not a finding.

**On §10.5 (the composition signal and its confound).** Worse: they found the signal,
quantified it better, and explained it.

| what | theirs | mine (§10.4–10.5) |
|---|---|---|
| strongest predictor of recovery | K+E α-helix ratio, **0.91** | DNA sequence entropy, 0.897 |
| DNA sequence entropy for recovery | **0.90** | 0.897 |
| K+E also enriched in functional hits | yes, Fig 5G, "to a lesser extent" | ratio_KE 0.726 pooled |
| why | **98.9%** of non-recovered designs came from ProteinMPNN/SolubleMPNN, which preferentially emit K and E; EE repeats and poly-A codons stall the ribosome | "a team fingerprint" |
| confound control | restricted analyses to MPNN-only subsets | within-team AUC |

Their explanation is better than mine and is mechanistic. My "team fingerprint" and their
"MPNN artefact" are the same confound seen from different sides — teams that used MPNN
produced K/E-rich helices — except they named the molecular cause and I named a proxy for
it. The within-team drop 0.726 → 0.576 is consistent with their account and adds a number
they did not publish, which is the most that can be claimed for it.

**On §10.3 (binding vs function), two corrections.**

First, the 100 SPR designs were **not a random sample**. The paper states they were
stratified by outcome in advance: the top ten, 20 more that enriched ≥2-fold, 10 that
depleted ≥2-fold, and 60 with no significant change. So the "36% of binders worked versus
28% of non-binders" framing in §10.3 reports rates that the sampling design fixed. **That
framing was wrong and is withdrawn.** The association test survives — this is a
case-control layout and Fisher's p = 0.53 on it is legitimate — but the percentages are
not population rates and must not be read as any.

Turned the right way round, the same 100 designs say:

| CAR-T outcome | fraction with detectable SPR binding |
|---|---|
| enriched (Up) | 5 of 30 — **17%** |
| no change | 7 of 60 — 12% |
| depleted (Down) | 2 of 10 — 20% |

Second, the paper already reports this qualitatively for its top ten: only three showed
appreciable binding, and the single quantifiable one came in at **K<sub>D</sub> = 643 nM**
against **1.75 nM** for the known scFv positive control. Extending it to all 100 is a
modest addition, not a discovery.

### 10.9 What survives, and what this costs us

**Survives as ours:**

- The cross-dataset confound comparison. Pooling across targets inflates AUC by **0.036**
  (§9.3); pooling across designers inflates it by **0.15** (§10.5). Neither paper makes
  that comparison, because neither has the other three datasets open beside this one.
- A harness that reproduces published numbers from raw data — DNA entropy 0.897 against
  their 0.90 — which is the reason to trust the rest of `benchmark/`.

**Does not survive:** both headline claims from the previous section. §10.4 is a
replication; §10.5 is a weaker version of a published result.

**What this costs strategically.** The "nobody can predict function" gap is **not
unclaimed**. At least two groups are working this dataset: its own authors, and the April
2026 biology-informed-features group. That is three times in two weeks that a direction
this project identified was already in print. The pattern is consistent and it is worth
stating as a conclusion rather than an accident: **the obvious analyses on the public
datasets are being done, quickly, by people closer to the data.** Any project that is
"re-analyse a public benchmark" is a race we start behind.

Of the gaps in §6 and §10.1, the one that has now been looked at hardest and still has
nobody on it is **target-conditional thresholds** — the DTU authors say in print that they
could not solve it, and no paper found in this pass claims to have. That is not a
recommendation; it is the one door that was still shut after this pass.

---

## 11. Fifth pass — the decade question, and where it is blocked

Sections 1–10 worked on binders, because that is what this project built. This pass asked a
larger question: what is the field's next decade actually about, and what is stopping it?
Roughly 120 papers and two authoritative roadmaps were surveyed; everything below is marked
**[verified]** where the source was read, **[reported]** where only an abstract was.

### 11.1 Two roadmaps, written independently, say the same thing

**[verified, abstract]** *The past, present and future of de novo protein design*, Yang,
Wang, Lee, Zhang … Baker, **Nature**, April 2026:

> "the long-standing challenges of designing new protein structures, assemblies and protein
> binders are **close to being solved**. The key current questions in these areas are not how
> to design, but **what to design** … Current challenges for methods development include
> design of catalysts for reactions with high energy barriers and, more generally, **design
> of switches and nanomachines that integrate binding, conformational change and catalysis**."

**[verified, full text]** *Programming biological behavior with AI-guided de novo protein
design*, Huettemann & Kortemme, **Molecular Cell**, October 2026. Its thesis is that the
field must move from bimolecular recognition to **multi-component reaction networks** —
signal amplification, switch-like responses, kinetic proofreading, error correction, and
directionality driven by energy expenditure. It names the blocker twice:

> "A key bottleneck is **datasets that inform on protein conformational landscapes at a scale
> that is sufficient to develop models that generalize** across proteins and mechanisms."

> "Critical to improved models are large-scale datasets of both informative parameters that
> quantify protein function (enzymatic activities, binding rates, affinities, **protein
> dynamics**) and their dependency on changes in atom-level interactions and cellular context."

Two groups that do not agree about much agree on this: **binders are finished; the decade is
about proteins that move; and the thing standing in the way is data about how proteins move.**

This project spent a month on binders — the problem both roadmaps call close to solved.

### 11.2 The blocker is real, and three independent sources measure it

**CASP has a category for this and it is starving. [verified]** The Critical Assessment of
Structure Prediction added *Macromolecular Conformational Ensembles* in 2022. The CASP16
assessment (Dube, Ramelot, … Moult, Kryshtafovych, Montelione, 2025; PubMed, [DOI](https://doi.org/10.1101/2025.09.02.673835)):

- **Ten** ensembles were released as community challenges. Single-structure CASP16 had 300 targets.
- "predictors **generally failed to capture key structural details distinguishing the states**"
- "accuracy was **significantly lower than for single-state targets**"
- CASP's own call for targets: "**Success of the macromolecular ensembles category depends on
  the availability of targets with multiple experimentally observed conformations**", and asks
  the community for cryo-EM, NMR, X-ray, SAXS or cross-link data to assess ensembles with.

**The models disagree with each other and nobody can score them. [reported]** An independent
comparison of BioEmu, AFSample2 and ESMFlow (Otten, Leung, Chong, Zuckerman, May 2026):

> "each tool produces clearly distinct conformational distributions. Thus, **it is impossible
> that all the tools generate ensembles of high biophysical quality** … our current analysis
> **cannot evaluate the absolute quality of an ensemble, and hence cannot identify a
> best-performing AI tool**."

A separate assessment of BioEmu (March 2026) finds it "does not reproduce a right
Boltzmann-weighted ensemble" and "fails to predict a mutation-induced shift in conformational
distribution" — which is precisely the operation a designer needs.

**The constraint is data, not architecture. [reported]** An ablation study published
**3 November 2026** — ten days before this was written — strips a flow-matching ensemble
generator down and varies architecture, conditioning, alignment and training source against
experimental SAXS and NMR:

> "This points to **the training data, rather than the architecture, as the main constraint**."

### 11.3 How thin the experimental ground truth is

| resource | what it holds | **[verified]** |
|---|---|---|
| SASBDB (small-angle scattering) | 5,716 experimental datasets, 4,508 protein — but only **298 ensemble models** | read from the site, 13 Nov 2026 |
| CASP16 ensemble category | **10** targets | assessment paper |
| ProteinConformers (Zhang group, Feb 2026) | 2.7M conformations — **generated by MD**, not measured | abstract |
| ATLAS / mdCATH | MD trajectories — simulation, not observation | abstract |
| HDX-MS | ~20 analysis packages found, **no standardised public repository** | absence of evidence, see caveat |

The pattern: nearly every ensemble model in this field is trained on molecular dynamics and
evaluated against molecular dynamics. MD is a model with a force field, not a measurement.
The field is, for the second time in this dossier, **grading its own homework** — and this
time at the exact point where its next decade depends on the answer.

### 11.4 The opening — and the group that got there nine months ago

Everything above concerns **natural** proteins. CASP targets are natural. ATLAS is natural.
Fold-switch benchmarks are natural. But the application both roadmaps name is **designed**
proteins — switches, receptors, machines — and designed proteins are the one class these
models have the least right to work on, because they have no evolutionary history and most of
these methods lean on multiple sequence alignments.

Three papers establish that this matters, and that the work has started:

- **[reported]** *Fold-switching proteins push the boundaries of conformational ensemble
  prediction* (Jan 2026): "DL models often predict conformational ensembles **by association
  with training-set structures**, limiting generalizability."
- **[reported]** *Mutational Robustness Predicts Protein Dynamics* (Mar 2026) covers ~2,000
  natural proteins and **~400 de novo designs**, and reports the "largest gains on designed
  proteins", with the signal "equally strong on de novo designs that lack any evolutionary
  history, pointing to a biophysical effect rather than a proxy for sequence conservation."
- **[reported], and the most important paper in this pass:** *Large-scale exploration of
  protein space by automated NMR* (Feb 2026). A pipeline producing **hundreds of isotopically
  labelled proteins per week** with one operator. They ran **384 de novo designs**, got usable
  2D NMR for **239 (62%)**, and report:

  > "NMR characterization confirmed that the designed proteins adopt their intended folds, and
  > **revealed unexpected local dynamics that are not captured by current computational
  > models**."

So the thesis — *the field cannot predict the dynamics of the designed proteins it now wants
to build* — is **not unclaimed**. One group has measured it on 384 designs and said so.

But note the scale. In the binder-benchmark space there were 4,280-row answer keys and a
dozen groups mining them, and this project was scooped three times in two weeks. Here there
are **three papers and one dataset of 239 spectra**. This is a field at its beginning, not a
race already run.

### 11.5 What this pass establishes

1. **The decade's question is designed dynamics** — proteins that switch, transduce, compute
   and expend energy, not proteins that stick. Two independent authoritative roadmaps.
2. **The blocker is measured ground truth about motion**, named by Kortemme, by CASP, and by a
   November 2026 ablation that puts it on data rather than architecture.
3. **The field currently validates motion models against simulations of motion.** That is the
   same category of error this project documented for binder scoring in §8 and §10 — one level
   up, and far more consequential.
4. **Designed proteins are the hardest and most relevant case**, and the one with almost no
   data. 239 NMR spectra exist. That is the entire public corpus I could find.

### 11.6 What was NOT done, and what would change the conclusion

- The Nature roadmap is **paywalled**; only its abstract was read. Its own list of challenges
  may be more specific than quoted here.
- The automated-NMR paper was read as an abstract only. **Whether its 384 designs and 239
  spectra are publicly downloadable has not been checked**, and that single fact changes what
  is buildable.
- "No standardised public HDX-MS repository" is an inference from not finding one across two
  searches. It is the weakest claim in this section and should be checked directly before
  anyone relies on it.
- No survey of cryo-EM conformational heterogeneity resources (cryoDRGN outputs, EMDB
  heterogeneity depositions), which is a plausible large and growing source of measured
  ensembles and was not examined at all.
- No assessment of who is funding this area, which would say how fast it is about to get
  crowded.
- Smaller frontier topics surfaced and were not pursued: de novo enzymes for high-barrier
  reactions (Nature, Oct 2026: de novo cysteine proteases still short of native activity),
  nanopores, protein materials, and the "de novo binderome" idea — a designed binder for every
  protein in a proteome, with BindCraft already reporting high-confidence design models against
  over 40% of the human proteome.
