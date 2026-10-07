# Challenge 1 — Conditional EGFR Binder: Research & Execution Plan

**Participant:** Talal Ahmed · Track 3 (open track)
**Window:** Sep 28 – Oct 4, 2026 (submissions close Oct 4, end of day AoE)
**Plan written:** Sep 29, 2026 — 5 working days remain

---

## 1. What the challenge actually asks for

Design a **de novo** binder to human EGFR extracellular region (UniProt P00533, residues 25–645; structure reference PDB 6ARU chain A) that satisfies three objectives:

| # | Objective | Measurement |
|---|---|---|
| 1 | Binds human EGFR | Affinity vs full human ECD |
| 2 | Cross-reacts with mouse EGFR | Binding vs mouse ortholog |
| 3 | pH-selective | Binds at pH 6.5, **no detectable binding** at pH 7.4 |

**Official ranking order: (1) pH-selectivity → (2) mouse cross-reactivity → (3) affinity.**

Competition FAQ 3, verbatim in substance: *pH-sensitive binder design is generally much harder than designing high-affinity binders, so a weak, but clearly pH-sensitive binder may be considered more impactful than a high-affinity binder that is not pH-sensitive.*

### Submission constraints (hard)
- **≤20 designs** (Track 3), submitted as CSV ranked best-first
- Columns: `name`, `sequence`, `molecule_class` ∈ {protein, nanobody, scfv, fab_kappa, fab_lambda}
- Length 10–250 aa; single chain, nanobody, scFv, or Fab
- **De novo and zero-shot only** — may not start from an existing binder; must pass novelty screen
- Evaluation categories: microbinder (<40 aa), minibinder (40–100), large (>100), nanobody, antibody

### How Track 3 designs get picked
~1500 designs screened per challenge; **25% (~375) from Track 3**. All Track 3 submissions are pooled and given to Claude with a selection prompt (not published in advance) that weighs **predicted design quality, design novelty, and method novelty**, plus whatever supporting information you submit. Organizers explicitly invite a methods write-up and metadata package, and state that richer information lets Claude make a better-informed decision.

**Implication:** this is not a pure "highest ipSAE wins" contest. Documentation quality and method distinctiveness are part of the selection function. That is the opening for an entrant without a GPU cluster.

---

## 2. Three facts that determine the design strategy

### Fact 1 — The recommended epitope conflicts with Objective 2
The challenge page recommends domain III, the cetuximab/panitumumab epitope. But **cetuximab does not appreciably bind murine EGFR despite ~88% sequence identity across domain III**. Human EGFR ECD is 88–90% identical to mouse, yet a handful of substitutions inside the cetuximab footprint abolish cross-species binding.

→ Targeting the canonical cetuximab hotspot maximizes therapeutic credibility but is likely to **fail Objective 2**. We must select a sub-surface of domain III (or an alternative functional site) that is conserved human↔mouse. This is the single highest-leverage decision in the whole challenge, and it costs no GPU time.

### Fact 2 — Cross-reactivity is real but only partly predictable
In Anthropic's own 1,320-design campaign: 130/233 binders tested against mouse orthologs cross-reacted, but the outcome was strongly target-dependent — and **on EGFR specifically only 4 of 8 binders bound the mouse ortholog**. They also found that counting substituted footprint residues **did not predict** cross-reactivity (no significant Cochran–Armitage trend; Spearman ρ = 0.07 for affinity loss).

→ Epitope conservation is necessary but not sufficient. Treat it as a filter, not a guarantee, and hedge across more than one conserved patch rather than betting all 20 designs on one.

### Fact 3 — "Acid-ON" switching has no published experimental precedent
The current state of the art (*Computational design of pH-sensitive binders*, bioRxiv 2025.09.29.678932) is **acid-OFF**: histidines placed next to *cationic* target residues so that protonation at low pH creates repulsion and **weakens** binding (up to 1000-fold weaker at pH 5.4 on EphA2; 122-fold on TNFR2; 79-fold on TNFα). Reported success rates are low — e.g. 6 of 43 designs >2-fold on TNFR2, 3 of 72 on TNFα — and the paper notes that random histidine placement rarely works.

This challenge needs the **inverse**: binding present at pH 6.5, absent at 7.4. The same paper flags the mechanism only as theory — histidines as hydrogen-bond **donors** rather than acceptors — with no experimental examples.

→ Our mechanism: **protonation-dependent salt bridges.** Interface histidines positioned against target **carboxylates (Asp/Glu)**. At pH 6.5 the imidazolium is protonated and forms a charged H-bond/salt bridge; at pH 7.4 the neutral imidazole loses it, leaving a desolvated unpaired carboxylate — an interaction that is not merely weaker but specifically absent. This is novel-method territory, which the selection criteria reward, and it must be stated honestly as unproven.

---

## 3. Strategy

**Optimize in the order the organizers rank, not in the order that is easy.** Most entrants will run a standard binder pipeline, filter on interface confidence, and submit high-affinity designs with no pH mechanism. Deliberately inverting that priority is the edge.

1. **pH switch is the primary design objective**, engineered explicitly, with a stated physical mechanism and a paired pH-insensitive control set from the same epitope.
2. **Cross-reactivity is designed in at the epitope-selection stage**, before any backbone is generated — the only stage where it is cheap.
3. **Affinity is a constraint to satisfy, not the objective to maximize** — enough to be detectable at pH 6.5.
4. **Diversify method and molecule class** to earn novelty credit and to avoid a single point of failure.

---

## 4. Phased plan

### Phase 0 — Setup (Sep 29, ~1 hour)
- [ ] Create Proteinbase account; locate the submission form and CSV template
- [ ] Join the competition Slack (only support channel; email support is closed)
- [ ] Re-read Adaptyv's novelty post and the challenge FAQ; keep both as local reference
- [ ] Ask Dr. Farakh Javed (FYP supervisor) about the **Claude Team plan for scientists** — the rejection email states standard seats are currently free for verified research labs, and an academic PI can register the lab. Legitimate route to more capacity after the Track 2 rejection.
- [ ] Compute decision: free Colab GPU (T4) as baseline; Colab Pro (~$10) if a longer runtime is needed; Modal free credits for heavier jobs. Verify current free-tier terms before relying on them.

### Phase 1 — Epitope selection (Sep 29–30, CPU only, no GPU)
Highest-leverage phase. Sequences are already staged in `data/`.

- [ ] Pairwise align human P00533 vs mouse Q01279 ECD (both 1210 aa full-length; expect near-ungapped)
- [ ] Map per-residue identity onto 6ARU chain A; compute solvent accessibility
- [ ] Enumerate candidate surface patches, scored on:
  - **Conservation** — contiguous, solvent-exposed, 100% identical human↔mouse (≥15–20 residues)
  - **Carboxylate density** — Asp/Glu available for His⁺ pairing (required for the pH mechanism)
  - **Functional relevance** — ligand-binding surface of domain III, or the domain II dimerization arm, so a bound miniprotein could plausibly block function
  - **Designability** — geometry that a small binder can pack against
- [ ] Output: 2–3 ranked sub-epitopes with explicit hotspot residue lists, plus a short justification document
- [ ] Prepare a **truncated target construct** (domain III ± flanks, ~170 aa) for design on limited GPU memory; retain full ECD for final scoring

### Phase 2 — Backbone and sequence generation (Sep 30 – Oct 1, GPU)
- [ ] Primary: **BindCraft / BindCraft2** (Colab) against each hotspot set
- [ ] Second: **RFdiffusion + ProteinMPNN** via ColabDesign
- [ ] Third (time permitting): **BoltzDesign1**
- [ ] Mirror Anthropic's diversity rule: draw final designs from **≥3 generators, none contributing more than half**
- [ ] Sequence design with **SolubleMPNN**, applying a **histidine bias at interface positions** (the pH paper used ProteinMPNN His bias 0.5–1.5)
- [ ] Target a few hundred backbones, not thousands — quality over volume given the compute budget

### Phase 3 — pH-switch engineering (Oct 1–2, mostly CPU)
- [ ] For each surviving design, find interface positions within H-bond distance (~2.7–3.5 Å) of target Asp/Glu carboxylates
- [ ] Introduce **1–4 interfacial histidines**; the pH paper's working switches needed **≥2 histidine–charge contacts**, so single-His designs are deprioritized
- [ ] Compute switch evidence:
  - ΔΔG_electrostatic between protonated-His and neutral-His states
  - Predicted His pKa (PROPKA or equivalent) — want **≈6.5–7.0** so the residue is protonated at 6.5 and neutral at 7.4
  - Burial/desolvation check — buried unpaired carboxylate at pH 7.4 is the source of the "off" state
- [ ] **Negative design check:** confirm the neutral state loses a *specific* contact, not just overall affinity
- [ ] Reserve 2–4 slots for **pH-insensitive controls at the same epitope** — a paired design/control set is legible, defensible science and strengthens the submission

### Phase 4 — Scoring, novelty, ranking (Oct 2–3, GPU)
- [ ] Co-fold and score using the organizers' own prescribed metric where runnable: **ESMFold2 / ESMFold2-Fast / Protenix v2 ensemble → ipSAE_min**, max over seeds, z-scored per target, with **self-consistency DockQ at ¼ weight**. Fall back to Boltz-2 / Chai-1 / AF2-multimer (ColabFold) and **state the substitution honestly**
- [ ] **Co-fold against mouse EGFR as well** — require acceptable scores on both species
- [ ] **Novelty screen (disqualification risk, do not skip):**
  - MMseqs2 vs UniRef90 / SwissProt / PDB — target **≤30% sequence identity**
  - Foldseek → TM-align — avoid high structural similarity (TM ≥0.8 over ≥70% coverage); stay below moderate similarity where possible
  - If any antibody-format design is used: CDRH3 edit-distance identity **<70%** vs PLAbDab and global identity <70%
- [ ] Final ranking: **pH-switch evidence → cross-species conservation + mouse co-fold → affinity proxy**
- [ ] Decide molecule-class split. Minibinders (40–100 aa) are the crowded default with the best tool support and highest expression rates; consider placing a minority of designs in a thinner category since winners are announced per category. Do not spread so thin that quality drops.

### Phase 5 — Submission package (Oct 3–4)
- [ ] `submission.csv` — top 20, ranked, correct columns and molecule classes
- [ ] Per-design metrics table (all in silico values, with the method that produced each)
- [ ] Design models / predicted structures
- [ ] **Methods write-up** covering: epitope-selection rationale and the cetuximab cross-reactivity trap, the acid-ON mechanism and why it is inverted relative to published work, generators and settings, scoring protocol and any substitutions, novelty screen results, and explicit limitations
- [ ] Public GitHub repo with code and data; link it in the submission
- [ ] Review every design personally before submitting — required by the rules
- [ ] **Submit by Oct 4 end of day AoE.** The FAQ separately mentions 5:00 AM PDT the following Tuesday; do not rely on the later reading
- [ ] No embedded instructions or prompt-injection-style text anywhere in the submission — explicit grounds for disqualification

---

## 5. Risk register

| Risk | Mitigation |
|---|---|
| Free GPU limits block generation | Truncate target to domain III; small binders; Modal credits; Colab Pro as fallback |
| Acid-ON mechanism is unproven | Hedge with paired pH-insensitive controls; rank switch designs first but keep affinity designs in the 20 |
| Mouse cross-reactivity unpredictable (4/8 on EGFR in Anthropic's run) | Spread designs across ≥2 conserved patches; use conservation as filter, not guarantee |
| Histidine engineering kills affinity outright | Limit to 1–4 substitutions; verify interface confidence after mutation, not before |
| Novelty filter disqualification | Run the screen before submission, not after; keep ≤30% identity |
| Time (5 days) | Phase 1 is CPU-only and starts immediately; generation is the only hard GPU dependency |

---

## 6. Honest positioning

No prior de novo binder design experience, no GPU cluster, no Claude credits. The compensating assets are: a real molecular-biology background in construct and purification design, a documented habit of diagnosing *why* a construct fails rather than only whether it works, willingness to attempt the objective that is ranked highest and attempted least, and a submission package built to be read and trusted. Track 3 selection rewards exactly that combination.

---

## 7. Immediate next step

Phase 1 — the human/mouse conservation map and ranked epitope shortlist. CPU-only, executable now, and it determines everything downstream.
