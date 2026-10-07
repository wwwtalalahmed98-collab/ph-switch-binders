# Challenge 1 — conditional EGFR binder: methods

**Track 3** · Talal Ahmed, PAF-IAST, Pakistan
**Status: DRAFT — design-specific numbers marked `[TBD]` pending generation**
**Code and data:** all scripts, target constructs and analysis outputs in the linked repository

---

## Summary

We designed de novo minibinders against a surface of human EGFR chosen so that the two
hardest objectives — mouse cross-reactivity and acid-ON pH selectivity — are satisfied by
the *target surface* before any backbone is generated, rather than being optimised for
afterwards.

The central decision is that we **did not target the cetuximab epitope**, despite it being
the epitope the challenge page recommends. Three independent lines of structural evidence
say it cannot satisfy objectives 2 and 3 simultaneously.

---

## 1. Epitope selection

### 1.1 Why not the recommended epitope

We computed the cetuximab footprint directly from PDB 6ARU as EGFR residues with any atom
within 5 Å of the Fab: 26 residues (P373, V374, R377, L406, Q408, Q432, H433, Q435, F436,
A439, V441, S442, I462, S464, G465, N466, K467, K489, I490, I491, S492, N493, R494, G495,
E496, N497).

1. **Cross-species.** 7 of those 26 are substituted in mouse, and they are the closest
   contacts: S492N (2.82 Å from the Fab), I491M (3.05 Å), N497K (3.08 Å), S442G (3.26 Å),
   K467R (3.39 Å), G495A (3.44 Å), R377K (3.61 Å). This is a structural account of the
   known failure of cetuximab to bind murine EGFR despite ~88% domain III identity, and a
   direct warning for objective 2.
2. **No chemistry for the switch.** The footprint contains exactly **one carboxylate
   (E496) in 26 residues**. An acid-ON histidine switch requires target Asp/Glu to pair
   with.
3. **Geometrically incapable.** Running our switch-placement engine on the cetuximab–EGFR
   interface returns no viable histidine position: the closest approach from any Fab Cβ to
   an EGFR carboxylate is 8.21 Å, outside the salt-bridge window entirely.

*Method check:* the same footprint calculation overlaps the EGF ligand site at 17/26
residues, consistent with cetuximab's known mechanism as a ligand blocker.

### 1.2 Conservation mapping

Human (P00533) and mouse (Q01279) EGFR were aligned and per-residue identity mapped onto
6ARU chain A, with solvent accessibility computed on EGFR alone (computing it with the Fab
bound artificially buries the cetuximab face and would bias patch selection).

| Region | Human/mouse identity |
|---|---|
| Extracellular region (25–645) | 88.7% |
| Domain I (L1) | 85.4% |
| Domain II (CR1) | 93.2% |
| Domain III (L2) | 87.7% |
| Domain IV (CR2) | 89.3% |

Of 284 candidate surface patches (13 Å radius, ≥12 exposed residues), **48 are 100%
identical** between species.

### 1.3 Functional sites

Computed empirically from PDB 1IVO (EGF-bound, dimerised): the EGF footprint (43 residues,
spanning domains I and III) and the receptor–receptor dimer interface (33 residues,
including the domain II dimerisation arm P266–M277).

### 1.4 Selected epitope

**Epitope A — domain I ligand-binding face:** F44, E45, D46, L49, Q52, Y74, D75.

- 100% identical human↔mouse across the patch
- 6–7 residues inside the EGF footprint, so a bound minibinder plausibly blocks ligand
  binding — the functional relevance the competition asks for
- Conserved, solvent-exposed carboxylates available for histidine pairing: E45 (rSASA
  0.46), D46 (0.36), D75 (0.47), E97 (0.34)
- **Accessible in the assayed state.** Relative solvent accessibility of all seven
  hotspots is identical in our trimmed domain I construct and in the intact tethered ECD
  (e.g. 0.37/0.37, 0.46/0.46), so binders designed against the trimmed target engage the
  same surface on the full receptor.

**Epitope B — domain II dimerisation arm (hedge):** T274, M277, D278, V279, P281, E282,
G283. 100% conserved, 6–7 residues at the dimer interface. **Y275 was deliberately
excluded**: its relative accessibility rises from 0.30 to 0.63 on trimming, showing it is
half-occluded by domain IV in the intact tethered receptor, so targeting it would produce
binders that clash with domain IV on the protein actually assayed.

---

## 2. Design generation

Target constructs: 6ARU chain A with the Fab removed, trimmed to domain I (PDB 1–164) and
domain II (PDB 165–310).

- Generators: `[TBD — RFdiffusion / BindCraft, counts per method]`
- Binder length sampled 55–95 aa, keeping all designs in the **minibinder** category
- Sequence design: ProteinMPNN **soluble** weights, cysteines omitted
- `[TBD]` designs generated, `[TBD]` passing filters

Where multiple generators contributed, no single method accounts for more than half of the
submitted designs, following the diversity constraint used in Anthropic's own campaign
protocol to avoid method-specific failure modes.

---

## 3. Acid-ON pH switch

### 3.1 Mechanism, and how it differs from published work

The state of the art (Marchand et al., bioRxiv 2025.09.29.678932) designs **acid-OFF**
binders: histidines adjacent to *cationic* target residues, so protonation at low pH
creates repulsion and weakens binding. That paper reports the acid-ON direction as theory
only — histidines as hydrogen-bond donors rather than acceptors — with no experimental
examples.

Challenge 1 requires acid-ON. Our mechanism:

- Histidines on the binder are positioned against **conserved target carboxylates**.
- At **pH 6.5** the imidazolium is protonated and forms a charged hydrogen bond / salt
  bridge → binding.
- At **pH 7.4** the neutral imidazole loses that contact, leaving a carboxylate that the
  binder has desolvated but no longer pairs with. The "off" state is a *specific* lost
  contact plus a desolvation penalty, not a diffuse weakening.

We note plainly that this direction is experimentally unproven. That is the reason we
consider it worth testing, and the reason each switch design is paired with a
pH-insensitive control at the same epitope.

### 3.2 Placement rules

Carried from the published work: deliberate placement only (random histidine incorporation
rarely switches); **≥2 independent His–carboxylate contacts** required; buried core
positions never substituted.

### 3.3 Geometric calibration

Every threshold was measured rather than assumed, from **23 genuine His–carboxylate salt
bridges** (N–O ≤ 3.3 Å) across 4HHB, 1IGY, 1BRS, 1A4Y, 2PTC, 1MBN, 5PTI, 1AKI, 3HHB, 1LZ1,
6ARU and 1IVO:

| Quantity | Observed | Used |
|---|---|---|
| Cβ → carboxylate O | 3.55–6.77 Å (median 6.22, mean 5.79) | window 3.5–7.0 Å, optimum 5.8 Å |
| Cα→Cβ direction cosine | −0.33 to 0.81 (p5 −0.26) | threshold ≥ −0.30 |
| N–O contact distance | 2.52–3.27 Å (median 2.84) | — |

This overturned two of our own starting assumptions:

1. A "the side chain must point at the target" filter (cos > 0.2) **rejects 48% of genuine
   salt bridges**, because χ1/χ2 rotation lets the imidazole swing far off the Cα→Cβ axis.
   The real barnase His102–barstar Asp39 bridge sits at cos 0.07.
2. Once a position is reachable, the contact forms at ~2.84 Å **regardless** of Cβ
   distance. Contact strength is therefore near-constant per contact, so we score the
   number and geometric quality of independent contacts rather than reporting a
   distance-dependent pseudo-ΔΔG in kcal/mol it has not earned.

### 3.4 Validation

- **Positive control:** given barnase/barstar (1BRS) with no knowledge of the answer, the
  engine proposes His102→Asp39 — the salt bridge that exists in nature — at rank 2 of 5.
- **Negative control:** the cetuximab–EGFR interface yields no viable position, as section
  1.1 predicts.

---

## 4. Scoring, novelty and ranking

- **Novelty (hard filter):** BLAST against SwissProt; designs retained only at ≤30%
  maximum identity to any known protein, matching the Proteinbase de novo threshold.
  `[TBD — identity range across submitted designs]`
- **Structure/affinity proxies:** `[TBD — pae_interaction, pLDDT per design]`
- **Ranking**, following the challenge's own stated priority rather than ours:
  pH-switch evidence → epitope conservation / cross-species → affinity proxy. The FAQ
  states that a weak but clearly pH-sensitive binder may be more impactful than a
  high-affinity binder that is not pH-sensitive, so switch evidence outranks predicted
  affinity in our ordering.

---

## 5. Submission composition

`[TBD]` designs: `[TBD]` switch designs and `[TBD]` pH-insensitive controls at the same
epitope.

The controls are deliberate. Anthropic's campaign found that only 4 of 8 EGFR binders
cross-reacted with the mouse ortholog, and that footprint substitution counts did not
predict cross-reactivity. Pairing switch designs with pH-insensitive controls at an
identical epitope makes the experiment interpretable whichever way the measurements fall:
a difference between the pair is attributable to the histidines rather than to the
epitope or the scaffold.

---

## 6. Limitations

- The acid-ON mechanism has no published experimental precedent; this is a test of it, not
  a validated design.
- Co-folding models are pH-blind, so no structure-prediction metric here evaluates the
  switch. Switch evidence is geometric and electrostatic, and is reported as such.
- Cross-reactivity was designed in through epitope conservation, which Anthropic's data
  show is necessary but not sufficient — their own analysis found no significant trend
  between footprint substitutions and ortholog binding.
- Compute was limited to free-tier GPU, so the number of trajectories is far below what
  these pipelines normally use. Designs were selected from a smaller pool than is typical,
  and filter thresholds are reported honestly rather than tuned to admit more designs.
- Domain I was trimmed from the receptor for design. We verified hotspot accessibility is
  unchanged in the intact ECD, but the construct remains a model of the assayed protein.
