# Phase 1 — Epitope selection for a mouse cross-reactive, pH-selective EGFR binder

**Date:** 2026-10-02 · **Challenge 1**, Anthropic × Adaptyv Protein Design Competition
**Inputs:** UniProt P00533 (human EGFR), Q01279 (mouse Egfr), PDB 6ARU (ECD + cetuximab Fab), PDB 1IVO (ECD + EGF, dimerized)
**Code:** `scripts/phase1_conservation.py`, `scripts/phase1_epitopes.py`, `scripts/phase1_functional.py`

All numbering is UniProt P00533 (mature numbering = UniProt − 24; PDB numbering in both
6ARU and 1IVO = UniProt − 24, verified at 99.7% and 100% sequence match respectively).

---

## 1. Summary

Objective 2 (mouse cross-reactivity) and Objective 3 (acid-ON pH switching) impose
constraints on the *target surface* that must be satisfied before any backbone is
generated. We mapped both constraints onto the structure and identified three candidate
epitopes, none of which is the epitope the challenge page recommends.

**Decision: do not target the cetuximab epitope.** Target the domain I ligand-binding
face (primary), with the domain II dimerization arm as a diversification hedge.

---

## 2. Cross-species conservation

| Region | Human/mouse identity |
|---|---|
| Full length (1210 aa) | 90.5% |
| Extracellular region (25–645) | 88.7% |
| Domain I (L1) | 85.4% |
| Domain II (CR1, dimerization arm) | **93.2%** |
| Domain III (L2, cetuximab epitope) | 87.7% |
| Domain IV (CR2) | 89.3% |

Global identity is high, but binding is decided locally: 48 of 284 candidate surface
patches are **100% identical** between human and mouse, and many others are not.

## 3. The cetuximab epitope is the wrong target

Computed empirically as EGFR residues with any atom within 5 Å of the cetuximab Fab in 6ARU.

**Footprint (26 residues):** P373, V374, R377, L406, Q408, Q432, H433, Q435, F436, A439,
V441, S442, I462, S464, G465, N466, K467, K489, I490, I491, S492, N493, R494, G495, E496, N497

**7 of 26 are substituted in mouse**, and they sit at the closest contacts:

| Substitution | Distance to Fab |
|---|---|
| S492N | 2.82 Å |
| I491M | 3.05 Å |
| N497K | 3.08 Å |
| S442G | 3.26 Å |
| K467R | 3.39 Å |
| G495A | 3.44 Å |
| R377K | 3.61 Å |

This is a structural explanation for the known fact that cetuximab does not appreciably
bind murine EGFR despite ~88% domain III identity — and it is a direct warning that a
binder designed against this surface is likely to fail Objective 2.

**A second, independent disqualifier:** the cetuximab footprint contains exactly **one
carboxylate (E496) in 26 residues**. The acid-ON switch mechanism requires target Asp/Glu
for protonated histidines to pair with. This surface cannot support the mechanism.

*Method validation:* the same footprint overlaps the EGF ligand site at 17/26 residues,
consistent with cetuximab's known mechanism as a ligand-blocking antibody.

## 4. Functional sites, computed empirically

From 1IVO (EGF-bound, dimerized):

- **EGF ligand footprint (43 residues):** S35, N36, K37, L38, T39, Q40, L41, G42, T43, D46,
  R53, Y69, A92, L93, Y113, E114, L122, S123, Y125, L349, H370, L372, P373, V374, D379,
  S380, F381, T382, L406, Q408, Q432, H433, Q435, F436, A439, V441, S442, I462, S464, G465,
  K489, I491, S492 — spanning domain I and domain III, as expected for the bivalent site.
- **Receptor–receptor dimer interface (33 residues):** includes P266, M268, Y270, P272, T273,
  T274, Y275, Q276, M277 — the domain II dimerization arm. (Y275/M277 here are Y251/M253 in
  mature numbering, the canonical arm residues — an independent check that the mapping is correct.)

## 5. Candidate epitopes

Criteria: 100% human/mouse identity across the patch, ≥3 conserved solvent-exposed
carboxylates for histidine pairing, and overlap with a functional site.

### Epitope A — domain I ligand-binding face **(PRIMARY)**
- **Center residues:** F44, E45, D46, L49, Q52, Y74, D75
- **Conservation:** 100% human/mouse across the patch
- **Functional:** 6–7 residues in the EGF footprint → direct ligand blockade
- **His-pairing partners:** E45 (rSASA 0.46), D46 (0.36), D75 (0.47), E97 (0.34)
- **Accessibility:** the domain I ligand surface is solvent-exposed in the tethered
  (autoinhibited) conformation, which is the state the isolated ECD predominantly adopts
  in solution — i.e. the state that will be on Adaptyv's chip.
- **Hotspots (PDB numbering, −24):** 20, 21, 22, 25, 28, 50, 51

### Epitope B — domain II dimerization arm **(SECONDARY / hedge)**
- **Center residues:** T274, Y275, M277, D278, V279, P281, E282, G283
- **Conservation:** 100%; domain II is the most conserved domain (93.2%)
- **Functional:** 6–7 residues in the receptor–receptor dimer interface → blocks activation
- **His-pairing partners:** D278 (0.20), E282 (0.62), D321 (0.79), E602 (0.43)
- **Risk:** the arm is partly occluded by domain IV in the tethered conformation — our own
  patch analysis picks up domain IV residues (600, 602, 603, 604) as spatial neighbours,
  confirming the packing. Binding to the isolated ECD may therefore be weak.
- **Hotspots (PDB numbering, −24):** 250, 251, 253, 254, 255, 257, 258, 259

### Epitope C — domain II K294 face **(BACKUP)**
- **Center residues:** K294 patch (257, 267, 282, 283, 286, 290, 293, 294, 296–298, 308, 319–322, 324)
- **Best carboxylate density of any conserved patch:** E257 (0.68), E282 (0.62), E319 (0.56),
  E320 (0.32), D321 (0.79) — five conserved, well-exposed His-pairing partners
- **Weakness:** only 2 residues at the dimer interface, so the functional claim is thin
- Use only if A and B both fail to produce designs with a credible switch

### Rejected
- **Cetuximab/domain III face** — 7/26 mouse substitutions at the closest contacts; 1 carboxylate in 26 residues
- **Domain III K454/S452 patches** — 100% conserved and carboxylate-rich, but they straddle the
  domain III/IV interface (members include 513, 521, 524, 525, 533), making them
  conformation-dependent and therapeutically harder to justify

---

## 6. Design allocation (20 submission slots)

| Epitope | Designs | Rationale |
|---|---|---|
| A — domain I ligand face | 11–12 | Primary: conserved, functional, accessible |
| B — domain II dimer arm | 5–6 | Diversification against conformational risk |
| Controls | 3–4 | pH-insensitive binders at epitope A, no interfacial His |

Rationale for controls: Anthropic found that on EGFR only 4 of 8 binders cross-reacted with
mouse, and that footprint substitution counts did not predict cross-reactivity. A paired
switch/control set at the same epitope makes the experiment interpretable whichever way it
goes, and interpretability is part of what the Track 3 selection prompt rewards.

---

## 7. Carried into Phase 2

1. Build target constructs from **6ARU chain A with the Fab stripped** (tethered ECD, the
   assay-relevant state). Trim to a design construct around each epitope to fit free-tier GPU memory.
2. Verify in the stripped construct that the epitope A hotspots remain solvent-exposed.
3. Hotspot lists above go straight into BindCraft / RFdiffusion as target hotspot residues.
4. Carry the carboxylate lists (E45, D46, D75, E97 for A; D278, E282, D321, E602 for B) into
   Phase 3 as the His-pairing partners for the acid-ON switch.
