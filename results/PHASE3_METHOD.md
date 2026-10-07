# Phase 3 — acid-ON pH switch engineering: method and validation

**Status:** pipeline built and validated against ground truth. Awaiting Phase 2 designs.
**Code:** `scripts/phase3_ph_switch.py` (`--selftest` runs both controls)

## Mechanism

Inverted relative to the published pH-switch literature. Marchand et al.
(bioRxiv 2025.09.29.678932) design **acid-OFF** binders: histidines next to *cationic*
target residues, so protonation at low pH creates repulsion and weakens binding. They
report the acid-ON direction only as theory, with no experimental examples.

Challenge 1 requires acid-ON. Our mechanism:

- Place histidines on the binder against **conserved target carboxylates** (Asp/Glu).
- At **pH 6.5** the imidazolium is protonated and forms a charged hydrogen bond / salt
  bridge → binding.
- At **pH 7.4** the neutral imidazole loses that contact, leaving a carboxylate that the
  binder has desolvated but no longer pairs with → the "off" state is a *specific* lost
  contact plus a desolvation penalty, not a diffuse weakening.

Design rules carried over from that paper: deliberate placement only (random histidine
incorporation rarely switches); require **≥2 independent His–carboxylate contacts**; never
substitute buried core positions.

## Geometric calibration

Every threshold was measured rather than assumed, from **23 genuine His–carboxylate salt
bridges** (N–O ≤ 3.3 Å) across 4HHB, 1IGY, 1BRS, 1A4Y, 2PTC, 1MBN, 5PTI, 1AKI, 3HHB, 1LZ1,
6ARU and 1IVO:

| Quantity | Observed range | Used |
|---|---|---|
| Cβ → carboxylate O | 3.55–6.77 Å (p5 3.83, median 6.22, mean 5.79) | window 3.5–7.0 Å, optimum 5.8 Å |
| Cα→Cβ direction cosine | −0.33 to 0.81 (p5 −0.26, median 0.27) | threshold ≥ −0.30 |
| N–O contact distance | 2.52–3.27 Å (median 2.84) | — |

**Two assumptions this overturned:**

1. *"The side chain must point toward the target."* A cos > 0.2 filter — which is what
   intuition suggests — **rejects 48% of genuine salt bridges**, because χ1/χ2 rotation
   lets the imidazole swing far off the Cα→Cβ axis. The real barnase His102–barstar Asp39
   bridge has cos = 0.07. The threshold is now set at the measured p5.
2. *"Score contacts with a distance-dependent electrostatic energy."* The calibration shows
   that once a position is reachable, the contact forms at ~2.84 Å **regardless** of how
   far the Cβ sits. Contact strength is therefore near-constant per contact, and a
   distance-dependent pseudo-ΔΔG would have been fitting noise and reporting it in
   kcal/mol it had not earned. Scoring is now the number and geometric quality of
   independent contacts, weighted by how much the binder desolvates each carboxylate —
   which is also the variable the published work found determinant.

## Validation

**Positive control — barnase/barstar (1BRS).** Barnase His102 forms a salt bridge with
barstar Asp39 (N–O 2.81 Å). Given the complex with no knowledge of the answer, the engine
proposes H102→D39 at rank 2 of 5 (contact quality 0.78). **Recovered: YES.**

**Negative control — cetuximab Fab/EGFR (6ARU).** No viable positions, as predicted in
Phase 1: the cetuximab epitope carries one carboxylate in 26 residues, and the closest
approach from any Fab Cβ to an EGFR carboxylate is 8.21 Å — outside the salt-bridge window
entirely.

That negative result is a third independent line of evidence against the epitope the
challenge page recommends, alongside its 7/26 mouse substitutions and its lack of
carboxylates. The surface is not merely a poor choice for cross-reactivity; it is
geometrically incapable of supporting this switch mechanism.

## What runs when designs arrive

```bash
python scripts/phase3_ph_switch.py <complex.pdb> --binder B --target A --carboxylates 21,22,51,73
```

(`--carboxylates` in PDB numbering: E45, D46, D75, E97 at epitope A → 21, 22, 51, 73.)

Output per design: ranked histidine substitution proposals, the count of independent
reachable carboxylates, and a PASS/FAIL on the ≥2-contact rule. Designs that fail that rule
are poor switch candidates regardless of their predicted affinity, and should be kept only
as pH-insensitive controls.
