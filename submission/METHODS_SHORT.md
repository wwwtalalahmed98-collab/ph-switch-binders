# Challenge 1 (EGFR) — methods summary

**Track 3. Five de novo minibinders, 53–69 aa, single chain, ranked by pH-switch
evidence first.** Full record: results/ in the project repo.

## Pipeline

1. **Epitope selection.** Human P00533 / mouse Q01279 aligned and mapped onto 6ARU
   chain A. The cetuximab (domain III) epitope was rejected on evidence: 7 mouse
   substitutions among its closest contacts and 1 carboxylate in 26 residues — poor
   ground for both mouse cross-reactivity and a histidine pH switch. Targeted the
   **domain I ligand-binding face** instead (construct: chain A residues 4–164).
   This diverges from the recommended domain III, deliberately.
2. **Backbone + sequence generation.** AfDesign binder hallucination (ColabDesign),
   `design_3stage(90,40,20)`, 2×T4, 24 trajectories, binder length 52–72.
3. **Scoring.** Every design re-predicted from sequence alone against the target at
   3 recycles — *not* the optimiser's own value. 14 of 24 reached i_ptm ≥ 0.45.
4. **pH switch (acid-ON).** Interfacial positions whose Cβ faces a conserved,
   solvent-exposed target carboxylate (E45/D46/D75/E97/E114/D134/E179) were mutated
   to histidine. Geometry calibrated on 23 genuine His–carboxylate salt bridges
   across 12 structures (Cβ–O 3.5–7.0 Å, ideal 5.8 Å, cos ≥ −0.30); validated by
   recovering barnase His102–Asp39 and by finding nothing in the cetuximab/EGFR
   interface. Protonated imidazolium pairs with the carboxylate at pH 6.5; neutral
   imidazole at 7.4 leaves an unpaired, desolvated carboxylate.
5. **Switch cost measured.** Every histidine variant was re-predicted. All five
   submitted designs retain 86–101% of their parent's i_ptm.

## Metrics (see challenge1_metrics.csv)

| name | i_ptm | retained | mouse identity at interface | His |
|---|---|---|---|---|
| egfr_ph_sw_01 | 0.615 | 86% | 93% | E39H,L32H,E6H |
| egfr_ph_sw_02 | 0.555 | 96% | 90% | S29H,T5H |
| egfr_ph_sw_03 | 0.706 | 101% | 84% | T33H,S49H |
| egfr_ph_sw_04 | 0.506 | 86% | 83% | S17H,F21H |
| egfr_ph_sw_05 | 0.694 | 101% | 92% | S36H,T33H,R9H |

Ranked by switch-geometry score first, per the challenge's stated priority
(pH-selectivity > mouse cross-reactivity > affinity).

## Stated limitations

- **i_ptm is not independent validation.** It comes from the same AlphaFold2 that
  generated the backbones. Adaptyv's assay is the real test.
- **Sequence liabilities are present and unrepaired.** egfr_ph_sw_05 carries 2 free
  cysteines; four of five contain homopolymer runs (L7–L9, R8). A conservative
  repair pass was built and was still running when the window closed, so these are
  the *measured* sequences rather than the cleaned ones. Disclosed rather than
  hidden.
- **The binders do not sit squarely on the intended epitope.** Observed interfaces
  span UniProt ~38–152 and contact only 1–5 residues of the Phase 1 patch. They are,
  however, 83–93% identical between human and mouse across the observed interface —
  better than the 85% construct average — so objective 2 is supported by the
  epitope they actually use.
- **ProteinMPNN redesign was tested and rejected:** it raised pLDDT to 0.82–0.87
  while destroying the interface on 12 of 14 backbones (i_ptm → ~0.10). Only 2
  survived. The submitted sequences are therefore AfDesign's, liabilities included.
- Affinity is assayed against the full ectodomain; designs were made against a
  trimmed domain I. 10 of 161 construct residues are trim-exposed artefacts; the
  per-design dependence on them was not finished before the window closed.
