# De novo pH-switchable protein binders

Work for the **Anthropic × Adaptyv Protein Design Competition 2026**, Track 3.

Two challenges, both asking for a binder that works in one pH environment and not
another, in opposite directions:

| | Challenge 1 — EGFR | Challenge 2 — TNF-α |
|---|---|---|
| Target | EGFR extracellular region (PDB 6ARU) | TNF-α soluble trimer (PDB 1TNF) |
| pH requirement | binds at **6.5**, not at 7.4 (acid-ON) | binds at **7.4**, not at 6.0 (acid-OFF) |
| Cross-species | mouse EGFR (Q01279) | mouse TNF-α (P06804) |
| Status | 5 designs submitted | in progress |

## What is here

```
scripts/      epitope selection, pH-switch placement, novelty screen, submission assembly
kaggle/       the notebook cells that ran on GPU (generation, scoring, switch, repair)
results/      reports, measurements, and the record of what went wrong
tnf/          Challenge 2: target analysis and epitope selection
targets/      trimmed target constructs used for design
submission/   the ranked CSV and metrics table
```

Structures and sequences are downloaded by the scripts rather than committed.

## Method, in brief

1. **Epitope selection before design.** Cross-species identity and the chemistry
   needed for a pH switch are mapped onto the target first, because they constrain
   the *surface*, not the binder. For EGFR this meant declining the recommended
   epitope: the cetuximab footprint carries 7 mouse substitutions among its closest
   contacts and one carboxylate in 26 residues.
2. **Backbone and sequence generation** with AfDesign binder hallucination
   (ColabDesign), two GPUs, 150 optimisation steps per trajectory.
3. **Independent scoring.** Every design is re-predicted from sequence alone at 3
   recycles. The optimiser's own score is not treated as a measurement — see below.
4. **pH-switch placement** on calibrated geometry, with every histidine variant
   re-predicted so the cost of the switch is measured.
5. **Targeted liability repair** — free cysteines and homopolymer runs substituted
   conservatively, switch histidines protected, each repair re-scored.

## Findings worth reusing

- **A design tool's own score is not a measurement.** Round-1 designs reporting
  i_ptm 0.725 fell to 0.100 when re-predicted independently. The optimiser was
  reporting the objective it was maximising, on sequences that had not converged.
- **pLDDT does not predict interface quality** here (r = 0.048). Ranking on fold
  confidence would have been meaningless.
- **A confident interface is not necessarily the intended one.** Designs scored
  i_ptm 0.69–0.79 while contacting only 1–5 of the 7 targeted residues.
- **ProteinMPNN redesign destroyed the interface on 12 of 14 backbones** (i_ptm to
  ~0.10) while raising pLDDT to 0.82–0.87. Its objective is monomer quality; nothing
  in it preserves an interface. Targeted repair kept 84–103% of binding instead.
- **Silent failures survive.** A residue-numbering mismatch between the reference
  construct (numbered from 4) and the design output (renumbered from 1) returned an
  empty result rather than an error, and went unnoticed for two rounds.

Details, including the measurements behind each, are in `results/`.

## Licence

Code under MIT. Design sequences and experimental results released through
Proteinbase are under ODC-BY per the competition terms.
