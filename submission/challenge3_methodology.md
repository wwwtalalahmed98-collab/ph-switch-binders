## Summary

20 single-chain de novo designs against human GLP1R (PDB 6X18, chain R). Four are
agonist candidates. The other sixteen are binders that, by mechanism, are blockers — and
this submission says so rather than presenting them as agonists.

**No design has been tested for activation.** Nothing in our pipeline can test it. What we
screened is the geometry that the known activation mechanism requires.

## The problem as we understood it

GLP-1 does two things at once. Its C-terminal half grips the extracellular domain groove;
its first ~8 residues thread ~33 Å down into the 7TM bundle, and that threading is what
opens the receptor. Grip without thread gives a blocker — which already exists as a drug,
exendin(9-39), so reproducing it wins nothing. Thread without grip has nothing holding it
on.

Our first approach split these into two molecule classes, groove-binders and
channel-threaders. That was wrong, and the record of why is the most useful thing we
produced.

## The mistake that taught us the most

Channel-threaders failed: 1 of 11 usable. The cause was our own target. To make the
receptor affordable to compute against we had cropped it to a shell around the binding
sites, which removed the floor of the transmembrane bundle. The designs were diving into
the hole we had made. At the depth they reached, the real receptor has **211 atoms** and
our crop had **4**.

We rebuilt the target with the channel floor sealed, re-docked the native hormone into it
first to confirm it reproduced published contacts, and re-ran the identical job as a
controlled experiment — same method, one variable changed. Channel-only designs failed
again, 4 of 7 touching zero channel residues. So the method was at fault, not only the
broken crop. That class received zero reserved slots here.

The design that works is one chain doing both jobs, which is what the hormone does.

## The finding we would submit even if every design fails

**On this target, AlphaFold2's interface confidence is anti-correlated with landing on
the right residues.**

- −0.963 across 11 designs (first target)
- −0.651 across 15 designs (sealed target)
- −0.804 among two-domain designs alone

Two targets, two design tasks, the same inverted sign. In the final batch the
highest-confidence design touches **1 of 30** channel residues; the lowest touches **19**.

Ranking by i_ptm is the obvious move and here it selects precisely the wrong designs. We
therefore applied **no confidence floor** to the channel-reaching populations and selected
on direct geometric measurement instead. The four agonist candidates are consequently the
four lowest-confidence designs in this submission (i_ptm 0.34–0.48). That is deliberate.
If they are filtered out on confidence, the submission loses the only designs that could
possibly activate anything.

## Method

AlphaFold2 hallucination (ColabDesign/AfDesign), three-stage logits → softmax → hard. No
design is derived from, threaded onto, or initialised from any natural sequence.

Every reported number comes from an **independent AF2 re-prediction at 3 recycles** of the
design–target complex, run after design. The optimiser's own score was recorded separately
and never used to select.

Epitope: two contact sets read off 6X18 rather than chosen by the optimiser — 22 groove
(ECD) residues and 30 channel (TM) residues. Cross-checked against 3IOL: 16 of its 17
epitope residues are in our set, Jaccard 0.70, zero residue-identity mismatches. Our set
adds the ECD C-terminal stalk (P137, E138, L141, L144, Y145), which does not exist in
3IOL and so is invisible to anyone working from that structure.

Native GLP-1 measured in the same target as a reference: reach 33.2 Å, 26 on-site
contacts, i_ptm 0.761.

## Selection

| population | n | rule | when written |
|---|---|---|---|
| landing-pad (groove) | 13 | i_ptm ≥ 0.50, groove ≥ 8 of 22, hydropathy ≤ 0.0 | before any design existed |
| two-domain (agonist candidates) | 4 | groove ≥ 4, channel ≥ 5, N-terminus leading, hydropathy ≤ 0.5 | after seeing the designs |
| dual-site reversed | 3 | groove ≥ 5, channel ≥ 15, C-terminus leading | after seeing the designs |
| all three | — | buildability: ≤ 6 cysteines, no run > 6 identical residues | after seeing the designs |

Only the first rule is pre-registered. The other three were set after seeing the data and
are weaker for it; we are marking them rather than presenting them as equivalent.

The two-domain pool filled 4 of 5 slots. Rather than pad it with a design that failed its
own filter, the spare went to the groove pool. `glp1r_dual_02` originated in the
channel-only class but qualified on its measurements, so it is included on geometry, not
on class.

## What this submission does not establish

- **Activation is unscreened.** The central limitation.
- No experimental result exists for any design. Every number is a prediction or a
  geometric measurement on a predicted structure.
- Foldability is not filtered. `glp1r_dual_03` matches the native hormone's 33.2 Å depth
  exactly and has helicity 0.14, pLDDT 0.458. It may not fold. We kept it because cutting
  on foldability removes the best geometric match to the native mechanism.
- The buildability gate was a problem we had not screened for at all until the final day;
  it is not applied to our earlier EGFR and TNF-α submissions, which are closed.
- 16 of the 20 are blockers by mechanism. They are submitted as binders, not as agonists.
