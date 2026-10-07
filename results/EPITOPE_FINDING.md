# The binders do not land on the Phase 1 epitope

The numbering fix worked - the log confirms the carboxylates are now resolved
correctly:

    source target: 161 residues, 4-164; carboxylates [21,22,51,73]
    at chain positions [17,18,47,69] (['GLU','ASP','ASP','GLU'])

and the switch placement still returned almost nothing: 0 switch, 14 control. The
reason is the epitope report, which is the most consequential result so far.

| design | i_ptm | observed epitope (UniProt) | n res | on the Phase 1 patch |
|---|---|---|---|---|
| r2_w0_002 | 0.786 | 38-152 | 19 | 1  (F44) |
| r2_w0_013 | 0.748 | 32-151 | 25 | 1  (F44) |
| r2_w0_007 | 0.718 | 36-152 | 27 | **5** (44,45,49,74,75) |
| r2_w0_006 | 0.701 | 38-152 | 19 | 1  (F44) |
| r2_w0_009 | 0.695 | 38-152 | 18 | 1  (F44) |
| r2_w0_010 | 0.684 | 38-181 | 25 | 2  (44,46) |
| r2_w0_003 | 0.678 | 38-181 | 19 | 1  (F44) |
| r2_w1_004 | 0.592 | 38-184 | 24 | 1  (F44) |
| r2_w1_001 | 0.580 | 38-152 | 21 | 2  (44,45) |
| r2_w0_005 | 0.560 | 28-62  | 13 | 0 |
| r2_w1_009 | 0.560 | 28-108 | 15 | 0 |
| r2_w1_006 | 0.557 | 38-152 | 18 | 1  (F44) |
| r2_w0_001 | 0.501 | 29-129 | 18 | **5** (44,45,49,52,74) |
| r2_w0_004 | 0.498 | 69-114 | 3  | 0 |

The hotspot specification was partly honoured - F44 (hotspot 20) is contacted by 11
of 14 - but the interfaces then sprawl across a far wider surface than the ligand
face, and the specific carboxylates E45/D46/D75 are mostly not engaged. A high i_ptm
says the interface is confident, not that it is the interface that was asked for.
Only r2_w0_007 and r2_w0_001 bind substantially where Phase 1 intended.

## This breaks the pH-switch plan as written, but not the science behind it

Two measurements decide whether this is recoverable.

**Mouse cross-reactivity survives.** The observed contact span (UniProt 38-152) is
105/115 = **91.3% identical** between human P00533 and mouse Q01279 - *better* than
the 85.1% average over the whole construct. Requirement 2 is not damaged by the
binders having moved; they moved onto surface that is, if anything, more conserved.

**There are more carboxylates than the four that were whitelisted.** Applying the
original selection rule to the whole construct - Asp or Glu, identical in mouse,
rel_SASA >= 0.20 - gives seven, not four:

    UniProt  45, 46, 75, 97, 114, 134, 179
    PDB      21, 22, 51, 73,  90, 110, 155

E114, D134 and E179 were excluded only because they sit outside the Phase 1 patch.
Designs whose interfaces reach to 152 and 181 had no eligible carboxylate in the
whitelist, which is the direct cause of "1 candidate contact, 1 carboxylate" on the
better designs and zero on the rest.

## What stage 4 changes

Nothing about the method, only the scope of its inputs: all seven conserved exposed
carboxylates are now eligible, and each design additionally reports the conservation
of *its own* interface with the specific mouse substitutions named. The epitope was
always a means to an end - conserved surface, carboxylates to protonate against -
and the end is what the challenge actually scores.
