# The switch stage found zero carboxylates - and it was a numbering bug

## Symptom

The first stage-3 run reported, for all 14 designs with i_ptm >= 0.45:

    r2_w0_002: 0 candidate contacts, 0 carboxylates -> control only

Zero *candidate contacts*, not merely too few carboxylates. Round 1 produced the
identical symptom on its own 13 backbones and it was written off then as a
consequence of the designs not binding. That explanation was wrong, and the same
bug has now been present through two rounds.

## Why it could not be biology

Three of the four target carboxylates - E45, D46, D75 - are hotspot residues the
binders were explicitly designed against (`hotspot="20,21,22,25,28,50,51"`). A
design with i_ptm 0.786 against those hotspots cannot have *no* residue within
3.5-7.0 A of any of them. The result was impossible, which is what made it a bug
rather than a finding.

## Cause

The target construct is chain A of 6ARU, residues 4-164 - 161 residues whose
numbering **starts at 4**, confirmed locally:

    targets/EGFR_A_domainI_ligand.pdb  chain A  n=161  range 4-164
    at 21,22,51,73: GLU, ASP, ASP, GLU        <- correct in the source

AfDesign renumbers each chain of the complex it writes from 1. In the saved file,
residue 21 is therefore the source's residue 24, not E45 - a three-residue shift.
The placement code looked for ASP/GLU at 21, 22, 51, 73 of the *saved* file, found
none, built an empty `acids` list, and so generated an empty candidate list. Every
design then fell through to "control only".

## Fix

Stop trusting the numbers. Map by position in the chain:

    _src_num = [r.id[1] for r in source_chain]      # [4, 5, ..., 164]
    ACID_IDX = [_src_num.index(n) for n in ACIDS]   # positions, not numbers
    acids    = [tgt_res[i] for i in ACID_IDX]       # same order in the saved file

with an assertion that each mapped residue really is ASP or GLU, and a length check
that the saved target chain matches the source. A silent empty list was the thing
that let this survive two rounds; it now raises.

## Added at the same time

A high i_ptm says the interface is confident, not that it is the interface we asked
for, and the cross-reactivity and pH-switch rationale both rest on it being the
domain I ligand face. Stage 3 now reports, per design, which target residues lie
within 5 A of the binder (in UniProt numbering) and how many fall on the Phase 1
patch F44/E45/D46/L49/Q52/Y74/D75. If the binders are docking somewhere else, that
needs to be known before anything is submitted, not after.
