# Where to bind TfR1, and why

Computed 2026-10-26, before any design work. Reproduce with:

    python tfr1/conservation_map.py
    python tfr1/find_patches.py

## The question

Challenge 4 needs one binder that works on human **and** mouse TfR1. The organizers
recommend the apical domain. They also state the apical patch differs at six of nine
positions while the ectodomains are 78% identical. So the recommended site is the worst
place on the molecule to meet the cross-species requirement, and we needed to know what
the alternatives actually look like rather than assume.

## Two checks that the pipeline is reading the structures correctly

Both numbers below were produced by our code and then compared against figures the
organizers publish. Neither was tuned to match.

| our result | their stated figure |
|---|---|
| ectodomain 121-760 is **78.1%** identical human/mouse | 78% |
| apical epitope from 11RN is **3 of 9** conserved | "differs at six of nine positions" |

Residue numbering was checked against UniProt P02786 across 1278 residues in 11QC with
**zero mismatches** before any geometry was used.

## A bug found and fixed before it misled us

The first run collapsed the homodimer by keeping whichever copy of each residue was more
exposed. That mixed chain A and chain B coordinates into one coordinate set, so residues
50-100 A apart on opposite subunits landed in the same "patch". Every spatial number from
that run was meaningless; it was caught because the patch counts came out absurdly low
(6 patches from 253 exposed residues).

Fixed: exposure is still computed on the full dimer, so the subunit-subunit interface is
correctly buried, but all geometry now comes from chain A alone.

## What the surface actually looks like

Exposed residues on one subunit: 247. Patch conservation across 213 candidate patches
(14 A radius, >=12 exposed members): median 67%, 90th percentile 86%, max 100%.
Background ectodomain identity is 78%.

Sites that must be avoided were read off real complexes rather than remembered:
1SUV (transferrin-receptor) and 1DE4 (HFE-receptor) give a combined 45-residue footprint.

Height is distance from residue 122, the membrane-proximal end of the supplied construct.
TfR1 is a type II membrane protein, so low numbers sit near the membrane.

| site | residues | height A | conserved | clear of Tf site |
|---|---|---|---|---|
| apical tip (VP1u, 11RN) - **the recommended one** | 9 | 73 | **33%** | yes |
| transferrin / HFE footprint | 45 | 33 | 82% | n/a, off limits |
| **S706 patch** | 17 | 50 | **94%** | 17 A |
| **T729 / L730 / R732 patch** | 12-15 | 55-57 | 86-92% | 11-16 A |
| **R325 / E272 / A269 / S327 patch** | 11-14 | 50-58 | 80-86% | 18-30 A |
| K495 patch | 14 | 34 | 93% | 9 A - tight |
| G751 patch | 15 | 33 | 93% | 9 A - tight |

## The decision

**Primary site: the conserved shoulder of the apical domain** (the R325 / E272 / A269 /
S327 cluster, roughly residues 240-331).

This is the result worth having. The apical domain spans about 189-383. Its *tip* - the
nine residues parvovirus VP1u touches - is only 33% conserved. Its *shoulder*, 10-20 A
away, is 80-86% conserved. So we do not have to choose between "the functional epitope the
organizers recommend" and "a site where mouse binding is plausible". We can take the
conserved part of the recommended domain: still apical, still leaves the transferrin site
free, still the region brain-shuttle antibodies engage, but roughly two and a half times
more conserved than the exact tip.

**Secondary site: the S706 patch** (94% conserved, 17 residues, 50 A up, 17 A clear of the
transferrin footprint). The most conserved usable surface on the molecule. Less obviously
"functional" than the apical domain, which costs us on objective 3, but it is the strongest
possible bet on objective 2.

**A few designs at the apical tip itself**, as declared long shots. They are the best bet
for affinity and functional relevance and the worst bet for mouse cross-reactivity, and
the selection rubric explicitly values epitope diversity.

Rejected: K495 and G751, both ~93% conserved, both only 9 A from the transferrin
footprint. A miniprotein interface is wider than 9 A, so a binder centred there would
likely overlap the transferrin site.

## What this does not establish

- **Conservation is necessary, not sufficient.** Identical residues do not guarantee
  identical backbone or side-chain conformations in the mouse protein, and the residues
  just outside each patch still differ. A 94% conserved patch can still present
  differently.
- **There is no mouse TfR1 structure in this analysis.** Every mouse statement here comes
  from sequence alignment onto a human structure. A predicted mouse model, and a check
  that the patch geometry survives, is still to do.
- **Exposure is computed on a crystal form**, 11QC, without glycans. TfR1 is glycosylated;
  a patch that looks open here may be shadowed by sugar in the real assay. Not checked.
- Nothing here says anything about the pH switch. That is a separate piece of work and it
  is not written yet.
