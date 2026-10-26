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

---

# REVISION, same day: the mouse check moved the primary site

Everything above ranked patches on human/mouse **sequence identity**. This section adds the
structural check and **changes the recommendation**. Reproduce with:

    python tfr1/mouse_check.py

## Method

Mouse structure is the AlphaFold DB model for Q62351 (v6, the current version - the URL
was read from the API rather than guessed, and v4 and v5 both 404). Human is 11QC chain A.
Global superposition over 638 ectodomain CA pairs gives **RMSD 1.41 A**, so the two
proteins are closely similar overall and per-patch differences mean something.

Each patch is then superposed **on its own CAs alone**, so the RMSD reported measures
whether the patch is the same *shape*, independent of where it sits.

Substitutions are classified by chemistry (+ / - / polar / hydrophobic / special). A swap
inside a class, such as Lys to Arg, keeps the chemistry a binder sees; a swap across
classes does not. "Functional conservation" counts only the across-class ones as changes.

Mouse model confidence (pLDDT) is reported per patch. Everything below is 86-94, so the
model is confident at every site discussed and none of these rows is weak for that reason.

## Results

| patch | n | identity | functional | patch RMSD | pLDDT | substitutions |
|---|---|---|---|---|---|---|
| **R732** | 14 | 86% | **100%** | 0.63 | 94 | K633R, E634D - both inside class |
| **S706** | 17 | 94% | 94% | 0.56 | 92 | P710S |
| E559 | 13 | 85% | 92% | **0.48** | 91 | T561A, R579K |
| G751 | 15 | 93% | 93% | 1.14 | 88 | D755N |
| T729 | 12 | 92% | 92% | 0.93 | 93 | G724T |
| R325 (apical shoulder) | 14 | 86% | 93% | 1.18 | 86 | **R325Q**, K693R |
| E714 | 19 | 84% | 84% | 0.79 | 93 | P710S, K720Q, Q721K |
| apical tip (recommended) | 9 | 33% | **33%** | 1.17 | 91 | D204Q, R208N, Y211D, L212P, N348K, K371Q |

None of these patches is occluded by the other subunit: dimer-context and monomer exposure
agree to within 0.05 everywhere, so the dimer does not bury any of them.

## What changed, and why

**The apical shoulder (R325) is demoted from primary.** It looked good on identity alone
(86%) but the structural check is unkind to it: the centre residue itself is **R325Q**, a
loss of positive charge at the middle of the patch, and it has the worst shape agreement
(1.18 A) and the lowest mouse-model confidence (86) of any conserved candidate. A binder
designed around an arginine at the centre of its epitope has nothing to grip in the mouse.

**The new primary is the R732 patch**: every substitution preserves chemistry class
(Lys to Arg, Glu to Asp), the shape is reproduced to 0.63 A, and the mouse model is
confident there. On what a binder actually contacts, this patch is effectively identical
between the two species.

**S706 stays strong** and is the largest of the good patches at 17 residues, with the
single substitution P710S.

**The apical tip is confirmed as a bad bet for objective 2** and not merely a weak one. Its
six substitutions include two losses of positive charge (R208N, K371Q) and a polar-to-acid
flip (Y211D) - a charge-reversing surface, not a conservative drift.

## Revised plan

| rank | site | why |
|---|---|---|
| primary | **R732 patch** | 100% functional conservation, 0.63 A shape agreement |
| secondary | **S706 patch** | 94%, largest usable patch, 0.56 A |
| third | **E559 patch** | best shape agreement anywhere (0.48 A), 92% |
| diversity | R325 apical shoulder | keeps an apical-domain epitope in the submission; rubric rewards epitope diversity and functional sites |
| declared long shots | apical tip | best functional story, worst cross-species odds; submitted as such |

The first three sit in one broad conserved region of the helical domain, 46 residues
spanning 35 x 32 x 53 A and 85% conserved. That is a large target area, which gives the
design runs room to find their own sub-sites rather than forcing one fixed footprint.

## What this still does not establish

- The mouse structure is **a prediction**. Confident (pLDDT 86-94) but not an experiment.
- **Glycans are still ignored.** TfR1 is glycosylated and sequons differ between species;
  a patch that is open in both models could be shadowed in one and not the other. Not
  checked, and it is the most likely way this analysis is wrong.
- Conservation says a binder *could* cross-react. It does not say one *will*.
- Still nothing here about the pH switch.
