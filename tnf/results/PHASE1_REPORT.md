# TNF-α Challenge 2 — Phase 1: where to bind

**Validation of the pipeline:** our alignment reproduces the challenge page's stated
figure exactly — soluble domains 77–233 are **124/157 = 79.0% identical**, and the
1TNF PDB→UniProt offset of +76 matches 151 of 152 resolved residues in chain A.

## The three candidate epitopes, measured

| Epitope | residues | human/mouse identity | conserved cationic anchors | mouse substitutions |
|---|---|---|---|---|
| TNFR2 receptor site (3ALQ) | 32 | **71.9%** | 4 — R107, R108, K166, K204 | P96H, R107Q, E129D, S147D, T148–, H149Y, V161I, T165E, I173V |
| Adalimumab (3WD5) | 19 | **84.2%** | 1 — K141 | P96H, A187L, D216K |
| Infliximab (4G3Y) | 28 | **71.4%** | 1 — R214 | G100E, S147D, H149Y, I173V, E180D, A187L, R214L, D216K |

**The recommended epitope is the least conserved of the three.** The receptor-binding
site sits at 71.9%, seven points *below* the 79% domain average, and carries a charge
reversal (R107Q) plus eight other substitutions. Objective 2 is measured against mouse
TNF-α directly, so this is a real cost, not a theoretical one.

## Why cationic residues decide the switch

Challenge 2 inverts Challenge 1: bind at pH 7.4, release at pH 6.0. That is the
**acid-OFF** direction, which unlike acid-ON has published precedent (Marchand et al.).
The mechanism needs a histidine that is neutral at 7.4 and protonated at 6.0, placed
where the new positive charge collides with an existing one — so the epitope must
carry conserved, exposed **Arg or Lys**.

Only **five** such residues exist on the whole human/mouse-conserved surface:
K87, R108, K141, K166, K204 (plus K188, exposed in chains B and C).

Spatial pairs close enough for one binder footprint (≤20 Å):

- **K141–K188: 10.4 Å** — the tightest pair, within one protomer
- K166–K204: 15.7 Å · K87–K204: 13.7 Å · R108–K166: 16.2 Å (across protomers)

## Patch scan

55 surface patches (12 Å radius) scored on conservation, cationic content and overlap
with the receptor footprint:

| centre | n | identity | cationic | TNFR2 overlap | adalimumab overlap | anchors |
|---|---|---|---|---|---|---|
| **K141** | 20 | **80.0%** | K141 + K188 | 7/32 | 11/19 | K141, K188 |
| **S223** | 15 | **86.7%** | 1 | 8/32 | 5/19 | R108 |
| I159 | 15 | 73.3% | 2 | 9/32 | 0 | K166, K204 |
| G116 | 16 | 81.2% | 2 | 2/32 | 0 | K87, K204 |

**No patch is simultaneously >80% conserved, two-anchored and strongly
receptor-blocking.** The trade-off is real and has to be made deliberately:

- **K141 patch** — the best balance. 80% conserved (above domain average), carries the
  tightest cationic pair in the whole structure, and overlaps the adalimumab epitope
  (11 of 19 residues), which is a clinically validated TNF-blocking site.
- **S223 patch** — the most conserved surface (86.7%) with solid receptor overlap, but
  only one cationic anchor, so the switch would rest on a single contact.
- **I159 patch** — the best receptor blocker with two anchors, but 73.3% conservation
  makes objective 2 a gamble.

## Recommendation

Split the 20 designs across the two best patches rather than betting on one:
**K141 (primary)** and **S223 (secondary)**. Challenge 1 showed binders drift from the
intended epitope regardless, so diversity of starting point is cheap insurance, and
the challenge scores pH-selectivity and cross-reactivity above raw affinity.
