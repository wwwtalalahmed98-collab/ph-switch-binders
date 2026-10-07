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

---

# Target construct decision

Adaptyv assays the **trimer**, but designing against all three chains (456 residues)
costs roughly five times as much AlphaFold compute per trajectory as Challenge 1's
target, which would leave about nine designs inside the available GPU budget. Two
chains (304 residues) costs about 2.7x and allows roughly eighteen. Dropping a chain
exposes surface that is buried in the real trimer, so the question is whether the
chosen epitope depends on any of it.

Measured residue by residue, rel_SASA in each construct vs the intact trimer:

| construct | residues | artificially exposed (Δ rel_SASA > 0.25) |
|---|---|---|
| A + B (drop C) | 304 | 16 |
| A alone | 152 | 18 |

Against the candidate patches:

| patch | size | artefacts in A+B | artefacts in A alone |
|---|---|---|---|
| I159 | 15 | 0 | 0 |
| S223 | 15 | 1 | 1 |
| **K141** | 20 | **4** | 5 |

**This overturns the K141 recommendation.** A fifth of that patch only exists because
the trimer was cut, and the binder would be free to exploit it — surface that is not
there in the assay. The earlier recommendation was made on conservation and switch
chemistry alone, before the construct was chosen, and it does not survive the
construct.

## The site actually chosen

Scanning composite sites centred between conserved cationic residues on *different*
protomers — the geometry the challenge recommends — gives a clear winner:

**K166(chain A) + R108(chain B), 16 Å apart**

| property | value |
|---|---|
| residues | 25 |
| human/mouse identity | 76.0% |
| conserved cationic anchors | 2 (K166, R108) |
| trim artefacts in A+B | **0** |
| TNFR2 footprint overlap | **13 of 32** |

It beats every single-chain patch on receptor overlap (13 vs 7–9), carries the two
switch anchors the acid-OFF mechanism requires, is free of construct artefacts, and
sits at 76.0% identity — below the 79% domain average but above the 71.9% of the
receptor site as a whole. It is also precisely the inter-protomer receptor site the
challenge page recommends, so no divergence needs defending this time.

Site residues (UniProt): 95, 105–112, 155, 159, 161–168, 201, 203, 207, 220, 221, 223.

Construct: `targets/TNF_AB_receptor_site.pdb`, chains A and B, 152 residues each.
Hotspots: `A86,A87,A88,A90,A91,A92,B29,B32,B33,B34,B35,B36` (PDB numbering = UniProt − 76).
