# ProteinMPNN redesign: mostly destroys the interface

Stage 5 ran ProteinMPNN (soluble weights, cysteines forbidden) on all 14 round-2
backbones, 8 sequences each, every one re-predicted at 3 recycles.

| backbone | AfDesign i_ptm | MPNN best i_ptm | MPNN pLDDT |
|---|---|---|---|
| r2_w0_007 | 0.718 | **0.803** | 0.801 |
| r2_w0_009 | 0.695 | **0.766** | 0.667 |
| r2_w0_013 | 0.748 | 0.217 | 0.587 |
| r2_w0_002 | 0.786 | 0.136 | 0.827 |
| r2_w0_004 | 0.498 | 0.129 | 0.720 |
| r2_w0_001 | 0.501 | 0.123 | 0.659 |
| r2_w0_003 | 0.678 | 0.116 | 0.868 |
| r2_w1_004 | 0.592 | 0.115 | 0.791 |
| r2_w0_006 | 0.701 | 0.111 | 0.550 |
| r2_w1_009 | 0.560 | 0.109 | 0.744 |
| r2_w0_005 | 0.560 | 0.106 | 0.592 |
| r2_w1_006 | 0.557 | 0.104 | 0.865 |
| r2_w0_010 | 0.684 | 0.102 | 0.823 |
| r2_w1_001 | 0.580 | 0.097 | 0.845 |

**12 of 14 collapsed.** Two improved.

The pLDDT column is the tell. MPNN raised it to 0.82-0.87 on several backbones whose
interface it had just destroyed - it is producing confidently-folded proteins that do
not bind. That is MPNN behaving exactly as designed: it optimises sequence-structure
compatibility for the monomer, and nothing in its objective preserves an interface it
was not told about. The AfDesign sequence, for all its ugliness, was selected *for*
the interface.

Placing histidines on the surviving MPNN sequences then failed too: r2_w0_007 went
0.803 -> 0.091 (11%). The other four "KEEP" verdicts in that run are artefacts of a
ratio taken against a parent that was already ~0.10 - 222% of nothing is nothing.

## Consequence for the submission

Wholesale redesign is the wrong instrument here. What the AfDesign sequences need is
not to be rewritten but to have a short list of specific local liabilities fixed:
free cysteines (9 of 14) and homopolymer runs (11 of 14). Stage 6 does exactly that -
conservative single-residue substitutions, switch histidines protected, every repair
re-predicted so its cost is measured rather than assumed.

Two MPNN sequences are worth keeping on their own merits as high-affinity,
liability-free entries: **r2_w0_007 (i_ptm 0.803, pLDDT 0.801, no cysteines, no runs)**
and **r2_w0_009 (0.766)**. They carry no pH switch, so they rank below working switch
designs under the challenge's own priority order, but they are the cleanest binders
in the set.
