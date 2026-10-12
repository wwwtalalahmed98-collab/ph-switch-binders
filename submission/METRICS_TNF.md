# Challenge 2 — per-design metrics

Twenty designs in submitted order. Rows 1-10 are pH-switch candidates;
rows 11-20 are their matched histidine-free controls, each in its partner's position.

`i_ptm` is AlphaFold2 interface pTM from an independent re-prediction of the
discrete sequence at 3 recycles, not the design optimiser's own score.
`switch_strength` sums `(8.8 - d)/(8.8 - 3.0)` over histidines within reach,
where d is the measured C-beta to cation-nitrogen distance in the confirmed structure.
`mouse_functional_switch` is *no* where every histidine faces R107, which is
glutamine in mouse TNF-alpha and therefore presents no cationic partner.

## Switch candidates

| name | pair | len | i_ptm | His | functional | switch str. | mouse switch | epitope cons. | His contacts (Cb, A) |
|---|---|---|---|---|---|---|---|---|---|
| tnf_sw_01 | tnf_ct_01 | 73 | 0.622 | 2 | 2 | 1.45 | yes | 82% | R108B 5.16 / R108B 4.01 |
| tnf_sw_02 | tnf_ct_02 | 73 | 0.725 | 3 | 3 | 1.51 | partial | 81% | R108B 4.92 / K166A 8.30 / R107B 4.41 |
| tnf_sw_03 | tnf_ct_03 | 70 | 0.487 | 2 | 2 | 1.26 | yes | 79% | R108B 4.74 / R108B 5.58 |
| tnf_sw_04 | tnf_ct_04 | 79 | 0.600 | 3 | 3 | 1.87 | no | 73% | R107B 4.35 / R107B 8.11 / R107B 3.11 |
| tnf_sw_05 | tnf_ct_05 | 75 | 0.487 | 2 | 2 | 0.95 | yes | 83% | R108B 7.03 / R108B 5.06 |
| tnf_sw_06 | tnf_ct_06 | 63 | 0.544 | 3 | 2 | 1.03 | partial | 79% | R108B 6.54 / R108B 10.20 (out of reach) / R107B 5.10 |
| tnf_sw_07 | tnf_ct_07 | 64 | 0.553 | 2 | 2 | 0.75 | yes | 79% | R108B 6.31 / R108B 6.93 |
| tnf_sw_08 | tnf_ct_08 | 80 | 0.489 | 2 | 2 | 0.62 | yes | 62% | R108B 8.30 / R108B 5.72 |
| tnf_sw_09 | tnf_ct_09 | 64 | 0.598 | 1 | 1 | 0.52 | yes | 82% | R108B 5.79 |
| tnf_sw_10 | tnf_ct_10 | 75 | 0.469 | 2 | 2 | 0.44 | yes | 81% | R108B 7.99 / R108B 7.03 |

## Matched controls

| name | pair | len | epitope residues | epitope cons. | note |
|---|---|---|---|---|---|
| tnf_ct_01 | tnf_sw_01 | 73 | 28 | 82% | histidine-free twin, identical repairs |
| tnf_ct_02 | tnf_sw_02 | 73 | 31 | 81% | histidine-free twin, identical repairs |
| tnf_ct_03 | tnf_sw_03 | 70 | 24 | 79% | histidine-free twin, identical repairs |
| tnf_ct_04 | tnf_sw_04 | 79 | 26 | 73% | histidine-free twin, identical repairs |
| tnf_ct_05 | tnf_sw_05 | 75 | 23 | 83% | histidine-free twin, identical repairs |
| tnf_ct_06 | tnf_sw_06 | 63 | 24 | 79% | histidine-free twin, identical repairs |
| tnf_ct_07 | tnf_sw_07 | 64 | 28 | 79% | histidine-free twin, identical repairs |
| tnf_ct_08 | tnf_sw_08 | 80 | 13 | 62% | histidine-free twin, identical repairs |
| tnf_ct_09 | tnf_sw_09 | 64 | 34 | 82% | histidine-free twin, identical repairs |
| tnf_ct_10 | tnf_sw_10 | 75 | 26 | 81% | histidine-free twin, identical repairs |

All 20 sequences are free of cysteines and contain no homopolymer run longer than 4.
