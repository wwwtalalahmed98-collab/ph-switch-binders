# Round 3 — switch placement at a histidine cap of 2, and confirmation

Source: `tnf_finish3` on Kaggle (GPU T4 x2), run from `kaggle/cellTNF15_switch3.py`
and `kaggle/cellTNF16_confirm3.py`. Inputs: the 23 round-3 designs that cleared both
bars in `round3_summary.csv`.

## What came out

30 records: 14 switch/control pairs plus two designs (`tnf3_w0_005`, `tnf3_w0_009`)
where no histidine site satisfied the geometric filter, so only the control twin
exists. Five switches cleared the retention and i_ptm bars and went to confirmation.

Every i_ptm below is from an independent re-prediction at 3 recycles, never the
optimiser's own score.

| design | len | His | confirmed i_ptm | control i_ptm | on-site | epitope kept | His geometry (Cb -> cation N) |
|---|---|---|---|---|---|---|---|
| tnf3_w0_012_sw | 76 | E54H, W57H | 0.722 | 0.702 | 16 | 96% | 54 -> R108B 7.21 A; 57 -> R108B 8.03 A |
| tnf3_w1_014_sw | 76 | L22H       | 0.721 | 0.720 |  8 | 88% | 22 -> R108B 7.05 A |
| tnf3_w1_001_sw | 68 | L45H, E49H | 0.616 | 0.638 |  6 | 91% | 45 -> R108B 4.39 A; 49 -> R108B 11.34 A |
| tnf3_w0_022_sw | 74 | D18H, I21H | 0.581 | 0.671 | 12 | 85% | 18 -> R108B 4.57 A; 21 -> R108B 6.01 A |
| tnf3_w0_017_sw | 56 | Y18H, E19H | 0.561 | 0.597 |  3 | 74% | 18 -> R107B 2.60 A; 19 -> R107B 8.41 A |

## What is actually good here

`tnf3_w0_012_sw` is the strongest switch produced in any round: confirmed i_ptm 0.722,
16 on-site contacts, 96% of the parent epitope retained, and both histidines inside the
3.0-8.8 A window facing R108 — which is arginine in mouse too (human R108 aligns to
mouse R111, one of the two anchors). Its control sits at 0.702 with 14 on-site contacts,
so the pair is well matched and a measured pH response would be interpretable.

`tnf3_w1_014_sw` carries a single histidine and binds at 0.721 against a control at
0.720. One histidine buys a smaller predicted pH swing, but it costs essentially nothing
in binding, which is the trade the cap-2 experiment was looking for.

## Three things the table shows that are not improvements

- **`tnf3_w1_001_sw` has one working histidine, not two.** His49 measures 11.34 A from
  the cation in the re-predicted structure, well outside the 3.0-8.8 A window. Placement
  used the parent structure's geometry; the backbone moved. It should be read as a
  one-histidine design.
- **`tnf3_w0_017_sw` is the weakest of the five and is aimed at the wrong residue.** Its
  partner is R107, which is glutamine in mouse — the same cross-reactivity liability
  found in round 1. It also has only 3 on-site contacts, 74% epitope retention, and a
  2.60 A Cb-to-cation distance that is below the lower bound of the intended window.
- **The round-3 hit rate is not evidence of a better method.** 20 of 46 designs reached
  i_ptm >= 0.50 against 9 of 27 across rounds 1 and 2. Settings were identical apart from
  the random seed; the run was simply longer. Those counts are not separable at these
  sample sizes.

## Still to do

Mouse binding for every cap-2 design from all three rounds is being measured in
`tnf_mouse2`, which replaces the sequence-conservation proxy with a predicted mouse
i_ptm. Re-ranking waits on that.
