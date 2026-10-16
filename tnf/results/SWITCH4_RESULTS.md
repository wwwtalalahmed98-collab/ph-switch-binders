# Mouse-aware histidine placement

Source: `tnf_switch4_mouse_aware`, GPU T4 x2, 30.2 min. 20 parents - every lineage the
screen measured as binding mouse - offered up to four histidine placements each, every
placement that held human then measured against mouse.

**20 parents, 68 placements tried, 30 hold both species, 12 parents yield at least one
usable switch.**

## The result that shows the mechanism

`tnf2_w0_007` is the same parent, the same engine, the same constants. Only the chosen
position differs.

| placement | human | mouse | mouse-site contacts |
|---|---|---|---|
| L30H,D27H (what every previous run built) | 0.487 | **0.150** | **0** |
| D27H alone | 0.502 | 0.146 | 0 |
| L30H alone | 0.585 | 0.656 | 15 |
| **E26H (kept)** | **0.566** | **0.678** | **15** |

A parent that measured 0.719 against mouse was being handed to the assay at 0.150 with
no mouse-site contact at all, because the position that scored best geometrically
happened to be one mouse could not tolerate. Three residues away is fine. Nothing about
the lineage was wrong.

`tnf_w0_012` shows the same thing more mildly: the old pair Q23H,D27H gives mouse 0.592;
A28H gives 0.638 at a HIGHER human i_ptm (0.710 against 0.691 for the parent itself).

## Best placement per parent

| parent | his | human | kept | mouse | kept | mouse-site |
|---|---|---|---|---|---|---|
| tnf3_w0_012 | E54H,W57H | 0.722 | 100% | 0.669 | 101% | 16 |
| tnf_w0_012  | A28H | 0.710 | 103% | 0.638 |  98% | 13 |
| tnf_w0_011  | Q19H | 0.651 | 102% | 0.638 |  91% |  5 |
| tnf3_w0_017 | Y18H | 0.606 | 100% | 0.517 |  98% |  8 |
| tnf2_w0_008 | E2H  | 0.598 |  93% | 0.628 |  92% | 10 |
| tnf2_w1_001 | L13H | 0.576 |  97% | 0.576 | 112% | 11 |
| tnf2_w0_007 | E26H | 0.566 |  97% | 0.678 |  94% | 15 |
| tnf3_w1_009 | I22H | 0.558 |  98% | 0.589 |  93% | 16 |
| tnf_w1_004  | Q45H | 0.650 | 103% | 0.425 |  93% | 12 |
| tnf_w0_007  | Q46H | 0.486 |  96% | 0.498 | 107% | 12 |
| tnf3_w0_016 | L9H  | 0.476 |  80% | 0.452 |  72% |  8 |
| tnf3_w0_004 | Q32H | 0.442 |  77% | 0.421 |  77% | 10 |

## How much of this is the new criterion

On the identical absolute bar used before (human i_ptm >= 0.50, mouse >= 0.45, 5+
mouse-site contacts), the old switch stages produced **5** usable switches and this one
produces **8**. That is the honest comparison; the "12 parents" headline uses a
retention-based bar that the earlier runs were never scored against.

5 to 8 is a modest count. Two changes were made at once - parents selected on measured
mouse binding, and several placements tried per parent instead of one - so the count
cannot be attributed to the acceptance criterion alone. The `tnf2_w0_007` rescue is the
stronger evidence, because there the parent, the engine and the constants are all held
fixed and only the position changes.

## Four things this run makes visible that are not improvements

- **Eleven of the twelve kept placements are single histidines.** That is the criterion
  doing what it was told: one histidine disturbs less, so it retains more. But one
  histidine is a weaker switch. At pKa 6.0 a single His goes from 3.8% protonated at pH
  7.4 to 50% at pH 6.0; two independent histidines give 7.5% to 75% for at least one
  charged. Retention was optimised and switch magnitude was not, so the designs are
  likelier to bind at both pH values than the two-histidine designs were. This is a real
  trade the run created and did not measure.

- **`tnf_w1_012` cannot be switched at all.** Its parent binds mouse at 0.701 with 13
  mouse-site contacts, and all four placements drop it to 0.206-0.327 with **zero**
  mouse-site contacts. It is the highest-affinity design in the current submission
  (tnf_sw_02, human 0.744). Some lineages appear incompatible with a histidine switch
  that survives in mouse, and no amount of reposition fixed this one.

- **The liability repair is now visibly expensive on some backbones.** `tnf3_w1_019`
  lost a third of its human i_ptm to repair alone (0.714 to 0.472 from C41S plus four
  E->Q forced by the homopolymer-run rule) and its repaired control's mouse i_ptm
  collapsed to 0.193, 28% of the parent. The switch stage never got a fair chance at
  that lineage. Same pattern, smaller, on `tnf3_w1_006`, `tnf3_w0_016` and `tnf3_w0_004`.
  This is a repair problem, not a switch problem, and nothing has been done about it.

- **`tnf3_w0_005` has no acid-OFF site at all** - no position in the binder passes the
  geometric filter against any exposed cation. It is a good mouse binder (0.652) with no
  way to make it pH-sensitive by this method.
