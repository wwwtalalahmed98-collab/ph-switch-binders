# Switch run 5: the repair fix, the widened parent pool, and the final set

Run: `tnf_switch5_repair_fix`, Kaggle, GPU T4 x2, 41m 24s, completed (both workers
reported done; all 29 parents appear in the summary table, so nothing was cut by the
time budget). Script: [kaggle/cellTNF19_switch_repair_fix.py](../../kaggle/cellTNF19_switch_repair_fix.py).

## What this run changed

Two fixes over `tnf_switch4_mouse_aware`:

1. **The liability repair became revertible and charge-preserving.** Free cysteines are
   still always replaced, but breaking up repetitive stretches is now tried, scored, and
   **reverted if it costs more than 15% of the parent's human i_ptm**. Substitutions
   preserve charge (E↔D, K↔R) instead of neutralising it (E→Q).
2. **The parent pool widened from human i_ptm ≥ 0.50 to ≥ 0.45**, adding 9 lineages the
   earlier screen had never examined, and each unscreened parent is measured against
   mouse inline before its variants are judged.

A third change was a correctness fix, not an improvement: the repair plan is now decided
**once per parent, on the control**, and the identical positional plan is applied to
every switch in that lineage, with an assertion that each switch differs from its control
at the histidine positions and nowhere else. Previously the plan was re-derived per
sequence, so a histidine landing inside a broken repetitive stretch could silently
desynchronise the switch/control twins.

## Result

> **29 parents offered histidines; 35 placements hold BOTH species; 14 parents yield at
> least one usable switch.**

Against switch4's 12 parents on the same bar, that is 12 → 14. The three fully new
lineages all came from the widened pool: `tnf3_w1_010`, `tnf3_w0_007`, `tnf_w1_004`.

## The repair fix did not do what it was built to do

The stated reason for this run was to rescue `tnf3_w1_019`, which lost a third of its
human i_ptm to the repair in switch4. It did not work:

| step | human i_ptm | note |
|---|---|---|
| parent | 0.714 | |
| after `C41S` alone | 0.463 | 65% of parent — the cysteine fix alone did the damage |
| after run-breaking too | 0.462 | reverted by the new rule, correctly, but it changed nothing |

My assumption had been that the forced `E→Q` substitutions were the expensive part. They
were not. A single free-cysteine replacement cost 35% of the binding, and that
substitution is not optional — a lone cysteine is a real manufacturing liability. All
four of `tnf3_w1_019`'s histidine placements then fell below the 70% human retention bar,
so the lineage is out. **The fix was correct and the hypothesis behind it was wrong.**

`tnf3_w1_016` went the same way for the same reason, more severely: six free cysteines,
human i_ptm 0.469 → 0.177 (38% of parent).

The revertible repair did earn its place elsewhere — it held for `tnf_w0_008`,
`tnf2_w1_001`, `tnf2_w0_001`, `tnf3_w0_012`, `tnf3_w0_007`, `tnf3_w1_010` and
`tnf_w1_004`, where run-breaking was free or better than the parent.

## `tnf2_w0_001`, resolved

This lineage sat in the proposed set on numbers from an older run and had never been
measured against mouse. It is now measured: **parent mouse i_ptm 0.679 with 12 on-site
contacts**, and its best placement `Q40H,Y37H` gives human 0.484 / mouse 0.652 / 12
on-site. It keeps its slot on verified numbers.

## The selection rule used for the final set

Objectives are ranked switch > mouse > affinity, so the rule is applied in that order:

1. **Prefer two histidines over one.** At histidine pKa 6.0, one His leaves 50% of
   molecules uncharged and still fully binding at pH 6.0; two leaves 25%. The brief asks
   for *no detectable binding* at pH 6.0.
2. **Prefer more histidines facing a conserved anchor** — a cation present in both human
   and mouse (K166A/R108B human, K168A/R111B mouse). A histidine facing a human-only
   cation cannot push off the mouse target, so the switch would be human-only.
3. **Prefer the higher worse-species i_ptm**, `min(human, mouse)`.

Then the ten highest-ranked parents take the ten slots. Built and verified by
[tnf/build_submission_v3.py](../build_submission_v3.py), which refuses any switch that
differs from its control anywhere other than the named histidine positions.

## The proposed set

| slot | lineage | histidines | nH | conserved anchors | human | mouse | worse | mouse-site |
|---|---|---|---|---|---|---|---|---|
| sw_01 | tnf3_w0_012 | E54H,W57H | 2 | 2 | 0.721 | 0.669 | 0.669 | 16 |
| sw_02 | tnf_w0_012 | Q23H,D27H | 2 | 2 | 0.622 | 0.591 | 0.591 | 12 |
| sw_03 | tnf2_w1_001 | L13H,D10H | 2 | 2 | 0.553 | 0.563 | 0.553 | 9 |
| sw_04 | tnf2_w0_001 | Q40H,Y37H | 2 | 2 | 0.484 | 0.652 | 0.484 | 12 |
| sw_05 | tnf3_w1_009 | I22H,L18H | 2 | 2 | 0.503 | 0.480 | 0.480 | 15 |
| sw_06 | tnf3_w0_017 | Y18H,E19H | 2 | 1 | 0.560 | 0.499 | 0.499 | 7 |
| sw_07 | tnf2_w0_008 | E2H | 1 | 1 | 0.598 | 0.628 | 0.598 | 10 |
| sw_08 | tnf2_w0_007 | L30H | 1 | 1 | 0.585 | 0.655 | 0.585 | 15 |
| sw_09 | tnf_w0_011 | L18H | 1 | 1 | 0.558 | 0.627 | 0.558 | 6 |
| sw_10 | tnf3_w0_007 | E55H | 1 | 1 | 0.505 | 0.495 | 0.495 | 11 |

Six of ten carry two histidines, against five in the previous proposal. Every switch
binds mouse at i_ptm ≥ 0.480 with at least 6 on-site mouse contacts.

Verified but outside the ten slots: `tnf_w0_007` Q46H (0.486/0.498/12), `tnf_w1_004`
Q45H (0.663/0.481/11), `tnf3_w0_016` L9H (0.479/0.452/8), `tnf3_w1_010` A30H
(0.452/0.521/12).

## Against what is currently banked

The banked submission has **six of ten switches failing the mouse objective**, which the
organisers rank above affinity. This set replaces five lineages outright, fixes one
placement, and keeps four.

| banked lineage | banked mouse i_ptm | mouse-site | outcome |
|---|---|---|---|
| tnf_w1_012 | 0.326 | 0 | dropped — no placement holds mouse |
| tnf_w0_003 | 0.353 | 7 | dropped — no placement holds mouse |
| tnf_w0_008 | 0.284 | 9 | dropped — no placement holds mouse |
| tnf_w0_004 | 0.189 | 6 | dropped — parent does not bind mouse |
| tnf_w1_004 | 0.251 | 7 | pair dropped; single `Q45H` works (0.481) but ranks 12th |
| tnf2_w0_007 | 0.148 | 0 | **kept, placement changed** to `L30H` → 0.655, 15 on-site |
| tnf_w0_012 | 0.591 | 12 | kept |
| tnf2_w1_001 | 0.563 | 9 | kept |
| tnf2_w0_008 | 0.628 | 10 | kept |
| tnf2_w0_001 | 0.652 | 12 | kept, now on measured mouse numbers |

## What is not done and not verified

- **Nothing has been physically tested.** Every number here is an AlphaFold2 prediction
  against a static crystal structure. i_ptm is a confidence score, not an affinity.
- **The pH switch itself has never been measured by anyone here.** It rests on a geometric
  argument — histidine placed within reach of a cation — plus textbook histidine
  chemistry. No calculation in this project models pH.
- **Two of the ten controls carry a native histidine** (`tnf2_w0_008` sw_07 and
  `tnf3_w0_007` sw_10). For those two pairs, a pH response in the lab cannot be
  attributed to our added histidine alone.
- **`sw_09` (`tnf_w0_011` L18H) has only 6 on-site mouse contacts**, the weakest in the
  set, and its parent's best-binding placement (`Q19H`, mouse 0.639) was passed over
  because that histidine faces a non-conserved cation and so would not switch in mouse.
  That is a judgement about an unmeasured mechanism, not a measurement.
- **The selection rule trades measured binding for an unvalidated switch argument.**
  Choosing pairs over singles cost real predicted binding in several lineages — `sw_05`
  gives up 0.561→0.503 human and 0.587→0.480 mouse versus its own single `I22H`. If the
  two-histidine reasoning is wrong, those trades were losses.
- **`switch_score`, the geometric switch heuristic, was not used as a selection
  criterion** beyond the conserved-anchor count. Several chosen variants score low on it
  (`sw_01` 0.351, `sw_10` 0.206) because the highest-scoring placements usually broke
  binding. The set is therefore optimised for binding-while-carrying-histidines, not for
  switch strength directly.
- **Kaggle's output preview silently truncates.** The result file is 113,975 bytes and
  the preview returns 102,400, closing the JSON brackets so it parses as valid. Reading
  it that way drops 3 of 29 records and corrupts a 4th without any error. The three lost
  records and the corrupted one all have zero passing placements, confirmed against the
  log's own complete table, so this set is unaffected — but any future read of a Kaggle
  output file has to check `totalBytes` against `previewBytes`.
