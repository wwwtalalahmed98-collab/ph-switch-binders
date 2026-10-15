# The parents bind mouse. The histidines are what break it.

Source: `tnf_mouse_parents2`, GPU T4 x2, 13.9 min. 89 parent designs across all three
rounds; the 31 that clear human i_ptm >= 0.50 with 5+ on-site contacts were predicted
against the mouse target, before any histidine is placed.

**20 of 31 parents bind mouse** at i_ptm >= 0.45 with 5+ mouse-site contacts. That is a
far larger pool than the 7 switches that passed the same bar, and it says the design
process was never the problem.

## The switch stage is what loses mouse binding

For the 17 lineages that have both a parent and a switch measured:

| lineage | parent h/m/site | switch h/m/site | mouse kept |
|---|---|---|---|
| tnf3_w0_012 | 0.720/0.661/15 | 0.722/0.669/16 | **101%** |
| tnf2_w1_001 | 0.595/0.513/11 | 0.553/0.566/9  | **110%** |
| tnf3_w0_017 | 0.608/0.529/10 | 0.561/0.499/7  | **94%** |
| tnf_w0_011  | 0.638/0.701/6  | 0.651/0.644/4  | **92%** |
| tnf2_w0_008 | 0.642/0.683/11 | 0.598/0.628/10 | **92%** |
| tnf_w0_012  | 0.691/0.653/13 | 0.620/0.592/12 | **91%** |
| tnf3_w0_016 | 0.599/0.631/9  | 0.289/0.441/6  | 70% |
| tnf_w0_008  | 0.682/0.560/8  | 0.710/0.284/9  | 51% |
| tnf_w1_004  | 0.632/0.455/12 | 0.493/0.222/7  | 49% |
| tnf_w0_003  | 0.717/0.726/12 | 0.650/0.351/7  | 48% |
| tnf_w1_012  | 0.758/0.701/13 | 0.744/0.326/0  | 47% |
| tnf3_w0_022 | 0.688/0.577/12 | 0.581/0.262/8  | 45% |
| tnf_w0_004  | 0.597/0.441/15 | 0.487/0.189/6  | 43% |
| tnf3_w1_006 | 0.658/0.609/6  | 0.202/0.245/2  | 40% |
| tnf_w0_007  | 0.504/0.465/12 | 0.284/0.174/7  | 37% |
| tnf3_w0_023 | 0.644/0.480/10 | 0.118/0.144/0  | 30% |
| tnf2_w0_007 | 0.586/0.719/14 | 0.487/0.148/0  | 21% |

The summary number:

    the same histidine substitutions keep a median 90% of HUMAN binding
    and a median 49% of MOUSE binding

That is not a coincidence, it is the selection criterion showing through. Histidine
positions were chosen on the human structure and kept or rejected on human i_ptm
retention alone. Mouse was never consulted, so positions that are tolerable in human and
destructive in mouse survive the filter by construction.

`tnf2_w0_007` is the clearest case: the parent is one of the best mouse binders measured
(0.719, 14 mouse-site contacts) and the two histidines take it to 0.148 with zero
mouse-site contacts, while keeping 83% of its human binding. It was ranked first in the
submission by the conservation proxy.

## Four mouse-compatible lineages have never been through a switch stage

`tnf3_w0_004`, `tnf3_w0_005`, `tnf3_w1_009`, `tnf3_w1_019` - all round 3, all clearing
both bars as parents, none ever offered a histidine.

`tnf3_w1_019` is the strongest of them: human 0.714, mouse 0.699, 98% retained.

## What this does not say

The six lineages whose histidines preserved mouse binding were not selected for that -
nothing in the pipeline was looking. Six out of seventeen is consistent with chance, so
the right reading is that mouse retention is currently unmanaged rather than that these
six designs are special. The next switch stage should place on human and then keep or
reject on retention against both species, which is a criterion no run so far has applied.
