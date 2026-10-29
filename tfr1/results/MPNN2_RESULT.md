# Rewriting run 2's well-aimed backbones (MPNN pass 2)

Kaggle `notebookdc9683b1fc` version 2, CPU only, 321 s, no GPU quota spent.
Script: [kaggle/cellTFR1j_mpnn2.py](../../kaggle/cellTFR1j_mpnn2.py).
Candidates: [mpnn2_candidates.tsv](mpnn2_candidates.tsv) (41 rows, copied from the log).

## What it did

10 of 36 run-2 backbones reach 3 or more epitope contacts. Each was handed to
ProteinMPNN with the backbone fixed: 24 sequences at temperature 0.1 and 24 at 0.2,
480 sequences in all. Selection was on composition, not on MPNN's own score.

Nine of the ten yielded something; `e559_005` yielded nothing (0 good, 0 marginal
out of 48). 41 candidates survived, 12 of them tier "good".

## The headline the run printed, and why it is wrong

> `12 sequences are BOTH on a well-aimed backbone and composition-good.`
> `That is the first time in this project.`

**That number does not survive a third screen.** The two-sided screen checks charge,
hydrophobicity, Kyte-Doolittle and identical-residue runs. It does not look at
proline. Natural globular proteins are about 5% proline; four of the twelve "good"
sequences are 25-32% proline, which is a polyproline rod, not a folded protein.

Re-screened with a proline/glycine gate ([rescreen_proline.py](../rescreen_proline.py),
P <= 12%, G <= 12%, P+G <= 20%):

| | count |
|---|---|
| candidates | 41 |
| survive the proline/glycine gate | 3 |
| of those, also tier "good" | **1** |

The survivor is **`tfr1b_e559_003_r02`**: 8% proline, 6% glycine, 35% charged,
40% hydrophobic, Kyte-Doolittle -0.65, on a backbone with 4/13 on-site contacts.

```
SPTKTYVSIEKTKSKGGKGPKPKPKGRKKLVEKLAERLAEREAARAAAEAAARAAALAAIPA
```

The two marginals that also pass are `tfr1b_e559_009_r05` and `tfr1b_e559_006_r01`,
both with Kyte-Doolittle near zero, i.e. on the oily edge.

Several candidates sit at 11-13% proline, a hair either side of the 12% line. That
line is my choice, not a measured threshold, so the count of 3 is soft in the
direction of being slightly generous or slightly harsh by one or two rows.

## This is the same mistake twice

Yesterday the screen had a ceiling on charge and a floor on hydrophobicity and
nothing on the other side, so it called a 64% hydrophobic aggregate "good". Today it
had no proline term at all, so it called a polyproline rod "good". Both times the
screen was written to catch the failure mode in front of me and silently passed the
next one. The gate is now three-sided; it is still not a fold prediction.

## What is confirmed, and what is not

Confirmed by running:

- MPNN can put a plausible composition onto a well-aimed backbone. One sequence out
  of 480 sampled, which is a low yield but not zero.
- The transcription of the candidate table is sound: `rescreen_proline.py` recomputes
  Kyte-Doolittle from the sequences and matches the Kaggle log on every row.

Not verified, and this is the important part:

- **No rewrite has been folded.** `on_site_parent` is the *parent backbone's* contact
  count. Whether `003_r02` still binds where its parent bound is unmeasured, and MPNN
  routinely produces sequences that fold to something other than the input backbone.
- The parent `tfr1b_e559_003` has i_ptm 0.559 and pLDDT 0.499. That is weak evidence
  it is a real fold at all.
- Run 2's designs have not been checked against the transferrin/HFE site. That check
  was run on run 1 only, where 3 of 47 blocked it.
- No acid-release switch has been measured on anything.
- Nothing is submitted.
